#' Run all statistics defined by the YAML `stats_report:` block
#'
#' Iterates over `comparisons:`, `correlations:`, `linear_models:` and
#' `ion_ratios:` and returns tidy data frames suitable for writing to xlsx
#' tabs.
#'
#' Each comparison names its own test, because the choice has to be made
#' before the p-values are seen:
#' \itemize{
#'   \item 2 levels  -> `welch` (default), `student`, or `wilcoxon`
#'   \item 3+ levels -> `anova` (default, with Tukey HSD post-hoc) or
#'     `kruskal` (with Holm-adjusted pairwise Wilcoxon post-hoc)
#' }
#' A comparison with no `method:` falls back to the default and says so, so
#' existing configurations keep running.
#'
#' Add `paired: true` together with `subject:` naming the `sample_meta` column
#' that links the two observations of each subject - pre/post, matched tissues,
#' the same animal at two timepoints. Only subjects present in both groups
#' contribute, and the effect size becomes Cohen's dz (or the matched-pairs
#' rank-biserial correlation for `wilcoxon`), which is *not* comparable to the
#' between-group values and is labelled distinctly so the two are never pooled.
#' The method label gains a `_paired` suffix for the same reason.
#'
#' Multiple-testing correction is applied *within families*, not across the
#' whole table: one family per comparison and method, per correlation analysis
#' and subset and method, and per linear-model term. The family used for each
#' row is recorded in its `adjust_family` column. Adjusting across unrelated
#' analyses would make the size of the correction depend on how many other
#' things happened to be configured in the same YAML.
#'
#' @param de A `struct::DatasetExperiment` (from `runMRManalyzeR()` or `readRDS()`).
#' @param st_params The `stats_report:` block from the YAML (a list).
#' @return A list with elements `stats`, `correlations`, `linear_models`,
#'   each a data frame (possibly with 0 rows if the corresponding YAML
#'   block was empty).
#' @examples
#' m  <- data.frame(A = c(10, 12, 11, 20, 24, 22), B = c(5, 6, 5, 11, 13, 12),
#'                  row.names = paste0("S", 1:6))
#' fm <- data.frame(Compound = c("A", "B"), row.names = c("A", "B"))
#' sm <- data.frame(Sample_ID = paste0("S", 1:6),
#'                  Group = rep(c("ctrl", "trt"), each = 3),
#'                  row.names = paste0("S", 1:6))
#' de <- struct::DatasetExperiment(data = m, sample_meta = sm, variable_meta = fm)
#' params <- list(comparisons = list(enabled = TRUE, entries = list(
#'   list(name = "ctrl_vs_trt", method = "welch",
#'        compare = list(factor = "Group", levels = c("ctrl", "trt"))))))
#' runStats(de, params)$stats
#' @family stats
#' @export
runStats = function(de, st_params){

  comparisons   = .section_entries(st_params$comparisons,   "comparisons")
  correlations  = .section_entries(st_params$correlations,  "correlations")
  linear_models = .section_entries(st_params$linear_models, "linear_models")

  # p-value adjustment is hard-coded to BH (= FDR). Bonferroni / others
  # were dropped from the YAML for simplicity.
  p_adjust      = "BH"
  sig_threshold = st_params$sig_threshold %||% 0.05

  ion_entries = .ion_entries(st_params$ion_ratios)

  list(
    stats         = .run_comparisons(de, comparisons,   p_adjust, sig_threshold),
    correlations  = .run_correlations(de, correlations, p_adjust, sig_threshold),
    linear_models = .run_linear_models(de, linear_models, p_adjust, sig_threshold),
    ion_ratios    = .run_ion_ratios(de, ion_entries, comparisons, p_adjust)
  )
}


# ---- Comparisons -----------------------------------------------------------

#' Adjust p-values within families rather than across the whole table
#'
#' Multiple-testing correction only means something relative to a stated
#' family. Adjusting one table that mixes several comparisons, both a
#' parametric and a non-parametric test of the same hypothesis, and global
#' tests alongside their own post-hoc results, controls nothing anybody asked
#' about: the size of the correction then depends on how many unrelated
#' analyses happened to be configured in the same YAML.
#'
#' So the family is made explicit and recorded in the returned table. Each
#' family is corrected across features independently.
#'
#' @param df A results table with a `p_value` column.
#' @param keys Column names whose combination defines a family.
#' @param p_adjust Method passed to [stats::p.adjust()].
#' @param sig_threshold Threshold applied to the adjusted value.
#' @return `df` with `adjust_family`, `p_adj` and `significant` filled in.
#' @keywords internal
#' @noRd
.adjust_within = function(df, keys, p_adjust, sig_threshold){
  if(nrow(df) == 0) return(df)
  fam = do.call(paste, c(df[keys], sep = " | "))
  df$adjust_family = fam
  df$p_adj = stats::ave(df$p_value, fam, FUN = function(p)
    stats::p.adjust(p, method = p_adjust))
  df$significant = df$p_adj < sig_threshold
  df
}

#' @keywords internal
#' @noRd
.run_comparisons = function(de, comparisons, p_adjust, sig_threshold){

  if(length(comparisons) == 0)
    return(.empty_stats_df())

  rows = lapply(comparisons, function(comp){
    .run_one_comparison(de, comp, p_adjust, sig_threshold)
  })

  out = dplyr::bind_rows(rows)
  # For `tukey` rows p_value is ALREADY adjusted across the pairwise
  # comparisons within a feature, by TukeyHSD() itself. p_adj then controls
  # the false discovery rate across features. Those are two different
  # families, so this is not a double correction of one family - but it is
  # conservative, because BH assumes p-values uniform under the null and
  # Tukey-adjusted values are not. Read p_adj on post-hoc rows accordingly.
  out = .adjust_within(out, c("comparison", "method"), p_adjust,
                       sig_threshold)
  out
}

#' Resolve the test method for one comparison
#'
#' The method must be decided before the p-values are seen. Running a
#' parametric and a non-parametric test of the same hypothesis and reporting
#' both lets the reader choose whichever came out lower, which is a forking
#' path rather than two pieces of evidence - so a comparison names one.
#'
#' `"welch"` is the default two-group test because [stats::t.test()] does not
#' assume equal variances and group sizes in a targeted study are usually
#' unequal; `"student"` requests the equal-variance form explicitly.
#'
#' @param comp A comparison entry; `comp$method` is read.
#' @param ngrp Number of non-empty levels.
#' @param comp_name Used in messages.
#' @return One of `"welch"`, `"student"`, `"wilcoxon"` (two groups) or
#'   `"anova"`, `"kruskal"` (three or more).
#' Hedges' g for two independent groups
#'
#' Cohen's d on the pooled standard deviation, times the small-sample
#' correction J. Reported for both t-test variants: the correction matters at
#' the group sizes a targeted study runs, where d is biased upwards.
#'
#' Sign follows `log2FC`: positive means group 2 is higher.
#'
#' @param x1,x2 Numeric vectors, group 1 and group 2.
#' @return Numeric, or `NA` when either group has fewer than two values.
#' @keywords internal
#' @noRd
.hedges_g = function(x1, x2){
  x1 = x1[!is.na(x1)]; x2 = x2[!is.na(x2)]
  n1 = length(x1); n2 = length(x2)
  if(n1 < 2 || n2 < 2) return(NA_real_)
  s_pool = sqrt(((n1 - 1) * stats::var(x1) + (n2 - 1) * stats::var(x2)) /
                (n1 + n2 - 2))
  if(!is.finite(s_pool) || s_pool == 0) return(NA_real_)
  d = (mean(x2) - mean(x1)) / s_pool
  d * (1 - 3 / (4 * (n1 + n2) - 9))          # Hedges' correction
}

#' Cohen's dz for paired observations
#'
#' The mean of the within-subject differences over their standard deviation.
#' Deliberately not comparable to Hedges' g: dz depends on how correlated the
#' pairs are, so a paired and an unpaired effect size must not be pooled.
#'
#' @param x1,x2 Paired numeric vectors, aligned subject by subject.
#' @return Numeric, or `NA`.
#' @keywords internal
#' @noRd
.cohens_dz = function(x1, x2){
  d = x2 - x1
  d = d[!is.na(d)]
  if(length(d) < 2) return(NA_real_)
  sdd = stats::sd(d)
  if(!is.finite(sdd) || sdd == 0) return(NA_real_)
  mean(d) / sdd
}

#' Matched-pairs rank-biserial correlation
#'
#' From the signed-rank statistic V: the difference between the proportion of
#' rank mass favouring each direction. Signed like `log2FC`.
#'
#' @param V The `statistic` from a paired [stats::wilcox.test()].
#' @param n Number of pairs.
#' @return Numeric, or `NA`.
#' @keywords internal
#' @noRd
.rank_biserial_paired = function(V, n){
  if(is.null(V) || is.na(V) || n < 1) return(NA_real_)
  1 - 4 * unname(V) / (n * (n + 1))
}

#' Rank-biserial correlation from a Wilcoxon rank-sum statistic
#'
#' The non-parametric companion to Hedges' g: the difference between the
#' proportion of cross-group pairs favouring each group, bounded -1 to 1.
#' Signed like `log2FC`, so positive means group 2 is higher.
#'
#' @param W The `statistic` from [stats::wilcox.test()] (group 1 first).
#' @param n1,n2 Group sizes used in the test.
#' @return Numeric, or `NA`.
#' @keywords internal
#' @noRd
.rank_biserial = function(W, n1, n2){
  if(is.null(W) || is.na(W) || n1 < 1 || n2 < 1) return(NA_real_)
  1 - 2 * unname(W) / (n1 * n2)
}

#' @keywords internal
#' @noRd
.resolve_method = function(comp, ngrp, comp_name){
  allowed = if(ngrp == 2) c("welch", "student", "wilcoxon")
            else          c("anova", "kruskal")
  m = comp$method

  if(is.null(m)){
    m = allowed[1]
    message(sprintf(
      "[runStats] comparison '%s' has no method:; defaulting to '%s'. Set method: in the YAML to pre-specify the test.",
      comp_name, m))
    return(m)
  }

  m = tolower(as.character(m)[1])
  if(identical(m, "t_test") || identical(m, "t-test")) m = "welch"
  if(!m %in% allowed)
    stop(sprintf(
      "[runStats] comparison '%s': method '%s' is not valid for %d groups. Use one of: %s.",
      comp_name, m, ngrp, paste(allowed, collapse = ", ")))
  m
}

#' @keywords internal
#' @noRd
.run_one_comparison = function(de, comp, p_adjust, sig_threshold){

  comp_name = comp$name %||% "unnamed_comparison"
  factor_name = comp$compare$factor
  levels_keep = comp$compare$levels
  if(is.null(factor_name) || is.null(levels_keep))
    stop(sprintf("Comparison '%s': compare.factor and compare.levels required.", comp_name))

  # Build a single combined subset: user-supplied conditions + the level
  # restriction. subsetDataset() supports vector values via %in%, so we can pass
  # `levels_keep` directly. Doing it in one call avoids mutating the DE in
  # place (which trips SummarizedExperiment's assay-replacement check).
  cond = comp$subset
  cond = if(is.null(cond)) list() else as.list(cond)

  if(!factor_name %in% colnames(de$sample_meta))
    stop(sprintf("Comparison '%s': factor '%s' not in sample_meta.",
                 comp_name, factor_name))
  cond[[factor_name]] = levels_keep
  de_sub = subsetDataset(de, conditions = cond)

  group = factor(de_sub$sample_meta[[factor_name]], levels = levels_keep)
  if(nlevels(droplevels(group)) < 2){
    warning(sprintf("Comparison '%s': fewer than 2 non-empty levels - skipped.",
                    comp_name))
    return(.empty_stats_df())
  }

  data = de_sub$data
  feats = colnames(data)
  ngrp  = nlevels(droplevels(group))

  method = .resolve_method(comp, ngrp, comp_name)
  paired = isTRUE(comp$paired)

  subject = NULL
  if(paired){
    if(ngrp != 2)
      stop(sprintf("Comparison '%s': paired analysis needs exactly 2 levels.",
                   comp_name))
    sub_col = comp$subject %||% stop(sprintf(
      "Comparison '%s': paired: true also needs subject: naming the sample_meta column that pairs the observations.",
      comp_name))
    if(!sub_col %in% colnames(de_sub$sample_meta))
      stop(sprintf("Comparison '%s': subject column '%s' not in sample_meta.",
                   comp_name, sub_col))
    subject = as.character(de_sub$sample_meta[[sub_col]])
  }

  if(ngrp == 2){
    .compare_two(comp_name, data, group, feats, levels_keep, method,
                 subject = subject)
  } else {
    .compare_many(comp_name, data, group, feats, levels_keep, method)
  }
}

#' @keywords internal
#' @noRd
.compare_two = function(comp_name, data, group, feats, levels_keep,
                        method = "welch", subject = NULL){
  g1 = levels_keep[1]; g2 = levels_keep[2]
  i1 = which(group == g1); i2 = which(group == g2)

  # Paired analysis needs the two vectors aligned subject by subject, and only
  # subjects measured in both conditions contribute. Dropping the unmatched
  # ones is the point of pairing: an unpaired test on the same data throws away
  # the very structure that makes it powerful.
  paired = !is.null(subject)
  if(paired){
    s1 = subject[i1]; s2 = subject[i2]
    if(anyDuplicated(s1) || anyDuplicated(s2))
      stop(sprintf("Comparison '%s': subject values are not unique within a group, so observations cannot be paired.",
                   comp_name))
    common = intersect(s1, s2)
    if(length(common) < 2)
      stop(sprintf("Comparison '%s': only %d subject(s) present in both groups.",
                   comp_name, length(common)))
    if(length(common) < length(s1) || length(common) < length(s2))
      message(sprintf(
        "[runStats] comparison '%s': %d subject(s) measured in both groups; unmatched observations dropped.",
        comp_name, length(common)))
    i1 = i1[match(common, s1)]
    i2 = i2[match(common, s2)]
  }

  rows = lapply(feats, function(f){
    x  = data[[f]]
    x1 = x[i1]; x2 = x[i2]
    n1 = sum(!is.na(x1)); n2 = sum(!is.na(x2))
    m1 = mean(x1, na.rm = TRUE); m2 = mean(x2, na.rm = TRUE)
    fc     = suppressWarnings(m2 / m1)
    log2fc = suppressWarnings(log2(fc))

    if(method == "wilcoxon"){
      # conf.int gives the Hodges-Lehmann location shift and its interval -
      # the non-parametric answer to "by how much", which a p-value alone
      # does not provide.
      wt = tryCatch(stats::wilcox.test(x1, x2, exact = FALSE, conf.int = TRUE,
                                       paired = paired),
                    error = function(e) NULL, warning = function(w) NULL)
      return(list(.row(
        comp_name, "wilcoxon", f, 2, g1, g2, m1, m2, fc, log2fc,
        statistic = if(is.null(wt)) NA_real_ else unname(wt$statistic),
        df        = NA_character_,
        p_value   = if(is.null(wt)) NA_real_ else wt$p.value,
        n_g1 = n1, n_g2 = n2,
        n_na_g1 = length(x1) - n1, n_na_g2 = length(x2) - n2,
        mean_diff = m2 - m1,
        conf_low  = if(is.null(wt$conf.int)) NA_real_ else -wt$conf.int[2],
        conf_high = if(is.null(wt$conf.int)) NA_real_ else -wt$conf.int[1],
        effect_size = if(paired)
          .rank_biserial_paired(if(is.null(wt)) NA else wt$statistic,
                                min(n1, n2))
        else
          .rank_biserial(if(is.null(wt)) NA else wt$statistic, n1, n2),
        effect_type = if(paired) "rank_biserial_paired"
                      else       "rank_biserial")))
    }

    tt = tryCatch(stats::t.test(x1, x2, var.equal = method == "student",
                                paired = paired),
                  error = function(e) NULL)
    list(.row(
      comp_name, if(paired) paste0(method, "_paired") else method,
      f, 2, g1, g2, m1, m2, fc, log2fc,
      statistic = if(is.null(tt)) NA_real_ else unname(tt$statistic),
      df        = if(is.null(tt)) NA_character_ else as.character(round(tt$parameter, 2)),
      p_value   = if(is.null(tt)) NA_real_ else tt$p.value,
      n_g1 = n1, n_g2 = n2,
      n_na_g1 = length(x1) - n1, n_na_g2 = length(x2) - n2,
      mean_diff = m2 - m1,
      # t.test() forms the interval as group1 - group2; negate so it is
      # signed like mean_diff and log2FC.
      conf_low  = if(is.null(tt)) NA_real_ else -tt$conf.int[2],
      conf_high = if(is.null(tt)) NA_real_ else -tt$conf.int[1],
      # Paired designs use Cohen's dz - the mean difference over the SD of the
      # differences - which is not comparable to the between-group g and is
      # named differently so the two are never averaged or ranked together.
      effect_size = if(paired) .cohens_dz(x1, x2) else .hedges_g(x1, x2),
      effect_type = if(paired) "cohens_dz" else "hedges_g"))
  })

  dplyr::bind_rows(unlist(rows, recursive = FALSE))
}

#' @keywords internal
#' @noRd
.compare_many = function(comp_name, data, group, feats, levels_keep,
                         method = "anova"){

  # Shared by both post-hoc paths: the per-pair descriptives.
  pair_stats = function(x, g_a, g_b){
    xa = x[group == g_a]; xb = x[group == g_b]
    list(xa = xa, xb = xb,
         na = sum(!is.na(xa)), nb = sum(!is.na(xb)),
         ma = mean(xa, na.rm = TRUE), mb = mean(xb, na.rm = TRUE))
  }

  rows = lapply(feats, function(f){
    x = data[[f]]
    df_in = data.frame(value = x, grp = group)
    n_tot = sum(!is.na(x))
    k     = length(levels_keep)

    if(method == "kruskal"){
      kw = tryCatch(stats::kruskal.test(value ~ grp, data = df_in),
                    error = function(e) NULL)
      # Epsilon-squared: H / (n - 1), the rank analogue of eta-squared.
      eps = if(is.null(kw) || n_tot < 2) NA_real_
            else unname(kw$statistic) / (n_tot - 1)

      global = list(.row(
        comp_name, "kruskal", f, k, NA, NA, NA, NA, NA, NA,
        statistic = if(is.null(kw)) NA_real_ else unname(kw$statistic),
        df        = if(is.null(kw)) NA_character_ else as.character(kw$parameter),
        p_value   = if(is.null(kw)) NA_real_ else kw$p.value,
        effect_size = eps, effect_type = "epsilon_squared"))

      # Post-hoc must match the global test, so pairwise Wilcoxon rather than
      # Tukey (which is built on the ANOVA fit). Done pair by pair rather than
      # via pairwise.wilcox.test() so each row can carry the same statistics,
      # interval and effect size as the two-group path; the p-values are then
      # Holm-adjusted across the pairs of this feature, exactly as
      # pairwise.wilcox.test() would.
      pairs = utils::combn(levels_keep, 2, simplify = FALSE)
      raw   = numeric(length(pairs))
      parts = vector("list", length(pairs))

      for(i in seq_along(pairs)){
        g_a = pairs[[i]][1]; g_b = pairs[[i]][2]
        ps  = pair_stats(x, g_a, g_b)
        wt  = tryCatch(stats::wilcox.test(ps$xa, ps$xb, exact = FALSE,
                                          conf.int = TRUE),
                       error = function(e) NULL, warning = function(w) NULL)
        raw[i]   = if(is.null(wt)) NA_real_ else wt$p.value
        parts[[i]] = list(g_a = g_a, g_b = g_b, ps = ps, wt = wt)
      }
      adj = stats::p.adjust(raw, method = "holm")

      ph_rows = lapply(seq_along(parts), function(i){
        q = parts[[i]]; ps = q$ps; wt = q$wt
        fc_pair = suppressWarnings(ps$mb / ps$ma)
        .row(comp_name, "pairwise_wilcoxon", f, 2, q$g_a, q$g_b,
             ps$ma, ps$mb, fc_pair, suppressWarnings(log2(fc_pair)),
             statistic = if(is.null(wt)) NA_real_ else unname(wt$statistic),
             df = NA_character_,
             p_value = adj[i],
             n_g1 = ps$na, n_g2 = ps$nb,
             n_na_g1 = length(ps$xa) - ps$na,
             n_na_g2 = length(ps$xb) - ps$nb,
             mean_diff = ps$mb - ps$ma,
             conf_low  = if(is.null(wt$conf.int)) NA_real_ else -wt$conf.int[2],
             conf_high = if(is.null(wt$conf.int)) NA_real_ else -wt$conf.int[1],
             effect_size = .rank_biserial(if(is.null(wt)) NA else wt$statistic,
                                          ps$na, ps$nb),
             effect_type = "rank_biserial")
      })
      return(c(global, ph_rows))
    }

    av = tryCatch(stats::aov(value ~ grp, data = df_in), error = function(e) NULL)
    sm = if(is.null(av)) NULL else summary(av)[[1]]
    av_p = if(is.null(sm)) NA_real_ else sm[["Pr(>F)"]][1]
    av_f = if(is.null(sm)) NA_real_ else sm[["F value"]][1]
    av_df = if(is.null(sm)) NA_character_ else
      sprintf("%d,%d", sm[["Df"]][1], sm[["Df"]][2])
    # Eta-squared: the share of total variance the grouping accounts for.
    eta = if(is.null(sm)) NA_real_ else
      sm[["Sum Sq"]][1] / sum(sm[["Sum Sq"]])

    global = list(.row(
      comp_name, "anova", f, k, NA, NA, NA, NA, NA, NA,
      statistic = av_f, df = av_df, p_value = av_p,
      effect_size = eta, effect_type = "eta_squared"))

    # Tukey HSD post-hoc on the aov fit. Its p-values are already adjusted
    # across the pairs within this feature, and it supplies the interval.
    tk_rows = list()
    if(!is.null(av)){
      tk = tryCatch(stats::TukeyHSD(av), error = function(e) NULL)
      if(!is.null(tk)){
        tk_df = as.data.frame(tk$grp)
        tk_df$pair = rownames(tk_df)
        tk_rows = lapply(seq_len(nrow(tk_df)), function(i){
          pr = strsplit(tk_df$pair[i], "-", fixed = TRUE)[[1]]
          g_a = pr[2]; g_b = pr[1]                  # TukeyHSD label = "B-A"
          ps  = pair_stats(x, g_a, g_b)
          fc_pair = suppressWarnings(ps$mb / ps$ma)
          .row(comp_name, "tukey", f, 2, g_a, g_b, ps$ma, ps$mb,
               fc_pair, suppressWarnings(log2(fc_pair)),
               statistic = tk_df$diff[i],
               df = NA_character_,
               p_value = tk_df$`p adj`[i],
               n_g1 = ps$na, n_g2 = ps$nb,
               n_na_g1 = length(ps$xa) - ps$na,
               n_na_g2 = length(ps$xb) - ps$nb,
               mean_diff = tk_df$diff[i],
               conf_low = tk_df$lwr[i], conf_high = tk_df$upr[i],
               effect_size = .hedges_g(ps$xa, ps$xb),
               effect_type = "hedges_g")
        })
      }
    }

    c(global, tk_rows)
  })

  dplyr::bind_rows(unlist(rows, recursive = FALSE))
}

# ---- Correlations (pairwise feature x feature) -----------------------------

#' @keywords internal
#' @noRd
.run_correlations = function(de, correlations, p_adjust, sig_threshold){
  if(length(correlations) == 0)
    return(.empty_corr_df())

  rows = lapply(correlations, function(corr){
    name    = corr$name    %||% "unnamed_correlation"
    methods = corr$methods %||% c("pearson", "spearman")

    # Accept either `subsets:` (plural list) or `subset:` (singular, legacy).
    subsets = corr$subsets
    if(is.null(subsets)){
      if(is.null(corr$subset)) subsets = list(list())   # one analysis on all samples
      else                     subsets = list(as.list(corr$subset))
    }

    include_self = isTRUE(corr$include_self %||% FALSE)

    per_subset = lapply(subsets, function(sub){
      sub_lab = if(length(sub) == 0) "all"
                else paste(sprintf("%s=%s", names(sub), unlist(sub)),
                           collapse = ", ")

      de_sub = subsetDataset(de, conditions = if(length(sub)) as.list(sub) else NULL)
      X = as.matrix(de_sub$data)

      if(ncol(X) < 2 || nrow(X) < 3){
        warning(sprintf("Correlation '%s' / %s: insufficient data (n=%d, features=%d).",
                        name, sub_lab, nrow(X), ncol(X)))
        return(.empty_corr_df())
      }

      per_method = lapply(methods, function(m){
        .pairwise_cor_long(X, method = m, name = name, sub_lab = sub_lab,
                           include_self = include_self)
      })
      dplyr::bind_rows(per_method)
    })

    dplyr::bind_rows(per_subset)
  })

  # One family per named analysis x subset x method: a Pearson and a Spearman
  # run over the same pairs are two answers to one question, not twice the
  # multiple testing.
  out = .adjust_within(dplyr::bind_rows(rows),
                       c("correlation", "subset", "method"),
                       p_adjust, sig_threshold)
  out
}

#' Pairwise feature correlations in long format
#'
#' Returns the lower triangle (plus the diagonal when `include_self = TRUE`),
#' one row per feature pair.
#'
#' Pearson p-values are computed vectorised, from the usual
#' \eqn{t = r\sqrt{(n-2)/(1-r^2)}}, which is exactly what [stats::cor.test()]
#' reports.
#'
#' Spearman is NOT computed that way. Applying the same formula to a rank
#' correlation gives the asymptotic approximation that `cor.test()` uses only
#' when `exact = FALSE`; for the small, untied datasets typical of a targeted
#' study `cor.test()` defaults to an exact permutation p-value instead, and the
#' two can differ materially at the sample sizes that matter here. So each
#' Spearman pair is passed to [stats::cor.test()], which picks exact or
#' asymptotic by its own rules. That costs one call per pair - a few seconds
#' for a panel of a few hundred features - and is worth it for p-values that
#' match what a reader gets when they check one by hand.
#'
#' @param X Numeric matrix, samples x features.
#' @param method `"pearson"` or `"spearman"`.
#' @param name,sub_lab Labels identifying the analysis and subset.
#' @param include_self Keep the diagonal (self-correlations).
#' @return A long data frame: `correlation`, `subset`, `method`, `feature_a`,
#'   `feature_b`, `n`, `estimate`, `statistic`, `p_value`. `statistic` is the
#'   t-statistic for Pearson and Spearman's S for Spearman.
#' @keywords internal
#' @noRd
.pairwise_cor_long = function(X, method, name, sub_lab,
                              include_self = FALSE){

  r_mat = suppressWarnings(stats::cor(X, use = "pairwise.complete.obs",
                                      method = method))

  # Per-pair sample count (handles the NA pattern)
  ok    = !is.na(X)
  n_mat = crossprod(ok)

  feats = colnames(X)
  idx   = which(lower.tri(r_mat, diag = include_self), arr.ind = TRUE)

  if(identical(method, "pearson")){
    denom = 1 - r_mat^2
    denom[denom <= 0] = NA_real_
    t_mat = r_mat * sqrt((n_mat - 2) / denom)
    p_mat = 2 * stats::pt(-abs(t_mat), df = n_mat - 2)
    diag(p_mat) = NA_real_
    stat = t_mat[idx]
    pval = p_mat[idx]
  } else {
    res = vapply(seq_len(nrow(idx)), function(k){
      i = idx[k, 1]; j = idx[k, 2]
      if(i == j) return(c(NA_real_, NA_real_))
      keep = !is.na(X[, i]) & !is.na(X[, j])
      if(sum(keep) < 3) return(c(NA_real_, NA_real_))
      ct = tryCatch(suppressWarnings(
             stats::cor.test(X[keep, i], X[keep, j], method = method)),
             error = function(e) NULL)
      if(is.null(ct)) return(c(NA_real_, NA_real_))
      c(unname(ct$statistic), ct$p.value)
    }, numeric(2))
    stat = res[1, ]
    pval = res[2, ]
  }

  # With fewer than three shared samples a correlation is +/-1 or undefined
  # by construction, not an estimate: left out rather than drawn as perfect.
  thin = n_mat[idx] < 3
  est  = r_mat[idx]
  est[thin] = NA_real_; stat[thin] = NA_real_; pval[thin] = NA_real_

  data.frame(
    correlation = name,
    subset      = sub_lab,
    method      = method,
    feature_a   = feats[idx[, 1]],
    feature_b   = feats[idx[, 2]],
    n           = n_mat[idx],
    estimate    = est,
    statistic   = stat,
    p_value     = pval,
    stringsAsFactors = FALSE
  )
}

# ---- Linear models ---------------------------------------------------------

#' @keywords internal
#' @noRd
.run_linear_models = function(de, linear_models, p_adjust, sig_threshold){
  if(length(linear_models) == 0)
    return(.empty_lm_df())

  rows = lapply(linear_models, function(mod){
    name = mod$name %||% "unnamed_model"
    rhs  = mod$formula %||% stop(sprintf("Linear model '%s': `formula` required.", name))
    cond = mod$subset; if(!is.null(cond)) cond = as.list(cond)
    de_sub = subsetDataset(de, conditions = cond)

    feats = colnames(de_sub$data)
    smeta = de_sub$sample_meta

    inner = lapply(feats, function(f){
      df_in = cbind(.feature = de_sub$data[[f]], smeta)
      fm = stats::as.formula(paste(".feature", rhs))
      fit = tryCatch(stats::lm(fm, data = df_in), error = function(e) NULL)
      if(is.null(fit)) return(NULL)
      co = tryCatch(stats::coef(summary(fit)), error = function(e) NULL)
      if(is.null(co)) return(NULL)

      # ---- friendlier term labels (S9) ----
      # Default lm() names terms like "(Intercept)" and "TreatmentHDM_Curdlan".
      # Rename to make the reference level explicit:
      #   "(Intercept)"            -> "Intercept (reference: Treatment=Ctrl, Sex=Female)"
      #   "TreatmentHDM_Curdlan"   -> "Treatment: HDM_Curdlan (vs Ctrl)"
      raw_terms = rownames(co)
      term_labels = .pretty_lm_terms(raw_terms, fit, smeta)

      data.frame(
        model      = name,
        formula    = rhs,
        feature    = f,
        term       = term_labels,
        term_raw   = raw_terms,
        estimate   = co[, "Estimate"],
        std_error  = co[, "Std. Error"],
        statistic  = co[, "t value"],
        p_value    = co[, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(inner)
  })

  # One family per model term, across features: the question asked of a
  # treatment coefficient is not the question asked of the intercept.
  out = .adjust_within(dplyr::bind_rows(rows), c("model", "term"),
                       p_adjust, sig_threshold)
  out
}


# ---- Empty-frame templates (so xlsx writer always sees the schema) --------

#' @keywords internal
#' @noRd
.empty_stats_df = function() data.frame(
  comparison = character(), method = character(), feature = character(),
  n_groups = integer(), group1 = character(), group2 = character(),
  n_g1 = integer(), n_g2 = integer(),
  n_na_g1 = integer(), n_na_g2 = integer(),
  mean_g1 = numeric(), mean_g2 = numeric(), mean_diff = numeric(),
  FC = numeric(), log2FC = numeric(),
  statistic = numeric(), df = character(),
  p_value = numeric(), conf_low = numeric(), conf_high = numeric(),
  effect_size = numeric(), effect_type = character(),
  p_adj = numeric(), significant = logical(), adjust_family = character(),
  stringsAsFactors = FALSE
)

#' @keywords internal
#' @noRd
.empty_corr_df = function() data.frame(
  correlation = character(), subset = character(), method = character(),
  feature_a = character(), feature_b = character(), n = integer(),
  estimate = numeric(), statistic = numeric(),
  p_value = numeric(), p_adj = numeric(), significant = logical(),
  adjust_family = character(),
  stringsAsFactors = FALSE
)

#' @keywords internal
#' @noRd
.empty_lm_df = function() data.frame(
  model = character(), formula = character(), feature = character(),
  term = character(), term_raw = character(),
  estimate = numeric(), std_error = numeric(), statistic = numeric(),
  p_value = numeric(), p_adj = numeric(), significant = logical(),
  adjust_family = character(),
  stringsAsFactors = FALSE
)

#' Translate raw lm() coefficient names into something readable.
#'
#' Examples (with Treatment  in  {Ctrl, HDM_Curdlan}, Sex  in  {Female, Male}):
#'  "(Intercept)"           -> "Intercept (ref: Treatment=Ctrl; Sex=Female)"
#'  "TreatmentHDM_Curdlan"  -> "Treatment: HDM_Curdlan (vs Ctrl)"
#'  "SexMale:TreatmentHDM_Curdlan" -> kept as-is (interaction); doc'd elsewhere
#' @keywords internal
#' @noRd
.pretty_lm_terms = function(raw_terms, fit, smeta){
  # Pull factor reference levels from the model frame xlevels
  xlev = tryCatch(fit$xlevels, error = function(e) list())
  ref_levels = vapply(names(xlev), function(v) xlev[[v]][1], character(1))
  names(ref_levels) = names(xlev)

  vapply(raw_terms, function(t){
    if(t == "(Intercept)"){
      if(length(ref_levels) == 0) return("Intercept")
      ref_str = paste(sprintf("%s=%s", names(ref_levels), unname(ref_levels)),
                      collapse = "; ")
      return(sprintf("Intercept (ref: %s)", ref_str))
    }
    # Try to peel "<factor><level>" off the start of the term
    for(v in names(ref_levels)){
      if(startsWith(t, v)){
        lvl = sub(paste0("^", v), "", t)
        # Skip pure numeric coefficients (continuous predictors)
        if(nzchar(lvl) && lvl %in% xlev[[v]])
          return(sprintf("%s: %s (vs %s)", v, lvl, ref_levels[[v]]))
      }
    }
    t
  }, character(1), USE.NAMES = FALSE)
}

#' Extract the entry list from a stats section, honouring its section-level
#' `enabled:` toggle.
#'
#' Schema: each of `comparisons:`, `correlations:`, `linear_models:` is a
#' block with a section-level `enabled:` switch and an `entries:` list:
#'   comparisons:
#'     enabled: True
#'     entries:
#'       - name: ...
#' Returns `list()` when the section is absent or `enabled: False`; otherwise
#' the `entries:` list, minus any entry switched off with its own
#' `enabled: False` - so runStats, the report and the validators all skip it
#' alike. A bare top-level list errors loudly so an un-migrated config is
#' caught rather than silently producing no results.
#' @keywords internal
#' @noRd
.section_entries = function(section, section_name){
  if(is.null(section)) return(list())
  if(!is.null(section$entries) || !is.null(section$enabled)){
    if(isFALSE(section$enabled)) return(list())
    return(Filter(function(e) !(is.list(e) && isFALSE(e$enabled)),
                  section$entries %||% list()))
  }
  stop(sprintf(
    "[runStats] '%s:' must be a block with 'enabled:' and 'entries:' (got a bare list). Wrap the list under 'entries:' and add 'enabled: True'.",
    section_name))
}

#' Entries of the `ion_ratios:` block
#'
#' Same contract as [.section_entries()], except this block names its list
#' `ratios:` rather than `entries:`.
#'
#' @param section The `ion_ratios:` block, or `NULL`.
#' @return A list of ratio specs; empty when absent or disabled.
#' @keywords internal
#' @noRd
.ion_entries = function(section){
  if(is.null(section)) return(list())
  if(isFALSE(section$enabled)) return(list())
  section$ratios %||% list()
}

#' @keywords internal
#' @noRd
.row = function(comp_name, method, feature, n_groups,
                group1, group2, mean_g1, mean_g2, FC, log2FC,
                statistic, df, p_value,
                n_g1 = NA_integer_, n_g2 = NA_integer_,
                n_na_g1 = NA_integer_, n_na_g2 = NA_integer_,
                mean_diff = NA_real_,
                conf_low = NA_real_, conf_high = NA_real_,
                effect_size = NA_real_, effect_type = NA_character_){
  data.frame(
    comparison = comp_name, method = method, feature = feature,
    n_groups = as.integer(n_groups),
    group1 = as.character(group1), group2 = as.character(group2),
    n_g1 = as.integer(n_g1), n_g2 = as.integer(n_g2),
    n_na_g1 = as.integer(n_na_g1), n_na_g2 = as.integer(n_na_g2),
    mean_g1 = as.numeric(mean_g1), mean_g2 = as.numeric(mean_g2),
    mean_diff = as.numeric(mean_diff),
    FC = as.numeric(FC),
    log2FC = as.numeric(log2FC),
    statistic = as.numeric(statistic),
    df = as.character(df),
    p_value = as.numeric(p_value),
    conf_low = as.numeric(conf_low), conf_high = as.numeric(conf_high),
    effect_size = as.numeric(effect_size),
    effect_type = as.character(effect_type),
    p_adj = NA_real_, significant = NA, adjust_family = NA_character_,
    stringsAsFactors = FALSE
  )
}


# ---- Ion ratios ------------------------------------------------------------

#' Per-sample ion ratios, with a group test per comparison
#'
#' Ion ratios are enzyme-activity surrogates: the ratio of a product to its
#' precursor tracks flux through a pathway more directly than either compound
#' alone, because it cancels the shared variation in how much substrate was
#' there to begin with.
#'
#' Returns one row per sample per ratio. Rows carrying a `comparison` are that
#' comparison's samples, with the group label and the group test attached
#' (Wilcoxon for two levels, Kruskal-Wallis for more); rows with `comparison`
#' `NA` are every sample, unscoped, for descriptive views.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param ratios The `ion_ratios$ratios` entries from the YAML - each a list
#'   with `num`, `den` and optionally `label`.
#' @param comparisons The comparison entries, as unpacked by
#'   `.section_entries()`.
#' @param p_adjust Adjustment method passed to [stats::p.adjust()].
#' @return A long data frame, or a zero-row frame when no ratio is computable.
#' @keywords internal
#' @noRd
.run_ion_ratios = function(de, ratios, comparisons, p_adjust = "BH"){

  empty = data.frame(ratio = character(0), numerator = character(0),
                     denominator = character(0), comparison = character(0),
                     sample = character(0), group = character(0),
                     value = numeric(0), method = character(0),
                     statistic = numeric(0), p_value = numeric(0),
                     p_adj = numeric(0), stringsAsFactors = FALSE)
  if(length(ratios) == 0) return(empty)

  dm  = as.data.frame(de$data)
  out = list()

  for(cfg in ratios){
    num = as.character(cfg$num)
    den = as.character(cfg$den)
    lbl = as.character(cfg$label %||% paste0(num, " / ", den))

    absent = setdiff(c(num, den), colnames(dm))
    if(length(absent)){
      # Named, because a silent skip leaves an empty ion_ratios sheet with no
      # hint that a feature name - often a combine prefix - did not match.
      warning(sprintf(
        "[runStats] metabolite ratio '%s' skipped: %s not in the data.",
        lbl, paste(sprintf("'%s'", absent), collapse = " and ")),
        call. = FALSE)
      next
    }

    # Zero denominators give Inf; treat as missing rather than as a huge ratio.
    val = as.numeric(dm[[num]]) / as.numeric(dm[[den]])
    val[!is.finite(val)] = NA_real_
    names(val) = rownames(dm)

    out[[length(out) + 1]] = data.frame(
      ratio = lbl, numerator = num, denominator = den,
      comparison = NA_character_, sample = rownames(dm),
      group = NA_character_, value = unname(val),
      method = NA_character_, statistic = NA_real_,
      p_value = NA_real_, p_adj = NA_real_, stringsAsFactors = FALSE)

    for(comp in comparisons){
      if(!isTRUE(comp$enabled %||% TRUE)) next
      cond = if(is.null(comp$subset)) list() else as.list(comp$subset)
      cond[[comp$compare$factor]] = comp$compare$levels
      de_sub = tryCatch(subsetDataset(de, conditions = cond),
                        error = function(e) NULL)
      if(is.null(de_sub) || nrow(as.data.frame(de_sub$sample_meta)) == 0) next

      rn  = rownames(as.data.frame(de_sub$data))
      v   = val[rn]
      grp = factor(as.data.frame(de_sub$sample_meta)[[comp$compare$factor]],
                   levels = comp$compare$levels)

      ngrp = length(comp$compare$levels)
      tst = tryCatch({
        if(ngrp == 2)
          stats::wilcox.test(v[grp == comp$compare$levels[1]],
                             v[grp == comp$compare$levels[2]], exact = FALSE)
        else
          stats::kruskal.test(v, grp)
      }, error = function(e) NULL)

      out[[length(out) + 1]] = data.frame(
        ratio = lbl, numerator = num, denominator = den,
        comparison = comp$name %||% "unnamed_comparison",
        sample = rn, group = as.character(grp), value = unname(v),
        method = if(ngrp == 2) "wilcoxon" else "kruskal",
        statistic = if(is.null(tst)) NA_real_ else unname(tst$statistic),
        p_value = if(is.null(tst)) NA_real_ else tst$p.value,
        p_adj = NA_real_, stringsAsFactors = FALSE)
    }
  }

  if(length(out) == 0) return(empty)
  res = dplyr::bind_rows(out)

  # One adjustment across the distinct ratio x comparison tests, not across
  # the repeated per-sample rows that carry them.
  keys = unique(res[!is.na(res$p_value), c("ratio", "comparison", "p_value")])
  if(nrow(keys)){
    keys$adj = stats::p.adjust(keys$p_value, method = p_adjust)
    idx = match(paste(res$ratio, res$comparison),
                paste(keys$ratio, keys$comparison))
    res$p_adj = keys$adj[idx]
  }
  res
}
