#' Append/replace stats sheets in the results xlsx
#'
#' Adds (or overwrites) the following tabs in the existing results xlsx:
#' \itemize{
#'   \item `key`           - what every column heading in this workbook means.
#'   \item `stats`         - long table of comparisons + post-hoc.
#'   \item `correlations`  - long table of feature x feature correlations.
#'   \item `linear_models` - long table of `lm()` coefficients.
#'   \item `ion_ratios`    - per-sample ion ratios with their group test.
#'   \item `summary`       - n / mean / SD / SE per group per compound.
#'   \item `cor_<name>_<subset>_<method>` - one wide correlation matrix tab
#'         per (correlation x subset x method). Sheet names are sanitised
#'         and truncated to 31 characters (Excel limit).
#' }
#' If the workbook does not exist (e.g. results dir was cleaned), a new one
#' is written containing only the stats tabs.
#'
#' @param out_xlsx Path to the results xlsx (from `runMRManalyzeR()`).
#' @param stats_tables A list with elements `stats`, `correlations`,
#'   `linear_models` and `ion_ratios` - as returned by [`runStats()`].
#' @param group_summary Optional data frame from [`summariseGroups()`],
#'   written as the `summary` tab. `NULL` skips it.
#' @return Invisibly, the path written.
#' @keywords internal
append_stats_xlsx <- function(out_xlsx, stats_tables,
                             group_summary = NULL){

  if(file.exists(out_xlsx)){
    wb <- openxlsx::loadWorkbook(out_xlsx)
  } else {
    wb <- openxlsx::createWorkbook()
  }

  # 1. Long-format tabs
  # Zero-row tables are skipped rather than written as a headers-only sheet:
  # an empty tab says nothing except that the block was disabled, and it is
  # easy to mistake for a run that failed.
  for(sheet_name in c("stats", "correlations", "linear_models",
                      "ion_ratios")){
    df <- stats_tables[[sheet_name]]
    if(is.null(df) || nrow(df) == 0) next
    .write_or_replace(wb, sheet_name, df)
  }

  # 1b. Per-comparison stats tabs
  if(!is.null(stats_tables$stats) && nrow(stats_tables$stats) > 0){
    for(cn in unique(stats_tables$stats$comparison)){
      sub <- stats_tables$stats[stats_tables$stats$comparison == cn, , drop = FALSE]
      .write_or_replace(wb, .safe_sheet_name(sprintf("stats_%s", cn)), sub)
    }
  }

  # 1c. Per-model linear model tabs
  if(!is.null(stats_tables$linear_models) && nrow(stats_tables$linear_models) > 0){
    for(mn in unique(stats_tables$linear_models$model)){
      sub <- stats_tables$linear_models[stats_tables$linear_models$model == mn, , drop = FALSE]
      .write_or_replace(wb, .safe_sheet_name(sprintf("lm_%s", mn)), sub)
    }
  }

  # 1d. Group summary - the numbers behind the boxplots / barplots.
  if(!is.null(group_summary) && nrow(group_summary) > 0)
    .write_or_replace(wb, "summary", group_summary)

  # 2. Wide-matrix tabs for correlations (one per subset x method x correlation)
  if(!is.null(stats_tables$correlations) && nrow(stats_tables$correlations) > 0){
    .write_corr_matrix_tabs(wb, stats_tables$correlations)
  }

  # 3. Column key, written last so it can describe only the sheets that were
  # actually produced, then moved to the front - a reader opening the
  # workbook should meet the glossary before the abbreviations.
  .write_key_sheet(wb)

  openxlsx::saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
  invisible(out_xlsx)
}


# ---- Internal helpers ------------------------------------------------------

#' Write the column key and move it to the front of the workbook
#'
#' Filtered to the sheets that were actually written, so a workbook produced
#' with `correlations: enabled: False` does not carry a page of definitions
#' for columns nobody will see.
#'
#' @keywords internal
#' @noRd
.write_key_sheet <- function(wb){
  present <- openxlsx::sheets(wb)
  # Every stats block disabled: a workbook holding nothing but a glossary of
  # columns it does not contain would be worse than no workbook.
  if(length(present) == 0) return(invisible(NULL))
  key <- .stats_key_df()

  # A trailing "*" marks a sheet family - `stats_<comparison>`, `lm_<model>`
  # and `cor_<...>` are generated one tab per entry, so the key describes the
  # family once rather than repeating itself per comparison.
  keep <- vapply(key$Sheet, function(s){
    if(s == "(all sheets)")             return(TRUE)
    if(endsWith(s, "*"))
      return(any(startsWith(present, sub("\\*$", "", s))))
    s %in% present
  }, logical(1), USE.NAMES = FALSE)

  key <- key[keep, , drop = FALSE]
  if(nrow(key) == 0) return(invisible(NULL))

  .write_or_replace(wb, "key", key)
  openxlsx::setColWidths(wb, "key", cols = seq_len(3),
                         widths = c(20, 22, 105))
  openxlsx::addStyle(wb, "key",
                     openxlsx::createStyle(wrapText = TRUE, valign = "top"),
                     rows = seq_len(nrow(key)) + 1L, cols = 3,
                     gridExpand = TRUE)
  openxlsx::freezePane(wb, "key", firstRow = TRUE)

  # worksheetOrder indexes the creation order that sheets() reports, and is
  # applied when the workbook is written - sheets() itself does not change.
  idx <- match("key", openxlsx::sheets(wb))
  openxlsx::`worksheetOrder<-`(
    wb, value = c(idx, setdiff(seq_along(openxlsx::sheets(wb)), idx)))
  invisible(NULL)
}

#' Definitions of every column heading the stats workbook can contain
#'
#' The workbook is the artefact that leaves the laboratory - it is read months
#' later, by people who did not run the pipeline, in Excel rather than in R.
#' `n_na_g1`, `epsilon_squared` and `adjust_family` are not self-explanatory,
#' and a column whose meaning has to be guessed is a column that will be
#' misread.
#'
#' @return A data frame of `Sheet` / `Column` / `Meaning`.
#' @keywords internal
#' @noRd
.stats_key_df <- function(){

  k <- function(sheet, ...){
    rows <- list(...)
    data.frame(
      Sheet   = sheet,
      Column  = vapply(rows, `[`, character(1), 1),
      Meaning = vapply(rows, `[`, character(1), 2),
      stringsAsFactors = FALSE)
  }

  rbind(
    k("(all sheets)",
      c("p_value", "Unadjusted p-value from the test named in `method`. Two exceptions, both flagged in `method`: `tukey` p-values are already adjusted across the pairs of that feature by Tukey HSD, and `pairwise_wilcoxon` p-values are already Holm-adjusted across those same pairs."),
      c("p_adj", "`p_value` after multiple-testing correction (Benjamini-Hochberg / FDR) applied within `adjust_family`."),
      c("significant", "TRUE when `p_adj` is below the `sig_threshold` set in the YAML config (0.05 unless changed)."),
      c("adjust_family", "The set of tests `p_adj` was computed over. Correction is applied within a family, never across the whole workbook: comparisons are corrected per comparison x method, correlations per correlation x subset x method, and linear models per model x term."),
      c("statistic", "The test statistic: t (Welch / Student), W (Wilcoxon), F (ANOVA), H (Kruskal-Wallis), the mean difference (Tukey), S (Spearman)."),
      c("feature", "Compound name, matching the `matrix` and `feature_metadata` sheets of the main results workbook."),
      c("(blank cell)", "Not applicable to that row rather than a failed calculation - e.g. `group1` / `group2` are blank on a global ANOVA or Kruskal-Wallis row, which tests all groups at once.")),

    k("stats",
      c("comparison", "Name of the comparison, from the `comparisons:` block of the YAML config."),
      c("method", "Test used. Two groups: `welch` (unequal-variance t-test), `student` (equal-variance t-test), `wilcoxon`; a `_paired` suffix means a paired design. Three or more groups: `anova` or `kruskal` for the global test, then `tukey` or `pairwise_wilcoxon` for the post-hoc pairs."),
      c("n_groups", "Number of groups this row's test covers - 2 for a pair, k for a global test across k groups."),
      c("group1", "First level of the comparison, the reference."),
      c("group2", "Second level. `FC`, `log2FC`, `mean_diff` and the confidence interval are all signed as group2 relative to group1."),
      c("n_g1", "Number of non-missing values in group1."),
      c("n_g2", "Number of non-missing values in group2."),
      c("n_na_g1", "Number of missing values in group1. Read alongside `n_g1`: a large fold-change resting on two surviving injections is not the same finding as one resting on ten."),
      c("n_na_g2", "Number of missing values in group2."),
      c("mean_g1", "Mean of group1, in the units of the matrix (peak area, ratio or concentration)."),
      c("mean_g2", "Mean of group2, same units."),
      c("mean_diff", "mean_g2 - mean_g1. For a Wilcoxon row this is still the difference in means, while the confidence interval is the Hodges-Lehmann location shift."),
      c("FC", "Fold change, mean_g2 / mean_g1."),
      c("log2FC", "log2(FC). Zero means no change; +1 is a doubling in group2, -1 a halving."),
      c("df", "Degrees of freedom. Text rather than a number because Welch's df is fractional and ANOVA has two (between,within)."),
      c("conf_low", "Lower bound of the 95% confidence interval on the difference, signed like `mean_diff`."),
      c("conf_high", "Upper bound of the same interval. An interval spanning zero is the same statement as p > 0.05."),
      c("effect_size", "Size of the difference, independent of n - the number to compare across compounds. See `effect_type` for which measure."),
      c("effect_type", "Which effect size: `hedges_g` (standardised mean difference, small-sample corrected; ~0.2 small, 0.5 medium, 0.8 large), `cohens_dz` (the paired equivalent, not comparable with g), `rank_biserial` / `rank_biserial_paired` (Wilcoxon, -1 to +1), `eta_squared` (ANOVA, share of total variance explained), `epsilon_squared` (its Kruskal-Wallis rank analogue).")),

    k("stats_*",
      c("(sheet)", "The `stats` rows for one comparison, split out for convenience. Same columns, same numbers - filtering the `stats` sheet gives the identical table.")),

    k("correlations",
      c("correlation", "Name of the correlation entry, from the `correlations:` block of the YAML config."),
      c("subset", "Samples the correlation was computed on, e.g. one treatment group. `all` means no subsetting."),
      c("method", "`pearson` (linear association, sensitive to outliers) or `spearman` (monotonic association on ranks)."),
      c("feature_a", "First compound of the pair."),
      c("feature_b", "Second compound of the pair. Order carries no meaning - correlation is symmetric."),
      c("n", "Number of samples with both compounds measured. Pairs are dropped one at a time, so `n` varies between rows."),
      c("estimate", "The correlation coefficient: r for Pearson, rho for Spearman. Runs -1 to +1.")),

    k("cor_*",
      c("(sheet)", "The same correlations as a square compound x compound matrix, one tab per correlation x subset x method, holding `estimate`. The diagonal is blank unless `include_self: True` was set.")),

    k("linear_models",
      c("model", "Name of the model, from the `linear_models:` block of the YAML config."),
      c("formula", "Right-hand side fitted for every compound, e.g. `~ Treatment * Sex` (main effects plus their interaction)."),
      c("term", "Readable name of the coefficient, including the reference level it is measured against."),
      c("term_raw", "The coefficient name exactly as `lm()` produced it, for anyone reproducing the fit in R."),
      c("estimate", "The coefficient: the change in the compound's value per unit of that term, holding the others fixed. For a factor, the difference from the reference level."),
      c("std_error", "Standard error of `estimate` - roughly, `estimate` +/- 2 x `std_error` is a 95% interval.")),

    k("lm_*",
      c("(sheet)", "The `linear_models` rows for one model, split out for convenience. Same columns, same numbers.")),

    k("ion_ratios",
      c("ratio", "Label for the ratio, from the `ion_ratios:` block of the YAML config."),
      c("numerator", "Compound on top of the ratio."),
      c("denominator", "Compound underneath. A zero denominator is recorded as missing rather than as an infinite ratio."),
      c("comparison", "Comparison the group test was run within, from the `comparisons:` block of the YAML config."),
      c("sample", "Sample the ratio was computed for - one row per sample, so the spread can be inspected rather than only the group test."),
      c("group", "Group the sample belongs to within `comparison`."),
      c("value", "numerator / denominator for that sample. Unitless, so it is comparable across batches and acquisitions in a way that either compound alone is not."),
      c("method", "Test used on the group values, chosen the same way as in `stats`.")),

    k("summary",
      c("Compound", "Compound name."),
      c("group", "Group the summary is over."),
      c("n", "Number of non-missing values in that group."),
      c("mean", "Arithmetic mean, in the units of the matrix."),
      c("sd", "Standard deviation - how spread out the samples are."),
      c("se", "Standard error of the mean, sd / sqrt(n) - how well the mean itself is pinned down. These are the exact numbers the report's bar charts and error bars are drawn from."))
  )
}

#' @keywords internal
#' @noRd
.write_or_replace <- function(wb, sheet_name, df){
  if(sheet_name %in% openxlsx::sheets(wb))
    openxlsx::removeWorksheet(wb, sheet_name)
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df,
                      colNames = TRUE, rowNames = FALSE, keepNA = FALSE)
}

#' Write one wide correlation matrix per (correlation x subset x method)
#' @keywords internal
#' @noRd
.write_corr_matrix_tabs <- function(wb, corr_long){
  combos <- unique(corr_long[, c("correlation", "subset", "method")])

  for(i in seq_len(nrow(combos))){
    cn <- combos$correlation[i]
    sb <- combos$subset[i]
    mt <- combos$method[i]

    sub <- corr_long[corr_long$correlation == cn &
                    corr_long$subset      == sb &
                    corr_long$method      == mt, , drop = FALSE]
    if(nrow(sub) == 0) next

    # Preserve the import order produced by .pairwise_cor_long() (column-major
    # walk over the lower triangle => first-seen order == colnames(X) order).
    feats <- unique(c(sub$feature_a, sub$feature_b))
    M <- matrix(NA_real_, nrow = length(feats), ncol = length(feats),
               dimnames = list(feats, feats))
    for(k in seq_len(nrow(sub))){
      a <- sub$feature_a[k]; b <- sub$feature_b[k]; r <- sub$estimate[k]
      M[a, b] <- r; M[b, a] <- r
    }
    # Self-self pairs are dropped from `correlations` (long table); leave the
    # diagonal as NA in the wide matrix so the empty cells make that obvious.
    # If the long table did include self-pairs (when include_self: True), the
    # loop above already filled the diagonal.

    sheet <- .safe_sheet_name(sprintf("cor_%s_%s_%s", cn, sb, mt))
    .write_or_replace(wb, sheet, as.data.frame(M, check.names = FALSE))
    # Re-write so row names appear (writeData with rowNames=TRUE)
    openxlsx::removeWorksheet(wb, sheet)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, as.data.frame(M, check.names = FALSE),
                        colNames = TRUE, rowNames = TRUE, keepNA = FALSE)
  }
}

#' Sanitise and truncate to Excel's 31-char sheet name limit.
#' @keywords internal
#' @noRd
.safe_sheet_name <- function(x){
  x <- gsub("[\\\\/?*\\[\\]:]", "_", x, perl = TRUE)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_|_$", "", x)
  if(nchar(x) > 31) x <- substr(x, 1, 31)
  if(!nzchar(x))    x <- "sheet"
  x
}
