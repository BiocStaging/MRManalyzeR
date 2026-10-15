#' Filter features by detection frequency
#'
#' Drops features that are missing from too much of the data. A feature's
#' detection fraction is the proportion of samples in which it has a non-`NA`
#' value; `method` decides which samples that proportion is computed over.
#'
#' `"across"` pools every study sample into one fraction. It is the strictest
#' option and the easiest to misuse: a feature present throughout one treatment
#' group and absent from the others scores roughly `1 / n_groups`, so a
#' threshold of 0.5 deletes it *because* it is treatment-specific. Use it only
#' when the question is "is this measured everywhere".
#'
#' `"within"` (the default) computes the fraction separately within each level
#' of `group_col` and keeps a feature that clears `min_frac` in **any one**
#' level, so the treatment-specific feature above survives on its own group.
#' Small groups quantise the threshold: with `n = 3` a fraction can only be 0,
#' 1/3, 2/3 or 1, so every `min_frac` in (1/3, 2/3] behaves identically.
#'
#' `"QC"` computes the fraction over QC injections only, which answers "is this
#' measured reproducibly" rather than "is this present in the biology". QCs are
#' pooled, so a feature confined to one treatment is diluted in the pool and can
#' fall below detection there - this mode also removes treatment-specific
#' features. With no QCs present the function warns and falls back.
#'
#' QC and blank injections never count towards `"within"` or `"across"`. They
#' remain in the returned dataset, since the blank filter and [correctBatch()]
#' still need them, but they are excluded from the denominator so a compound
#' seen only in the QCs cannot keep itself alive.
#'
#' The fraction actually used is recorded per feature in `variable_meta` as
#' `detected_frac`, so the decision is inspectable rather than only visible as
#' absent columns.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param min_frac Numeric in `[0, 1]`. A feature is kept when its detection
#'   fraction is greater than or equal to this.
#' @param method `"within"` (default), `"across"` or `"QC"`; see Details.
#' @param group_col `sample_meta` column whose levels are tested separately
#'   when `method = "within"`. Required for that method.
#' @param sample_col `sample_meta` column identifying the sample type.
#' @param sample_labels Values of `sample_col` marking study samples. Only
#'   these count towards `"within"` and `"across"`.
#' @param qc_label Value of `sample_col` marking QC injections, used by
#'   `method = "QC"`.
#' @return `de` with failing features removed from `data` and `variable_meta`,
#'   and a `detected_frac` column added to `variable_meta`.
#' @family peak-matrix processing
#' @seealso [filterSamples()], which drops injections rather than features and
#'   is meant to run after this one.
#' @examples
#' de <- struct::DatasetExperiment(
#'   data = data.frame(
#'     good   = c(10, 12, 11, 9, 10, 11),
#'     npOnly = c(NA, NA, NA, 20, 22, 21),
#'     sparse = c(5, NA, NA, NA, NA, 6),
#'     row.names = paste0("S", 1:6)),
#'   sample_meta = data.frame(
#'     Sample_type = rep("Sample", 6),
#'     Treatment   = rep(c("Control", "NP"), each = 3),
#'     row.names   = paste0("S", 1:6)),
#'   variable_meta = data.frame(
#'     Compound  = c("good", "npOnly", "sparse"),
#'     row.names = c("good", "npOnly", "sparse")))
#'
#' # 'within' keeps npOnly on the strength of its own group; 'across' would
#' # delete it for being absent from two thirds of the samples.
#' filterFeatures(de, min_frac = 0.5, method = "within", group_col = "Treatment")
#' @export
filterFeatures <- function(de, min_frac = 0.5,
                          method = c("within", "across", "QC"),
                          group_col = NULL,
                          sample_col = "Sample_type",
                          sample_labels = "Sample",
                          qc_label = "QC"){
  .check_de(de, "filterFeatures")

  method <- match.arg(method)
  if(!is.numeric(min_frac) || length(min_frac) != 1L ||
     is.na(min_frac) || min_frac < 0 || min_frac > 1)
    stop("filterFeatures: min_frac must be a single number in [0, 1].")

  X    <- as.data.frame(de$data)
  smeta <- as.data.frame(de$sample_meta)
  if(!ncol(X)) return(de)

  # Which rows are study samples. A missing column is a config error worth
  # saying out loud, but not worth stopping a run over - fall back to every
  # row and report what that means.
  study <- rep(TRUE, nrow(X))
  if(method %in% c("within", "across") && !is.null(sample_col)){
    if(!sample_col %in% colnames(smeta)){
      fmt <- paste0("filterFeatures: '%s' not in sample_meta - every ",
                    "injection counts towards detection, QCs and ",
                    "blanks included.")
      warning(sprintf(fmt, sample_col))
    } else {
      study <- as.character(smeta[[sample_col]]) %in% sample_labels
      if(!any(study))
        stop(sprintf("filterFeatures: no rows where %s is one of: %s.",
                     sample_col, paste(sample_labels, collapse = ", ")))
    }
  }

  # QC mode needs QCs; without them it cannot answer its own question, and a
  # copied config block should not take the whole run down.
  if(identical(method, "QC")){
    has_qc <- sample_col %in% colnames(smeta) &&
             any(as.character(smeta[[sample_col]]) %in% qc_label)
    if(!has_qc){
      method <- if(is.null(group_col)) "across" else "within"
      fmt <- paste0("filterFeatures: no QC injections (%s == '%s') - ",
                    "falling back to method = '%s'.")
      warning(sprintf(fmt, sample_col, qc_label, method))
      if(sample_col %in% colnames(smeta))
        study <- as.character(smeta[[sample_col]]) %in% sample_labels
    }
  }

  frac_of <- function(rows) colMeans(!is.na(X[rows, , drop = FALSE]))

  if(identical(method, "QC")){
    best <- frac_of(as.character(smeta[[sample_col]]) %in% qc_label)

  } else if(identical(method, "across")){
    best <- frac_of(study)

  } else {
    if(is.null(group_col) || !group_col %in% colnames(smeta)){
      fmt <- paste0("filterFeatures: method 'within' needs group_col to ",
                    "name a sample_meta column (got %s).")
      got <- if(is.null(group_col)) "NULL" else sprintf("'%s'", group_col)
      stop(sprintf(fmt, got))
    }
    Xs <- X[study, , drop = FALSE]
    g <- as.character(smeta[[group_col]])[study]
    g[is.na(g) | !nzchar(g)] <- "_unassigned_"
    # Best group wins: a feature that is real in one treatment survives the
    # groups it is genuinely absent from.
    per <- vapply(split(seq_len(nrow(Xs)), g),
                 function(i) colMeans(!is.na(Xs[i, , drop = FALSE])),
                 numeric(ncol(Xs)))
    best <- if(is.matrix(per)) apply(per, 1L, max) else per
    names(best) <- colnames(Xs)
  }

  keep <- best >= min_frac
  drop <- colnames(X)[!keep]

  vm <- as.data.frame(de$variable_meta)
  .warn_overwrite(vm, "detected_frac", "filterFeatures")
  vm$detected_frac <- round(unname(best[match(rownames(vm), names(best))]), 4)
  de$variable_meta <- vm

  if(length(drop)){
    message(sprintf(
      "[filterFeatures] %s: dropped %d of %d features below min_frac = %g%s",
      method, length(drop), ncol(X), min_frac,
      if(identical(method, "within")) sprintf(" (best of %s)", group_col) else ""))
    message("   ", paste(utils::head(drop, 10), collapse = ", "),
            if(length(drop) > 10)
              sprintf(" ... and %d more", length(drop) - 10) else "")
  } else {
    message(sprintf("[filterFeatures] %s: all %d features retained.",
                    method, ncol(X)))
  }

  if(!any(keep))
    stop("filterFeatures: min_frac removed every feature.")

  subsetDataset(de, features = colnames(X)[keep],
                drop_empty_features = FALSE)
}
