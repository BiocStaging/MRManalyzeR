# runStats() is the scientific core: the test named in the YAML is the test
# run, the numbers are the ones stats:: would give for that test, and the
# multiple-testing correction is applied within the family each row names.
# Every expected value below is computed from stats:: on the same vectors, so
# a change in the engine has to be a change in the arithmetic to pass.

.F1  <- c(10, 11, 12, 13,  20, 21, 22, 23,  30, 31, 32, 33)  # A < B < C
.F2  <- c(5, 7, 6, 8,      9, 7, 10, 8,      6, 8, 7, 9)     # modest shift
.F3  <- c(2, NA, 3, 4,     9, 8, NA, 10,     3, 4, 5, NA)    # with gaps
.F4  <- c(1, 2, 3, 40,     8, 9, 10, 11,     5, 6, 7, 8)     # unequal spread
.GRP <- rep(c("A", "B", "C"), each = 4)
.SEX <- rep(c("F", "M"), 6)

.de_stats <- function() {
  ids <- paste0("S", 1:12)
  struct::DatasetExperiment(
    name = "stats_fixture",
    data = data.frame(F1 = .F1, F2 = .F2, F3 = .F3, F4 = .F4,
                      row.names = ids),
    sample_meta = data.frame(
      Sample_ID = ids,
      Treatment = .GRP,
      Sex       = .SEX,
      Subject   = rep(paste0("P", 1:4), 3),   # the same four, three times
      row.names = ids),
    variable_meta = data.frame(Compound = c("F1", "F2", "F3", "F4"),
                               row.names = c("F1", "F2", "F3", "F4")))
}

# One comparison, its arguments named as the YAML names them.
.cmp <- function(...) {
  suppressMessages(runStats(.de_stats(), list(
    comparisons = list(enabled = TRUE, entries = list(list(...))))))$stats
}

.AB  <- list(factor = "Treatment", levels = c("A", "B"))
.AC  <- list(factor = "Treatment", levels = c("A", "C"))
.ABC <- list(factor = "Treatment", levels = c("A", "B", "C"))

test_that("welch reports t.test's own numbers, signed towards group 2", {
  st <- .cmp(name = "AvB", method = "welch", compare = .AB)
  expect_equal(sort(unique(st$feature)), c("F1", "F2", "F3", "F4"))
  expect_true(all(st$method == "welch"))
  expect_true(all(st$n_groups == 2))
  expect_true(all(st$group1 == "A" & st$group2 == "B"))

  a <- .F1[.GRP == "A"]; b <- .F1[.GRP == "B"]
  tt <- stats::t.test(a, b)
  r  <- st[st$feature == "F1", ]

  expect_equal(r$p_value, tt$p.value)
  expect_equal(r$statistic, unname(tt$statistic))
  expect_equal(r$mean_g1, mean(a))
  expect_equal(r$mean_g2, mean(b))
  expect_equal(r$mean_diff, mean(b) - mean(a))
  expect_equal(r$FC, mean(b) / mean(a))
  expect_equal(r$log2FC, log2(mean(b) / mean(a)))
  expect_equal(r$n_g1, 4L)
  expect_equal(r$n_g2, 4L)

  # t.test() forms its interval as group1 - group2. The table reports it the
  # other way round, so the interval and the fold change point the same way.
  expect_equal(r$conf_low, -tt$conf.int[2])
  expect_equal(r$conf_high, -tt$conf.int[1])
  expect_gt(r$conf_low, 0)
  expect_gt(r$effect_size, 0)
  expect_equal(r$effect_type, "hedges_g")

  # Missing values are counted, not silently dropped.
  f3 <- st[st$feature == "F3", ]
  expect_equal(f3$n_g1, 3L)
  expect_equal(f3$n_na_g1, 1L)
  expect_equal(f3$n_g2, 3L)
  expect_equal(f3$n_na_g2, 1L)
  expect_equal(f3$mean_g1, mean(.F3[.GRP == "A"], na.rm = TRUE))
})

test_that("student is the equal-variance test, and says so", {
  st <- .cmp(name = "AvB", method = "student", compare = .AB)
  a  <- .F1[.GRP == "A"]; b <- .F1[.GRP == "B"]
  r  <- st[st$feature == "F1", ]

  expect_true(all(st$method == "student"))
  expect_equal(r$p_value, stats::t.test(a, b, var.equal = TRUE)$p.value)
  expect_equal(r$df, "6")               # n1 + n2 - 2, not Welch's fraction

  # On F4, where the group variances differ by two orders of magnitude, the
  # two tests part company - which is the whole reason for naming one.
  welch <- .cmp(name = "AvB", method = "welch", compare = .AB)
  f4 <- function(d) d[d$feature == "F4", ]
  expect_equal(f4(st)$df, "6")
  expect_equal(f4(welch)$df, as.character(round(
    stats::t.test(.F4[.GRP == "A"], .F4[.GRP == "B"])$parameter, 2)))
  expect_false(identical(f4(welch)$df, f4(st)$df))
  expect_false(isTRUE(all.equal(f4(welch)$p_value, f4(st)$p_value)))
})

test_that("wilcoxon reports the rank test and its own effect size", {
  st <- .cmp(name = "AvB", method = "wilcoxon", compare = .AB)
  a  <- .F1[.GRP == "A"]; b <- .F1[.GRP == "B"]
  wt <- suppressWarnings(
    stats::wilcox.test(a, b, exact = FALSE, conf.int = TRUE))
  r  <- st[st$feature == "F1", ]

  expect_true(all(st$method == "wilcoxon"))
  expect_equal(r$p_value, wt$p.value)
  expect_equal(r$statistic, unname(wt$statistic))
  expect_true(is.na(r$df))
  expect_equal(r$effect_type, "rank_biserial")
  expect_equal(r$effect_size, 1 - 2 * unname(wt$statistic) / 16)

  # The Hodges-Lehmann interval, negated like the parametric one.
  expect_equal(r$conf_low, -wt$conf.int[2])
  expect_equal(r$conf_high, -wt$conf.int[1])
})

test_that("a comparison with no method: gets the default, and is told so", {
  de <- .de_stats()

  expect_message(
    two <- runStats(de, list(comparisons = list(enabled = TRUE, entries = list(
      list(name = "AvB", compare = .AB))))),
    "defaulting to 'welch'")
  expect_true(all(two$stats$method == "welch"))

  expect_message(
    many <- runStats(de, list(comparisons = list(enabled = TRUE, entries = list(
      list(name = "ABC", compare = .ABC))))),
    "defaulting to 'anova'")
  expect_true("anova" %in% many$stats$method)
})

test_that("a method that does not fit the design is refused", {
  expect_error(.cmp(name = "AvB", method = "kruskal", compare = .AB),
               "not valid for 2 groups")
  expect_error(.cmp(name = "ABC", method = "welch", compare = .ABC),
               "not valid for 3 groups")
  expect_error(.cmp(name = "AvB", method = "ttest", compare = .AB),
               "Use one of: welch, student, wilcoxon")

  # The two spellings configurations used before the methods were named.
  expect_true(all(.cmp(name = "AvB", method = "t_test",
                       compare = .AB)$method == "welch"))
  expect_true(all(.cmp(name = "AvB", method = "t-test",
                       compare = .AB)$method == "welch"))
})

test_that("a missing factor or levels is refused by name", {
  expect_error(.cmp(name = "x", compare = list(levels = c("A", "B"))),
               "compare.factor and compare.levels required")
  expect_error(.cmp(name = "x", compare = list(factor = "Nope",
                                               levels = c("A", "B"))),
               "factor 'Nope' not in sample_meta")
})

test_that("three groups gives the global test and its Tukey post-hoc", {
  st <- .cmp(name = "ABC", method = "anova", compare = .ABC)
  d  <- data.frame(value = .F1, grp = factor(.GRP))
  av <- stats::aov(value ~ grp, data = d)
  sm <- summary(av)[[1]]

  glob <- st[st$method == "anova", ]
  expect_equal(nrow(glob), 4)                  # one per feature
  expect_true(all(glob$n_groups == 3))
  expect_true(all(glob$effect_type == "eta_squared"))
  expect_equal(glob$p_value[glob$feature == "F1"], sm[["Pr(>F)"]][1])
  expect_equal(glob$statistic[glob$feature == "F1"], sm[["F value"]][1])
  expect_equal(glob$df[glob$feature == "F1"], "2,9")

  # The global row describes no pair, so it carries no group means.
  expect_true(all(is.na(glob$mean_g1)))

  ph <- st[st$method == "tukey" & st$feature == "F1", ]
  expect_equal(nrow(ph), 3)
  expect_true(all(ph$n_groups == 2))

  # TukeyHSD labels a pair "B-A". Read the wrong way round, every post-hoc
  # fold change would carry the wrong sign.
  tk  <- stats::TukeyHSD(av)$grp
  row <- ph[ph$group1 == "A" & ph$group2 == "B", ]
  expect_equal(nrow(row), 1)
  expect_equal(row$mean_g1, mean(.F1[.GRP == "A"]))
  expect_equal(row$mean_g2, mean(.F1[.GRP == "B"]))
  expect_equal(row$mean_diff, unname(tk["B-A", "diff"]))
  expect_gt(row$mean_diff, 0)
  expect_equal(row$p_value, unname(tk["B-A", "p adj"]))
  expect_equal(row$conf_low, unname(tk["B-A", "lwr"]))
  expect_equal(row$conf_high, unname(tk["B-A", "upr"]))
})

test_that("kruskal gets a rank post-hoc, Holm-adjusted within the feature", {
  st <- .cmp(name = "ABC", method = "kruskal", compare = .ABC)
  d  <- data.frame(value = .F1, grp = factor(.GRP))
  kw <- stats::kruskal.test(value ~ grp, data = d)

  g <- st[st$method == "kruskal" & st$feature == "F1", ]
  expect_equal(g$p_value, kw$p.value)
  expect_equal(g$statistic, unname(kw$statistic))
  expect_equal(g$effect_type, "epsilon_squared")
  expect_equal(g$effect_size, unname(kw$statistic) / (12 - 1))

  # The post-hoc matches the global test - pairwise Wilcoxon, not Tukey, which
  # would be read off an ANOVA fit the global test never used.
  ph <- st[st$method == "pairwise_wilcoxon" & st$feature == "F1", ]
  expect_equal(nrow(ph), 3)
  expect_false("tukey" %in% st$method)

  pairs <- list(c("A", "B"), c("A", "C"), c("B", "C"))
  raw <- vapply(pairs, function(p) suppressWarnings(stats::wilcox.test(
    .F1[.GRP == p[1]], .F1[.GRP == p[2]], exact = FALSE,
    conf.int = TRUE)$p.value), numeric(1))
  expect_equal(sort(ph$p_value), sort(stats::p.adjust(raw, method = "holm")))
})

test_that("a paired comparison pairs by subject and names its own effect", {
  st <- .cmp(name = "AvB", method = "welch", paired = TRUE,
             subject = "Subject", compare = .AB)
  a  <- .F2[.GRP == "A"]; b <- .F2[.GRP == "B"]
  r  <- st[st$feature == "F2", ]

  expect_equal(r$method, "welch_paired")
  expect_equal(r$p_value, stats::t.test(a, b, paired = TRUE)$p.value)
  expect_equal(r$effect_type, "cohens_dz")
  expect_equal(r$effect_size, mean(b - a) / stats::sd(b - a))

  # A paired effect size is not comparable to the unpaired one, so the two
  # never come back under the same label.
  unp <- .cmp(name = "AvB", method = "welch", compare = .AB)
  expect_equal(unique(unp$effect_type), "hedges_g")
})

test_that("only subjects measured in both groups are paired", {
  de <- .de_stats()
  sm <- as.data.frame(de$sample_meta)
  sm$Subject[8] <- "P9"          # the fourth B sample is someone else's
  de$sample_meta <- sm

  params <- list(comparisons = list(enabled = TRUE, entries = list(
    list(name = "AvB", method = "welch", paired = TRUE, subject = "Subject",
         compare = .AB))))
  expect_message(out <- runStats(de, params),
                 "3 subject\\(s\\) measured in both groups")

  r <- out$stats[out$stats$feature == "F2", ]
  expect_equal(r$n_g1, 3L)
  expect_equal(r$n_g2, 3L)
  expect_equal(r$p_value,
               stats::t.test(.F2[1:3], .F2[5:7], paired = TRUE)$p.value)
})

test_that("pairing that cannot be done is refused, not approximated", {
  de  <- .de_stats()
  ent <- function(...) list(comparisons = list(enabled = TRUE,
                                               entries = list(list(...))))

  expect_error(runStats(de, ent(name = "x", method = "anova", paired = TRUE,
                                subject = "Subject", compare = .ABC)),
               "needs exactly 2 levels")
  expect_error(runStats(de, ent(name = "x", method = "welch", paired = TRUE,
                                compare = .AB)),
               "also needs subject:")
  expect_error(runStats(de, ent(name = "x", method = "welch", paired = TRUE,
                                subject = "Animal", compare = .AB)),
               "subject column 'Animal' not in sample_meta")

  sm <- as.data.frame(de$sample_meta)
  sm$Subject[1:4] <- "P1"        # one subject, four observations
  de$sample_meta <- sm
  expect_error(runStats(de, ent(name = "x", method = "welch", paired = TRUE,
                                subject = "Subject", compare = .AB)),
               "not unique within a group")
})

test_that("p-values are adjusted within the family each row names", {
  out <- suppressMessages(runStats(.de_stats(), list(
    comparisons = list(enabled = TRUE, entries = list(
      list(name = "AvB", method = "welch",    compare = .AB),
      list(name = "AvB", method = "wilcoxon", compare = .AB),
      list(name = "AvC", method = "welch",    compare = .AC))))))
  st <- out$stats

  # One family per comparison and method: a parametric and a rank test of the
  # same hypothesis are two answers to one question, not six tests.
  expect_equal(sort(unique(st$adjust_family)),
               c("AvB | welch", "AvB | wilcoxon", "AvC | welch"))

  for(fam in unique(st$adjust_family)){
    rows <- st[st$adjust_family == fam, ]
    expect_equal(nrow(rows), 4)          # four features
    expect_equal(rows$p_adj, stats::p.adjust(rows$p_value, method = "BH"))
    expect_equal(rows$significant, rows$p_adj < 0.05)
  }

  # The threshold is the one configured, not a constant. F2 lands between
  # the two, so it is the row that changes its answer.
  thr <- function(x) suppressMessages(runStats(.de_stats(), list(
    sig_threshold = x,
    comparisons = list(enabled = TRUE, entries = list(
      list(name = "AvB", method = "welch", compare = .AB))))))$stats
  strict <- thr(0.05); loose <- thr(0.5)

  expect_equal(loose$significant, loose$p_adj < 0.5)
  expect_equal(strict$significant, strict$p_adj < 0.05)
  expect_true(loose$significant[loose$feature == "F2"])
  expect_false(strict$significant[strict$feature == "F2"])
})

test_that("correlations come back as the lower triangle, both methods", {
  out <- runStats(.de_stats(), list(correlations = list(
    enabled = TRUE, entries = list(
      list(name = "pairs", methods = c("pearson", "spearman"))))))
  cr <- out$correlations

  expect_equal(sort(unique(cr$method)), c("pearson", "spearman"))
  expect_true(all(cr$correlation == "pairs"))
  expect_true(all(cr$subset == "all"))

  # Four features give six pairs, and nothing is correlated with itself.
  expect_equal(nrow(cr[cr$method == "pearson", ]), 6)   # 4 features
  expect_false(any(cr$feature_a == cr$feature_b))

  pair_of <- function(d, m, f) d[d$method == m & d$feature_a %in% f &
                                 d$feature_b %in% f, ]
  pe <- pair_of(cr, "pearson", c("F1", "F2"))
  expect_equal(pe$n, 12L)
  expect_equal(pe$estimate, stats::cor(.F1, .F2))
  expect_equal(pe$p_value, stats::cor.test(.F1, .F2)$p.value)

  # Spearman goes through cor.test() so the p-value matches what a reader
  # gets checking one pair by hand, exact or asymptotic as cor.test decides.
  sp <- pair_of(cr, "spearman", c("F1", "F2"))
  expect_equal(sp$estimate, stats::cor(.F1, .F2, method = "spearman"))
  expect_equal(sp$p_value, suppressWarnings(
    stats::cor.test(.F1, .F2, method = "spearman")$p.value))

  # A pair with gaps is counted on the samples it has.
  expect_equal(pair_of(cr, "pearson", c("F1", "F3"))$n, sum(!is.na(.F3)))
})

test_that("correlations take subsets, and the diagonal only when asked", {
  out <- runStats(.de_stats(), list(correlations = list(
    enabled = TRUE, entries = list(
      list(name = "by_group", methods = "pearson", include_self = TRUE,
           subsets = list(list(Treatment = "A"), list(Treatment = "B")))))))
  cr <- out$correlations

  expect_equal(sort(unique(cr$subset)), c("Treatment=A", "Treatment=B"))
  expect_true(any(cr$feature_a == cr$feature_b))            # include_self
  expect_equal(nrow(cr[cr$subset == "Treatment=A", ]), 10)  # 6 pairs + 4 self

  # Each subset and method is its own family.
  expect_equal(sort(unique(cr$adjust_family)),
               c("by_group | Treatment=A | pearson",
                 "by_group | Treatment=B | pearson"))
})

test_that("linear models report the fit, with the reference level named", {
  out <- runStats(.de_stats(), list(linear_models = list(
    enabled = TRUE, entries = list(
      list(name = "adjusted", formula = "~ Treatment + Sex")))))
  lmt <- out$linear_models

  expect_equal(sort(unique(lmt$term_raw)),
               c("(Intercept)", "SexM", "TreatmentB", "TreatmentC"))
  expect_true(all(lmt$model == "adjusted"))
  expect_true(all(lmt$formula == "~ Treatment + Sex"))

  d   <- data.frame(F1 = .F1, Treatment = .GRP, Sex = .SEX)
  co  <- stats::coef(summary(stats::lm(F1 ~ Treatment + Sex, data = d)))
  row <- lmt[lmt$feature == "F1" & lmt$term_raw == "TreatmentB", ]

  expect_equal(row$estimate, unname(co["TreatmentB", "Estimate"]))
  expect_equal(row$std_error, unname(co["TreatmentB", "Std. Error"]))
  expect_equal(row$statistic, unname(co["TreatmentB", "t value"]))
  expect_equal(row$p_value, unname(co["TreatmentB", "Pr(>|t|)"]))

  # A coefficient means nothing without the level it is measured against.
  expect_match(row$term, "vs A")
  expect_match(lmt$term[lmt$term_raw == "(Intercept)"][1],
               "ref: Treatment=A; Sex=F", fixed = TRUE)

  # One family per term: the question asked of a treatment coefficient is not
  # the question asked of the intercept.
  expect_equal(length(unique(lmt$adjust_family)), 4)
  fam <- lmt[lmt$term_raw == "TreatmentB", ]
  expect_equal(fam$p_adj, stats::p.adjust(fam$p_value, method = "BH"))
})

test_that("a model without a formula is refused", {
  expect_error(runStats(.de_stats(), list(linear_models = list(
    enabled = TRUE, entries = list(list(name = "nf"))))),
    "formula` required")
})

test_that("nothing configured gives empty tables with the full schema", {
  out <- runStats(.de_stats(), list())

  expect_equal(nrow(out$stats), 0)
  expect_equal(nrow(out$correlations), 0)
  expect_equal(nrow(out$linear_models), 0)

  # The xlsx writer and both reports read these columns by name, so an empty
  # run has to carry them rather than hand back a bare frame.
  expect_true(all(c("comparison", "method", "feature", "p_value", "p_adj",
                    "significant", "adjust_family", "effect_size") %in%
                  colnames(out$stats)))
  expect_true(all(c("correlation", "subset", "method", "feature_a",
                    "feature_b", "n", "estimate", "p_adj") %in%
                  colnames(out$correlations)))
  expect_true(all(c("model", "formula", "feature", "term", "term_raw",
                    "estimate", "std_error", "p_adj") %in%
                  colnames(out$linear_models)))
})

test_that("a section switched off is skipped, and an unmigrated one caught", {
  de <- .de_stats()

  off <- runStats(de, list(comparisons = list(enabled = FALSE, entries = list(
    list(name = "AvB", method = "welch", compare = .AB)))))
  expect_equal(nrow(off$stats), 0)

  # An entry can be switched off on its own.
  one <- suppressMessages(runStats(de, list(comparisons = list(
    enabled = TRUE, entries = list(
      list(name = "AvB", method = "welch", compare = .AB),
      list(name = "AvC", enabled = FALSE, method = "welch",
           compare = .AC))))))
  expect_equal(unique(one$stats$comparison), "AvB")

  # The pre-entries: layout, where the section was a bare list.
  expect_error(runStats(de, list(comparisons = list(
    list(name = "AvB", method = "welch", compare = .AB)))),
    "must be a block with 'enabled:' and 'entries:'")
})

test_that("runStats says what it was handed", {
  expect_error(runStats(data.frame(a = 1), list()), "DatasetExperiment")
  expect_error(runStats(data.frame(a = 1), list()), "runStats")
})
