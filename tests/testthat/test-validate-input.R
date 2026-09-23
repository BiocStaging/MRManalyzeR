fm_ok = function()
  data.frame(Compound        = c("A", "B"),
             Processing_name = c("a", "b"),
             Report          = c("YES", "YES"),
             stringsAsFactors = FALSE)

sm_ok = function()
  data.frame(Name    = c("S1", "S2"),
             Include = c("YES", "YES"),
             stringsAsFactors = FALSE)

test_that("a well-formed pair of sheets passes", {
  v = validateInput(fm_ok(), sm_ok())
  expect_length(v$errors, 0)
  expect_length(v$warnings, 0)
})

test_that("report_col is optional, but its absence is announced", {
  # A feature table curated outside a TargetLynx workflow - the usual case for
  # data_source = "matrix" - has no Report column. That must not be an error,
  # because there is no exclusion information to honour.
  fm = fm_ok()[, c("Compound", "Processing_name")]

  v = validateInput(fm, sm_ok())
  expect_length(v$errors, 0)
  expect_match(v$warnings, "no 'Report' column", all = FALSE)

  # Silence would be wrong too: a typo in report_col looks identical from here
  # and would quietly reinstate features that were meant to be dropped.
  expect_true(any(v$summary$status == "warning"))
})

test_that("duplicate compound names are still caught without report_col", {
  fm = data.frame(Compound        = c("A", "A"),
                  Processing_name = c("a", "b"),
                  stringsAsFactors = FALSE)
  v = validateInput(fm, sm_ok())
  expect_match(v$errors, "duplicate", all = FALSE)
})

test_that("compound and processing-name columns are still required", {
  fm = data.frame(Something = c("A", "B"), stringsAsFactors = FALSE)
  v  = validateInput(fm, sm_ok())
  expect_match(v$errors, "missing required column", all = FALSE)
})

test_that("an included injection's order has to be a number", {
  sm <- sm_ok()
  sm$Injection_order <- c("_48", "2")
  v <- validateInput(fm_ok(), sm, injection_order_head = "Injection_order")
  expect_match(v$errors, "non-numeric Injection_order", all = FALSE)
  expect_match(v$errors, "'_48'", all = FALSE)

  sm$Injection_order <- c(1, 2)
  v <- validateInput(fm_ok(), sm, injection_order_head = "Injection_order")
  expect_length(v$errors, 0)

  # Excluded injections are not checked.
  sm$Injection_order <- c("x", "2")
  sm$Include <- c("NO", "YES")
  v <- validateInput(fm_ok(), sm, injection_order_head = "Injection_order")
  expect_length(v$errors, 0)
})
