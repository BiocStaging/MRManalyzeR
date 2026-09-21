# ---------------------------------------------------------------------------
# Build the synthetic example dataset
#
#   inst/extdata/example_synthetic.xlsx  - Skyline-layout workbook
#   inst/extdata/example_synthetic.RDS   - the processed DatasetExperiment
#
# NOTHING HERE IS REAL DATA. Every value is simulated from a fixed seed. The
# bundled TargetLynx workbook (example_data.xlsx) is the real dataset; this one
# exists because that one cannot exercise the whole workflow: it is a single
# chromatographic batch, one reconstitution volume for every vial, and has no
# per-sample amount to normalise to, so normaliseMatrix(), correctBatch() and
# adjustConcentration() would all be no-ops on it.
#
# What is deliberately planted, so the documentation can show each step
# recovering something known:
#
#   * a TREATMENT EFFECT on 5 of the 20 analytes (recorded in
#     feature_metadata$Simulated_log2FC; the other 15 are flat), so the stats
#     stage can be checked against the truth;
#   * a per-analyte RESPONSE SHIFT between the two chromatographic batches,
#     so correctBatch() has something to remove and a PCA coloured by batch
#     separates before correction and not after;
#   * QC injections drawn with ~8 percent noise against ~35 percent between
#     animals, so CV_QC, CV_sample and their ratio are realistic;
#   * signal proportional to protein_ug (divided out by normaliseMatrix) and
#     inversely proportional to the reconstitution volume, which differs
#     between batches;
#   * a per-analyte LOD that masks the bottom of each distribution, with two
#     analytes deliberately sparse.
#
# The design is BALANCED ACROSS BATCHES (5 animals per treatment x genotype
# cell in each batch) - otherwise the batch effect and the treatment effect
# would be confounded and neither could be demonstrated.
#
# Run from the package root, after devtools::load_all():
#   source("data-raw/make_example_synthetic.R")
# ---------------------------------------------------------------------------

set.seed(20260721)

out_xlsx <- "inst/extdata/example_synthetic.xlsx"
out_rds  <- "inst/extdata/example_synthetic.RDS"

# --- Analytes ---------------------------------------------------------------
n_comp   <- 20L
analyte  <- sprintf("Analyte_%02d", seq_len(n_comp))
pathway  <- rep(c("COX", "LOX", "CYP"), length.out = n_comp)
precursor <- rep(c("AA", "EPA", "DHA", "LA"), length.out = n_comp)

# Baseline abundance spans three orders of magnitude, as a real panel does.
base <- 10^stats::runif(n_comp, 1.5, 4.5)

# Planted treatment effect: 3 up, 2 down, the rest flat.
effect_idx  <- c(2L, 5L, 9L, 13L, 18L)
effect_fc   <- c(2.5, 2.0, 1.8, 0.50, 0.60)
fold_change <- rep(1, n_comp)
fold_change[effect_idx] <- effect_fc

# Two analytes sit near the detection limit and go missing in many samples.
sparse_idx <- c(7L, 16L)

# --- Injection sequence -----------------------------------------------------
# One batch = a blank, then five blocks of (4 samples + 1 QC), then a blank.
one_batch <- c("Blank", rep(c(rep("Sample", 4), "QC"), 5), "Blank")
type      <- rep(one_batch, 2)
chrom     <- rep(c("B1", "B2"), each = length(one_batch))
n_inj     <- length(type)

is_sample <- type == "Sample"
is_qc     <- type == "QC"
is_blank  <- type == "Blank"

# --- Study design, balanced within each batch -------------------------------
cells <- expand.grid(Treatment = c("PBS", "HDM"), Type = c("WT", "KO"),
                     stringsAsFactors = FALSE)
one_design <- cells[rep(seq_len(nrow(cells)), each = 5), ]   # 20 per batch
design <- rbind(one_design[sample(nrow(one_design)), ],
                one_design[sample(nrow(one_design)), ])

treatment <- rep(NA_character_, n_inj)
genotype  <- rep(NA_character_, n_inj)
treatment[is_sample] <- design$Treatment
genotype[is_sample]  <- design$Type

# --- Per-injection metadata -------------------------------------------------
# Protein drives the measured signal and is divided out by normaliseMatrix().
protein <- rep(250, n_inj)
protein[is_sample] <- round(stats::rnorm(sum(is_sample), 250, 45))

# Batch 2 was reconstituted in twice the volume, so its extracts are twice as
# dilute before any correction.
recon_vol <- ifelse(chrom == "B1", 50, 100)

sample_meta <- data.frame(
  Name             = sprintf("SYN_%03d", seq_len(n_inj)),
  Sample_ID        = ifelse(is_sample,
                            sprintf("animal_%02d", cumsum(is_sample)),
                            sprintf("%s_%02d", type, seq_len(n_inj))),
  Injection_order  = seq_len(n_inj),
  Treatment        = treatment,
  Type             = genotype,
  Sex              = ifelse(is_sample,
                            rep(c("F", "M"), length.out = n_inj), NA),
  protein_ug       = protein,
  starting_vol_uL  = 200,
  sample_volume_uL = recon_vol,
  cal_vol_uL       = 50,
  sample_IS_vol_uL = 10,
  cal_IS_vol_uL    = 10,
  Chrom_Batch      = chrom,
  Extraction_Batch = ifelse(chrom == "B1", "E1", "E2"),
  Include          = "YES",
  Sample_type      = type,
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

# --- Simulate the measured signal -------------------------------------------
areas <- matrix(NA_real_, nrow = n_inj, ncol = n_comp,
                dimnames = list(sample_meta$Name, analyte))

# Between-batch response shift, one factor per analyte.
batch_shift <- exp(stats::rnorm(n_comp, 0, 0.25))

for(j in seq_len(n_comp)){
  v <- numeric(n_inj)

  # Biological samples: treatment effect (HDM only), animal-to-animal spread,
  # and signal proportional to the protein loaded.
  fc <- ifelse(treatment[is_sample] == "HDM", fold_change[j], 1)
  v[is_sample] <- base[j] * fc *
    exp(stats::rnorm(sum(is_sample), 0, 0.35)) *
    (protein[is_sample] / 250)

  # Pooled QCs: identical material, so only instrument noise.
  v[is_qc] <- base[j] * exp(stats::rnorm(sum(is_qc), 0, 0.08))

  # Blanks: carryover at a fraction of a percent.
  v[is_blank] <- base[j] * 0.003 * exp(stats::rnorm(sum(is_blank), 0, 0.5))

  # Dilution by reconstitution volume, then the batch response shift.
  v <- v * (50 / recon_vol)
  v[chrom == "B2"] <- v[chrom == "B2"] * batch_shift[j]

  areas[, j] <- v
}

# --- Per-analyte LOD --------------------------------------------------------
# Masks the bottom few percent of each analyte; the two sparse analytes get a
# much higher floor so they go missing in roughly a third of samples.
probs <- rep(0.04, n_comp)
probs[sparse_idx] <- 0.35
lod <- vapply(seq_len(n_comp), function(j)
  stats::quantile(areas[is_sample, j], probs = probs[j], na.rm = TRUE),
  numeric(1))

areas <- round(areas, 1)
lod   <- round(lod, 1)

# --- Feature metadata -------------------------------------------------------
feature_meta <- data.frame(
  Processing_name   = analyte,      # Skyline molecules are named canonically
  Compound          = analyte,
  Report            = "YES",
  Comment           = NA_character_,
  LOD               = lod,
  Enzymatic_pathway = pathway,
  PUFA              = precursor,
  RT                = round(stats::runif(n_comp, 2, 12), 2),
  Simulated_log2FC  = round(log2(fold_change), 3),
  stringsAsFactors  = FALSE,
  check.names       = FALSE
)

# --- Skyline layout: molecules down the rows, samples across the columns ----
skyline_data <- cbind(
  data.frame(Molecule = analyte, check.names = FALSE,
             stringsAsFactors = FALSE),
  as.data.frame(t(areas), check.names = FALSE))

openxlsx::write.xlsx(
  list(feature_metadata = feature_meta,
       sample_metadata  = sample_meta,
       skyline_data     = skyline_data),
  file = out_xlsx, overwrite = TRUE)

message("[make_example_synthetic] wrote ", out_xlsx, " - ",
        n_inj, " injections x ", n_comp, " analytes")

# --- Processed DatasetExperiment, for the man-page examples ------------------
# Shipping the processed object means every plot example is three lines and
# needs no workbook parsing during R CMD check.
de <- processDataset(
  feature_meta, sample_meta,
  xlsx_path        = out_xlsx,
  data_source      = "skyline",
  data_tab_names   = "skyline_data",
  signal_filter    = "LOD",
  blank_filter     = 3,
  normalize        = "protein_ug",
  adjust_conc      = TRUE,
  starting_vol_col = "starting_vol_uL",
  replace_MVs      = 0.5,
  batch_correction = TRUE,
  bc_qc_label      = "QC",
  bc_factor_name   = "Sample_type",
  bc_header        = "Chrom_Batch",
  blank_head       = "Sample_type")[[1]]

de <- addCVMetrics(de, sample_type_head = "Sample_type",
                     qc_label = "QC", sample_labels = "Sample")

saveRDS(de, out_rds)

message("[make_example_synthetic] wrote ", out_rds, " - ",
        nrow(as.data.frame(de$data)), " x ",
        ncol(as.data.frame(de$data)))

# --- What was planted, for the record ---------------------------------------
message("[make_example_synthetic] planted treatment effects:")
print(feature_meta[feature_meta$Simulated_log2FC != 0,
                   c("Compound", "Enzymatic_pathway", "Simulated_log2FC")])
