# MRManalyzeR <img src="man/figures/logo.png" align="right" height="139" alt="" />

Reproducible post-acquisition processing of targeted LC-MS/MS metabolomics
results exported from Waters TargetLynx, Skyline and plain matrix formats
(extendable to other vendors). It reads the quantitative result tables,
integrates them with sample and feature metadata, applies configurable
quality-control and matrix-processing steps, and returns an analysis-ready
`struct::DatasetExperiment` together with self-contained data-quality and
statistical reports.

Each processing and analysis step is available as an exported R function
operating on that object. Complete workflows can additionally be run from a
single YAML configuration, without writing a driver script.

MRManalyzeR takes vendor software outputs. It does **not** perform
chromatographic peak detection, peak integration or calibration-curve fitting
from raw mass-spectrometry data.

## Scope

Built for **targeted quantitative assays** — methods of roughly >20 MRM
transitions, where the analytes are known compounds carrying pathway or class
annotation, and the experiment is structured around defined target lists,
whether testing a specific hypothesis or conducting targeted exploratory
profiling.

Two consequences follow, and they are deliberate:

* **Reported values are commonly internal-standard normalised.** Normalisation
  and batch correction are therefore intentionally limited. Signal-drift
  modelling is available in untargeted metabolomics tools.
* **Feature metadata carries biology.** Compound class and enzymatic pathway
  are first-class: they group and colour heatmap annotations, PCA loadings
  and correlation blocks throughout.

The primary input is TargetLynx. Skyline exports and plain sample-by-analyte
matrices from any other software are supported through the same workflow, and
the metadata schema is configurable for other laboratories.

### Where MRManalyzeR fits

| Task | Tool |
|---|---|
| Raw untargeted LC-MS feature detection and alignment | `xcms` |
| General metabolomics peak-matrix preprocessing | `pmp` |
| Metabolite annotation-table workflows | `MetMashR` |
| Targeted assay processing, QC and reporting | **MRManalyzeR** |

---

## Installation

### Bioconductor dependency (do this first)

MRManalyzeR depends on `struct`, which is distributed through Bioconductor
rather than CRAN, so install it with `BiocManager` before installing
MRManalyzeR from GitHub:

```r
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("struct")
```

(`structToolbox` and `pmp` are **not** required. Earlier versions of this
README listed them; they are no longer used.)

### Latest release

```r
# install.packages("remotes")
remotes::install_github("MJS-708/MRManalyzeR", ref = "main")
```

Development branch:

```r
remotes::install_github("MJS-708/MRManalyzeR", ref = "dev")
```

Local install (from a cloned working copy):

```r
devtools::install("path/to/MRManalyzeR")
```

---

## Quick start — bundled example

Run the workflow against a small bundled oxylipin dataset to verify the
install and produce both reports. No paths to configure:

```r
res <- MRManalyzeR::runExample()
res$output_directory   # where the xlsx + html reports were written
res$datasetExperiment  # the processed dataset
res$stats_tables       # comparisons, correlations, linear models
```

`runExample()` copies the bundled YAML, points it at
`inst/extdata/example_data.xlsx` and a fresh `tempdir()`, runs
`runMRManalyzeR()`, and (in interactive sessions) opens the rendered
HTML reports in the RStudio viewer / default browser.

### The returned object

Every entry point returns a list whose `datasetExperiment` element is a
`struct::DatasetExperiment` holding four aligned pieces:

| Slot | Contents |
|---|---|
| `data` | the sample × feature matrix, in the reported datatype |
| `sample_meta` | one row per injection, plus everything the YAML added |
| `variable_meta` | one row per compound: class, pathway, CV metrics, exclusion reasons |
| `name` / `description` | provenance carried through to the reports |

Because it is a `DatasetExperiment`, the result drops straight into other
Bioconductor tooling rather than being a private format.

### Example dataset

A subset of the BAL fluid oxylipin panel from **Kolmert *et al.* (2018),
*Prostaglandins & Other Lipid Mediators* 137, 11–18.**
DOI: <https://doi.org/10.1016/j.prostaglandins.2018.05.005>.

The bundled xlsx (`inst/extdata/example_data.xlsx`) contains a
TargetLynx-style workbook with `feature_metadata`, `sample_metadata` and
`lcms_data_*` sheets; the bundled YAML (`inst/extdata/example_config.yml`)
is annotated and ready to copy as a starting point for your own studies.

A fully synthetic two-batch dataset (`example_synthetic.RDS`) is also
bundled, and is what every function's `@examples` block runs on.

---

## Step-by-step API

The YAML route below is a convenience layer over exported functions. Any
step can be called directly on a `DatasetExperiment`:

```r
library(MRManalyzeR)

de <- assembleDataset(readTargetLynx("study.xlsx", datatype = "Area"),
                       feature_meta, sample_meta)
de <- filterBlanks(de, blank_filter = 3)
de <- correctBatch(de, batch_head = "Chrom_Batch")
de <- normaliseMatrix(de, column = "protein_mg")
de <- imputeMissing(de, scalar = 0.2)

qc  <- addCVMetrics(de)
pca <- runPCA(de, transform = "log2", scale = TRUE)
plotPCA(pca, de, colour_by = "Treatment")
```

See `vignette("MRManalyzeR")` for the full walkthrough, and `?MRManalyzeR`
for the four stages the functions are grouped into: **data parse →
peak-matrix processing → QC check → stats**.

---

## Usage on your own data

### Step 1 — prepare the data xlsx

The input workbook needs at least three sheets (all sheet names are
configurable):

* **`feature_metadata`** — one row per compound. Required: `Processing_name`
  (the vendor's identifier, used to join to the source table) and `Compound`
  (display name). Optional `Report` column, `YES` / `NO`, indicating whether
  the feature passed prior chromatographic review in TargetLynx or Skyline
  and should be included downstream; if the column is absent every feature is
  reported. MRManalyzeR does not itself inspect chromatographic peak shape —
  that judgement stays in the vendor software.
* **`sample_metadata`** — one row per acquisition. The acquisition name must
  match the vendor export exactly. Add any study metadata you want available
  for QC plots and statistics (e.g. `Injection_order`, `Chrom_Batch`, `Sex`,
  `Treatment`, `Sample_type`).
* **`lcms_data_1`** (and optionally `lcms_data_2`, …) — the TargetLynx
  summary table. All `lcms_data_*` sheets must share the same feature set.
  Skyline and generic-matrix inputs use their own sheets instead.

### Step 2 — copy and edit the YAML template

```r
file.copy(
  system.file("extdata", "example_config.yml", package = "MRManalyzeR"),
  "config.yml"
)
```

Then edit `config.yml`. The key blocks are:

| Block | Purpose |
|---|---|
| `paths:` | input xlsx location, output directory, output filename stem |
| `PeakMatrixProcessing:` | data source and datatype, signal / blank / detection filters, batch correction, normalisation, IS and volume adjustment, missing-value imputation |
| `data_quality_report:` | per-compound measurement vs injection order, Q-Q normality, QC PCA |
| `stats_report:` | configured group comparisons and multigroup tests, boxplots, correlations, linear models, PCA, volcano, ion ratios, heatmaps |

The statistical method is chosen explicitly per comparison (`welch`,
`student`, `wilcoxon` for two groups; `anova` or `kruskal` for three or
more) — nothing is selected automatically on your behalf.

Every option is documented inline in the bundled YAML — read it once
before editing.

### Step 3 — run

```r
res <- MRManalyzeR::runMRManalyzeR("config.yml")
```

This writes, into `paths.result_dir`:

* `<fn>_<datatype><suffix>.xlsx` — processed matrix, feature metadata and
  sample metadata; with `replace_MVs` set the matrix is imputed and a
  `matrix_measured` sheet holds the values as measured
* `<fn>_<datatype><suffix>_stats.xlsx` — a `key` sheet defining every column,
  then `stats`, `correlations`, `linear_models`, `ion_ratios` and `summary`
* `<fn>_<datatype><suffix>_data_quality_report.html`
* `<fn>_<datatype><suffix>_stats_report.html`
* `<fn>_<datatype><suffix>.RDS` and `.txt` — the dataset and the resolved
  parameters, plus `_measured.RDS` when `replace_MVs` is set

The reports and statistics use the measured values; set `use_imputed: True`
in a report's block to give it the imputed matrix instead. The PCAs use the
imputed values when `replace_MVs` is set, and impute for themselves otherwise.

The bundled config sets `suffix: _`, which is why the shipped examples read
`example_data_ng_mL__stats.xlsx` with two underscores. Set `suffix: ""` for
single ones.

To re-render reports from an already-processed RDS (skipping the xlsx
ingest), set `PeakMatrixProcessing.execute: False` and re-run — `fn`,
`datatype` and `suffix` must then match the stored file. The reports read the
`_measured.RDS` beside it when there is one.

---

## Combining multiple acquisitions

To merge separately-acquired panels (e.g. GOM, cysLT, SPM, or PL-pos /
PL-neg / SL) sharing biological samples, use the combine workflow:

```r
file.copy(
  system.file("extdata", "example_combine_config.yml", package = "MRManalyzeR"),
  "combine_config.yml"
)
res <- MRManalyzeR::runMRManalyzeRCombine("combine_config.yml")
```

The combine YAML lists per-panel paths (`.RDS` or `.xlsx`), optional QC
remappings (e.g. SL `QC1..5` ↔ PL `QC1, QC2, QC5, QC7, QC8`), feature
prefixing, and the same `stats_report:` block as a regular run.
`PeakMatrixProcessing` and the data quality report are skipped — the
inputs are already-processed matrices.

---

## Reports

Both reports are self-contained HTML with tabset navigation and (by
default) plotly hover tooltips that surface `Sample_ID` on points and
boxplot jitter. Toggle off with `interactive_plots: False` in either
report block of the YAML.

| Report | Sections |
|---|---|
| `data_quality_report` | sample summary; per-compound measurement vs injection order, in the selected datatype; per-compound Q-Q normality; QC PCA (one tab per colouring) |
| `stats_report` | global summaries; per-comparison statistical tables; per-compound boxplots (faceted across comparisons + by free-form factor); per-comparison PCA (scores + loadings); volcano; ion ratios (e.g. EpOME/DiHOME for sEH activity); samples × features heatmap; correlations; linear models |

Every figure in both reports is drawn by an exported `plot_*()` function, so
any panel can be reproduced, restyled or subset in one line without editing
an R Markdown file.

---

## License

See `LICENSE`.
