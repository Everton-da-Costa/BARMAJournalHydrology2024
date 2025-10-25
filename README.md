betaARMA: A Package for Beta Autoregressive Moving Average Models
================
Everton da Costa

# BARMAJournalHydrology2024

[![DOI](https://img.shields.io/badge/DOI-10.1016/j.jhydrol.2024.131489-blue.svg)](https://doi.org/10.1016/j.jhydrol.2024.131489)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/Everton-da-Costa/BARMAJournalHydrology2024/actions/workflows/R.yml/badge.svg)](https://github.com/Everton-da-Costa/BARMAJournalHydrology2024/actions/workflows/R.yml)

This repository contains the R package and associated data for the
scientific article:

**“Test inferences and link function selection in dynamic beta modeling
of seasonal hydro-environmental time series with temporary abnormal
regimes”** by Costa, E., Cribari-Neto, F., and Scher, V. T.  
Published in the *Journal of Hydrology*, Volume 638, 2024, 131489.

[**View Article on
ScienceDirect**](https://doi.org/10.1016/j.jhydrol.2024.131489)

------------------------------------------------------------------------

## 📚 Table of Contents

- [🎯 Project Motivation](#-project-motivation)
- [✨ Key Features](#-key-features)
- [📂 Repository Structure](#-repository-structure)
- [🛠️ Installation](#️-installation)
- [🚀 Getting Started & Examples
  (Vignettes)](#-getting-started--examples-vignettes)
- [🎓 Citation](#-citation)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [📬 Contact](#-contact)

------------------------------------------------------------------------

## 🎯 Project Motivation

Climate change has increased the frequency of extreme weather events
like prolonged droughts. This poses a significant challenge for managing
water resources, particularly for hydroelectric power plants. This
project addresses the need for accurate modeling and forecasting of
reservoir water levels, which are time series data naturally bounded
between 0 and 1.

We developed and validated a dynamic beta model ($\beta$ARMA) to capture
the unique characteristics of these time series: seasonality and
abnormally dry periods. The analysis focuses on the useful volume of the
Itaparica reservoir in Brazil.

This work provides statistical tools for hydrologists and data
scientists to:

1.  Perform accurate hypothesis tests for model validation.
2.  Select the best model configuration (i.e., the link function).
3.  Generate reliable in-sample predictions and out-of-sample forecasts.

------------------------------------------------------------------------

## ✨ Key Features

This package provides a full toolchain for dynamic beta modeling. The
key technical components include:

- **Dynamic Beta ($\beta$ARMA) Model:** The core model is implemented in
  `R/barma.R`, providing the main function for fitting Beta
  Autoregressive Moving Average models.
- **Classical Hypothesis Tests:** A key contribution of the paper,
  implemented in `R/barma_classical_tests.R`, providing functions for
  Wald, Score, and Likelihood-Ratio tests.
- **Core Estimation Engine:** The mathematical foundation is implemented
  in a series of functions for computing the log-likelihood
  (`loglik_terms_ar.R`, `loglik_terms_ma.R`), score vector
  (`score_vector_arma.R`), and information matrix (`inf_matrix_arma.R`).
- **Vignettes as Case Studies:** Two detailed vignettes serve as
  practical guides and complete applications of the methodology:
  - An **end-to-end application** forecasting reservoir levels.
  - A **technical deep-dive** into the implementation of statistical
    hypothesis tests.

------------------------------------------------------------------------

## 📂 Repository Structure

The repository is structured as a standard R package for clarity and
reproducibility.

``` plaintext
.
├── R/                  # Source code for all R functions.
├── data/               # Processed data included in the package (.rda).
├── data-raw/           # Raw data and scripts used to process it.
├── vignettes/          # Detailed tutorials and case studies (.Rmd).
├── man/                # R package documentation files for functions.
├── inst/               # Additional files (e.g., CITATION, REFERENCES.bib).
├── reports/            # Pre-rendered HTML reports for archival and reproducibility.
├── DESCRIPTION         # Package metadata and dependencies.
├── NAMESPACE           # Manages the package's namespace.
├── LICENSE             # MIT License file.
└── README.md           # This file.
```

------------------------------------------------------------------------

## Code of Conduct

Please note that the BARMAJournalHydrology2024 project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

## 🛠️ Installation

This research compendium can be installed as an R package directly from
GitHub. This is the recommended method as it handles all dependencies
automatically.

First, ensure you have the remotes package. If not, install it from
CRAN:

``` r
if (!require("remotes")) {
  install.packages("remotes")
}
```

Then, install the package from GitHub:

``` r
remotes::install_github("everton-da-costa/BARMAJournalHydrology2024", 
                        dependencies = TRUE,
                        build_vignettes = TRUE)
```

Last Tested Environment The scripts were last successfully tested on:

- **R version:** 4.4.2 (“Pile of Leaves”)

- **Platform:** x86_64-pc-linux-gnu (64-bit)

------------------------------------------------------------------------

## 🚀 Getting Started & Examples (Vignettes)

The best way to understand and replicate the analysis is through the
package vignettes, which provide detailed, narrated code examples.

### 1. View Pre-Rendered Reports (Recommended)

This is the fastest way to see the full analysis. You can view the pre-rendered HTML reports directly in your browser without installing the package.

* **Portfolio Case Study: [View `reservoir_itaparica` Report (HTML)](reports/reservoir_itaparica_2025-10-20.html)**

    > An end-to-end data science project demonstrating how to forecast the Itaparica reservoir water levels. It covers exploratory data analysis, feature engineering, model fitting, and comparison against benchmarks.

* **Technical Deep-Dive: [View `simulated_ts_classical_tests` Report (HTML)](reports/simulated_ts_classical_tests_2025-10-20.html)**

    > A detailed walkthrough of the implementation and validation of the Wald, Score, and Likelihood Ratio tests, replicating the simulation study from the original paper.

### 2. Run Locally (After Installation)

After installation, you can also see all available vignettes and run them interactively from your R console:
    
```R
# 1. Lists all tutorials for this package
browseVignettes("BARMAJournalHydrology2024")

# 2. Open a specific vignette
vignette("reservoir_itaparica", package = "BARMAJournalHydrology2024")
vignette("simulated_ts_classical_tests", package = "BARMAJournalHydrology2024")
```

------------------------------------------------------------------------

## 🎓 Citation

If you use this code or data in your research, please cite the original
article:

``` bibtex
@article{Costa+Cribari+Scher_2024,
  title     = {Test inferences and link function selection in dynamic beta modeling of seasonal hydro-environmental time series with temporary abnormal regimes},
  author    = {Costa, E. and Cribari-Neto, F. and Scher, V. T.},
  journal   = {Journal of Hydrology},
  volume    = {638},
  pages     = {131489}, 
  year      = {2024},
  doi       = {10.1016/j.jhydrol.2024.131489}
}
```

## 🤝 Contributing

Contributions are welcome! If you find any issues or have suggestions
for improvements, please open an issue or submit a pull request.

## 📄 License

This project is licensed under the MIT License. See the LICENSE file for
details.

## 📬 Contact

For questions, suggestions, or issues related to the code, please
contact:

Everton da Costa

📧 <everto.cost@gmail.com>
