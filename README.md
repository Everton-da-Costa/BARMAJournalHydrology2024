# BARMAJournalHydrology2024

[![DOI](https://img.shields.io/badge/DOI-10.1016/j.jhydrol.2024.131489-blue.svg)](https://doi.org/10.1016/j.jhydrol.2024.131489)


This repository contains the R package and associated data for the scientific article:

**"Test inferences and link function selection in dynamic beta modeling of seasonal hydro-environmental time series with temporary abnormal regimes"**  
by Costa, E., Cribari-Neto, F., and Scher, V. T.  
Published in the *Journal of Hydrology*, Volume 638, 2024, 131489.

[**View Article on ScienceDirect**](https://doi.org/10.1016/j.jhydrol.2024.131489)

---

## 📚 Table of Contents

- [🎯 Project Motivation](#-project-motivation)
- [✨ Key Features](#-key-features)
- [📂 Repository Structure](#-repository-structure)
- [🛠️ Installation](#️-installation)
- [🚀 How to Replicate the Analysis](#-how-to-replicate-the-analysis)
- [🎓 Citation](#-citation)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [📬 Contact](#-contact)

---

## 🎯 Project Motivation

Climate change has increased the frequency of extreme weather events like prolonged droughts. This poses a significant challenge for managing water resources, particularly for hydroelectric power plants. This project addresses the need for accurate modeling and forecasting of reservoir water levels, which are time series data naturally bounded between 0 and 1 .

We developed and validated a dynamic beta model (BARMA) to capture the unique characteristics of these time series: seasonality and abnormally dry periods. The analysis focuses on the useful volume of three major Brazilian reservoirs: Itaparica, Sobradinho, and Três Marias.

This work provides robust statistical tools for hydrologists and data scientists to:

1.  Perform accurate hypothesis tests for model validation.

2.  Select the best model configuration (i.e., the link function).

3.  Generate reliable in-sample predictions and out-of-sample forecasts.

---

## ✨ Key Features

This package provides a full toolchain for dynamic beta modeling. The key technical components include:

*   **Dynamic Beta (BARMA) Model:** The core model is implemented in `R/barma.R`, providing the main function for fitting Beta Autoregressive Moving Average models to doubly-bounded time series.

*   **Classical Hypothesis Tests:** A key contribution of the paper, implemented in `R/barma_classical_tests.R`. This script provides functions for performing Wald, Score, and Likelihood-Ratio tests, which are essential for model validation and inference.

*   **Core Estimation Engine:** The mathematical foundation of the model is implemented in a series of functions, including:
    *   `R/loglik_terms_ar.R` & `R/loglik_terms_ma.R`: For computing the log-likelihood terms.
    *   `R/score_vector_arma.R`: For computing the score vector (gradient).
    *   `R/inf_matrix_arma.R`: For computing the information matrix, crucial for standard errors and statistical tests.

*   **Tools for Link Function Selection:** Provides code to analyze and compare the performance of different link functions (e.g., logit, probit, cloglog) for the model.

*   **Full Reproducibility:** The `analysis/` directory contains scripts to replicate all key findings from the paper.

---

## 📂 Repository Structure

The repository is structured as a standard R package for clarity and reproducibility. This shows the main components; see the "Key Features" section above for details on the most important files.

```plaintext
.
├── R/                  # Source code for all R functions in the package.
├── data/               # Raw data files (.csv) for the reservoirs.
├── analysis/           # R scripts to replicate the paper's analysis.
├── man/                # R package documentation files.
├── inst/               # Additional files (e.g., CITATION, .bib).
├── output/             # Generated outputs (e.g., PDFs from R Markdown).
│
├── DESCRIPTION         # Package metadata and dependencies.
├── NAMESPACE           # Manages the package's namespace.
├── BARMAJournalHydrology2024.Rproj # RStudio project file.
├── LICENSE             # MIT License file.
└── README.md           # This file.
```

---

## 🛠️ Installation
This research compendium can be installed as an R package directly from GitHub. This is the recommended method as it handles all dependencies automatically.

First, ensure you have the remotes package. If not, install it from CRAN:

```R
if (!require("remotes")) {
  install.packages("remotes")
}
```
```R
remotes::install_github("everton-da-costa/BARMAJournalHydrology2024", dependencies = TRUE)
```

Last Tested Environment
The scripts were last successfully tested on:

*   **R version:** 4.4.2 ("Pile of Leaves")

*   **Platform:** x86_64-pc-linux-gnu (64-bit)

---

## 🚀 How to Replicate the Analysis

After installing the package, you can replicate the analyses presented in the paper.

1. Clone the Repository

First, clone this repository to get access to the analysis scripts.

```R
git clone https://github.com/everton-da-costa/BARMAJournalHydrology2024.git
cd BARMAJournalHydrology2024
```

2. Load the Package
Start a new R session within the project directory and load the package:

```R
library(BARMAJournalHydrology2024)
```

3. Run the Analysis Scripts
The analysis scripts are located in the analysis/ directory. You can run them directly in your R console.

For example, to run the analysis for the Sobradinho reservoir:

```R
# This script applies the model to the Sobradinho data
source("analysis/application_sobradinho.R")
```
Similarly, you can run the analysis for the other reservoirs:

```R
# source("analysis/application_itaparica.R")
# source("analysis/application_tres_marias.R")
```

4. Generate the PDF Report


The paper includes an example of classical tests. You can generate the PDF report for this example by running the rendering script:

```R
# This will execute the .Rmd file and save the PDF in the output/ directory
source("analysis/render_classical_tests_example.R")
```

---

## 🎓 Citation

If you use this code or data in your research, please cite the original article:

```bibtex
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
Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

## 📄 License
This project is licensed under the MIT License. See the LICENSE file for details.

## 📬 Contact
For questions, suggestions, or issues related to the code, please contact:

Everton da Costa

📧 everto.cost@gmail.com
