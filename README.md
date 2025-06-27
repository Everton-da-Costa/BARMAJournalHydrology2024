# BARMAJournalHydrology2024

This repository contains the R code and associated data for the scientific article:

**"Test inferences and link function selection in dynamic beta modeling of seasonal hydro-environmental time series with temporary abnormal regimes"**  
by Costa, E., Cribari-Neto, F., and Scher, V. T.  
Published in the *Journal of Hydrology*, Volume 638, 2024, 131489.

[**View Article on ScienceDirect**](https://doi.org/10.1016/j.jhydrol.2024.131489)

---

## 📚 Table of Contents

- [📄 Project Overview](#-project-overview)
- [📂 Repository Structure](#-repository-structure)
- [📦 Repository Contents](#-repository-contents)
- [🛠️ Setup and Usage](#️-setup-and-usage)
- [🎓 Citation](#-citation)
- [📬 Contact](#contact)

---

## 📄 Project Overview

This repository provides R scripts and data to replicate the time series analysis on the useful volume of three water reservoirs: Itaparica, Sobradinho, and Três Marias (ONS, 2024). The code implements dynamic beta modeling for seasonal hydro-environmental time series, addressing test inferences and link function selection in the presence of temporary abnormal regimes.

---

## 📂 Repository Structure

```plaintext
.
├── R/                  # Core R functions and scripts
│   ├── barma.R
│   ├── barma3ClassicalTests.R
│   ├── barma3ClassicalTestsAuxFun.R
│   ├── makeLinkStructure.R
│   └── simuBarma.R
│
├── data/               # Raw data files
│   ├── itaparica.csv
│   ├── sobradinho.csv
│   └── tres_marias.csv
│
├── analysis/           # Application scripts and Rmd sources
│   ├── application_itaparica.R
│   ├── application_sobradinho.R
│   ├── application_tres_marias.R
│   └── classical_tests_example.Rmd
│
├── output/             # Generated outputs (PDFs)
│   └── classical_tests_example.pdf
│
├── .gitignore          # Files/folders to ignore in Git
├── BARMAJournalHydrology2024.Rproj # RStudio project file
└── README.md           # Project overview
```

---

## 📦 Repository Contents

This section details the contents of the repository, organized by directory and file purpose.

*   `R/`
    *   `barma.R`: Contains core functions for the BARMA model, used by the main application scripts.
    *   `barma3ClassicalTests.R`: Main script implementing the three classical tests for misspecification.
    *   `barma3ClassicalTestsAuxFun.R`: Auxiliary/helper functions supporting the classical tests.
    *   `makeLinkStructure.R`: Functions for creating the link function structure for the models.
    *   `simuBarma.R`: Functions for simulating BARMA time series.
*   `data/`
    *   `itaparica.csv`: Useful volume data for the Itaparica reservoir.
    *   `sobradinho.csv`: Useful volume data for the Sobradinho reservoir.
    *   `tres_marias.csv`: Useful volume data for the Três Marias reservoir.
*   `analysis/`
    *   `application_itaparica.R`: Analysis script for replicating results for the Itaparica reservoir.
    *   `application_sobradinho.R`: Analysis script for replicating results for the Sobradinho reservoir.
    *   `application_tres_marias.R`: Analysis script for replicating results for the Três Marias reservoir.
    *   `classical_tests_example.Rmd`: R Markdown source file for the numerical example of the three classical tests.
*   `output/`
    *   `classical_tests_example.pdf`: The rendered PDF report generated from `classical_tests_example.Rmd`.
*   `.gitignore`: Specifies intentionally untracked files to ignore.
*   `BARMAJournalHydrology2024.Rproj`: RStudio project file for easy project management.
*   `README.md`: This overview document for the repository.

---

## 🛠️ Setup and Usage

To run the R scripts in this repository, you will need **R** installed.

### 1. Dependencies

This code requires several R packages. You can install all dependencies by running the following command in your R console:

```R
install.packages(c("tseries", "forecast", "zoo", "lbfgs", "moments", "rmarkdown"))
```

### 2. Last Tested Environment
The scripts were last successfully tested on:
*   **R version:** 4.4.2 ("Pile of Leaves")
*   **Platform:** x86_64-pc-linux-gnu (64-bit)

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

## 📬 Contact
For questions, suggestions, or issues related to the code, please contact:

Everton da Costa

📧 everto.cost@gmail.com.br
