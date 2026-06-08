# Geopolitical Risk and Gold Return Predictability across Quantile States

This repository provides the replication code for the study:

**Geopolitical Risk and Gold Return Predictability across Quantile States: Quantile-on-Quantile Regression with Block Bootstrap and Scenario Forecasts**

The workflow estimates baseline lagged quantile regressions, main and interactive multivariate quantile-on-quantile regression (mQQR) models, bootstrap and bandwidth robustness checks, and real-time/scenario-based forecasts.

## Repository structure

```text
.
├── R/
│   └── 01_reproduce_empirica_mqqr.R
├── data/
│   └── data_template.csv
├── outputs/
├── README.md
└── .gitignore
```

## Data format

Place the study dataset in the `data/` folder using one of these names:

- `data/rd.xlsx` with a sheet named `data`, or
- `data/rd.csv`

The dataset must contain the following variables:

| Required column | Description |
|---|---|
| `date` | Monthly date; accepted formats include `YYYY-MM-DD`, `YYYY-MM`, `YYYYMmm`, or Excel dates |
| `Gold Price` | Monthly gold price level |
| `Gold Return` | Monthly gold return in percent |
| `GPR` | Geopolitical risk index |
| `Dollar Index` | Dollar index |
| `VIX` | Volatility index |

The script converts the variables into standardized internal names and constructs lagged predictors.

## Main output folders

Running the script creates:

```text
outputs/01_Manuscript_Tables_Figures
outputs/02_Appendix_Tables_Figures
outputs/03_Model_Objects
outputs/04_Data_Audit
```

## How to run

Open R or RStudio from the repository root and run:

```r
source("R/01_reproduce_empirica_mqqr.R")
```

## Required R packages

The script checks and installs missing CRAN packages automatically. The main packages are:

- `readxl`, `openxlsx`, `dplyr`, `tidyr`, `stringr`, `lubridate`, `zoo`
- `quantreg`, `mqqr`, `QuantileOnQuantile`
- `ggplot2`, `patchwork`, `ragg`, `plot3D`, `scales`

## Notes

- The script uses lagged predictors to support the predictability design.
- If realized observations after the last available month are not included, the forecasting section is implemented as real-time and scenario-based forecasting rather than ex-post forecast evaluation.
- Full cell-level surfaces are exported to Excel files and the most important figures are saved at 600 dpi.
