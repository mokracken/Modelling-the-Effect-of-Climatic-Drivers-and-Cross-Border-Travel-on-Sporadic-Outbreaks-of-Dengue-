# Modelling the Effect of Climatic Drivers and Cross-Border Travel on Sporadic Outbreaks of Dengue in Mauritius

## Overview

This repository contains the analytical scripts and selected model outputs used to examine associations between climatic drivers, cross-border travel, and dengue incidence in Mauritius.

The study applies Distributed Lag Non-linear Models (DLNMs) within a Generalised Linear Mixed Model (GLMM) framework to investigate delayed and non-linear associations between environmental exposures, international travel, and dengue risk in a non-endemic setting.

## Repository Contents

### 1. Climate Data Extraction

**`ERA-5_Temp_extraction.js`**  
Google Earth Engine (GEE) script for extracting temperature data from ERA5-Land.

**`CHIRPS_extraction.js`**  
GEE script for extracting rainfall data from the CHIRPS dataset.

These scripts can be run in the Google Earth Engine Code Editor. Users must configure the geographical boundaries and study period before executing the extraction.

### 2. DLNM–GLMM Analysis (R)

The analysis is organised into three sequential scripts.

**Script 1 — `DLNM_GLMM_Final_Grid_script.R`**

Performs the DLNM grid search to identify candidate combinations of maximum lag periods and degrees of freedom for temperature, rainfall, vapour pressure deficit (VPD), and resident arrivals.

- Preprocesses epidemiological, climate, and travel data.
- Constructs DLNM cross-basis functions.
- Fits negative-binomial GLMMs using `glmmTMB`.
- Evaluates model convergence, AIC, BIC, and Pearson dispersion.
- Exports grid-search results and ranked candidate models.

**Script 2 — `DLNM_Model_Selection_and_Diagnostics.R`**

Screens the grid-search results and evaluates potentially competitive models.

- Automatically identifies Tmin, Tmean, or Tmax from the selected grid-results CSV.
- Verifies grid specifications and model identifiers.
- Refits competitive models and checks reproducibility against grid-search estimates.
- Evaluates model convergence, covariance matrices, relative risks, and confidence intervals.
- Performs DHARMa and temporal residual diagnostics.
- Exports model-selection and diagnostic results.

**Script 3 — `DLNM_GLMM_Final_Script_Rev.R`**

Fits and evaluates the selected DLNM–GLMM specification and generates the final analytical outputs.

- Fits the selected temperature, rainfall, VPD, and resident-arrival specification.
- Estimates cumulative and lag-specific relative risks (RRs) with 95% confidence intervals.
- Generates exposure–response and lag–response curves.
- Produces model summaries, figures, and tables for reporting.

### 3. Grid-Search Results and Model Outputs

The repository includes:

- `DLNM_tmin_GRID_31104_ALL.csv` — minimum-temperature grid-search results.
- `DLNM_tmean_GRID_46656_ALL.csv` — mean-temperature grid-search results.
- `DLNM_TMAX_GRID_46656_ALL.csv` — maximum-temperature grid-search results.
- `FINAL_MODEL_SELECTION_TABLE_TMIN.csv`, `FINAL_MODEL_SELECTION_TABLE_TMEAN.csv`, and `FINAL_MODEL_SELECTION_TABLE_TMAX.csv` — model-selection summaries.
- `overall_RR_curve_primary_reference.csv` — cumulative exposure–response estimates.
- `lagwise_RR_primary_reference.csv` — lag-specific relative risk estimates.
- `dlnm_synthetic_data.csv` — synthetic dataset illustrating the required input structure.

## Usage

### Required Data

Prepare a weekly dataset containing the variables required by the R scripts, including:

- Local and imported dengue cases.
- Temperature (Tmin, Tmean, or Tmax).
- Rainfall and VPD.
- Resident international arrivals.
- District identifiers, epidemiological week, and year.
- Population for the model offset.

The synthetic CSV provides an example of the input structure but cannot reproduce the study results.

### Running the Analysis

**Step 1 — Grid Search**

Run `DLNM_GLMM_Final_Grid_script.R` in R or RStudio.

Select the prepared weekly dataset when prompted. Ensure the temperature variable, maximum lag specifications, and degrees of freedom are configured for the intended grid search.

The script generates model-fitting results that can be used for subsequent model selection.

**Step 2 — Model Selection and Diagnostics**

Run `DLNM_Model_Selection_and_Diagnostics.R`.

When prompted:

1. Select the grid-search results CSV (Tmin, Tmean, or Tmax).
2. Select the original weekly analysis dataset.

The script automatically detects the temperature variable from the grid-results CSV, verifies model specifications, refits competitive models, and performs diagnostic evaluations.

**Step 3 — Final Model Estimation**

Run `DLNM_GLMM_Final_Script_Rev.R`.

Select the prepared weekly analysis dataset and set the temperature variable and model specifications corresponding to the selected model.

The script generates final relative risk estimates, model summaries, and exposure–response and lag–response figures.

## Methods Summary

- **DLNM:** Models non-linear and delayed associations between exposures and dengue incidence.
- **GLMM:** Negative-binomial regression implemented using `glmmTMB`, accounting for population and random effects.
- **Model selection:** AIC, BIC, model convergence, parameter stability, and diagnostic evaluations.
- **Model diagnostics:** Including Pearson dispersion, DHARMa residual assessments, and temporal residual checks.
- **Reference values:** Temperature, VPD, and resident arrivals are centred at their respective medians, while rainfall is centred at the local minimum of its cumulative exposure–response curve.
- **Study period:** 2014–2024, excluding the specified COVID-19 disruption period.

## Data Availability

The original dengue surveillance data were provided by the Ministry of Health and Wellness, Mauritius, and are subject to institutional and legal data-governance restrictions.

This repository does not contain restricted individual or district-level surveillance records. Instead, it provides synthetic data illustrating the required input structure, analytical code, and selected aggregate model outputs.

**Note on synthetic data:** `dlnm_synthetic_data.csv` is provided solely as a structural template. It does not represent the original study dataset. Models fitted using these synthetic records may produce different estimates or convergence behaviour, and such results should not be interpreted as reproductions of the study findings.

Access to the original surveillance data may be considered upon reasonable request and with permission from the Ministry of Health and Wellness, Mauritius, subject to applicable data-sharing agreements.

## Contact

Amar Teeluck  
The University of Queensland, Australia  
Email: m.teeluck@student.uq.edu.au
