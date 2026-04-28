# Modelling the Effect of Climatic Drivers and Cross-Border Travel on Sporadic Outbreaks of Dengue in Mauritius

## Overview  
This repository contains the analytical scripts used to examine the association between climatic drivers, cross-border travel, and dengue incidence in Mauritius.

The study applies a **Distributed Lag Non-linear Model (DLNM)** within a **Generalized Linear Mixed Model (GLMM)** framework to capture delayed and non-linear relationships between environmental exposures and dengue risk in a **non-endemic setting**.

---

## Repository Contents  

- **Google Earth Engine (GEE) Script**  
  - Extraction and aggregation of climate variables (temperature, rainfall, dew temperature)

- **DLNM–GLMM Script (R)**  
  - Data preprocessing and lag creation  
  - Crossbasis construction for DLNM  
  - Model fitting using `glmmTMB`  
  - Model evaluation (QAIC, QBIC, overdispersion, autocorrelation, R²)  
  - Relative risk (RR) estimation using `crosspred()`  
  - Generation of overall RR and lag-response plots  

---

## Data Availability  

The data used in this study are derived from routine dengue surveillance records provided by the Ministry of Health and Wellness, Mauritius, and are subject to institutional and legal data governance restrictions.

- Raw data are **not publicly available**  
- This repository provides **code only** for transparency and reproducibility  

---

## Usage  

To run the analysis:

1. Prepare a dataset with the following variables:
   - Dengue cases (outcome)  
   - Climate variables (temperature, rainfall, humidity/VPD)  
   - Cross-border travel indicators (e.g., imported cases, resident arrivals)  
   - Time variables (week, year)  
   - Population (for offset)
   - dummy_data_csv provided as template. (Local and imported cases data in spreadsheet are generated randomly. Data are not publicly available due to institutional and legal data-governance restrictions)

2. Run:
   - GEE script (for climate data extraction)  
   - DLNM–GLMM R script (for modelling and analysis)

---

## Methods Summary  

- DLNM used to model **non-linear and lagged effects**  
- GLMM fitted using `glmmTMB` with random effects  
- Model selection based on QAIC, QBIC, and diagnostic metrics  
- Relative risks estimated using **minimum-risk value (MRV) centering**

---

## Contact  

Amar Teeluck  
Email: m.teeluck@student.uq.edu.au 
