# Workflow

## Input preparation

Prepare predictors and training samples following the study methods and data sources. Use matching predictor band names, physical units, soil depth and target definitions for training and prediction. Configure a common valid vegetation mask where requested. Public Earth Engine imports remain in the relevant scripts.

Future predictors are prepared for SSP1-2.6, SSP2-4.5 and SSP5-8.5 at five-year intervals. The ensemble script averages five models after the temporal aggregation and grid alignment have been prepared.

## Carbon and adsorption

Use the biomass scripts for historical or future prediction. For SOC and SIC, fit each pool separately using the prepared 2001–2022 mean predictors and the appropriate sample field. Carbon model outputs retain the sample-target units; do not apply a second biomass-to-carbon conversion to an input already expressed as carbon.

PM2.5 adsorption is calculated as `Q = P * Vd * LAI * T / 1e6`, with P in micrograms/m3, Vd in m/s, LAI in m2/m2 and T in seconds/year. Q is in g/m2/year. Historical dry duration uses precipitation below 0.1 mm/day. Future Vd and T are fixed to 2022 while P and LAI vary.

## Trends, totals and attribution

Use OLS for the supplied trend workflows. The MATLAB script defaults to a ten-year scaling; the temporal attribution script calculates slopes using the supplied year coordinates. Check units before combining outputs with other scripts or source tables.

Carbon totals use density multiplied by pixel area and are converted to Pg C. Spatial correlations and attribution use paired valid data. Temporal classification uses the dominant driver, the sign of Q–C correlation and the sign of the dominant driver's trend.

## Figures

Use the source-data index to locate tables and vectors. The two included Python plotting scripts require external rasters. Confirm that input trend scaling matches the map class breaks and legend: annual and decadal slopes differ by a factor of ten.

This upload package preserves the current source files. Packaging checks verify file integrity and script syntax; they do not rerun the Earth Engine models or establish full numerical reproduction of every manuscript panel.
