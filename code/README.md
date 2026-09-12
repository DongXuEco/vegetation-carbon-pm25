# Code index

Supporting code for **Global vegetation carbon stocks and PM2.5 interception decouple across future pollution pathways**.

Private assets and local paths must be configured before execution. Output units follow the input target where indicated.

| Directory | Script | Purpose |
|---|---|---|
| Data_generation | Biomass_carbon_stock.js | Fit stacked biomass models, validate and predict annual inputs |
| Data_generation | Biomass_carbon_SSP_prediction.js | Train biomass models and predict future SSP inputs |
| Data_generation | Soil_carbon_stock.js | Fit SOC or SIC using 2001–2022 mean predictors, then predict annual inputs |
| Data_generation | PM25_air_purification.js | Historical PM2.5 adsorption using annual Vd bands and dry-day duration |
| Data_generation | PM25_SSP_adsorption.js | SSP adsorption using scenario PM2.5 and LAI with 2022 Vd and dry duration |
| Data_generation | Code for cluster analysis (SOM).R | Standardize four variables, train 16 SOM units and group them into four clusters |
| Data_processing_and_analysis | Multi_model_ensemble_mean.py | Equal-weight mean over five models with common valid pixels |
| Data_processing_and_analysis | Global_carbon_stock_total.py | Area-integrated carbon stock in Pg C |
| Data_processing_and_analysis | OLS_trend_analysis.m | Pixel-wise OLS slopes; configurable annual or decadal units |
| Data_processing_and_analysis | Carbon_PM25_spatial_correlation.js | Local Pearson correlation and P-values in 19 × 19 windows |
| Data_processing_and_analysis | Carbon_adsorption_spatial_attribution.js | Four-class attribution of spatial correlation in 9 × 9 windows |
| Data_processing_and_analysis | Carbon_adsorption_temporal_attribution.py | Eight-class attribution of temporal correlation |
| Data_processing_and_analysis | Adsorption_temporal_driver_attribution.py | Relative LAI and PM2.5 contributions to adsorption trends |
| Figure_generation | Fig4.py | Bivariate trend maps; corresponds to the Fig. 5 trend data |
| Figure_generation | Fig5.py | Bivariate mean-value maps; corresponds to the Fig. 4 mean data |

The two spatial scripts perform different analyses and intentionally have different window sizes. Both currently use a nominal 10 km grid; ground distances depend on projection and latitude.
