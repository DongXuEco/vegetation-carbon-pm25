# Spatial drivers of carbon–PM2.5 adsorption correlation

Categorical rasters supplied for the spatial attribution of the correlation between carbon stock (C) and PM2.5 adsorption (Q).

| Raster | Scenario |
|---|---|
| HIS_Spatial_Attribution_9x9.tif | Historical |
| 126_Spatial_Attribution_9x9.tif | SSP1-2.6 |
| 245_Spatial_Attribution_9x9.tif | SSP2-4.5 |
| 585_Spatial_Attribution_9x9.tif | SSP5-8.5 |

The filenames identify a 9 × 9 pixel attribution window. The separate local correlation script uses a different window; do not infer its settings from these filenames.

The four-class convention used by the spatial attribution code is:

| Code | Dominant driver | Q–C correlation |
|---|---|---|
| 11 | LAI | Positive |
| 12 | LAI | Negative |
| 21 | PM2.5 concentration | Positive |
| 22 | PM2.5 concentration | Negative |

These signs refer to spatial correlation, not temporal trends. See the [data dictionary](../../docs/data_dictionary.md).

Original rasters, attribute tables, encoding files and overviews are retained. Download the entire directory and use nearest-neighbour resampling for categorical values. Temporary GIS lock files are excluded.

The supplied outputs are copied without reclassification or recalculation. Checksums establish copy integrity; model settings and manuscript statistics have not been independently reconstructed from these files.
