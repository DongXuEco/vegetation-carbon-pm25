# Temporal drivers of PM2.5 adsorption

Categorical rasters supplied for the temporal driver attribution of vegetation PM2.5 adsorption (Q).

| Raster | Scenario |
|---|---|
| Attribution_HIS.tif | Historical |
| Attribution_SSP126.tif | SSP1-2.6 |
| Attribution_SSP245.tif | SSP2-4.5 |
| Attribution_SSP585.tif | SSP5-8.5 |

The companion raster attribute tables contain classes 11, 12, 21 and 22. The adsorption temporal driver code uses the following convention:

| Code | Dominant driver | Contribution sign |
|---|---|---|
| 11 | LAI | Positive |
| 12 | LAI | Negative |
| 21 | PM2.5 concentration | Positive |
| 22 | PM2.5 concentration | Negative |

These codes describe contributions to Q trends, not Q–C correlation signs. They are distinct from the temporal and spatial correlation attribution datasets. See the [data dictionary](../../docs/data_dictionary.md).

Download the entire directory to retain raster attribute tables, encoding files, overviews and available auxiliary metadata. Attribute-table Count fields are pixel counts, not area-weighted percentages. Use nearest-neighbour resampling for categorical data.

Regional area proportions must be calculated from the rasters using pixel-area weights.
