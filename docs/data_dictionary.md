# Units and classification codes

## Variables

| Variable | Meaning | Typical unit |
|---|---|---|
| C | Carbon stock density | t C/ha |
| Q | Annual PM2.5 adsorption | g/m2/year |
| P | PM2.5 concentration | micrograms/m3 |
| LAI | Leaf area index | m2/m2 |
| Vd | Deposition velocity | m/s |
| T | Annual dry duration | seconds/year |

Training-target fields determine the carbon-model output units. Check whether biomass inputs represent dry biomass or carbon before integrating pools.

## Fig. 4 and Fig. 5 CSV tables

Both tables contain `Scenario`, `Indicator` and `Global_Land_Weighted_Mean`. `Historical` denotes the historical period; `SSP126`, `SSP245` and `SSP585` denote the three future pathways. `Carbon` and `PM25_Adsorption` identify the two indicators.

The Fig. 4 table contains mean levels in t C/ha and g/m2/year. The Fig. 5 table contains annual OLS rates, in t C/ha/year and (g/m2/year)/year. Multiply the Fig. 5 values by 10 to reproduce decadal bar values; do not apply this conversion to Fig. 4. The legacy column name is retained in both tables. Raster-path fields are excluded.

## Latitude profiles

The Fig. 1 workbook contains annual latitude-band sums, the across-year mean, sample standard deviation and 1.96 times that standard deviation. These are sums of raster values, not area-integrated Pg C or Mt totals. Mean ± 1.96 SD describes temporal spread; it is not a confidence interval for the mean. Any separate display multiplier must be checked against the final plotted axes.

## Eight-class temporal attribution

The hundreds digit is the dominant driver (1 = LAI, 2 = PM2.5 concentration), the tens digit is the Q–C correlation sign (1 = positive, 2 = negative), and the units digit is the dominant driver's OLS trend sign (1 = positive, 2 = negative).

| Code | Interpretation |
|---|---|
| 111 | LAI-synergy (+) |
| 112 | LAI-synergy (-) |
| 121 | LAI-tradeoff (+) |
| 122 | LAI-tradeoff (-) |
| 211 | PM2.5-synergy (+) |
| 212 | PM2.5-synergy (-) |
| 221 | PM2.5-tradeoff (+) |
| 222 | PM2.5-tradeoff (-) |

Here synergy/tradeoff denotes the sign of Q–C correlation. The parenthetical sign denotes the driver trend. Zero denotes invalid or unclassified boundary cases in the eight-class raster.

Spatial attribution instead uses 11/12 for LAI-dominated positive/negative correlation and 21/22 for PM2.5-dominated positive/negative correlation. The separate adsorption-trend driver script uses signs of contributions, with 30 for equal nonzero magnitudes and 40 for both zero; its codes must not be interpreted as correlation classes.

SOM cluster numbers are categorical labels. Interpret each cluster from its variable profile, rather than assuming numeric IDs carry a fixed ecological meaning across separate runs.
