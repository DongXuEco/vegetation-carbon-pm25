# Global vegetation carbon stocks and PM2.5 interception decouple across future pollution pathways

This repository contains the supporting materials, source data and analysis code for our manuscript, **Global vegetation carbon stocks and PM2.5 interception decouple across future pollution pathways**, prepared for submission to *Nature Ecology & Evolution*.

The study examines vegetation carbon stocks and PM2.5 interception under historical conditions (2001–2022) and three future pathways: SSP1-2.6, SSP2-4.5 and SSP5-8.5 (2025–2100 at five-year intervals).

The repository includes carbon-stock modelling, PM2.5 adsorption calculations, temporal and spatial analyses, and selected map-generation scripts. Figure source data are provided as spreadsheets, CSV tables, vector datasets and temporal and spatial attribution rasters.

## Repository structure

```text
vegetation-carbon-pm25/
├── README.md
├── code/
│   ├── Data_generation/               # Earth Engine models and R clustering
│   ├── Data_processing_and_analysis/  # Trends, correlations and attribution
│   └── Figure_generation/             # Python bivariate maps
├── Data/                             # Figure source tables, vectors and rasters
├── docs/
│   ├── workflow.md
│   ├── data_dictionary.md
│   └── reuse.md
├── requirements.txt
└── manifest.csv                      # Source-file sizes and SHA-256 hashes
```

## Getting started

1. Read the [code index](code/README.md) and [workflow](docs/workflow.md).
2. Prepare the inputs described in the scripts and in the study methods. Private Earth Engine asset IDs and local input/output paths are intentionally blank. Public Earth Engine dataset imports are retained.
3. Set the target fields, predictor bands, years, scenarios and export settings for the selected analysis. Some defaults cover only example years and do not execute the complete study period.
4. Run the relevant script in Google Earth Engine, Python, R or MATLAB.

The source tables can be inspected independently of the modelling inputs. This repository does not contain the complete set of processed predictor rasters, training samples or private Earth Engine assets; those inputs must be prepared before running the modelling workflows.

## Software

- **Google Earth Engine:** run `.js` files in the Code Editor with an enabled account and access to the configured assets.
- **Python:** NumPy, Rasterio, GeoPandas, Matplotlib and Cartopy. Install with `python -m pip install -r requirements.txt` in a suitable environment.
- **R:** dplyr, readxl and kohonen (`install.packages(c("dplyr", "readxl", "kohonen"))`).
- **MATLAB:** Mapping Toolbox is required for the GeoTIFF workflow.

Python dependencies are listed without version pins; this package is not a locked computational environment. Cartopy may download Natural Earth resources when first plotting.

## Source data

See the [data index](Data/README.md) for figure-to-file mappings and the [data dictionary](docs/data_dictionary.md) for units and classification codes. Keep all companion files together when using a shapefile.

The plotting filenames follow an earlier numbering: `Fig4.py` plots trends corresponding to the Fig. 5 trend source table, while `Fig5.py` plots levels corresponding to the Fig. 4 mean source table.

## Reuse

Please cite the associated study when reusing its methods or data. Publication metadata and a reuse licence have not yet been supplied for this repository; see [reuse information](docs/reuse.md).
