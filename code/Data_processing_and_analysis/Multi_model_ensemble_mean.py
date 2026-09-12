"""
Calculate equal-weight means of five models for each year.
Inputs must share the grid, units and averaging period; prepare time means first.
A pixel is retained only when all five models have valid values.
"""

import os

import numpy as np
import rasterio


# 1. Input and output settings
input_dir = ""   # Folder containing the model GeoTIFF files.
output_dir = ""  # Folder for the ensemble mean GeoTIFF files.
variable = "lai"
table = "Lmon"
scenario = "ssp126"  # ssp126, ssp245 or ssp585.
file_template = "{variable}_{table}_{model}_{scenario}_{year}.tif"
models = [
    "CESM2-WACCM", "CNRM-ESM2-1", "GFDL-ESM4",
    "MPI-ESM1-2-HR", "UKESM1-0-LL"
]
years = range(2025, 2101, 5)


def strict_ensemble_mean(arrays, model_count):
    """Average equally weighted models over their common valid pixels."""
    total = None
    common_valid = None
    count = 0
    for array in arrays:
        values = np.asarray(np.ma.getdata(array), dtype=np.float64)
        valid = ~np.ma.getmaskarray(array) & np.isfinite(values)
        if total is None:
            total = np.zeros(values.shape, dtype=np.float64)
            common_valid = np.ones(values.shape, dtype=bool)
        if values.shape != total.shape:
            raise ValueError("All model arrays must have the same shape.")
        total += np.where(valid, values, 0.0)
        common_valid &= valid
        count += 1
    if count != model_count or count == 0:
        raise ValueError("Exactly one array per model is required.")
    if not np.any(common_valid):
        raise ValueError("No pixels are valid in all models.")
    mean = np.full(total.shape, np.nan, dtype=np.float32)
    values = total[common_valid] / count
    if np.any(~np.isfinite(values)) or np.any(np.abs(values) > np.finfo(np.float32).max):
        raise ValueError("Ensemble values exceed the finite Float32 range.")
    mean[common_valid] = values
    return mean


def read_models(file_list, reference):
    """Read each model using its own NoData value and raster mask."""
    for filename in file_list:
        with rasterio.open(filename) as src:
            if src.count != 1:
                raise ValueError(f"A single-band raster is required: {filename}")
            if src.crs is None or (
                src.crs != reference["crs"]
                or src.transform != reference["transform"]
                or src.width != reference["width"]
                or src.height != reference["height"]
            ):
                raise ValueError(f"Model grids or coordinate systems differ: {filename}")
            if src.scales[0] != 1 or src.offsets[0] != 0:
                raise ValueError(f"Convert scaled values to physical units first: {filename}")
            yield src.read(1, masked=True)


def main():
    if not input_dir or not output_dir:
        raise ValueError("Set input_dir and output_dir before running the script.")
    if len(models) != 5 or len(set(models)) != 5:
        raise ValueError("Specify five distinct model names.")

    # 2. Match exactly one filename per model and year
    files_by_year = {}
    for year in years:
        file_list = [
            os.path.join(input_dir, file_template.format(
                variable=variable, table=table, model=model,
                scenario=scenario, year=year
            )) for model in models
        ]
        if len(set(file_list)) != len(models):
            raise ValueError("The filename template must distinguish every model.")
        for filename in file_list:
            if not os.path.isfile(filename):
                raise FileNotFoundError(f"Required model file not found: {filename}")
        files_by_year[year] = file_list

    os.makedirs(output_dir, exist_ok=True)
    print(f"Calculating five-model ensemble means for {scenario}...")
    for year, file_list in files_by_year.items():
        # 3. Read reference metadata and calculate the ensemble mean
        with rasterio.open(file_list[0]) as src:
            meta = src.meta.copy()
        mean_array = strict_ensemble_mean(read_models(file_list, meta), len(models))

        # 4. Export using the shared spatial grid
        meta.update(driver="GTiff", dtype="float32", count=1, nodata=np.nan, compress="lzw")
        output_filename = f"{variable}_{table}_EnsembleMean_{scenario}_{year}.tif"
        output_path = os.path.join(output_dir, output_filename)
        with rasterio.open(output_path, "w", **meta) as dst:
            dst.write(mean_array, 1)
        print(f"{year}: saved {output_filename}")

    print("All requested ensemble means have been saved.")


if __name__ == "__main__":
    main()
