"""
Sum carbon density (t C/ha) * pixel area (ha), and convert to Pg C.
Use north-up EPSG:4326 rasters masked to the valid vegetation domain.
Areas use a 6371 km sphere; boundary pixels are not fractionally weighted.
"""

import csv
import os

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.mask import mask


# 1. Input and output settings
input_folder = ""  # Folder containing the carbon density GeoTIFF files.
shp_path = ""      # Land boundary vector file with a defined CRS.
output_csv = ""    # Full path of the output CSV file.
file_template = "Total_126_{year}.tif"
years = range(2025, 2101, 5)  # Use range(2001, 2023) for annual historical data.
R = 6371.0  # Mean Earth radius in km.


def pixel_row_areas_ha(transform, rows):
    """Return the spherical area of one pixel in each latitude row."""
    if transform.b != 0 or transform.d != 0 or transform.a <= 0 or transform.e >= 0:
        raise ValueError("A north-up, unrotated geographic grid is required.")
    if transform.a > 360:
        raise ValueError("Pixel longitude width must not exceed 360 degrees.")
    latitude_edges = transform.f + np.arange(rows + 1) * transform.e
    if np.any(latitude_edges < -90 - 1e-8) or np.any(latitude_edges > 90 + 1e-8):
        raise ValueError("Raster latitude boundaries extend beyond the poles.")
    latitude_edges = np.radians(np.clip(latitude_edges, -90, 90))
    return R**2 * np.radians(transform.a) * np.abs(
        np.sin(latitude_edges[:-1]) - np.sin(latitude_edges[1:])
    ) * 100


def summarize_carbon(data, row_areas_ha):
    """Sum density times area over valid pixels without a full area grid."""
    values = np.asarray(data.data, dtype=np.float64)
    valid = ~np.ma.getmaskarray(data) & np.isfinite(values) & (values >= 0)
    if not np.any(valid):
        raise ValueError("No valid carbon density pixels remain inside the boundary.")
    row_totals = np.sum(np.where(valid, values, 0.0), axis=1, dtype=np.float64)
    total_carbon_pg = np.sum(row_totals * row_areas_ha, dtype=np.float64) / 1e9
    valid_area_ha = np.sum(valid.sum(axis=1) * row_areas_ha, dtype=np.float64)
    return float(total_carbon_pg), float(valid_area_ha), int(valid.sum())


def main():
    if not input_folder or not shp_path or not output_csv:
        raise ValueError("Set input_folder, shp_path and output_csv before running.")

    # 2. Load the land boundary
    land_gdf = gpd.read_file(shp_path)
    if land_gdf.crs is None:
        raise ValueError("The land boundary must have a defined CRS.")
    land_gdf = land_gdf.loc[land_gdf.geometry.notna() & ~land_gdf.geometry.is_empty]
    if land_gdf.empty or not land_gdf.geometry.is_valid.all():
        raise ValueError("Provide a nonempty land boundary with valid geometries.")
    land_gdf = land_gdf.to_crs("EPSG:4326")

    results = []
    print(f"Calculating global carbon stock for {len(years)} years...")
    for year in years:
        file_name = file_template.format(year=year)
        file_path = os.path.join(input_folder, file_name)
        record = {
            "Year": year, "Total_Carbon_PgC": "", "Valid_Area_ha": "",
            "Valid_Pixels": "", "Status": "", "Message": ""
        }
        if not os.path.isfile(file_path):
            record.update(Status="missing_file", Message=f"File not found: {file_name}")
        else:
            try:
                with rasterio.open(file_path) as src:
                    if src.crs != rasterio.crs.CRS.from_epsg(4326):
                        raise ValueError("Input rasters must use EPSG:4326.")
                    if src.count != 1:
                        raise ValueError("Input rasters must contain one carbon density band.")

                    # 3. Clip while retaining raster and boundary masks
                    data, out_transform = mask(
                        src, land_gdf.geometry, crop=True, indexes=1,
                        filled=False, all_touched=False
                    )

                    # 4. Calculate pixel areas and aggregate carbon stock
                    row_areas_ha = pixel_row_areas_ha(out_transform, data.shape[0])
                    total_pg, valid_area, valid_pixels = summarize_carbon(data, row_areas_ha)
                    record.update(
                        Total_Carbon_PgC=total_pg, Valid_Area_ha=valid_area,
                        Valid_Pixels=valid_pixels, Status="ok"
                    )
            except Exception as exc:
                record.update(Status="error", Message=str(exc))

        results.append(record)
        if record["Status"] == "ok":
            print(f"{year}: {record['Total_Carbon_PgC']:.4f} Pg C")
        else:
            print(f"{year}: {record['Status']} - {record['Message']}")

    # 5. Write every requested year, including missing or failed records
    fields = ["Year", "Total_Carbon_PgC", "Valid_Area_ha", "Valid_Pixels", "Status", "Message"]
    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    success_count = sum(row["Status"] == "ok" for row in results)
    print(f"Saved {len(results)} records ({success_count} successful) to: {output_csv}")


if __name__ == "__main__":
    main()
