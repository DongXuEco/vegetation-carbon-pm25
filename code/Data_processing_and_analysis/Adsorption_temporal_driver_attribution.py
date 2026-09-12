"""
Compare LAI trend * mean(PM2.5) with PM2.5 trend * mean(LAI).
Use OLS trend rasters with the same period, units and valid vegetation domain; P is concentration.
11/12: LAI dominant, positive/negative contribution; 21/22: PM2.5 dominant, positive/negative.
30: equal nonzero magnitudes; 40: both zero; 0: NoData.
First-order approximation with positive, time-constant Vd*T; signs refer to contributions.
"""

import os

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling


# 1. Input and output settings
base_dir = ""       # Folder containing trend and period-mean rasters.
output_folder = ""  # Folder for the attribution GeoTIFF files.
scenarios = ["126", "245", "585"]
lai_trend_template = "cj_OLS_LAI_{ssp}.tif"
pm25_trend_template = "cj_OLS_pm25_{ssp}.tif"
lai_mean_template = "lai_ssp{ssp}.tif"
pm25_mean_template = "mmrpm2p5_ssp{ssp}.tif"


def read_values(src):
    """Read physical values, preserving the source mask and NoData."""
    if src.count != 1 or src.crs is None:
        raise ValueError(f"A single-band raster with a defined CRS is required: {src.name}")
    if src.scales[0] != 1 or src.offsets[0] != 0:
        raise ValueError(f"Apply scale and offset before using this raster: {src.name}")
    data = src.read(1, masked=True).astype(np.float64).filled(np.nan)
    data[~np.isfinite(data)] = np.nan
    return data


def read_and_align(target_path, reference, nonnegative=False):
    """Align continuous values to the reference grid using bilinear resampling."""
    with rasterio.open(target_path) as src:
        data = read_values(src)
        if nonnegative:
            data[data < 0] = np.nan
        shape = (reference["height"], reference["width"])
        if (src.crs == reference["crs"] and src.transform == reference["transform"]
                and data.shape == shape):
            return data
        destination = np.full(shape, np.nan, dtype=np.float64)
        reproject(
            source=data,
            destination=destination,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=np.nan,
            dst_transform=reference["transform"],
            dst_crs=reference["crs"],
            dst_nodata=np.nan,
            init_dest_nodata=True,
            resampling=Resampling.bilinear
        )
        return destination


def classify_drivers(s_lai, s_pm25, m_lai, m_pm25):
    """Classify contributions on the common valid domain after alignment."""
    if not (s_lai.shape == s_pm25.shape == m_lai.shape == m_pm25.shape):
        raise ValueError("All input arrays must have the same shape.")
    valid = (np.isfinite(s_lai) & np.isfinite(s_pm25)
             & np.isfinite(m_lai) & np.isfinite(m_pm25)
             & (m_lai >= 0) & (m_pm25 >= 0))
    with np.errstate(over="ignore", invalid="ignore"):
        con_lai = s_lai * m_pm25
        con_pm25 = s_pm25 * m_lai
    valid &= np.isfinite(con_lai) & np.isfinite(con_pm25)
    if not np.any(valid):
        raise ValueError("No common valid pixels remain after alignment.")

    attr_map = np.zeros(s_lai.shape, dtype=np.int16)
    lai_dominant = np.abs(con_lai) > np.abs(con_pm25)
    pm25_dominant = np.abs(con_pm25) > np.abs(con_lai)
    both_zero = (con_lai == 0) & (con_pm25 == 0)
    equal_magnitude = np.abs(con_lai) == np.abs(con_pm25)

    attr_map[valid & lai_dominant & (con_lai > 0)] = 11
    attr_map[valid & lai_dominant & (con_lai < 0)] = 12
    attr_map[valid & pm25_dominant & (con_pm25 > 0)] = 21
    attr_map[valid & pm25_dominant & (con_pm25 < 0)] = 22
    attr_map[valid & equal_magnitude & ~both_zero] = 30
    attr_map[valid & both_zero] = 40
    return attr_map


def align_and_calculate():
    if not base_dir or not output_folder:
        raise ValueError("Set base_dir and output_folder before running the script.")
    os.makedirs(output_folder, exist_ok=True)

    for ssp in scenarios:
        print(f"Processing SSP{ssp}...")
        path_s_lai = os.path.join(base_dir, lai_trend_template.format(ssp=ssp))
        path_s_pm = os.path.join(base_dir, pm25_trend_template.format(ssp=ssp))
        path_m_lai = os.path.join(base_dir, lai_mean_template.format(ssp=ssp))
        path_m_pm = os.path.join(base_dir, pm25_mean_template.format(ssp=ssp))

        # 2. Read the LAI trend reference and align the other inputs
        with rasterio.open(path_s_lai) as ref:
            s_lai = read_values(ref)
            ref_meta = ref.meta.copy()
        s_pm25 = read_and_align(path_s_pm, ref_meta)
        m_lai = read_and_align(path_m_lai, ref_meta, nonnegative=True)
        m_pm25 = read_and_align(path_m_pm, ref_meta, nonnegative=True)

        # 3. Calculate contributions and identify the dominant driver
        attr_map = classify_drivers(s_lai, s_pm25, m_lai, m_pm25)

        # 4. Export the categorical attribution map
        ref_meta.update(driver="GTiff", dtype="int16", count=1, nodata=0, compress="lzw")
        output_path = os.path.join(output_folder, f"Attribution_SSP{ssp}.tif")
        with rasterio.open(output_path, "w", **ref_meta) as dst:
            dst.write(attr_map, 1)
            dst.update_tags(
                CLASS_11="LAI dominant, positive", CLASS_12="LAI dominant, negative",
                CLASS_21="PM2.5 dominant, positive", CLASS_22="PM2.5 dominant, negative",
                CLASS_30="Equal nonzero absolute contributions",
                CLASS_40="Both contributions zero"
            )
        print(f"Saved: {output_path}")


if __name__ == "__main__":
    align_and_calculate()
