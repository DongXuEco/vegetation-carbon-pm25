"""
Identify drivers of Q-C temporal correlation from complete, matching yearly rasters.
Use physical values and PM2.5 concentration, not mass mixing ratio.
Class digits: driver (1 LAI, 2 PM2.5), R sign, dominant-driver OLS trend sign.
Signs: 1 positive, 2 negative; 0 marks invalid or unclassified boundary cases.
First-order covariance attribution; assumes positive, time-constant Vd*T and omits interaction.
"""

import os

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling


# 1. Input and output settings
data_dirs = {"carbon": "", "absorption": "", "lai": "", "pm25": ""}
# Enter one filename template per variable, including {year}.
# Templates may also contain {scenario}; e.g. Total_{scenario}_{year}.tif.
file_templates = {"carbon": "", "absorption": "", "lai": "", "pm25": ""}
scenario = "126"
years = np.arange(2025, 2101, 5)
output_file = ""    # Full path of the eight-class GeoTIFF.
r_output_file = ""  # Full path of the Pearson correlation GeoTIFF.
pixel_chunk_size = 32768  # Number of pixels processed per statistical block.


def load_tif_stack_aligned(files, ref_meta=None):
    """Load an explicitly ordered series and align every raster to one grid."""
    stack = None
    for index, filename in enumerate(files):
        with rasterio.open(filename) as src:
            if src.count != 1 or src.crs is None:
                raise ValueError(f"A single-band raster with a CRS is required: {filename}")
            if src.scales[0] != 1 or src.offsets[0] != 0:
                raise ValueError(f"Apply scale and offset before analysis: {filename}")
            data = src.read(1, masked=True).astype(np.float32).filled(np.nan)
            data[~np.isfinite(data) | (data < 0)] = np.nan
            if ref_meta is None:
                ref_meta = src.meta.copy()
            shape = (ref_meta["height"], ref_meta["width"])
            if stack is None:
                stack = np.empty((len(files), *shape), dtype=np.float32)
            if (data.shape != shape or src.crs != ref_meta["crs"]
                    or src.transform != ref_meta["transform"]):
                destination = np.full(shape, np.nan, dtype=np.float32)
                reproject(
                    source=data, destination=destination,
                    src_transform=src.transform, src_crs=src.crs, src_nodata=np.nan,
                    dst_transform=ref_meta["transform"], dst_crs=ref_meta["crs"],
                    dst_nodata=np.nan, init_dest_nodata=True,
                    resampling=Resampling.bilinear
                )
                data = destination
            stack[index] = data
    return stack, ref_meta


def ols_slope(values, time_years):
    """Return annual OLS slopes with an intercept for complete time-by-pixel arrays."""
    time_years = np.asarray(time_years, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    centered_years = time_years - time_years.mean()
    centered_values = values - values.mean(axis=0)
    return np.sum(centered_years[:, None] * centered_values, axis=0) / np.sum(centered_years**2)


def classify_eight(con_lai, con_pm25, r, trend_lai, trend_pm25):
    """Append trend class to the original two-digit attribution code."""
    result = np.zeros(r.shape, dtype=np.int16)
    valid = (np.isfinite(con_lai) & np.isfinite(con_pm25) & np.isfinite(r)
             & np.isfinite(trend_lai) & np.isfinite(trend_pm25))
    for driver, dominant, trend in [
        (1, np.abs(con_lai) > np.abs(con_pm25), trend_lai),
        (2, np.abs(con_pm25) > np.abs(con_lai), trend_pm25)
    ]:
        for r_code, r_mask in [(1, r > 0), (2, r < 0)]:
            original_class = driver * 10 + r_code
            for trend_code, trend_mask in [(1, trend > 0), (2, trend < 0)]:
                result[valid & dominant & r_mask & trend_mask] = original_class * 10 + trend_code
    return result


def calculate_attribution(c_stack, q_stack, l_stack, p_stack, time_years, chunk_size=32768):
    """Calculate covariance, correlation and OLS trends on identical samples."""
    stacks = [c_stack, q_stack, l_stack, p_stack]
    if any(stack.shape != c_stack.shape for stack in stacks):
        raise ValueError("All four time-series stacks must have identical shapes.")
    time_years = np.asarray(time_years, dtype=np.float64)
    if (time_years.ndim != 1 or len(time_years) < 3
            or len(time_years) != c_stack.shape[0]
            or not np.all(np.isfinite(time_years)) or np.any(np.diff(time_years) <= 0)):
        raise ValueError("Provide at least three matching, strictly increasing years.")
    if not isinstance(chunk_size, int) or chunk_size < 1:
        raise ValueError("chunk_size must be a positive integer.")
    shape = c_stack.shape[1:]
    flat = [stack.reshape(len(time_years), -1) for stack in stacks]
    attr = np.zeros(flat[0].shape[1], dtype=np.int16)
    r_map = np.full(attr.shape, np.nan, dtype=np.float32)
    for start in range(0, attr.size, chunk_size):
        stop = min(start + chunk_size, attr.size)
        blocks = [array[:, start:stop] for array in flat]
        valid = np.logical_and.reduce([
            np.all(np.isfinite(block) & (block >= 0), axis=0) for block in blocks
        ])
        if not np.any(valid):
            continue
        c, q, l, p = [block[:, valid].astype(np.float64) for block in blocks]
        mean_l, mean_p = l.mean(axis=0), p.mean(axis=0)
        dc, dq = c - c.mean(axis=0), q - q.mean(axis=0)
        cov_qc = np.mean(dq * dc, axis=0)
        con_lai = mean_p * np.mean((l - mean_l) * dc, axis=0)
        con_pm25 = mean_l * np.mean((p - mean_p) * dc, axis=0)
        std_c = np.sqrt(np.mean(dc**2, axis=0))
        std_q = np.sqrt(np.mean(dq**2, axis=0))
        r = np.full(std_c.shape, np.nan)
        defined = (std_c > 0) & (std_q > 0)
        np.divide(cov_qc, std_c * std_q, out=r, where=defined)
        r = np.clip(r, -1, 1)
        trend_lai = ols_slope(l, time_years)
        trend_pm25 = ols_slope(p, time_years)
        indices = np.flatnonzero(valid) + start
        attr[indices] = classify_eight(con_lai, con_pm25, r, trend_lai, trend_pm25)
        r_map[indices] = r
    return attr.reshape(shape), r_map.reshape(shape)


def main():
    if (not all(data_dirs.values()) or not all(file_templates.values())
            or not output_file or not r_output_file):
        raise ValueError("Set all input folders, filename templates and output paths.")
    if os.path.normcase(os.path.abspath(output_file)) == os.path.normcase(os.path.abspath(r_output_file)):
        raise ValueError("Classification and correlation outputs must use different paths.")
    files_by_variable = {}
    for variable, folder in data_dirs.items():
        files = [os.path.join(folder, file_templates[variable].format(
            year=int(year), scenario=scenario)) for year in years]
        if len(set(files)) != len(years):
            raise ValueError(f"The filename template must distinguish every year: {variable}")
        for filename in files:
            if not os.path.isfile(filename):
                raise FileNotFoundError(filename)
        files_by_variable[variable] = files

    # 2. Load all variables on the carbon reference grid
    print("Reading and aligning the four time series...")
    c_stack, meta = load_tif_stack_aligned(files_by_variable["carbon"])
    q_stack, _ = load_tif_stack_aligned(files_by_variable["absorption"], meta)
    l_stack, _ = load_tif_stack_aligned(files_by_variable["lai"], meta)
    p_stack, _ = load_tif_stack_aligned(files_by_variable["pm25"], meta)

    # 3. Calculate the eight-class attribution and Pearson correlation
    attr_map, r_map = calculate_attribution(
        c_stack, q_stack, l_stack, p_stack, years, pixel_chunk_size
    )
    print(f"Classified pixels: {np.count_nonzero(attr_map)}")
    print(f"Pixels with defined correlation: {np.count_nonzero(np.isfinite(r_map))}")

    # 4. Export both rasters with the reference spatial metadata
    for filename, data, dtype, nodata in [
        (output_file, attr_map, "int16", 0),
        (r_output_file, r_map, "float32", np.nan)
    ]:
        os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
        profile = meta.copy()
        profile.update(driver="GTiff", dtype=dtype, count=1, nodata=nodata, compress="lzw")
        with rasterio.open(filename, "w", **profile) as dst:
            dst.write(data, 1)
            if dtype == "int16":
                dst.update_tags(
                    class_encoding="driver * 100 + correlation_sign * 10 + trend_sign",
                    driver="1=LAI; 2=PM2.5 concentration",
                    signs="1=positive; 2=negative",
                    zero="NoData or unclassified boundary case",
                    trend_method="OLS slope with intercept using actual years"
                )
        print(f"Saved: {filename}")


if __name__ == "__main__":
    main()
