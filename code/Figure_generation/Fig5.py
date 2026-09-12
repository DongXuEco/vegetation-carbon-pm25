# -*- coding: utf-8 -*-
"""
Plot carbon density and PM2.5 adsorption as bivariate maps.
Use geographic rasters with consistent periods and units across scenarios.
"""

import os
import glob
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
from cartopy.mpl.ticker import LatitudeFormatter
import matplotlib.ticker as mticker
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import reproject, Resampling

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LatitudeFormatter
warnings.filterwarnings("ignore")

# =========================================================
# 1. Path settings
# =========================================================
base_dir = ""  # Directory containing the input rasters.

output_png = ""  # Full output path for the PNG figure.
output_pdf = ""  # Full output path for the PDF figure.

# Match input filenames
def find_one(pattern):
    files = glob.glob(os.path.join(base_dir, pattern))
    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {pattern}")
    if len(files) > 1:
        print(f"[Warning] Multiple files match {pattern}; using the first:\n{files[0]}")
    return files[0]

scenario_files = {
    "historical": {
        "carbon": find_one("Total_HIS*.tif"),
        "pm25":   find_one("Adsorption_HIS*.tif"),
        "title":  "Historical"
    },
    "ssp126": {
        "carbon": find_one("Total_126*.tif"),
        "pm25":   find_one("Adsorption_126*.tif"),
        "title":  "SSP1-2.6"
    },
    "ssp245": {
        "carbon": find_one("Total_245*.tif"),
        "pm25":   find_one("Adsorption_245*.tif"),
        "title":  "SSP2-4.5"
    },
    "ssp585": {
        "carbon": find_one("Total_585*.tif"),
        "pm25":   find_one("Adsorption_585*.tif"),
        "title":  "SSP5-8.5"
    }
}

scenario_order = ["historical", "ssp126", "ssp245", "ssp585"]

# =========================================================
# 2. Parameter settings
# =========================================================
CARBON_MIN, CARBON_MAX = 0.0, 280.0
PM_MIN, PM_MAX = 0.0, 30.0
NCLASS = 10

LAT_MIN, LAT_MAX = -60, 83
LON_MIN, LON_MAX = -180, 180

DPI = 600

# =========================================================
# 3. Global plot style
# =========================================================
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Roboto", "Arial", "Helvetica", "DejaVu Sans"]
mpl.rcParams["font.size"] = 10
mpl.rcParams["axes.titlesize"] = 12
mpl.rcParams["axes.labelsize"] = 11
mpl.rcParams["figure.facecolor"] = "white"
mpl.rcParams["savefig.facecolor"] = "white"
mpl.rcParams["axes.unicode_minus"] = False

# Background and boundaries
OCEAN_COLOR = "#FFFFFF"
COAST_COLOR = "#AEB5BA"
BORDER_COLOR = "#D5DADD"

# =========================================================
# 4. 10×10 bivariate palette
#    Top = high PM2.5
#    Bottom = low PM2.5
#    Left = low carbon
#    Right = high carbon
# =========================================================
palette_rgb = np.array([

[(114, 47,118),(155, 66,136),(180, 83,151),(200,106,163),(216,122,175),(230,143,178),(236,130,147),(241,118,116),(246,106, 84),(252, 93, 53)],

[(139, 72,150),(174, 95,166),(199,118,180),(215,135,189),(229,158,197),(235,168,194),(239,156,161),(242,145,129),(246,134, 96),(250,123, 64)],

[(151,100,175),(181,127,189),(201,151,201),(223,169,208),(234,186,213),(240,192,210),(242,182,176),(244,173,142),(246,163,108),(248,154, 74)],

[(161,128,194),(191,159,210),(214,180,219),(229,195,224),(240,208,229),(245,216,226),(245,208,191),(245,200,156),(245,192,120),(245,184, 85)],

[(167,152,214),(200,183,225),(221,207,235),(234,216,235),(243,231,241),(250,241,242),(248,234,206),(246,228,169),(245,221,132),(243,214, 96)],

[(152,163,221),(180,190,231),(205,213,239),(224,226,243),(238,239,246),(249,246,232),(250,246,214),(245,240,180),(238,233,147),(222,212, 83)],

[(130,148,212),(159,177,225),(191,205,236),(211,219,241),(228,234,244),(244,244,225),(242,241,192),(232,233,170),(218,223,134),(190,203, 74)],

[(106,131,204),(134,160,217),(165,188,230),(186,205,236),(211,224,241),(231,237,222),(227,232,194),(216,223,162),(197,211,126),(165,191, 68)],

[( 90,113,193),(112,143,209),(139,170,225),(158,190,230),(184,207,236),(213,224,214),(204,219,180),(189,208,146),(169,197,112),(128,169, 63)],

[( 72, 91,179),( 96,124,199),(115,147,214),(128,167,222),(157,188,229),(188,208,202),(177,201,160),(163,191,128),(142,177, 97),( 88,140, 52)]

], dtype=float) / 255.0

# =========================================================
# 5. Utility functions
# =========================================================
def read_raster_crop(path, lon_min, lat_min, lon_max, lat_max):
    with rasterio.open(path) as src:
        win = from_bounds(lon_min, lat_min, lon_max, lat_max, src.transform)
        win = win.round_offsets().round_lengths()

        data = src.read(1, window=win, masked=True)
        transform = src.window_transform(win)

        arr = data.filled(np.nan).astype(np.float32)
        if src.nodata is not None:
            arr[arr == src.nodata] = np.nan
        arr[~np.isfinite(arr)] = np.nan

        meta = {
            "shape": arr.shape,
            "transform": transform,
            "crs": src.crs,
            "width": arr.shape[1],
            "height": arr.shape[0]
        }
    return arr, meta


def same_grid(meta1, meta2):
    return (
        meta1["shape"] == meta2["shape"] and
        meta1["crs"] == meta2["crs"] and
        tuple(meta1["transform"]) == tuple(meta2["transform"])
    )


def reproject_to_match(src_arr, src_meta, dst_meta, resampling=Resampling.bilinear):
    src_fill = -9999.0
    dst_fill = -9999.0

    src_tmp = np.where(np.isfinite(src_arr), src_arr, src_fill).astype(np.float32)
    dst_arr = np.full(dst_meta["shape"], dst_fill, dtype=np.float32)

    reproject(
        source=src_tmp,
        destination=dst_arr,
        src_transform=src_meta["transform"],
        src_crs=src_meta["crs"],
        src_nodata=src_fill,
        dst_transform=dst_meta["transform"],
        dst_crs=dst_meta["crs"],
        dst_nodata=dst_fill,
        resampling=resampling
    )

    dst_arr[dst_arr == dst_fill] = np.nan
    return dst_arr


def classify_fixed(arr, vmin, vmax, nclass):
    edges = np.linspace(vmin, vmax, nclass + 1)
    arr_clip = np.clip(arr, vmin, vmax)

    idx = np.digitize(arr_clip, bins=edges[1:-1], right=False)
    idx = np.clip(idx, 0, nclass - 1)

    idx[~np.isfinite(arr)] = -1
    return idx, edges


def make_rgba_from_classes(x_idx, y_idx, palette):
    """
    x_idx: carbon classes (0~9), left low -> right high
    y_idx: pm25 classes   (0~9), bottom low -> top high
    palette row 0 = top(high PM), row 9 = bottom(low PM)
    """
    h, w = x_idx.shape
    rgba = np.zeros((h, w, 4), dtype=np.float32)

    valid = (x_idx >= 0) & (y_idx >= 0)

    # Reverse the y class indices to obtain palette row indices.
    rows = (NCLASS - 1 - y_idx[valid]).astype(int)
    cols = x_idx[valid].astype(int)

    rgba[valid, :3] = palette[rows, cols]
    rgba[valid, 3] = 1.0
    return rgba


def add_bivariate_legend(fig, palette, box=(0.81, 0.18, 0.16, 0.34)):
    """
    Draw the bivariate legend on the right.
    """
    ax_leg = fig.add_axes(box)

    # Keep the palette row order:
    # Top = high PM
    # Bottom = low PM
    ax_leg.imshow(
        palette,
        origin="upper",
        extent=[0, NCLASS, 0, NCLASS],
        interpolation="nearest"
    )

    ax_leg.set_aspect("equal", adjustable="box")
    ax_leg.set_xlim(0, NCLASS)
    ax_leg.set_ylim(0, NCLASS)

    # Grid lines
    for k in range(NCLASS + 1):
        ax_leg.axhline(k, color="white", lw=0.38, alpha=0.85)
        ax_leg.axvline(k, color="white", lw=0.38, alpha=0.85)

    # Outer frame
    for spine in ax_leg.spines.values():
        spine.set_linewidth(0.7)
        spine.set_edgecolor("#AEB5BA")

    # Tick labels in data units
    ax_leg.set_xticks([0, 2.5, 5, 7.5, 10])
    ax_leg.set_xticklabels(["0", "70", "140", "210", "280"],
                           fontsize=15, color="#4A4F54")

    # With this extent and origin='upper', the top is y=10 and the bottom is y=0.
    # Labels increase from bottom to top.
    ax_leg.set_yticks([0, 2.5, 5, 7.5, 10])
    ax_leg.set_yticklabels(["0", "7.5", "15", "22.5", "30"],
                           fontsize=15, color="#4A4F54")

    ax_leg.tick_params(axis="both", which="major",
                       labelsize=15, width=0.6, length=3, colors="#4A4F54")

    ax_leg.set_xlabel("Carbon density (t/ha)", fontsize=15, labelpad=6, color="#2D3134")
    ax_leg.set_ylabel("PM$^{2.5}$ adsorption (g/m$^{2}$)", fontsize=15, labelpad=8, color="#2D3134")

    # Corner labels
    ax_leg.text(0.4, 9.6, "PM$_{2.5}$ ↑\nCarbon ↓",
                ha="left", va="top", fontsize=15, color="white")

    ax_leg.text(9.6, 9.6, "Co-high",
                ha="right", va="top", fontsize=15, color="white")

    ax_leg.text(0.4, 0.4, "Co-low",
                ha="left", va="bottom", fontsize=15, color="white")

    ax_leg.text(9.6, 0.4, "Carbon ↑\nPM$_{2.5}$ ↓",
                ha="right", va="bottom", fontsize=15, color="white")

    return ax_leg


def add_two_band_explanation(fig, palette, box=(0.78, 0.58, 0.22, 0.16)):
    """
    Two horizontal color strips:
    1) Carbon storage density: low on the left, high on the right
    2) PM2.5 adsorption: low on the left, high on the right
    Both strips have the same width.
    """
    ax_note = fig.add_axes(box)
    ax_note.set_xlim(0, 1)
    ax_note.set_ylim(0, 1)
    ax_note.axis("off")

    # Title
    ax_note.text(
        0.62, 1, "Value intensity",
        ha="center", va="top",
        fontsize=15, color="#333333"
    )

    # Shared horizontal bounds for both color strips
    band_x0, band_x1 = 0.41, 0.81
    band_h = 0.18

    # =========================
    # 1. Carbon color strip
    # =========================
    carbon_strip = palette[-1, :, :]   # Low on the left, high on the right

    y0, y1 = 0.66, 0.66 + band_h

    ax_note.imshow(
        carbon_strip[np.newaxis, :, :],
        extent=[band_x0, band_x1, y0, y1],
        origin="lower",
        aspect="auto",
        interpolation="nearest"
    )

    for k in range(1, carbon_strip.shape[0]):
        xx = band_x0 + (band_x1 - band_x0) * k / carbon_strip.shape[0]
        ax_note.plot([xx, xx], [y0, y1], color="white", lw=0.8, alpha=0.7)

    ax_note.text(
        0.01, (y0 + y1) / 2,
        "Carbon density",
        ha="left", va="center",
        fontsize=15, color="#444444"
    )

    # =========================
    # 2. PM2.5 color strip
    # =========================
    pm_strip = palette[::-1, 0, :]   # Low on the left, high on the right

    y0p, y1p = 0.30, 0.30 + band_h

    ax_note.imshow(
        pm_strip[np.newaxis, :, :],
        extent=[band_x0, band_x1, y0p, y1p],
        origin="lower",
        aspect="auto",
        interpolation="nearest"
    )

    for k in range(1, pm_strip.shape[0]):
        xx = band_x0 + (band_x1 - band_x0) * k / pm_strip.shape[0]
        ax_note.plot([xx, xx], [y0p, y1p], color="white", lw=0.8, alpha=0.7)

    ax_note.text(
        0.01, (y0p + y1p) / 2,
        "PM$_{2.5}$ adsorption",
        ha="left", va="center",
        fontsize=15, color="#444444"
    )

    # =========================
    # 3. Lower -> Higher arrow
    # =========================
    ax_note.annotate(
        "",
        xy=(0.82, 0.10), xytext=(0.41, 0.10),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", lw=1.3, color="#333333"),
        annotation_clip=False
    )

    ax_note.text(
        0.469, 0, "Lower",
        ha="center", va="top",
        fontsize=15, color="#333333",
        transform=ax_note.transAxes
    )
    ax_note.text(
        0.751, 0, "Higher",
        ha="center", va="top",
        fontsize=15, color="#333333",
        transform=ax_note.transAxes
    )

    return ax_note


# =========================================================
# 6. Check input files
# =========================================================
print("Base directory:", base_dir)
print("Directory exists:", os.path.exists(base_dir))

for scn in scenario_order:
    for key in ["carbon", "pm25"]:
        path = scenario_files[scn][key]
        print(path, "->", os.path.exists(path))
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")

# =========================================================
# 7. Read rasters and align grids
# =========================================================
ref_arr, ref_meta = read_raster_crop(
    scenario_files["historical"]["carbon"],
    LON_MIN, LAT_MIN, LON_MAX, LAT_MAX
)

data = {}

for scn in scenario_order:
    carbon_arr, carbon_meta = read_raster_crop(
        scenario_files[scn]["carbon"],
        LON_MIN, LAT_MIN, LON_MAX, LAT_MAX
    )
    pm_arr, pm_meta = read_raster_crop(
        scenario_files[scn]["pm25"],
        LON_MIN, LAT_MIN, LON_MAX, LAT_MAX
    )

    if not same_grid(carbon_meta, ref_meta):
        print(f"[Info] Resampling the {scn} carbon raster to the reference grid...")
        carbon_arr = reproject_to_match(carbon_arr, carbon_meta, ref_meta, resampling=Resampling.bilinear)
        carbon_meta = ref_meta

    if not same_grid(pm_meta, ref_meta):
        print(f"[Info] Resampling the {scn} PM2.5 raster to the reference grid...")
        pm_arr = reproject_to_match(pm_arr, pm_meta, ref_meta, resampling=Resampling.bilinear)
        pm_meta = ref_meta

    data[scn] = {
        "carbon": carbon_arr,
        "pm25": pm_arr,
        "meta": ref_meta,
        "title": scenario_files[scn]["title"]
    }

# =========================================================
# 8. Classify values and generate RGBA maps
# =========================================================
for scn in scenario_order:
    carbon_idx, _ = classify_fixed(data[scn]["carbon"], CARBON_MIN, CARBON_MAX, NCLASS)
    pm_idx, _ = classify_fixed(data[scn]["pm25"], PM_MIN, PM_MAX, NCLASS)

    rgba = make_rgba_from_classes(carbon_idx, pm_idx, palette_rgb)
    data[scn]["rgba"] = rgba

# =========================================================
# 9. Plot maps using a manual panel layout
# =========================================================
proj = ccrs.Robinson(central_longitude=0)
pc = ccrs.PlateCarree()

fig = plt.figure(figsize=(20.5, 8.4), facecolor="white")

# Four maps on the left and legends on the right
panel_pos = {
    "historical": [0.02, 0.55, 0.36, 0.38],
    "ssp126":     [0.40, 0.55, 0.36, 0.38],
    "ssp245":     [0.02, 0.08, 0.36, 0.38],
    "ssp585":     [0.40, 0.08, 0.36, 0.38],
}

panel_labels = ["A", "B", "C", "D"]

for i, scn in enumerate(scenario_order):
    ax = fig.add_axes(panel_pos[scn], projection=proj)

    ax.set_facecolor("white")
    ax.add_feature(cfeature.OCEAN.with_scale("110m"), facecolor=OCEAN_COLOR, zorder=0)

    rgba = data[scn]["rgba"]
    h, w = rgba.shape[:2]
    transform = data[scn]["meta"]["transform"]

    left = transform.c
    top = transform.f
    right = left + transform.a * w
    bottom = top + transform.e * h
    extent = [left, right, bottom, top]

    ax.imshow(
        rgba,
        origin="upper",
        extent=extent,
        transform=pc,
        interpolation="nearest",
        zorder=2
    )

    ax.add_feature(
        cfeature.COASTLINE.with_scale("110m"),
        linewidth=0.28,
        edgecolor=COAST_COLOR,
        zorder=3
    )
    ax.add_feature(
        cfeature.BORDERS.with_scale("110m"),
        linewidth=0.10,
        edgecolor=BORDER_COLOR,
        zorder=3
    )

    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=pc)

    # Add latitude labels to panels A and C only.
    # Show latitude labels on the left of panels A and C.
    if scn in ["historical", "ssp245"]:
        gl = ax.gridlines(
            crs=pc,
            draw_labels=True,
            linewidth=0.0,
            color="none",
            alpha=0
        )

        gl.top_labels = False
        gl.bottom_labels = False
        gl.right_labels = False
        gl.left_labels = True

        gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60])
        gl.yformatter = LatitudeFormatter()

        gl.ylabel_style = {
            "size": 15,
            "color": "#4A4F54"
        }

    else:
        # Omit latitude labels from panels B and D.
        gl = ax.gridlines(
            crs=pc,
            draw_labels=False,
            linewidth=0.0,
            color="none",
            alpha=0
        )

    # Set the GeoAxes frame when the geo spine is available.
    if "geo" in ax.spines:
        ax.spines["geo"].set_linewidth(1.15)
        ax.spines["geo"].set_edgecolor("black")

    # Panel labels A/B/C/D in the upper-left corner
    ax.text(
        0.018, 0.985, panel_labels[i],
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=22, fontweight="bold", color="#222222"
    )

    # Scenario titles in the lower-left part of each map
    ax.text(
        0.10, 0.12, data[scn]["title"],
        transform=ax.transAxes,
        ha="left", va="bottom",
        fontsize=15, color="#111111", fontweight="normal"
    )

# Legends on the right
add_two_band_explanation(fig, palette_rgb, box=(0.78, 0.58, 0.22, 0.16))
add_bivariate_legend(fig, palette_rgb, box=(0.81, 0.17, 0.16, 0.34))

# Save figures
fig.savefig(output_png, dpi=DPI, bbox_inches="tight", facecolor="white")
fig.savefig(output_pdf, dpi=DPI, bbox_inches="tight", facecolor="white")
plt.show()

print("\nPlotting complete.")
print("PNG:", output_png)
print("PDF:", output_pdf)
