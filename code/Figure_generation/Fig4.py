# -*- coding: utf-8 -*-
"""
Plot carbon and PM2.5 adsorption trends as bivariate maps.
Use geographic rasters; input trend units must match the class ranges and legend.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle

import rasterio
from rasterio.warp import reproject, Resampling
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# =========================================================
# 1. Path settings
# =========================================================

base_dir = ""  # Directory containing the input rasters.

scenario_files = {
    "historical": {
        "carbon": os.path.join(base_dir, "cj_OLS_carbonHIS.tif"),
        "pm25":   os.path.join(base_dir, "cj_OLS_PMHIS.tif"),
        "title":  "Historical"
    },
    "ssp126": {
        "carbon": os.path.join(base_dir, "cj_OLS_carbon126.tif"),
        "pm25":   os.path.join(base_dir, "cj_OLS_PM126.tif"),
        "title":  "SSP1-2.6"
    },
    "ssp245": {
        "carbon": os.path.join(base_dir, "cj_OLS_carbon245.tif"),
        "pm25":   os.path.join(base_dir, "cj_OLS_PM245.tif"),
        "title":  "SSP2-4.5"
    },
    "ssp585": {
        "carbon": os.path.join(base_dir, "cj_OLS_carbon585.tif"),
        "pm25":   os.path.join(base_dir, "cj_OLS_PM585.tif"),
        "title":  "SSP5-8.5"
    }
}

output_png = ""  # Full output path for the PNG figure.
output_pdf = ""  # Full output path for the PDF figure.

scenario_order = ["historical", "ssp126", "ssp245", "ssp585"]

# =========================================================
# 2. Fixed classification parameters
# =========================================================

# Fixed ranges
CARBON_MIN, CARBON_MAX = -1.0, 1.0
PM_MIN, PM_MAX = -0.4, 0.4

# Ten classes per variable
NCLASS = 10

# Map extent
LAT_MIN, LAT_MAX = -60, 83

# Output resolution
DPI = 600

# =========================================================
# 3. Quadrant colors
#    Lower left: red; lower right: green; upper left: purple; upper right: blue
# =========================================================

# =========================================================
# 3. Four quadrant palettes
# =========================================================

# Upper left: Adsorption ↑ / Carbon ↓
PURPLE_5 = ["#FBFAFC", "#DADAEB", "#9D99C7", "#6A52A3", "#40027D"]

# Lower right: Carbon ↑ / Adsorption ↓
GREEN_5  = ["#F5FAF2", "#C5E8BE", "#74C477", "#238C46", "#00451C"]

# Lower left: Co-decrease
RED_5    = ["#FFF5F0", "#FCBBA2", "#FA6948", "#CC181E", "#6B000C"]

# Upper right: Co-increase
BLUE_5   = ["#F7FBFF", "#C7DCF0", "#6BAFD6", "#2170B5", "#08326E"]



# Background and boundaries
OCEAN_COLOR = "#FFFFFF"
COAST_COLOR = "#AAB2B7"
BORDER_COLOR = "#D0D5D8"
FRAME_COLOR = "black"

# =========================================================
# 4. Global plot style
# =========================================================

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["font.size"] = 10
mpl.rcParams["axes.titlesize"] = 12
mpl.rcParams["axes.labelsize"] = 10
mpl.rcParams["figure.facecolor"] = "white"
mpl.rcParams["savefig.facecolor"] = "white"


# =========================================================
# 5. Utility functions
# =========================================================
def hex_list_to_rgb01(hex_list):
    out = []
    for c in hex_list:
        c = c.lstrip("#")
        out.append(np.array([int(c[i:i+2], 16) for i in (0, 2, 4)], dtype=float) / 255.0)
    return np.array(out)


def interpolate_palette(hex_list, n=10):
    """
    Interpolate a five-color palette to n colors.
    """
    base = hex_list_to_rgb01(hex_list)
    xp = np.linspace(0, 1, len(base))
    xnew = np.linspace(0, 1, n)

    out = np.zeros((n, 3), dtype=float)
    for k in range(3):
        out[:, k] = np.interp(xnew, xp, base[:, k])
    return out
def hex_to_rgb01(hex_color):
    hex_color = hex_color.lstrip("#")
    return np.array([int(hex_color[i:i+2], 16) for i in (0, 2, 4)], dtype=float) / 255.0


def read_raster(path):
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        nodata = src.nodata
        if nodata is not None:
            arr[arr == nodata] = np.nan
        arr[~np.isfinite(arr)] = np.nan

        meta = {
            "shape": src.shape,
            "transform": src.transform,
            "crs": src.crs,
            "width": src.width,
            "height": src.height,
            "bounds": src.bounds
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
    """
    Equal-interval classification over a fixed range:
    - Values below vmin are assigned to the first class.
    - Values above vmax are assigned to the last class.
    Returns:
        idx   -> 0 ~ nclass-1
        edges -> Class boundaries
    """
    edges = np.linspace(vmin, vmax, nclass + 1)
    arr_clip = np.clip(arr, vmin, vmax)

    # digitize with internal edges
    idx = np.digitize(arr_clip, bins=edges[1:-1], right=False)
    idx = np.clip(idx, 0, nclass - 1)

    idx[~np.isfinite(arr)] = -1
    return idx, edges


def build_bivariate_lut(nclass, x_min, x_max, y_min, y_max):
    """
    Build a 10×10 bivariate color lookup table.
    Use the four quadrant palettes:
    - Upper left purple
    - Upper right blue
    - Lower left red
    - Lower right green
    """

    purple = interpolate_palette(PURPLE_5, n=5)
    green  = interpolate_palette(GREEN_5,  n=5)
    red    = interpolate_palette(RED_5,    n=5)
    blue   = interpolate_palette(BLUE_5,   n=5)

    x_edges = np.linspace(x_min, x_max, nclass + 1)
    y_edges = np.linspace(y_min, y_max, nclass + 1)

    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2

    lut = np.zeros((nclass, nclass, 4), dtype=np.float32)

    for iy, y in enumerate(y_centers):
        for ix, x in enumerate(x_centers):

            # ---------- Determine the quadrant ----------
            if x <= 0 and y >= 0:
                palette = purple   # Upper left
            elif x > 0 and y >= 0:
                palette = blue     # Upper right
            elif x <= 0 and y < 0:
                palette = red      # Lower left
            else:
                palette = green    # Lower right

            # ---------- Determine the level from the center ----------
            # Map the maximum normalized distance along x and y to five levels.
            if x <= 0:
                dx = abs(x) / abs(x_min) if x_min != 0 else 0
            else:
                dx = abs(x) / abs(x_max) if x_max != 0 else 0

            if y <= 0:
                dy = abs(y) / abs(y_min) if y_min != 0 else 0
            else:
                dy = abs(y) / abs(y_max) if y_max != 0 else 0

            d = max(dx, dy)  # Nested squares extending outward from the center
            d = np.clip(d, 0, 1)

            # Map to five colors: lighter at the center, darker at the edges.
            level = int(np.floor(d * 5))
            if level == 5:
                level = 4

            rgb = palette[level]

            lut[iy, ix, :3] = rgb
            lut[iy, ix, 3] = 1.0

    return lut, x_edges, y_edges


def make_rgba_from_classes(x_idx, y_idx, lut):
    """
    Generate a discrete RGBA map from class indices.
    """
    h, w = x_idx.shape
    rgba = np.zeros((h, w, 4), dtype=np.float32)

    valid = (x_idx >= 0) & (y_idx >= 0)
    rgba[valid] = lut[y_idx[valid], x_idx[valid]]

    return rgba


def add_classified_legend_square(fig, lut, box=(0.82, 0.14, 0.15, 0.15)):
    """
    Draw a square legend.
    """
    ax_leg = fig.add_axes(box)

    ax_leg.imshow(
        lut,
        origin="lower",
        extent=[0, NCLASS, 0, NCLASS],
        interpolation="nearest"
    )

    ax_leg.set_aspect("equal", adjustable="box")
    ax_leg.set_xlim(0, NCLASS)
    ax_leg.set_ylim(0, NCLASS)

    for spine in ax_leg.spines.values():
        spine.set_linewidth(0.8)
        spine.set_edgecolor("#B3B9BD")

    # Tick labels in data units
    ax_leg.set_xticks([0, 2.5, 5, 7.5, 10])
    ax_leg.set_xticklabels(["-10", "-5", "0", "5", "10"],
                           fontsize=11, color="#5A5F63")

    ax_leg.set_yticks([0, 2.5, 5, 7.5, 10])
    ax_leg.set_yticklabels(["-5", "-2.5", "0", "2.5", "5"],
                           fontsize=11, color="#5A5F63")

    ax_leg.tick_params(axis="both", which="major",
                       labelsize=11, width=0.6, length=3, colors="#5A5F63")

    ax_leg.set_xlabel("Carbon stock trend (t/ha/decade)", fontsize=11, labelpad=4, color="#2D3134")
    ax_leg.set_ylabel("Adsorption trend (g/m$^{2}$/decade)",
                      fontsize=11, labelpad=4, color="#2D3134")

    # Nested white squares
    rect1 = Rectangle((4, 4), 2, 2, fill=False, lw=1.0, ec=(1, 1, 1, 0.95))
    rect2 = Rectangle((3, 3), 4, 4, fill=False, lw=1.0, ec=(1, 1, 1, 0.72))
    rect3 = Rectangle((2, 2), 6, 6, fill=False, lw=1.0, ec=(1, 1, 1, 0.48))

    ax_leg.add_patch(rect1)
    ax_leg.add_patch(rect2)
    ax_leg.add_patch(rect3)

    ax_leg.text(0.35, 9.65, "Adsorption ↑\nCarbon ↓",
                ha="left", va="top", fontsize=11, color="white")

    ax_leg.text(9.65, 9.65, "Co-increase",
                ha="right", va="top", fontsize=11, color="white")

    ax_leg.text(0.35, 0.35, "Co-decrease",
                ha="left", va="bottom", fontsize=11, color="white")

    ax_leg.text(9.65, 0.35, "Carbon ↑\nAdsorption ↓",
                ha="right", va="bottom", fontsize=11, color="white")

    return ax_leg

def add_gradient_explanation(fig, box=(0.80, 0.53, 0.24, 0.15)):
    """
    Draw palette strips above the main legend.
    Show change intensity within each quadrant.
    The main legend is drawn separately.
    """
    ax_note = fig.add_axes(box)
    ax_note.set_xlim(0, 1)
    ax_note.set_ylim(0, 1)
    ax_note.axis("off")

    # Title
    ax_note.text(
        0.44, 1.00, "Trend intensity",
        ha="center", va="top",
        fontsize=11, color="#333333"
    )

    labels = ["Ads. ↑ - Car. ↓", "Car. ↑ - Ads. ↓", "Co-decrease", "Co-increase"]
    palettes = [PURPLE_5, GREEN_5, RED_5, BLUE_5]

    # Row positions
    y_positions = [0.72, 0.52, 0.32, 0.12]

    # Color strip positions
    x0 = 0.29
    bar_w = 0.3
    bar_h = 0.13

    for label, palette, y in zip(labels, palettes, y_positions):
        # Labels on the left
        ax_note.text(
            0.03, y + bar_h / 2,
            label,
            ha="left", va="center",
            fontsize=11, color="#333333"
        )

        # Color strips (light to dark)
        colors = hex_list_to_rgb01(palette)[np.newaxis, :, :]
        ax_note.imshow(
            colors,
            extent=[x0, x0 + bar_w, y, y + bar_h],
            origin="lower",
            aspect="auto",
            interpolation="nearest"
        )

        # Separators between color levels
        for k in range(1, 5):
            xx = x0 + bar_w * k / 5
            ax_note.plot([xx, xx], [y, y + bar_h],
                         color="white", lw=0.8, alpha=0.75)

    # Arrow below the strips
    ax_note.annotate(
        "",
        xy=(0.60, -0.02), xytext=(0.29, -0.02),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", lw=1.2, color="#333333"),
        annotation_clip=False
    )

    ax_note.text(
        0.34, -0.08, "Lower",
        ha="center", va="top",
        fontsize=11, color="#333333",
        transform=ax_note.transAxes
    )
    ax_note.text(
        0.54, -0.08, "Higher",
        ha="center", va="top",
        fontsize=11, color="#333333",
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

ref_arr, ref_meta = read_raster(scenario_files["historical"]["carbon"])

data = {}

for scn in scenario_order:
    carbon_arr, carbon_meta = read_raster(scenario_files[scn]["carbon"])
    pm_arr, pm_meta = read_raster(scenario_files[scn]["pm25"])

    if not same_grid(carbon_meta, ref_meta):
        carbon_arr = reproject_to_match(carbon_arr, carbon_meta, ref_meta, resampling=Resampling.bilinear)
        carbon_meta = ref_meta

    if not same_grid(pm_meta, ref_meta):
        pm_arr = reproject_to_match(pm_arr, pm_meta, ref_meta, resampling=Resampling.bilinear)
        pm_meta = ref_meta

    data[scn] = {
        "carbon": carbon_arr,
        "pm25": pm_arr,
        "meta": ref_meta,
        "title": scenario_files[scn]["title"]
    }


# =========================================================
# 8. Build the 10×10 discrete bivariate color table
# =========================================================

lut, carbon_edges, pm_edges = build_bivariate_lut(
    NCLASS,
    CARBON_MIN, CARBON_MAX,
    PM_MIN, PM_MAX
)


# =========================================================
# 9. Classify values and generate RGBA maps
# =========================================================

for scn in scenario_order:
    carbon_idx, _ = classify_fixed(data[scn]["carbon"], CARBON_MIN, CARBON_MAX, NCLASS)
    pm_idx, _ = classify_fixed(data[scn]["pm25"], PM_MIN, PM_MAX, NCLASS)

    rgba = make_rgba_from_classes(carbon_idx, pm_idx, lut)
    data[scn]["rgba"] = rgba


# =========================================================
# =========================================================
# 10. Plot maps using a manual panel layout
# =========================================================

proj = ccrs.Robinson(central_longitude=0)
pc = ccrs.PlateCarree()

# Wide canvas with space for legends on the right
fig = plt.figure(figsize=(18.5, 8.1), facecolor="white")

# Position map panels to leave room for legends on the right.
panel_pos = {
    "historical": [0.03, 0.57, 0.36, 0.40],
    "ssp126":     [0.42, 0.57, 0.36, 0.40],
    "ssp245":     [0.03, 0.12, 0.36, 0.40],
    "ssp585":     [0.42, 0.12, 0.36, 0.40],
}



panel_labels = ["A", "B", "C", "D"]

for i, scn in enumerate(scenario_order):
    ax = fig.add_axes(panel_pos[scn], projection=proj)

    ax.set_facecolor("white")
    ax.add_feature(cfeature.OCEAN.with_scale("110m"), facecolor=OCEAN_COLOR, zorder=0)

    bounds = data[scn]["meta"]["bounds"]
    extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]

    ax.imshow(
        data[scn]["rgba"],
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

    ax.set_extent([-180, 180, LAT_MIN, LAT_MAX], crs=pc)

    try:
        ax.spines["geo"].set_linewidth(1.2)
        ax.spines["geo"].set_edgecolor("black")
    except Exception:
        pass

    ax.set_title(data[scn]["title"], pad=4, color="#222222", fontweight="normal")

    ax.text(
        0.018, 0.985, panel_labels[i],
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=18, fontweight="bold", color="#222222"
    )

# Place the legends to the right of the map panels.
# Legend positions and sizes use figure coordinates.
add_classified_legend_square(fig, lut, box=(0.77, 0.2, 0.24, 0.24))
add_gradient_explanation(fig, box=(0.80, 0.53, 0.24, 0.15))
fig.savefig(output_png, dpi=DPI, bbox_inches="tight", facecolor="white")
fig.savefig(output_pdf, dpi=DPI, bbox_inches="tight", facecolor="white")
plt.show()

print("\nPlotting complete.")
print("PNG:", output_png)
print("PDF:", output_pdf)