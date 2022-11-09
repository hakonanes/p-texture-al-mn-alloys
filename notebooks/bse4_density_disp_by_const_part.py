# coding: utf-8
#
# Density of dispersoids by constituent particles
# 
# Håkon Wiik Ånes (hakon.w.anes@ntnu.no)
# 
# Procedure:
# 1. Get label map for particles
# 2. Expand label map using `skimage.segmentation.expand_labels()`, which uses `scipy.ndimage.distance_transform_edt()`
# 3. Get masks of region inside precipitate-free zones (PFZs), excluding big particles, and outside
# 4. Calculate:
#     1. Density of particles inside and outside PFZs
#     2. Area weighted particle size inside and outside PFZs

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from skimage.measure import label
from skimage.segmentation import expand_labels

from mapregions import MapRegions


# Parse input parameters
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("sample")
parser.add_argument("dset")
parser.add_argument("ct")
parser.add_argument("tpfz")

args = parser.parse_args()

sample = str(args.sample)
dset = str(args.dset)
constituent_threshold = float(args.ct)
pfz_extent = float(args.tpfz)

step_size = dict(bse=1 / 39.2)

# Particle size thresholds in um for constituent particles (PSN) and
# dispersoids (Smith-Zener drag)
dispersoid_threshold_max = 0.24  # um
dispersoid_threshold_min = 0.03  # um

# Extent of PFZ out from constituent particles in um
pfz_extent_px = int(np.ceil(pfz_extent * step_size["bse"] * 1e3))

# Load data
dir_data = Path("/home/hakon/phd/data/p/prover") / sample / str(dset) / "partdet"
img = np.load(dir_data / "bse_labels_filled_filtered.npy")

# Get label map of particles
seg1 = label(img)
mreg = MapRegions(seg1, dx=step_size["bse"], dy=step_size["bse"], background_label=0, scan_unit="um")

# Get big particles
mask_big = 0.816 * 2 * mreg.equivalent_radius >= constituent_threshold
mreg_big = mreg[mask_big]

# Extend label map of big particles
seg2 = expand_labels(mreg_big.label_map, distance=pfz_extent_px)
mask = seg2 > 0

# Mask and area of precipitate-free zones (PFZs)
pfz = mask.astype(int) - (~mreg_big.is_background_map).astype(int)
pfz = pfz.astype(bool)
area_pfz = np.sum(pfz) * step_size["bse"] ** 2

# Mask and area of pixels outside PFZs (excluding constituent particles themselves)
not_pfz = ~pfz
not_pfz[~mreg_big.is_background_map] = False
area_not_pfz = np.sum(not_pfz) * step_size["bse"] ** 2

# Sanity checks (raise assertion error(s) if false)
# Mask sizes
assert np.sum([pfz, not_pfz, ~mreg_big.is_background_map]) == img.size
# Areas
assert np.allclose(
    np.sum([area_pfz, area_not_pfz, np.sum(~mreg_big.is_background_map) * step_size["bse"] ** 2]),
    img.size * step_size["bse"] ** 2,
    atol=1e-5
)

# Get labels of constituent particles, dispersoids inside PFZs and dispersoids outside
# PFZs
# Constituent particles
labels_big = mreg.label[mask_big]

# All particles inside PFZs (regardless of size), excluding constituent particles
labels_inside_all = mreg.label_map[pfz]
labels_inside = np.unique(labels_inside_all)[1:]

# Outside PFZs
mask_outside_all = ~np.isin(mreg.label, np.concatenate([labels_big, labels_inside]))
labels_outside = mreg.label[mask_outside_all]

# Sanity check (raises assertion error if false)
assert np.allclose(mreg.label, np.sort(np.concatenate([labels_big, labels_inside, labels_outside])))

# Size and area of dispersoids inside PFZ
r_all = 0.816 * 2 * mreg.equivalent_radius
mask_dispersoid_all = (r_all >= dispersoid_threshold_min) & (r_all <= dispersoid_threshold_max)
mask_inside_all = np.isin(mreg.label, labels_inside)
mreg_inside = mreg[mask_inside_all & mask_dispersoid_all]
n_inside = mreg_inside.label.size
size_inside = 0.816 * 2 * mreg_inside.equivalent_radius
area_inside = mreg_inside.area

# Size and area of dispersoids outside PFZ
mreg_outside = mreg[mask_outside_all & mask_dispersoid_all]
n_outside = mreg_outside.label.size
size_outside = 0.816 * 2 * mreg_outside.equivalent_radius
area_outside = mreg_outside.area

# Concatenate data and save to one file
columns = ["r", "area", "is_constituent", "is_dispersoid", "inside_pfz"]
particles = pd.DataFrame(
    np.zeros((mreg.label.size, len(columns))),
    columns=columns,
)

particles.loc[mask_big, "is_constituent"] = 1
particles.loc[mask_dispersoid_all, "is_dispersoid"] = 1
particles.loc[mask_inside_all, "inside_pfz"] = 1

particles["r"] = mreg.equivalent_radius
particles["area"] = mreg.area

# Sanity checks
assert np.allclose(
    np.sum(particles["inside_pfz"].astype(bool) & particles["is_dispersoid"].astype(bool)),
    n_inside
)
assert np.allclose(particles["is_constituent"].sum(), mreg_big.label.size)

# Write to file
np.savetxt(
    dir_data / f"pfz_stats_threshold{constituent_threshold:.2f}_pfz{pfz_extent:.2f}.csv",
    particles.values,
    fmt="%.5f,%.5f,%i,%i,%i",
    header=(
        "r,area,is_constituent,is_dispersoid,inside_pfz\n"
        f"area pfz: {area_pfz:.5f} \n"
        f"area not pfz: {area_not_pfz:.5f}"
    )
)
