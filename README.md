# Ice Floe Validation Dataset
The Ice Floe Validation Dataset is a collection of randomly sampled 100 km by 100 km MODIS scenes from the _Aqua_ and _Terra_ satellites. This repository contains the dataset itself along with the notebooks and scripts used to produce the dataset.

# Contents
## /data
The data directory has subdirectories
1. `/metadata` containing specification/config files
2. `/modis` containing truecolor, falsecolor, landmask, and cloudfraction images
3. `/masie` containing images from the MASIE 4 km data product, cropped and aligned with the MODIS projection via nearest-neighbors interpolation
4. `/validation_dataset` containing the binary PNG, labeled TIFF files, and Photoshop files. The subdirectories are
   * binary_floes
   * binary_landfast
   * binary_landmask
   * labeled_floes
   * labeling_psd_files
   * property_tables

The files all follow a strict, informative naming convention. 
1. MODIS and MASIE images: `NNN-region_name-100km-YYYYMMDD.satellite_name.image_type.250m.tiff` NNN is the zero-padded, three digit case number. Each case specifies a date and location, and contains two satellite images (Aqua and Terra). MODIS image types are `cloudfraction` (RGB), `truecolor` (Bands 1-4-3), `falsecolor` (Bands 7-2-1), and `landmask`. MASIE images have `masie` instead of satellite, since there's only one MASIE image each day, and have image types `seaice` and `landmask`.
2. Validation images: `NNN-region_name-YYYYMMDD-satellite-image_type.file_type`. Image types with file type `png` are `binary_floes`, `binary_landfast`, `binary_landmask`. Labeled floe masks are `tiff` with image type `labeled_floes` and have integer labels corresponding to the labels in the property tables.

Finally, the `/property_tables` contains the main metadata table with all case information and manual image assessments `ice_floe_validation_dataset.csv` as well as automatically generated property tables for each image with ice floes in it, and a table of hybrid automatic/manual floe matching for paired analyses.

## /notebooks
1. `bookkeeping.ipynb` Scripts and diagnostics to check the completeness of the datasets.
2. `sample_selection.ipynb` Code for iterative random sampling

## /scripts
1. setup_locations.py
2. calculate_sea_ice_fraction.py
3. binarize_landmasks.py
4. extract_features_and_pair_floes.py
5. plot_locations.py
6. plot_sif_climatology.py

## /figures
The figures directory has diagnostic and explanatory figures, including the sea ice fraction climatology, the map of analyzed cases, and a summary of the metadata.

# TBD options
1. Filetypes: the labeled images could be bundled together into a single h5 file per case.
