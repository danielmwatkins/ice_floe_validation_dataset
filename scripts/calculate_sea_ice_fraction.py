"""For the data selection step, we need to know what times of year there is sea ice present in each of our study areas. 
To that end, we use the NSIDC Sea Ice Concentration Climate Data Record (CDR). For each day between the start of March
and the end of September, we open the corresponding daily CDR, and for each study region in the region_definitions file,
we select the area inside the region bounding box, and identify the pixels with SIC > 0.15 and the pixels inside the land
and coast masks. We then calculate the sea ice fraction as SIC Pixels / (Total Pixels - Land Pixels). We then save these
time series to data/metadata/daily_sea_ice_fraction.csv."""
import numpy as np
import pandas as pd
import xarray as xr
import os

sic_dataloc = '/Users/dwatkin2/Documents/research/data/nsidc_daily_cdr_v4/' #check v5 variables for update
saveloc = '../data/metadata/'
regions = pd.read_csv('../data/metadata/region_definitions.csv', index_col=0)
print(regions)
date_index = []
results = []

def compute_sic(left_x, right_x, bottom_y, top_y, sic_data):
    """Computes the sea ice extent as a fraction of total area within the region bounded
    by <left_x>, <right_x>, <bottom_y>, and <top_y> using the netcdf file <sic_data>. Assumes
    that sic_data is the NSIDC SIC CDR."""

    # Note: this is in polar stereographic coordinates
    # We should update to an equal area method if possible
    
    x_idx = (sic_data.xgrid >= left_x) & (sic_data.xgrid <= right_x)
    y_idx = (sic_data.ygrid >= bottom_y) & (sic_data.ygrid <= top_y)
    
    with_ice = ((sic_data.sel(x=x_idx, y=y_idx)['cdr_seaice_conc'] > 0.15) & \
                (sic_data.sel(x=x_idx, y=y_idx)['cdr_seaice_conc'] <= 1))
    coast_mask = (sic_data.sel(x=x_idx, y=y_idx)['cdr_seaice_conc'] > 1).sum() 
    total_area_pixels = np.prod(with_ice.shape)
    sic_area_pixels = with_ice.sum().data
    return np.round(sic_area_pixels/(total_area_pixels - coast_mask.data), 3)

# TBD: Update to 2006 for the MASIE data coverage?
for year in range(2003, 2023):
    files = os.listdir(sic_dataloc + str(year))
    for file in files:
        if '.nc' in file:
            date = pd.to_datetime(file.split('_')[4], format='%Y%m%d')
            if (date.month >= 3) & (date.month <= 9):
                with xr.open_dataset(sic_dataloc + str(year) + '/' + file) as ds_sic:
                    date_index.append(date)
                    temp_results = []
                    for region in regions.index:
                        temp_results.append(
                            compute_sic(regions.loc[region, 'left_x'], regions.loc[region, 'right_x'], 
                            regions.loc[region, 'lower_y'], regions.loc[region, 'upper_y'], ds_sic))
                results.append(temp_results)
    
sic_timeseries = pd.DataFrame(results, index=date_index, columns=regions.index)
sic_timeseries = sic_timeseries.sort_index()
sic_timeseries.to_csv(saveloc + 'daily_sea_ice_fraction.csv')