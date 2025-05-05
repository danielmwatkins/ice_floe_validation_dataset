import pandas as pd
import proplot as pplt
import cartopy.crs as ccrs
import numpy as np
import warnings
import xarray as xr
import pyproj
warnings.simplefilter('ignore')

pplt.rc['reso'] = 'med'
pplt.rc['cartopy.circular'] = False

# Load data
# TBD: Update to version 5 of the NSIDC sea ice concentration climate data record
sic_data = []
for year in range(2000, 2022):
    ds = xr.open_dataset("~/Documents/research/data/nsidc_daily_cdr_v4/aggregate/seaice_conc_daily_nh_{y}_v04r00.nc".format(y=year))
    sic_data.append(ds)

ds_april = xr.concat([ds.sel(tdim=ds.time.dt.month == 4) for ds in sic_data], dim='tdim')
ds_sept = xr.concat([ds.sel(tdim=ds.time.dt.month == 9) for ds in sic_data], dim='tdim')


# Calculate mean
ds_april_mean = ds_april[['cdr_seaice_conc']].where(ds_april[['cdr_seaice_conc']] <= 1).mean(dim='tdim')
ds_sept_mean = ds_sept[['cdr_seaice_conc']].where(ds_sept[['cdr_seaice_conc']] <= 1).mean(dim='tdim')


# Setup the polar stereographic coordinate arrays
crs0 = pyproj.CRS('WGS84')
crs1 = pyproj.CRS('epsg:3413')
transformer_xy = pyproj.Transformer.from_crs(crs0, crs_to=crs1, always_xy=True)
transformer_ll = pyproj.Transformer.from_crs(crs1, crs_to=crs0, always_xy=True)
x, y = transformer_xy.transform(np.ravel(ds_april['longitude'].isel(tdim=0).data.squeeze()),
                                np.ravel(ds_april['latitude'].isel(tdim=0).data.squeeze()))
s = ds_april['longitude'].isel(tdim=0).data.squeeze().shape
X = np.reshape(x, s)
Y = np.reshape(y, s)

# Load region data
regions = pd.read_csv('../data/metadata/region_definitions.csv', index_col=0)
regions['print_title'] = [c.replace('_', ' ').title().replace('Of', 'of') for c in regions.index]


crs = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
fig, ax = pplt.subplots(width=4.5, proj='npstere', proj_kw={'lon_0': -45})
ax.format(land=True, color='k', boundinglat=52, landzorder=0, latmax=90, facecolor='w', lonlocator=10)
ax.set_extent([-3.5e6, 2.8e6, -4e6, 4.4e6], crs=crs)

for region in regions.index:
    xbox = np.array([regions.loc[region, coord] for coord in ['left_x', 'left_x', 'right_x', 'right_x', 'left_x']])
    ybox = np.array([regions.loc[region, coord] for coord in ['lower_y', 'upper_y', 'upper_y', 'lower_y', 'lower_y']])
    ax.plot(xbox, ybox, transform=crs, lw=2, color='pink5')

miz = (ds_april_mean['cdr_seaice_conc'] > 0.15) & (ds_april_mean['cdr_seaice_conc'] < 0.85)
april_miz = np.ma.masked_array(np.ones(miz.data.shape), mask=~miz)    

miz = (ds_sept_mean['cdr_seaice_conc'] > 0.15) & (ds_sept_mean['cdr_seaice_conc'] < 0.85)
sept_miz = np.ma.masked_array(np.ones(miz.data.shape), mask=~miz)    

apr_pack_ice = (ds_april_mean['cdr_seaice_conc'] >= 0.85) & (ds_april_mean['cdr_seaice_conc'] <= 1)
apr_pack_ice = np.ma.masked_array(np.ones(apr_pack_ice.data.shape), mask=~apr_pack_ice)    

sep_pack_ice = (ds_sept_mean['cdr_seaice_conc'] >= 0.85) & (ds_sept_mean['cdr_seaice_conc'] <= 1)
sep_pack_ice = np.ma.masked_array(np.ones(sep_pack_ice.data.shape), mask=~sep_pack_ice)    


ax.pcolormesh(X, Y, apr_pack_ice, vmin=0, vmax=1, color='blue1', alpha=1,
              transform=crs, label='')

ax.pcolormesh(X, Y, april_miz, vmin=0, vmax=1, color='blue4', alpha=1,
              transform=crs, label='')

ax.pcolormesh(X, Y, sep_pack_ice, vmin=0, vmax=1, color='orange1', alpha=1,
              transform=crs, label='')

ax.pcolormesh(X, Y, sept_miz, vmin=0, vmax=1, color='orange3', alpha=1,
              transform=crs, label='September MIZ')

h = [ax.plot([],[],marker='s', lw=0, color=color) for color in ['blue4', 'blue1', 'orange3', 'orange1']]
ax.legend(h, ['April MIZ', 'April Pack Ice', 'September MIZ', 'September Pack Ice'], loc='ur', ncols=1, alpha=1)
idx = 1
for region in regions.index:
    ax.text(regions.loc[region, 'left_x'] + 300e3,
            regions.loc[region, 'upper_y'] - 400e3, str(idx),
            transform=crs, bbox=True, bboxalpha=1,
            border=False, color='k', borderwidth=0,
            bboxstyle='circle', bboxcolor='w', zorder=10)
    idx += 1
ax.text(100, 69, 'Arctic Circle', color='lightgray', transform=ccrs.PlateCarree(), rotation=-40)

ax.plot(np.linspace(0, 360, 100), np.ones(100)*66.3, ls='-.', color='light gray')
fig.save('../figures/sample_locations_map.png', dpi=300)