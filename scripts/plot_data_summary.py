import proplot as pplt
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import os

df = pd.read_csv('../data/validation_dataset/validation_dataset.csv', parse_dates=['start_date'])
df['month'] = df['start_date'].dt.month
df = df.loc[df.satellite=='aqua']
df.index = [str(x).zfill(3) for x in df.case_number]

import pandas as pd
fsd_data = {}
for satellite in ['aqua', 'terra']:
    fsd_data[satellite] = {}
    for file in os.listdir('../data/validation_dataset/property_tables/' + satellite):
        if 'csv' in file:
            df_fsd = pd.read_csv('../data/validation_dataset/property_tables/' + satellite + '/' + file)
            if len(df_fsd) > 0:
                df_fsd['region'] = file.split('-')[1]
                df_fsd['date'] = pd.to_datetime(file.split('-')[2])
                fsd_data[satellite][file.split('-')[0]] = df_fsd


df_aqua_fsd = pd.concat(fsd_data['aqua']).reset_index()
df_terra_fsd = pd.concat(fsd_data['terra']).reset_index()
df_aqua_fsd.rename({'level_0': 'case_number'}, axis=1, inplace=True)
df_terra_fsd.rename({'level_0': 'case_number'}, axis=1, inplace=True)
df_aqua_fsd.head()


regions = pd.read_csv('../data/metadata/region_definitions.csv', index_col=0)

colors = {region: c['color'] for region, c in zip(
            regions.index,
            pplt.Cycle('dark2', len(regions)))}
linestyles = {region: ls for region, ls in zip(regions.index,
                        ['-', '-.', '--', '-', '-.', '--', '-.', '-', '--'])}

regions['print_title'] = [c.replace('_', ' ').title().replace('Of', 'of') for c in regions.index]
regions = regions.sort_values('center_lon')

for idx, row in regions.iterrows():
    if row.print_title == 'Barents Kara Seas':
        regions.loc[idx, 'print_title'] = 'Barents-Kara Seas'
    if row.print_title == 'Bering Chukchi Seas':
        regions.loc[idx, 'print_title'] = 'Bering-Chukchi Seas'

pplt.rc['cartopy.circular'] = False
pplt.rc['reso'] = 'med'
crs = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
fig, axs = pplt.subplots([[1, 3], [2, 3]], height=6, refwidth=4, refnum=3, share=False,
    proj={3: 'npstere'}, proj_kw={3: {'lon_0': -45}})

ax = axs[2]
ax.format(land=True, landzorder=0)

ice_floes_cases = df.loc[df.visible_floes == 'yes']
no_floes_cases = df.loc[df.visible_floes == 'no']

for idx, region, lat, lon in zip(range(len(regions)), regions.index, regions.center_lat, regions.center_lon):

    xbox = np.array(regions.loc[region, ['left_x', 'left_x', 'right_x', 'right_x', 'left_x']].astype(float))
    ybox = np.array(regions.loc[region, ['lower_y', 'upper_y', 'upper_y', 'lower_y', 'lower_y']].astype(float))
    
    ax.plot(xbox, ybox, transform=ccrs.CRS('epsg:3413'),
            label='({n}) {t}'.format(n=idx + 1, t=regions.loc[region, 'print_title']), 
               color=colors[region], ls=linestyles[region], m='', zorder=5, lw=1.5)
    
ax.set_extent([-3.5e6, 2.8e6, -4e6, 4.4e6], crs=crs)

for idx, case in df.iterrows():
    x0 = case.center_x
    y0 = case.center_y
    left = x0 - 50e3
    right = x0 + 50e3
    bottom = y0 - 50e3
    top = y0 + 50e3 
    region = case.region
    if str(case.case_number).zfill(3) in no_floes_cases.index:
        c = 'light gray'
        z = 1
        ax.plot(x0, y0, m='s', color=c, alpha=0.5, s=5, transform=crs, label='', zorder=0)
    else:
        c = colors[region]
        z = 20
    ax.plot([left, left, right, right, left],
        [bottom, top, top, bottom, bottom], label='',
            transform=crs, color=c, lw=1, zorder=z)
   
ax.plot([],[],m='s', color='gray', label='No visible floes', lw=0)
ax.format(title='Location of samples')
fig.legend(ncols=1, alpha=1, loc='r', order='F')

# Bar chart
ax = axs[0]

df = pd.read_csv('../data/validation_dataset/validation_dataset.csv', parse_dates=['start_date'])
df['month'] = df['start_date'].dt.month
df.index = [str(x).zfill(3) for x in df.case_number]

df.loc[:, 'visible_sea_ice'] = df.loc[:, 'visible_sea_ice'].where(df.loc[:, 'visible_sea_ice']=='yes')
df.loc[:, 'visible_landfast_ice'] = df.loc[:, 'visible_landfast_ice'].where(df.loc[:, 'visible_landfast_ice']=='yes')
df.loc[:, 'visible_floes'] = df.loc[:, 'visible_floes'].where(df.loc[:, 'visible_floes']=='yes')
df.rename({'cloud_fraction_manual': 'Number of images',
                  'visible_sea_ice': 'Visible ice',
                  'visible_floes': 'Visible floes',
                  'visible_landfast_ice': 'Visible landfast'}, axis=1, inplace=True)

ax.bar(df.loc[:, ['month', 'Number of images', 'Visible ice', 'Visible floes', 'Visible landfast']].groupby('month').count())
ax.legend(loc='b', ncols=2, fontsize=8)
ax.format(ylabel='Count', xlabel='', xlocator=np.arange(3, 10), title='Summary of metadata',
          xformatter=['Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'],
          xtickminor=False)

# FSD
ax = axs[1]
for data, label in zip([df_aqua_fsd, df_terra_fsd], ['Aqua', 'Terra']):
    bins = np.logspace(1, np.log(np.sqrt(100)), 25)*0.25
    x, xb = np.histogram(np.sqrt(data['area'])*0.25, bins)
    
    ax.plot(0.5*(bins[:-1] + bins[1:]), x, marker='.', label=label)
    ax.format(yscale='log', xscale='log', ylabel='Count', xlocator=(1, 5, 10, 20, 40),
              xlim=(2.5, 50), xlabel='Length scale (km)', title='Floe size distribution')
ax.legend(loc='ur', ncols=1)


fig.format(abc=True)
fig.save('../figures/fig_04_random_sample_summary.png', dpi=300)

