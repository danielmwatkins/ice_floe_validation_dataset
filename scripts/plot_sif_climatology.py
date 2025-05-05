"""Generate figure 2, the sea ice fraction climatology for the study regions."""
import pandas as pd
import proplot as pplt
import numpy as np
import warnings
warnings.simplefilter('ignore')

regions = pd.read_csv('../data/metadata/region_definitions.csv', index_col=0)
regions['print_title'] = [c.replace('_', ' ').title().replace('Of', 'of') for c in regions.index]

colors = {region: c['color'] for region, c in zip(
            regions.index,
            pplt.Cycle('colorblind10', len(regions)))}


idx = 1
for region in ['greenland_sea', 'barents_kara_seas', 'laptev_sea', 'sea_of_okhostk', 'east_siberian_sea', 'bering_chukchi_seas', 'beaufort_sea', 'hudson_bay', 'baffin_bay']:
    regions.loc[region, 'idx_num'] = idx
    idx += 1
regions = regions.sort_values('idx_num')

regions.loc['bering_chukchi_seas', 'print_title'] = 'Bering-Chukchi Seas'
regions.loc['barents_kara_seas', 'print_title'] = 'Barents-Kara Seas'

sic_timeseries = pd.read_csv('../data/metadata/daily_sea_ice_fraction.csv', index_col=0, parse_dates=True)

df['doy'] = df.index.dayofyear

p90 = df.groupby('doy').quantile(0.9)
for region in p90:
    print(region, p90[region][p90[region] > 0.05].index.max())

linestyles = {region: ls for region, ls in zip(regions.index,
                        ['-', '-.', '--', '-', '-.', '--', '-.', '-', '--'])}

fig, ax = pplt.subplots(width=6, height=4)
for region in regions.index:
    ax.plot(sic_timeseries[region].groupby(sic_timeseries[region].index.dayofyear).median(),
            shadedata=[sic_timeseries[region].groupby(sic_timeseries[region].index.dayofyear).quantile(0.25),
                       sic_timeseries[region].groupby(sic_timeseries[region].index.dayofyear).quantile(0.75)],
            fadedata=[sic_timeseries[region].groupby(sic_timeseries[region].index.dayofyear).quantile(0.1),
                       sic_timeseries[region].groupby(sic_timeseries[region].index.dayofyear).quantile(0.9)],
            c=colors[region], lw=2, ls=linestyles[region], 
            label='({n}) {t}'.format(n=int(regions.loc[region, 'idx_num']), t=regions.loc[region, 'print_title']))
ax.legend(loc='b', ncols=3,lw=2, order='F')
ax.format(ylabel='Sea Ice Fraction', xlabel='Day of Year', ylim=(0, 1),
         ylocator=np.arange(0.1, 0.91, 0.2), xlim=(60, 274))
fig.save('../figures/fig03_sea_ice_fraction.png', dpi=300)