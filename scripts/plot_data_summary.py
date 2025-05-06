import proplot as pplt
import pandas as pd
import numpy as np

fig, ax = pplt.subplots()

df = pd.read_csv('../data/validation_dataset/validation_dataset.csv', parse_dates=['start_date'])
df['month'] = df['start_date'].dt.month
df.loc[:, 'visible_sea_ice'] = df.loc[:, 'visible_sea_ice'].where(df.loc[:, 'visible_sea_ice']=='yes')
df.loc[:, 'visible_landfast_ice'] = df.loc[:, 'visible_landfast_ice'].where(df.loc[:, 'visible_landfast_ice']=='yes')
df.loc[:, 'visible_floes'] = df.loc[:, 'visible_floes'].where(df.loc[:, 'visible_floes']=='yes')
df.rename({'cloud_fraction_manual': 'Number of images',
                  'visible_sea_ice': 'Sea ice',
                  'visible_floes': 'Sea ice floes',
                  'visible_landfast_ice': 'Landfast ice'}, axis=1, inplace=True)

ax.bar(df.loc[:, ['month', 'Number of images', 'Sea ice', 'Sea ice floes', 'Landfast ice']].groupby('month').count())
ax.legend(loc='b', ncols=2)
ax.format(ylabel='Count', xlabel='', xlocator=np.arange(3, 10),
          xformatter=['Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'],
          xtickminor=False)
fig.save('../figures/random_sample_summary_chart.png', dpi=300)



     
