"""
Generate table with region definitions
"""

import pandas as pd
import pyproj
import numpy as np
import os
import warnings

warnings.simplefilter('ignore')

regions = {
    'greenland_sea': 
        {'left_x': 290374,
         'right_x': 1244051,
         'upper_y': -725683,
         'lower_y': -2234665
        },
    'barents_kara_seas':
        {'left_x': 700817,
         'right_x': 2398925,
         'upper_y': 1173623,
         'lower_y': -725683
        },
    'laptev_sea': 
        {'left_x': -152260, 
         'right_x': 700817,
         'upper_y': 2119251,
         'lower_y': 895970
        },
    'east_siberian_sea':
        {'left_x': -1472117,
         'right_x': -152260,
         'upper_y': 2038772,
         'lower_y': 753120
        },
    'bering_chukchi_seas':
        {'left_x': -1472117 - 1500e3,
         'right_x': -1472117,
         'upper_y': 753120 + 1500e3,
         'lower_y': 753120
        },
    'beaufort_sea': 
        {'left_x':-2383879,
         'right_x': -883879,
         'upper_y': 753120,
         'lower_y': -303168
        },
    'baffin_bay':
        {'left_x': -987231,
         'right_x': -97937,
         'upper_y': -818234,
         'lower_y': -818234 - 1500e3, 
        },
    'hudson_bay':
        {'left_x': -2795941,
         'right_x': -1655207,
         'upper_y': -1868686,
         'lower_y': -3635000
        }, 
    'sea_of_okhostk':
        {'left_x': -1549890,
         'right_x': -49890,
         'upper_y': 4214708,
         'lower_y': 2714708
        }
}

crs0 = pyproj.CRS('WGS84')
crs1 = pyproj.CRS('epsg:3413')
transformer_ll = pyproj.Transformer.from_crs(crs1, crs_to=crs0, always_xy=True)

for region in regions:

    x0 = 0.5*(regions[region]['left_x'] + regions[region]['right_x'])
    y0 = 0.5*(regions[region]['upper_y'] + regions[region]['lower_y'])
    lon0, lat0 = transformer_ll.transform(x0, y0)
    
    regions[region]['center_x'] = np.round(x0, 0)
    regions[region]['center_y'] = np.round(y0, 0)
    regions[region]['center_lon'] = np.round(lon0, 4)
    regions[region]['center_lat'] = np.round(lat0, 4)

location_df = pd.DataFrame(regions).T
print(location_df)


location_df.loc[:, ['center_lat', 'center_lon', 'center_x', 'center_y', 'left_x', 'right_x', 'lower_y', 'upper_y']].to_csv('../data/metadata/region_definitions.csv')