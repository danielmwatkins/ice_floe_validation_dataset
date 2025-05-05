"""Loops through all the images and overlays the validation dataset masks for verification"""
import numpy as np
import os
import pandas as pd
import proplot as pplt
import rasterio as rio
from rasterio.plot import reshape_as_image
from sklearn.model_selection import KFold
from scipy.interpolate import interp1d

# Load the list of cloud clearing evaluation cases
dataloc = '../../ice_floe_validation_dataset/'
df = pd.read_csv(dataloc + '/data/validation_dataset/validation_dataset.csv')
df['case_number'] = [str(cn).zfill(3) for cn in df['case_number']]
df.groupby('region').count()
df['start_date'] = pd.to_datetime(df['start_date'].values)
df.index = [cn + '_' + sat for cn, sat in zip(df.case_number, df.satellite)]

def fname(case_data, imtype='labeled_floes'):
    """Generates filenames from rows in the overview table. imtype can be 'labeled_floes', 
    'binary_floes', 'binary_landfast', or 'binary_landmask', 'truecolor', or 'falsecolor'.
    The imtype determines whether a 'png' or 'tiff' is returned.
    """

    cn = case_data['case_number']
    date = pd.to_datetime(case_data['start_date']).strftime('%Y%m%d')
    region = case_data['region']
    sat = case_data['satellite']
    if 'binary' in imtype:
        return  '-'.join([cn, region, date, sat, imtype + '.png'])
        
    elif imtype in ['truecolor', 'falsecolor', 'cloudfraction', 'labeled_floes',]:
        prefix = '-'.join([cn, region, '100km', date])
        return '.'.join([prefix, sat, imtype, '250m', 'tiff'])

    elif imtype in ['seaice', 'landmask',]:
        prefix = '-'.join([cn, region, '100km', date])
        return '.'.join([prefix, 'masie', imtype, '250m', 'tiff'])     

    elif imtype == 'cloudfraction_numeric':
        
        return '-'.join([cn, region, date, sat, 'cloudfraction.csv'])

# Load raster data and masks
fc_dataloc = dataloc + 'data/modis/falsecolor/'
tc_dataloc = dataloc + 'data/modis/truecolor/'
cf_dataloc = dataloc + 'data/modis/cloudfraction_numeric/'

lm_dataloc = dataloc + 'data/validation_dataset/binary_landmask/'
lb_dataloc = dataloc + 'data/validation_dataset/binary_floes/'
lf_dataloc = dataloc + 'data/validation_dataset/binary_landfast/'

masie_ice_loc = dataloc + 'data/masie/seaice/'
masie_land_loc = dataloc + 'data/masie/landmask/'

tc_images = {}
fc_images = {}
cf_images = {}
lb_images = {}
lf_images = {}
lm_images = {}
mi_images = {}
ml_images = {}

missing = []
for row, data in df.iterrows():
    for datadir, imtype, data_dict in zip([tc_dataloc, fc_dataloc, cf_dataloc,
                                           lb_dataloc, lf_dataloc, lm_dataloc,
                                           masie_ice_loc, masie_land_loc],
                                          ['truecolor', 'falsecolor', 'cloudfraction_numeric',
                                           'binary_floes', 'binary_landfast', 'binary_landmask',
                                           'seaice', 'landmask'],
                                          [tc_images, fc_images, cf_images,
                                           lb_images, lf_images, lm_images,
                                           mi_images, ml_images]):
        try:
            if imtype != 'cloudfraction_numeric':
                with rio.open(datadir + fname(df.loc[row,:], imtype)) as im:
                    data_dict[row] = im.read()
            else:
                data_dict[row] = pd.read_csv(datadir + fname(df.loc[row,:], imtype), index_col=0) 
                data_dict[row].index = data_dict[row].index.astype(int)
                data_dict[row].columns = data_dict[row].columns.astype(int)
                
        except:
            if imtype in ['falsecolor', 'cloudfraction_numeric', 'landmask']:
                print('Couldn\'t read', fname(df.loc[row,:], imtype), imtype)
            elif imtype == 'binary_floes':
                if df.loc[row, 'visible_floes'] == 'yes':
                    missing.append(fname(df.loc[row,:], imtype))
            elif imtype == 'binary_landfast':
                if df.loc[row, 'visible_landfast_ice'] == 'yes':
                    missing.append(fname(df.loc[row,:], imtype))
            elif imtype in ['seaice', 'landmask']: # masie images
                missing.append(fname(df.loc[row,:], imtype))
                

# Make the index of the df a unique case label
df.index = [x.case_number + '_' + x.satellite for row, x in df.iterrows()]

df['cloud_fraction_modis'] = np.nan
for case in cf_images:
    df.loc[case, 'cloud_fraction_modis'] = np.mean(cf_images[case]/100)

# Plotting images
# 3 by 2 layout
# Top row: TC, FC, Cloud Fraction
# Bottom row: MASIE, Manual Labels, Manual Labels + Clouds
for case in fc_images:
    region = df.loc[case,'region']
    cf_manual = (df.loc[case, 'cloud_fraction_manual']*100).round(0)
    cf_modis = (df.loc[case, 'cloud_fraction_modis']*100).round(0)
    
    fig, axs = pplt.subplots(ncols=3, nrows=2)

    for ax, data, title in zip(axs[0,:], 
                              [tc_images[case], fc_images[case], cf_images[case]],
                              ['TC Image', 'FC Image', 'Cloud Fraction (%)']):
        if title in ['TC Image', 'FC Image']:
            ax.imshow(reshape_as_image(data))
        elif title == 'Cloud Fraction (%)':
            c = ax.pcolormesh(data.values, vmin=0, vmax=100, N=17, cmap='Blues')
            ax.colorbar(c, label='Cloud Fraction (%)')
            ax.format(urtitle='Manual: {c0}%\nModis: {c1}'.format(c0=cf_manual, c1=cf_modis))
        ax.format(title=title, yreverse=True)
    
    # plot MASIE data
    ax = axs[1, 0]
    masie_ice = mi_images[case].squeeze()
    masie_land = ml_images[case].squeeze()
    ax.imshow(np.ma.masked_array(masie_land, mask=masie_land != 1), c='steelblue')
    ax.imshow(np.ma.masked_array(masie_ice, mask=masie_ice==0), c='w')
    ax.imshow(np.ma.masked_array(masie_land, mask=masie_land != 2), c='gray')
    ax.imshow(np.ma.masked_array(masie_land, mask=masie_land != 4), c='darkgray')
    h = []
    for c in ['steelblue', 'w', 'lightgray', 'darkgray']:
        h.append(ax.plot([],[],m='s', lw=0, c=c, edgecolor='k'))
    ax.legend(h, ['Water', 'Ice', 'Coast', 'Land'], loc='b', ncols=2)
    ax.format(title='MASIE')
    
    
    # plot manual label data on the next two axes
    binary_land = lm_images[case][0,:,:]
    clouds25 = (cf_images[case] >= 25).astype(int)
    clouds50 = (cf_images[case] >= 50).astype(int)
    clouds = (clouds25 + clouds50)/2

    for ax in axs[1, 1:]:
        ax.imshow(np.ma.masked_array(binary_land, mask=binary_land == 0), c='gray9')
        if case in lb_images:
            manual_ice = lb_images[case][0,:,:]
            ax.imshow(np.ma.masked_array(manual_ice, mask=manual_ice==0), c='red5')
        else:
            ax.format(ultitle='No ice mask')
        
        if case in lf_images:
            manual_landfast = lf_images[case][0,:,:]
            ax.imshow(np.ma.masked_array(manual_landfast, mask=manual_landfast == 0), c='yellow4')
        else:
            ax.format(urtitle='No landfast mask')

    # Overlay clouds on last pannel
    ax.imshow(np.ma.masked_array(clouds, mask=clouds==0), cmap='Blues', vmin=0, vmax=2, alpha=0.5)
        
    h = []
    for c in ['red5', 'yellow4', 'darkgray', 'lightblue']:
        h.append(ax.plot([],[],m='s', lw=0, c=c, edgecolor='k'))

    axs[1,1].legend(h[:-1], ['Sea Ice', 'Landfast Ice', 'Land'], loc='b', ncols=2)
    axs[1,1].format(title='Manual Labels')
    
        
    axs[1,2].legend(h, ['Sea Ice', 'Landfast Ice', 'Land', 'MODIS cloud'], loc='b', ncols=2)
    axs[1,2].format(title='Manual Labels + Cloud')
    
    axs.format(yreverse=True, suptitle=case.replace('_', ' ').title() + ' ' + region.replace('_', ' ').title())
    fig.save('../data/validation_dataset/quicklook_images/' + case + '_' + region + '_quicklook.png', dpi=300)
    pplt.close(fig)
