"""Use image-by-image thresholds to create"""

import pandas as pd
import numpy as np
import os
import rasterio as rio
import skimage

data_loc = '../data/modis/truecolor/'
save_loc = '../data/validation_dataset/binary_watermask_init/'

# replace with list of adapted images
grayscale_threshold = 0.4


files = os.listdir(data_loc)
files = [f for f in files if 'tiff' in f]

for file in files:
    info, satellite, imtype, dx, ftype = file.split('.')
    info = info.replace('-100km', '') # image dims not needed 
    
    new_fname = '-'.join([info, satellite, 'binary_water_init']) + '.png'

    im = rio.open(data_loc + file).read()

    im = np.mean(im[:,:,:], axis=0)/255 # Convert to grayscale
    im = (im > grayscale_threshold).astype(int) * 255
    im = im.astype(np.uint8)

    skimage.io.imsave(save_loc + new_fname, im)
    
