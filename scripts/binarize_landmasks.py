"""Import the RGB landmasks from the modis folder and convert to binary PNG files"""

import pandas as pd
import numpy as np
import os
import rasterio as rio
import skimage

data_loc = '../data/modis/landmask/'
save_loc = '../data/validation_images/binary_landmask/'
grayscale_threshold = 0.4
files = os.listdir(data_loc)
files = [f for f in files if 'tiff' in f]

for file in files:
    info, satellite, imtype, dx, ftype = file.split('.')
    info = info.replace('-100km', '') # image dims not needed 
    
    new_fname = '-'.join([info, satellite, 'binary_landmask']) + '.png'

    im = rio.open(data_loc + file).read()

    # Convert to grayscale with equal weighting. Coastlines blur into the ocean
    # in the images, so we set a threshold to divide land and water.
    im = np.mean(im[:,:,:], axis=0)/255 # Convert to grayscale
    im = (im > grayscale_threshold).astype(int) * 255
    im = im.astype(np.uint8)

    skimage.io.imsave(save_loc + new_fname, im)
    