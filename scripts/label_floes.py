# simple script: save the labeled images
import skimage as sk
import os

binary_dataloc = "../data/validation_dataset/binary_floes/"
labeled_dataloc = "../data/validation_dataset/labeled_floes/"

files = [f for f in os.listdir(binary_dataloc) if 'png' in f]
for file in files:
    im = sk.io.imread(binary_dataloc + file)[:,:,0]
    labeled_dataloc = binary_dataloc.replace("binary", "labeled")
    labeled_im = sk.measure.label(im).astype(float)
    sk.io.imsave(labeled_dataloc + file.replace("binary", "labeled").replace("png", "tiff"), 
             labeled_im, check_contrast=False)