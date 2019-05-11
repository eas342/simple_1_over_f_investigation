from astropy.io import fits, ascii

from astropy.io import ascii,fits
import os
import numpy as np
import pdb
import sys
from sys import argv
import glob
from copy import deepcopy
from scipy import ndimage

#redFile = 'proc/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.red.fits'
redFile = 'proc_red_smoothedRowKernel/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.red_smoothedRowKernel.fits'

#grpDir = 'grp_split_red_file'
grpDir = 'grp_split_red_smoothed_RowKernel_file'

overwrite= True

HDUListOrig = HDUList = fits.open(redFile)
origDat = HDUListOrig[0].data
headOrig = HDUListOrig[0].header

baseName = os.path.splitext(os.path.basename(redFile))[0]

if os.path.exists(grpDir) == False:
    os.mkdir(grpDir)

for oneGroup in np.arange(headOrig['NGROUP']):
    oneGroupDat = origDat[oneGroup]
    fitsName = baseName + "_grp_{:04d}.fits".format(oneGroup)
    imgHDU = fits.PrimaryHDU(oneGroupDat,headOrig)
    
    imgHDU.header["ON_GRP"] = (oneGroup, "Which group this is a slice of")
    outPath = os.path.join(grpDir,fitsName)
    imgHDU.writeto(outPath,overwrite=True)
    
