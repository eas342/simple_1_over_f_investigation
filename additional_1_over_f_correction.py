from astropy.io import ascii,fits
import os
import numpy as np
import pdb
import sys
from sys import argv
import glob
from copy import deepcopy
from scipy import ndimage

"""
Command line script to apply different forms of 1/f corrections to data
These are intended to go beyond traditional reference pixel correction methods

Usage
python additional_1_over_f_correction.py Method
where Method is 'rowColSub', 'smoothedRowKernel'
refAmpFlip# where # is 0 through 3.
If Method is not specified, it does rowSub or 'row-by-row subtraction'
"""

#redFile = 'proc/NRCNRCA3-DARK-72552035481_1_483_SE_2017-09-12T23h40m37.red.fits'
redFile = 'proc/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.red.fits'

overWrite = True

if len(argv) > 1:
    if argv[1] == 'rowColSub':
        outDir = 'proc_red_additional_rowcol_sub'
        correctionMode = 'rowColSub'
    elif argv[1] == 'smoothedRowKernel':
        windowSize = 40
        outDir = 'proc_red_smoothedRowKernel'
        correctionMode = 'smoothedRowKernel'
        kernel = np.ones([1,windowSize])/np.float(windowSize)
    elif 'refAmpFlip' in argv[1]:
        outDir = 'proc_red_{}'.format(argv[1])
        correctionMode = argv[1]
    else:
        print("unrecognized correction type")
        sys.exit()
else:
    outDir = 'proc_red_additional_rowSub'
    correctionMode = 'rowSub'

if os.path.exists(outDir) == False:
    os.mkdir(outDir)

HDUListOrig = HDUList = fits.open(redFile)
origDat = HDUListOrig[0].data
headOrig = HDUListOrig[0].header

baseName = os.path.splitext(os.path.basename(redFile))[0]
newDat = np.zeros_like(origDat)

for ind,oneGroup in enumerate(origDat):
    dat = np.array(oneGroup)
    head = deepcopy(headOrig)
    
    medianOfEachRow = np.median(dat,axis=1)
    correction2D = np.tile(medianOfEachRow,[head['NAXIS1'],1]).transpose()

    if correctionMode == 'rowColSub':
        medianOfEachColumn = np.median(dat,axis=0)
        correction2D_col = np.tile(medianOfEachColumn,[head['NAXIS1'],1])
        correctedDat = dat - correction2D - correction2D_col
    elif correctionMode == 'smoothedRowKernel':
        smoothedImage = ndimage.convolve(dat,kernel)
        correctedDat = dat - smoothedImage
    elif correctionMode == 'rowSub':
        correctedDat = dat - correction2D
    elif 'refAmpFlip' in correctionMode:
        try:
            refAmp = int(correctionMode[-1])
        except ValueError:
            print('No valid amplifier recognized in {}'.format(correctionMode))
            sys.exit()
        
        refAmpImg = dat[:,refAmp * 512:(refAmp + 1) *512]
        correctionImg = np.zeros_like(dat)
        ampDirections = [1,-1,1,-1]
        nAmp = 4
        for oneAmp in np.arange(nAmp):
            if ampDirections[oneAmp] == ampDirections[refAmp]:
                orientedAmpImage = refAmpImg
            else:
                orientedAmpImage = np.flip(refAmpImg,axis=1)
            correctionImg[:,oneAmp * 512:(oneAmp + 1) * 512] = orientedAmpImage
        correctedDat = dat - correctionImg
    else:
        print("unrecognized correction mode")
        sys.exit()
    
    newDat[ind] = correctedDat

imgHDU = fits.PrimaryHDU(newDat,headOrig)
HDUList = fits.HDUList([imgHDU])
fitsName = "{}_{}.fits".format(baseName,correctionMode)
if correctionMode == 'rowColSub':
    imgHDU.header["CCORRECT"] = (True, "Is a column-by-column median subtraction applied?")
else:
    imgHDU.header["CCORRECT"] = (False, "Is a column-by-column median subtraction applied?")

imgHDU.header["RCORRECT"] = (True, "Is a row-by-row median subtration applied to each group?")

outPath = os.path.join(outDir,fitsName)
HDUList.writeto(outPath,overwrite=overWrite)

HDUListOrig.close()
