from astropy.io import ascii,fits
import os
import numpy as np
import pdb
import sys
from sys import argv
import glob
from copy import deepcopy

redFile = 'proc/NRCNRCA3-DARK-72552035481_1_483_SE_2017-09-12T23h40m37.red.fits'

overWrite = True

if len(argv) > 1:
    if argv[1] == 'rowColSub':
        outDir = 'proc_pairwise_sub_red_additional_rowcol_sub'
        correctionMode = 'rowColSub'
    else:
        print("unrecognized correction type")
        sys.exit()
else:
    outDir = 'proc_pairwise_sub_red_additional_rowSub'
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
    else:
        correctedDat = dat - correction2D
    
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
