from astropy.io import ascii,fits
import os
import numpy as np
import pdb
import sys
from sys import argv
from copy import deepcopy

if len(argv) >= 2:
    if argv[1] == 'A3':
        fullDir = '/data1/tso_analysis/otis_long_dark/proc/NRCNRCA3-DARK-72552035481_1_483_SE_2017-09-12T23h40m37.fits'
        pairWiseDir = 'pairwise_sub_A3'
    elif argv[1] == 'redNoSide':
        print("Making subtrative pairs from reduced data w/ no side reference correction")
        fullDir = '/data1/tso_analysis/otis_long_dark/proc2/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.red.fits'
        pairWiseDir = 'pairwise_sub_red_no_side_refcor'
    elif argv[1] == 'red':
        print("Making subtractive pairs from reduced data...")
        fullDir = '/data1/tso_analysis/otis_long_dark/proc/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.red.fits'
        pairWiseDir = 'pairwise_sub_red'
    elif argv[1] == 'redSubThenPairwise':
        print("Making subtractive pairs from reduced median row-subtracted data...")
        fullDir = 'proc_pairwise_sub_red_additional_rowSub/NRCNRCA3-DARK-72552035481_1_483_SE_2017-09-12T23h40m37.red_rowSub.fits'
        pairWiseDir = 'pairwise_sub_red_rowSub'
    else:
        print("Unrecognized argument. Returing ... ")
        sys.exit()
else:
    fullDir = '/data1/External/OTIS_unzipped/NRCNRCALONG-DARK-7235074213_1_1_924_JW1_JLAB40_20170823T074323.606_20170823T080253.912/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.fits'
    pairWiseDir = 'pairwise_sub'

overWrite = True

fullDat = fits.open(fullDir)

dat = np.array(fullDat['PRIMARY'].data,dtype='int32')
head = fullDat['PRIMARY'].header

baseName = os.path.splitext(os.path.basename(fullDir))[0]



if os.path.exists(pairWiseDir) == False:
    os.mkdir(pairWiseDir)

for oneGroup in np.arange(head['NGROUP'] - 1):
    if np.mod(oneGroup,2) == 0:
        diff = dat[oneGroup+1,:,:] - dat[oneGroup,:,:]
        imgHDU = fits.PrimaryHDU(diff,head)
        HDUList = fits.HDUList([imgHDU])
        fitsName = baseName + "_int_{:03d}_{:03d}.fits".format(oneGroup+1,oneGroup)
        imgHDU.header["ON_NINT"] = (oneGroup, "Which integration is the subtrahend")
        imgHDU.header["INT_NEXT"] = (oneGroup+1, "Which integration is the minuend")
        imgHDU.header["INTTIME"] = (head['TFRAME'], "Total integration time for one fudged MULTIACCUM pairwise subtracted")
        outPath = os.path.join(pairWiseDir,fitsName)
        HDUList.writeto(outPath,overwrite=overWrite)

    
