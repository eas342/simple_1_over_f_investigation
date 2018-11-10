from astropy.io import ascii,fits
import os
import numpy as np
import pdb

fullDir = '/data1/External/OTIS_unzipped/NRCNRCALONG-DARK-7235074213_1_1_924_JW1_JLAB40_20170823T074323.606_20170823T080253.912/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.fits'
pairWiseDir = 'pairwise_sub'
overWrite = True

fullDat = fits.open(fullDir)

dat = fullDat['PRIMARY'].data
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
        outPath = os.path.join(pairWiseDir,fitsName)
        HDUList.writeto(outPath,overwrite=overWrite)

    
