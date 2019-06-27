from astropy.io import ascii,fits
import os
import numpy as np
import pdb
import sys
from sys import argv
import glob
from copy import deepcopy
from scipy import ndimage
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from copy import deepcopy
from astropy.convolution import convolve as astropyConvolve

"""
Command line script to apply different forms of 1/f corrections to data
These are intended to go beyond traditional reference pixel correction methods

Usage
python additional_1_over_f_correction.py Method
where Method is 'rowColSub', 'smoothedRowKernel'
refAmpFlip# where # is 0 through 3.
'pcaEach' - New PCA on each image
If Method is not specified, it does rowSub or 'row-by-row subtraction'
"""

#redFile = 'proc/NRCNRCA3-DARK-72552035481_1_483_SE_2017-09-12T23h40m37.red.fits'
redFile = 'proc/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.red.fits'

overWrite = True


def do_pca_model(inputArr,nComp=5):
    scaler = StandardScaler(with_std=False,with_mean=True)
    inputX = scaler.fit_transform(np.array(inputArr,dtype=np.float))
    pca = PCA(n_components=nComp)
    principalComponents = pca.fit_transform(inputX)
    modelPCA_unscaled = np.dot(principalComponents,pca.components_)
    modelPCA_2D = scaler.inverse_transform(modelPCA_unscaled)
    return modelPCA_2D

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
    elif ('pcaEach' in argv[1]) | ('pcaInd' in argv[1]):
        outDir = 'proc_red_{}'.format(argv[1])
        correctionMode = argv[1]
    elif 'rowKernelInterp' in argv[1]:
        windowSize = 161 # window for kernel bigger than source ap
        outDir = 'proc_red_rowKernelInterp'
        correctionMode = argv[1]
        kernel = np.ones([1,windowSize])/np.float(windowSize)
        
        ### Make a mask for apertures and bad pixels
        ## bad pixel file
        headRed = fits.getheader(redFile)
        if headRed['DETECTOR'] == 'NRCALONG':
            bpmFile = '/usr/local/nircamsuite/cal/BPM/NRCA5_17158_BPM_ISIMCV3_2016-01-21.fits'
        else:
            raise Exception("No BPM for {}".format(headRed['DETECTOR']))
        
        bpm = fits.getdata(bpmFile)
        mask = bpm > 0
        
        ## apertures
        aps = [[  768, 540],[ 1280, 540],[  768, 1250]]
        x, y = np.meshgrid(np.arange(bpm.shape[0]),np.arange(bpm.shape[0]))
        for oneAp in  aps:
            r = np.sqrt((x - oneAp[0])**2 + (y - oneAp[1])**2)
            mask = mask | (r < 70)
        
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
    elif correctionMode == 'rowKernelInterp':
        smoothedImage = astropyConvolve(dat,kernel,mask=mask)
        correctedDat = dat - smoothedImage
    elif correctionMode == 'rowSub':
        correctedDat = dat - correction2D
    elif ('pcaEach' in correctionMode) | ('pcaInd' in correctionMode):
        badP = np.abs(dat) > 200.
        cleanDat = deepcopy(dat)
        cleanDat[badP] = 0
        if correctionMode == 'pcaEach':
            nComp = 10
        elif correctionMode == 'pcaInd':
            nComp = 5
        else:
            nComp = int(correctionMode.split('pcaEach')[-1])

        if 'pcaEach' in correctionMode:
            modelPCA_2D = do_pca_model(cleanDat,nComp=nComp)
        elif 'pcaInd' in correctionMode:
            modelPCA_2D = np.zeros_like(cleanDat)
            for oneAmp in np.arange(4):
                xStart, xEnd = int(oneAmp) * 512, (int(oneAmp) + 1) * 512
                modelPCA_2D[:,xStart:xEnd] = do_pca_model(cleanDat[:,xStart:xEnd],nComp=nComp)
        else:
            Exception("PCA correction mode {} not understdood".format(correctionMode))
            
        
        correctedDat = dat - modelPCA_2D
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
