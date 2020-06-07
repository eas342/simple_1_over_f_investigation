from astropy.io import fits, ascii
import numpy as np
import glob
import os
from copy import deepcopy
import pdb
import os


exampleFile = 'pairwise_sub_red_eachAmpAvg/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.red_eachAmpAvg_int_055_054.fits'

def collect_darks(searchDir = 'pairwise_sub_red_eachAmpAvg',yStart=30,yEnd=35,xStart=4,xEnd=2044):
    fileL = glob.glob(os.path.join(searchDir,'NRCN*.fits'))
    nFile = len(fileL)
    nY = yEnd - yStart
    nX = xEnd - xStart
    dataKeep = np.zeros([nY,nX,nFile])
    for ind,oneFile in enumerate(fileL):
        dat = fits.getdata(oneFile)
        dataKeep[:,:,ind] = dat[yStart:yEnd,xStart:xEnd]
    
    return dataKeep

def get_median_cov_matrix(dataKeep,nCols=1):
    nX = np.int(np.floor(dataKeep.shape[1]/nCols))
    nY = dataKeep.shape[0]
    nFile = dataKeep.shape[2]
    all_cov = np.zeros([nY * nCols,nY * nCols,nX])
    for oneX in np.arange(nX):
        inData = dataKeep[:,oneX * nCols:oneX * nCols + nCols,:]
        reshaped = np.reshape(inData,[nY * nCols,nFile])
#        reshaped = np.reshape(np.transpose(inData,[2,1,0]),[nY * nCols,nFile])
        all_cov[:,:,oneX] = np.cov(reshaped)
    
    return np.median(all_cov,axis=2)

def show_cov_matrix(yStart=30,yEnd=35,dataKeep=None,nCols=1):
    if dataKeep is None:
        dataKeep = collect_darks(yStart=yStart,yEnd=yEnd)
    cov_matrix = get_median_cov_matrix(dataKeep,nCols=nCols)
    print(np.array_str(cov_matrix,suppress_small=True))

def several_cov_matrices(nCols=1):
    yStarts = np.array([25,30,35,40,45])
    yEnds = yStarts + 5
    for ind in np.arange(len(yStarts)):
        print('Cov centered at {}'.format(np.mean([yStarts[ind],yEnds[ind]])))
        show_cov_matrix(yStart=yStarts[ind],yEnd=yEnds[ind])

def show_spatial_covariance(whiteNoise=False):
    """ Showt a covariance matrix w/ 2 spatial and 2 spectral pixels"""
    if whiteNoise == True:
        dataKeep = np.random.randn(2,2040,54) * np.sqrt(40.)
    else:
        dataKeep = collect_darks(searchDir='pairwise_sub_red',yStart=30,yEnd=32)
    
    show_cov_matrix(yStart=30,yEnd=32,nCols=2,dataKeep=dataKeep)

def save_cov_matrix(searchDir = 'pairwise_sub_red'):
    """ Save a typical covariance matrix for a 0.1 um wavelength bin """
    spatialAp = 14
    nCols = 100
    dataKeep = collect_darks(yStart=30,yEnd=30 + spatialAp,searchDir=searchDir)
    cov_matrix = get_median_cov_matrix(dataKeep,nCols=nCols)
    
    ## get the original file list
    fileL = glob.glob(os.path.join(searchDir,'NRCN*.fits'))
    orig_head = fits.getheader(fileL[0])

    primHDU = fits.PrimaryHDU(cov_matrix)
    primHDU.header['NSPATIAL'] = (spatialAp, "number of pixels in the spatial direction (Y)")
    primHDU.header['NDISP'] = (nCols, "Number of pixels in the dispersion direction (X)")
    primHDU.header['PROCTYPE'] = (searchDir, "Processing that was appled")
    primHDU.header['ORIGFILE'] = (orig_head['FILENAME'], 'Original filename')
    primHDU.header['OBS_ID'] = (orig_head['OBS_ID'], 'Observation ID')
    primHDU.writeto('cov_matrices/read_noise_cov_matrix_spatial_spectral_{}.fits'.format(searchDir),overwrite=True)
    
def save_cov_eachAmpAvg():
    save_cov_matrix(searchDir='pairwise_sub_red_eachAmpAvg')

def save_cov_colRow():
    save_cov_matrix(searchDir='pairwise_sub_red_rowColSub')
