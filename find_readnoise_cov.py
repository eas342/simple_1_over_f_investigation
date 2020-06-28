from astropy.io import fits, ascii
import numpy as np
import glob
import os
from copy import deepcopy
import pdb
import os


exampleFile = 'pairwise_sub_red_eachAmpAvg/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.red_eachAmpAvg_int_055_054.fits'

def collect_darks(searchDir = 'pairwise_sub_red_eachAmpAvg',yStart=30,yEnd=35,xStart=4,xEnd=2044,
                  wildCard='NRCN*.fits'):
    fileL = np.sort(glob.glob(os.path.join(searchDir,wildCard)))
    nFile = len(fileL)
    nY = yEnd - yStart
    nX = xEnd - xStart
    dataKeep = np.zeros([nY,nX,nFile])
    for ind,oneFile in enumerate(fileL):
        dat = fits.getdata(oneFile)
        if len(dat.shape) == 3:
            ## for NCDHAS data cubes, grab the first plane for the slope
            useDat = dat[0]
        else:
            useDat = dat
        dataKeep[:,:,ind] = useDat[yStart:yEnd,xStart:xEnd]
    
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

def save_cov_matrix(searchDir = 'pairwise_sub_red',yStart=30,outName=None,
                    wildCard='NRCN*.fits',isSlope=False,gainDivide=False):
    """ Save a typical covariance matrix for a 0.1 um wavelength bin """
    spatialAp = 14
    nCols = 100
    dataKeep = collect_darks(yStart=yStart,yEnd=yStart + spatialAp,searchDir=searchDir,
                             wildCard=wildCard)
    cov_matrix = get_median_cov_matrix(dataKeep,nCols=nCols)
    
    ## get the original file list
    fileL = np.sort(glob.glob(os.path.join(searchDir,wildCard)))
    orig_head = fits.getheader(fileL[0])

    if isSlope == True:
        if ('EFFINTTM' in orig_head):
            itime = orig_head['EFFINTTM']
        elif 'ITIME' in orig_head:
            itime = orig_head['ITIME']
        else:
            print("Couldn't find integration time. Setting to 1.")
            itime = 1
        cov_matrix = cov_matrix * itime**2
    if gainDivide == True:
        cov_matrix = cov_matrix / 1.8**2
    

    primHDU = fits.PrimaryHDU(cov_matrix)
    primHDU.header['NSPATIAL'] = (spatialAp, "number of pixels in the spatial direction (Y)")
    primHDU.header['NDISP'] = (nCols, "Number of pixels in the dispersion direction (X)")
    primHDU.header['PROCTYPE'] = (searchDir, "Processing that was appled")
    primHDU.header['ORIGFILE'] = (orig_head['FILENAME'], 'Original filename')
    primHDU.header['OBS_ID'] = (orig_head['OBS_ID'], 'Observation ID')
    if outName == None:
        outName = os.path.basename(searchDir)
    primHDU.writeto('cov_matrices/read_noise_cov_matrix_spatial_spectral_{}.fits'.format(outName),overwrite=True)
    
def save_cov_eachAmpAvg():
    save_cov_matrix(searchDir='pairwise_sub_red_eachAmpAvg')

def save_cov_colRow():
    save_cov_matrix(searchDir='pairwise_sub_red_rowColSub')

def save_cov_mirage_sim():
    searchPath1 = '/surtrdata1/tso_analysis/sim_data/mirage_sim/nircam_011_proc/proc/raw_separated_MMM_refpix/'
    searchPath2 = 'nircam_011_nircam_011_jw88888001001_01101_00002-seg00?_nrca5_uncal/'
    searchDir = os.path.join(searchPath1,searchPath2)
    
    save_cov_matrix(searchDir=searchDir,yStart=83,outName='mirage_sim_001',wildCard='*.slp.fits',
                    isSlope=True)

def save_cov_mirage_sim(backsub=True):
    if backsub == True:
        searchPath1 = '/surtrdata1/tso_analysis/sim_data/mirage_sim/nircam_011_proc/proc/raw_separated_MMM_refpix/'
        searchPath2 = 'nircam_011_nircam_011_jw88888001001_01101_00002-seg001_nrca5_uncal/backsub_img'
        outName = 'mirage_sim_001_backsub'
        isSlope=False
        gainDivide=True
    else:
        searchPath1 = '/surtrdata1/tso_analysis/sim_data/mirage_sim/nircam_011_proc/proc/raw_separated_MMM_refpix/'
        searchPath2 = 'nircam_011_nircam_011_jw88888001001_01101_00002-seg00?_nrca5_uncal/'
        outName = 'mirage_sim_001'
        isSlope=True
        gainDivide=False
    
    searchDir = os.path.join(searchPath1,searchPath2)
    
    save_cov_matrix(searchDir=searchDir,yStart=83,outName=outName,wildCard='*.slp.fits',
                    isSlope=isSlope,gainDivide=gainDivide)
