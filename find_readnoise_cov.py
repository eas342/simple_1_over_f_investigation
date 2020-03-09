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
    dataKeep = np.zeros([nY,nFile,nX])
    for ind,oneFile in enumerate(fileL):
        dat = fits.getdata(oneFile)
        dataKeep[:,ind,:] = dat[yStart:yEnd,xStart:xEnd]
    
    return dataKeep

def get_median_cov_matrix(dataKeep):
    nX = dataKeep.shape[2]
    nY = dataKeep.shape[0]
    nFile = dataKeep.shape[1]
    all_cov = np.zeros([nY,nY,nX])
    for oneX in np.arange(nX):
        all_cov[:,:,oneX] = np.cov(dataKeep[:,:,oneX])
    
    return np.median(all_cov,axis=2)

def show_cov_matrix(yStart=30,yEnd=35):
    dataKeep = collect_darks(yStart=yStart,yEnd=yEnd)
    cov_matrix = get_median_cov_matrix(dataKeep)
    print(np.array_str(cov_matrix,suppress_small=True))

def several_cov_matrices():
    yStarts = np.array([25,30,35,40,45])
    yEnds = yStarts + 5
    for ind in np.arange(len(yStarts)):
        print('Cov centered at {}'.format(np.mean([yStarts[ind],yEnds[ind]])))
        show_cov_matrix(yStart=yStarts[ind],yEnd=yEnds[ind])

