from astropy.io import fits, ascii
import numpy as np
import glob
import os
from copy import deepcopy
import pdb
import os

def add_source(inputDir='pairwise_sub_red',
               sourceToAdd='slope_NRCALONG_FULL_F322W2_sim.fits',custText='',
               slopeFiles=False):
    
    if slopeFiles == True:
        searchName = 'NRC*.slp.fits'
    else:
        searchName = 'NRC*.fits'
    
    fileL = glob.glob(os.path.join(inputDir,searchName))
    
    sourceImgFile = os.path.join('sim_grism_slope',sourceToAdd)
    sourceImg = fits.getdata(sourceImgFile)
    sourceHdr = fits.getheader(sourceImgFile)

    outDir = inputDir + "_grism" + custText
    if os.path.exists(outDir) == False:
        os.mkdir(outDir)
    

    frameTime = 10.73677
    gain = 1.8
    
    if slopeFiles == True:
        sourceToAdd = sourceImg / gain ## convert to ADU
    else:
        sourceToAdd = sourceImg * frameTime / gain ## convert to ADU and add to counts
    
    for ind,oneFile in enumerate(fileL):
        if np.mod(ind,10) == 0:
            print("Working on {} of {}".format(ind+1,len(fileL)))
        HDUList = fits.open(oneFile)
            
        HDUList[0].data = HDUList[0].data + sourceToAdd
        
        outFileName = os.path.basename(oneFile)
        outFilePath = os.path.join(outDir,outFileName)
        HDUList[0].header['FILTER'] = (sourceHdr['FILTER'], sourceHdr.comments['FILTER'])
        HDUList.writeto(outFilePath,overwrite=True)
        HDUList.close()
        
if __name__ == "__main__":
    add_source()
    

def azlab_dark_sim():
    """ Simulate grism image using AZ Lab subarray darks """
    add_source(sourceToAdd='slope_NRCALONG_SUBGRISM256_F322W2_sim.fits',
               inputDir='/fenrirdata1/es_tso/AZLab05/all_azlab05/proc/raw_separated_MMM_refpix/NRCBSELINEDARKNCOLS2048_1_487_S_2019-08-03T07h36m37',
               slopeFiles=True,custText='_f322w2')
