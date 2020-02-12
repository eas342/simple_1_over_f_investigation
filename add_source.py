from astropy.io import fits, ascii
import numpy as np
import glob
import os
from copy import deepcopy
import pdb

def add_source(inputDir='pairwise_sub_red',
               sourceToAdd='slope_NRCALONG_FULL_F322W2_sim.fits'):
    fileL = glob.glob(os.path.join(inputDir,'NRC*.fits'))
    
    sourceImgFile = os.path.join('sim_grism_slope',sourceToAdd)
    sourceImg = fits.getdata(sourceImgFile)
    sourceHdr = fits.getheader(sourceImgFile)

    outDir = inputDir + "_grism"

    frameTime = 10.73677
    gain = 1.8
    
    for ind,oneFile in enumerate(fileL):
        if np.mod(ind,10) == 0:
            print("Working on {} of {}".format(ind+1,len(fileL)))
        HDUList = fits.open(oneFile)
        
        HDUList[0].data = HDUList[0].data + sourceImg * frameTime / gain ## convert to ADU
        
        outFileName = os.path.basename(oneFile)
        outFilePath = os.path.join(outDir,outFileName)
        HDUList[0].header['FILTER'] = (sourceHdr['FILTER'], sourceHdr.comments['FILTER'])
        HDUList.writeto(outFilePath,overwrite=True)
        HDUList.close()
        
if __name__ == "__main__":
    add_source()
    

