import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits, ascii
from astropy.table import Table
import glob
import os
import pdb
from tser_tools import phot_pipeline
import pixeltime
from astropy.stats import LombScargle
from scipy.stats import binned_statistic

filePath = 'grp_split_red_extra_bias_sub_file'
grpFileName = 'NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.red_extra_bias_sub_grp_{:04d}.fits'
freqGrid = np.linspace(1./5.,1./2e-5,65536)

def get_files():
    fileList = glob.glob(os.path.join(filePath,grpFileName.replace("{:04d}","*")))
    return len(fileList), fileList
    

def get_pixeltime_series(grp=53):
    """ 
    Get a time series by pixel from an example group
    """
    HDUList = fits.open(os.path.join(filePath,grpFileName.format(grp)))
    oneImg = HDUList[0].data
    outT = pixeltime.main.all_pixels_tser(oneImg) 
    HDUList.close()
    outT['t (sec)'] = outT['time'] * 1e-5
    
    return outT

def get_periodogram(t,amp=0,clip=80):

    time = t['t (sec)']
    y = t['Amp {} Val'.format(amp)]
    pts = np.abs(y) < clip
    
    psd = LombScargle(time[pts],y[pts]).power(frequency=freqGrid)
    return psd

def example_psd():
    t = get_pixeltime_series()
    psd = get_periodogram(t)
    fig, ax = plt.subplots()
    ax.loglog(freqGrid,psd)
    fig.savefig('power_spectra/plots/ex_psd.pdf')
    plt.close(fig)

def all_psd_one_amp(amp=0,testMode=False):
    ngroups, fileList = get_files()
    if testMode == True:
        ngroups = 5
        fileList = fileList[0:5]
        update = 1
    else:
        update = 10
    
    nFreq = len(freqGrid)
    allPSD = np.zeros([ngroups,nFreq])
    groupList = np.arange(ngroups)
    for oneGrp in groupList:
        if np.mod(oneGrp,update) == 0:
            print("Working on {} of {}".format(oneGrp,ngroups - 1))
        
        t = get_pixeltime_series(grp=oneGrp)
        allPSD[oneGrp,:] = get_periodogram(t,amp=amp)
    
    primHDU = fits.PrimaryHDU(allPSD)
    primHDU.name = 'POWER'
    
    freqHDU = fits.ImageHDU(freqGrid)
    freqHDU.name = 'FREQ'
    freqHDU.header['UNIT'] = ('freq (HZ)', 'Frequency Units')
    
    grpHDU = fits.ImageHDU(groupList)
    grpHDU.name = 'GROUPS'
    grpHDU.header['UNIT'] = ('Group No','Group number up the ramp')
    
    fileTable = Table()
    
    fileTable['Full_Path'] = fileList
    fileTable['Group_Num'] = groupList
    fileHDU = fits.BinTableHDU(fileTable)
    fileHDU.name = 'FILENAMES'
    
    HDUList = fits.HDUList([primHDU,freqHDU,grpHDU,fileHDU])
    HDUList.writeto('power_spectra/by_group/all_power_spec_amp_{}.fits'.format(amp),
                    overwrite=True)
    
    
def all_psd():
    for oneAmp in np.arange(4):
        print("Workin on Amp {} of {}".format(oneAmp+1,4))
        all_psd_one_amp(oneAmp)

if __name__ == "__main__":
    all_psd()
