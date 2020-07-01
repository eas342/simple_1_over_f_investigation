import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import pdb
from scipy import stats
import os

HDUList, head, dat, testDat = None, None, None, None

def get_data(detector='A3'):
    global HDUList
    global head
    global dat
    global testDat
    if detector == 'A3':
        HDUList = fits.open('proc_red_extra_bias_sub/NRCNRCA3-DARK-72552035481_1_483_SE_2017-09-12T23h40m37.red_extra_bias_sub.fits')
    elif detector == 'ALONG':
        HDUList = fits.open('proc_red_extra_bias_sub/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.red_extra_bias_sub.fits')
    else:
        raise Exception("No red file for detector {}".format(detector))
    
    head = HDUList[0].header
    dat = HDUList[0].data
    testDat = dat[:,0:300,0:200]

get_data()

def find_rtn(testMode=True,algorithm='multi-jumps'):
    if testMode == True:
        useDat = testDat
    else:
        useDat = dat
    
    
    medianAbsDev = np.median(np.abs(useDat - np.median(useDat)))


    if algorithm == 'mean-mode':
        ## round the data into bins slightly smaller than the stdev to find histogram peaks
        roundedDat = np.array(np.round(useDat/3.) * 3,dtype=np.int)
        ## find the median for each pixel
        medianMap = np.median(useDat,axis=0)
        
        #scipy.statsmode is very slow for this big data set
        #so cache the results
        if testMode == True:
            modeMap_name = 'mode_test_data.fits'
        else:
            modeMap_name = 'mode_{}'.format(head['FILENAME'])
        modeMap_path = os.path.join('modes_and_histograms','mode_maps',modeMap_name)
        
        if (os.path.exists(modeMap_path) == False):
            modeMap = stats.mode(roundedDat,axis=0,nan_policy='omit')
            hduModeMap = fits.PrimaryHDU(modeMap)
            hduModeMap.writeto(modeMap_path)
        else:
            modeMap = fits.getdata(modeMap_path)
        
        ## find points where mode and median are significantly differenty 
        modeOffsets = np.abs(medianMap - modeMap[0][0]) > 10.
        
        jumpMap = modeOffsets
        
    else:
        diffDat = np.diff(useDat,axis=0)
        diffDat = diffDat[1:] ## throw out the first diff, since it can be weird
        medianDiff = np.median(diffDat,axis=0)
        nPlanes = useDat.shape[0]
            #    medianAbsDev = np.median(np.abs(useDat - np.tile(np.median(useDat,axis=0),[nPlanes,1,1])),axis=0)
        
        if algorithm == 'big-jumps':
            ## Identify potential jumps from the deltas up the ramp
            potJumps = (np.abs(diffDat) > 100)
            jumpMap = np.sum(potJumps,axis=0)
        #    jumpMaxSizesMap = np.max(np.abs(potJumps * diffDat),axis=0)
        elif algorithm == 'multi-jumps':
            potJumps = (np.abs(diffDat) > 40) & (np.abs(diffDat) < 200)
            jumpCount = np.sum(potJumps,axis=0)
            jumpMap = (jumpCount > 1) & (jumpCount < 4)
        else:
            raise Exception("Unrecognized RTN algorithm")
        
        ## Make sure jumps stand out from rest of ramp (unlike RC say)
        #   outlierJumps = np.greater(jumpMaxSizesMap,medianAbsDev * 20.)
        
    ## Make sure the pixel doesn't have significant RC at the beginning
    earlySlope = np.median(useDat[5:10,:,:],axis=0) - np.median(useDat[0:5,:,:],axis=0)
    lateSlope = np.median(useDat[-10:-5,:,:],axis=0) - np.median(useDat[-5:,:,:],axis=0)
    linearSlopes = np.less(np.abs(lateSlope - earlySlope),20. * medianAbsDev)
    overallSlopes = np.mean(useDat[-25:,:,:],axis=0) - np.mean(useDat[0:25,:,:],axis=0)
    flatSlopes = np.less(overallSlopes,15)
    normalSlopes = flatSlopes & linearSlopes

    ## Find the X & Y locations
    ## look for multiple Jumps
    RTN = np.array(jumpMap & normalSlopes,dtype=np.uint8)
    
    whereJumps = np.where(RTN)
    yJumps, xJumps = whereJumps[0], whereJumps[1]
    
    return yJumps, xJumps, RTN

def find_rc(testMode=True,diffLogDiff=False):
    if testMode == True:
        useDat = testDat
    else:
        useDat = dat
    diffData = np.diff(useDat,axis=0)
    if diffLogDiff == True:
        diffLogDiff = np.diff(np.log10(np.diff(diffData,axis=0) + 65537),axis=0)
        medDiffLogDiff = np.median(diffLogDiff,axis=0)
        rcMap = medDiffLogDiff < -1e-5
    else:
        bigDiff = (diffData > 90)
        numBigDiff = np.sum(bigDiff,axis=0)
        rcMap = np.array(numBigDiff > 6,dtype=np.uint8)
    
    whereRC = np.where(rcMap)
    yRC, xRC = whereRC[0], whereRC[1]
    return yRC, xRC, rcMap

def save_map(thisMap,mapType='rtn'):
    primHDU = fits.PrimaryHDU(thisMap)
    primHDU.header['DETECTOR'] = head['DETECTOR']
    primHDU.writeto('maps/{}_map_{}.fits'.format(mapType,head['DETECTOR']),
                    overwrite=True)

def plot_ramps(yJumps,xJumps,plotType='rtn',startIndex=0,numRamps=10):
    fig, ax = plt.subplots()
    for ind,oneJump in enumerate(np.arange(numRamps) + startIndex):
        if plotType == 'rtn':
            offset = ind + 50 * ind
        else:
            offset = 0
        
        ax.plot(dat[:,yJumps[oneJump],xJumps[oneJump]] - offset)
        
    
    ax.set_xlabel('Frame Number')
    ax.set_ylabel('$\Delta$ DN')
    fig.savefig('pixel_ramps/example_{}_{}.pdf'.format(plotType,head['DETECTOR']),
                bbox_inches='tight')
    plt.close(fig)

    
def do_all():
    for oneDetector in ['A3','ALONG']:
        get_data(oneDetector)
        for pixType in ['rtn','rc']:
            if pixType == 'rtn':
                yPx, xPx, thisMap = find_rtn(testMode=False)
                numRamps=5
            else:
                yPx, xPx, thisMap = find_rc(testMode=False)
                numRamps=5
            
            print("Found {} {} pixels".format(len(yPx),pixType))
            plot_ramps(yPx,xPx,plotType=pixType,numRamps=numRamps)
            save_map(thisMap,mapType=pixType)
    
if __name__ == "__main__":
    do_all()
