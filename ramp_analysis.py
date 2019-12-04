import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import pdb


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
    testDat = dat[:,0:200,0:100]

get_data()

def find_rtn(testMode=True):
    if testMode == True:
        useDat = testDat
    else:
        useDat = dat
    
    diffDat = np.diff(useDat,axis=0)
    diffDat = diffDat[1:] ## throw out the first diff, since it can be weird
    medianDiff = np.median(diffDat,axis=0)
    nPlanes = useDat.shape[0]
    medianAbsDev = np.median(np.abs(useDat - np.median(useDat)))
#    medianAbsDev = np.median(np.abs(useDat - np.tile(np.median(useDat,axis=0),[nPlanes,1,1])),axis=0)
    
    ## Identify potential jumps from the deltas up the ramp
    potJumps = (np.abs(diffDat) > 200.)
    jumpMap = np.sum(potJumps,axis=0)
#    jumpMaxSizesMap = np.max(np.abs(potJumps * diffDat),axis=0)

    ## Make sure jumps stand out from rest of ramp (unlike RC say)
 #   outlierJumps = np.greater(jumpMaxSizesMap,medianAbsDev * 20.)

    ## Make sure the pixel doesn't have significant RC at the beginning
    earlySlope = np.median(useDat[5:10,:,:],axis=0) - np.median(useDat[0:5,:,:],axis=0)
    lateSlope = np.median(useDat[-10:-5,:,:],axis=0) - np.median(useDat[-5:,:,:],axis=0)
    normalSlopes = np.less(np.abs(lateSlope - earlySlope),20. * medianAbsDev)

    ## Find the X & Y locations
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
    for oneJump in np.arange(numRamps) + startIndex:
        ax.plot(dat[:,yJumps[oneJump],xJumps[oneJump]])
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
                numRamps=10
            else:
                yPx, xPx, thisMap = find_rc(testMode=False)
                numRamps=5
            
            plot_ramps(yPx,xPx,plotType=pixType,numRamps=numRamps)
            save_map(thisMap,mapType=pixType)
    
if __name__ == "__main__":
    do_all()
