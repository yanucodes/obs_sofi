import os
import glob
import numpy as np
import lsst.afw.image as afwImage
from astropy.io import fits
import matplotlib.pyplot as plt

def blkavg(arr,x1,x2,y1,y2):
    arr = arr[x1:x2,y1:y2]
    
    plt.imshow(arr)
    plt.show()
    
    arr = arr.reshape((arr.shape[0], -1, 1))
    arr = np.mean(arr, axis=1)

    return arr


def createFlat(flist):
    
    tempFlats = []
    
    for i in range(4):
        temp = []
        tmp1 = afwImage.ExposureF(flist[0+i])
        tmp2 = afwImage.ExposureF(flist[7-i])
        temp.append(tmp1.getMaskedImage().getImage().getArray())
        temp.append(tmp2.getMaskedImage().getImage().getArray())

        tempFlats.append(np.mean(temp, axis = 0))

    tempOff = tempFlats[0]
    tempOffMask = tempFlats[1]
    tempOnMask = tempFlats[2]
    tempOn = tempFlats[3]
    
    plt.imshow(tempOff)
    plt.show()

    tempOnA = blkavg(tempOn, 0, 1024, 500, 600)
    tempOnC = blkavg(tempOnMask, 0, 1024, 500, 600)
    tempOnB = blkavg(tempOnMask, 0, 1024, 50, 150)

    tempOnAC = tempOnA - tempOnC
    tempOnACB = tempOnAC + tempOnB
    
    print tempOnACB.shape

    tempOn2D = np.repeat(tempOnACB, 1024, axis=1)

    tempOnBias = tempOn - tempOn2D

    tempOffA = blkavg(tempOff, 0, 1024, 500, 600)
    tempOffC = blkavg(tempOffMask, 0, 1024, 500, 600)
    tempOffB = blkavg(tempOffMask, 0, 1024, 50, 150)

    tempOffAC = tempOffA - tempOffC
    tempOffACB = tempOffAC + tempOffB

    tempOff2D = np.repeat(tempOffACB, 1024, axis=1)

    tempOffBias = tempOff - tempOff2D

    flat = tempOnBias - tempOffBias

    norm = np.mean(flat)

    flat = flat/norm

    return flat



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Create a special flat for SOFI")
    parser.add_argument("--inputdir", default=".", help="Input directory")
    args = parser.parse_args()
    
    inputdir = args.inputdir

    flist = glob.glob(os.path.join(inputdir, "FLAT_06Feb*.fits"))

    specialFlat = createFlat(flist)

    hdu = fits.PrimaryHDU(specialFlat)

    hdu.writeto("flat.fits")
