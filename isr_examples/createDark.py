import os
import glob
import numpy as np
import lsst.afw.image as afwImage
from astropy.io import fits

def createDark(flist):
    
    darks = []
    
    for i in range(len(flist)):
        darkExposure = afwImage.ExposureF(flist[i])
        print "exp", i
        darks.append(darkExposure.getMaskedImage().getImage().getArray())
    
    print darks



    dark = np.median(darks, axis = 0)

    return dark



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Create a special flat for SOFI")
    parser.add_argument("--inputdir", default=".", help="Input directory")
    parser.add_argument("--outputdir", default=".", help="Output directory")
    parser.add_argument("--prefix", default="D_10", help="Prefix")
    args = parser.parse_args()
    
    inputdir = args.inputdir
    outputdir = args.outputdir
    filename = args.prefix +"*.fits.gz"

    flist = glob.glob(os.path.join(inputdir, filename))
    
    olddata, newheader = fits.getdata(flist[0], header=True)
    newheader.remove('ESO DET CHIP PXSPACE')

    dark = createDark(flist)

    hdu = fits.PrimaryHDU(dark)
    hdu.header = newheader
    
    fn = args.prefix+".fits"

    hdu.writeto(os.path.join(outputdir, fn), clobber=True)
