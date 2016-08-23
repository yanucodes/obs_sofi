import os
import glob
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.convolution import convolve, Box1DKernel
import lsst.utils
from lsst.daf.butlerUtils import ExposureIdInfo
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
from lsst.ip.isr.isr import darkCorrection

def createAndSubtractDark():
    
    imlist = glob.glob(os.path.join(outputdir, "df_F02_*.fits"))
    
    darkArray = np.zeros((5, 1024,1024))
    
    kbeg = 0
    
    for k in range(len(imlist)):
    
        for j in range(kbeg,kbeg+5):
            im = afwImage.ExposureF(imlist[j])
            maskedImage = im.getMaskedImage()
    
            imSlice = np.zeros((32, 1024,32))
                           
            for i in range(0,32):
                lowBound = (i+1)*32 - 32
                upperBound = (i+1)*32
                imSliceMaskedImage = maskedImage[lowBound:upperBound,0:1024]
                imSlice[i] = imSliceMaskedImage.getImage().getArray()
                    
            #imcombine (need to add sigmaclip)
            #imcombine (input="@lista.tmp",output= "dark_"//ima, combine="average", reject="sigclip", lsigma=5, hsigma=5)
            tempDark = np.median(imSlice, axis=0)
            tempDark = sigma_clip(tempDark, sigma=5, iters=5)

            #  further compress the image on the x axis...
            #  blkavg (input="dark_"//ima,output = "tmp1", b1 = 34, b2=1, option = "average")

            tempDark = tempDark.reshape((tempDark.shape[0], -1, 1))
            tempDark = np.median(tempDark, axis=1)
        
            # smooth it
            tempDark = tempDark.reshape(-1)
            tempDark = convolve(tempDark, Box1DKernel(10))
            tempDark = tempDark.reshape(1024,1)
        
            tempDark = np.repeat(tempDark, 1024, axis=1)
        
            imArray = maskedImage.getImage().getArray()
            imArray = sigma_clip(imArray, sigma=3, iters=5)
            im_median = np.median(imArray)
            darkArray[j-kbeg] = tempDark - im_median
    
        darkArr = np.median(darkArray, axis=0)
        darkArrClipped = np.ma.getdata(sigma_clip(darkArr, sigma=5, iters=5))
    
        hdu = fits.PrimaryHDU(darkArrClipped)
        hdu.writeto(os.path.join(calibdir, "dark.fits"), clobber = True)
    
        darkExposure = afwImage.ExposureF(os.path.join(calibdir, "dark.fits"))
        darkMaskedImage = darkExposure.getMaskedImage()
        
        
        rawExposure = afwImage.ExposureF(imlist[k])
        maskedImage = rawExposure.getMaskedImage()
        darkCorrection(maskedImage, darkMaskedImage, 1.0, 1.0)
        
        fn = imlist[k]
        name = "dd_" + str(fn[len(outputdir)+3:len(fn)])
        
        maskedImage.writeFits(os.path.join(outputdir,name))

        if (j >= 5 and j < (len(imlist)-4)):
            kbeg += 1


def runIsr():
    '''Run the task to do ISR on a ccd'''
    from lsst.obs.sofi.runIsrTask import SofiIsrTask
    
    #Create the isr task with modified config
    isrConfig = SofiIsrTask.ConfigClass()
    isrConfig.doBias = False
    isrConfig.doOverscan = False
    isrConfig.overscanFitType = 'POLY'
    isrConfig.doDark = True
    isrConfig.doFlat = True
    isrConfig.doAssembleCcd = False
    isrConfig.doFringe = False
    isrConfig.doLinearize = False
    isrConfig.assembleCcd.doRenorm = False #We'll take care of gain in the flats
    isrConfig.assembleCcd.setGain = False
    sofiIsrTask = SofiIsrTask(config = isrConfig)
    
    darkExposure = afwImage.ExposureF(os.path.join(calibdir, darkfn))
    flatExposure = afwImage.ExposureF(os.path.join(calibdir, flatfn))
    
    imlist = glob.glob(os.path.join(inputdir,"F02_*.fits.gz"))
    
    for fn in imlist:
        rawExposure = afwImage.ExposureF(fn)
        output = sofiIsrTask.runIsr(rawExposure, dark=darkExposure, flat=flatExposure)
        name = "df_" + str(fn[len(inputdir):len(fn)-3])
        exposure = output.exposure
        exposure.writeFits(os.path.join(outputdir, name))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of IsrTask and DetectAndMeasureTask")
    parser.add_argument('--inputdir', default=".", help="Input directory")
    parser.add_argument('--outputdir', default=".", help="Output directory")
    parser.add_argument('--calibdir', default=".", help="Calib directory")
    parser.add_argument('--flat', default="FLAT", help="Flat exposure name")
    parser.add_argument('--dark', default="DARK", help="Dark exposure name")
    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    parser.add_argument('--display', action="store_true", help="Display the output image and source catalog", default=False)
    args = parser.parse_args()
    inputdir = args.inputdir
    outputdir = args.outputdir
    calibdir = args.calibdir
    flatfn = args.flat
    darkfn = args.dark
    
    if args.debug:
        try:
            import debug
        except ImportError as e:
            print >> sys.stderr, e

    runIsr()

    createAndSubtractDark()
