import os
import glob
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.convolution import convolve, Box1DKernel
import lsst.utils
from lsst.daf.butlerUtils import ExposureIdInfo
from lsst.afw.detection import GaussianPsf
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
from lsst.meas.astrom import displayAstrometry
from lsst.meas.algorithms import estimateBackground
from lsst.pipe.tasks.calibrate import DetectAndMeasureTask
from lsst.pipe.tasks.repair import RepairTask
from lsst.ip.isr.isr import darkCorrection


np.random.seed(1)

FilterName = "Ks"

def loadData(psfSigma=1.5):
    """Prepare the data we need to run the example"""
    
    # Load sample input from disk
    imFile = os.path.join(inputdir, "postISRCCD.fits")

    exposure = afwImage.ExposureF(imFile)
    # add a filter
    afwImage.Filter.define(afwImage.FilterProperty(FilterName, 2162))
    exposure.setFilter(afwImage.Filter(FilterName))
    # add a simple Gaussian PSF model
    psfModel = GaussianPsf(11, 11, psfSigma)
    exposure.setPsf(psfModel)
    
    return exposure

def createAndSubtractDark():
    
    imlist = glob.glob(os.path.join(inputdir, "postISRSTD_*.fits"))
    
    darkArray = np.zeros((5, 1024,1024))
    
    k = 0
    
    for fitsfile in imlist:
        
        im = afwImage.ExposureF(fitsfile)
        maskedImage = im.getMaskedImage()
        
        imSlice = np.zeros((32, 1024,32))
                           
        for i in range(0,32):
            lowBound = (i+1)*32 - 32
            upperBound = (i+1)*32
            imSliceMaskedImage = maskedImage[lowBound:upperBound,0:1024]
            imSlice[i] = imSliceMaskedImage.getImage().getArray()

        #imcombine (need to add sigmaclip)
        #imcombine (input="@lista.tmp",output= "dark_"//ima, combine="average", reject="sigclip", lsigma=5, hsigma=5)
        tempDark = np.mean(imSlice, axis=0)
        tempDark = sigma_clip(tempDark, sigma=5, iters=5)

        #  further compress the image on the x axis...
        #  blkavg (input="dark_"//ima,output = "tmp1", b1 = 34, b2=1, option = "average")

        tempDark = tempDark.reshape((tempDark.shape[0], -1, 1))
        tempDark = np.mean(tempDark, axis=1)
        
        # smooth it
        tempDark = tempDark.reshape(-1)
        tempDark = convolve(tempDark, Box1DKernel(10))
        tempDark = tempDark.reshape(1024,1)
        
        tempDark = np.repeat(tempDark, 1024, axis=1)
        
        imArray = maskedImage.getImage().getArray()
        imArray = sigma_clip(imArray, sigma=3, iters=5)
        im_mean = np.mean(imArray)
        darkArray[k] = tempDark - im_mean
        k = k + 1
    
    darkArr = np.mean(darkArray, axis=0)
    print darkArr.shape
    darkArrClipped = np.ma.getdata(sigma_clip(darkArr, sigma=5, iters=5))
        
    hdu = fits.PrimaryHDU(darkArrClipped)
    hdu.writeto('dark.fits')
    
    darkExposure = afwImage.ExposureF(os.path.join(inputdir, "dark.fits"))
    darkMaskedImage = darkExposure.getMaskedImage()
    
    k = 1
    for fitsfile in imlist:
        rawExposure = afwImage.ExposureF(fitsfile)
        maskedImage = rawExposure.getMaskedImage()
        darkCorrection(maskedImage, darkMaskedImage, 1.0, 1.0)
        name = "postISR" + str(k) + ".fits"
        maskedImage.writeFits(name)
        k = k + 1

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
    isrConfig.assembleCcd.doRenorm = False #We'll take care of gain in the flats
    isrConfig.assembleCcd.setGain = False
    SofiIsrTask = SofiIsrTask(config = isrConfig)
    
    darkExposure = afwImage.ExposureF(os.path.join(inputdir, "DARK_1.2.fits.gz"))
    flatExposure = afwImage.ExposureF(os.path.join(inputdir,"Flat06Feb.fits.gz"))
    
    imlist = glob.glob(os.path.join(inputdir,"STD_9104_05feb_01*.fits"))
    k = 1
    for fitsfile in imlist:
        rawExposure = afwImage.ExposureF(fitsfile)
        output = SofiIsrTask.runIsr(rawExposure, dark=darkExposure, flat=flatExposure)
        name = "postISRSTD_" + str(k) + ".fits"
        exposure = output.exposure
        exposure.writeFits(name)
        k = k + 1

def runDam(display=False):
    """Subtract background, mask cosmic rays, then detect and measure
        """
    # Create the tasks; note that background estimation is performed by a function,
    # not a task, though it has a config
    repairConfig = RepairTask.ConfigClass()
    repairTask = RepairTask(config=repairConfig)
    
    backgroundConfig = estimateBackground.ConfigClass()
    
    damConfig = DetectAndMeasureTask.ConfigClass()
    damConfig.detection.thresholdValue = 5.0
    damConfig.detection.includeThresholdMultiplier = 1.0
    damConfig.measurement.doApplyApCorr = "yes"
    detectAndMeasureTask = DetectAndMeasureTask(config=damConfig)
    
    # load the data
    # Exposure ID and the number of bits required for exposure IDs are usually obtained from a data repo,
    # but here we pick reasonable values (there are 64 bits to share between exposure IDs and source IDs).
    exposure = loadData()
    exposureIdInfo = ExposureIdInfo(expId=1, expBits=32)
    
    #repair cosmic rays
    repairTask.run(exposure=exposure)
    
    # subtract an initial estimate of background level
    estBg, exposure = estimateBackground(
        exposure = exposure,
        backgroundConfig = backgroundConfig,
        subtract = True,
        )
        
    # detect and measure
    damRes = detectAndMeasureTask.run(exposure=exposure, exposureIdInfo=exposureIdInfo)
    if display:
        displayAstrometry(frame=2, exposure=damRes.exposure, sourceCat=damRes.sourceCat, pause=False)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of IsrTask and DetectAndMeasureTask")
    parser.add_argument("--inputdir", default=".", help="Input directory")
    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    parser.add_argument('--ds9', action="store_true", help="Display the result?", default=False)
    parser.add_argument('--write', '-w', action="store_true", help="Write the result?", default=False)
    parser.add_argument('--display', action="store_true", help="Display the output image and source catalog", default=False)
    args = parser.parse_args()
    inputdir = args.inputdir
    
    if args.debug:
        try:
            import debug
        except ImportError as e:
            print >> sys.stderr, e

    runIsr()
    
    if args.ds9:
        im = exposure.getMaskedImage().getImage()
        im_median = numpy.median(im.getArray())
        ds9.mtv(im)
        ds9.scale(min=im_median*0.90, max=im_median*1.1, type='SQRT')

    if args.write:
        exposure.writeFits("postISRCCD.fits")

    createAndSubtractDark()
