import os
import numpy as np
import lsst.utils
from lsst.daf.butlerUtils import ExposureIdInfo
from lsst.afw.detection import GaussianPsf
import lsst.afw.image as afwImage
from lsst.meas.astrom import displayAstrometry
from lsst.meas.algorithms import estimateBackground
from lsst.pipe.tasks.calibrate import DetectAndMeasureTask
from lsst.pipe.tasks.repair import RepairTask

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

def createDark():
    imlist = glob.iglob(os.path.join(inputdir, "postISRCCD*.fits"))
    for fits in imlist:
        im = afwImage.Exposure(fits)
        maskedImage = im.getMaskedImage()
        bbox1 = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(32, 1024))
        imslice = {}
        imslice[1] = maskedImage.Factory(maskedImage, bbox1)
        bbox2 = afwGeom.Box2I(afwGeom.Point2I(32, 0), afwGeom.Extent2I(64, 1024))
        imslice[2] = maskedImage.Factory(maskedImage, bbox2)
        bbox3 = afwGeom.Box2I(afwGeom.Point2I(64, 0), afwGeom.Extent2I(96, 1024))
        imslice[3] = maskedImage.Factory(maskedImage, bbox3)
        bbox4 = afwGeom.Box2I(afwGeom.Point2I(96, 0), afwGeom.Extent2I(128, 1024))
        imslice[4] = maskedImage.Factory(maskedImage, bbox4)
        bbox5 = afwGeom.Box2I(afwGeom.Point2I(128, 0), afwGeom.Extent2I(160, 1024))
        imslice[5] = maskedImage.Factory(maskedImage, bbox5)
        bbox6 = afwGeom.Box2I(afwGeom.Point2I(160, 0), afwGeom.Extent2I(192, 1024))
        imslice[6] = maskedImage.Factory(maskedImage, bbox6)
        bbox7 = afwGeom.Box2I(afwGeom.Point2I(192, 0), afwGeom.Extent2I(224, 1024))
        imslice[7] = maskedImage.Factory(maskedImage, bbox7)
        bbox8 = afwGeom.Box2I(afwGeom.Point2I(224, 0), afwGeom.Extent2I(256, 1024))
        imslice[8] = maskedImage.Factory(maskedImage, bbox8)
        bbox9 = afwGeom.Box2I(afwGeom.Point2I(256, 0), afwGeom.Extent2I(288, 1024))
        imslice[9] = maskedImage.Factory(maskedImage, bbox9)
        bbox10 = afwGeom.Box2I(afwGeom.Point2I(288, 0), afwGeom.Extent2I(320, 1024))
        imslice[10] = maskedImage.Factory(maskedImage, bbox10)
        bbox11 = afwGeom.Box2I(afwGeom.Point2I(320, 0), afwGeom.Extent2I(352, 1024))
        imslice[11] = maskedImage.Factory(maskedImage, bbox11)
        bbox12 = afwGeom.Box2I(afwGeom.Point2I(352, 0), afwGeom.Extent2I(384, 1024))
        imslice[12] = maskedImage.Factory(maskedImage, bbox12)
        bbox13 = afwGeom.Box2I(afwGeom.Point2I(384, 0), afwGeom.Extent2I(416, 1024))
        imslice[13] = maskedImage.Factory(maskedImage, bbox13)
        bbox14 = afwGeom.Box2I(afwGeom.Point2I(416, 0), afwGeom.Extent2I(448, 1024))
        imslice[14] = maskedImage.Factory(maskedImage, bbox14)
        bbox15 = afwGeom.Box2I(afwGeom.Point2I(448, 0), afwGeom.Extent2I(480, 1024))
        imslice[15] = maskedImage.Factory(maskedImage, bbox15)
        bbox16 = afwGeom.Box2I(afwGeom.Point2I(480, 0), afwGeom.Extent2I(512, 1024))
        imslice[16] = maskedImage.Factory(maskedImage, bbox16)
        bbox17 = afwGeom.Box2I(afwGeom.Point2I(512, 0), afwGeom.Extent2I(544, 1024))
        imslice[17] = maskedImage.Factory(maskedImage, bbox17)
        bbox18 = afwGeom.Box2I(afwGeom.Point2I(544, 0), afwGeom.Extent2I(576, 1024))
        imslice[18] = maskedImage.Factory(maskedImage, bbox18)
        bbox19 = afwGeom.Box2I(afwGeom.Point2I(576, 0), afwGeom.Extent2I(608, 1024))
        imslice[19] = maskedImage.Factory(maskedImage, bbox19)
        bbox20 = afwGeom.Box2I(afwGeom.Point2I(608, 0), afwGeom.Extent2I(640, 1024))
        imslice[20] = maskedImage.Factory(maskedImage, bbox20)
        bbox21 = afwGeom.Box2I(afwGeom.Point2I(640, 0), afwGeom.Extent2I(672, 1024))
        imslice[21] = maskedImage.Factory(maskedImage, bbox21)
        bbox22 = afwGeom.Box2I(afwGeom.Point2I(672, 0), afwGeom.Extent2I(704, 1024))
        imslice[22] = maskedImage.Factory(maskedImage, bbox22)
        bbox23 = afwGeom.Box2I(afwGeom.Point2I(704, 0), afwGeom.Extent2I(736, 1024))
        imslice[23] = maskedImage.Factory(maskedImage, bbox23)
        bbox24 = afwGeom.Box2I(afwGeom.Point2I(736, 0), afwGeom.Extent2I(768, 1024))
        imslice[24] = maskedImage.Factory(maskedImage, bbox24)
        bbox25 = afwGeom.Box2I(afwGeom.Point2I(768, 0), afwGeom.Extent2I(800, 1024))
        imslice[25] = maskedImage.Factory(maskedImage, bbox25)
        bbox26 = afwGeom.Box2I(afwGeom.Point2I(800, 0), afwGeom.Extent2I(832, 1024))
        imslice[26] = maskedImage.Factory(maskedImage, bbox26)
        bbox27 = afwGeom.Box2I(afwGeom.Point2I(832, 0), afwGeom.Extent2I(864, 1024))
        imslice[27] = maskedImage.Factory(maskedImage, bbox27)
        bbox28 = afwGeom.Box2I(afwGeom.Point2I(864, 0), afwGeom.Extent2I(896, 1024))
        imslice[28] = maskedImage.Factory(maskedImage, bbox28)
        bbox29 = afwGeom.Box2I(afwGeom.Point2I(896, 0), afwGeom.Extent2I(928, 1024))
        imslice[29] = maskedImage.Factory(maskedImage, bbox29)
        bbox30 = afwGeom.Box2I(afwGeom.Point2I(928, 0), afwGeom.Extent2I(960, 1024))
        imslice[30] = maskedImage.Factory(maskedImage, bbox30)
        bbox31 = afwGeom.Box2I(afwGeom.Point2I(960, 0), afwGeom.Extent2I(992, 1024))
        imslice[31] = maskedImage.Factory(maskedImage, bbox31)
        bbox32 = afwGeom.Box2I(afwGeom.Point2I(992, 0), afwGeom.Extent2I(1024, 1024))
        imslice[32] = maskedImage.Factory(maskedImage, bbox32)

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
    
    darkExposure = afwImage.ExposureF(os.path.join(inputdir, "DARK_10.fits.gz"))
    flatExposure = afwImage.ExposureF(os.path.join(inputdir,"Flat06Feb.fits.gz"))
    rawExposure = afwImage.ExposureF(os.path.join(inputdir,"F02_S22_10_030.fits.gz"))
    
    output = SofiIsrTask.runIsr(rawExposure, dark=darkExposure, flat=flatExposure)
    
    return output.exposure

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

    exposure = runIsr()
    
    if args.ds9:
        im = exposure.getMaskedImage().getImage()
        im_median = numpy.median(im.getArray())
        ds9.mtv(im)
        ds9.scale(min=im_median*0.90, max=im_median*1.1, type='SQRT')

    if args.write:
        exposure.writeFits("postISRCCD.fits")

#runDam(display=args.display)
