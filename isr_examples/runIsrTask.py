import os
from lsst.ip.isr import IsrTask
import lsst.afw.display.ds9 as ds9
from lsst.obs.sofi.makeCamera import makeCamera
import sys, numpy
import lsst.afw.image as afwImage
from lsst.afw.cameraGeom.fitsUtils import getByKey, setByKey, HeaderAmpMap, HeaderDetectorMap, DetectorBuilder

def runIsr():
    '''Run the task to do ISR on a ccd'''
    
    
    #Create the isr task with modified config
    isrConfig = IsrTask.ConfigClass()
    isrConfig.doBias = False #We didn't make a zero frame
    isrConfig.doDark = True
    isrConfig.doFlat = True
    isrConfig.doFringe = False #There is no fringe frame for this example
    
    isrConfig.assembleCcd.doRenorm = False #We'll take care of gain in the flats
    isrConfig.assembleCcd.setGain = False
    isrTask = IsrTask(config=isrConfig)
    
    darkExposure = afwImage.ExposureF(os.path.join(inputdir, "DARK_10.fits.gz"))
    flatExposure = afwImage.ExposureF(os.path.join(inputdir,"Flat06Feb.fits.gz"))
    rawExposure = afwImage.ExposureF(os.path.join(inputdir,"F02_S22_10_032.fits.gz"))
    camera, det = makeCamera()
    rawExposure.setDetector(det)
    
    output = isrTask.run(rawExposure, dark=darkExposure, flat=flatExposure)
    return output.exposure

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of IsrTask")
    parser.add_argument("--inputdir", default=".", help="Input directory")
    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    parser.add_argument('--ds9', action="store_true", help="Display the result?", default=False)
    parser.add_argument('--write', '-w', action="store_true", help="Write the result?", default=False)
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
