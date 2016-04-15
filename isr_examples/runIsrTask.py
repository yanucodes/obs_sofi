import os
import lsst.afw.image as afwImage

def runIsr():
    '''Run the task to do ISR on a ccd'''
    from lsst.obs.sofi.runIsrTask import SofiIsrTask
    
    #Create the isr task with modified config
    isrConfig = SofiIsrTask.ConfigClass()
    isrConfig.doBias = False
    isrConfig.doOverscan = True
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
