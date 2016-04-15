import os
from lsst.ip.isr import IsrTask, isr
import lsst.afw.display.ds9 as ds9
import makeCamera
import sys, numpy
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
from lsst.obs.sofi.makeCamera import createDetector
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class SofiIsrTaskConfig(IsrTask.ConfigClass):
    
    doOverscan = pexConfig.Field(
        dtype = bool,
        doc = "Apply overscan correction?",
        default = True,
        )
        
    def setDefaults(self):
        IsrTask.ConfigClass.setDefaults(self)


class SofiIsrTask(IsrTask):
    
    ConfigClass = SofiIsrTaskConfig
    _DefaultName = "isr"
    dataPrefix = ""
    
    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
    
    def readIsrData(self, dataRef, rawExposure):
        """!Retrieve necessary frames for instrument signature removal
            \param[in] dataRef -- a daf.persistence.butlerSubset.ButlerDataRef
            of the detector data to be processed
            \param[in] rawExposure -- a reference raw exposure that will later be
            corrected with the retrieved calibration data;
            should not be modified in this method.
            \return a pipeBase.Struct with fields containing kwargs expected by run()
            - dark: exposure of dark frame
            - flat: exposure of flat field
            """
        darkExposure = self.getIsrExposure(dataRef, "dark") if self.config.doDark else None
        flatExposure = self.getIsrExposure(dataRef, "flat") if self.config.doFlat else None

        return pipeBase.Struct(dark = darkExposure, flat = flatExposure)
    
    def doOverscanCorrection(self, exposure):
        lowerDataBBox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(1024, 1011))
        lowerOverscanBBox = afwGeom.Box2I(afwGeom.Point2I(0, 1012), afwGeom.Extent2I(1024, 12))
                                      
        maskedImage = exposure.getMaskedImage()
        lowerDataView = maskedImage.Factory(maskedImage, lowerDataBBox)
                                      
        expImage = exposure.getMaskedImage().getImage()
        lowerOverscanImage = expImage.Factory(expImage, lowerOverscanBBox)
        
        isr.overscanCorrection(
                               ampMaskedImage = lowerDataView,
                               overscanImage = lowerOverscanImage,
                               fitType = self.config.overscanFitType,
                               order = self.config.overscanOrder,
                               collapseRej = self.config.overscanRej,
                               )

    @pipeBase.timeMethod
    def runIsr(self, rawExposure, dark=None,  flat=None):
        '''Run the task to do ISR on a ccd'''

        det = createDetector(2, 2, 512, 512, 0, 0, 0, 0, False)
        rawExposure.setDetector(det['detector'])
        
        maskedExposure = rawExposure.getMaskedImage()
        
        if self.config.doOverscan:
            self.doOverscanCorrection(rawExposure)
        
        return IsrTask.run(self, rawExposure, dark = dark, flat = flat)
