#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
"""Demonstrate how to create a coadd by warping and adding.
"""
import os
import sys
import time
import traceback
import numpy as np
from math import ceil
import lsst.pex.config as pexConfig
import lsst.pex.logging as pexLog
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.coadd.utils as coaddUtils
import lsst.meas.algorithms as measAlg
from lsst.meas.algorithms import SubtractBackgroundTask
from SourceDetectionTask import run as sd
from correlation import offsets

def setPsf(exposure):
    psfConfig = measAlg.GaussianPsfFactory()
    psfConfig.defaultFwhm = 3.4
    psf = psfConfig.apply(3.4)
    
    exposure.setPsf(psf)
    
    im = exposure.getMaskedImage().getImage()
    im -= float(np.median(im.getArray()))

class WarpAndCoaddConfig(pexConfig.Config):
    saveDebugImages = pexConfig.Field(
        doc = "Save intermediate images?",
        dtype = bool,
        default = False,
    )
    bboxMin = pexConfig.ListField(
        doc = "Lower left corner of bounding box used to subframe to all input images",
        dtype = int,
        default = (0, 0),
        length = 2,
    )
    bboxSize = pexConfig.ListField(
        doc = "Size of bounding box used to subframe all input images; 0 0 for full input images",
        dtype = int,
        default = (0, 0),
        length = 2,
    )
    coaddZeroPoint = pexConfig.Field(
        dtype = float,
        doc = "Photometric zero point of coadd (mag).",
        default = 27.0,
    )
    coadd = pexConfig.ConfigField(dtype = coaddUtils.Coadd.ConfigClass, doc = "")
    warp = pexConfig.ConfigField(dtype = afwMath.Warper.ConfigClass, doc = "")

def warpAndCoadd(coaddPath, exposureListPath, xoffsets, yoffsets, borders, config):
    """Create a coadd by warping and psf-matching
    
    Inputs:
    - coaddPath: path to desired coadd; ovewritten if it exists
    - exposureListPath: a file containing a list of paths to input exposures;
        blank lines and lines that start with # are ignored
    - config: an instance of WarpAndCoaddConfig
    
    The first exposure in exposureListPath is used as the reference: all other exposures
    are warped to match to it.
    """
    weightPath = os.path.splitext(coaddPath)[0] + "_weight.fits"

    bbox = afwGeom.Box2I(
        afwGeom.Point2I(config.bboxMin[0], config.bboxMin[1]),
        afwGeom.Extent2I(config.bboxSize[0], config.bboxSize[1]),
    )
    print "SaveDebugImages =", config.saveDebugImages
    print "bbox =", bbox
    
    #zpScaler = coadd.ZeropointScaler(config.coaddZeroPoint)
    
    #bsconfig = SubtractBackgroundTask.ConfigClass()
    #bsconfig.statisticsProperty = "MEDIAN"
    #backgroundTask = SubtractBackgroundTask(config=bsconfig)

    # process exposures
    accumGoodTime = 0
    coadd = None
    expNum = 0
    numExposuresInCoadd = 0
    numExposuresFailed = 0
    
    
    xmin, xmax, ymin, ymax = borders
    
    if xmin<0 and xmax>0:
        NX = 1024+ceil(xmax)+ceil(-xmin)+10
        dxaux = 517.0 + ceil(-xmin)
    if ymin<0 and ymax>0:
        NY = 1024+ceil(ymax)+ceil(-ymin)+10
        dyaux = 517.0 + ceil(-ymin)
    if xmin>0 and xmax>0:
        NX = 1024+ceil(xmax)+10
        dxaux = 517.0
        print 'xmin>0'
    if ymin>0 and ymax>0:
        NY = 1024+ceil(ymax)+10
        dyaux = 517.0
        print 'ymin>0'
    if ymin<0 and ymax<0:
        NY = 1024+ceil(-ymin)+10
        dyaux = 517.0 + ceil(-ymin)
    if xmin<0 and xmax<0:
        NX = 1024+ceil(-xmin)+10
        dxaux = 517.0 + ceil(-xmin)


    auxMaskedImage = afwImage.makeMaskedImage(afwImage.makeImageFromArray(np.zeros((int(NX),int(NY)))))
    aux = afwImage.makeExposure(auxMaskedImage)
    
    bkgmean = afwImage.MaskedImageF('bkgmean.fits')
    
    with file(exposureListPath, "rU") as infile:
        for exposurePath in infile:
            exposurePath = exposurePath.strip()
            if not exposurePath or exposurePath.startswith("#"):
                continue
            expNum += 1

            try:
                print >> sys.stderr, "Processing exposure: %s" % (exposurePath,)
                startTime = time.time()
                exposure = afwImage.ExposureF(exposurePath)
                #backgroundTask.run(exposure=exposure)
                setPsf(exposure)
                sd(exposure, display=False, threshold=5.0)
                
                mi = exposure.getMaskedImage()
                mi+= bkgmean
                
                if config.saveDebugImages:
                    exposure.writeFits("exposure%s.fits" % (expNum,))
                
                if not coadd:
                    print >> sys.stderr, "Create warper and coadd with size and WCS matching the first/reference exposure"
                    
                    exposureRef = exposure
                    
                    cd11, cd12 = exposure.getWcs().getCDMatrix()[0]
                    cd21, cd22 = exposure.getWcs().getCDMatrix()[1]
                    
                    metadata = exposure.getWcs().getFitsMetadata()
                    metadata.setDouble("CRPIX1",     dxaux)
                    metadata.setDouble("CRPIX2",     dyaux)
        
                    auxWcs = afwImage.makeWcs(metadata)
                    aux.setWcs(auxWcs)
                    print aux.getWcs().getPixelOrigin(), aux.getWcs().getSkyOrigin()
                    print exposure.getWcs().getPixelOrigin(), exposure.getWcs().getSkyOrigin()
                    
                    warper = afwMath.Warper.fromConfig(config.warp)
                    coadd = coaddUtils.Coadd.fromConfig(
                        bbox = aux.getBBox(),
                        wcs = auxWcs,
                        config = config.coadd)
                    print "badPixelMask=", coadd.getBadPixelMask()
                    
                    warpedExposure = warper.warpExposure(
                                                         destWcs = auxWcs,
                                                         srcExposure = exposure,
                                                         maxBBox = aux.getBBox(),
                                                         )
                    
                    print >> sys.stderr, "Add reference exposure to coadd (without warping)"
                    coadd.addExposure(warpedExposure)
                else:
                    xoff = xoffsets[expNum-2]
                    yoff = yoffsets[expNum-2]
                    
                    metadata = exposure.getWcs().getFitsMetadata()
                    
                    metadata.setDouble("CRPIX1",     512.0 + xoff)
                    metadata.setDouble("CRPIX2",     512.0 + yoff)
                    
                    exposure.setWcs(afwImage.makeWcs(metadata))
                    
                    print >> sys.stderr, "Warp exposure"
                    warpedExposure = warper.warpExposure(
                        destWcs = coadd.getWcs(),
                        srcExposure = exposure,
                        maxBBox = coadd.getBBox(),
                    )
                    if config.saveDebugImages:
                        warpedExposure.writeFits("warped%s.fits" % (expNum,))
                    
                    #print >> sys.stderr, "Scale exposure to desired photometric zeropoint"
                    #zpScaler.scaleExposure(warpedExposure)

                    print >> sys.stderr, "Add warped exposure to coadd"
                    coadd.addExposure(warpedExposure)

                    # ignore time for first exposure since nothing happens to it
                    deltaTime = time.time() - startTime
                    print >> sys.stderr, "Elapsed time for processing exposure: %0.1f sec" % (deltaTime,)
                    accumGoodTime += deltaTime
                numExposuresInCoadd += 1
            except Exception, e:
                print >> sys.stderr, "Exposure %s failed: %s" % (exposurePath, e)
                traceback.print_exc(file=sys.stderr)
                numExposuresFailed += 1
                continue

    coaddExposure = coadd.getCoadd()
    arr = np.isnan(coaddExposure.getMaskedImage().getImage().getArray())
    for i in range(int(NX)):
        for k in range(int(NY)):
            if arr[i,k]:
                coaddExposure.getMaskedImage().getImage().set(k,i,0.0)

    coaddExposure.writeFits(coaddPath)
    print >> sys.stderr, "Wrote coadd: %s" % (coaddPath,)
    weightMap = coadd.getWeightMap()
    weightMap.writeFits(weightPath)
    print >> sys.stderr, "Wrote weightMap: %s" % (weightPath,)

    print >> sys.stderr, "Coadded %d exposures and failed %d" % (numExposuresInCoadd, numExposuresFailed)
    if numExposuresInCoadd > 1:
        timePerGoodExposure = accumGoodTime / float(numExposuresInCoadd - 1)
        print >> sys.stderr, "Processing speed: %.1f seconds/exposure (ignoring first and failed)" % \
            (timePerGoodExposure,)

if __name__ == "__main__":
    pexLog.Trace.setVerbosity('lsst.coadd', 3)
    helpStr = """Usage: warpAndCoadd.py coaddPath exposureListPath

where:
- coaddPath is the desired name or path of the output coadd
- exposureListPath is a file containing a list of:
    pathToExposure
  where:
  - pathToExposure is the path to an Exposure
  - the first exposure listed is taken to be the reference exposure,
    which determines the size and WCS of the coadd
  - empty lines and lines that start with # are ignored.
"""
    if len(sys.argv) != 4:
        print helpStr
        sys.exit(0)
    
    coaddPath = sys.argv[1]
    '''
    if os.path.exists(coaddPath):
        print >> sys.stderr, "Coadd file %s already exists" % (coaddPath,)
        sys.exit(1)
        '''
    
    exposureListPath = sys.argv[2]

    exposureDirPath = sys.argv[3]
    
    config = WarpAndCoaddConfig()
    
    xoffsets, yoffsets, borders = offsets(coaddPath, exposureListPath, exposureDirPath)
    
    warpAndCoadd(coaddPath, exposureListPath, xoffsets, yoffsets, borders, config)
