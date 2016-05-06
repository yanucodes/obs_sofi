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

"""
Tests for PSF code

Run with:
   python psf.py
or
   python
   >>> import psf; psf.run()
"""

import math
import numpy
import unittest
import lsst.utils.tests as utilsTests
import lsst.pex.logging as logging
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.display.ds9 as ds9
import lsst.daf.base as dafBase
import lsst.afw.display.utils as displayUtils
import lsst.meas.algorithms as measAlg
import lsst.meas.base as measBase
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
from lsst.meas.algorithms.detection import SourceDetectionTask
from lsst.meas.base import BaseMeasurementTask

try:
    type(verbose)
except NameError:
    verbose = 0
    logging.Trace.setVerbosity("meas.algorithms.Interp", verbose)
    logging.Trace.setVerbosity("afw.detection.Psf", verbose)
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def psfVal(ix, iy, x, y, sigma1, sigma2, b):
    """Return the value at (ix, iy) of a double Gaussian
       (N(0, sigma1^2) + b*N(0, sigma2^2))/(1 + b)
    centered at (x, y)
    """
    return (math.exp        (-0.5*((ix - x)**2 + (iy - y)**2)/sigma1**2) +
            b*math.exp        (-0.5*((ix - x)**2 + (iy - y)**2)/sigma2**2))/(1 + b)

def measure(exposure):
    
    psf = measAlg.SingleGaussianPsf(11, 11, 2)
    exposure.setPsf(psf)
    
    """Measure a set of Footprints, returning a SourceCatalog"""
    schema = afwTable.SourceTable.makeMinimalSchema()
    config = SourceDetectionTask.ConfigClass()
    config.thresholdPolarity = "both"
    config.background.isNanSafe = True
    config.thresholdValue = 3.0
    detectionTask = SourceDetectionTask(config=config, schema=schema)
    config = BaseMeasurementTask.ConfigClass()
    config.slots.psfFlux = "base_PsfFlux"
    config.slots.apFlux = "base_CircularApertureFlux_3_0"
        
    task = measBase.BaseMeasurementTask(schema, config=config)
    table = afwTable.SourceCatalog(schema)
    result = detectionTask.run(table, exposure)
    
    sources = result.sources
    print "Found %d sources (%d +ve, %d -ve)" % (len(sources), result.fpSets.numPos, result.fpSets.numNeg)
    
    # Then run the default SFM task.  Results not checked
    task.callMeasure(sources, exposure)
    
    return table

'''
def setUp(exposure):
    width, height = 110, 301
    
    mi = exposure.getMaskedImage()
    mi.set(0)
    sd = 3                          # standard deviation of image
    mi.getVariance().set(sd*sd)
    mi.getMask().addMaskPlane("DETECTED")

    FWHM = 5
    ksize = 31                      # size of desired kernel

    sigma1 = 1.75
    sigma2 = 2*sigma1

    exposure.setPsf(measAlg.DoubleGaussianPsf(ksize, ksize,
                                                    1.5*sigma1, 1, 0.1))
                                                    
                                                    #set the right detector
    exposure.setDetector(DetectorWrapper().detector)


    #
    # Make a kernel with the exactly correct basis functions.  Useful for debugging
    #
    basisKernelList = afwMath.KernelList()
    for sigma in (sigma1, sigma2):
        basisKernel = afwMath.AnalyticKernel(self.ksize, self.ksize,
                                                 afwMath.GaussianFunction2D(sigma, sigma))
        basisImage = afwImage.ImageD(basisKernel.getDimensions())
        basisKernel.computeImage(basisImage, True)
        basisImage /= numpy.sum(basisImage.getArray())

        if sigma == sigma1:
            basisImage0 = basisImage
        else:
            basisImage -= basisImage0

        basisKernelList.append(afwMath.FixedKernel(basisImage))

    order = 1                                # 1 => up to linear
    spFunc = afwMath.PolynomialFunction2D(order)

    exactKernel = afwMath.LinearCombinationKernel(basisKernelList, spFunc)
    exactKernel.setSpatialParameters([[1.0, 0,          0],
                                          [0.0, 0.5*1e-2, 0.2e-2]])
    exactPsf = measAlg.PcaPsf(exactKernel)

    rand = afwMath.Random()               # make these tests repeatable by setting seed
    
    catalog = measure(exposure)

    for source in catalog:
        try:
            cand = measAlg.makePsfCandidate(source, exposure)
            cellSet.insertCandidate(cand)

        except Exception, e:
            print e
            continue
'''


def setupDeterminer(exposure, nEigenComponents=3, starSelectorAlg="secondMoment"):
    """Setup the starSelector and psfDeterminer"""
    if starSelectorAlg == "secondMoment":
        starSelectorClass = measAlg.SecondMomentStarSelectorTask
        starSelectorConfig = starSelectorClass.ConfigClass()
        starSelectorConfig.clumpNSigma = 5.0
        starSelectorConfig.histSize = 14
        starSelectorConfig.badFlags = ["base_PixelFlags_flag_edge",
                                           "base_PixelFlags_flag_interpolatedCenter",
                                           "base_PixelFlags_flag_saturatedCenter",
                                           "base_PixelFlags_flag_crCenter",
                                           ]
    elif starSelectorAlg == "objectSize":
        starSelectorClass = measAlg.ObjectSizeStarSelectorTask
        starSelectorConfig = starSelectorClass.ConfigClass()
        starSelectorConfig.sourceFluxField = "base_GaussianFlux_flux"
        starSelectorConfig.badFlags = ["base_PixelFlags_flag_edge",
                                           "base_PixelFlags_flag_interpolatedCenter",
                                           "base_PixelFlags_flag_saturatedCenter",
                                           "base_PixelFlags_flag_crCenter",
                                           ]
        starSelectorConfig.widthStdAllowed = 0.5

    starSelector = starSelectorClass(config=starSelectorConfig)

    psfDeterminerFactory = measAlg.psfDeterminerRegistry["pca"]
    psfDeterminerConfig = psfDeterminerFactory.ConfigClass()
    width, height = exposure.getMaskedImage().getDimensions()
    psfDeterminerConfig.sizeCellX = width
    psfDeterminerConfig.sizeCellY = height//3
    psfDeterminerConfig.nEigenComponents = nEigenComponents
    psfDeterminerConfig.spatialOrder = 1
    psfDeterminerConfig.kernelSizeMin = 31
    psfDeterminerConfig.nStarPerCell = 0
    psfDeterminerConfig.nStarPerCellSpatialFit = 0 # unlimited
    psfDeterminer = psfDeterminerFactory(psfDeterminerConfig)

    return starSelector, psfDeterminer

def testPsfDeterminer(exposure):
    """Test the (PCA) psfDeterminer"""
    """
    for starSelectorAlg in ["secondMoment",
                                "objectSize",
                                ]:
        print "Using %s star selector" % (starSelectorAlg)
        """
    
    starSelector, psfDeterminer = setupDeterminer(exposure,
                                                        nEigenComponents=2, starSelectorAlg="objectSize")
    catalog = measure(exposure)
    
    psfCandidateList = starSelector.run(exposure, catalog).psfCandidates
        
    psf, cellSet = psfDeterminer.determinePsf(exposure, psfCandidateList, metadata)
    exposure.setPsf(psf)