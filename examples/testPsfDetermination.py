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

def measure(footprintSet, exposure):

    schema = afwTable.SourceTable.makeMinimalSchema()

    config = measBase.BaseMeasurementConfig()
    config.slots.psfFlux = "base_PsfFlux"
    measureTask = measBase.BaseMeasurementTask(schema, config=config)

    table = afwTable.SourceCatalog(schema)
    footprintSet.makeSources(table)
    measureTask.callMeasure(table, exposure)

    return table
'''
def measure(footprintSet, exposure):

    schema = afwTable.SourceTable.makeMinimalSchema()
    config = measBase.SingleFrameMeasurementConfig()
    config.algorithms.names = ["base_PixelFlags",
                 "base_SdssCentroid",
                 "base_GaussianFlux",
                 "base_SdssShape",
                 "base_CircularApertureFlux",
                 "base_PsfFlux",
                 ]
    config.algorithms["base_CircularApertureFlux"].radii = [3.0]
    config.slots.centroid = "base_SdssCentroid"
    config.slots.psfFlux = "base_PsfFlux"
    config.slots.apFlux = "base_CircularApertureFlux_3_0"
    config.slots.modelFlux = None
    config.slots.instFlux = None
    config.slots.calibFlux = None
    config.slots.shape = "base_SdssShape"
    measureTask = measBase.SingleFrameMeasurementTask(schema, config=config)
    table = afwTable.SourceCatalog(schema)
    footprintSet.makeSources(table)
    measureTask.run(table, exposure)

    return table
    '''

def setupDeterminer(exposure, nEigenComponents=3, starSelectorAlg="secondMoment"):
    """Setup the starSelector and psfDeterminer"""
    schema = afwTable.SourceTable.makeMinimalSchema()
    if starSelectorAlg == "secondMoment":
        starSelectorClass = measAlg.SecondMomentStarSelectorTask
        starSelectorConfig = starSelectorClass.ConfigClass()
        starSelectorConfig.clumpNSigma = 3.0
        starSelectorConfig.badFlags = [                                            ]
    elif starSelectorAlg == "objectSize":
        starSelectorClass = measAlg.ObjectSizeStarSelectorTask
        starSelectorConfig = starSelectorClass.ConfigClass()
        starSelectorConfig.badFlags = ["base_PixelFlags_flag_edge",
                                           "base_PixelFlags_flag_interpolatedCenter",
                                           "base_PixelFlags_flag_saturatedCenter",
                                           "base_PixelFlags_flag_crCenter",
                                           ]
        starSelectorConfig.widthStdAllowed = 0.5

    starSelector = starSelectorClass(config=starSelectorConfig,schema=schema)

    psfDeterminerFactory = measAlg.psfDeterminerRegistry["pca"]
    psfDeterminerConfig = psfDeterminerFactory.ConfigClass()
    width, height = exposure.getMaskedImage().getDimensions()
    psfDeterminerConfig.sizeCellX = width
    psfDeterminerConfig.sizeCellY = height//3
    psfDeterminerConfig.nEigenComponents = nEigenComponents
    psfDeterminerConfig.spatialOrder = 1
    psfDeterminerConfig.spatialReject = 6.0
    psfDeterminerConfig.kernelSizeMin = 31
    psfDeterminerConfig.nStarPerCell = 0
    psfDeterminerConfig.nStarPerCellSpatialFit = 0 # unlimited
    psfDeterminer = psfDeterminerFactory(psfDeterminerConfig)

    return starSelector, psfDeterminer

def testPsfDeterminer(exposure):
    """starSelectorAlg "secondMoment",
                                "objectSize"
        """

    starSelector, psfDeterminer = setupDeterminer(exposure,
                                                    nEigenComponents=2, starSelectorAlg="secondMoment")

    footprintSet = afwDetection.FootprintSet(exposure.getMaskedImage(), afwDetection.Threshold(100), "DETECTED")
    catalog = measure(footprintSet, exposure)

    psfCandidateList = starSelector.run(exposure, catalog).psfCandidates

    psf, cellSet = psfDeterminer.determinePsf(exposure, psfCandidateList, metadata)
    exposure.setPsf(psf)

if __name__ == '__main__':
    exposure = afwImage.ExposureF("dd_F02_S22_10_022.fits")
    psfConfig = measAlg.GaussianPsfFactory()
    psfConfig.defaultFwhm = 3.08
    psf = psfConfig.apply(3.08)
    exposure.setPsf(psf)
    testPsfDeterminer(exposure)
