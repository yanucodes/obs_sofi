#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2015 LSST Corporation.
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

import os
import sys
import numpy

import eups

import lsst.daf.base               as dafBase
import lsst.afw.table              as afwTable
import lsst.afw.image              as afwImage
import lsst.afw.display.ds9        as ds9
import lsst.meas.algorithms        as measAlg
from lsst.meas.algorithms.detection import SourceDetectionTask
from lsst.meas.base import SingleFrameMeasurementTask
from lsst.meas.base import BaseMeasurementTask

def loadData(inputdir, name):
    """Prepare the data we need to run the example"""

    # Load sample input from disk
    #inputdir = "/LSST/SOFI/FIELDS/postISR/"

    imFile = os.path.join(inputdir, name + ".fits")

    exposure = afwImage.ExposureF(imFile)
    
    #psf = measAlg.SingleGaussianPsf(11, 11, 2)
    psfConfig = measAlg.GaussianPsfFactory()
    psfConfig.defaultFwhm = 3.08
    psf = psfConfig.apply(3.08)
    
    exposure.setPsf(psf)

    im = exposure.getMaskedImage().getImage()
    im -= float(numpy.median(im.getArray()))

    return exposure

def run(exposure, display=False, framenumber=1, threshold=5.0):
    
    schema = afwTable.SourceTable.makeMinimalSchema()
    #
    # Create the detection task
    #
    config = SourceDetectionTask.ConfigClass()
    config.background.isNanSafe = True
    config.thresholdValue = threshold
    detectionTask = SourceDetectionTask(config=config, schema=schema)
    #
    # And the measurement Task
    #
    config = SingleFrameMeasurementTask.ConfigClass()
    config.plugins.names.clear()
    for plugin in ["base_SdssCentroid", "base_SdssShape", "base_CircularApertureFlux", "base_GaussianFlux"]:
        config.plugins.names.add(plugin)
    config.slots.psfFlux = None
    config.slots.apFlux = "base_CircularApertureFlux_3_0"

    measureTask = SingleFrameMeasurementTask(schema, config=config)

    #
    # Print the schema the configuration produced
    #
    #print schema

    #
    # Create the output table
    #
    tab = afwTable.SourceTable.make(schema)
    #
    # Process the data
    #
    result = detectionTask.run(tab, exposure)

    sources = result.sources

    print "Found %d sources (%d +ve, %d -ve)" % (len(sources), result.fpSets.numPos, result.fpSets.numNeg)

    measureTask.run(sources, exposure)
    
    #sources.writeFits('sources.fits')
    

    if display:                         # display on ds9 (see also --debug argparse option)
        frame = framenumber
        ds9.mtv(exposure, frame=frame)

        with ds9.Buffering():
            for s in sources:
                xy = s.getCentroid()
                t = not s.get("flags_negative")
                if t:
                    ds9.dot('*', *xy, ctype=ds9.RED if t else ds9.BLACK, frame=frame)
                #ds9.dot(s.getShape(), *xy, ctype=ds9.RED, frame=frame)
                #ds9.dot('o', *xy, size=config.plugins["base_CircularApertureFlux"].radii[0],
                #ctype=ds9.YELLOW, frame=frame)

    return sources, tab, result

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of Source{Detection,Measurement}Task")
    parser.add_argument("--dir", default="/LSST/SOFI/FIELDS/postISR/", help="Input directory")
    parser.add_argument("--fn", default="coadd", help="Filename")
    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    parser.add_argument('--ds9', action="store_true", help="Display sources on ds9", default=False)

    args = parser.parse_args()

    if args.debug:
        try:
            import debug
        except ImportError as e:
            print >> sys.stderr, e

    inputdir = args.dir

    exposure = loadData(inputdir, args.fn)

    run(exposure, display=args.ds9)
