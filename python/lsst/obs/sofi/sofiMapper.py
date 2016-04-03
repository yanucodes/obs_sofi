#
# LSST Data Management System
# Copyright 2012 LSST Corporation.
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

import pyfits

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils

from lsst.daf.butlerUtils import CameraMapper, exposureFromImage
import lsst.pex.policy as pexPolicy

# Solely to get boost serialization registrations for Measurement subclasses
import lsst.meas.algorithms

class SofiMapper(CameraMapper):
    packageName = "obs_sofi"

    def __init__(self, **kwargs):
        policyFile = pexPolicy.DefaultPolicyFile("obs_sofi", "SOFIMapper.paf", "policy")
        policy = pexPolicy.Policy(policyFile)
        super(SofiMapper, self).__init__(policy, policyFile.getRepositoryPath(), **kwargs)

        # The "ccd" provided by the user is translated through the registry into an extension name for the "raw"
        # template.  The template therefore doesn't include "ccd", so we need to ensure it's explicitly included
        # so the ArgumentParser can recognise and accept it.

        self.exposures['raw'].keyDict['ccd'] = int

        afwImageUtils.defineFilter('Ks', lambdaEff=2175, alias=['KS'])
        self.filterIdMap = dict(Ks=0)

        # Ensure each dataset type of interest knows about the full range of keys available from the registry
        keys = {'pointing': str,
                'dateObs': str,
                'expTime': float,
                }
        for name in ("raw", "calexp", "postISRCCD", "src", "icSrc", "icMatch"):
            self.mappings[name].keyDict.update(keys)

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.
        @param dataId (dict) Data identifier with visit, ccd
        """
        pathId = self._transformId(dataId)

        return 0
    
    def bypass_ccdExposureId(self, datasetType, pythonType, location, dataId):
        return self._computeCcdExposureId(dataId)
    
    def bypass_ccdExposureId_bits(self, datasetType, pythonType, location, dataId):
        return 32
    
    def _extractDetectorName(self, dataId):
        return "Hawaii HgCdTe"

    def _standardizeDetrend(self, detrend, image, dataId, filter=False):
        """Hack up detrend images to remove troublesome keyword"""
        md = image.getMetadata()
        removeKeyword(md, 'RADECSYS') # Irrelevant, and use of "GAPPT" breaks wcslib
        exp = exposureFromImage(image)
        return self._standardizeExposure(self.calibrations[detrend], exp, dataId, filter=filter, trimmed=False)

    def std_dark(self, image, dataId):
        return self._standardizeDetrend("dark", image, dataId, filter=False)

    def std_flat(self, image, dataId):
        return self._standardizeDetrend("flat", image, dataId, filter=True)


def removeKeyword(md, key):
    """Remove a keyword from a header without raising an exception if it doesn't exist"""
    if md.exists(key):
        md.remove(key)
