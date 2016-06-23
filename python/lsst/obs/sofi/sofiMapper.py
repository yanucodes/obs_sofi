#
# LSST Data Management System
# Copyright 2016 LSST Corporation.
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


import lsst.afw.image.utils as afwImageUtils
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.daf.butlerUtils import CameraMapper
import lsst.pex.policy as pexPolicy
from .monocam import Monocam
from .hack import getDatabase, fakeWcs

__all__ = ["MonocamMapper"]


class MonocamMapper(CameraMapper):
    packageName = 'obs_sofi'

    def __init__(self, inputPolicy=None, **kwargs):
        policyFile = pexPolicy.DefaultPolicyFile(self.packageName, "SOFIMapper.paf", "policy")
        policy = pexPolicy.Policy(policyFile)

        CameraMapper.__init__(self, policy, policyFile.getRepositoryPath(), **kwargs)

        getDatabase(kwargs["root"])

        # Ensure each dataset type of interest knows about the full range of keys available from the registry
        keys = {'expNum': int,
                'pointing': str,
                'filter': str,
                'dateObs': str,
                'expTime': float,
        }
        for name in ("raw", "raw_amp",
                     # processCcd outputs
                     "postISRCCD", "calexp", "postISRCCD", "src", "icSrc", "srcMatch",
                     ):
            self.mappings[name].keyDict.update(keys)


        self.filterIdMap = {'K': 0}

        # SOFI K-band filter
        afwImageUtils.defineFilter('K', lambdaEff=2175, alias=['Ks'])

    def _extractDetectorName(self, dataId):
        return "0"

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.

        @param dataId (dict) Data identifier with visit
        """
        expNum = dataId['expNum']
        return long(expNum)

    def bypass_ccdExposureId(self, datasetType, pythonType, location, dataId):
        return self._computeCcdExposureId(dataId)

    def bypass_ccdExposureId_bits(self, datasetType, pythonType, location, dataId):
        return 32

    def validate(self, dataId):
        expNum = dataId.get("expNum")
        if expNum is not None and not isinstance(expNum, int):
            dataId["expNum"] = int(expNum)
        return dataId

    def _setCcdExposureId(self, propertyList, dataId):
        propertyList.set("Computed_ccdExposureId", self._computeCcdExposureId(dataId))
        return propertyList

    def _makeCamera(self, policy, repositoryDir):
        """Make a camera (instance of lsst.afw.cameraGeom.Camera) describing the camera geometry
        """
        return Sofi()

    def bypass_defects(self, datasetType, pythonType, location, dataId):
        """ since we have no defects, return an empty list.  Fix this when defects exist """
        return []

    def _defectLookup(self, dataId):
        """ This function needs to return a non-None value otherwise the mapper gives up
        on trying to find the defects.  I wanted to be able to return a list of defects constructed
        in code rather than reconstituted from persisted files, so I return a dummy value.
        """
        return "hack"

    def bypass_raw(self, datasetType, pythonType, location, dataId):
        """Read raw image with hacked metadata"""
        filename = location.getLocations()[0]
        md = afwImage.readMetadata(filename, 1)
        removeKeyword(md, 'ESO DET CHIP PXSPACE')
        image = afwImage.DecoratedImageU(filename)
        image.setMetadata(md)
        return self.std_raw(image, dataId)

    bypass_raw_amp = bypass_raw


    def standardizeCalib(self, dataset, item, dataId):
        """Standardize a calibration image read in by the butler

        Some calibrations are stored on disk as Images instead of MaskedImages
        or Exposures.  Here, we convert it to an Exposure.

        @param dataset  Dataset type (e.g., "bias", "dark" or "flat")
        @param item  The item read by the butler
        @param dataId  The data identifier (unused, included for future flexibility)
        @return standardized Exposure
        """
        mapping = self.calibrations[dataset]
        if "MaskedImage" in mapping.python:
            exp = afwImage.makeExposure(item)
        elif "Image" in mapping.python:
            if hasattr(item, "getImage"): # For DecoratedImageX
                item = item.getImage()
            exp = afwImage.makeExposure(afwImage.makeMaskedImage(item))
        elif "Exposure" in mapping.python:
            exp = item
        else:
            raise RuntimeError("Unrecognised python type: %s" % mapping.python)

        if hasattr(CameraMapper, "std_" + dataset):
            return getattr(parent, "std_" + dataset)(self, exp, dataId)
        return self._standardizeExposure(mapping, exp, dataId)

    def _extractDetectorName(self, dataId):
        return "Hawaii HgCdTe"

    def std_dark(self, item, dataId):
        exp = self.standardizeCalib("dark", item, dataId)
        exp.getCalib().setExptime(1.0)
        return exp

    def std_flat(self, item, dataId):
        return self.standardizeCalib("flat", item, dataId)

def removeKeyword(md, key):
    """Remove a keyword from a header without raising an exception if it doesn't exist"""
    if md.exists(key):
        md.remove(key)

