#based on ip_isr/examples/exampleUtils.py and lsst.obs.sdss.makeCamera.py and obs.monocam
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
import numpy
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.geom as afwGeom
from lsst.afw.table import AmpInfoCatalog, AmpInfoTable, LL
from lsst.afw.cameraGeom.cameraFactory import makeDetector


class Sofi(cameraGeom.Camera):
    """The monocam Camera

    There is one ccd with name "0"
    It has sixteen amplifiers with names like "00".."07" and "10".."17"

    Standard keys are:
    amp: amplifier name: one of 00, 01, 02, 03, 04, 05, 06, 07, 10, 11, 12, 13, 14, 15, 16, 17
    ccd: ccd name: always 0
    visit: exposure number; this will be provided by the DAQ
    """
    # Taken from fit4_20160413-154303.pdf
    gain = {(0, 0): 5.271,
            (1, 0): 5.271,
            (0, 1): 5.271,
            (1, 1): 5.271}

    readNoise = {(0, 0): 2.239,
                 (1, 0): 2.092,
                 (0, 1): 2.117,
                 (1, 1): 2.092}

    def __init__(self):
        """Construct a TestCamera
        """
        plateScale = afwGeom.Angle((0.288/0.0185), afwGeom.arcseconds)  # plate scale, in angle on sky/mm
        radialDistortion = 0.  # radial distortion in mm/rad^2
        radialCoeff = numpy.array((0.0, 1.0, 0.0, radialDistortion)) / plateScale.asRadians()
        focalPlaneToPupil = afwGeom.RadialXYTransform(radialCoeff)
        pupilToFocalPlane = afwGeom.InvertedXYTransform(focalPlaneToPupil)
        cameraTransformMap = cameraGeom.CameraTransformMap(cameraGeom.FOCAL_PLANE,
                                                           {cameraGeom.PUPIL: pupilToFocalPlane})
        detectorList = self._makeDetectorList(pupilToFocalPlane, plateScale)
        cameraGeom.Camera.__init__(self, "sofi", detectorList, cameraTransformMap)

    def _makeDetectorList(self, focalPlaneToPupil, plateScale):
        """!Make a list of detectors

        @param[in] focalPlaneToPupil  lsst.afw.geom.XYTransform from FOCAL_PLANE to PUPIL coordinates
        @param[in] plateScale  plate scale, in angle on sky/mm
        @return a list of detectors (lsst.afw.cameraGeom.Detector)
        """
        detectorList = []
        detectorConfigList = self._makeDetectorConfigList()
        for detectorConfig in detectorConfigList:
            ampInfoCatalog = self._makeAmpInfoCatalog()
            detector = makeDetector(detectorConfig, ampInfoCatalog, focalPlaneToPupil,
                                    plateScale.asArcseconds())
            detectorList.append(detector)
        return detectorList

    def _makeDetectorConfigList(self):
        """!Make a list of detector configs

        @return a list of detector configs (lsst.afw.cameraGeom.DetectorConfig)
        """
        # There is only a single detector assumed perfectly centered and aligned.
        detConfig = cameraGeom.DetectorConfig()
        detConfig.name = 'HawaiiHgCdTe'
        detConfig.id = 0
        detConfig.serial = '0'
        detConfig.detectorType = 0
        # This is the orientation we need to put the serial direciton along the x-axis
        detConfig.bbox_x0 = 0
        detConfig.bbox_x1 = 1023
        detConfig.bbox_y0 = 0
        detConfig.bbox_y1 = 1023
        detConfig.pixelSize_x = 0.0185  # in mm
        detConfig.pixelSize_y = 0.0185  # in mm
        detConfig.transformDict.nativeSys = 'Pixels'
        detConfig.transformDict.transforms = None
        detConfig.refpos_x = 511.5
        detConfig.refpos_y = 511.5
        detConfig.offset_x = 0.0
        detConfig.offset_y = 0.0
        detConfig.transposeDetector = False
        detConfig.pitchDeg = 0.0
        detConfig.yawDeg = 0.0  # this is where chip rotation goes in.
        detConfig.rollDeg = 0.0
        return [detConfig]

    def _makeAmpInfoCatalog(self):
        """Construct an amplifier info catalog
        """

        schema = afwTable.AmpInfoTable.makeMinimalSchema()
        ampCatalog = AmpInfoCatalog(schema)
        
        iy = 0
        for ix in range(nAmpX):
            record = ampCatalog.addNew()
            self.populateAmpBoxes(nPixX, nPixY, pre, hOscan, 0, ext, ix, iy,
                         isPerAmp, record)
        iy = 1
            for ix in range(nAmpX):
            record = ampCatalog.addNew()
            self.populateAmpBoxes(nPixX, nPixY-vOscan, pre, hOscan, vOscan, ext, ix, iy,
                         isPerAmp, record)
        
        return ampCatalog


    def populateAmpBoxes(nx, ny, nprescan, nhoverscan, nvoverscan, nextended, ix, iy,
                      isPerAmp, record):
        '''!Fill ampInfo tables
        \param[in] isPerAmp -- If True, return a dictionary of amp exposures keyed by amp name.
                           If False, return a single exposure with amps mosaiced preserving non-science pixels
                           (e.g. overscan)
        \param[in] nx -- number of pixels in the serial register
        \param[in] ny -- number of rows in the parallel direction
        \param[in] nprescan -- number of prescan rows
        \param[in] nhoverscan -- number of horizonatal overscan columns
        \param[in] nvoverscan -- number of vertical overscan rows
        \param[in] nextended -- number of pixels in the extended register
        \param[in] ix -- index in x direction of the amp in the chip
        \param[in] iy -- index in y direction of the amp in the chip
        \param[in] isPerAmp -- are the raw data per amp or assembled into a mosaiced image
        \param[in, out] record -- record to add this amp to
        '''
    
        readNoise = {}
        readNoise[0,0] = 2.239
        readNoise[1,0] = 2.092
        readNoise[0,1] = 2.117
        readNoise[1,1] = 2.092
    
        def makeBbox(x0, y0, x_extent, y_extent):
            return afwGeom.BoxI(afwGeom.PointI(x0, y0), afwGeom.ExtentI(x_extent, y_extent))

        bbox = makeBbox(0, 0, nx, ny)

        dataBox = makeBbox(0, 0, nx, ny)
        dataBox.shift(afwGeom.ExtentI(nextended, nprescan))

        allBox = afwGeom.BoxI()

        preBox = makeBbox(0, 0, nx, nprescan)
        preBox.shift(afwGeom.ExtentI(nextended, 0))

        extBox = makeBbox(0, 0, nextended, ny)
        extBox.shift(afwGeom.ExtentI(0, nprescan))

        hOscanBox = makeBbox(0, 0, nhoverscan, ny)
        hOscanBox.shift(afwGeom.ExtentI(nextended+nx, nprescan))

        vOscanBox = makeBbox(0, 0, nx, nvoverscan)
        vOscanBox.shift(afwGeom.ExtentI(nextended, nprescan+ny))

        allBox.include(dataBox)
        allBox.include(preBox)
        allBox.include(extBox)
        allBox.include(hOscanBox)
        allBox.include(vOscanBox)

        bbox.shift(afwGeom.ExtentI(ix*512, iy*512))
        xtot = allBox.getDimensions().getX()
        ytot = allBox.getDimensions().getY()
        rShiftExt = afwGeom.ExtentI(ix*xtot, iy*ytot)
        #Set read corner in assembled coordinates
        record.setReadoutCorner(afwTable.LL)

        if not isPerAmp:

            allBox.shift(rShiftExt)
            dataBox.shift(rShiftExt)
            preBox.shift(rShiftExt)
            extBox.shift(rShiftExt)
            hOscanBox.shift(rShiftExt)
            vOscanBox.shift(rShiftExt)
            rawXoff = 0
            rawYoff = 0

        else:
            #We assume that single amp images have the first pixel read in the
            #lower left and that the pixels are arrange such that the
            #serial is along the x-axis.
            rawXoff = rShiftExt.getX()
            rawYoff = rShiftExt.getY()

        record.setBBox(bbox)

        record.setName("A:%i,%i"%(ix, iy))
        record.setGain(5.271)
        record.setReadNoise(readNoise[ix,iy])
        record.setSaturation(10000)
        record.setLinearityCoeffs((2.407751840712499E-4, -0.2882800562336101, -0.2883747432059056, -1.875847827836264E-4)) #linear astrometric coefficients taken from http://www.ls.eso.org/sci/facilities/lasilla/instruments/sofi/tecdoc/distorion_LF.ps.gz
        record.setLinearityType('Polynomial')
        record.setHasRawInfo(True)
        record.setRawFlipX(False)
        record.setRawFlipY(False)
        record.setRawBBox(allBox)
        record.setRawXYOffset(afwGeom.ExtentI(rawXoff, rawYoff))
        record.setRawDataBBox(dataBox)
        record.setRawHorizontalOverscanBBox(hOscanBox)
        record.setRawVerticalOverscanBBox(vOscanBox)
        record.setRawPrescanBBox(preBox)
