#based on ip_isr/examples/exampleUtils.py and lsst.obs.sdss.makeCamera.py

import numpy

from lsst.afw.cameraGeom import makeCameraFromCatalogs, CameraConfig, DetectorConfig,\
    SCIENCE, PIXELS, PUPIL, FOCAL_PLANE
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage

#
# Make an Amp
#

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

    bbox.shift(afwGeom.ExtentI(ix*nx, iy*ny))
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

#
# Make a Ccd out of 4 Amps
#

def createDetector(nAmpX, nAmpY, nPixX, nPixY, pre, hOscan, vOscan, ext, isPerAmp):
    '''!Fill ampInfo tables
    \param[in] nAmpX -- Number of amps in the x direction
    \param[in] nAmpY -- Number of amps in the y direction
    \param[in] nPixX -- Number of pixels in the amp in the x direction
    \param[in] nPixY -- Number of pixels in the amp in the y direction
    \param[in] pre -- Number of prescan rows
    \param[in] hOscan -- Number of horizontal overscan columns
    \param[in] vOscan -- Number of vertical overscan rows
    \param[in] ext -- Number of pixels in the extended register
    \param[in] isPerAmp -- Are the raw amp data in separate images?
    \return an lsst.afw.cameraGeom.Detector object
    '''
    
    schema = afwTable.AmpInfoTable.makeMinimalSchema()
    ampCatalog = afwTable.AmpInfoCatalog(schema)
    for iy in range(nAmpY):
        for ix in range(nAmpX):
            record = ampCatalog.addNew()
            populateAmpBoxes(nPixX, nPixY, pre, hOscan, vOscan, ext, ix, iy,
                              isPerAmp, record)

    detConfig = DetectorConfig()
    detConfig.name = 'HawaiiHgCdTe'
    detConfig.id = 0
    detConfig.bbox_x0 = 0
    detConfig.bbox_y0 = 0
    detConfig.bbox_x1 = nAmpX*nPixX - 1
    detConfig.bbox_y1 = nAmpY*nPixY - 1
    detConfig.detectorType = 0 #Science type
    detConfig.serial = '0'
    detConfig.offset_x = 0.
    detConfig.offset_y = 0.
    detConfig.refpos_x = nAmpX*nPixX*0.5 - 0.5
    detConfig.refpos_y = nAmpY*nPixY*0.5 - 0.5
    detConfig.yawDeg = 0.
    detConfig.pitchDeg = 0.
    detConfig.rollDeg = 0.
    detConfig.pixelSize_x = 18.5/1000. #in mm
    detConfig.pixelSize_y = 18.5/1000. #in mm
    detConfig.transposeDetector = False
    detConfig.transformDict.nativeSys = PIXELS.getSysName()
    
    return {'ccdConfig':detConfig, 'ampInfo':ampCatalog}

#
# Make a Camera
#
def makeCamera(name="SOFI"):
    """Make a camera
        @param name: name of the camera
        @param outputDir: If not None, write the objects used to make the camera to this location
        @return a camera object
        """
    camConfig = CameraConfig()
    camConfig.name = name
    camConfig.detectorList = {}
    camConfig.plateScale = 0.288/0.0185 # arcsec/mm
    pScaleRad = afwGeom.arcsecToRad(camConfig.plateScale)
    
    '''The following part probably needs to be corrected'''
    radialDistortCoeffs = [0.0, 1.0/pScaleRad]
    tConfig = afwGeom.TransformConfig()
    tConfig.transform.name = 'inverted'
    radialClass = afwGeom.xyTransformRegistry['radial']
    tConfig.transform.active.transform.retarget(radialClass)
    tConfig.transform.active.transform.coeffs = radialDistortCoeffs
    tmc = afwGeom.TransformMapConfig()
    tmc.nativeSys = FOCAL_PLANE.getSysName()
    tmc.transforms = {PUPIL.getSysName():tConfig}
    camConfig.transformDict = tmc
    
    ccdId = 0
    ampInfoCatDict = {}
    detName = "HawaiiHgCdTe"
    det = createDetector(2, 2, 512, 512, 0, 0, 0, 0, False)
    ampInfoCatDict[detName] = det['ampInfo']
    camConfig.detectorList[ccdId] = det['ccdConfig']

    return makeCameraFromCatalogs(camConfig, ampInfoCatDict)


#
# Print a Camera
#
def printCamera(title, camera):
    """Print information about a camera
        @param title: title for camera output
        @param camera: Camera object to use to print the information
        """
    print title, "Camera:", camera.getName()
    
    for det in camera:
        print "%s %dx%d centre (mm): %s" % \
            (det.getName(),
             det.getBBox().getWidth(), det.getBBox().getHeight(),
             det.getCenter(FOCAL_PLANE).getPoint())

#************************************************************************************************************

def main():
    camera = makeCamera("SOFI")
    
    print
    printCamera("", camera)

if __name__ == "__main__":
    main()
