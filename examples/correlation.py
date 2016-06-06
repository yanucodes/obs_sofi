import os
import sys
import warnings
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.meas.algorithms as measAlg
from SourceDetectionTask import run as sd
import traceback
import glob
from scipy.optimize import curve_fit

class GetOutOfLoop(Exception):
    pass

def parabola(c, center_x, center_y, a, b):
    a = float(a)
    b = float(b)
    return lambda x,y: c*( ((center_x-x)**2)/(a**2) + ((center_y-y)**2)/(b**2) )

def fit_parabola(arr, ncc, x):
    print arr[x-1], arr[x], arr[x+1]
    d = arr[x+1] - (2.0*arr[x]+arr[x-1])
    print d
    center = float(x) + 0.5 - (arr[x+1] - arr[x]) / d
    return center

def setPsf(exposure):
    psfConfig = measAlg.GaussianPsfFactory()
    psfConfig.defaultFwhm = 3.4
    psf = psfConfig.apply(3.4)
    
    exposure.setPsf(psf)
    
    im = exposure.getMaskedImage().getImage()
    im -= float(np.median(im.getArray()))

def makeccimg(hwid, p, x, y, image, ixoff, iyoff):
    ncc = 2*hwid + 1
    ccimg = np.zeros((ncc,ncc))
    N = len(p)
    
    for k in range (-hwid,hwid+1):
        for l in range (-hwid,hwid+1):
            for i in range(N):
                xx=x[i]-int(ixoff)+k
                yy=y[i]-int(iyoff)+l
                if (xx>=0 and xx<1024 and yy>=0 and yy<1024):
                    if (p[i]>0.0 and image[xx,yy]>0.0):
                        ccimg[k+hwid][l+hwid]=ccimg[k+hwid][l+hwid] + p[i]*image[xx,yy]
    return ccimg



def findoffsets(exposure1, exposure2, sources1, sources2, hwid, disp = False):
    
    nx = 1024
    ny = 1024
    ncc = 2*hwid + 1
    nlist = nx*ny
    x = []
    y = []
    p = []
    j = 0
    maxnlist = 99999
    image1 = np.zeros((nx,ny))
    image2 = np.zeros((nx,ny))
    
    ixoff, iyoff = exposure1.getWcs().getPixelOrigin() - exposure2.getWcs().skyToPixel(exposure1.getWcs().getSkyOrigin())
    
    img1 = exposure1.getMaskedImage().getImage().getArray()
    img2 = exposure2.getMaskedImage().getImage().getArray()
    imag2 = exposure2.getMaskedImage().getImage().getArray().ravel()
    
    s2 = sources2
    
    sum = 0.0
    for i in range(len(sources1)):
        x0, y0 = sources1[i].getFootprint().getBBox().getBegin()
        xN, yN = sources1[i].getFootprint().getBBox().getEnd()
        for xx in range(x0,xN):
            for yy in range(y0,yN):
                if (sources1[i].getFootprint().contains(afwGeom.Point2I(xx,yy)) and img1[xx,yy]>0.0):
                    if (j<maxnlist):
                        x.append(int(xx))
                        y.append(int(yy))
                        p.append(img1[xx,yy])
                        image1[xx,yy]=img1[xx,yy]
                        j = j +1
        sum = 0.0


    print "npix", j
    N = j
    
    sum = 0.0
    for i in range(len(sources2)):
        x0, y0 = sources2[i].getFootprint().getBBox().getBegin()
        xN, yN = sources2[i].getFootprint().getBBox().getEnd()
        for xx in range(x0,xN):
            for yy in range(y0,yN):
                if (sources2[i].getFootprint().contains(afwGeom.Point2I(xx,yy)) and img2[xx,yy]>0.0):
                    image2[xx,yy]=img2[xx,yy]
        sum = 0.0



    ccimg = makeccimg(hwid, p, x, y, image2, ixoff, iyoff)
    ix0, iy0 = np.where(ccimg == ccimg.max())

    print ccimg.max()

    xprof = ccimg[:,iy0].ravel()
    yprof = ccimg[ix0,:].ravel()

    print ix0, iy0
    print xprof
    
    xoff = (float(fit_parabola(xprof, ncc, ix0)) - float(hwid)) + ixoff - int(ixoff)
    yoff = (float(fit_parabola(yprof, ncc, iy0)) - float(hwid)) + iyoff - int(iyoff)
    
    print "offsets", ixoff + xoff, iyoff + yoff
    
    
    print "xoff, yoff"
    print xoff, yoff
    
    return xoff, yoff


def offsets(coaddPath, exposureListPath, exposureDirPath):
    
    expNum = 0
    refNum = 0
    xoffRef = 0.0
    yoffRef = 0.0
    exposureRef = None
    changedRef = False
    
    hwid = 20
    ncc = 2*hwid + 1
    ccimgfull = np.zeros((ncc,ncc))
    
    
    f = open(os.path.join(exposureDirPath, exposureListPath), 'w')
    
    
    infile = glob.glob(os.path.join(exposureDirPath, "bs_dd_F02_S22_*.fits"))
    
    flist = np.zeros(len(infile))
    s = []
        
    for j in range(len(infile)):
        exposure = afwImage.ExposureF(infile[j])
        setPsf(exposure)
        sources, table, result = sd(exposure, display=False, threshold=5.0)

        i = 0
        while (i<len(sources)):
            if (int(sources[i].getFootprint().getNpix())>100 or sources[j].get("flags_negative")):
                del sources[i]
            else:
                i+=1
        flist[j]=len(sources)
        s.append(sources)
        
    indices = np.where(flist == flist.max())[0]
    print indices[0]
    i = int(indices[0])

    expRef = afwImage.ExposureF(infile[i])
    exposureRef = expRef
    f.write('%s\n' % infile[i])
    del infile[i]
    refSs = s[i]
    refSources = refSs
    del s[i]
    xoffsets = [0.0 for i in range(len(infile))]
    yoffsets = [0.0 for i in range(len(infile))]
    files=[]
    
    #print infile
    i = 0
    Ende = False
    
    while (len(infile)>0):
        exposurePath = infile[0]
        try:
            print >> sys.stderr, "Processing exposure: %d %s" % (i, exposurePath,)
            exposure = afwImage.ExposureF(exposurePath)
            xoff, yoff = findoffsets(exposureRef,exposure,refSources,s[i],hwid)
            xoffsets[i]=(xoff+xoffRef)
            yoffsets[i]=(yoff+yoffRef)
            print >> sys.stderr, "Offsets for exposure %s: %f %f" % (exposurePath, xoff, yoff,)
            i+=1
            files.append(infile[0])
            del infile[0]
            #refNum = 0
            xoffRef = 0.0
            yoffRef = 0.0
        except Exception, e:
            print >> sys.stderr, "Exposure %s failed: %s" % (exposurePath, e)
            traceback.print_exc(file=sys.stderr)
            continue


    for i in range(len(files)):
        f.write('%s\n' % files[i])
    
    f.close()
    
    
    return xoffsets, yoffsets
