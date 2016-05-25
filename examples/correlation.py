import os
import sys
import warnings
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.meas.algorithms as measAlg
from SourceDetectionTask import run as sd
import cv2
import matplotlib.pyplot as plt
from scipy import optimize
from pylab import *
import traceback
import glob


class GetOutOfLoop(Exception):
    pass

class TooSmallSlist(Exception):
    pass

def gaussian(height, center_x, center_y, width_x, width_y, rotation):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    
    rotation = np.deg2rad(rotation)
    center_x = center_x * np.cos(rotation) - center_y * np.sin(rotation)
    center_y = center_x * np.sin(rotation) + center_y * np.cos(rotation)
        
    def rotgauss(x,y):
        xp = x * np.cos(rotation) - y * np.sin(rotation)
        yp = x * np.sin(rotation) + y * np.cos(rotation)
        g = height*np.exp(-(((center_x-xp)/width_x)**2+((center_y-yp)/width_y)**2)/2.)
        return g
        
    return rotgauss

def parabola(c, center_x, center_y, a, b):
    a = float(a)
    b = float(b)
    return lambda x,y: c*( ((center_x-x)**2)/(a**2) + ((center_y-y)**2)/(b**2) )

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y, 0.0

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

def fitparabola(data):
    """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(parabola(*p)(*indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

def fit_parabola(arr, ncc, x):
    d = arr[x+1] - (2.0*arr[x]+arr[x-1])
    center = float(x) + 0.5 - (arr[x+1] - arr[x]) / d
    return center

def loadexp(name):
    """Prepare the data"""
    
    # Load sample input from disk
    inputdir = "."
    
    imFile = os.path.join(inputdir, name)
    
    exposure = afwImage.ExposureF(imFile)
    
    #psf = measAlg.SingleGaussianPsf(11, 11, 2)
    
    return exposure

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

def findoffsets(exposure1, exposure2, hwid, disp = False):
    
    setPsf(exposure1)
    setPsf(exposure2)
    
    sources1, table1, result1 = sd(exposure1, display=disp)
    sources2, table2, result2 = sd(exposure2, display=disp, framenumber=2)

    nx = 1024
    ny = 1024
    ncc = 2*hwid + 1
    nlist = nx*ny
    x = []
    y = []
    p = []
    j = 0
    maxnlist = 99999
    image1 = np.zeros((1024,1024))
    image2 = np.zeros((1024,1024))
    
    ixoff, iyoff = exposure1.getWcs().getPixelOrigin() - exposure2.getWcs().skyToPixel(exposure1.getWcs().getSkyOrigin())
    
    
    print "ixoff, iyoff"
    print ixoff, iyoff
    
    img1 = exposure1.getMaskedImage().getImage().getArray()
    img2 = exposure2.getMaskedImage().getImage().getArray()
    imag2 = exposure2.getMaskedImage().getImage().getArray().ravel()
    
    s2 = sources2
    slist1 = []
    slist2 = []
    
    for j in range(len(sources1)):
        if (int(sources1[j].getFootprint().getNpix())<2000):
            xc, yc = sources1[j].getCentroid()
            xc = int(xc-ixoff)
            yc = int(yc-iyoff)
            for i in range(len(s2)):
                try:
                    for k in range(xc-3,xc+4):
                        for l in range(yc-3,yc+4):
                            if (s2[i].getFootprint().contains(afwGeom.Point2I(k,l)) and s2[i].getFootprint().getNpix()<2000):
                                slist1.append(j)
                                slist2.append(i)
                                raise GetOutOfLoop
                except GetOutOfLoop:
                    pass
        
    #if (len(slist1)<10):
    #raise TooSmallSlist
        


    print "len slist", len(slist1)
    j = 0
    
    for i in slist1:
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

    print "npix", j
    N = j
    
    for i in slist2:
        x0, y0 = sources2[i].getFootprint().getBBox().getBegin()
        xN, yN = sources2[i].getFootprint().getBBox().getEnd()
        for xx in range(x0,xN):
            for yy in range(y0,yN):
                if (sources2[i].getFootprint().contains(afwGeom.Point2I(xx,yy)) and img2[xx,yy]>0.0):
                    image2[xx,yy]=img2[xx,yy]

    #endofuncomment
        
    hwid = 10
    
   

    ccimg = makeccimg(hwid, p, x, y, image2, ixoff, iyoff)

    ix0, iy0 = np.where(ccimg == ccimg.max())

    #print "ccimg", ccimg

    #check whether ix0, iy0 are on the edge

    xprof = ccimg[ix0,:].ravel()
    yprof = ccimg[:,iy0].ravel()
    
    xoff = (float(fit_parabola(xprof, ncc, iy0)) - float(hwid))
    yoff = (float(fit_parabola(yprof, ncc, ix0)) - float(hwid))
    
    print "xoff, yoff"
    print xoff, yoff
    """
    plt.imshow(ccimg)
    plt.show()
    """
    return xoff, yoff


def offsets(coaddPath, exposureListPath, exposureDirPath):
    
    expNum = 0
    refNum = 0
    xoffRef = 0.0
    yoffRef = 0.0
    xoffsets = []
    yoffsets = []
    #xoffsets.append(0.0)
    #yoffsets.append(0.0)
    exposureRef = None
    changedRef = False
    
    hwid = 10
    ncc = 2*hwid + 1
    ccimgfull = np.zeros((ncc,ncc))
    
    
    f = open(os.path.join(exposureDirPath, exposureListPath), 'w')
    
    infile = glob.glob(os.path.join(exposureDirPath, "*.fits"))
    #print infile
    i = 0
    while (i < len(infile) or (i==len(infile) and changedRef)):
        if changedRef:
            i = i - 1
            print infile[i]
        exposurePath = infile[i]
        if not exposurePath or exposurePath.startswith("#"):
            continue
        expNum += 1
            
        try:
            print >> sys.stderr, "Processing exposure: %s" % (exposurePath,)
            exposure = afwImage.ExposureF(exposurePath)
            if not exposureRef:
                exposureRef = exposure
                xoffsets.append(0.0)
                yoffsets.append(0.0)
                f.write('%s\n' % exposurePath)
            else:
                try:
                    xoff, yoff = findoffsets(exposureRef,exposure,hwid)
                    xoffsets.append(xoff+xoffRef)
                    yoffsets.append(yoff+yoffRef)
                    print >> sys.stderr, "Offsets for exposure %s: %f %f" % (exposurePath, xoff, yoff,)
                    f.write('%s\n' % exposurePath)
                    changedRef = False
                except TooSmallSlist:
                    refNum += 1
                    if (refNum<i):
                        exposurePath = infile[refNum]
                        print >> sys.stderr, "Too small slist, using exposure %s as reference" % (exposurePath)
                        exposureRef = afwImage.ExposureF(os.path.join(exposureDirPath, exposurePath))
                        xoffRef = xoffsets[refNum]
                        yoffRef = yoffsets[refNum]
                        changedRef = True
                    else:
                        print >> sys.stderr, "Exposure %s failed: %s" % (exposurePath, e)
                        print >> sys.stderr, "Not enough objects to correlate"
        except Exception, e:
            print >> sys.stderr, "Exposure %s failed: %s" % (exposurePath, e)
            traceback.print_exc(file=sys.stderr)
            xoffsets.append(0.0)
            yoffsets.append(0.0)
            continue
        i = i+1
                
    f.close()
    return xoffsets, yoffsets


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="New correlation procedure")
    
    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    parser.add_argument('--ds9', action="store_true", help="Display sources on ds9", default=False)
    
    args = parser.parse_args()
    
    if args.debug:
        try:
            import debug
        except ImportError as e:
            print >> sys.stderr, e

    name = "dd_F02_S22_10_021.fits"
    exposure1 = loadexp(name)
    name = "dd_F02_S22_10_023skysub.fits"
    exposure2 = loadexp(name)

    findoffsets( exposure1, exposure2, disp=args.ds9)
