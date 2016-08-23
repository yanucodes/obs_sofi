#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
"""Example use of SubtractBackgroundTask
"""
from __future__ import absolute_import, division, print_function
import os.path
import glob
import lsst.utils
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.meas.algorithms import SubtractBackgroundTask
import numpy as np
from astropy.io import fits

def modeFromArray(arr):
    min = arr.min()
    max = arr.max()
    
    bins = np.zeros(int(max)-int(min)+1)
    
    arr = arr.ravel()
    
    for val in arr:
        bins[int(val)-int(min)]+=1
    
    mode = np.where(bins == bins.max())
    mode = mode[0]
    #print(mode)

    bins = np.arange(min,max,1.0)
    mode = bins[int(mode[0])]

    return mode

def mode(exposure):
    
    image = exposure.getMaskedImage().getImage()
    s = afwMath.makeStatistics(image, afwMath.MIN | afwMath.MAX)
    min = s.getValue(afwMath.MIN)
    max = s.getValue(afwMath.MAX)
    #print(min,max)

    bins = np.zeros(int(max)-int(min)+1)

    arr = image.getArray().ravel()
    #print(np.mean(arr))

    for val in arr:
        bins[int(val)-int(min)]+=1

    mode = np.where(bins == bins.max())
    mode = mode[0]
    bins = np.arange(min,max,1.0)
    mode = bins[int(mode)]

    #print(mode)
    return mode

def loadExp(mypath):
    """Load the data we need to run the example"""
    imFile = os.path.join(mypath, "dd_F02_S22_10_021.fits")
    return afwImage.ExposureF(imFile)

def subtractBkg(exposure, binSize = 256, algorithm = "NATURAL_SPLINE"):
    # create the task
    config = SubtractBackgroundTask.ConfigClass()
    config.algorithm = algorithm
    config.statisticsProperty = "MEANCLIP"
    config.binSize = binSize
    backgroundTask = SubtractBackgroundTask(config=config)
    
    # subtract an initial estimate of background level
    return backgroundTask.run(exposure=exposure)


def calculateBkg(exposure, binSize = 256, algorithm = "NATURAL_SPLINE"):
    """Subtract background
    """
    
    bgRes = subtractBkg(exposure, binSize = binSize, algorithm = algorithm)

    background = bgRes.background

    # compute mean and variance of the background
    backgroundImage = background.getImage()
    bkgArr = backgroundImage.getArray()
    s = afwMath.makeStatistics(backgroundImage, afwMath.MEAN | afwMath.VARIANCE)
    bgmean = s.getValue(afwMath.MEAN)
    bgvar = s.getValue(afwMath.VARIANCE)
    print("background mean=%0.1f; variance=%0.1f" % (bgmean, bgvar))
    
    """
    backgroundImage.writeFits("background.fits")
    backgroundSubImage = exposure.getMaskedImage()
    #backgroundSubImage -= background.getImage()
    backgroundSubImage.getImage().writeFits("dd_bs_F02_S22_10_021.fits")
    """

    return bgmean, bgvar


def scale(inputdir,pointing):
    
    expPath = inputdir[0:14]+"/FIELDS/"+ pointing
    infiles = glob.glob(os.path.join(expPath, "offsets.*.txt"))
    f = open(infiles[0], 'r')
    infile = f.readlines()
    
    avgscale = 0.0
    scale = []
    data = []
    bkgs = []
    for fil in infile:
        file = fil.split()[0]
        file = str(os.path.join(inputdir,file[0:len(file)-7]))
        exposure = afwImage.ExposureF(file)
        data.append(exposure.getMaskedImage().getImage().getArray())
        exposure1 = afwImage.ExposureF(file)
        bkg, sig = calculateBkg(exposure1)
        scale.append(bkg)
        avgscale +=bkg
        bkgs.append(bkg)

    avgscale = avgscale/len(infile)
    print(avgscale)
    
    bkgmean = np.zeros((1024,1024))

    for i in range(1024):
        for k in range(1024):
            bkgmean[i,k]=avgscale

    for i in range(len(scale)):
        scale[i]=avgscale-scale[i]

    return scale, data, bkgs, bkgmean

def kselect(a, n):
    
    k = int(n/2) + (n%2 - 1)
    
    a.sort()

    return a[k]



def medplane(data, scale):
    nx, ny = data[0].shape
    buf = [0.0 for i in range(len(scale))]
    plane = np.zeros((nx,ny))
    for x in range(nx):
        for y in range(ny):
            for i in range(len(scale)):
                d = data[i]
                buf[i] = (d[x,y]) + scale[i]
            
            plane[x,y] = kselect(buf, len(scale))

    return plane

def gainMap(flat, calibdir, pointing):
    
    #m = np.mode(flat)
    m = np.median(flat.getMaskedImage().getImage().getArray())
    
    #mean, var = calculateBkg(flat, binSize = 2048, algorithm = "CONSTANT")
    mi = flat.getMaskedImage()
    #scale = 1.0/mean
    #mi *=scale
    
    gain = mi.getImage().getArray()
    gain = gain/m
    
    nbad=0
    for x in range(1024):
        for y in range(1024):
            if (gain[x,y]<0.3 or gain[x,y]>1.7):
                gain[x,y]=0.0
                nbad+=1


    dev = np.zeros((1024,1024))

    for i in range(64):
        for j in range(64):
            buf = []
            for x in range(i*16,(i+1)*16):
                for y in range(j*16,(j+1)*16):
                    if (gain[x,y]>0.0):
                        buf.append(gain[x,y])
            med = np.median(buf)
            for x in range(i*16,(i+1)*16):
                for y in range(j*16,(j+1)*16):
                    if (gain[x,y]>0.0):
                        dev[x,y]=gain[x,y]-med

    med = np.median(dev)
    sig = np.median(np.abs(dev - med))/0.6745
    l = med - 5.0*sig
    h = med + 5.0*sig

    for x in range(1024):
        for y in range(1024):
            if (dev[x,y]<l or dev[x,y]>h):
                gain[x,y]=0.0
                nbad+=1

    print(nbad)



    hdu = fits.PrimaryHDU(gain)
    hdu.writeto(os.path.join(calibdir,'gain'+pointing+'.fits'), clobber = True)

    return gain


def fixmodes(arr, bkg):
    
    fixmode = bkg
    
    for i in range(512):
        if (arr[i]>0):
            fixmode = arr[i]
            break

    for i in range(512):
        if (arr[i]==0):
            arr[i] = fixmode
        else:
            fixmode = arr[i]

    return arr

def destriperows(img, bpm, bkg):
    
    QUADSIZE = 512
    
    rowMode = np.zeros(QUADSIZE)

    for k in range(4):
        XX = int ((k%2)*QUADSIZE)
        YY = int (int(k/2)*QUADSIZE)
        buf = []
        print(XX, XX+QUADSIZE)
        print(YY, YY+QUADSIZE)
        for y in range(YY,YY+QUADSIZE):
            for x in range(XX,XX+QUADSIZE):
                if (bpm[x,y]>0):
                    buf.append(img[x,y])

            if (len(buf)>QUADSIZE/2):
                rowMode[y-YY]=modeFromArray(buf)
            else:
                rowMode[y-YY] = 0.0

        rowMode = fixmodes(rowMode,bkg)

        for y in range(YY, YY + QUADSIZE):
            mode = rowMode[y - YY]
            for x in range(XX,XX+QUADSIZE):
                img[x,y]+= (bkg - mode)

    return img


def skysub(data, bkg, gainmap, sky):
    
    nx, ny = data.shape
    
    #skybkg = modeFromArray(sky)
    skybkg = np.median(sky)
    
    imgout = np.zeros((nx,ny))

    for x in range(nx):
        for y in range(ny):
            imgout[x,y]=data[x,y] + skybkg - sky[x,y]

    #imgout = destriperows(imgout, gainmap, bkg)

    for x in range(nx):
        for y in range(ny):
            if ( gainmap[x,y]<=0.0):
                imgout[x,y]=bkg

    return imgout

def skyfilter(data, bkgs, gain, inputdir, pointing):
    
    gainmap = gain
    hwid = 3

    skybeg = 0
    
    expPath = inputdir[0:14]+"/FIELDS/"+ pointing
    print(expPath)
    
    expPath = inputdir[0:14]+"/FIELDS/"+ pointing
    infiles = glob.glob(os.path.join(expPath, "offsets.*.txt"))
    f = open(infiles[0], 'r')
    infile = f.readlines()
    
    """
    for i in range(len(data)):
        bkgs[i] = modeFromArray(data[i])
        """

    for i in range(len(data)):
        print(i)
        avgscale = 0.0
        skyend = skybeg + 2*hwid + 1
        scale = []
        buf = []
        
        for j in range(skybeg, skyend):
            if (j!=i):
                buf.append(data[j])
                avgscale+=bkgs[j]
                scale.append(bkgs[j])
                    
        avgscale = avgscale / float(len(buf))

        for j in range(len(buf)):
            scale[j]=avgscale - scale[j]
            
        sky = medplane(buf, scale)

        fimg = skysub(data[i], bkgs[i], gainmap, sky)

        #exposure.writeFits(os.path.join(inputdir, name))
        fn = infile[i].split()[0]
        print(expPath)
        name = "bs_" + str(fn[3:len(fn)-7])
        print(name)
        olddata, newheader = fits.getdata(os.path.join(expPath, fn[3:len(fn)-7]+'.gz'), header=True)
        
        newheader.remove('ESO DET CHIP PXSPACE')
        
        hdu = fits.PrimaryHDU(fimg)
        hdu.header = newheader
        hdu.writeto(os.path.join(inputdir,name), clobber = True)
        print("done")

        if (i >= hwid and i < (len(data)-hwid - 1)):
            skybeg += 1




if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description="Calculate mean plane of the images")
    parser.add_argument("--inputdir", default=".", help="Input directory")
    parser.add_argument('--flat', action="store_true", help="Calculate flat?", default=False)
    parser.add_argument('--gain', action="store_true", help="Creat gain?", default=False)
    parser.add_argument('--pointing', default="D", help="Pointing")
    args = parser.parse_args()
    
    inputdir = args.inputdir
    
    calibdir = inputdir[0:14] + "/CALIBS/"
    pointing = args.pointing
    
    scale, data, bkgs, bkgmean = scale(inputdir,pointing)
    
    if args.flat:
        
        plane = medplane(data, scale)
        hdu = fits.PrimaryHDU(plane)
        hdu.writeto(os.path.join(calibdir,'flat'+pointing+'.fits'), clobber = True)
        print("created flat")
    
    hdu = fits.PrimaryHDU(bkgmean)
    hdu.writeto(os.path.join(calibdir,'bkgmean'+pointing+'.fits'), clobber = True)


    flatExp = afwImage.ExposureF(os.path.join(calibdir,'flat'+pointing+'.fits'))
    flat = flatExp.getMaskedImage()
    if args.gain:
        gain = gainMap(flatExp, calibdir, pointing)
        print("created gainmap")
    else:
        exp1 = afwImage.ExposureF(os.path.join(calibdir,'gain'+pointing+'.fits'))
        gain = exp1.getMaskedImage().getImage().getArray()



    skyfilter(data, bkgs, gain, inputdir, pointing)
