import numpy as np
import lsst.afw.image as afwImage

def fit_parabola(arr, ncc, x):
    d = arr[x+1] - (2.0*arr[x]+arr[x-1])
    center = float(x) + 0.5 - (arr[x+1] - arr[x]) / d
    return center

def find_offsets(coadd, exposure):
    
    print " "
    print " "
    print "Find offsets"
    print " "
    
    nx = 1024
    ny = 1024
    
    nlist = nx*ny
    
    x = {}
    y = {}
    
    for i in range(nlist):
        x[i] = i%nx
        y[i] = i/nx
    
    ixoff, iyoff = exposure.getWcs().skyToPixel(coadd.getWcs().getSkyOrigin()) - coadd.getWcs().getPixelOrigin()
    
    img = exposure.getMaskedImage().getImage().getArray().ravel()
    p = coadd.getMaskedImage().getImage().getArray().ravel()

    hwid = 40
    ncc = 2*hwid + 1

    ccimg = np.zeros((ncc,ncc))
    
    for i in range(nlist):
        ipix = int((y[i]+iyoff)*nx + x[i] + ixoff)
        for k in range(-hwid,hwid+1):
            m = ipix + k*nx
            for l in range(-hwid,hwid+1):
                n = m + l
                if (n>=0 and n<nlist):
                    ccimg[k+hwid][l+hwid]+=p[i]*img[n]
    
    """
    for j in range(nx):
        for i in range(ny):
            xcoff = i + ixoff
            ycoff = j + iyoff
            for l in range(-hwid,hwid):
                for k in range(-hwid,hwid):
                    m = xcoff + k
                    n = ycoff + l
                    if ( m>=0 and n>=0 and m<nx and n<ny):
                        ccimg[k+hwid][l+hwid]+=p[i][j]*img[m][n]
"""
    ix0, iy0 = np.where(ccimg == ccimg.max())

    #check whether ix0, iy0 are on the edge

    xprof = ccimg[ix0,:].ravel()
    yprof = ccimg[:,iy0].ravel()

    print ccimg
    print ix0, iy0
    print xprof
    print yprof

    xoff = float(fit_parabola(xprof, ncc, iy0))
    yoff = float(fit_parabola(yprof, ncc, ix0))
    
    print "xoff, yoff"
    print xoff, yoff

    return xoff, yoff