import lsst.afw.image as afwImage


im1 = afwImage.ExposureF("dd_F02_S22_10_022.fits")
im2 = afwImage.ExposureF("fpC-004192-g4-0300.fit.gz")

print im1.getMaskedImage().getVariance().getArray()
print im2.getMaskedImage().getVariance().getArray()