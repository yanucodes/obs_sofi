import os
import glob
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser(description="Add extenstion type 'IMAGE' to fits files")
parser.add_argument("--inputdir", default=".", help="Input directory")
args = parser.parse_args()
inputdir = args.inputdir

list = glob.glob(os.path.join(inputdir, "*.fits")) + glob.glob(os.path.join(inputdir, "*.fits.gz")) + glob.glob(os.path.join(inputdir, "*.fit.gz")) + glob.glob(os.path.join(inputdir, "*.fit"))

for fitsfile in list:
    f = fits.open(fitsfile, mode='update')
    f[0].header.set('EXTTYPE', 'IMAGE   ')
    f.flush()
    f.close()