#!/usr/bin/env python2
# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012, 2013 LSST Corporation.
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
import argparse
import glob
import os
import re
try:
    import sqlite3
except ImportError:
    # try external pysqlite package; deprecated
    import sqlite as sqlite3
import sys
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage

parser = argparse.ArgumentParser()
parser.add_argument("--create", default=False, action="store_true", help="Create new registry (clobber old)?")
parser.add_argument("--root", default=".", help="Root directory")
args = parser.parse_args()

root = args.root
files = glob.glob(os.path.join(root, "*.fits.gz"))
sys.stderr.write('processing %d files...\n' % (len(files)))

registryName = "calibregistry.sqlite3"
if os.path.exists(registryName) and args.create:
    os.unlink(registryName)

makeTables = not os.path.exists(registryName)
conn = sqlite.connect(registryName)
if makeTables:
    cmd = "create table calib (id integer primary key autoincrement"
    cmd += ", mjd float, filename text, calibtype text)"
    conn.execute(cmd)
    conn.commit()

for fits in files:
    matches = re.search(r'calib/(DARK|FLAT)(\w{5}|\d{2}|\d{1}\.d{1}).fits.gz', fits)
    if not matches:
        print >>sys.stderr, "Warning: skipping unrecognized filename:", fits
        continue

    sys.stderr.write("Processing %s\n" % (fits))
    
    # Extract information from header
    im = afwImage.ExposureF(fits)
    h = im.getMetadata()
    mjd = h.get('MJD-OBS')
    filename = h.get('ARCFILE')
    objtype = h.get('OBJECT')
    
    try:
        conn.execute("INSERT INTO raw VALUES (NULL, ?, ?, ?)",
                     (mjd, filename, objtype))
    
    except Exception, e:
        print "skipping botched %s: %s" % (fits, e)
        continue

conn.commit()
conn.close()
