#taken from obs_monocam, adapted for obs_sofi
import os
import re
from lsst.pipe.tasks.ingest import IngestTask, ParseTask, RegisterTask, assertCanCopy
from lsst.pipe.tasks.ingestCalibs import CalibsParseTask

# Lookup table for lab data
# XXX This is completely made up and needs to be updated with the real values when they become available.
filters = {
    0: 'Ks',
    "UNK": "UNKNOWN",
}

EXTENSIONS = ["fits", "gz"]  # Filename extensions to strip off

class SofiParseTask(ParseTask):
    """Parser suitable for lab data"""
    _counter = -1  # Visit counter; negative so it doesn't overlap with the 'id' field, which autoincrements

    def getInfo(self, filename):
        # Grab the basename
        phuInfo, infoList = ParseTask.getInfo(self, filename)
        basename = os.path.basename(filename)
        while any(basename.endswith("." + ext) for ext in EXTENSIONS):
            basename = basename[:basename.rfind('.')]
        phuInfo['basename'] = basename
        return phuInfo, infoList

    def assign_expNum(self, md):
        expNum = self._counter
        self._counter -= 1
        return expNum

    def translate_ccd(self, md):
        return 0  # There's only one

    def translate_filter(self, md):
        return filters[md.get("FILTER")]


class SofiRegisterTask(RegisterTask):
    """Put data in the registry
    """
    def addExpNums(self, conn, dryrun=False, table=None):
        """Set the exposure numbers to match the 'id' field"""
        if table is None:
            table = self.config.table
        sql = "UPDATE %s SET expNum = id" % table
        if dryrun:
            print "Would execute: %s" % sql
        else:
            conn.execute(sql)
        return RegisterTask.addVisits(self, conn, dryrun=dryrun, table=table)


class SofiIngestTask(IngestTask):
    """Ingest data
    """
    def run(self, args):
        """Open the database"""
        #getDatabase(args.butler.repository._mapper.root)
        #(see lsst.obs.monocam.hack.getDatabase)
        return IngestTask.run(self, args)


##############################################################################################################

class SofiCalibsParseTask(CalibsParseTask):
    """Parser for calibs"""
    def _translateFromCalibId(self, field, md):
        """Get a value from the CALIB_ID written by constructCalibs"""
        data = md.get("CALIB_ID")
        match = re.search(".*%s=(\S+)" % field, data)
        return match.groups()[0]

    def translate_ccd(self, md):
        return self._translateFromCalibId("ccd", md)

    def translate_filter(self, md):
        return self._translateFromCalibId("filter", md)
