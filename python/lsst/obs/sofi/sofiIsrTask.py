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
import lsst.ip.isr as ip_isr
import lsst.pipe.base as pipe_base

class MonocamIsrTask(ip_isr.IsrTask):
    @pipe_base.timeMethod
    def run(self, ccdExposure, bias=None, dark=None,  flat=None, defects=None, fringes=None, bfKernel=None):
        """!Perform instrument signature removal on an exposure

        Steps include:
        - Detect saturation, apply overscan correction, bias, dark and flat
        - Perform CCD assembly
        - Interpolate over defects, saturated pixels and all NaNs

        \param[in] ccdExposure  -- lsst.afw.image.exposure of detector data
        \param[in] bias -- exposure of bias frame
        \param[in] dark -- exposure of dark frame
        \param[in] flat -- exposure of flatfield
        \param[in] defects -- list of detects
        \param[in] fringes -- a pipe_base.Struct with field fringes containing
                              exposure of fringe frame or list of fringe exposure
        \param[in] bfKernel -- kernel for brighter-fatter correction

        \return a pipe_base.Struct with field:
         - exposure
        """

        #Validate Input
        if self.config.doBias and bias is None:
            raise RuntimeError("Must supply a bias exposure if config.doBias True")
        if self.config.doDark and dark is None:
            raise RuntimeError("Must supply a dark exposure if config.doDark True")
        if self.config.doFlat and flat is None:
            raise RuntimeError("Must supply a flat exposure if config.doFlat True")
        if self.config.doBrighterFatter and bfKernel is None:
            raise RuntimeError("Must supply a kernel if config.doBrighterFatter True")
        if fringes is None:
            fringes = pipe_base.Struct(fringes=None)
        if self.config.doFringe and not isinstance(fringes, pipe_base.Struct):
            raise RuntimeError("Must supply fringe exposure as a pipe_base.Struct")

        defects = [] if defects is None else defects

        ccd = ccdExposure.getDetector()

        if self.config.doBias:
            self.biasCorrection(ccdExposure, bias)

        if self.config.doBrighterFatter:
            self.brighterFatterCorrection(ccdExposure, bfKernel,
                                          self.config.brighterFatterMaxIter,
                                          self.config.brighterFatterThreshold,
                                          self.config.brighterFatterApplyGain,
                                          )

        if self.config.doDark:
            self.darkCorrection(ccdExposure, dark)

        for amp in ccd:
            #if ccdExposure is one amp, check for coverage to prevent performing ops multiple times
            if ccdExposure.getBBox().contains(amp.getBBox()):
                ampExposure = ccdExposure.Factory(ccdExposure, amp.getBBox())
                self.updateVariance(ampExposure, amp)

        if self.config.doFringe and not self.config.fringeAfterFlat:
            self.fringe.run(ccdExposure, **fringes.getDict())

        if self.config.doFlat:
            self.flatCorrection(ccdExposure, flat)

        self.maskAndInterpDefect(ccdExposure, defects)

        self.saturationInterpolation(ccdExposure)

        self.maskAndInterpNan(ccdExposure)

        if self.config.doFringe and self.config.fringeAfterFlat:
            self.fringe.run(ccdExposure, **fringes.getDict())

        ccdExposure.getCalib().setFluxMag0(self.config.fluxMag0T1 * ccdExposure.getCalib().getExptime())

        return pipe_base.Struct(
            exposure = ccdExposure,
        )

    @pipe_base.timeMethod
    def runDataRef(self, sensorRef):
        """!Perform instrument signature removal on a ButlerDataRef of a Sensor

        - Read in necessary detrending/isr/calibration data
        - Process raw exposure in run()
        - Persist the ISR-corrected exposure as "postISRCCD" if config.doWrite is True

        \param[in] sensorRef -- daf.persistence.butlerSubset.ButlerDataRef of the
                                detector data to be processed
        \return a pipe_base.Struct with fields:
        - exposure: the exposure after application of ISR
        """
        self.log.info("Performing ISR on sensor %s" % (sensorRef.dataId))
        # We should probably loop over this using the butler.
        ampDict = {}
        for channel in range(16):
            sensorRef.dataId['channel'] = channel+1 # to get the correct channel
            ampExposure = sensorRef.get('raw_amp', immediate=True)
            ampExposure = self.convertIntToFloat(ampExposure)
            # assumes amps are in order of the channels
            amp = ampExposure.getDetector()[channel]

            self.saturationDetection(ampExposure, amp)
            self.overscanCorrection(ampExposure, amp)
            ampDict[amp.getName()] = ampExposure

        ccdExposure = self.assembleCcd.assembleCcd(ampDict)
        isrData = self.readIsrData(sensorRef, ccdExposure)
        result = self.run(ccdExposure, **isrData.getDict())


        if self.config.doWrite:
            sensorRef.put(result.exposure, "postISRCCD")

        return result

