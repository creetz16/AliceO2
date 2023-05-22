// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_MID_CHANNELCALIBRATORPARAM_H
#define O2_MID_CHANNELCALIBRATORPARAM_H

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2
{
namespace mid
{

/**
 * @class ChannelCalibratorParam
 * @brief Configurable parameters for the Bad Channel Calibrator
 */
struct ChannelCalibratorParam : public o2::conf::ConfigurableParamHelper<ChannelCalibratorParam> {

  float maxNoise = 10000.f;                  ///< maximum allowed noise value (Hz)
  float maxDead = 0.9f;                      ///< maximum fraction of time a strip was not responding to FET
  unsigned long int nCalibTriggers = 115000; ///< Number of calibration triggers before sending
  bool onlyAtEndOfStream = false;            ///< Run only at end of stream

  O2ParamDef(ChannelCalibratorParam, "MIDChannelCalibratorParam");
};
} // namespace mid
} // namespace o2

#endif
