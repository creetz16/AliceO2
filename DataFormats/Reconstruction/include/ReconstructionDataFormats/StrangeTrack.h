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

/// \file StrangeTrack.h
/// \brief
///

#ifndef _ALICEO2_STRANGETRACK_
#define _ALICEO2_STRANGETRACK_

#include <array>
#include "ReconstructionDataFormats/Track.h"

namespace o2
{
namespace dataformats
{

enum kPartType { kStrkV0,
                 kStrkCascade,
                 kStrkThreeBody };

struct StrangeTrack {
  kPartType mPartType;
  o2::track::TrackParCovF mMother;
  unsigned int mITSRef = -1;
  unsigned int mDecayRef = -1;
  std::array<float, 3> mDecayVtx;
  std::array<float, 3> mDecayMom;
  std::array<float, 2> mMasses; // V0: hypertriton and hyperhydrongen4, cascade: Xi and Omega.
  float mITSClusSize;
  float mMatchChi2;
  float mTopoChi2;

  std::array<float, 3> mIUMom;
  float mDecayPt;
  float mGeoChi2;
  float mGeoChi2KFcreation;
  float mMotherZ;
  float mMotherZkf;
  float mMotherR;
  float mMotherRkf;
  float mMotherErrZ;
  float mMotherErrZkf;

  float mMassCascV0;
  float mChargeCascV0;
  float mMassInit;
  float mErrX;
  float mErrY;
  float mErrZ;
  float mErrPx;
  float mErrPy;
  float mErrPz;
  float mErrXV0;
  float mErrYV0;
  float mErrZV0;
  float mErrPxV0;
  float mErrPyV0;
  float mErrPzV0;
};

} // namespace dataformats
} // namespace o2

#endif // _ALICEO2_STRANGETRACK_
