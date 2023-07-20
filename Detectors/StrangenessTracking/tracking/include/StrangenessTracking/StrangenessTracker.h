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

/// \file StrangenessTracker.h
/// \brief
///

#ifndef _ALICEO2_STRANGENESS_TRACKER_
#define _ALICEO2_STRANGENESS_TRACKER_

#include <gsl/gsl>

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "StrangenessTracking/IndexTableUtils.h"
#include "StrangenessTracking/StrangenessTrackingConfigParam.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/StrangeTrack.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"

#include "DataFormatsITS/TrackITS.h"
#include "ITSBase/GeometryTGeo.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"

#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"

#ifndef HomogeneousField
#define HomogeneousField
#endif

/// includes KFParticle
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

namespace o2
{
namespace strangeness_tracking
{

enum DauType : int {
  kV0DauPos = 0,
  kV0DauNeg = 1,
  kBach = 2
};

struct ClusAttachments {

  std::array<unsigned int, 7> arr;
};

class StrangenessTracker
{
 public:
  using StrangeTrack = o2::dataformats::StrangeTrack;
  using PID = o2::track::PID;
  using TrackITS = o2::its::TrackITS;
  using ITSCluster = o2::BaseCluster<float>;
  using V0 = o2::dataformats::V0;
  using Cascade = o2::dataformats::Cascade;
  using GIndex = o2::dataformats::VtxTrackIndex;
  using DCAFitter2 = o2::vertexing::DCAFitterN<2>;
  using DCAFitter3 = o2::vertexing::DCAFitterN<3>;
  using MCLabContCl = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
  using MCLabSpan = gsl::span<const o2::MCCompLabel>;
  using VBracket = o2::math_utils::Bracket<int>;

  StrangenessTracker() = default;
  ~StrangenessTracker() = default;

  bool loadData(const o2::globaltracking::RecoContainer& recoData);
  bool matchDecayToITStrack(float decayR);
  void prepareITStracks();
  void process();
  bool updateTrack(const ITSCluster& clus, o2::track::TrackParCov& track);

  std::vector<ClusAttachments>& getClusAttachments() { return mClusAttachments; };
  std::vector<StrangeTrack>& getStrangeTrackVec() { return mStrangeTrackVec; };
  std::vector<o2::MCCompLabel>& getStrangeTrackLabels() { return mStrangeTrackLabels; };

  float getBz() const { return mBz; }
  void setBz(float d) { mBz = d; }
  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDict = d; }
  void setCorrType(const o2::base::PropagatorImpl<float>::MatCorrType& type) { mCorrType = type; }
  void setConfigParams(const StrangenessTrackingParamConfig* params) { mStrParams = params; }
  void setMCTruthOn(bool v) { mMCTruthON = v; }

  void clear()
  {
    mDaughterTracks.clear();
    mClusAttachments.clear();
    mStrangeTrackVec.clear();
    mTracksIdxTable.clear();
    mSortedITStracks.clear();
    mSortedITSindexes.clear();
    mITSvtxBrackets.clear();
    mInputITSclusters.clear();
    mInputClusterSizes.clear();
    if (mMCTruthON) {
      mStrangeTrackLabels.clear();
    }
  }

  void setupFitters()
  {
    mFitterV0.setBz(mBz);
    mFitter3Body.setBz(mBz);
    mFitterV0.setUseAbsDCA(true);
    mFitter3Body.setUseAbsDCA(true);
  }

  void setPID()
  {
    if (mStrParams->pidV0 == 0) {
      pidV0 = PID::HyperTriton;
      pidV0comp = PID::Hyperhydrog4;
    } else if (mStrParams->pidV0 == 1) {
      pidV0 = PID::Hyperhydrog4;
      pidV0comp = PID::HyperTriton;
    }
    if (mStrParams->pidCasc == 0) {
      pidCasc = PID::XiMinus;
      pidCascComp = PID::OmegaMinus;
    } else if (mStrParams->pidCasc == 1) {
      pidCasc = PID::OmegaMinus;
      pidCascComp = PID::XiMinus;
    }
  }

  double calcV0alpha(const V0& v0)
  {
    std::array<float, 3> momT, momP, momN;
    v0.getProng(0).getPxPyPzGlo(momP);
    v0.getProng(1).getPxPyPzGlo(momN);
    v0.getPxPyPzGlo(momT);
    float qNeg = momN[0] * momT[0] + momN[1] * momT[1] + momN[2] * momT[2];
    float qPos = momP[0] * momT[0] + momP[1] * momT[1] + momP[2] * momT[2];
    return (qPos - qNeg) / (qPos + qNeg);
  };

  double calcMotherMass(const std::array<float, 3>& pDauFirst, const std::array<float, 3>& pDauSecond, PID pidDauFirst, PID pidDauSecond)
  {
    double m2DauFirst = PID::getMass2(pidDauFirst);
    double m2DauSecond = PID::getMass2(pidDauSecond);
    double p2DauFirst = (pDauFirst[0] * pDauFirst[0]) + (pDauFirst[1] * pDauFirst[1]) + (pDauFirst[2] * pDauFirst[2]);
    double p2DauSecond = (pDauSecond[0] * pDauSecond[0]) + (pDauSecond[1] * pDauSecond[1]) + (pDauSecond[2] * pDauSecond[2]);
    float ePos = std::sqrt(p2DauFirst + m2DauFirst), eNeg = std::sqrt(p2DauSecond + m2DauSecond);

    double e2Mother = (ePos + eNeg) * (ePos + eNeg);
    double pxMother = (pDauFirst[0] + pDauSecond[0]);
    double pyMother = (pDauFirst[1] + pDauSecond[1]);
    double pzMother = (pDauFirst[2] + pDauSecond[2]);
    double p2Mother = (pxMother * pxMother) + (pyMother * pyMother) + (pzMother * pzMother);
    return std::sqrt(e2Mother - p2Mother);
  }

  bool createKFV0(const o2::track::TrackParCov& posTrack, const o2::track::TrackParCov& negTrack, PID pidMother)
  {

    float massPosDaughter, massNegDaughter;
    if (posTrack.getAbsCharge() == 2) { // if charge of positive daughter is two, it is 3He or 4He
      (pidMother == PID::HyperTriton) ? massPosDaughter = o2::constants::physics::MassHelium3 : massPosDaughter = o2::constants::physics::MassAlpha;
      massNegDaughter = o2::constants::physics::MassPionCharged;
    } else {
      massPosDaughter = o2::constants::physics::MassPionCharged;
      (pidMother == PID::HyperTriton) ? massNegDaughter = o2::constants::physics::MassHelium3 : massNegDaughter = o2::constants::physics::MassAlpha;
    }

    int nV0Daughters = 2;
    // create KFParticle objects from trackParCovs
    KFParticle kfpDauPos = createKFParticleFromTrackParCov(posTrack, posTrack.getCharge(), massPosDaughter);
    KFParticle kfpDauNeg = createKFParticleFromTrackParCov(negTrack, negTrack.getCharge(), massNegDaughter);
    const KFParticle* V0Daughters[2] = {&kfpDauPos, &kfpDauNeg};

    // construct mother
    KFParticle KFV0;
    KFV0.SetConstructMethod(mStrParams->kfConstructMethod);
    try {
      KFV0.Construct(V0Daughters, nV0Daughters);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to create KFParticle V0 from daughter tracks." << e.what();
      return false;
    }
    LOG(debug) << "KFParticle V0 position before TransportToDecayVertex(): (" <<  KFV0.GetX() << ", " << KFV0.GetY() << ", " << KFV0.GetZ() << ")";
    KFV0.TransportToDecayVertex();
    LOG(debug) << "KFParticle V0 position after TransportToDecayVertex(): (" <<  KFV0.GetX() << ", " << KFV0.GetY() << ", " << KFV0.GetZ() << ")";
    kfpMother = KFV0;

    return true;
  };

  bool createKFCascade(const o2::track::TrackParCov& posTrack, const o2::track::TrackParCov& negTrack, const o2::track::TrackParCov& bachTrack, PID pidMother)
    {

      float massPosDaughter, massNegDaughter;
      if (bachTrack.getCharge() < 0) { // if bachelor charge negative, cascade is a Xi- --> Lam + pi- --> (p+ + pi-) + pi-   OR   Omega- --> Lam + K- --> (p+ + pi-) + K-
        massPosDaughter = o2::constants::physics::MassProton;
        massNegDaughter = o2::constants::physics::MassPionCharged;
      }  else { // if bachelor charge positive, cascade is a Xi+ --> Lam + pi+ --> (p- + pi+) + pi+  OR   Omega+ --> Lam + K+ --> (p- + pi+) + K+
        massPosDaughter = o2::constants::physics::MassPionCharged;
        massNegDaughter = o2::constants::physics::MassProton;
      }

      int nV0Daughters = 2;
      // create KFParticle objects from trackParCovs
      KFParticle kfpDauPos = createKFParticleFromTrackParCov(posTrack, posTrack.getCharge(), massPosDaughter);
      KFParticle kfpDauNeg = createKFParticleFromTrackParCov(negTrack, negTrack.getCharge(), massNegDaughter);
      const KFParticle* V0Daughters[2] = {&kfpDauPos, &kfpDauNeg};

      // construct V0
      KFParticle KFV0;
      KFV0.SetConstructMethod(mStrParams->kfConstructMethod);
      try {
        KFV0.Construct(V0Daughters, nV0Daughters);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct cascade V0 from daughter tracks: " << e.what();
        return false;
      }
      if (mStrParams->kfMassConst) {
        KFV0.SetNonlinearMassConstraint(o2::constants::physics::MassLambda);
      }

      int nCascDaughters = 2;
      KFParticle kfpBach;
      // create KFParticle objects from trackParCovs
      if (pidMother == PID::XiMinus) {
        kfpBach = createKFParticleFromTrackParCov(bachTrack, bachTrack.getCharge(), o2::constants::physics::MassPionCharged);
      } 
      else if (pidMother == PID::OmegaMinus) {
        kfpBach = createKFParticleFromTrackParCov(bachTrack, bachTrack.getCharge(), o2::constants::physics::MassKaonCharged);
      }
      const KFParticle* CascDaugthers[2] = {&kfpBach, &KFV0};

      // construct mother
      KFParticle KFCasc;
      KFCasc.SetConstructMethod(mStrParams->kfConstructMethod);
      try {
        KFCasc.Construct(CascDaugthers, nCascDaughters);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct cascade from V0 and bachelor track: " << e.what();
        return false;
      }
      KFCasc.TransportToDecayVertex();
      kfpMother = KFCasc;

      return true;
    };

  std::vector<ITSCluster> getTrackClusters()
  {
    std::vector<ITSCluster> outVec;
    auto firstClus = mITStrack.getFirstClusterEntry();
    auto ncl = mITStrack.getNumberOfClusters();
    for (int icl = 0; icl < ncl; icl++) {
      outVec.push_back(mInputITSclusters[mInputITSidxs[firstClus + icl]]);
    }
    return outVec;
  };

  std::vector<int> getTrackClusterSizes()
  {
    std::vector<int> outVec;
    auto firstClus = mITStrack.getFirstClusterEntry();
    auto ncl = mITStrack.getNumberOfClusters();
    for (int icl = 0; icl < ncl; icl++) {
      outVec.push_back(mInputClusterSizes[mInputITSidxs[firstClus + icl]]);
    }
    return outVec;
  };

  void getClusterSizes(std::vector<int>& clusSizeVec, const gsl::span<const o2::itsmft::CompClusterExt> ITSclus, gsl::span<const unsigned char>::iterator& pattIt, const o2::itsmft::TopologyDictionary* mdict)
  {
    for (unsigned int iClus{0}; iClus < ITSclus.size(); ++iClus) {
      auto& clus = ITSclus[iClus];
      auto pattID = clus.getPatternID();
      int npix;
      o2::itsmft::ClusterPattern patt;

      if (pattID == o2::itsmft::CompCluster::InvalidPatternID || mdict->isGroup(pattID)) {
        patt.acquirePattern(pattIt);
        npix = patt.getNPixels();
      } else {

        npix = mdict->getNpixels(pattID);
        patt = mdict->getPattern(pattID);
      }
      clusSizeVec[iClus] = npix;
    }
    // LOG(info) << " Patt Npixel: " << pattVec[0].getNPixels();
  }

  float getMatchingChi2(o2::track::TrackParCovF v0, const TrackITS& itsTrack)
  {
    if (v0.rotate(itsTrack.getParamOut().getAlpha()) && v0.propagateTo(itsTrack.getParamOut().getX(), mBz)) {
      return v0.getPredictedChi2(itsTrack.getParamOut());
    }
    return -100;
  };

  o2::MCCompLabel getStrangeTrackLabel() // ITS label with fake flag recomputed
  {
    bool isFake = false;
    auto itsTrkLab = mITSTrkLabels[mStrangeTrack.mITSRef];
    for (unsigned int iLay = 0; iLay < 7; iLay++) {
      if (mITStrack.hasHitOnLayer(iLay) && mITStrack.isFakeOnLayer(iLay) && mStructClus.arr[iLay] == 0) {
        isFake = true;
        break;
      }
    }
    itsTrkLab.setFakeFlag(isFake);
    return itsTrkLab;
  }

  template <typename T>
  KFParticle createKFParticleFromTrackParCov(const o2::track::TrackParametrizationWithError<T>& trackparCov, int charge, float mass) 
  {
    std::array<T, 3> xyz, pxpypz;
    float xyzpxpypz[6];
    trackparCov.getPxPyPzGlo(pxpypz);
    trackparCov.getXYZGlo(xyz);
    for (int i{0}; i < 3; ++i) {
      xyzpxpypz[i] = xyz[i];
      xyzpxpypz[i + 3] = pxpypz[i];
    }

    LOG(debug) << "Glo position of TrackParCov before KFParticle creation: (" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ")";

    std::array<float, 21> cv;
    try {
      trackparCov.getCovXYZPxPyPzGlo(cv);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to get cov matrix from TrackParCov" << e.what();
    }

    KFParticle kfPart;
    float Mini, SigmaMini, M, SigmaM;
    kfPart.GetMass(Mini, SigmaMini);
    LOG(debug) << "Daughter KFParticle mass before creation: " << Mini << " +- " << SigmaMini;
    
    try {
      kfPart.Create(xyzpxpypz, cv.data(), charge, mass);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to create KFParticle from daughter TrackParCov" << e.what();
    }

    kfPart.GetMass(M, SigmaM);
    LOG(debug) << "Daughter KFParticle mass after creation: " << M << " +- " << SigmaM;

    LOG(debug) << "Position of KFParticle after creation from TrackParCov: (" << kfPart.GetX() << ", " << kfPart.GetY() << ", " << kfPart.GetZ() << ")";

    return kfPart;
  }

  bool getTrackParCovFromKFP(const KFParticle& kfParticle, const PID pid, const int sign, o2::track::TrackParCovF& track)
  {

    // position check
    LOG(debug) << "Position of kfmother before transformation: (" << kfParticle.GetX() << ", " << kfParticle.GetY() << ", " << kfParticle.GetZ() << ")";
    // LOG(info) << "Momentum of kfmother before transformation: (" << kfParticle.GetPx() << ", " << kfParticle.GetPy() << ", " << kfParticle.GetPz() << ")";
    
    o2::gpu::gpustd::array<float, 3> xyz, pxpypz;
    o2::gpu::gpustd::array<float, 21> cv;

    // get parameters from kfParticle
    xyz[0] = kfParticle.GetX();
    xyz[1] = kfParticle.GetY();
    xyz[2] = kfParticle.GetZ();
    pxpypz[0] = kfParticle.GetPx();
    pxpypz[1] = kfParticle.GetPy();
    pxpypz[2] = kfParticle.GetPz();
    
    // set covariance matrix elements (lower triangle)
    for (int i =0; i < 21; i++) {
      cv[i] = kfParticle.GetCovariance(i);
    }

    // create TrackParCov track
    track = o2::track::TrackParCovF(xyz, pxpypz, cv, sign, true, pid);

    // position check
    std::array<float, 3> xyzGlo, pxpypzGlo;
    track.getXYZGlo(xyzGlo);
    track.getPxPyPzGlo(pxpypzGlo);
    LOG(debug) << "Position of TrackParCov after transformation: (" << track.getX() << ", " << track.getY() << ", " << track.getZ() << ")";
    LOG(debug) << "Glo position of TrackParCov after transformation: (" << xyzGlo[0] << ", " << xyzGlo[1] << ", " << xyzGlo[2] << ")";
    // LOG(info) << "Glo momentum of TrackParCov after transformation: (" << pxpypzGlo[0] << ", " << pxpypzGlo[1] << ", " << pxpypzGlo[2] << ")";

    return true;
  }

  bool propagateToOuterParam(const TrackITS& itsTrack, o2::track::TrackParCovF& track)
  {
    auto geom = o2::its::GeometryTGeo::Instance();
    auto propInstance = o2::base::Propagator::Instance();
    float x = itsTrack.getParamOut().getX();
    int layer{geom->getLayer(itsTrack.getFirstClusterLayer())};
    if (!track.rotate(itsTrack.getParamOut().getAlpha())) {
        return false;
    }

    // transport to OuterParam
    if (!propInstance->propagateToX(track, x, getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mCorrType)) {
      return false;
    }

    // apply material correction 
    if (mCorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE) {
      float thick = layer < 3 ? 0.005 : 0.01;
      constexpr float radl = 9.36f; // Radiation length of Si [cm]
      constexpr float rho = 2.33f;  // Density of Si [g/cm^3]
      if (!track.correctForMaterial(thick, thick * rho * radl)) {
        return false;
      }
    }

    return true;
  }

  bool fillKFinfo(const KFParticle& kfParticle) {
    float M, SigmaM;
    kfParticle.GetMass(M, SigmaM);
    mStrangeTrack.mMasses[0] = M;
    mStrangeTrack.mDecayMom[0] = kfParticle.GetPx();
    mStrangeTrack.mDecayMom[1] = kfParticle.GetPy();
    mStrangeTrack.mDecayMom[2] = kfParticle.GetPz();
    mStrangeTrack.mDecayPt = kfParticle.GetPt();
    mStrangeTrack.mDecayVtx[0] = kfParticle.GetX();
    mStrangeTrack.mDecayVtx[1] = kfParticle.GetY();
    mStrangeTrack.mDecayVtx[2] = kfParticle.GetZ();
    mStrangeTrack.mGeoChi2 = kfParticle.GetChi2();
    return true;
  }


 protected: 
  bool mMCTruthON = false;                      /// flag availability of MC truth
  gsl::span<const TrackITS> mInputITStracks;    // input ITS tracks
  std::vector<VBracket> mITSvtxBrackets;        // time brackets for ITS tracks
  std::vector<int> mTracksIdxTable;             // index table for ITS tracks
  std::vector<int> mInputClusterSizes;          // input cluster sizes
  std::vector<ITSCluster> mInputITSclusters;    // input ITS clusters
  gsl::span<const int> mInputITSidxs;           // input ITS track-cluster indexes
  gsl::span<const V0> mInputV0tracks;           // input V0 of decay daughters
  gsl::span<const Cascade> mInputCascadeTracks; // input V0 of decay daughters
  const MCLabContCl* mITSClsLabels = nullptr;   /// input ITS Cluster MC labels
  MCLabSpan mITSTrkLabels;                      /// input ITS Track MC labels

  std::vector<o2::its::TrackITS> mSortedITStracks; // sorted ITS tracks
  std::vector<int> mSortedITSindexes;              // indexes of sorted ITS tracks
  IndexTableUtils mUtils;                          // structure for computing eta/phi matching selections

  std::vector<StrangeTrack> mStrangeTrackVec;       // structure containing updated mother and daughter tracks
  std::vector<ClusAttachments> mClusAttachments;    // # of attached tracks, -1 not attached, 0 for the mother, > 0 for the daughters
  std::vector<o2::MCCompLabel> mStrangeTrackLabels; // vector of MC labels for mother track

  const StrangenessTrackingParamConfig* mStrParams = nullptr;
  float mBz = -5; // Magnetic field
  const o2::itsmft::TopologyDictionary* mDict = nullptr;

  DCAFitter2 mFitterV0;    // optional DCA Fitter for recreating V0 with hypertriton mass hypothesis
  DCAFitter3 mFitter3Body; // optional DCA Fitter for final 3 Body refit

  o2::base::PropagatorImpl<float>::MatCorrType mCorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrNONE; // use mat correction

  std::vector<o2::track::TrackParCovF> mDaughterTracks;  // vector of daughter tracks
  StrangeTrack mStrangeTrack;                            // structure containing updated mother and daughter track refs
  ClusAttachments mStructClus;                           // # of attached tracks, 1 for mother, 2 for daughter
  o2::its::TrackITS mITStrack;                           // ITS track
  std::array<GIndex, 2> mV0dauIDs;                       // V0 daughter IDs
  KFParticle kfpMother;                                  // mother KFParticle
  o2::track::PID pidV0;                                  // PID hypothesis for the V0 fitting
  o2::track::PID pidV0comp;                                  // commpeting PID hypothesis for the V0 fitting
  o2::track::PID pidCasc;                                // PID hypothesis for the cascade fitting
  o2::track::PID pidCascComp;                                // competing PID hypothesis for the cascade fitting

  ClassDefNV(StrangenessTracker, 1);
};

} // namespace strangeness_tracking
} // namespace o2

#endif //  _ALICEO2_STRANGENESS_TRACKER_