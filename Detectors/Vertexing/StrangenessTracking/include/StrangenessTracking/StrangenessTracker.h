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
#include "ReconstructionDataFormats/Decay3Body.h"
#include "ReconstructionDataFormats/DecayNBodyIndex.h"
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
  using V0Index = o2::dataformats::V0Index;
  using Cascade = o2::dataformats::Cascade;
  using CascadeIndex = o2::dataformats::CascadeIndex;
  using Decay3Body = o2::dataformats::Decay3Body;
  using Decay3BodyIndex = o2::dataformats::Decay3BodyIndex;
  using GIndex = o2::dataformats::VtxTrackIndex;
  using DCAFitter2 = o2::vertexing::DCAFitterN<2>;
  using DCAFitter3 = o2::vertexing::DCAFitterN<3>;
  using DCAFitter4 = o2::vertexing::DCAFitterN<4>;
  using MCLabContCl = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
  using MCLabSpan = gsl::span<const o2::MCCompLabel>;
  using VBracket = o2::math_utils::Bracket<int>;

  StrangenessTracker() = default;
  ~StrangenessTracker() = default;

  bool loadData(const o2::globaltracking::RecoContainer& recoData);
  bool matchDecayToITStrack(float decayR, StrangeTrack& strangeTrack, ClusAttachments& structClus, const TrackITS& itsTrack, std::vector<o2::track::TrackParCovF>& daughterTracks, int iThread = 0);
  void prepareITStracks();
  void process();
  void processV0(int iv0, const V0& v0, const V0Index& v0Idx, int iThread = 0);
  void processCascade(int icasc, const Cascade& casc, const CascadeIndex& cascIdx, const V0& cascV0, int iThread = 0);
  void process3Body(int i3body, const Decay3Body& dec3body, const Decay3BodyIndex& dec3bodyIdx, int iThread = 0);
  bool updateTrack(const ITSCluster& clus, o2::track::TrackParCov& track);

  std::vector<ClusAttachments>& getClusAttachments(int iThread = 0) { return mClusAttachments[iThread]; };
  std::vector<StrangeTrack>& getStrangeTrackVec(int iThread = 0) { return mStrangeTrackVec[iThread]; };
  std::vector<o2::MCCompLabel>& getStrangeTrackLabels(int iThread = 0) { return mStrangeTrackLabels[iThread]; };
  size_t getNTracks(int ithread = 0) const { return ithread < (int)mStrangeTrackVec.size() ? mStrangeTrackVec[ithread].size() : 0; }

  float getBz() const { return mBz; }
  void setBz(float d) { mBz = d; }
  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDict = d; }
  void setCorrType(const o2::base::PropagatorImpl<float>::MatCorrType& type) { mCorrType = type; }
  void setConfigParams(const StrangenessTrackingParamConfig* params) { mStrParams = params; }
  void setMCTruthOn(bool v) { mMCTruthON = v; }
  bool getMCTruthOn() const { return mMCTruthON; }

  void clear()
  {
    for (int i = 0; i < mNThreads; i++) {
      mDaughterTracks[i].clear();
      mClusAttachments[i].clear();
      mStrangeTrackVec[i].clear();
      if (mMCTruthON) {
        mStrangeTrackLabels[i].clear();
      }
    }
    mTracksIdxTable.clear();
    mSortedITStracks.clear();
    mSortedITSindexes.clear();
    mITSvtxBrackets.clear();
    mInputITSclusters.clear();
    mInputClusterSizes.clear();
  }

  void setupThreads(int nThreads = 1)
  {
    mNThreads = nThreads;
    mFitterV0.resize(nThreads);
    mFitter3Body.resize(nThreads);
    mFitter4Body.resize(nThreads);
    mStrangeTrackVec.resize(nThreads);
    mClusAttachments.resize(nThreads);
    mStrangeTrackLabels.resize(nThreads);
    mDaughterTracks.resize(nThreads);
  }

  void setupFitters()
  {
    for (auto& fitter : mFitterV0) {
      fitter.setBz(mBz);
      fitter.setUseAbsDCA(true);
    }
    for (auto& fitter : mFitter3Body) {
      fitter.setBz(mBz);
      fitter.setUseAbsDCA(true);
    }
    for (auto& fitter : mFitter4Body) {
      fitter.setBz(mBz);
      fitter.setUseAbsDCA(true);
    }
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

  double calcMotherMass3body(const std::array<float, 3>& pDauFirst, const std::array<float, 3>& pDauSecond, const std::array<float, 3>& pDauThird, PID pidDauFirst, PID pidDauSecond, PID pidDauThird)
  {
    double m2DauFirst = PID::getMass2(pidDauFirst);
    double m2DauSecond = PID::getMass2(pidDauSecond);
    double m2DauThird = PID::getMass2(pidDauThird);
    double p2DauFirst = (pDauFirst[0] * pDauFirst[0]) + (pDauFirst[1] * pDauFirst[1]) + (pDauFirst[2] * pDauFirst[2]);
    double p2DauSecond = (pDauSecond[0] * pDauSecond[0]) + (pDauSecond[1] * pDauSecond[1]) + (pDauSecond[2] * pDauSecond[2]);
    double p2DauThird = (pDauThird[0] * pDauThird[0]) + (pDauThird[1] * pDauThird[1]) + (pDauThird[2] * pDauThird[2]);
    float eFirst = std::sqrt(p2DauFirst + m2DauFirst), eSecond = std::sqrt(p2DauSecond + m2DauSecond), eThird = std::sqrt(p2DauThird + m2DauThird);

    double e2Mother = (eFirst + eSecond + eThird) * (eFirst + eSecond + eThird);
    double pxMother = (pDauFirst[0] + pDauSecond[0] + pDauThird[0]);
    double pyMother = (pDauFirst[1] + pDauSecond[1] + pDauThird[1]);
    double pzMother = (pDauFirst[2] + pDauSecond[2] + pDauThird[2]);
    double p2Mother = (pxMother * pxMother) + (pyMother * pyMother) + (pzMother * pzMother);
    return std::sqrt(e2Mother - p2Mother);
  }

  bool recreateV0(const o2::track::TrackParCov& posTrack, const o2::track::TrackParCov& negTrack, V0& newV0, int iThread = 0)
  {
    int nCand;
    try {
      nCand = mFitterV0[iThread].process(posTrack, negTrack);
    } catch (std::runtime_error& e) {
      return false;
    }
    if (!nCand || !mFitterV0[iThread].propagateTracksToVertex()) {
      return false;
    }

    const auto& v0XYZ = mFitterV0[iThread].getPCACandidatePos();

    auto& propPos = mFitterV0[iThread].getTrack(0, 0);
    auto& propNeg = mFitterV0[iThread].getTrack(1, 0);

    std::array<float, 3> pP, pN;
    propPos.getPxPyPzGlo(pP);
    propNeg.getPxPyPzGlo(pN);
    std::array<float, 3> pV0 = {pP[0] + pN[0], pP[1] + pN[1], pP[2] + pN[2]};
    newV0 = V0(v0XYZ, pV0, mFitterV0[iThread].calcPCACovMatrixFlat(0), propPos, propNeg, PID::HyperTriton);
    return true;
  };

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

    std::array<float, 21> cv;
    try {
      trackparCov.getCovXYZPxPyPzGlo(cv);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to get cov matrix from TrackParCov" << e.what();
    }
    KFParticle kfPart;
    try {
      kfPart.Create(xyzpxpypz, cv.data(), charge, mass);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to create KFParticle from daughter TrackParCov" << e.what();
    }

    return kfPart;
  }

  bool getTrackParCovFromKFP(const KFParticle& kfParticle, const PID pid, const int sign, o2::track::TrackParCovF& track)
  {
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

  bool fillKFinfo(const KFParticle& kfParticle, StrangeTrack& strangeTrack) {
    float M, SigmaM;
    kfParticle.GetMass(M, SigmaM);
    strangeTrack.mMasses[0] = M;
    for (int i = 0; i < 8; i++) {
      strangeTrack.mDecayParams[i] = kfParticle.GetParameter(i);
    }
    for (int i = 0; i < 36; i++) {
      strangeTrack.mDecayCov[i] = kfParticle.GetCovariance(i); // lower triangular form of covariance matrix
    }
    strangeTrack.mGeoChi2 = kfParticle.GetChi2();
    return true;
  }

  std::vector<ITSCluster> getTrackClusters(const TrackITS& itsTrack)
  {
    std::vector<ITSCluster> outVec;
    outVec.reserve(7);
    auto firstClus = itsTrack.getFirstClusterEntry();
    auto ncl = itsTrack.getNumberOfClusters();
    for (int icl = 0; icl < ncl; icl++) {
      outVec.push_back(mInputITSclusters[mInputITSidxs[firstClus + icl]]);
    }
    return outVec;
  };

  std::vector<int> getTrackClusterSizes(const TrackITS& itsTrack)
  {
    std::vector<int> outVec;
    outVec.reserve(7);
    auto firstClus = itsTrack.getFirstClusterEntry();
    auto ncl = itsTrack.getNumberOfClusters();
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

  o2::MCCompLabel getStrangeTrackLabel(const TrackITS& itsTrack, const StrangeTrack& strangeTrack, const ClusAttachments& structClus) // ITS label with fake flag recomputed
  {
    bool isFake = false;
    auto itsTrkLab = mITSTrkLabels[strangeTrack.mITSRef];
    for (unsigned int iLay = 0; iLay < 7; iLay++) {
      if (itsTrack.hasHitOnLayer(iLay) && itsTrack.isFakeOnLayer(iLay) && structClus.arr[iLay] == 0) {
        isFake = true;
        break;
      }
    }
    itsTrkLab.setFakeFlag(isFake);
    return itsTrkLab;
  }

 protected:
  bool mMCTruthON = false;                             /// flag availability of MC truth
  int mNThreads = 1;                                   /// number of threads (externally driven)
  gsl::span<const TrackITS> mInputITStracks;           // input ITS tracks
  std::vector<VBracket> mITSvtxBrackets;               // time brackets for ITS tracks
  std::vector<int> mTracksIdxTable;                    // index table for ITS tracks
  std::vector<int> mInputClusterSizes;                 // input cluster sizes
  std::vector<ITSCluster> mInputITSclusters;           // input ITS clusters
  gsl::span<const int> mInputITSidxs;                  // input ITS track-cluster indexes
  gsl::span<const V0> mInputV0tracks;                  // input V0 of decay daughters
  gsl::span<const V0Index> mInputV0Indices;            // input V0 indices of decay daughters
  gsl::span<const Cascade> mInputCascadeTracks;        // input cascade of decay daughters
  gsl::span<const CascadeIndex> mInputCascadeIndices;  // input cascade indices of decay daughters
  gsl::span<const Decay3Body> mInput3BodyTracks;       // input decay3body of decay daughters
  gsl::span<const Decay3BodyIndex> mInput3BodyIndices; // input decay3body indices of decay daughters
  const MCLabContCl* mITSClsLabels = nullptr;          /// input ITS Cluster MC labels
  MCLabSpan mITSTrkLabels;                             /// input ITS Track MC labels

  std::vector<o2::its::TrackITS> mSortedITStracks; // sorted ITS tracks
  std::vector<int> mSortedITSindexes;              // indexes of sorted ITS tracks
  IndexTableUtils mUtils;                          // structure for computing eta/phi matching selections

  std::vector<std::vector<StrangeTrack>> mStrangeTrackVec;       // structure containing updated mother and daughter tracks (per thread)
  std::vector<std::vector<ClusAttachments>> mClusAttachments;    // # of attached tracks, -1 not attached, 0 for the mother, > 0 for the daughters (per thread)
  std::vector<std::vector<o2::MCCompLabel>> mStrangeTrackLabels; // vector of MC labels for mother track (per thread)

  const StrangenessTrackingParamConfig* mStrParams = nullptr;
  float mBz = -5; // Magnetic field
  const o2::itsmft::TopologyDictionary* mDict = nullptr;

  std::vector<DCAFitter2> mFitterV0;    // optional DCA Fitter for recreating V0 with hypertriton mass hypothesis (per thread)
  std::vector<DCAFitter3> mFitter3Body; // optional DCA Fitter for final 3 Body refit (per thread)
  std::vector<DCAFitter4> mFitter4Body; // optional DCA Fitter for final 4 Body refit (per thread)

  o2::base::PropagatorImpl<float>::MatCorrType mCorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrNONE; // use mat correction

  std::vector<std::vector<o2::track::TrackParCovF>> mDaughterTracks; // vector of daughter tracks (per thread)

  ClusAttachments mStructClus;                          // # of attached tracks, 1 for mother, 2 for daughter
  KFParticle kfpMother;                                  // mother KFParticle
  o2::track::PID pidV0;                                  // PID hypothesis for the V0 fitting
  o2::track::PID pidV0comp;                              // commpeting PID hypothesis for the V0 fitting
  o2::track::PID pidCasc;                                // PID hypothesis for the cascade fitting
  o2::track::PID pidCascComp;                            // competing PID hypothesis for the cascade fitting

  ClassDefNV(StrangenessTracker, 1);
};

} // namespace strangeness_tracking
} // namespace o2

#endif //  _ALICEO2_STRANGENESS_TRACKER_
