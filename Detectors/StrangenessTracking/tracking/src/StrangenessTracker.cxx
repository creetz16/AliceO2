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
/// \file StrangenessTracker.cxx
/// \brief

#include <numeric>
#include "StrangenessTracking/StrangenessTracker.h"
#include "ITStracking/IOUtils.h"

namespace o2
{
namespace strangeness_tracking
{

template <typename T>
KFParticle createKFParticleFromTrackParCov(const o2::track::TrackParametrizationWithError<T>& trackparCov, int charge, T mass) {
  T xyzpxpypz[6];
  o2::gpu::gpustd::array<T, 3> xyz, mom;
  trackparCov.getPxPyPzGlo(mom);
  trackparCov.getXYZGlo(xyz);
  for (int i{0}; i < 3; ++i) {
    xyzpxpypz[i] = xyz[i];
    xyzpxpypz[i + 3] = mom[i];
  }

  o2::gpu::gpustd::array<T, 21> cv;
  //trackparCov.getCovXYZPxPyPzGlo(cv);
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
  //kfPart.Create(xyzpxpypz, cv.data(), charge, mass);
  return kfPart;
}

bool StrangenessTracker::loadData(const o2::globaltracking::RecoContainer& recoData)
{
  clear();

  mInputV0tracks = recoData.getV0s();
  mInputCascadeTracks = recoData.getCascades();
  mInputITStracks = recoData.getITSTracks();
  mInputITSidxs = recoData.getITSTracksClusterRefs();

  auto compClus = recoData.getITSClusters();
  auto clusPatt = recoData.getITSClustersPatterns();
  auto pattIt = clusPatt.begin();
  mInputITSclusters.reserve(compClus.size());
  mInputClusterSizes.resize(compClus.size());
  o2::its::ioutils::convertCompactClusters(compClus, pattIt, mInputITSclusters, mDict);
  auto pattIt2 = clusPatt.begin();
  getClusterSizes(mInputClusterSizes, compClus, pattIt2, mDict);

  mITSvtxBrackets.resize(mInputITStracks.size());
  for (int i = 0; i < mInputITStracks.size(); i++) {
    mITSvtxBrackets[i] = {-1, -1};
  }

  // build time bracket for each ITS track
  auto trackIndex = recoData.getPrimaryVertexMatchedTracks(); // Global ID's for associated tracks
  auto vtxRefs = recoData.getPrimaryVertexMatchedTrackRefs(); // references from vertex to these track IDs

  if (mStrParams->mVertexMatching) {
    int nv = vtxRefs.size();
    for (int iv = 0; iv < nv; iv++) {
      const auto& vtref = vtxRefs[iv];
      int it = vtref.getFirstEntry(), itLim = it + vtref.getEntries();
      for (; it < itLim; it++) {
        auto tvid = trackIndex[it];
        if (!recoData.isTrackSourceLoaded(tvid.getSource()) || tvid.getSource() != GIndex::ITS) {
          continue;
        }
        if (mITSvtxBrackets[tvid.getIndex()].getMin() == -1) {
          mITSvtxBrackets[tvid.getIndex()].setMin(iv);
          mITSvtxBrackets[tvid.getIndex()].setMax(iv);
        } else {
          mITSvtxBrackets[tvid.getIndex()].setMax(iv);
        }
      }
    }
  }

  if (mMCTruthON) {
    mITSClsLabels = recoData.mcITSClusters.get();
    mITSTrkLabels = recoData.getITSTracksMCLabels();
  }

  LOG(debug) << "V0 tracks size: " << mInputV0tracks.size();
  LOG(debug) << "Cascade tracks size: " << mInputCascadeTracks.size();
  LOG(debug) << "ITS tracks size: " << mInputITStracks.size();
  LOG(debug) << "ITS idxs size: " << mInputITSidxs.size();
  LOG(debug) << "ITS clusters size: " << mInputITSclusters.size();
  LOG(debug) << "VtxRefs size: " << vtxRefs.size();

  return true;
}

void StrangenessTracker::prepareITStracks() // sort tracks by eta and phi and select only tracks with vertex matching
{

  for (int iTrack{0}; iTrack < mInputITStracks.size(); iTrack++) {
    if (mStrParams->mVertexMatching && mITSvtxBrackets[iTrack].getMin() == -1) {
      continue;
    }
    mSortedITStracks.push_back(mInputITStracks[iTrack]);
    mSortedITSindexes.push_back(iTrack);
  }

  mTracksIdxTable.resize(mUtils.mPhiBins * mUtils.mEtaBins + 1);
  std::sort(mSortedITStracks.begin(), mSortedITStracks.end(), [&](o2::its::TrackITS& a, o2::its::TrackITS& b) { return mUtils.getBinIndex(a.getEta(), a.getPhi()) < mUtils.getBinIndex(b.getEta(), b.getPhi()); });
  std::sort(mSortedITSindexes.begin(), mSortedITSindexes.end(), [&](int i, int j) { return mUtils.getBinIndex(mInputITStracks[i].getEta(), mInputITStracks[i].getPhi()) < mUtils.getBinIndex(mInputITStracks[j].getEta(), mInputITStracks[j].getPhi()); });

  for (auto& track : mSortedITStracks) {
    mTracksIdxTable[mUtils.getBinIndex(track.getEta(), track.getPhi())]++;
  }
  std::exclusive_scan(mTracksIdxTable.begin(), mTracksIdxTable.begin() + mUtils.mPhiBins * mUtils.mEtaBins, mTracksIdxTable.begin(), 0);
  mTracksIdxTable[mUtils.mPhiBins * mUtils.mEtaBins] = mSortedITStracks.size();
}

void StrangenessTracker::process()
{
  #ifdef HomogeneousField
    KFParticle::SetField(mBz);
  #endif

  // Loop over V0s
  mDaughterTracks.resize(2); // resize to 2 prongs: first positive second negative

  for (int iV0{0}; iV0 < mInputV0tracks.size(); iV0++) {
    LOG(debug) << "Analysing V0: " << iV0 + 1 << "/" << mInputV0tracks.size();
    auto& v0 = mInputV0tracks[iV0];
    mV0dauIDs[kV0DauPos] = v0.getProngID(kV0DauPos), mV0dauIDs[kV0DauNeg] = v0.getProngID(kV0DauNeg);
    auto posTrack = v0.getProng(kV0DauPos);
    auto negTrack = v0.getProng(kV0DauNeg);
    auto alphaV0 = calcV0alpha(v0);
    alphaV0 > 0 ? posTrack.setAbsCharge(2) : negTrack.setAbsCharge(2);
    V0 correctedV0; // recompute V0 for Hypertriton

    if (!recreateV0(posTrack, negTrack, correctedV0)) {
      continue;
    }

    mStrangeTrack.mPartType = dataformats::kStrkV0;

    auto v0R2 = v0.calcR2();
    auto iBinsV0 = mUtils.getBinRect(correctedV0.getEta(), correctedV0.getPhi(), mStrParams->mEtaBinSize, mStrParams->mPhiBinSize);
    for (int& iBinV0 : iBinsV0) {
      for (int iTrack{mTracksIdxTable[iBinV0]}; iTrack < TMath::Min(mTracksIdxTable[iBinV0 + 1], int(mSortedITStracks.size())); iTrack++) {
        mStrangeTrack.mMother = (o2::track::TrackParCovF)correctedV0;
        mDaughterTracks[kV0DauPos] = correctedV0.getProng(kV0DauPos);
        mDaughterTracks[kV0DauNeg] = correctedV0.getProng(kV0DauNeg);
        mITStrack = mSortedITStracks[iTrack];
        auto& ITSindexRef = mSortedITSindexes[iTrack];

        if (mStrParams->mVertexMatching && (mITSvtxBrackets[ITSindexRef].getMin() > v0.getVertexID() ||
                                            mITSvtxBrackets[ITSindexRef].getMax() < v0.getVertexID())) {
          continue;
        }

        if (matchDecayToITStrack(sqrt(v0R2))) {

          auto propInstance = o2::base::Propagator::Instance();
          o2::track::TrackParCov decayVtxTrackClone = mStrangeTrack.mMother; // clone track and propagate to decay vertex
          if (!propInstance->propagateToX(decayVtxTrackClone, mStrangeTrack.mDecayVtx[0], getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mCorrType)) {
            LOG(debug) << "Mother propagation to decay vertex failed";
            continue;
          }

          decayVtxTrackClone.getPxPyPzGlo(mStrangeTrack.mDecayMom);
          std::array<float, 3> momPos, momNeg;
          mFitter3Body.getTrack(kV0DauPos).getPxPyPzGlo(momPos);
          mFitter3Body.getTrack(kV0DauNeg).getPxPyPzGlo(momNeg);
          if (alphaV0 > 0) {
            mStrangeTrack.mMasses[0] = calcMotherMass(momPos, momNeg, PID::Helium3, PID::Pion); // Hypertriton invariant mass at decay vertex
            mStrangeTrack.mMasses[1] = calcMotherMass(momPos, momNeg, PID::Alpha, PID::Pion);   // Hyperhydrogen4Lam invariant mass at decay vertex
          } else {
            mStrangeTrack.mMasses[0] = calcMotherMass(momPos, momNeg, PID::Helium3, PID::Pion); // Anti-Hypertriton invariant mass at decay vertex
            mStrangeTrack.mMasses[1] = calcMotherMass(momPos, momNeg, PID::Alpha, PID::Pion);   // Anti-Hyperhydrogen4Lam invariant mass at decay vertex
          }

          LOG(debug) << "ITS Track matched with a V0 decay topology ....";
          LOG(debug) << "Number of ITS track clusters attached: " << mITStrack.getNumberOfClusters();
          mStrangeTrack.mDecayRef = iV0;
          mStrangeTrack.mITSRef = mSortedITSindexes[iTrack];
          mStrangeTrackVec.push_back(mStrangeTrack);
          mClusAttachments.push_back(mStructClus);
          if (mMCTruthON) {
            auto lab = getStrangeTrackLabel();
            mStrangeTrackLabels.push_back(lab);
          }
        }
      }
    }
  }

  // Loop over Cascades
  mDaughterTracks.resize(3); // resize to 3 prongs: first V0 pos, second V0 neg, third bachelor

  for (int iCasc{0}; iCasc < mInputCascadeTracks.size(); iCasc++) {
    LOG(debug) << "Analysing Cascade: " << iCasc + 1 << "/" << mInputCascadeTracks.size();

    auto& casc = mInputCascadeTracks[iCasc];
    auto& cascV0 = mInputV0tracks[casc.getV0ID()];
    mV0dauIDs[kV0DauPos] = cascV0.getProngID(kV0DauPos);
    mV0dauIDs[kV0DauNeg] = cascV0.getProngID(kV0DauNeg);

    mStrangeTrack.mPartType = dataformats::kStrkCascade;

    // first: bachelor, second: V0 pos, third: V0 neg
    auto cascR2 = casc.calcR2();
    auto iBinsCasc = mUtils.getBinRect(casc.getEta(), casc.getPhi(), mStrParams->mEtaBinSize, mStrParams->mPhiBinSize);
    for (int& iBinCasc : iBinsCasc) {
      for (int iTrack{mTracksIdxTable[iBinCasc]}; iTrack < TMath::Min(mTracksIdxTable[iBinCasc + 1], int(mSortedITStracks.size())); iTrack++) {
        mStrangeTrack.mMother = (o2::track::TrackParCovF)casc;
        mDaughterTracks[kV0DauPos] = cascV0.getProng(kV0DauPos);
        mDaughterTracks[kV0DauNeg] = cascV0.getProng(kV0DauNeg);
        mDaughterTracks[kBach] = casc.getBachelorTrack();
        mITStrack = mSortedITStracks[iTrack];
        auto& ITSindexRef = mSortedITSindexes[iTrack];
        LOG(debug) << "----------------------";
        LOG(debug) << "CascV0: " << casc.getV0ID() << ", Bach ID: " << casc.getBachelorID() << ", ITS track ref: " << mSortedITSindexes[iTrack];

        if (mStrParams->mVertexMatching && (mITSvtxBrackets[ITSindexRef].getMin() > casc.getVertexID() ||
                                            mITSvtxBrackets[ITSindexRef].getMax() < casc.getVertexID())) {
          LOG(debug) << "Vertex ID mismatch: " << mITSvtxBrackets[ITSindexRef].getMin() << " < " << casc.getVertexID() << " < " << mITSvtxBrackets[ITSindexRef].getMax();
          continue;
        }

        if (matchDecayToITStrack(sqrt(cascR2))) {

          auto propInstance = o2::base::Propagator::Instance();
          o2::track::TrackParCov decayVtxTrackClone = mStrangeTrack.mMother; // clone track and propagate to decay vertex
          if (!propInstance->propagateToX(decayVtxTrackClone, mStrangeTrack.mDecayVtx[0], getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mCorrType)) {
            LOG(debug) << "Mother propagation to decay vertex failed";
            continue;
          }
          decayVtxTrackClone.getPxPyPzGlo(mStrangeTrack.mDecayMom);
          std::array<float, 3> momV0, mombach;
          mFitter3Body.getTrack(0).getPxPyPzGlo(momV0);                                      // V0 momentum at decay vertex
          mFitter3Body.getTrack(1).getPxPyPzGlo(mombach);                                    // bachelor momentum at decay vertex
          mStrangeTrack.mMasses[0] = calcMotherMass(momV0, mombach, PID::Lambda, PID::Pion); // Xi invariant mass at decay vertex
          mStrangeTrack.mMasses[1] = calcMotherMass(momV0, mombach, PID::Lambda, PID::Kaon); // Omega invariant mass at decay vertex

          LOG(debug) << "ITS Track matched with a Cascade decay topology ....";
          LOG(debug) << "Number of ITS track clusters attached: " << mITStrack.getNumberOfClusters();

          mStrangeTrack.mDecayRef = iCasc;
          mStrangeTrack.mITSRef = mSortedITSindexes[iTrack];
          mStrangeTrackVec.push_back(mStrangeTrack);
          mClusAttachments.push_back(mStructClus);
          if (mMCTruthON) {
            auto lab = getStrangeTrackLabel();
            mStrangeTrackLabels.push_back(lab);
          }
        }
      }
    }
  }
}

bool StrangenessTracker::matchDecayToITStrack(float decayR)
{
  auto geom = o2::its::GeometryTGeo::Instance();
  auto trackClusters = getTrackClusters();
  auto trackClusSizes = getTrackClusterSizes();
  auto& lastClus = trackClusters[0];
  mStrangeTrack.mMatchChi2 = getMatchingChi2(mStrangeTrack.mMother, mITStrack);

  auto radTol = decayR < 4 ? mStrParams->mRadiusTolIB : mStrParams->mRadiusTolOB;
  auto nMinClusMother = trackClusters.size() < 4 ? 2 : mStrParams->mMinMotherClus;

  std::vector<ITSCluster> motherClusters;
  std::vector<int> motherClusSizes;
  std::array<unsigned int, 7> nAttachments;
  nAttachments.fill(-1); // fill arr with -1

  int nUpdates = 0;
  bool isMotherUpdated = false;

  for (int iClus{0}; iClus < trackClusters.size(); iClus++) {
    auto& clus = trackClusters[iClus];
    auto& compClus = trackClusSizes[iClus];
    int nUpdOld = nUpdates;
    double clusRad = sqrt(clus.getX() * clus.getX() - clus.getY() * clus.getY());
    auto diffR = decayR - clusRad;
    auto relDiffR = diffR / decayR;
    // Look for the Mother if the Decay radius allows for it, within a tolerance
    LOG(debug) << "decayR: " << decayR << ", diffR: " << diffR << ", clus rad: " << clusRad << ", radTol: " << radTol;
    if (relDiffR > -radTol) {
      LOG(debug) << "Try to attach cluster to Mother, layer: " << geom->getLayer(clus.getSensorID());
      if (updateTrack(clus, mStrangeTrack.mMother)) {
        motherClusters.push_back(clus);
        motherClusSizes.push_back(compClus);
        nAttachments[geom->getLayer(clus.getSensorID())] = 0;
        isMotherUpdated = true;
        nUpdates++;
        LOG(debug) << "Cluster attached to Mother";
        continue; // if the cluster is attached to the mother, skip the rest of the loop
      }
    }

    // if Mother is not found, check for V0 daughters compatibility
    if (relDiffR < radTol && !isMotherUpdated) {
      bool isDauUpdated = false;
      LOG(debug) << "Try to attach cluster to Daughters, layer: " << geom->getLayer(clus.getSensorID());
      for (int iDau{0}; iDau < mDaughterTracks.size(); iDau++) {
        auto& dauTrack = mDaughterTracks[iDau];
        if (updateTrack(clus, dauTrack)) {
          nAttachments[geom->getLayer(clus.getSensorID())] = iDau + 1;
          isDauUpdated = true;
          break;
        }
      }
      if (!isDauUpdated) {
        break; // no daughter track updated, stop the loop
      }
      nUpdates++;
    }
    if (nUpdates == nUpdOld) {
      break; // no track updated, stop the loop
    }
  }

  if (nUpdates < trackClusters.size() || motherClusters.size() < nMinClusMother) {
    return false;
  }

  o2::track::TrackParCov motherTrackClone = mStrangeTrack.mMother; // clone and reset covariance for final topology refit
  motherTrackClone.resetCovariance();

  LOG(debug) << "Clusters attached, starting inward-outward refit";

  std::reverse(motherClusters.begin(), motherClusters.end());

  for (auto& clus : motherClusters) {
    if (!updateTrack(clus, motherTrackClone)) {
      break;
    }
  }

  // compute mother average cluster size
  mStrangeTrack.mITSClusSize = float(std::accumulate(motherClusSizes.begin(), motherClusSizes.end(), 0)) / motherClusSizes.size();

  LOG(debug) << "Inward-outward refit finished, starting final topology refit";
  // final Topology refit

  int cand = 0; // best V0 candidate
  int nCand;
  int nDauPos;

  // ========== refit cascade ============
  if (mStrangeTrack.mPartType == dataformats::kStrkCascade) {
    // Do cascade with V0 from KF --> auto KFParticle cascade with create function
    // print here the vertex I get from the DCA fitter and the one fomr the KF!!!
    // How it is done in DCAFitter: GetPCA and then get XYZ from that (or similar)
    // create cascade from Lam and bachelor (with mass from pion) --> check if it is far away from Xi mass --> if yes, SetMass of daughter track to kaon and do mother again --> if it is omega, store is as omega
    // in principle I can also constrain Lambda to be a Lambda with mass constraint (= one of the next steps)

    /// ======= DCA fitter reconstruction =======
    // refit cascade
    V0 cascV0Upd;
    if (!recreateV0(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], cascV0Upd)) {
      LOG(debug) << "Cascade V0 refit failed";
      return false;
    }
    cascV0Upd.setAbsCharge(0);
    try {
      nCand = mFitter3Body.process(cascV0Upd, mDaughterTracks[kBach], motherTrackClone);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Fitter3Body failed: " << e.what();
      return false;
    }
    if (!nCand || !mFitter3Body.propagateTracksToVertex()) {
      LOG(debug) << "Fitter3Body failed: propagation to vertex failed";
      return false;
    }
    int nCasc{0};
    try {
      nCasc = mFitterV0.process(cascV0Upd, mDaughterTracks[kBach]);  // refit V0
    } catch (std::runtime_error& e) {
      LOG(debug) << "Cascade refit failed: " << e.what();
    }
    if (nCasc) {
      mResettedMotherTrack = mFitterV0.createParentTrackParCov();
    }

    /// ======== KF reconstruction =========
    float massProton = o2::constants::physics::MassProton;
    float massPion = o2::constants::physics::MassPionCharged;
    float massKaon = o2::constants::physics::MassKaonCharged;
    float massLambda = o2::constants::physics::MassLambda;
    float massXi = o2::constants::physics::MassXiMinus;
    // set daughter masses
    float massPosDaughter, massNegDaughter;
    if (mDaughterTracks[kBach].getCharge() < 0) { // if bachelor charge negative, cascade is a Xi- --> Lam + pi- --> (p+ + pi-) + pi-   OR   Omega- --> Lam + K- --> (p+ + pi-) + K-
      massPosDaughter = massProton;
      massNegDaughter = massPion;
    }  else { // if bachelor charge positive, cascade is a Xi+ --> Lam + pi+ --> (p- + pi+) + pi+  OR   Omega+ --> Lam + K+ --> (p- + pi+) + K+
      massPosDaughter = massPion;
      massNegDaughter = massProton;
    }
    /// V0 reconstruction
    KFParticle cascV0KF;
    int nV0Daughters = 2;
    // create KFParticle objects from trackParCovs
    KFParticle kfpDauPos = createKFParticleFromTrackParCov(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauPos].getCharge(), massPosDaughter); // positive prong
    KFParticle kfpDauNeg = createKFParticleFromTrackParCov(mDaughterTracks[kV0DauNeg], mDaughterTracks[kV0DauNeg].getCharge(), massNegDaughter); // negative prong
    const KFParticle* V0Daughters[2] = {&kfpDauPos, &kfpDauNeg};
    mStrangeTrack.kfpDauPosMass = kfpDauPos.GetMass();
    mStrangeTrack.kfpDauNegMass = kfpDauNeg.GetMass();
    if (kfpDauPos.GetMass()==0 || kfpDauNeg.GetMass()==0) {
      LOG(debug) << "Daughter mass 0. Excluding candidate.";
      return false;
    }
    // construct mother
    cascV0KF.SetConstructMethod(2);
    //cascV0KF.Construct(V0Daughters, nV0Daughters);
    try {
      cascV0KF.Construct(V0Daughters, nV0Daughters);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to construct cascade V0 from daughter tracks: " << e.what();
      return false;
    }

    /// cascade reconstruction: Xi
    KFParticle cascKFXi;
    KFParticle cascKFXi_wMassConstLam;
    int nCascDaughters = 2;
    // set mass constraint to V0
    KFParticle cascV0KF_wMassConst = cascV0KF;
    // cascV0KF_wMassConst.SetNonlinearMassConstraint(massLambda);
    try {
      cascV0KF_wMassConst.SetNonlinearMassConstraint(massLambda);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to set non-linear mass constraint on V0 from cascade: " << e.what();
      return false;
    }
    // create KFParticle objects from trackParCovs
    KFParticle kfpBachelorXi = createKFParticleFromTrackParCov(mDaughterTracks[kBach], mDaughterTracks[kBach].getCharge(), massPion); // bachelor
    if (kfpBachelorXi.GetMass()==0) {
      LOG(debug) << "Bachelor mass 0. Excluding candidate."; 
      return false;
    }
    const KFParticle* CascDaugthersXi[2] = {&kfpBachelorXi, &cascV0KF};
    const KFParticle* CascDaugthersXi_wMassConstLambda[2] = {&kfpBachelorXi, &cascV0KF_wMassConst};
    mStrangeTrack.kfpBachelorMass = kfpBachelorXi.GetMass();
    mStrangeTrack.kfpCascV0Mass = cascV0KF.GetMass();
    mStrangeTrack.kfpCascV0MassConst = cascV0KF_wMassConst.GetMass();
    // construct mother
    cascKFXi.SetConstructMethod(2);
    cascKFXi.Construct(CascDaugthersXi, nCascDaughters);
    cascKFXi_wMassConstLam.SetConstructMethod(2);
    cascKFXi_wMassConstLam.Construct(CascDaugthersXi_wMassConstLambda, nCascDaughters);
    // transport mother to decay Vtx
    KFParticle cascKFXidecayVtx = cascKFXi;
    KFParticle cascKFXidecayVtx_wMassConstLam = cascKFXi_wMassConstLam;
    cascKFXidecayVtx.TransportToDecayVertex();
    cascKFXidecayVtx_wMassConstLam.TransportToDecayVertex();
    mResettedMotherTrackKF = cascKFXidecayVtx;
    mResettedMotherTrackKF_wMassConstLam = cascKFXidecayVtx_wMassConstLam;
    // set mass of strange track
    mStrangeTrack.mMassesKF[0] = mResettedMotherTrackKF.GetMass(); // cascades with Xi mass hypothesis

    /// cascade reconstruction: Omega
    KFParticle cascKFOm;
    KFParticle cascKFOm_wMassConstLam;
    // create KFParticle objects from trackParCovs
    KFParticle kfpBachelorOm = createKFParticleFromTrackParCov(mDaughterTracks[kBach], mDaughterTracks[kBach].getCharge(), massKaon); // bachelor
    if (kfpBachelorOm.GetMass()==0) {
      LOG(debug) << "Bachelor mass 0. Excluding candidate."; 
      return false;
    }
    const KFParticle* CascDaugthersOm[2] = {&kfpBachelorOm, &cascV0KF};
    const KFParticle* CascDaugthersOm_wMassConstLambda[2] = {&kfpBachelorOm, &cascV0KF_wMassConst};
    // construct mother
    cascKFOm.SetConstructMethod(2);
    cascKFOm.Construct(CascDaugthersOm, nCascDaughters);
    cascKFOm_wMassConstLam.SetConstructMethod(2);
    cascKFOm_wMassConstLam.Construct(CascDaugthersOm_wMassConstLambda, nCascDaughters);
    // transport mother to decay Vtx
    KFParticle cascKFOmdecayVtx = cascKFOm;
    KFParticle cascKFOmdecayVtx_wMassConstLam = cascKFOm_wMassConstLam;
    cascKFOmdecayVtx.TransportToDecayVertex();
    cascKFOmdecayVtx_wMassConstLam.TransportToDecayVertex();
    mResettedMotherTrackKF = cascKFOmdecayVtx;
    mResettedMotherTrackKF_wMassConstLam = cascKFOmdecayVtx_wMassConstLam;
    // set mass of strange track
    mStrangeTrack.mMassesKF[1] = cascKFOmdecayVtx.GetMass(); // cascades with mass of Omega
    
  }


  // ========== refit V0 ==========
  else if (mStrangeTrack.mPartType == dataformats::kStrkV0) {
    /// ======= DCA fitter reconstruction =======
    try {
      nCand = mFitter3Body.process(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], motherTrackClone);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Fitter3Body failed: " << e.what();
      return false;
    }
    if (!nCand || !mFitter3Body.propagateTracksToVertex()) {
      LOG(debug) << "Fitter3Body failed: propagation to vertex failed";
      return false;
    }
    int nCasc{0};
    try {
      nCasc = mFitterV0.process(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg]);  // refit V0
    } catch (std::runtime_error& e) {
      LOG(debug) << "V0 refit failed " << e.what();
    }
    if (nCasc) {
      mResettedMotherTrack = mFitterV0.createParentTrackParCov();
    }

    /// ======== KF reconstruction =========
    float massHelium3 = o2::constants::physics::MassHelium3;
    float massAlpha = o2::constants::physics::MassAlpha;
    float massPion = o2::constants::physics::MassPionCharged;
    float massHyperTriton = o2::constants::physics::MassHyperTriton;
    // set daughter masses
    float massPosDaughter, massNegDaughter;
    if (mDaughterTracks[kV0DauPos].getAbsCharge() == 2) { // if charge of positive daughter is two, it is 3He (or 4He)
      massPosDaughter = massHelium3;
      massNegDaughter = massPion;
    } else {
      massPosDaughter = massPion;
      massNegDaughter = massHelium3;
    }
    // V0 reconstruction: hypertriton
    KFParticle V0KFtr;
    int nV0Daughters = 2;
    // create KFParticle objects from trackParCovs
    KFParticle kfpDauPos = createKFParticleFromTrackParCov(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauPos].getCharge(), massPosDaughter); // positive prong
    KFParticle kfpDauNeg = createKFParticleFromTrackParCov(mDaughterTracks[kV0DauNeg], mDaughterTracks[kV0DauNeg].getCharge(), massNegDaughter); // negative prong
    const KFParticle* V0Daughters_tr[2] = {&kfpDauPos, &kfpDauNeg};
    mStrangeTrack.kfpDauPosMass = kfpDauPos.GetMass();
    mStrangeTrack.kfpDauNegMass = kfpDauNeg.GetMass();
    if (kfpDauPos.GetMass()==0 || kfpDauNeg.GetMass()==0) {
      LOG(debug) << "Daughter mass 0. Excluding candidate."; 
      return false;
    }
    // construct mother
    V0KFtr.SetConstructMethod(2);
    V0KFtr.Construct(V0Daughters_tr, nV0Daughters);
    KFParticle V0KFtrdecayVtx = V0KFtr;
    V0KFtrdecayVtx.TransportToDecayVertex();
    mResettedMotherTrackKF = V0KFtrdecayVtx;
    // set mass of strange track
    mStrangeTrack.mMassesKF[0] = mResettedMotherTrackKF.GetMass(); // V0s with hypertriton invariant mass hypothesis

    // V0 reconstruction: hyperhydrogen
    KFParticle V0KFhydro;
    // set daughter masses
    if (mDaughterTracks[kV0DauPos].getAbsCharge() == 2) { // if charge of positive daughter is 2, it is 4He
      massPosDaughter = massAlpha;
      massNegDaughter = massPion;
    }  else {
      massPosDaughter = massPion;
      massNegDaughter = massAlpha;
    }
    // create KFParticle objects from trackParCovs
    KFParticle kfpDauPos_hydro = createKFParticleFromTrackParCov(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauPos].getCharge(), massPosDaughter); // positive prong
    KFParticle kfpDauNeg_hydro = createKFParticleFromTrackParCov(mDaughterTracks[kV0DauNeg], mDaughterTracks[kV0DauNeg].getCharge(), massNegDaughter); // negative prong
    const KFParticle* V0Daughters_hydro[2] = {&kfpDauPos_hydro, &kfpDauNeg_hydro};
    if (kfpDauPos_hydro.GetMass()==0 || kfpDauNeg_hydro.GetMass()==0) {
      LOG(debug) << "Daughter mass 0. Excluding candidate."; 
      return false;
    }
    // construct mother
    V0KFhydro.SetConstructMethod(2);
    V0KFhydro.Construct(V0Daughters_hydro, nV0Daughters);
    KFParticle V0KFhydrodecayVtx = V0KFhydro;
    V0KFhydrodecayVtx.TransportToDecayVertex();
    // set mass of strange track
    mStrangeTrack.mMassesKF[1] = V0KFhydrodecayVtx.GetMass(); // V0s with mass of hyperhydrogen4

    mStrangeTrack.kfpBachelorMass = -999;
    mStrangeTrack.kfpCascV0Mass = -999;
    mStrangeTrack.kfpCascV0MassConst = -999;
  }

  // get vertex position and chi2 of refitted track
  mStrangeTrack.mDecayVtx = mFitter3Body.getPCACandidatePos();
  mStrangeTrack.decayVtxX = mStrangeTrack.mDecayVtx[0];
  mStrangeTrack.decayVtxY = mStrangeTrack.mDecayVtx[1];
  mStrangeTrack.decayVtxZ = mStrangeTrack.mDecayVtx[2];
  mStrangeTrack.mTopoChi2 = mFitter3Body.getChi2AtPCACandidate();
  // // get vertex properties of refitted KF track
  mStrangeTrack.decayVtxXKF = mResettedMotherTrackKF.GetX();
  mStrangeTrack.decayVtxYKF = mResettedMotherTrackKF.GetY();
  mStrangeTrack.decayVtxZKF = mResettedMotherTrackKF.GetZ();
  mStrangeTrack.mGeoChi2KF = mResettedMotherTrackKF.GetChi2();
  mStrangeTrack.mPtKF = mResettedMotherTrackKF.GetPt();

  mStructClus.arr = nAttachments;

  return true;
}

bool StrangenessTracker::updateTrack(const ITSCluster& clus, o2::track::TrackParCov& track)
{
  auto geom = o2::its::GeometryTGeo::Instance();
  auto propInstance = o2::base::Propagator::Instance();
  float alpha = geom->getSensorRefAlpha(clus.getSensorID()), x = clus.getX();
  int layer{geom->getLayer(clus.getSensorID())};

  if (!track.rotate(alpha)) {
    return false;
  }

  if (!propInstance->propagateToX(track, x, getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mCorrType)) {
    return false;
  }

  if (mCorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE) {
    float thick = layer < 3 ? 0.005 : 0.01;
    constexpr float radl = 9.36f; // Radiation length of Si [cm]
    constexpr float rho = 2.33f;  // Density of Si [g/cm^3]
    if (!track.correctForMaterial(thick, thick * rho * radl)) {
      return false;
    }
  }
  auto chi2 = std::abs(track.getPredictedChi2(clus)); // abs to be understood
  LOG(debug) << "Chi2: " << chi2;
  if (chi2 > mStrParams->mMaxChi2 || chi2 < 0) {
    return false;
  }

  if (!track.update(clus)) {
    return false;
  }

  return true;
}

} // namespace strangeness_tracking
} // namespace o2