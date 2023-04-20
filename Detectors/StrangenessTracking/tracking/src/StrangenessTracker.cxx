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
  trackparCov.getCovXYZPxPyPzGlo(cv);

  KFParticle kfPart;
  kfPart.Create(xyzpxpypz, cv.data(), charge, mass);
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
    auto& DecIndexRef = iV0;
    auto& v0 = mInputV0tracks[iV0];
    mV0dauIDs[0] = v0.getProngID(0), mV0dauIDs[1] = v0.getProngID(1);
    auto posTrack = v0.getProng(0);
    auto negTrack = v0.getProng(1);
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
        mDaughterTracks[0] = correctedV0.getProng(0);
        mDaughterTracks[1] = correctedV0.getProng(1);
        mITStrack = mSortedITStracks[iTrack];
        auto& ITSindexRef = mSortedITSindexes[iTrack];
        LOG(debug) << "V0 pos: " << v0.getProngID(0) << " V0 neg: " << v0.getProngID(1) << ", ITS track ref: " << mSortedITSindexes[iTrack];
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
          auto p2mom = decayVtxTrackClone.getP2();
          auto p2pos = mFitter3Body.getTrack(0).getP2(); // positive V0 daughter
          auto p2neg = mFitter3Body.getTrack(1).getP2(); // negative V0 daughter
          if (alphaV0 > 0) {
            mStrangeTrack.mMasses[0] = calcMotherMass(p2mom, p2pos, p2neg, PID::Helium3, PID::Pion); // Hypertriton invariant mass at decay vertex
            mStrangeTrack.mMasses[1] = calcMotherMass(p2mom, p2pos, p2neg, PID::Alpha, PID::Pion);   // Hyperhydrogen4Lam invariant mass at decay vertex
          } else {
            mStrangeTrack.mMasses[0] = calcMotherMass(p2mom, p2neg, p2pos, PID::Helium3, PID::Pion); // Anti-Hypertriton invariant mass at decay vertex
            mStrangeTrack.mMasses[1] = calcMotherMass(p2mom, p2neg, p2pos, PID::Alpha, PID::Pion);   // Anti-Hyperhydrogen4Lam invariant mass at decay vertex
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
  mDaughterTracks.resize(3); // resize to 3 prongs: first bachelor, second V0 pos, third V0 neg

  for (int iCasc{0}; iCasc < mInputCascadeTracks.size(); iCasc++) {
    LOG(debug) << "Analysing Cascade: " << iCasc + 1 << "/" << mInputCascadeTracks.size();
    auto& DecIndexRef = iCasc;
    auto& casc = mInputCascadeTracks[iCasc];
    auto& cascV0 = mInputV0tracks[casc.getV0ID()];
    mV0dauIDs[0] = cascV0.getProngID(0), mV0dauIDs[1] = cascV0.getProngID(1);

    mStrangeTrack.mPartType = dataformats::kStrkCascade;
    // first: bachelor, second: V0 pos, third: V0 neg
    auto cascR2 = casc.calcR2();
    auto iBinsCasc = mUtils.getBinRect(casc.getEta(), casc.getPhi(), mStrParams->mEtaBinSize, mStrParams->mPhiBinSize);
    for (int& iBinCasc : iBinsCasc) {
      for (int iTrack{mTracksIdxTable[iBinCasc]}; iTrack < TMath::Min(mTracksIdxTable[iBinCasc + 1], int(mSortedITStracks.size())); iTrack++) {
        mStrangeTrack.mMother = (o2::track::TrackParCovF)casc;
        mDaughterTracks[0] = casc.getBachelorTrack(), mDaughterTracks[1] = cascV0.getProng(0), mDaughterTracks[2] = cascV0.getProng(1);
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
          auto p2mom = decayVtxTrackClone.getP2();
          auto p2V0 = mFitter3Body.getTrack(0).getP2();
          auto p2bach = mFitter3Body.getTrack(1).getP2();
          mStrangeTrack.mMasses[0] = calcMotherMass(p2mom, p2V0, p2bach, PID::Lambda, PID::Pion); // Xi invariant mass at decay vertex
          mStrangeTrack.mMasses[1] = calcMotherMass(p2mom, p2V0, p2bach, PID::Lambda, PID::Kaon); // Omega invariant mass at decay vertex

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
    if (!recreateV0(mDaughterTracks[1], mDaughterTracks[2], cascV0Upd)) {
      LOG(debug) << "Cascade V0 refit failed";
      return false;
    }
    cascV0Upd.setAbsCharge(0);
    try {
      nCand = mFitter3Body.process(cascV0Upd, mDaughterTracks[0], motherTrackClone);
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
      nCasc = mFitterV0.process(cascV0Upd, mDaughterTracks[0]);  // refit V0
    } catch (std::runtime_error& e) {
      LOG(debug) << "Cascade refit failed " << e.what();
    }
    if (nCasc) {
      mResettedMotherTrack = mFitterV0.createParentTrackParCov();
    }


    /// ======== KF reconstruction =========
    // check mass of daughters
    float massPosDaughter, massNegDaughter;
    float massBachelor = o2::constants::physics::MassPionCharged;
    if (mDaughterTracks[0].getSign() < 0) { // check charge of bachelor
      massPosDaughter = o2::constants::physics::MassProton;
      massNegDaughter = o2::constants::physics::MassPionCharged;
    }  else {
      massPosDaughter = o2::constants::physics::MassPionCharged;
      massNegDaughter = o2::constants::physics::MassProton;
    }
    /// V0 reconstruction
    KFParticle cascV0KF;
    int nV0Daughters = 2;
    // create KFParticle objects from trackParCovs
    KFParticle kfpDaughter1 = createKFParticleFromTrackParCov(mDaughterTracks[1], mDaughterTracks[1].getSign(), massPosDaughter); // prong 1 (pos)
    KFParticle kfpDaughter2 = createKFParticleFromTrackParCov(mDaughterTracks[2], mDaughterTracks[2].getSign(), massNegDaughter); // prong 2 (neg)
    const KFParticle* V0Daughters[2] = {&kfpDaughter1, &kfpDaughter2};
    // construct mother
    cascV0KF.SetConstructMethod(2);
    cascV0KF.Construct(V0Daughters, nV0Daughters);

    /// cascade reconstruction
    KFParticle cascKF;
    int nCascDaughters = 2;
    // create KFParticle objects from trackParCovs --> CONTINUE HERE!!!
    KFParticle kfpDaughter0 = createKFParticleFromTrackParCov(mDaughterTracks[0], mDaughterTracks[0].getSign(), massBachelor); // bachelor
    const KFParticle* CascDaugthers[2] = {&kfpDaughter0, &cascV0KF};
    // construct mother
    cascKF.SetConstructMethod(2);
    cascKF.Construct(CascDaugthers, nCascDaughters);
    mResettedMotherTrackKF = cascKF;

  }

  // ========== refit V0 ==========
  else if (mStrangeTrack.mPartType == dataformats::kStrkV0) {
    /// ======= DCA fitter reconstruction =======
    try {
      nCand = mFitter3Body.process(mDaughterTracks[0], mDaughterTracks[1], motherTrackClone);
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
      nCasc = mFitterV0.process(mDaughterTracks[0], mDaughterTracks[1]);  // refit V0
    } catch (std::runtime_error& e) {
      LOG(debug) << "Cascade refit failed " << e.what();
    }
    if (nCasc) {
      mResettedMotherTrack = mFitterV0.createParentTrackParCov();
    }

    /// ======== KF reconstruction =========
    // check mass of daughters
    float massPosDaughter, massNegDaughter;
    if (mDaughterTracks[0].getAbsCharge() == 2) { // check charges
      massPosDaughter = o2::constants::physics::MassHelium3;
      massNegDaughter = o2::constants::physics::MassPionCharged;
    }  else {
      massPosDaughter = o2::constants::physics::MassPionCharged;
      massNegDaughter = o2::constants::physics::MassHelium3;
    }
    // V0 reconstruction
    KFParticle V0KF;
    int nV0Daughters = 2;
    // create KFParticle objects from trackParCovs
    KFParticle kfpDaughter0 = createKFParticleFromTrackParCov(mDaughterTracks[0], mDaughterTracks[0].getSign(), massPosDaughter); // prong 1 (pos)
    KFParticle kfpDaughter1 = createKFParticleFromTrackParCov(mDaughterTracks[1], mDaughterTracks[1].getSign(), massNegDaughter); // prong 2 (neg)
    const KFParticle* V0Daughters[2] = {&kfpDaughter0, &kfpDaughter1};
    // construct mother
    V0KF.SetConstructMethod(2);
    V0KF.Construct(V0Daughters, nV0Daughters);
    mResettedMotherTrackKF = V0KF;
  }

  // get vertex position and chi2 of refitted track
  mStrangeTrack.decayVtxKFx = mResettedMotherTrackKF.GetX();
  mStrangeTrack.decayVtxKFy = mResettedMotherTrackKF.GetY();
  mStrangeTrack.decayVtxKFz = mResettedMotherTrackKF.GetZ();
  mStrangeTrack.mGeoChi2KF = mResettedMotherTrackKF.GetChi2();
  mStructClus.arr = nAttachments;
  LOG(info) << "Chi2 DCA fitter: " << mStrangeTrack.mTopoChi2;
  LOG(info) << "Chi2 KF: " << mStrangeTrack.mGeoChi2KF;
  LOG(info) << "Vtx X KF: " << mStrangeTrack.decayVtxKFx;
  LOG(info) << "Vtx Y KF: " << mStrangeTrack.decayVtxKFy;
  LOG(info) << "Vtx Z KF: " << mStrangeTrack.decayVtxKFz;

  mStrangeTrack.mDecayVtx = mFitter3Body.getPCACandidatePos();
  mStrangeTrack.mTopoChi2 = mFitter3Body.getChi2AtPCACandidate();
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