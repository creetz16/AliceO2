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
    mV0dauIDs[kV0DauPos] = v0.getProngID(kV0DauPos);
    mV0dauIDs[kV0DauNeg] = v0.getProngID(kV0DauNeg);
    auto posTrack = v0.getProng(kV0DauPos);
    auto negTrack = v0.getProng(kV0DauNeg);
    auto alphaV0 = calcV0alpha(v0);
    alphaV0 > 0 ? posTrack.setAbsCharge(2) : negTrack.setAbsCharge(2);

    if (!createKFV0(posTrack, negTrack, PID::HyperTriton)) { // reconstruct V0 with KF using Hypertriton PID
      continue;
    }

    o2::track::TrackParCovF correctedV0;
    if (!getTrackParCovFromKFP(kfpMother, PID::HyperTriton, alphaV0 > 0 ? 1 : -1, correctedV0)) { // convert KFParticle V0 to TrackParCov object
        continue;
    }

    mStrangeTrack.mPartType = dataformats::kStrkV0;

    auto v0R2 = v0.calcR2();
    auto iBinsV0 = mUtils.getBinRect(correctedV0.getEta(), correctedV0.getPhi(), mStrParams->mEtaBinSize, mStrParams->mPhiBinSize);
    for (int& iBinV0 : iBinsV0) {
      for (int iTrack{mTracksIdxTable[iBinV0]}; iTrack < TMath::Min(mTracksIdxTable[iBinV0 + 1], int(mSortedITStracks.size())); iTrack++) {
        mDaughterTracks[kV0DauPos] = posTrack;
        mDaughterTracks[kV0DauNeg] = negTrack;
        mITStrack = mSortedITStracks[iTrack];
        auto& ITSindexRef = mSortedITSindexes[iTrack];

        if (mStrParams->mVertexMatching && (mITSvtxBrackets[ITSindexRef].getMin() > v0.getVertexID() ||
                                            mITSvtxBrackets[ITSindexRef].getMax() < v0.getVertexID())) {
          continue;
        }

        if (!propagateToOuterParam(mITStrack, correctedV0)) { // propagate V0 to OuterParam of ITS track
            continue;
        }
        mStrangeTrack.mMother = correctedV0; // V0 at OuterParam of ITS track

        if (matchDecayToITStrack(sqrt(v0R2))) { // now mother at IU

          LOG(debug) << "ITS Track matched with a V0 decay topology ....";
          LOG(debug) << "Number of ITS track clusters attached: " << mITStrack.getNumberOfClusters();

          // get parameters at IU
          mStrangeTrack.mMother.getPxPyPzGlo(mStrangeTrack.mIUMom);

          // get parameters at decay vertex
          /// TODO: function --> fill KF info
          mStrangeTrack.mDecayMom[0] = kfpMother.GetPx();
          mStrangeTrack.mDecayMom[1] = kfpMother.GetPy();
          mStrangeTrack.mDecayMom[2] = kfpMother.GetPz();
          mStrangeTrack.mDecayPt = kfpMother.GetPt();
          mStrangeTrack.mDecayVtx[0] = kfpMother.GetX();
          mStrangeTrack.mDecayVtx[1] = kfpMother.GetY();
          mStrangeTrack.mDecayVtx[2] = kfpMother.GetZ();
          mStrangeTrack.mGeoChi2 = kfpMother.GetChi2();

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
    auto posTrack = cascV0.getProng(kV0DauPos);
    auto negTrack = cascV0.getProng(kV0DauNeg);
    auto bachTrack = casc.getBachelorTrack();

    if (!createKFCascade(posTrack, negTrack, bachTrack, PID::XiMinus)) { // reconstruct cascade with KF using XiMinus PID
      continue;
    }

    o2::track::TrackParCovF cascade;
    if (!getTrackParCovFromKFP(kfpMother, PID::XiMinus, bachTrack.getCharge()<0 ? -1 : 1, cascade)) { // convert KFParticle cascade to TrackParCov object
        continue;
    }

    mStrangeTrack.mPartType = dataformats::kStrkCascade;

    auto cascR2 = casc.calcR2();
    auto iBinsCasc = mUtils.getBinRect(casc.getEta(), casc.getPhi(), mStrParams->mEtaBinSize, mStrParams->mPhiBinSize);
    for (int& iBinCasc : iBinsCasc) {
      for (int iTrack{mTracksIdxTable[iBinCasc]}; iTrack < TMath::Min(mTracksIdxTable[iBinCasc + 1], int(mSortedITStracks.size())); iTrack++) {
        mDaughterTracks[kV0DauPos] = posTrack;
        mDaughterTracks[kV0DauNeg] = negTrack;
        mDaughterTracks[kBach] = bachTrack;
        mITStrack = mSortedITStracks[iTrack];
        auto& ITSindexRef = mSortedITSindexes[iTrack];
        LOG(debug) << "----------------------";
        LOG(debug) << "CascV0: " << casc.getV0ID() << ", Bach ID: " << casc.getBachelorID() << ", ITS track ref: " << mSortedITSindexes[iTrack];

        if (mStrParams->mVertexMatching && (mITSvtxBrackets[ITSindexRef].getMin() > casc.getVertexID() ||
                                            mITSvtxBrackets[ITSindexRef].getMax() < casc.getVertexID())) {
          LOG(debug) << "Vertex ID mismatch: " << mITSvtxBrackets[ITSindexRef].getMin() << " < " << casc.getVertexID() << " < " << mITSvtxBrackets[ITSindexRef].getMax();
          continue;
        }

        if (!propagateToOuterParam(mITStrack, cascade)) { // propagate cascade to OuterParam of ITS track
            continue;
        }
        mStrangeTrack.mMother = cascade; // cascade at OuterParam of ITS track

        if (matchDecayToITStrack(sqrt(cascR2))) // now mother is at IU
        { 

          LOG(debug) << "ITS Track matched with a Cascade decay topology ....";
          LOG(debug) << "Number of ITS track clusters attached: " << mITStrack.getNumberOfClusters();

          // get parameters at IU
          mStrangeTrack.mMother.getPxPyPzGlo(mStrangeTrack.mIUMom);

          // get parameters at decay vertex
          /// TODO: put mass here
          mStrangeTrack.mDecayMom[0] = kfpMother.GetPx();
          mStrangeTrack.mDecayMom[1] = kfpMother.GetPy();
          mStrangeTrack.mDecayMom[2] = kfpMother.GetPz();
          mStrangeTrack.mDecayPt = kfpMother.GetPt();
          mStrangeTrack.mDecayVtx[0] = kfpMother.GetX();
          mStrangeTrack.mDecayVtx[1] = kfpMother.GetY();
          mStrangeTrack.mDecayVtx[2] = kfpMother.GetZ();
          mStrangeTrack.mGeoChi2 = kfpMother.GetChi2();

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

    // construct updated decay vertex with KF
    std::array<float, 3> xyz_tr, pxpypz_tr;
    std::array<float, 3> xyz_hyd, pxpypz_hyd;
    std::array<float, 3> deltaXYZ, deltaPxPyPz;
    std::array<float, 21> cv_tr, cv_hyd;
    std::array<float, 21> deltaCV;

    if (mStrangeTrack.mPartType == dataformats::kStrkV0) {
      if (!createKFV0(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], PID::Hyperhydrog4)) {
        return false;
      }
      if (!getTrackParCovFromKFP(kfpMother, PID::HyperTriton, mDaughterTracks[kV0DauPos].getAbsCharge()==2 ? 1 : -1, mStrangeTrack.mMother)) {
        return false;
      }
      mStrangeTrack.mMother.getPxPyPzGlo(pxpypz_hyd);
      mStrangeTrack.mMother.getXYZGlo(xyz_hyd);
      mStrangeTrack.mMother.getCovXYZPxPyPzGlo(cv_hyd);

      if (!createKFV0(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], PID::HyperTriton)) {
        return false;
      }
      if (!getTrackParCovFromKFP(kfpMother, PID::HyperTriton, mDaughterTracks[kV0DauPos].getAbsCharge()==2 ? 1 : -1, mStrangeTrack.mMother)) {
        return false;
      }
      mStrangeTrack.mMother.getPxPyPzGlo(pxpypz_tr);
      mStrangeTrack.mMother.getXYZGlo(xyz_tr);
      mStrangeTrack.mMother.getCovXYZPxPyPzGlo(cv_tr);

      for (int i=0; i<3; i++) {
        deltaXYZ[i] = xyz_tr[i] - xyz_hyd[i];
        deltaPxPyPz[i] = pxpypz_tr[i] - pxpypz_hyd[i];
      }

      for (int i=0; i<21; i++) {
        deltaCV[i] = cv_tr[i] - cv_hyd[i];
      }

      // LOG(info) << "################################################";
      // LOG(info) << "Delta x        " << deltaXYZ[0];
      // LOG(info) << "Delta y        " << deltaXYZ[1];
      // LOG(info) << "Delta z        " << deltaXYZ[2];
      // LOG(info) << "Delta Px       " << deltaPxPyPz[0];
      // LOG(info) << "Delta Py       " << deltaPxPyPz[1];
      // LOG(info) << "Delta Pz       " << deltaPxPyPz[2];
      if (deltaCV[0]>1e-5 || deltaCV[1]>1e-5 || deltaCV[2]>1e-5 || deltaCV[3]>1e-5 || deltaCV[4]>1e-5 || deltaCV[5]>1e-5 || deltaCV[6]>1e-5 || deltaCV[7]>1e-5 || deltaCV[8]>1e-5 || 
      deltaCV[9]>1e-5 || deltaCV[10]>1e-5 || deltaCV[11]>1e-5 || deltaCV[12]>1e-5 || deltaCV[13]>1e-5 || deltaCV[14]>1e-5 || deltaCV[15]>1e-5 || deltaCV[16]>1e-5 || deltaCV[17]>1e-5 || 
      deltaCV[18]>1e-5 || deltaCV[19]>1e-5 || deltaCV[20]>1e-5) {
        for (int i=0; i<21; i++) {
          LOG(info) << "Delta cov[" << i << "]   " << deltaCV[i];
        }
        LOG(info) << "################################################";
      }
    }

    else if (mStrangeTrack.mPartType == dataformats::kStrkCascade) {
      if (!createKFCascade(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], mDaughterTracks[kBach], PID::XiMinus)) {
        return false;
      }
      if (!getTrackParCovFromKFP(kfpMother, PID::XiMinus, mDaughterTracks[kBach].getCharge()<0 ? -1 : 1, mStrangeTrack.mMother)) {
        return false;
      }
    }
  }

  if (nUpdates < trackClusters.size() || motherClusters.size() < nMinClusMother) {
    return false;
  }

  // reconstruct mother at decay vertex with hypertriton/hyperhydrogen and Xi/Omega mass hypotheses
  float M, SigmaM;
  if (mStrangeTrack.mPartType == dataformats::kStrkV0) {
    // hyperhydrogen
    if (!createKFV0(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], PID::Hyperhydrog4)) {
      return false;
    }
    // get invariant mass at decay vertex
    kfpMother.GetMass(M, SigmaM);
    mStrangeTrack.mMasses[1] = M;

    // hypertriton
    if (!createKFV0(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], PID::HyperTriton)) {
      return false;
    }
    // get invariant mass at decay vertex
    kfpMother.GetMass(M, SigmaM);
    mStrangeTrack.mMasses[0] = M;
  }
  else if (mStrangeTrack.mPartType == dataformats::kStrkCascade) {
    // OmegaMinus
    if (!createKFCascade(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], mDaughterTracks[kBach], PID::OmegaMinus)) {
      return false;
    }
    // get invariant mass at decay vertex
    kfpMother.GetMass(M, SigmaM);
    mStrangeTrack.mMasses[1] = M;

    // XiMinus
    if (!createKFCascade(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], mDaughterTracks[kBach], PID::XiMinus)) {
      return false;
    }
    // get invariant mass at decay vertex
    kfpMother.GetMass(M, SigmaM);
    mStrangeTrack.mMasses[0] = M;
  }

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