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

    std::array<float, 3> xyzGloV0;
    std::array<float, 21> covV0;
    v0.getXYZGlo(xyzGloV0);
    v0.getCovXYZPxPyPzGlo(covV0);
    LOG(debug) << "Position of V0 before KFParticle creation: (" << v0.getX() << ", " << v0.getY() << ", " << v0.getZ() << ")";
    LOG(debug) << "Glo position of V0 before KFParticle creation: (" << xyzGloV0[0] << ", " << xyzGloV0[1] << ", " << xyzGloV0[2] << ")";
    // fill Glo position in histogram
    float motherZ = xyzGloV0[2];
    float motherR = sqrt(xyzGloV0[0] * xyzGloV0[0] + xyzGloV0[1] * xyzGloV0[1]);
    float motherErrZ = sqrt(fabs(v0.getSigmaZ2()));
    
    if (!createKFV0(posTrack, negTrack, pidV0)) { // reconstruct V0 with KF using Hypertriton PID // PID::HyperTriton
      continue;
    }
    
    float motherZkf = kfpMother.GetZ();
    float motherRkf = sqrt(kfpMother.GetX() * kfpMother.GetX() + kfpMother.GetY() * kfpMother.GetY());
    float motherErrZkf = kfpMother.GetErrZ();

    float M, SigmaM;
    kfpMother.GetMass(M, SigmaM);
    mStrangeTrack.mMassInit = M;

    o2::track::TrackParCovF correctedV0;
    if (!getTrackParCovFromKFP(kfpMother, pidV0, alphaV0 > 0 ? 1 : -1, correctedV0)) { // convert KFParticle V0 to TrackParCov object
        continue;
    }

    mStrangeTrack.mPartType = dataformats::kStrkV0;

    mStrangeTrack.mMotherZ = motherZ;
    mStrangeTrack.mMotherZkf = motherZkf;
    mStrangeTrack.mMotherR = motherR;
    mStrangeTrack.mMotherRkf = motherRkf;
    mStrangeTrack.mMotherErrZ = motherErrZ;
    mStrangeTrack.mMotherErrZkf = motherErrZkf;

    mStrangeTrack.mGeoChi2KFcreation = kfpMother.GetChi2();

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

        std::array<float, 3> OuterParamXYZ;
        mITStrack.getParamOut().getXYZGlo(OuterParamXYZ);
        LOG(info) << "Outer param radius of ITS track: " << sqrt(OuterParamXYZ[0] * OuterParamXYZ[0] + OuterParamXYZ[1] * OuterParamXYZ[1]);

        if (matchDecayToITStrack(sqrt(v0R2))) { // now mother at IU

          LOG(debug) << "ITS Track matched with a V0 decay topology ....";
          LOG(debug) << "Number of ITS track clusters attached: " << mITStrack.getNumberOfClusters();

          // get parameters at IU
          mStrangeTrack.mMother.getPxPyPzGlo(mStrangeTrack.mIUMom);

          // get parameters at decay vertex
          if (!fillKFinfo(kfpMother)) {
            continue;
          }

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
    
    std::array<float, 3> xyzGloCasc;
    std::array<float, 21> covCasc;
    casc.getXYZGlo(xyzGloCasc);
    casc.getCovXYZPxPyPzGlo(covCasc);
    LOG(debug) << "Position of cascade before KFParticle creation: (" << casc.getX() << ", " << casc.getY() << ", " << casc.getZ() << ")";
    LOG(debug) << "Glo position of cascade before KFParticle creation: (" << xyzGloCasc[0] << ", " << xyzGloCasc[1] << ", " << xyzGloCasc[2] << ")";
    // get Glo z position and radius before KFParticle creation
    float motherZ = xyzGloCasc[2];
    float motherR = sqrt(xyzGloCasc[0] * xyzGloCasc[0] + xyzGloCasc[1] * xyzGloCasc[1]);
    float motherErrZ = sqrt(fabs(casc.getSigmaZ2()));
    if (!createKFCascade(posTrack, negTrack, bachTrack, pidCasc, false)) { // reconstruct cascade with KF using XiMinus PID // PID::XiMinus
      continue;
    }

    float motherZkf = kfpMother.GetZ();
    float motherRkf = sqrt(kfpMother.GetX() * kfpMother.GetX() + kfpMother.GetY() * kfpMother.GetY());
    float motherErrZkf = kfpMother.GetErrZ();

    float M, SigmaM;
    kfpMother.GetMass(M, SigmaM);
    mStrangeTrack.mMassInit = M;

    o2::track::TrackParCovF cascade;
    if (!getTrackParCovFromKFP(kfpMother, pidCasc, bachTrack.getCharge()<0 ? -1 : 1, cascade)) { // convert KFParticle cascade to TrackParCov object // PID::XiMinus
        continue;
    }

    mStrangeTrack.mPartType = dataformats::kStrkCascade;

    mStrangeTrack.mMotherZ = motherZ;
    mStrangeTrack.mMotherZkf = motherZkf;
    mStrangeTrack.mMotherR = motherR;
    mStrangeTrack.mMotherRkf = motherRkf;
    mStrangeTrack.mMotherErrZ = motherErrZ;
    mStrangeTrack.mMotherErrZkf = motherErrZkf;

    mStrangeTrack.mGeoChi2KFcreation = kfpMother.GetChi2();

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

        std::array<float, 3> OuterParamXYZ;
        mITStrack.getParamOut().getXYZGlo(OuterParamXYZ);
        LOG(info) << "Outer param radius of ITS track: " << sqrt(OuterParamXYZ[0] * OuterParamXYZ[0] + OuterParamXYZ[1] * OuterParamXYZ[1]);

        if (matchDecayToITStrack(sqrt(cascR2))) // now mother is at IU
        { 

          LOG(debug) << "ITS Track matched with a Cascade decay topology ....";
          LOG(debug) << "Number of ITS track clusters attached: " << mITStrack.getNumberOfClusters();

          // get parameters at IU
          mStrangeTrack.mMother.getPxPyPzGlo(mStrangeTrack.mIUMom);

          // get parameters at decay vertex
          if (!fillKFinfo(kfpMother)) {
            continue;
          }

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

  LOG(info) << "##########################";
  int counterBad = 0;
  std::array<float, 3> mom;
  std::array<float, 21> cv;
  mStrangeTrack.mMother.getCovXYZPxPyPzGlo(cv);
  mStrangeTrack.mMother.getPxPyPzGlo(mom);
  LOG(debug) << "Mother momentum before update = (" << mom[0] << ", " << mom[1] << ", " << mom[2] << ")";
  LOG(debug) << "Mother squared unc on 1/pT before update = " << cv[20];
  float deltaCV = cv[20];
  for (int iClus{0}; iClus < trackClusters.size(); iClus++) {
    auto& clus = trackClusters[iClus];
    auto& compClus = trackClusSizes[iClus];
    int nUpdOld = nUpdates;
    double clusRad = sqrt(clus.getX() * clus.getX() + clus.getY() * clus.getY());
    auto decayRupd = sqrt(kfpMother.GetX() * kfpMother.GetX() + kfpMother.GetY() * kfpMother.GetY());
    auto diffR = decayRupd - clusRad;
    auto relDiffR = diffR / decayRupd;
    // print position of mother
    std::array<float, 3> xyzMother;
    mStrangeTrack.mMother.getXYZGlo(xyzMother);
    double radius2Mother = xyzMother[0] * xyzMother[0] + xyzMother[1] * xyzMother[1];
    LOG(info) << "Glo mother radius from kfpMother TrackParCov (initially at outer param of ITS track): " << sqrt(radius2Mother);
    // Look for the Mother if the Decay radius allows for it, within a tolerance
    LOG(info) << "decayR: " << decayR << ", decayRupd: " << decayRupd << ", clus rad: " << clusRad << ", diffR: " << diffR << ", relDiffR: " << relDiffR << ", radTol: " << radTol;

    if (relDiffR > -radTol) {
      LOG(debug) << "Try to attach cluster to Mother, layer: " << geom->getLayer(clus.getSensorID());

      if (updateTrack(clus, mStrangeTrack.mMother)) {
        motherClusters.push_back(clus);
        motherClusSizes.push_back(compClus);
        nAttachments[geom->getLayer(clus.getSensorID())] = 0;
        isMotherUpdated = true;
        nUpdates++;
        LOG(info) << "Cluster attached to Mother";

        // print momentumm change
        std::array<float, 3> momup;
        std::array<float, 21> cvup;
        mStrangeTrack.mMother.getCovXYZPxPyPzGlo(cvup);
        mStrangeTrack.mMother.getPxPyPzGlo(momup);
        LOG(debug) << "Updated mother momentum = (" << momup[0] << ", " << momup[1] << ", " << momup[2] << ")";
        deltaCV = abs(deltaCV - cvup[20]);
        LOG(debug) << "Difference in squared unc on 1/pT = " << deltaCV;
        if (deltaCV > 1e-2) counterBad++;

        continue; // if the cluster is attached to the mother, skip the rest of the loop
      }
    }

    // if Mother is not found, check for daughters compatibility
    if (relDiffR < radTol && !isMotherUpdated) {
      bool isDauUpdated = false;
      LOG(debug) << "Try to attach cluster to Daughters, layer: " << geom->getLayer(clus.getSensorID());
      for (int iDau{0}; iDau < mDaughterTracks.size(); iDau++) {
        auto& dauTrack = mDaughterTracks[iDau];

        if (updateTrack(clus, dauTrack)) {
          nAttachments[geom->getLayer(clus.getSensorID())] = iDau + 1;
          isDauUpdated = true;
          LOG(info) << "Cluster attached to Daughter";
          break;
        }
      }
      if (!isDauUpdated) {
        LOG(info) << "No daughter track updated.";
        break; // no daughter track updated, stop the loop
      }
      nUpdates++;
      
    }
    if (nUpdates == nUpdOld) {
      LOG(info) << "No track updated.";
      break; // no track updated, stop the loop
    }


    if (mStrangeTrack.mPartType == dataformats::kStrkV0) {

      // create hyperhydrogen4 V0 and fill mass
      if (!createKFV0(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], pidV0comp)) {
        return false;
      }
      float M, SigmaM;
      kfpMother.GetMass(M, SigmaM);
      mStrangeTrack.mMasses[1] = M;

      // recreate hypertriton V0
      if (!createKFV0(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], pidV0)) {
        return false;
      }
      // print mommentum change
      LOG(debug) << "Updated mother momentum = (" << kfpMother.GetPx() << ", " << kfpMother.GetPy() << ", " << kfpMother.GetPz() << ")";
      float cv20up = kfpMother.GetCovariance(20);
      deltaCV = abs(deltaCV - cv20up);
      if (deltaCV > 1e-2) counterBad++;
      LOG(debug) << "Difference in squared unc on 1/pT = " << deltaCV;

      // get TrackParCov from KFParticle mother
      if (!getTrackParCovFromKFP(kfpMother, pidV0, mDaughterTracks[kV0DauPos].getAbsCharge()==2 ? 1 : -1, mStrangeTrack.mMother)) {
        return false;
      }

    }

    else if (mStrangeTrack.mPartType == dataformats::kStrkCascade) {

      // create Omega cascade and fill mass
      if (!createKFCascade(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], mDaughterTracks[kBach], pidCascComp, false)) {
        return false;
      }
      float M, SigmaM;
      kfpMother.GetMass(M, SigmaM);
      mStrangeTrack.mMasses[1] = M;

      // recreate XiMinus cascade
      if (!createKFCascade(mDaughterTracks[kV0DauPos], mDaughterTracks[kV0DauNeg], mDaughterTracks[kBach], pidCasc, true)) {
        return false;
      }

      // get TrackParCov from KFParticle mother
      if (!getTrackParCovFromKFP(kfpMother, pidCasc, mDaughterTracks[kBach].getCharge()<0 ? -1 : 1, mStrangeTrack.mMother)) {
        return false;
      }
    }
  } // end of cluster loop
  LOG(info) << "##########################";

  if (nUpdates < trackClusters.size() || motherClusters.size() < nMinClusMother) {
    return false;
  }

  mStructClus.arr = nAttachments;

  LOG(info) << "Number of updates: " << nUpdates;
  LOG(debug) << "Number of times deltaCV > 1e-2: " << counterBad;

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

  // print position
  double StartX = track.getX();
  double StartY = track.getY();
  double StartZ = track.getZ();
  double radius2 = StartX * StartX + StartY * StartY;
  if (sqrt(radius2) > 40) {
    LOG(info) << "Radius outside ITS: " << sqrt(radius2);
  }
  if (abs(StartZ) > 75) {
    LOG(info) << "z-position outside ITS: " << StartZ;
  }

  if (!propInstance->propagateToX(track, x, getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mCorrType)) {
    LOG(info) << "Propagation failed.";
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