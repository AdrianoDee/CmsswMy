#ifndef CANtupleHLTGenerator_H
#define CANtupleHLTGenerator_H

#include "CombinedHitQuadrupletGeneratorCA.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoUtilities.h"
#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPairCA.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

//#include "RecoPixelVertexing/PixelTriplets/interface/HitQuadrupletGenerator.h"
#include "RecoPixelVertexing/PixelTriplets/interface/OrderedHitSeeds.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoPixelVertexing/PixelTriplets/interface/CANtupleGenerator.h"

#include <utility>
#include <vector>

class SeedComparitor;

class CANtupleHLTGenerator : public CANtupleGenerator {

public:
  CANtupleHLTGenerator( const edm::ParameterSet& cfg, edm::ConsumesCollector& iC);

  virtual ~CANtupleHLTGenerator();

    virtual void getNTuplets(const TrackingRegion& region, OrderedHitSeeds & ntuplets,
                                const edm::Event & ev, const edm::EventSetup& es,
                                SeedingLayerSetsHits::SeedingLayerSet fourLayers
                                ) override;

private:
  const float extraHitRZtolerance;
  const float extraHitRPhitolerance;
  const QuantityDependsPt extraPhiTolerance;
  const QuantityDependsPt maxChi2;
  const bool fitFastCircle;
  const bool fitFastCircleChi2Cut;
  const bool useBendingCorrection;
  const bool useFixedPreFiltering;
  const bool useMScat;
  const bool useBend;
  const float dphi;
  std::unique_ptr<SeedComparitor> theComparitor;

};
#endif


