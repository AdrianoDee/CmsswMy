#ifndef CombinedHitQuadrupletGeneratorCA_H
#define CombinedHitQuadrupletGeneratorCA_H

#include <vector>
#include "RecoPixelVertexing/PixelTriplets/interface/HitQuadrupletGenerator.h"
#include "RecoTracker/TkHitPairs/interface/LayerFKDTreeCache.h"
#include "RecoTracker/TkHitPairs/interface/LayerDoubletsCache.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include <string>
#include <memory>

class TrackingRegion;
class CANtupleGenerator;
class SeedingLayerSetsHits;

namespace edm { class Event; class EventSetup; }

class CombinedHitQuadrupletGeneratorCA : public HitQuadrupletGenerator {
public:

  CombinedHitQuadrupletGeneratorCA( const edm::ParameterSet& cfg, edm::ConsumesCollector& iC);
    
  virtual ~CombinedHitQuadrupletGeneratorCA();
    
  virtual void hitQuadruplets( const TrackingRegion& reg, OrderedHitSeeds & triplets,
                                const edm::Event & ev,  const edm::EventSetup& es);
    
  virtual void clear() final;

private:
    
  edm::EDGetTokenT<SeedingLayerSetsHits> theSeedingLayerToken;
    
  std::unique_ptr<CANtupleGenerator> theGenerator;
    
  LayerFKDTreeCache           theKDTReeCache;
  LayerDoubletsCache          theDoubletsCache; 
    
};

#endif
