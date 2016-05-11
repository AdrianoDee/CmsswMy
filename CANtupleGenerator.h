#ifndef CANtupleGenerator_H
#define CANtupleGenerator_H

//#include "RecoPixelVertexing/PixelTriplets/interface/OrderedHitQuadruplets.h"
//#include "RecoPixelVertexing/PixelTriplets/interface/OrderedHitTriplets.h"
#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPairCA.h"
#include "RecoTracker/TkSeedingLayers/interface/OrderedSeedingHits.h"
#include "TrackingTools/TransientTrackingRecHit/interface/SeedingLayerSetsHits.h"
#include "RecoTracker/TkHitPairs/interface/LayerDoubletsCache.h"
#include "RecoTracker/TkHitPairs/interface/LayerFKDTreeCache.h"
#include <vector>

namespace edm { class ParameterSet; class Event; class EventSetup; class ConsumesCollector; }
class TrackingRegion;
class HitPairGeneratorFromLayerPairCA;

class CANtupleGenerator {

public:

  explicit CANtupleGenerator(unsigned int maxElement=0);
  explicit CANtupleGenerator(const edm::ParameterSet& pset);
  virtual ~CANtupleGenerator();

  void init(std::unique_ptr<HitPairGeneratorFromLayerPairCA>&& doubletsGenerator,LayerFKDTreeCache *kdTReeCache,LayerDoubletsCache *doubletsCache);

  virtual void getNTuplets(const TrackingRegion& region, SeedingHitSet & ntuplets,
                             const edm::Event & ev, const edm::EventSetup& es,
                             SeedingLayerSetsHits::SeedingLayerSet fourLayers
                             ) = 0;
protected:
    
  LayerFKDTreeCache *theKDTreeCache;
  LayerDoubletsCache *theDoubletsCache;
    
  std::unique_ptr<HitPairGeneratorFromLayerPairCA>&& theDoubletsGeneratorCA;
  
  const unsigned int theMaxElement;
};
#endif


