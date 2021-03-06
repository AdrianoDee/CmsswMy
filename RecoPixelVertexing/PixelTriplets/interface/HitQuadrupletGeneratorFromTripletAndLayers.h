#ifndef RecoPixelVertexing_PixelTriplets_HitQuadrupletGeneratorFromTripletAndLayers_h
#define RecoPixelVertexing_PixelTriplets_HitQuadrupletGeneratorFromTripletAndLayers_h

/** A HitQuadrupletGenerator from HitTripletGenerator and vector of
    Layers. The HitTripletGenerator provides a set of hit triplets.
    For each triplet the search for compatible hit(s) is done among
    provided Layers
 */

#include "RecoPixelVertexing/PixelTriplets/interface/OrderedHitSeeds.h"
#include <vector>
#include "TrackingTools/TransientTrackingRecHit/interface/SeedingLayerSetsHits.h"
#include "RecoTracker/TkHitPairs/interface/LayerHitMapCache.h"
#include "RecoTracker/TkHitPairs/interface/LayerFKDTreeCache.h"
#include "RecoTracker/TkHitPairs/interface/LayerDoubletsCache.h"

class HitTripletGeneratorFromPairAndLayers;

class HitQuadrupletGeneratorFromTripletAndLayers {

public:
    typedef LayerHitMapCache    LayerCacheType;
    typedef LayerFKDTreeCache   LayerTreeCacheType;
    typedef LayerDoubletsCache  LayerDoubletsCacheType;

  HitQuadrupletGeneratorFromTripletAndLayers();
  virtual ~HitQuadrupletGeneratorFromTripletAndLayers();

  void init( std::unique_ptr<HitTripletGeneratorFromPairAndLayers>&& tripletGenerator, LayerCacheType* layerCache);
    
  void init( std::unique_ptr<HitTripletGeneratorFromPairAndLayers>&& tripletGenerator, LayerCacheType* layerCache, LayerTreeCacheType* treeCache, LayerDoubletsCacheType* doubletsCache);

  virtual void hitQuadruplets( const TrackingRegion& region, OrderedHitSeeds& result,
                               const edm::Event& ev, const edm::EventSetup& es,
                               const SeedingLayerSetsHits::SeedingLayerSet& tripletLayers,
                               const std::vector<SeedingLayerSetsHits::SeedingLayer>& fourthLayers) = 0;

  virtual void hitQuadruplets( const TrackingRegion& region, OrderedHitSeeds& result,
                                const edm::Event& ev, const edm::EventSetup& es,
                                const SeedingLayerSetsHits::SeedingLayerSet& fourLayers) = 0;
    

protected:
  std::unique_ptr<HitTripletGeneratorFromPairAndLayers> theTripletGenerator;

  LayerCacheType            *theLayerCache;
  LayerTreeCacheType        *theKDTreeCache;
  LayerDoubletsCacheType    *theDoubletsCache;

  
  
};
#endif

