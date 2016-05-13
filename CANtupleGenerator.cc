#include "RecoPixelVertexing/PixelTriplets/interface/CANtupleGenerator.h"
#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPairCA.h"

CANtupleGenerator::CANtupleGenerator():
theKDTreeCache(nullptr),
theDoubletsCache(nullptr),
{}

CANtupleGenerator::~CANtupleGenerator() {}

void CANtupleGenerator::init(std::unique_ptr<HitPairGeneratorFromLayerPairCA>&& doubletsGenerator,LayerFKDTreeCache *kdTReeCache,LayerDoubletsCache *doubletsCache) {
    theDoubletsGeneratorCA =  std::move(doubletsGenerator);
    theKDTreeCache = kdTReeCache;
    theDoubletsCache = doubletsCache;
}
