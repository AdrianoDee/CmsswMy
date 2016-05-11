#include "RecoPixelVertexing/PixelTriplets/interface/CANtupleGenerator.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

CANtupleGenerator::CANtupleGenerator(unsigned int maxElement):
theKDTreeCache(nullptr),
theDoubletsCache(nullptr),
theMaxElement(maxElement)
{}

CANtupleGenerator::CANtupleGenerator(const edm::ParameterSet& pset):
CANtupleGenerator(pset.getParameter<unsigned int>("maxElement"))
{}

CANtupleGenerator::~CANtupleGenerator() {}

void CANtupleGenerator::init(std::unique_ptr<HitPairGeneratorFromLayerPairCA>&& doubletsGenerator,LayerFKDTreeCache *kdTReeCache,LayerDoubletsCache *doubletsCache) {
    theDoubletsGeneratorCA =  doubletsGenerator;
    theKDTreeCache = kdTReeCache;
    theDoubletsCache = doubletsCache;
}
