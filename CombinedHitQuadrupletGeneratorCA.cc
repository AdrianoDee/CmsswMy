#include "CombinedHitQuadrupletGeneratorCA.h"

#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPairCA.h"
#include "RecoPixelVertexing/PixelTriplets/interface/CANtupleGenerator.h"
#include "RecoPixelVertexing/PixelTriplets/interface/CANtupleGeneratorFactory.h"
#include "LayerQuadruplets.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

using namespace std;
using namespace ctfseeding;

CombinedHitQuadrupletGeneratorCA::CombinedHitQuadrupletGeneratorCA(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC):
theSeedingLayerToken(iC.consumes<SeedingLayerSetsHits>(cfg.getParameter<edm::InputTag>("SeedingLayers")))
{
    edm::ParameterSet generatorPSet = cfg.getParameter<edm::ParameterSet>("GeneratorPSet");
    std::string       generatorName = generatorPSet.getParameter<std::string>("ComponentName");
    edm::ParameterSet cantupleGeneratorPSet = cfg.getParameter<edm::ParameterSet>("CAGeneratorPSet");
    std::string cantupleGeneratorName = cantupleGeneratorPSet.getParameter<std::string>("ComponentName");
    
    theGenerator.reset(CANtupleGeneratorFactory::get()->create(generatorName, generatorPSet, iC));
    theGenerator->init(std::make_unique<HitPairGeneratorFromLayerPairCA>(0, 1), &theKDTReeCache, &theDoubletsCache);
}

CombinedHitQuadrupletGeneratorCA::~CombinedHitQuadrupletGeneratorCA() {}

void CombinedHitQuadrupletGeneratorCA::hitQuadruplets(
                                                    const TrackingRegion& region, OrderedHitSeeds & result,
                                                    const edm::Event& ev, const edm::EventSetup& es)
{
    edm::Handle<SeedingLayerSetsHits> hlayers;
    ev.getByToken(theSeedingLayerToken, hlayers);
    const SeedingLayerSetsHits& layers = *hlayers;
    if(layers.numberOfLayersInSet() != 4)
        throw cms::Exception("Configuration") << "CombinedHitQuadrupletsGenerator expects SeedingLayerSetsHits::numberOfLayersInSet() to be 4, got " << layers.numberOfLayersInSet();
    
    for(int j=0; j<(int)layers->size();j++) {
        theGenerator->getNTuplets(region, result, ev, es, layers[j]);
    }
    
    theCACellsCache.clear();
    theKDTReeCache.clear();

}


