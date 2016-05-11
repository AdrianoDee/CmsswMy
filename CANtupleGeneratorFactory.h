#ifndef CANtupleGeneratorFactory_H
#define CANtupleGeneratorFactory_H

#include "RecoPixelVertexing/PixelTriplets/interface/HitQuadrupletGenerator.h"
#include "RecoPixelVertexing/PixelTriplets/interface/CANtupleGenerator.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"

namespace edm {class ParameterSet; class ConsumesCollector;}

typedef edmplugin::PluginFactory<CANtupleGenerator *(const edm::ParameterSet &, edm::ConsumesCollector&)>
	CANtupleGeneratorFactory;
 
#endif
