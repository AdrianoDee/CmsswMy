#ifndef RecoTracker_TkTrackingRegions_SeedingHitsGeneratorFactory_H
#define RecoTracker_TkTrackingRegions_SeedingHitsGeneratorFactory_H

#include "RecoTracker/TkTrackingRegions/interface/SeedingHitsGenerator.h"
namespace edm {class ParameterSet; class ConsumesCollector;}

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SeedingHitsGenerator *(const edm::ParameterSet &, edm::ConsumesCollector& iC)> SeedingHitsGeneratorFactory;

#endif

