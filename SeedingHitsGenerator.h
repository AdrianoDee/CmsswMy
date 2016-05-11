#ifndef TkTrackingRegions_SeedingHitsGenerator_H
#define TkTrackingRegions_SeedingHitsGenerator_H

#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"
#include <vector>

class TrackingRegion;
namespace edm { class Event; class EventSetup; class ConsumesCollector;}

class SeedingHitsGenerator {
public:
  SeedingHitsGenerator() : theMaxElement(0){}
  virtual ~SeedingHitsGenerator() {}

  virtual const SeedingHitSet & run(
      const TrackingRegion& reg, const edm::Event & ev, const edm::EventSetup& es ) = 0;

  virtual void clear() { }  //fixme: should be purely virtual!

  unsigned int theMaxElement;
};

#endif
