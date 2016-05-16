import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEParmError_cfi import *
from RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi import *
from RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cfi import *
from TrackingTools.MaterialEffects.MaterialPropagator_cfi import *
from RecoTracker.TkSeedingLayers.TTRHBuilderWithoutAngle4MixedTriplets_cfi import *
from RecoTracker.TkSeedingLayers.TTRHBuilderWithoutAngle4MixedPairs_cfi import *
#from RecoTracker.TkSeedingLayers.MixedLayerTriplets_cfi import *
from RecoTracker.TkSeedingLayers.PixelLayerQuadrupletsCA_cfi import *
from RecoPixelVertexing.PixelTriplets.CANtupleHLTGenerator_cfi import *
#from RecoPixelVertexing.PixelTriplets.PixelTripletLargeTipGenerator_cfi import *

import RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi
globalSeedsFromNTuplets = RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi.seedGeneratorFromRegionHitsEDProducer.clone(
    OrderedHitsFactoryPSet = cms.PSet(
      ComponentName = cms.string('StandardCANtupletGenerator'),
      SeedingLayers = cms.InputTag('PixelLayerQuadrupletsCA'),
      GeneratorPSet = cms.PSet(CANtupleHLTGenerator_cfi.clone(maxElement = cms.uint32(1000000)))
# this one uses an exact helix extrapolation and can deal correctly with
# arbitrarily large D0 and generally exhibits a smaller fake rate:
#     GeneratorPSet = cms.PSet(PixelTripletLargeTipGenerator)
    )
)


