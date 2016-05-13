import FWCore.ParameterSet.Config as cms

from RecoTracker.TkSeedingLayers.seedingLayersEDProducer_cfi import *

PixelLayerQuadrupletsCA = seedingLayersEDProducer.clone()
PixelLayerQuadrupletsCA.layerList = cms.vstring('BPix1+BPix2+BPix3+BPix4',
                                           'BPix1+BPix2+BPix3+FPix1_pos',
                                           'BPix1+BPix2+BPix3+FPix1_neg',
                                           'BPix1+BPix2+FPix1_pos+FPix2_pos',
                                           'BPix1+BPix2+FPix1_neg+FPix2_neg',
                                           'BPix1+FPix1_pos+FPix2_pos+FPix3_pos',
                                           'BPix1+FPix1_neg+FPix2_neg+FPix3_neg',
                                           'BPix2+FPix1_pos+FPix2_pos+FPix3_pos',
                                           'BPix2+FPix1_neg+FPix2_neg+FPix3_neg',
                                           'BPix1+BPix2+FPix2_pos+FPix3_pos',
                                           'BPix1+BPix2+FPix2_neg+FPix3_neg',
                                           'BPix1+BPix2+FPix1_pos+FPix3_pos',
                                           'BPix1+BPix2+FPix1_neg+FPix3_neg')
PixelLayerQuadrupletsCA.BPix = cms.PSet(
    TTRHBuilder = cms.string('WithTrackAngle'),
    HitProducer = cms.string('siPixelRecHits'),
)
PixelLayerQuadrupletsCA.FPix = cms.PSet(
    TTRHBuilder = cms.string('WithTrackAngle'),
    HitProducer = cms.string('siPixelRecHits'),
)


