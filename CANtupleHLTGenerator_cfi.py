import FWCore.ParameterSet.Config as cms

# moving to the block.  Will delete the PSet once transition is done
CANtupleHLTGenerator = cms.PSet(
    maxElement = cms.uint32(100000),
    useBending = cms.bool(True),
    useFixedPreFiltering = cms.bool(False),
    ComponentName = cms.string('CANtupleHLTGenerator'),
    extraHitRPhitolerance = cms.double(0.032),
    extraHitRZtolerance = cms.double(0.037),
    useMScat = cms.bool(True),
    phiPreFiltering = cms.double(0.3),
    useBendingCorrection = True,
    SeedComparitorPSet = cms.PSet(
     ComponentName = cms.string('none')
     )
)

import RecoPixelVertexing.PixelLowPtUtilities.LowPtClusterShapeSeedComparitor_cfi
CANtupleHLTGeneratorWithFilter = PixelTripletHLTGenerator.clone()
CANtupleHLTGeneratorWithFilter.SeedComparitorPSet = RecoPixelVertexing.PixelLowPtUtilities.LowPtClusterShapeSeedComparitor_cfi.LowPtClusterShapeSeedComparitor


