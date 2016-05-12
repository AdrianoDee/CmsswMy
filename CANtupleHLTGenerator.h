#ifndef CANtupleHLTGenerator_H
#define CANtupleHLTGenerator_H

#include "CombinedHitQuadrupletGeneratorCA.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoUtilities.h"
#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPairCA.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

//#include "RecoPixelVertexing/PixelTriplets/interface/HitQuadrupletGenerator.h"
#include "RecoPixelVertexing/PixelTriplets/interface/OrderedHitSeeds.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoPixelVertexing/PixelTriplets/interface/CANtupleGenerator.h"

#include <utility>
#include <vector>

class SeedComparitor;

class CANtupleHLTGenerator : public CANtupleGenerator {

public:
  CANtupleHLTGenerator( const edm::ParameterSet& cfg, edm::ConsumesCollector& iC);

  virtual ~CANtupleHLTGenerator();

    virtual void getNTuplets(const TrackingRegion& region, OrderedHitSeeds & ntuplets,
                                const edm::Event & ev, const edm::EventSetup& es,
                                SeedingLayerSetsHits::SeedingLayerSet fourLayers
                                ) override;

private:

  std::unique_ptr<SeedComparitor> theComparitor;
    
    class QuantityDependsPtEval {
    public:
        QuantityDependsPtEval(float v1, float v2, float c1, float c2):
        value1_(v1), value2_(v2), curvature1_(c1), curvature2_(c2)
        {}
        
        float value(float curvature) const {
            if(value1_ == value2_) // not enabled
                return value1_;
            
            if(curvature1_ < curvature)
                return value1_;
            if(curvature2_ < curvature && curvature <= curvature1_)
                return value2_ + (curvature-curvature2_)/(curvature1_-curvature2_) * (value1_-value2_);
            return value2_;
        }
        
    private:
        const float value1_;
        const float value2_;
        const float curvature1_;
        const float curvature2_;
    };
    
    // Linear interpolation (in curvature) between value1 at pt1 and
    // value2 at pt2. If disabled, value2 is given (the point is to
    // allow larger/smaller values of the quantity at low pt, so it
    // makes more sense to have the high-pt value as the default).
    class QuantityDependsPt {
    public:
        explicit QuantityDependsPt(const edm::ParameterSet& pset):
        value1_(pset.getParameter<double>("value1")),
        value2_(pset.getParameter<double>("value2")),
        pt1_(pset.getParameter<double>("pt1")),
        pt2_(pset.getParameter<double>("pt2")),
        enabled_(pset.getParameter<bool>("enabled"))
        {
            if(enabled_ && pt1_ >= pt2_)
                throw cms::Exception("Configuration") << "PixelQuadrupletGenerator::QuantityDependsPt: pt1 (" << pt1_ << ") needs to be smaller than pt2 (" << pt2_ << ")";
            if(pt1_ <= 0)
                throw cms::Exception("Configuration") << "PixelQuadrupletGenerator::QuantityDependsPt: pt1 needs to be > 0; is " << pt1_;
            if(pt2_ <= 0)
                throw cms::Exception("Configuration") << "PixelQuadrupletGenerator::QuantityDependsPt: pt2 needs to be > 0; is " << pt2_;
        }
        
        QuantityDependsPtEval evaluator(const edm::EventSetup& es) const {
            if(enabled_) {
                return QuantityDependsPtEval(value1_, value2_,
                                             PixelRecoUtilities::curvature(1.f/pt1_, es),
                                             PixelRecoUtilities::curvature(1.f/pt2_, es));
            }
            return QuantityDependsPtEval(value2_, value2_, 0.f, 0.f);
        }
        
    private:
        const float value1_;
        const float value2_;
        const float pt1_;
        const float pt2_;
        const bool enabled_;
    };

    const float extraHitRZtolerance;
    const float extraHitRPhitolerance;
    const QuantityDependsPt extraPhiTolerance;
    const QuantityDependsPt maxChi2;
    const bool fitFastCircle;
    const bool fitFastCircleChi2Cut;
    const bool useBendingCorrection;
    const bool useFixedPreFiltering;
    const bool useMScat;
    const bool useBend;
    const float dphi;
};
#endif


