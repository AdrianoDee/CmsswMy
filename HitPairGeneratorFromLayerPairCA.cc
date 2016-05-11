#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPairCA.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"

#include "RecoTracker/TkTrackingRegions/interface/HitRZCompatibility.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionBase.h"
#include "RecoTracker/TkHitPairs/interface/OrderedHitPairs.h"
#include "RecoTracker/TkHitPairs/src/InnerDeltaPhi.h"

#include "FWCore/Framework/interface/Event.h"

#define greatz 5E3

using namespace GeomDetEnumerators;
using namespace std;

typedef PixelRecoRange<float> Range;

namespace {
  template<class T> inline T sqr( T t) {return t*t;}
}


#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/FKDPoint.h"


HitPairGeneratorFromLayerPairCA::~HitPairGeneratorFromLayerPairCA() {}

// devirtualizer
#include<tuple>
namespace {

  template<typename Algo>
  struct Kernel {
    using  Base = HitRZCompatibility;
    void set(Base const * a) {
      assert( a->algo()==Algo::me);
      checkRZ=reinterpret_cast<Algo const *>(a);
    }
    
      void operator()(LayerTree& innerTree,const SeedingLayerSetsHits::SeedingLayer& innerLayer,const PixelRecoRange<float>& phiRange,std::vector<unsigned int>& foundHits) const {
          
          constexpr float nSigmaRZ = 3.46410161514f; // std::sqrt(12.f);
          
          const BarrelDetLayer& layerBarrelGeometry = static_cast<const BarrelDetLayer&>(*innerLayer.detLayer());
          
          float vErr = 0.0;
          
          for(auto hit : innerLayer.hits()){
              auto const & gs = static_cast<BaseTrackerRecHit const &>(*hit).globalState();
              auto dv = innerLayer.detLayer()->isBarrel() ? gs.errorZ : gs.errorR;
              auto max = std::max(vErr,dv);
              vErr = max;
          }
          
          vErr *= nSigmaRZ;
          
          float rmax,rmin,zmax,zmin;
          
          auto thickness = innerLayer.detLayer()->surface().bounds().thickness();
          auto u = innerLayer.detLayer()->isBarrel() ? layerBarrelGeometry.specificSurface().radius() : innerLayer.detLayer()->position().z(); //BARREL? Raggio //FWD? z
          Range upperRange = checkRZ->range(u+thickness);
          Range lowerRange = checkRZ->range(u-thickness);
          
          if(innerLayer.detLayer()->isBarrel()){
              
              zmax = std::max(upperRange.max(),lowerRange.max());
              zmin = -std::max(-upperRange.min(),-lowerRange.min());
              rmax = u+thickness+vErr;
              rmin = u-thickness-vErr;
              
          }else{
             
              rmax = std::max(upperRange.max(),lowerRange.max());
              rmin = -std::max(-upperRange.min(),-lowerRange.min());
              zmin = u+thickness+vErr;
              zmax = u-thickness-vErr;
              
          }
          
          LayerPoint minPoint(phiRange.min(),zmin,rmin,1000);
          LayerPoint maxPoint(phiRange.max(),zmax,rmax,1001);
          
          innerTree.LayerTree::search_in_the_box(minPoint,maxPoint,foundHits);
          
      }
    
    Algo const * checkRZ;
    
  };


  template<typename ... Args> using Kernels = std::tuple<Kernel<Args>...>;

}

/*
void HitPairGeneratorFromLayerPair::hitPairs(
					     const TrackingRegion & region, OrderedHitPairs & result,
					     const edm::Event& iEvent, const edm::EventSetup& iSetup, Layers layers) {

  auto const & ds = doublets(region, iEvent, iSetup, layers);
  for (std::size_t i=0; i!=ds.size(); ++i) {
    result.push_back( OrderedHitPair( ds.hit(i,HitDoublets::inner),ds.hit(i,HitDoublets::outer) ));
  }
  if (theMaxElement!=0 && result.size() >= theMaxElement){
     result.clear();
    edm::LogError("TooManyPairs")<<"number of pairs exceed maximum, no pairs produced";
  }
}*/


HitDoubletsCA HitPairGeneratorFromLayerPairCA::doublets (const TrackingRegion& reg,
                                                         const edm::Event & ev,  const edm::EventSetup& es,const SeedingLayerSetsHits::SeedingLayer& innerLayer,
                                                         const SeedingLayerSetsHits::SeedingLayer& outerLayer, LayerTree & innerTree) {

  HitDoubletsCA result(innerLayer,outerLayer);

  InnerDeltaPhi deltaPhi(*outerLayer.detLayer(),*innerLayer.detLayer(), reg, es);

  // std::cout << "layers " << theInnerLayer.detLayer()->seqNum()  << " " << outerLayer.detLayer()->seqNum() << std::endl;

  constexpr float nSigmaPhi = 3.f;
  for (int io = 0; io!=int(outerLayer.hits().size()); ++io) {
      
    Hit const & ohit = outerLayer.hits()[io];
    auto const & gs = static_cast<BaseTrackerRecHit const &>(*ohit).globalState();
    auto loc = gs.position-reg.origin().basicVector();
      
    float oX = gs.position.x();
    float oY = gs.position.y();
    float oZ = gs.position.z();
    float oRv = loc.perp();
      
    float oDrphi = gs.errorRPhi;
    float oDr = gs.errorR;
    float oDz = gs.errorZ;
      
    if (!deltaPhi.prefilter(oX,oY)) continue;
      
    PixelRecoRange<float> phiRange = deltaPhi(oX,oY,oZ,nSigmaPhi*oDrphi);

    if (phiRange.empty()) continue;

    const HitRZCompatibility *checkRZ = reg.checkRZ(innerLayer.detLayer(), ohit, es, outerLayer.detLayer(), oRv, oZ, oDr, oDz);
    if(!checkRZ) continue;

    Kernels<HitZCheck,HitRCheck,HitEtaCheck> kernels;
      
    std::vector<unsigned int> foundHitsInRange;

      switch (checkRZ->algo()) {
          case (HitRZCompatibility::zAlgo) :
              std::get<0>(kernels).set(checkRZ);
              std::get<0>(kernels)(innerTree,innerLayer,phiRange,foundHitsInRange);
              break;
          case (HitRZCompatibility::rAlgo) :
              std::get<1>(kernels).set(checkRZ);
              std::get<0>(kernels)(innerTree,innerLayer,phiRange,foundHitsInRange);
              break;
          case (HitRZCompatibility::etaAlgo) :
              break;
      }
      
      for (auto i=0; i!=(int)foundHitsInRange.size(); ++i) {

          if (theMaxElement!=0 && result.size() >= theMaxElement){
              result.clear();
              edm::LogError("TooManyPairs")<<"number of pairs exceed maximum, no pairs produced";
              delete checkRZ;
              return result;
          }
          result.add(foundHitsInRange[i],io);
      }
      delete checkRZ;
  
  }
  LogDebug("HitPairGeneratorFromLayerPairCA")<<" total number of pairs provided back: "<<result.size();
  result.shrink_to_fit();
  return result;
}
