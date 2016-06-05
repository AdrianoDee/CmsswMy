#include "RecoPixelVertexing/PixelTriplets/plugins/PixelTripletHLTGenerator.h"
#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPair.h"

#include "ThirdHitPredictionFromInvParabola.h"
#include "ThirdHitRZPrediction.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoUtilities.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "ThirdHitCorrection.h"
#include "RecoTracker/TkHitPairs/interface/RecHitsSortedInPhi.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>

#include "RecoTracker/TkSeedingLayers/interface/SeedComparitorFactory.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedComparitor.h"

#include "DataFormats/GeometryVector/interface/Pi.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/KDTreeLinkerAlgo.h" //amend to point at your copy...
#include "RecoPixelVertexing/PixelTriplets/plugins/KDTreeLinkerTools.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/FKDTree.h"
#include "RecoTracker/TkHitPairs/interface/HitDoubletsCA.h"
#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPairCA.h"

#include "CommonTools/Utils/interface/DynArray.h"

#include "DataFormats/Math/interface/normalizedPhi.h"

#include<cstdio>

#include <fstream>
#include <sstream>
#include <iostream>

using LayerTree = FKDTree<float,3>;

using pixelrecoutilities::LongitudinalBendingCorrection;
using Range=PixelRecoRange<float>;

using namespace std;
using namespace ctfseeding;


PixelTripletHLTGenerator:: PixelTripletHLTGenerator(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC)
  : HitTripletGeneratorFromPairAndLayers(cfg),
    useFixedPreFiltering(cfg.getParameter<bool>("useFixedPreFiltering")),
    extraHitRZtolerance(cfg.getParameter<double>("extraHitRZtolerance")),
    extraHitRPhitolerance(cfg.getParameter<double>("extraHitRPhitolerance")),
    useMScat(cfg.getParameter<bool>("useMultScattering")),
    useBend(cfg.getParameter<bool>("useBending")),
    dphi(useFixedPreFiltering ?  cfg.getParameter<double>("phiPreFiltering") : 0)
{
  edm::ParameterSet comparitorPSet =
    cfg.getParameter<edm::ParameterSet>("SeedComparitorPSet");
  std::string comparitorName = comparitorPSet.getParameter<std::string>("ComponentName");
  if(comparitorName != "none") {
    theComparitor.reset( SeedComparitorFactory::get()->create( comparitorName, comparitorPSet, iC) );
  }
}

PixelTripletHLTGenerator::~PixelTripletHLTGenerator() {}

void PixelTripletHLTGenerator::hitTriplets(const TrackingRegion& region,
					   OrderedHitTriplets & result,
					   const edm::Event & ev,
					   const edm::EventSetup& es,
					   const SeedingLayerSetsHits::SeedingLayerSet& pairLayers,
					   const std::vector<SeedingLayerSetsHits::SeedingLayer>& thirdLayers)
{

  std::cout<<"PixelTripletsHLT : in!"<<std::endl;
  std::ofstream legacy("Txts/legacydoublets.txt",std::ofstream::app);
  /*
  //FeliceKDTree!
  LayerTree alberoFuori;
  alberoFuori.FKDTree<float,3>::make_FKDTreeFromRegionLayer(pairLayers[1],region,ev,es);
  //alberoFuori->FKDTree<float,3>::build();
  std::cout<<"Built?"<<std::endl;
  bool corretto = alberoFuori.FKDTree<float,3>::test_correct_build();
  if(corretto) std::cout<<"Tree Correctly Built"<<std::endl;
  HitPairGeneratorFromLayerPairCA caDoubletsGenerator(0,1,10000);
  */
  if (theComparitor) theComparitor->init(ev, es);


  //std::cout<<"Thickness :  " <<pairLayers[1].detLayer()->surface().bounds().thickness()<<std::endl;


  auto const & doublets = thePairGenerator->doublets(region,ev,es, pairLayers);
  std::size_t n(0);

  const RecHitsSortedInPhi* innerHitsMap = &(*theLayerCache)(pairLayers[0],region,ev,es);
  const RecHitsSortedInPhi* outerHitsMap = &(*theLayerCache)(pairLayers[1],region,ev,es);
  std::vector<Hit> innerHits = innerHitsMap->hits();
  std::vector<Hit> innerRegionHits = region.hits(ev,es,pairLayers[0]);
  std::vector<Hit> outerHits = innerHitsMap->hits();
  std::vector<Hit> outerRegionHits = region.hits(ev,es,pairLayers[1]);


  std::vector<RecHitsSortedInPhi::HitWithPhi> innerRegionHitsPhi;
  for (auto const & hp : innerRegionHits) innerRegionHitsPhi.emplace_back(hp);
  std::vector<int> sortedIndecesInner(innerRegionHits.size());
  std::generate(std::begin(sortedIndecesInner), std::end(sortedIndecesInner), [&]{ return n++; });
  std::sort( std::begin(sortedIndecesInner),std::end(sortedIndecesInner),[&](int i1, int i2) { return innerRegionHitsPhi[i1].phi()<innerRegionHitsPhi[i2].phi(); } );
  std::sort( innerRegionHitsPhi.begin(), innerRegionHitsPhi.end(), RecHitsSortedInPhi::HitLessPhi()); n=0;

  std::vector<RecHitsSortedInPhi::HitWithPhi> outerRegionHitsPhi;
  for (auto const & hp : outerRegionHits) outerRegionHitsPhi.emplace_back(hp);
  std::vector<int> sortedIndecesOuter(outerRegionHits.size());
  std::generate(std::begin(sortedIndecesOuter), std::end(sortedIndecesOuter), [&]{ return n++; });
  std::sort( std::begin(sortedIndecesOuter),std::end(sortedIndecesOuter),[&](int i1, int i2) { return outerRegionHitsPhi[i1].phi()<outerRegionHitsPhi[i2].phi(); } );
  std::sort( outerRegionHitsPhi.begin(), outerRegionHitsPhi.end(), RecHitsSortedInPhi::HitLessPhi()); n=0;


  //std::cout<<"INNER LAYER :  " <<pairLayers[0].name()<<"    "<<"OUTER LAYER :  " <<pairLayers[1].name()<<std::endl;
  //std::cout<<"Legacy Doublets : done!"<<std::endl;
  //std::cout<<doublets.size()<<" doublets found!"<<std::endl;

  legacy<<"==========================["<<pairLayers[0].name()<<" - "<<pairLayers[1].name()<<"]=========================="<<std::endl;
  legacy<<"====== "<<doublets.size()<<" Legacy doublets found!"<<std::endl;
    for(int j=0;j <(int)doublets.size();j++){
        std::cout<<"No phi ------"<<std::endl;
        std::cout<<" [ "<<doublets.innerHitId(j) <<" - "<<doublets.outerHitId(j)<<" ]  ";
        std::cout<<"Inner hit "<<innerHitsMap->x[doublets.innerHitId(j)]<<" - "<<innerHitsMap->y[doublets.innerHitId(j)]<<" - "<<innerHitsMap->z[doublets.innerHitId(j)]<<std::endl;

        Hit const & innerHit = innerRegionHits[doublets.innerHitId(j)].hit();
        Hit const & outerHit = outerRegionHits[doublets.outerHitId(j)].hit();
        auto const & gsInner = innerHit->globalState();
        auto const & gsOuter = outerHit->globalState();
        //auto locInner = gsInner.position-region.origin().basicVector();
        //auto locOuter = gsOuter.position-region.origin().basicVector();

        auto zI = gsInner.position.z(); auto zO = gsOuter.position.z();
        auto xI = gsInner.position.x(); auto xO = gsOuter.position.x();
        auto yI = gsInner.position.y(); auto yO = gsOuter.position.y();
        std::cout<<"Yes phi ------"<<std::endl;
        int innerPhi = sortedIndecesInner[doublets.innerHitId(j)];
        int outerPhi = sortedIndecesOuter[doublets.outerHitId(j)];
        std::cout<<" [ "<<innerPhi <<" - "<<outerPhi<<" ]  ";
        std::cout<<"Inner hit "<<innerHitsMap->x[doublets.innerHitId(j)]<<" - "<<innerHitsMap->y[doublets.innerHitId(j)]<<" - "<<innerHitsMap->z[doublets.innerHitId(j)]<<std::endl;

        Hit const & innerHitPhi = innerRegionHitsPhi[doublets.innerHitId(j)]->hit();
        Hit const & outerHitPhi = outerRegionHitsPhi[doublets.outerHitId(j)]->hit();
        auto const & gsInnerPhi = innerHitPhi->globalState();
        auto const & gsOuterPhi = outerHitPhi->globalState();
        //auto locInner = gsInner.position-region.origin().basicVector();
        //auto locOuter = gsOuter.position-region.origin().basicVector();

        auto zIPhi = gsInnerPhi.position.z(); auto zOPhi = gsOuterPhi.position.z();
        auto xIPhi = gsInnerPhi.position.x(); auto xOPhi = gsOuterPhi.position.x();
        auto yIPhi = gsInnerPhi.position.y(); auto yOPhi = gsOuterPhi.position.y();


          legacy<<" [ "<<doublets.innerHitId(j) <<" - "<<doublets.outerHitId(j)<<" ]  ";
          legacy<<"[ ("<<xI<<" ; "<<yI<<" ; "<<zI<<")"<<"("<<xO<<" ; "<<yO<<" ; "<<zO<<") ]"<<std::endl;
    }

  /*
  auto const & CADoublets = caDoubletsGenerator.doublets(region,ev,es, pairLayers[0],pairLayers[1],alberoFuori);
  std::cout<<"CA Doublets : done!"<<std::endl;
  std::cout<<CADoublets.size()<<" CA doublets found!"<<std::endl;
    for(int j=0;j <(int)CADoublets.size();j++){
        std::cout<<" [ "<<CADoublets.innerHitId(j) <<" - "<<CADoublets.outerHitId(j)<<" ]  ";
    }
  */
  if (doublets.empty()) return;

  auto outSeq =  doublets.detLayer(HitDoublets::outer)->seqNum();


  // std::cout << "pairs " << doublets.size() << std::endl;

  float regOffset = region.origin().perp(); //try to take account of non-centrality (?)
  int size = thirdLayers.size();

  declareDynArray(ThirdHitRZPrediction<PixelRecoLineRZ>, size, preds);
  declareDynArray(ThirdHitCorrection, size, corrections);

  const RecHitsSortedInPhi * thirdHitMap[size];
  typedef RecHitsSortedInPhi::Hit Hit;

  using NodeInfo = KDTreeNodeInfo<unsigned int>;
  std::vector<NodeInfo > layerTree; // re-used throughout
  std::vector<unsigned int> foundNodes; // re-used thoughout
  foundNodes.reserve(100);

  declareDynArray(KDTreeLinkerAlgo<unsigned int>,size, hitTree);
  float rzError[size]; //save maximum errors


  const float maxDelphi = region.ptMin() < 0.3f ? float(M_PI)/4.f : float(M_PI)/8.f; // FIXME move to config??
  const float maxphi = M_PI+maxDelphi, minphi = -maxphi; // increase to cater for any range
  const float safePhi = M_PI-maxDelphi; // sideband


  // fill the prediction vector
  for (int il=0; il<size; ++il) {
    thirdHitMap[il] = &(*theLayerCache)(thirdLayers[il], region, ev, es);
    auto const & hits = *thirdHitMap[il];
    ThirdHitRZPrediction<PixelRecoLineRZ> & pred = preds[il];
    pred.initLayer(thirdLayers[il].detLayer());
    pred.initTolerance(extraHitRZtolerance);

    corrections[il].init(es, region.ptMin(), *doublets.detLayer(HitDoublets::inner), *doublets.detLayer(HitDoublets::outer),
                         *thirdLayers[il].detLayer(), useMScat, useBend);

    layerTree.clear();
    float minv=999999.0f, maxv= -minv; // Initialise to extreme values in case no hits
    float maxErr=0.0f;
    for (unsigned int i=0; i!=hits.size(); ++i) {
      auto angle = hits.phi(i);
      auto v =  hits.gv(i);
      //use (phi,r) for endcaps rather than (phi,z)
      minv = std::min(minv,v);  maxv = std::max(maxv,v);
      float myerr = hits.dv[i];
      maxErr = std::max(maxErr,myerr);
      layerTree.emplace_back(i, angle, v); // save it
      // populate side-bands
      if (angle>safePhi) layerTree.emplace_back(i, angle-Geom::ftwoPi(), v);
      else if (angle<-safePhi) layerTree.emplace_back(i, angle+Geom::ftwoPi(), v);
    }
    KDTreeBox phiZ(minphi, maxphi, minv-0.01f, maxv+0.01f);  // declare our bounds
    //add fudge factors in case only one hit and also for floating-point inaccuracy
    hitTree[il].build(layerTree, phiZ); // make KDtree
    rzError[il] = maxErr; //save error
    // std::cout << "layer " << thirdLayers[il].detLayer()->seqNum() << " " << layerTree.size() << std::endl;
  }

  float imppar = region.originRBound();
  float imppartmp = region.originRBound()+region.origin().perp();
  float curv = PixelRecoUtilities::curvature(1.f/region.ptMin(), es);

  for (std::size_t ip =0;  ip!=doublets.size(); ip++) {
    auto xi = doublets.x(ip,HitDoublets::inner);
    auto yi = doublets.y(ip,HitDoublets::inner);
    auto zi = doublets.z(ip,HitDoublets::inner);
    auto rvi = doublets.rv(ip,HitDoublets::inner);
    auto xo = doublets.x(ip,HitDoublets::outer);
    auto yo = doublets.y(ip,HitDoublets::outer);
    auto zo = doublets.z(ip,HitDoublets::outer);
    auto rvo = doublets.rv(ip,HitDoublets::outer);

    auto toPos = std::signbit(zo-zi);

    PixelRecoPointRZ point1(rvi, zi);
    PixelRecoPointRZ point2(rvo, zo);
    PixelRecoLineRZ  line(point1, point2);
    ThirdHitPredictionFromInvParabola predictionRPhi(xi-region.origin().x(),yi-region.origin().y(),
						     xo-region.origin().x(),yo-region.origin().y(),
						     imppar,curv,extraHitRPhitolerance);

    ThirdHitPredictionFromInvParabola predictionRPhitmp(xi,yi,xo,yo,imppartmp,curv,extraHitRPhitolerance);

    // printf("++Constr %f %f %f %f %f %f %f\n",xi,yi,xo,yo,imppartmp,curv,extraHitRPhitolerance);

    // std::cout << ip << ": " << point1.r() << ","<< point1.z() << " "
    //                        << point2.r() << ","<< point2.z() <<std::endl;

    for (int il=0; il!=size; ++il) {
      const DetLayer * layer = thirdLayers[il].detLayer();
      auto barrelLayer = layer->isBarrel();

      if ( (!barrelLayer) & (toPos != std::signbit(layer->position().z())) ) continue;

      if (hitTree[il].empty()) continue; // Don't bother if no hits

      auto const & hits = *thirdHitMap[il];

      auto & correction = corrections[il];

      correction.init(line, point2, outSeq);

      auto & predictionRZ =  preds[il];

      predictionRZ.initPropagator(&line);
      Range rzRange = predictionRZ();
      correction.correctRZRange(rzRange);

      Range phiRange;
      if (useFixedPreFiltering) {
	float phi0 = doublets.phi(ip,HitDoublets::outer);
	phiRange = Range(phi0-dphi,phi0+dphi);
      }
      else {
	Range radius;
	if (barrelLayer) {
	  radius =  predictionRZ.detRange();
	} else {
	  radius = Range(max(rzRange.min(), predictionRZ.detSize().min()),
			 min(rzRange.max(), predictionRZ.detSize().max()) );
	}
	if (radius.empty()) continue;

	// std::cout << "++R " << radius.min() << " " << radius.max()  << std::endl;

	auto rPhi1 = predictionRPhitmp(radius.max());
        bool ok1 = !rPhi1.empty();
        if (ok1) {
          correction.correctRPhiRange(rPhi1);
          rPhi1.first  /= radius.max();
          rPhi1.second /= radius.max();
        }
	auto rPhi2 = predictionRPhitmp(radius.min());
        bool ok2 = !rPhi2.empty();
       	if (ok2) {
	  correction.correctRPhiRange(rPhi2);
	  rPhi2.first  /= radius.min();
	  rPhi2.second /= radius.min();
        }

        if (ok1) {
          rPhi1.first = normalizedPhi(rPhi1.first);
          rPhi1.second = proxim(rPhi1.second,rPhi1.first);
          if(ok2) {
            rPhi2.first = proxim(rPhi2.first,rPhi1.first);
            rPhi2.second = proxim(rPhi2.second,rPhi1.first);
            phiRange = rPhi1.sum(rPhi2);
          } else phiRange=rPhi1;
        } else if(ok2) {
          rPhi2.first = normalizedPhi(rPhi2.first);
          rPhi2.second = proxim(rPhi2.second,rPhi2.first);
          phiRange=rPhi2;
        } else continue;

      }

      constexpr float nSigmaRZ = 3.46410161514f; // std::sqrt(12.f); // ...and continue as before
      constexpr float nSigmaPhi = 3.f;

      foundNodes.clear(); // Now recover hits in bounding box...
      float prmin=phiRange.min(), prmax=phiRange.max();


      if (prmax-prmin>maxDelphi) {
        auto prm = phiRange.mean();
        prmin = prm - 0.5f*maxDelphi;
        prmax = prm + 0.5f*maxDelphi;
      }

      if (barrelLayer)
	{
	  Range regMax = predictionRZ.detRange();
	  Range regMin = predictionRZ(regMax.min()-regOffset);
	  regMax = predictionRZ(regMax.max()+regOffset);
	  correction.correctRZRange(regMin);
	  correction.correctRZRange(regMax);
	  if (regMax.min() < regMin.min()) { swap(regMax, regMin);}
	  KDTreeBox phiZ(prmin, prmax, regMin.min()-nSigmaRZ*rzError[il], regMax.max()+nSigmaRZ*rzError[il]);
	  hitTree[il].search(phiZ, foundNodes);
	}
      else
	{
	  KDTreeBox phiZ(prmin, prmax,
			 rzRange.min()-regOffset-nSigmaRZ*rzError[il],
			 rzRange.max()+regOffset+nSigmaRZ*rzError[il]);
	  hitTree[il].search(phiZ, foundNodes);
	}

      // std::cout << ip << ": " << thirdLayers[il].detLayer()->seqNum() << " " << foundNodes.size() << " " << prmin << " " << prmax << std::endl;


      // int kk=0;
      for (auto KDdata : foundNodes) {

	if (theMaxElement!=0 && result.size() >= theMaxElement){
	  result.clear();
	  edm::LogError("TooManyTriplets")<<" number of triples exceeds maximum. no triplets produced.";
	  return;
	}

	float p3_u = hits.u[KDdata];
	float p3_v =  hits.v[KDdata];
	float p3_phi =  hits.lphi[KDdata];

       //if ((kk++)%100==0)
       //std::cout << kk << ": " << p3_u << " " << p3_v << " " << p3_phi << std::endl;


	Range allowed = predictionRZ(p3_u);
	correction.correctRZRange(allowed);
	float vErr = nSigmaRZ *hits.dv[KDdata];
	Range hitRange(p3_v-vErr, p3_v+vErr);
	Range crossingRange = allowed.intersection(hitRange);
	if (crossingRange.empty())  continue;

	float ir = 1.f/hits.rv(KDdata);
        // limit error to 90 degree
        constexpr float maxPhiErr = 0.5*M_PI;
	float phiErr = nSigmaPhi * hits.drphi[KDdata]*ir;
        phiErr = std::min(maxPhiErr, phiErr);
        bool nook=true;
	for (int icharge=-1; icharge <=1; icharge+=2) {
	  Range rangeRPhi = predictionRPhi(hits.rv(KDdata), icharge);
          if(rangeRPhi.first>rangeRPhi.second) continue; // range is empty
	  correction.correctRPhiRange(rangeRPhi);
	  if (checkPhiInRange(p3_phi, rangeRPhi.first*ir-phiErr, rangeRPhi.second*ir+phiErr,maxPhiErr)) {
	    // insert here check with comparitor
	    OrderedHitTriplet hittriplet( doublets.hit(ip,HitDoublets::inner), doublets.hit(ip,HitDoublets::outer), hits.theHits[KDdata].hit());
	    if (!theComparitor  || theComparitor->compatible(hittriplet,region) ) {
	      result.push_back( hittriplet );
	    } else {
	      LogDebug("RejectedTriplet") << "rejected triplet from comparitor ";
	    }
	    nook=false; break;
	  }
	}
        if (nook) LogDebug("RejectedTriplet") << "rejected triplet from second phicheck " << p3_phi;
      }
    }
  }
  // std::cout << "triplets " << result.size() << std::endl;
}
