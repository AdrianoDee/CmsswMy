#include "PixelQuadrupletGenerator.h"
#include "ThirdHitRZPrediction.h"
#include "RecoPixelVertexing/PixelTriplets/interface/ThirdHitPredictionFromCircle.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoLineRZ.h"

#include "RecoPixelVertexing/PixelTriplets/plugins/KDTreeLinkerAlgo.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/KDTreeLinkerTools.h"
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"

#include "RecoPixelVertexing/PixelTrackFitting/src/RZLine.h"
#include "RecoTracker/TkSeedGenerator/interface/FastCircleFit.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoUtilities.h"
#include "RecoTracker/TkMSParametrization/interface/LongitudinalBendingCorrection.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "RecoTracker/TkSeedingLayers/interface/SeedComparitorFactory.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedComparitor.h"

#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPairCA.h"
#include "RecoTracker/TkHitPairs/interface/HitDoubletsCA.h"

#include "CommonTools/Utils/interface/DynArray.h"

#include "FWCore/Utilities/interface/isFinite.h"

#include <fstream>
#include <sstream>
#include <iostream>

namespace {
    template <typename T>
    T sqr(T x) {
        return x*x;
    }
}

typedef PixelRecoRange<float> Range;

PixelQuadrupletGenerator::PixelQuadrupletGenerator(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC):
extraHitRZtolerance(cfg.getParameter<double>("extraHitRZtolerance")), //extra window in ThirdHitRZPrediction range
extraHitRPhitolerance(cfg.getParameter<double>("extraHitRPhitolerance")), //extra window in ThirdHitPredictionFromCircle range (divide by R to get phi)
extraPhiTolerance(cfg.getParameter<edm::ParameterSet>("extraPhiTolerance")),
maxChi2(cfg.getParameter<edm::ParameterSet>("maxChi2")),
fitFastCircle(cfg.getParameter<bool>("fitFastCircle")),
fitFastCircleChi2Cut(cfg.getParameter<bool>("fitFastCircleChi2Cut")),
useBendingCorrection(cfg.getParameter<bool>("useBendingCorrection"))
{
    if(cfg.exists("SeedComparitorPSet")) {
        edm::ParameterSet comparitorPSet =
        cfg.getParameter<edm::ParameterSet>("SeedComparitorPSet");
        std::string comparitorName = comparitorPSet.getParameter<std::string>("ComponentName");
        if(comparitorName != "none") {
            theComparitor.reset(SeedComparitorFactory::get()->create(comparitorName, comparitorPSet, iC));
        }
    }
}

PixelQuadrupletGenerator::~PixelQuadrupletGenerator() {}


void PixelQuadrupletGenerator::hitQuadruplets(const TrackingRegion& region, OrderedHitSeeds& result,
                                              const edm::Event& ev, const edm::EventSetup& es,
                                              const SeedingLayerSetsHits::SeedingLayerSet& tripletLayers,
                                              const std::vector<SeedingLayerSetsHits::SeedingLayer>& fourthLayers)
{
    if (theComparitor) theComparitor->init(ev, es);

    OrderedHitTriplets triplets;
    theTripletGenerator->hitTriplets(region, triplets, ev, es,
                                     tripletLayers, // pair generator picks the correct two layers from these
                                     std::vector<SeedingLayerSetsHits::SeedingLayer>{tripletLayers[2]});
    if(triplets.empty()) return;

    const size_t size = fourthLayers.size();

    const RecHitsSortedInPhi *fourthHitMap[size];
    typedef RecHitsSortedInPhi::Hit Hit;

    using NodeInfo = KDTreeNodeInfo<unsigned int>;
    std::vector<NodeInfo > layerTree; // re-used throughout
    std::vector<unsigned int> foundNodes; // re-used thoughout

    declareDynArray(KDTreeLinkerAlgo<unsigned int>, size, hitTree);
    float rzError[size]; //save maximum errors

    declareDynArray(ThirdHitRZPrediction<PixelRecoLineRZ>, size, preds);

    // Build KDtrees
    for(size_t il=0; il!=size; ++il) {
        fourthHitMap[il] = &(*theLayerCache)(fourthLayers[il], region, ev, es);
        auto const& hits = *fourthHitMap[il];

        ThirdHitRZPrediction<PixelRecoLineRZ> & pred = preds[il];
        pred.initLayer(fourthLayers[il].detLayer());
        pred.initTolerance(extraHitRZtolerance);

        layerTree.clear();
        float maxphi = Geom::ftwoPi(), minphi = -maxphi; // increase to cater for any range
        float minv=999999.0, maxv= -999999.0; // Initialise to extreme values in case no hits
        float maxErr=0.0f;
        for (unsigned int i=0; i!=hits.size(); ++i) {
            auto angle = hits.phi(i);
            auto v =  hits.gv(i);
            //use (phi,r) for endcaps rather than (phi,z)
            minv = std::min(minv,v);  maxv = std::max(maxv,v);
            float myerr = hits.dv[i];
            maxErr = std::max(maxErr,myerr);
            layerTree.emplace_back(i, angle, v); // save it
            if (angle < 0)  // wrap all points in phi
            { layerTree.emplace_back(i, angle+Geom::ftwoPi(), v);}
            else
            { layerTree.emplace_back(i, angle-Geom::ftwoPi(), v);}
        }
        KDTreeBox phiZ(minphi, maxphi, minv-0.01f, maxv+0.01f);  // declare our bounds
        //add fudge factors in case only one hit and also for floating-point inaccuracy
        hitTree[il].build(layerTree, phiZ); // make KDtree
        rzError[il] = maxErr; // save error
    }

    const QuantityDependsPtEval maxChi2Eval = maxChi2.evaluator(es);
    const QuantityDependsPtEval extraPhiToleranceEval = extraPhiTolerance.evaluator(es);

    // re-used thoughout, need to be vectors because of RZLine interface
    std::vector<float> bc_r(4), bc_z(4), bc_errZ(4);

    // Loop over triplets
    for(const auto& triplet: triplets) {
        GlobalPoint gp0 = triplet.inner()->globalPosition();
        GlobalPoint gp1 = triplet.middle()->globalPosition();
        GlobalPoint gp2 = triplet.outer()->globalPosition();

        PixelRecoLineRZ line(gp0, gp2);
        ThirdHitPredictionFromCircle predictionRPhi(gp0, gp2, extraHitRPhitolerance);

        const double curvature = predictionRPhi.curvature(ThirdHitPredictionFromCircle::Vector2D(gp1.x(), gp1.y()));

        const float abscurv = std::abs(curvature);
        const float thisMaxChi2 = maxChi2Eval.value(abscurv);
        const float thisExtraPhiTolerance = extraPhiToleranceEval.value(abscurv);

        constexpr float nSigmaRZ = 3.46410161514f; // std::sqrt(12.f); // ...and continue as before

        auto isBarrel = [](const unsigned id) -> bool {
            return id == PixelSubdetector::PixelBarrel;
        };

        declareDynArray(GlobalPoint,4, gps);
        declareDynArray(GlobalError,4, ges);
        declareDynArray(bool,4, barrels);
        gps[0] = triplet.inner()->globalPosition();
        ges[0] = triplet.inner()->globalPositionError();
        barrels[0] = isBarrel(triplet.inner()->geographicalId().subdetId());

        gps[1] = triplet.middle()->globalPosition();
        ges[1] = triplet.middle()->globalPositionError();
        barrels[1] = isBarrel(triplet.middle()->geographicalId().subdetId());

        gps[2] = triplet.outer()->globalPosition();
        ges[2] = triplet.outer()->globalPositionError();
        barrels[2] = isBarrel(triplet.outer()->geographicalId().subdetId());

        for(size_t il=0; il!=size; ++il) {
            if(hitTree[il].empty()) continue; // Don't bother if no hits

            auto const& hits = *fourthHitMap[il];
            const DetLayer *layer = fourthLayers[il].detLayer();
            bool barrelLayer = layer->isBarrel();

            auto& predictionRZ = preds[il];
            predictionRZ.initPropagator(&line);
            Range rzRange = predictionRZ(); // z in barrel, r in endcap

            // construct search box and search
            Range phiRange;
            if(barrelLayer) {
                auto radius = static_cast<const BarrelDetLayer *>(layer)->specificSurface().radius();
                double phi = predictionRPhi.phi(curvature, radius);
                phiRange = Range(phi-thisExtraPhiTolerance, phi+thisExtraPhiTolerance);
            }
            else {
                double phi1 = predictionRPhi.phi(curvature, rzRange.min());
                double phi2 = predictionRPhi.phi(curvature, rzRange.max());
                phiRange = Range(std::min(phi1, phi2)-thisExtraPhiTolerance, std::max(phi1, phi2)+thisExtraPhiTolerance);
            }

            KDTreeBox phiZ(phiRange.min(), phiRange.max(),
                           rzRange.min()-nSigmaRZ*rzError[il],
                           rzRange.max()+nSigmaRZ*rzError[il]);

            foundNodes.clear();
            hitTree[il].search(phiZ, foundNodes);

            if(foundNodes.empty()) {
                continue;
            }

            SeedingHitSet::ConstRecHitPointer selectedHit = nullptr;
            float selectedChi2 = std::numeric_limits<float>::max();
            for(auto hitIndex: foundNodes) {
                const auto& hit = hits.theHits[hitIndex].hit();

                // Reject comparitor. For now, because of technical
                // limitations, pass three hits to the comparitor
                // TODO: switch to using hits from 2-3-4 instead of 1-3-4?
                // Eventually we should fix LowPtClusterShapeSeedComparitor to
                // accept quadruplets.
                if(theComparitor) {
                    SeedingHitSet tmpTriplet(triplet.inner(), triplet.outer(), hit);
                    if(!theComparitor->compatible(tmpTriplet, region)) {
                        continue;
                    }
                }

                gps[3] = hit->globalPosition();
                ges[3] = hit->globalPositionError();
                barrels[3] = isBarrel(hit->geographicalId().subdetId());

                float chi2 = std::numeric_limits<float>::quiet_NaN();
                // TODO: Do we have any use case to not use bending correction?
                if(useBendingCorrection) {
                    // Following PixelFitterByConformalMappingAndLine
                    const float simpleCot = ( gps.back().z()-gps.front().z() )/ (gps.back().perp() - gps.front().perp() );
                    const float pt = 1/PixelRecoUtilities::inversePt(abscurv, es);
                    for (int i=0; i< 4; ++i) {
                        const GlobalPoint & point = gps[i];
                        const GlobalError & error = ges[i];
                        bc_r[i] = sqrt( sqr(point.x()-region.origin().x()) + sqr(point.y()-region.origin().y()) );
                        bc_r[i] += pixelrecoutilities::LongitudinalBendingCorrection(pt,es)(bc_r[i]);
                        bc_z[i] = point.z()-region.origin().z();
                        bc_errZ[i] =  (barrels[i]) ? sqrt(error.czz()) : sqrt( error.rerr(point) )*simpleCot;
                    }
                    RZLine rzLine(bc_r,bc_z,bc_errZ);
                    float      cottheta, intercept, covss, covii, covsi;
                    rzLine.fit(cottheta, intercept, covss, covii, covsi);
                    chi2 = rzLine.chi2(cottheta, intercept);
                }
                else {
                    RZLine rzLine(gps, ges, barrels);
                    float  cottheta, intercept, covss, covii, covsi;
                    rzLine.fit(cottheta, intercept, covss, covii, covsi);
                    chi2 = rzLine.chi2(cottheta, intercept);
                }
                if(edm::isNotFinite(chi2) || chi2 > thisMaxChi2) {
                    continue;
                }
                // TODO: Do we have any use case to not use circle fit? Maybe
                // HLT where low-pT inefficiency is not a problem?
                if(fitFastCircle) {
                    FastCircleFit c(gps, ges);
                    chi2 += c.chi2();
                    if(edm::isNotFinite(chi2))
                    continue;
                    if(fitFastCircleChi2Cut && chi2 > thisMaxChi2)
                    continue;
                }


                if(chi2 < selectedChi2) {
                    selectedChi2 = chi2;
                    selectedHit = hit;
                }
            }
            if(selectedHit)
            result.emplace_back(triplet.inner(), triplet.middle(), triplet.outer(), selectedHit);
        }
    }
}


void PixelQuadrupletGenerator::hitQuadruplets( const TrackingRegion& region, OrderedHitSeeds& result,
                                              const edm::Event& ev, const edm::EventSetup& es,
                                              const SeedingLayerSetsHits::SeedingLayerSet& fourLayers)
{
    std::cout<<std::endl<<"PixelQuadruplets CA : in!"<<std::endl;
    std::ofstream cadoublets("Txts/cadoublets.txt",std::ofstream::app);
    if (theComparitor) theComparitor->init(ev, es);

    HitPairGeneratorFromLayerPairCA caDoubletsGenerator(0,1,10000);

    //std::vector<FKDTree<float,3>*> layersHitsTree;

    /*
    bool treesFlags[3];
    for (int j=0; j<3; j++) {
        std::cout<<"Checking cache for id : "<<fourLayers[j].index()<<std::endl;
        if(theKDTreeCache->checkCache(fourLayers[j].index())){
            std::cout<<"Tree already done "<<std::endl;
            theKDTreeCache->getCache(fourLayers[j].index(),&trees[j]);
            std::cout<<"Tree copied "<<std::endl;
            //layersHitsTree.push_back(&trees[j]);
        } else {
            std::cout<<"Tree to be created "<<fourLayers[j].index()<<std::endl;
            trees[j].LayerTree::make_FKDTreeFromRegionLayer(fourLayers[j],region,ev,es);
            std::cout<<"Tree created "<<std::endl;
            theKDTreeCache->writeCache(fourLayers[j].index(),&trees[j]);
            std::cout<<"Tree cached "<<std::endl;
            //layersHitsTree.push_back(&trees[j]);
        }
    }*/

    LayerTree & innerTree = theKDTreeCache->getTree(fourLayers[0],region,ev,es);
    LayerTree & middleTree = theKDTreeCache->getTree(fourLayers[1],region,ev,es);
    LayerTree & outerTree = theKDTreeCache->getTree(fourLayers[2],region,ev,es);

    //for (int j=0; j<4; j++) {
    //    if(fourLayers[j]) std::cout<<fourLayers[j].hits()[61]->globalState().r;
    //}

    //TESTING
    std::vector<unsigned int> idS = innerTree.getIdVector();

    if(innerTree.empty()) std::cout<<"Tree Empty"<<std::endl;
        else{
            for (int j=0; j<(int)idS.size(); j++) {
                std::cout<<" "<<idS[j]<<" ";
            }
        }
    /*
    LayerTree treeFourth; treeFourth.FKDTree<float,3>::make_FKDTreeFromRegionLayer(fourLayers[3],region,ev,es);
    layersHitsTree.push_back(&treeFourth);*/


    /*
    auto const & albero1 = (theKDTreeCache)(fourLayers[0],region,ev,es); layersHitsTree.push_back(&albero1);
    auto const & albero2 = (theKDTreeCache)(fourLayers[1],region,ev,es); layersHitsTree.push_back(&albero2);
    auto const & albero3 = (theKDTreeCache)(fourLayers[2],region,ev,es); layersHitsTree.push_back(&albero3);
    auto const & albero4 = (theKDTreeCache)(fourLayers[3],region,ev,es); layersHitsTree.push_back(&albero4);
    */

    //if (&innerTree == nullptr) std::cout<<"CAZZOOOOOOO!"<<std::endl;

    std::size_t n(0);
    std::vector<Hit> hits0 = region.hits(ev,es,fourLayers[0]);
    std::vector<RecHitsSortedInPhi::HitWithPhi> hit0phi;
    for (auto const & hp : hits0) hit0phi.emplace_back(hp);
    std::vector<int> sortedIndeces0(hits0.size());
    std::generate(std::begin(sortedIndeces0), std::end(sortedIndeces0), [&]{ return n++; });
    std::sort( std::begin(sortedIndeces0),std::end(sortedIndeces0),[&](int i1, int i2) { return hit0phi[i1].phi()<hit0phi[i2].phi(); } );
    std::sort( hit0phi.begin(), hit0phi.end(), RecHitsSortedInPhi::HitLessPhi());
    for(int j=0;j <(int)n;j++){
      Hit const & hit = hits0[j]->hit();
      Hit const & phiHit = hit0phi[j].hit();
      auto const & gsHit = hit->globalState();
      auto const & gsPhi = phiHit->globalState();
      //auto locInner = gsInner.position-region.origin().basicVector();
      //auto locOuter = gsOuter.position-region.origin().basicVector();

      auto zI = gsHit.position.z(); auto zO = gsPhi.position.z();
      auto xI = gsHit.position.x(); auto xO = gsPhi.position.x();
      auto yI = gsHit.position.y(); auto yO = gsPhi.position.y();

      std::cout<<"[ NoPhi : ("<<xI<<" ; "<<yI<<" ; "<<zI<<")"<<" Phi : ("<<xO<<" ; "<<yO<<" ; "<<zO<<") ]"<<std::endl;

    }

    std::vector<Hit> hits1 = region.hits(ev,es,fourLayers[1]);
    std::vector<RecHitsSortedInPhi::HitWithPhi> hit1phi;
    for (auto const & hp : hits1) hit1phi.emplace_back(hp);
    std::vector<int> sortedIndeces1(hits0.size());
    std::vector<Hit> hits2 = region.hits(ev,es,fourLayers[2]);
    std::vector<RecHitsSortedInPhi::HitWithPhi> hit2phi;
    for (auto const & hp : hits2) hit2phi.emplace_back(hp);
    std::vector<int> sortedIndeces2(hits0.size());
    std::vector<Hit> hits3 = region.hits(ev,es,fourLayers[3]);
    std::vector<RecHitsSortedInPhi::HitWithPhi> hit3phi;
    for (auto const & hp : hits3) hit3phi.emplace_back(hp);
    std::vector<int> sortedIndeces3(hits0.size());

    auto const & doublets1 = caDoubletsGenerator.doublets(region,ev,es,fourLayers[0],fourLayers[1],&innerTree);
    //std::cout<<"INNER LAYER :  " <<fourLayers[0].name()<<"    "<<"OUTER LAYER :  " <<fourLayers[1].name()<<std::endl;
    //std::cout<<"CA Doublets : done!"<<std::endl;
    //std::cout<<doublets1.size()<<" CA doublets found!"<<std::endl;

    cadoublets<<"==========================["<<fourLayers[0].name()<<" - "<<fourLayers[1].name()<<"]=========================="<<std::endl;
    cadoublets<<"====== "<<doublets1.size()<<" CA doublets found!"<<std::endl;

    for(int j=0;j <(int)doublets1.size();j++){
      Hit const & innerHit = hits0[doublets1.innerHitId(j)]->hit();
      Hit const & outerHit = hits1[doublets1.outerHitId(j)]->hit();
      auto const & gsInner = innerHit->globalState();
      auto const & gsOuter = outerHit->globalState();
      //auto locInner = gsInner.position-region.origin().basicVector();
      //auto locOuter = gsOuter.position-region.origin().basicVector();

      auto zI = gsInner.position.z(); auto zO = gsOuter.position.z();
      auto xI = gsInner.position.x(); auto xO = gsOuter.position.x();
      auto yI = gsInner.position.y(); auto yO = gsOuter.position.y();

        cadoublets<<" [ "<<doublets1.innerHitId(j) <<" - "<<doublets1.outerHitId(j)<<" ]  ";
        cadoublets<<"[ ("<<xI<<" ; "<<yI<<" ; "<<zI<<")"<<"("<<xO<<" ; "<<yO<<" ; "<<zO<<") ]"<<std::endl;
    }



    auto const & doublets2 = caDoubletsGenerator.doublets(region,ev,es,fourLayers[1],fourLayers[2],&middleTree);
    std::cout<<"INNER LAYER :  " <<fourLayers[1].name()<<"    "<<"OUTER LAYER :  " <<fourLayers[2].name()<<std::endl;
    std::cout<<"CA Doublets : done!"<<std::endl;
    std::cout<<doublets2.size()<<" CA doublets found!"<<std::endl;
    for(int j=0;j <(int)doublets2.size();j++){
        std::cout<<" [ "<<doublets2.innerHitId(j) <<" - "<<doublets2.outerHitId(j)<<" ]  ";
    }

    auto const & doublets3 = caDoubletsGenerator.doublets(region,ev,es,fourLayers[2],fourLayers[3],&outerTree);
    std::cout<<"INNER LAYER :  " <<fourLayers[2].name()<<"    "<<"OUTER LAYER :  " <<fourLayers[3].name()<<std::endl;
    std::cout<<"CA Doublets : done!"<<std::endl;
    std::cout<<doublets3.size()<<" CA doublets found!"<<std::endl;
    for(int j=0;j <(int)doublets3.size();j++){
        std::cout<<" [ "<<doublets3.innerHitId(j) <<" - "<<doublets3.outerHitId(j)<<" ]  ";
    }

    cadoublets<<"==========================["<<fourLayers[1].name()<<" - "<<fourLayers[2].name()<<"]=========================="<<std::endl;
    cadoublets<<"====== "<<doublets2.size()<<" CA doublets found!"<<std::endl;

    for(int j=0;j <(int)doublets2.size();j++){
      Hit const & innerHit = hits1[doublets2.innerHitId(j)]->hit();
      Hit const & outerHit = hits2[doublets2.outerHitId(j)]->hit();
      auto const & gsInner = innerHit->globalState();
      auto const & gsOuter = outerHit->globalState();
      //auto locInner = gsInner.position-region.origin().basicVector();
      //auto locOuter = gsOuter.position-region.origin().basicVector();

      auto zI = gsInner.position.z(); auto zO = gsOuter.position.z();
      auto xI = gsInner.position.x(); auto xO = gsOuter.position.x();
      auto yI = gsInner.position.y(); auto yO = gsOuter.position.y();

        cadoublets<<" [ "<<doublets1.innerHitId(j) <<" - "<<doublets1.outerHitId(j)<<" ]  ";
        cadoublets<<"[ ("<<xI<<" ; "<<yI<<" ; "<<zI<<")"<<"("<<xO<<" ; "<<yO<<" ; "<<zO<<") ]"<<std::endl;
    }

    cadoublets<<"==========================["<<fourLayers[2].name()<<" - "<<fourLayers[3].name()<<"]=========================="<<std::endl;
    cadoublets<<"====== "<<doublets3.size()<<" CA doublets found!"<<std::endl;

    for(int j=0;j <(int)doublets3.size();j++){
      Hit const & innerHit = hits2[doublets3.innerHitId(j)]->hit();
      Hit const & outerHit = hits3[doublets3.outerHitId(j)]->hit();
      auto const & gsInner = innerHit->globalState();
      auto const & gsOuter = outerHit->globalState();
      //auto locInner = gsInner.position-region.origin().basicVector();
      //auto locOuter = gsOuter.position-region.origin().basicVector();

      auto zI = gsInner.position.z(); auto zO = gsOuter.position.z();
      auto xI = gsInner.position.x(); auto xO = gsOuter.position.x();
      auto yI = gsInner.position.y(); auto yO = gsOuter.position.y();

        cadoublets<<" [ "<<doublets1.innerHitId(j) <<" - "<<doublets1.outerHitId(j)<<" ]  ";
        cadoublets<<"[ ("<<xI<<" ; "<<yI<<" ; "<<zI<<")"<<"("<<xO<<" ; "<<yO<<" ; "<<zO<<") ]"<<std::endl;
    }

    //std::cout<<"ZETA QUI -----"<<std::endl;
    //std::cout<<doublets1.z(0,HitDoubletsCA::inner)<<std::endl;
    //std::cout<<doublets1.z(0,HitDoubletsCA::inner)<<std::endl;
    //std::cout<<doublets3.z(0,HitDoubletsCA::inner)<<std::endl;

    //std::vector<HitDoubletsCA> layersDoublets;

    //std::vector<CACell::CAntuplet> foundQuadruplets;
    std::vector<unsigned int> indexOfFirstCellOfLayer;
    std::vector<unsigned int> numberOfCellsPerLayer;

    /* CON LE CACHEs : NON TESTATO
    for (auto layer : fourLayers)
    {
        layersHitsTree.push_back(&(*theKDTreeCache)(layer,region,ev,es));
    }

    for (int j=0;j<(int)layersHitsTree.size()-1;j++)
    {
        layersDoublets.push_back(&(*theDoubletsCache)(fourLayers[j],fourLayers[j+1],(*layersHitsTree[j]),region,ev,es)); //Passa vector
    }
    */

    //QUI SOTTO CI VA IL CELLULAR AUTOMATON
/*
    OrderedHitTriplets triplets;
    theTripletGenerator->hitTriplets(region, triplets, ev, es,
                                     tripletLayers, // pair generator picks the correct two layers from these
                                     std::vector<SeedingLayerSetsHits::SeedingLayer>{tripletLayers[2]});
    if(triplets.empty()) return;

    const size_t size = fourthLayers.size();

    const RecHitsSortedInPhi *fourthHitMap[size];
    typedef RecHitsSortedInPhi::Hit Hit;

    using NodeInfo = KDTreeNodeInfo<unsigned int>;
    std::vector<NodeInfo > layerTree; // re-used throughout
    std::vector<unsigned int> foundNodes; // re-used thoughout

    declareDynArray(KDTreeLinkerAlgo<unsigned int>, size, hitTree);
    float rzError[size]; //save maximum errors

    declareDynArray(ThirdHitRZPrediction<PixelRecoLineRZ>, size, preds);

    // Build KDtrees
    for(size_t il=0; il!=size; ++il) {
        fourthHitMap[il] = &(*theLayerCache)(fourthLayers[il], region, ev, es);
        auto const& hits = *fourthHitMap[il];

        ThirdHitRZPrediction<PixelRecoLineRZ> & pred = preds[il];
        pred.initLayer(fourthLayers[il].detLayer());
        pred.initTolerance(extraHitRZtolerance);

        layerTree.clear();
        float maxphi = Geom::ftwoPi(), minphi = -maxphi; // increase to cater for any range
        float minv=999999.0, maxv= -999999.0; // Initialise to extreme values in case no hits
        float maxErr=0.0f;
        for (unsigned int i=0; i!=hits.size(); ++i) {
            auto angle = hits.phi(i);
            auto v =  hits.gv(i);
            //use (phi,r) for endcaps rather than (phi,z)
            minv = std::min(minv,v);  maxv = std::max(maxv,v);
            float myerr = hits.dv[i];
            maxErr = std::max(maxErr,myerr);
            layerTree.emplace_back(i, angle, v); // save it
            if (angle < 0)  // wrap all points in phi
            { layerTree.emplace_back(i, angle+Geom::ftwoPi(), v);}
            else
            { layerTree.emplace_back(i, angle-Geom::ftwoPi(), v);}
        }
        KDTreeBox phiZ(minphi, maxphi, minv-0.01f, maxv+0.01f);  // declare our bounds
        //add fudge factors in case only one hit and also for floating-point inaccuracy
        hitTree[il].build(layerTree, phiZ); // make KDtree
        rzError[il] = maxErr; // save error
    }
*/

    /////////////
    //FITTING
    ////////////

/*
    const QuantityDependsPtEval maxChi2Eval = maxChi2.evaluator(es);
    const QuantityDependsPtEval extraPhiToleranceEval = extraPhiTolerance.evaluator(es);

    // re-used thoughout, need to be vectors because of RZLine interface
    std::vector<float> bc_r(4), bc_z(4), bc_errZ(4);

    // Loop over triplets
    for(const auto& triplet: triplets) {
        GlobalPoint gp0 = triplet.inner()->globalPosition();
        GlobalPoint gp1 = triplet.middle()->globalPosition();
        GlobalPoint gp2 = triplet.outer()->globalPosition();

        PixelRecoLineRZ line(gp0, gp2);
        ThirdHitPredictionFromCircle predictionRPhi(gp0, gp2, extraHitRPhitolerance);

        const double curvature = predictionRPhi.curvature(ThirdHitPredictionFromCircle::Vector2D(gp1.x(), gp1.y()));

        const float abscurv = std::abs(curvature);
        const float thisMaxChi2 = maxChi2Eval.value(abscurv);
        const float thisExtraPhiTolerance = extraPhiToleranceEval.value(abscurv);

        constexpr float nSigmaRZ = 3.46410161514f; // std::sqrt(12.f); // ...and continue as before

        auto isBarrel = [](const unsigned id) -> bool {
            return id == PixelSubdetector::PixelBarrel;
        };

        declareDynArray(GlobalPoint,4, gps);
        declareDynArray(GlobalError,4, ges);
        declareDynArray(bool,4, barrels);
        gps[0] = triplet.inner()->globalPosition();
        ges[0] = triplet.inner()->globalPositionError();
        barrels[0] = isBarrel(triplet.inner()->geographicalId().subdetId());

        gps[1] = triplet.middle()->globalPosition();
        ges[1] = triplet.middle()->globalPositionError();
        barrels[1] = isBarrel(triplet.middle()->geographicalId().subdetId());

        gps[2] = triplet.outer()->globalPosition();
        ges[2] = triplet.outer()->globalPositionError();
        barrels[2] = isBarrel(triplet.outer()->geographicalId().subdetId());

        for(size_t il=0; il!=size; ++il) {
            if(hitTree[il].empty()) continue; // Don't bother if no hits

            auto const& hits = *fourthHitMap[il];
            const DetLayer *layer = fourthLayers[il].detLayer();
            bool barrelLayer = layer->isBarrel();

            auto& predictionRZ = preds[il];
            predictionRZ.initPropagator(&line);
            Range rzRange = predictionRZ(); // z in barrel, r in endcap

            // construct search box and search
            Range phiRange;
            if(barrelLayer) {
                auto radius = static_cast<const BarrelDetLayer *>(layer)->specificSurface().radius();
                double phi = predictionRPhi.phi(curvature, radius);
                phiRange = Range(phi-thisExtraPhiTolerance, phi+thisExtraPhiTolerance);
            }
            else {
                double phi1 = predictionRPhi.phi(curvature, rzRange.min());
                double phi2 = predictionRPhi.phi(curvature, rzRange.max());
                phiRange = Range(std::min(phi1, phi2)-thisExtraPhiTolerance, std::max(phi1, phi2)+thisExtraPhiTolerance);
            }

            KDTreeBox phiZ(phiRange.min(), phiRange.max(),
                           rzRange.min()-nSigmaRZ*rzError[il],
                           rzRange.max()+nSigmaRZ*rzError[il]);

            foundNodes.clear();
            hitTree[il].search(phiZ, foundNodes);

            if(foundNodes.empty()) {
                continue;
            }

            SeedingHitSet::ConstRecHitPointer selectedHit = nullptr;
            float selectedChi2 = std::numeric_limits<float>::max();
            for(auto hitIndex: foundNodes) {
                const auto& hit = hits.theHits[hitIndex].hit();

                // Reject comparitor. For now, because of technical
                // limitations, pass three hits to the comparitor
                // TODO: switch to using hits from 2-3-4 instead of 1-3-4?
                // Eventually we should fix LowPtClusterShapeSeedComparitor to
                // accept quadruplets.
                if(theComparitor) {
                    SeedingHitSet tmpTriplet(triplet.inner(), triplet.outer(), hit);
                    if(!theComparitor->compatible(tmpTriplet, region)) {
                        continue;
                    }
                }

                gps[3] = hit->globalPosition();
                ges[3] = hit->globalPositionError();
                barrels[3] = isBarrel(hit->geographicalId().subdetId());

                float chi2 = std::numeric_limits<float>::quiet_NaN();
                // TODO: Do we have any use case to not use bending correction?
                if(useBendingCorrection) {
                    // Following PixelFitterByConformalMappingAndLine
                    const float simpleCot = ( gps.back().z()-gps.front().z() )/ (gps.back().perp() - gps.front().perp() );
                    const float pt = 1/PixelRecoUtilities::inversePt(abscurv, es);
                    for (int i=0; i< 4; ++i) {
                        const GlobalPoint & point = gps[i];
                        const GlobalError & error = ges[i];
                        bc_r[i] = sqrt( sqr(point.x()-region.origin().x()) + sqr(point.y()-region.origin().y()) );
                        bc_r[i] += pixelrecoutilities::LongitudinalBendingCorrection(pt,es)(bc_r[i]);
                        bc_z[i] = point.z()-region.origin().z();
                        bc_errZ[i] =  (barrels[i]) ? sqrt(error.czz()) : sqrt( error.rerr(point) )*simpleCot;
                    }
                    RZLine rzLine(bc_r,bc_z,bc_errZ);
                    float      cottheta, intercept, covss, covii, covsi;
                    rzLine.fit(cottheta, intercept, covss, covii, covsi);
                    chi2 = rzLine.chi2(cottheta, intercept);
                }
                else {
                    RZLine rzLine(gps, ges, barrels);
                    float  cottheta, intercept, covss, covii, covsi;
                    rzLine.fit(cottheta, intercept, covss, covii, covsi);
                    chi2 = rzLine.chi2(cottheta, intercept);
                }
                if(edm::isNotFinite(chi2) || chi2 > thisMaxChi2) {
                    continue;
                }
                // TODO: Do we have any use case to not use circle fit? Maybe
                // HLT where low-pT inefficiency is not a problem?
                if(fitFastCircle) {
                    FastCircleFit c(gps, ges);
                    chi2 += c.chi2();
                    if(edm::isNotFinite(chi2))
                    continue;
                    if(fitFastCircleChi2Cut && chi2 > thisMaxChi2)
                    continue;
                }


                if(chi2 < selectedChi2) {
                    selectedChi2 = chi2;
                    selectedHit = hit;
                }
            }
            if(selectedHit)
            result.emplace_back(triplet.inner(), triplet.middle(), triplet.outer(), selectedHit);
        }
    }
 */
}
