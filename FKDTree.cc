#include "FKDTree.h"

#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <utility>
#include <iostream>
#include <cassert>
#include <deque>
#include "FKDPoint.h"
#include "FQueue.h"
#include "DataFormats/GeometryVector/interface/Pi.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"

FKDTree::FKDTree(const std::vector<Hit>& hits, GlobalPoint const & origin, DetLayer const * il) :

    theNumberOfPoints(2*hits.size()),
    theDepth(FLOOR_LOG2(theNumberOfPoints)),
    theDimensions[0](theNumberOfPoints,0),
    theDimensions[1](theNumberOfPoints,0),
    theDimensions[2](theNumberOfPoints,0),
    theIntervalLength(theNumberOfPoints, 0),
    theIntervalMin(theNumberOfPoints, 0),
    theIds(theNumberOfPoints,0)
    thePoints(theNumberOfPoints)

{
    unsigned int pointID = 0;
    for (int i = 0; i!=(int)hits.size(); i++) {
        Hit const & hit = hits[i]->hit();
        auto const & gs = hit->globalState();
        auto loc = gs.position-origin.basicVector();
        
        auto phi = gs.position.barePhi();
        auto z = gs.position.z();
        float r = loc.perp();
        
        std::cout<<"Make Point? - Phi = "<<phi<<" Z = "<<z<<" R = "<<r<<std::endl;
        thePoints.push_back(make_FKDPoint(phi,z,r,pointID)); pointID++;
        std::cout<<"Made!"<<std::endl;
        std::cout<<"Point : "<<points[i][0]<<" - "<<points[i][0]<<" - "<<points[i][1]<<" - "<<points[i][2]<<std::endl;
        if (phi>safePhi) {thePoints.push_back(make_FKDPoint(phi-Geom::ftwoPi(),z,r,pointID)); pointID++;}
        else if (phi<-safePhi) {thePoints.push_back(make_FKDPoint(phi+Geom::ftwoPi(),z,r,pointID));pointID++;}
        
    }
    theNumberOfPoints = pointID;
    theIntervalMin.shrink_to_fit();
    theIntervalMin.shrink_to_fit();
    theDimensions[0].shrink_to_fit();
    theDimensions[1].shrink_to_fit();
    theDimensions[2].shrink_to_fit();
    theDepth = FLOOR_LOG2(theNumberOfPoints);
    thePoints.shrink_to_fit();
    // standard region have origin as 0,0,z (not true!!!!0
    // cosmic region never used here
    // assert(origin.x()==0 && origin.y()==0);

    for (unsigned int i=0; i!=hits.size(); ++i) {
        auto const & h = *hits[i].hit();
        auto const & gs = static_cast<BaseTrackerRecHit const &>(h).globalState();
        auto loc = gs.position-origin.basicVector();
        float lr = loc.perp();
        // float lr = gs.position.perp();
        float lz = gs.position.z();
        float dr = gs.errorR;
        float dz = gs.errorZ;
        // r[i] = gs.position.perp();
        // phi[i] = gs.position.barePhi();
        x[i] = gs.position.x();
        y[i] = gs.position.y();
        z[i] = lz;
        drphi[i] = gs.errorRPhi;
        u[i] = isBarrel ? lr : lz;
        v[i] = isBarrel ? lz : lr;
        du[i] = isBarrel ? dr : dz;
        dv[i] = isBarrel ? dz : dr;
        lphi[i] = loc.barePhi();
    }
    
}


