//
//  ParticleSwarmOptimiser.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 17/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#ifndef ParticleSwarmOptimiser_hpp
#define ParticleSwarmOptimiser_hpp

#include <iomanip>
#include "MatrixManip.hpp"

class ParticleSwarmOptimiser : public MatrixManip {
public:
    struct PSOparams {
        datatype phi1, phi2, phi, chi, w, wdamp, c1, c2;
        datatype velocityRangeScale;
        M1 velocityMin, velocityMax;
        int convreps = 8, verbose = 1;
        datatype convthreshold = (0.01)/100;
    };
    struct particleStruct {
        M1 position;
        M1 velocity;
        datatype cost;
        int* neighbour = nullptr;
        struct {
            M1 position;
            M1 velocity;
            datatype cost;
        } best, neighbourhoodbest;
    };
    
    int paramcount;
    particleStruct* particle = nullptr;
    particleStruct globalbest, localbest;
    Prior* prior;
    
    int particlecount;
    int neighbours = 5;
    int maxiters = 1e3;
    PSOparams psoparams;
    int bestcount = 0;
    datatype bestcost = INFINITY, lastbestcost = INFINITY, globalbestcost = INFINITY;
    
    ParticleSwarmOptimiser( int paramcount, Prior* prior );
    ~ParticleSwarmOptimiser();
    void initialise();
    void run();                  // Returns vectors of optimised parameter values
    virtual datatype costFunction( M1 state, Prior* prior, int n ) { return 0; };
    datatype getBestCost();
    M1 getBestPos();
};

#endif /* ParticleSwarmOptimiser_hpp */
