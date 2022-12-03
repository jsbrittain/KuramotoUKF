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
        datatype velocityRangeScale, *velocityMin, *velocityMax;
        int convreps = 8, verbose = 1;
        datatype convthreshold = (0.01)/100;
    };
    struct particleStruct {
        datatype* position = nullptr;
        datatype* velocity = nullptr;
        datatype cost;
        int* neighbour = nullptr;
        struct {
            datatype* position = nullptr;
            datatype* velocity = nullptr;
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
    virtual datatype costFunction( datatype* state, Prior* prior, int n ) { return 0; };
    datatype getBestCost();
    datatype* getBestPos();
};

#endif /* ParticleSwarmOptimiser_hpp */
