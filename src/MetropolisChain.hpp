//
//  MetropolisChain
//  mcmc
//
//  Created by John-Stuart Brittain on 14/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#ifndef MetropolisChain_hpp
#define MetropolisChain_hpp

#include <thread>
#include <iomanip>

#include "MatrixManip.hpp"

class MetropolisChain : public MatrixManip {
public:
    Prior* prior = nullptr;
    
    struct MetropolisChainParams {
        M1 state0;
        M2 covarProposal;
        bool tuneProposal = false;
        int tuningiters = 0;
        Prior* prior;
        int statedim = 0;
        int randseed = 0;
        int burnin = 200;
        int chainlength = 3000;
        int verbose = 0;
    };

    int randseed;
    int burnin;
    int chainlength;
    int statedim;
    int chaint;
    int verbose = 0;
    M1 state0;
    M2 stateCurrent;       // Allocate new memory block (may be running chains in parallel)
    M1 stateNegLogLikeli;
    M1 stateProposal;       //  |
    M2 covarProposal;
    bool tuneProposal = false;
    int tuningiters = 0;
    datatype logPcurrent, logPproposal, logPcrit;
    
    MetropolisChain( MetropolisChainParams mcmcparams );
    MetropolisChain( M1 state0, M2 covarProposal, bool tuneProposal, int tuningiters, Prior* prior, int statedim, int randseed, int burnin, int chainlength, int verbose );
    ~MetropolisChain();
    void initialise();
    void run( );
    datatype steps( int maxcount );
    bool step( );
    void setVerbose( int val );
    datatype negLogLikeliTestFcn( M1 state );
    void saveChain( const std::string filename );
    datatype getMAPcost( );
    M1 getMAPstate( );
    
    // Log-likelihood function --- needs to be implemented by children
    virtual datatype negLogLikeliFcn( M1 state, Prior* prior, int n ) { return NAN; };
};

#endif /* MetropolisChain */
