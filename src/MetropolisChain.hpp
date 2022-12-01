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
    Prior* prior = NULL;
    
    struct MetropolisChainParams {
        datatype* state0 = NULL;
        datatype** covarProposal = NULL;
        bool tuneProposal = false;
        int tuningiters = NAN;
        Prior* prior;
        int statedim = NAN;
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
    datatype* state0 = NULL;
    datatype** stateCurrent = NULL;       // Allocate new memory block (may be running chains in parallel)
    datatype* stateNegLogLikeli = NULL;
    datatype* stateProposal = NULL;       //  |
    datatype** covarProposal = NULL;
    bool tuneProposal = false;
    int tuningiters = 0;
    datatype logPcurrent, logPproposal, logPcrit;
    
    MetropolisChain( MetropolisChainParams mcmcparams );
    MetropolisChain( datatype* state0, datatype** covarProposal, bool tuneProposal, int tuningiters, Prior* prior, int statedim, int randseed, int burnin, int chainlength, int verbose );
    ~MetropolisChain();
    void initialise();
    void run( );
    datatype steps( int maxcount );
    bool step( );
    void setVerbose( int val );
    datatype negLogLikeliTestFcn( datatype* state );
    void saveChain( const std::string filename );
    datatype getMAPcost( );
    datatype* getMAPstate( );
    
    // Log-likelihood function --- needs to be implemented by children
    virtual datatype negLogLikeliFcn( datatype* state, Prior* prior, int n ) { return NAN; };
};

#endif /* MetropolisChain */
