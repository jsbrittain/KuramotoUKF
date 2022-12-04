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
        datatype* state0 = nullptr;
        datatype** covarProposal = nullptr;
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
    datatype* state0 = nullptr;
    datatype** stateCurrent = nullptr;       // Allocate new memory block (may be running chains in parallel)
    datatype* stateNegLogLikeli = nullptr;
    datatype* stateProposal = nullptr;       //  |
    datatype** covarProposal = nullptr;
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
