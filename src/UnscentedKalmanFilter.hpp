//
//  UnscentedKalmanFilter.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 14/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#ifndef UnscentedKalmanFilter_hpp
#define UnscentedKalmanFilter_hpp

#include <iostream>
#include <iomanip>
#include <fstream>

#include <cmath>
#include <cassert>
#include <cstdlib>

#include "MatrixManip.hpp"

#define TRIGRES 65536
static const datatype dph = twopi/TRIGRES;
static const datatype phase2lookup = TRIGRES/twopi;

class UnscentedKalmanFilter : public MatrixManip {
public:
    int n_statevars, n_obs, N, t;	// N = sample count
    datatype **state = NULL;		// State variables [t,k]
    datatype **statepred = NULL;
    datatype ***stateP = NULL;      // State-transition covariance [t,k,k]
    datatype **ypred = NULL;		// Calculated (predicted) measurements
    datatype **y = NULL;			// Observed data (fixed)
    datatype costable[TRIGRES], sintable[TRIGRES];
    
    bool** statePmask=NULL;
    bool** obsPmask=NULL;
    bool** crossPmask=NULL;
    
    int n_sigmavecs;
    datatype **sigmaX = NULL;		// Sigma vectors [L',k]
    datatype **sigmaXpred = NULL;
    datatype **sigmaXpredMS = NULL;
    datatype **sigmasqrtP = NULL;
    datatype *sigmaWc = NULL;
    datatype *sigmaWm = NULL;
    datatype **sigmaY = NULL;
    datatype **sigmaYpredMS = NULL;
    datatype ***sigmaPcovarX = NULL;
    datatype ***sigmaPcovarY = NULL;
    datatype ***sigmaPcovarXY = NULL;
    datatype ***statePXXpred = NULL;
    datatype ***statePXYpred = NULL;
    datatype ***statePYYpred = NULL;
    datatype **stateNoise = NULL;
    datatype **obsNoise = NULL;
    datatype **PyyInv = NULL;
    datatype *ydiff = NULL;
    datatype **K = NULL;
    datatype **KT = NULL;
    datatype **PyyKT = NULL;
    datatype **KPyyKT = NULL;
    datatype **L = NULL;
    datatype **Linv = NULL;
    datatype **LinvT = NULL;
    datatype clambda, calpha, cbeta, ckappa;
    datatype sqrtPscaling;
    
    UnscentedKalmanFilter();
    ~UnscentedKalmanFilter();
    void reset();
    virtual void stateTransitionFunction( datatype* x, datatype* xpred ) {};
    virtual void observationFunction( datatype* x, datatype* y ) {};
    void genSigmaX( datatype* x, datatype** P, int dim, datatype** sigmaX );
    void genSigmaW( );
    void sqrtPscaled(datatype** L);
    void predict();
    void update();
    void initialise();
    void checkCholesky();
    void checkMatMultiply();
    void checkMatInv();
    void mathCheck();
    void saveStates(      const std::string filename );
    void saveStatesCovar( const std::string filename );
    void saveObs(         const std::string filename );
    void saveObsPred(     const std::string filename );
    void saveObsCovar(    const std::string filename );
    
    // These routines are used for generation, so don't have to be optimised
    void generate( int steps );
};

#endif /* UnscentedKalmanFilter_hpp */
