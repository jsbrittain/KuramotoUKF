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
    datatype **state = nullptr;		// State variables [t,k]
    datatype **statepred = nullptr;
    datatype ***stateP = nullptr;      // State-transition covariance [t,k,k]
    datatype **ypred = nullptr;		// Calculated (predicted) measurements
    datatype **y = nullptr;			// Observed data (fixed)
    datatype costable[TRIGRES], sintable[TRIGRES];
    
    bool** statePmask=nullptr;
    bool** obsPmask=nullptr;
    bool** crossPmask=nullptr;
    
    int n_sigmavecs;
    datatype **sigmaX = nullptr;		// Sigma vectors [L',k]
    datatype **sigmaXpred = nullptr;
    datatype **sigmaXpredMS = nullptr;
    datatype **sigmasqrtP = nullptr;
    datatype *sigmaWc = nullptr;
    datatype *sigmaWm = nullptr;
    datatype **sigmaY = nullptr;
    datatype **sigmaYpredMS = nullptr;
    datatype ***sigmaPcovarX = nullptr;
    datatype ***sigmaPcovarY = nullptr;
    datatype ***sigmaPcovarXY = nullptr;
    datatype ***statePXXpred = nullptr;
    datatype ***statePXYpred = nullptr;
    datatype ***statePYYpred = nullptr;
    datatype **stateNoise = nullptr;
    datatype **obsNoise = nullptr;
    datatype **PyyInv = nullptr;
    datatype *ydiff = nullptr;
    datatype **K = nullptr;
    datatype **KT = nullptr;
    datatype **PyyKT = nullptr;
    datatype **KPyyKT = nullptr;
    datatype **L = nullptr;
    datatype **Linv = nullptr;
    datatype **LinvT = nullptr;
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
