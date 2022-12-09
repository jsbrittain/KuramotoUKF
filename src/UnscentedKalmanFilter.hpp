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
    M2 state;		// State variables [t,k]
    M2 statepred;
    M3 stateP;      // State-transition covariance [t,k,k]
    M2 ypred;		// Calculated (predicted) measurements
    M2 y;			// Observed data (fixed)
    datatype costable[TRIGRES], sintable[TRIGRES];
    
    std::vector<std::vector<bool>> statePmask;
    std::vector<std::vector<bool>> obsPmask;
    std::vector<std::vector<bool>> crossPmask;
    
    int n_sigmavecs;
    M2 sigmaX;		// Sigma vectors [L',k]
    M2 sigmaXpred;
    M2 sigmaXpredMS;
    M2 sigmasqrtP;
    M1 sigmaWc;
    M1 sigmaWm;
    M2 sigmaY;
    M2 sigmaYpredMS;
    M3 sigmaPcovarX;
    M3 sigmaPcovarY;
    M3 sigmaPcovarXY;
    M3 statePXXpred;
    M3 statePXYpred;
    M3 statePYYpred;
    M2 stateNoise;
    M2 obsNoise;
    M2 PyyInv;
    M1 ydiff;
    M2 K;
    M2 KT;
    M2 PyyKT;
    M2 KPyyKT;
    M2 L;
    M2 Linv;
    M2 LinvT;
    datatype clambda, calpha, cbeta, ckappa;
    datatype sqrtPscaling;
    
    UnscentedKalmanFilter();
    ~UnscentedKalmanFilter();
    void reset();
    virtual void stateTransitionFunction( M1 x, M1& xpred ) = 0;
    virtual void observationFunction( M1 x, M1& y ) = 0;
    void genSigmaX( M1 x, M2 P, int dim, M2 sigmaX );
    void genSigmaW( );
    void sqrtPscaled(M2 L);
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
