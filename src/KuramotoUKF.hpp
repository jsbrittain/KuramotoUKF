//
//  KuramotoUKF.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 14/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#ifndef KuramotoUKF_hpp
#define KuramotoUKF_hpp

#include <iostream>
#include <fstream>

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>

#include "UnscentedKalmanFilter.hpp"
#include "MatrixManip.hpp"

using namespace std;

class KuramotoUKF : public UnscentedKalmanFilter {
public:
    int nodecount;     // -1 is no connection
    int feedbacklag = 0, feedback_obs_chan = 0;
    datatype feedbackStrength = 0.0, feedbacklag_secs = 0.0, feedbacklag_secs_min = 0.0;
    datatype feedback_lagmag, feedback_lagphase, lag_xx, lag_yy;
    datatype nodecountf;
    datatype dt, ycos;
    datatype ascaling;
    datatype* yGivenX = nullptr;
    int* paramPriorList = nullptr;
    Prior* prior = nullptr;
    int n_params, n_priors;
    
    struct {
        int* phase0;
        int* instfreq;
        int connK = -1;
        int feedbackStrength = -1;
        int feedbackLagSecs = -1;
        int ascaling = -1;
        int* phase0Sigma;
        int* instfreqSigma;
        int connKSigma = -1;
        int feedbackStrengthSigma = -1;
        int* phaseSigma;
        int* instfreq0Sigma;
        int connK0Sigma = -1;
        int feedbackStrength0Sigma = -1;
        int ySigma = -1;
    } paramindex;
    
    struct KurfParams {
        // Provide defaults
        int    nodecount   = 3;
        int    n_statevars = 7;
        datatype alpha       = 1e-3;
        datatype beta        = 2.0;
        datatype kappa       = 0.0;
        datatype dt          = 1;
    };
    
    enum ParamMode { none, fixed, parameter, variable } mode;
    enum Masktype { zero, diagonal, full, block, defaultmask };
    struct ModelVar {
        ParamMode mode;             // Using default values here makes the structure non-aggregate
        Masktype covarMask;
        datatype value, P, P0;      //  preventing easy initialisation per instance below
    };
    
    struct ModelParamsSimple {
        KurfParams kurfparams;
        
        enum GenFromPriorsMode { means, noisy } gen_from_priors_mode;
        int nodecount      = 3;
        int n_obs          = 1;
        Masktype covarMask = Masktype::full;
        datatype feedbacklag_secs_min = 0.050;
        int      feedback_obs_chan    = 0;
        
        ModelVar theta = {
            .mode=variable, covarMask=Masktype::defaultmask,
            .value=NAN, .P=0.00001, .P0=0.00001
        };
        ModelVar omega = {
            .mode=variable, covarMask=Masktype::defaultmask,
            .value=2*pi*5.0, .P=0.00001, .P0=0.00001
        };
        ModelVar connK = {
            .mode=variable, covarMask=Masktype::defaultmask,
            .value=1, .P=0.00001, .P0=0.00001
        };
        ModelVar feedback_strength = {
            .mode=parameter, covarMask=Masktype::defaultmask,
            .value=0.1, .P=0.00001, .P0=0.00001,
        };
        ModelVar feedbacklag_secs = {
            .mode=parameter, covarMask=Masktype::defaultmask,
            .value=0.051, .P=0.00001, .P0=0.00001,
        };
        ModelVar ascaling = {
            .mode=parameter, covarMask=Masktype::defaultmask,
            .value=0.0014, .P=0.00001, .P0=0.00001
        };
        ModelVar sigmay = {
            .mode=parameter, covarMask=Masktype::defaultmask,
            .value=0.000001, .P=0.00001, .P0=0.00001
        };
    };
    
    struct stateConditions {
        datatype* x;
        datatype** P;
        datatype** stateNoise;
        datatype** obsNoise;
        datatype feedbackStrength = 0.0;
        datatype feedbacklag_secs = 0.0;
        datatype ascaling = 1.0;
    };
    
    KuramotoUKF( int nodecount, int n_statevars, datatype calpha, datatype cbeta, datatype ckappa, datatype dt );
    KuramotoUKF( KurfParams kurfparams );
    ~KuramotoUKF();
    void initialise();
    void initialise( int n_obs, int samplecount );
    void stateTransitionFunction( datatype* x, datatype* xpred ) override;
    void observationFunction( datatype* x, datatype* y ) override;
    datatype cosT( datatype x );
    datatype sinT( datatype x );
    void readMeasurementFile( const string filename, int n_obs );
    void predict() { UnscentedKalmanFilter::predict(); }
    void update() { UnscentedKalmanFilter::update(); }
    void generate( datatype T ) { UnscentedKalmanFilter::generate( (int) T/dt ); };
    void generate( int N ) { UnscentedKalmanFilter::generate( N ); };
    void applymask( datatype** X, int dim1, int dim2, int** mask );
    void printState( );
    void printPrediction( );
    void printSigmaX();
    void printK();
    void setInitialConditionsFromPriors( );
    void setInitialConditions( stateConditions* initcond );
    void setInitialConditions( datatype* x0, datatype** P0, datatype** stateNoise, datatype** obsNoise, datatype ascaling, datatype feedbackStrength, datatype feedbacklag_secs );
    void checkCholesky() { UnscentedKalmanFilter::checkCholesky(); }
    void testTranspose( ) ;
    void mathCheck( );
    datatype run( );
    void saveStates( const string filename ) { UnscentedKalmanFilter::saveStates( filename ); };
    void saveObs( const string filename ) { UnscentedKalmanFilter::saveObs( filename ); };
    void saveObsPred( const string filename ) { UnscentedKalmanFilter::saveObsPred( filename ); };
    void testTrigonometricApproximations( const string filename );
    void testRandomNumberGenerators( const string filestem ) { UnscentedKalmanFilter::testRandomNumberGenerators( filestem ); };
    int getSampleCount() { return UnscentedKalmanFilter::N; };
    void reset() { UnscentedKalmanFilter::reset(); };
    datatype negLogLikeliFcn( datatype* x, Prior* prior, int n );
    datatype negLogLikeliFcn( datatype* x );
    void formPriors( KuramotoUKF::ModelParamsSimple params );
    bool** initCovarMask( int dim1, int dim2, Masktype masktype );
    void updateCovarMasks_StateVar( int k, Masktype masktype );
    void updateCovarMasks_StateVar( int k, Masktype masktype, int* blockindices, int blockcount );
    bool isMaskFull( bool** mask, int dim1, int dim2 );
    void setPriors( int* paramPriorList, int n_priors, int n_params );
    Prior* getPriors( );
    int* getParamPriorList( );
    void printStateCond( stateConditions* initcond );
    datatype* priorsToPriorVec( Prior* prior, int* n_priors );
    stateConditions* unpackParamVec( datatype* paramvec );
    datatype* priorVecToParamVec( datatype* priorvec, int* paramPriorList );
    datatype* rmsy();
    bool isVariable( ModelVar var );
    bool isParameter( ModelVar var );
};

#endif /* KuramotoUKF_hpp */
