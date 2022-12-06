//
//  KuramotoUKF.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 14/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#include "KuramotoUKF.hpp"

KuramotoUKF::KuramotoUKF( int nodecount, int n_statevars, datatype alpha, datatype beta, datatype kappa, datatype dt ) : nodecount(nodecount), dt(dt) {
    calpha = alpha;
    cbeta = beta;
    ckappa = kappa;
    UnscentedKalmanFilter::n_statevars = n_statevars;
}
KuramotoUKF::KuramotoUKF( KurfParams kurfparams ) : nodecount(kurfparams.nodecount), dt(kurfparams.dt) {
    calpha = kurfparams.alpha;
    cbeta = kurfparams.beta;
    ckappa = kurfparams.kappa;
    UnscentedKalmanFilter::n_statevars = kurfparams.n_statevars;
}
KuramotoUKF::~KuramotoUKF() {
    //deallocMatrix(y,N,n_obs);
    //deallocMatrix(yGivenX,n_obs);
}
void KuramotoUKF::initialise( int n_obs, int samplecount ) {
    UnscentedKalmanFilter::n_obs = n_obs;
    N = samplecount;
    initialise();
    y = allocMatrix(N,n_obs);
}
void KuramotoUKF::initialise( ) {
    assert(N>0);
    // Allocate memory
    //n_statevars = 2*nodecount+1;		// These are (theta,omega,K)
    //stateindexK = n_statevars-1;
    nodecountf = (datatype) nodecount;
    feedbacklag = (int)(feedbacklag_secs/dt);
    UnscentedKalmanFilter::initialise();
    yGivenX = allocMatrix(n_obs);
}
void KuramotoUKF::stateTransitionFunction( M1 x, M1 xpred ) {
    // Check for sensory feedback
    if ( paramindex.feedbackStrength != -1 ) {
        // Calculate average phase and resultant magnitude from # steps back
        feedback_lagmag = 0; feedback_lagphase = 0;
        if ( feedbacklag < 0 ) {
            for ( int k = nodecount; k < n_statevars; k++ )
                xpred[k] = NAN;
            return;
        }
        if ( t > feedbacklag ) {
            // Determine resultant phase vector
            lag_xx = 0; lag_yy = 0;
            for (int k=0; k<nodecount; k++) {
                lag_xx += cosT(state[t-feedbacklag][k]);
                lag_yy += sinT(state[t-feedbacklag][k]);
            }
            lag_xx /= nodecountf; lag_yy /= nodecountf;
            // Calculate magnitude and phase of resultant vector
            feedback_lagmag = sqrt( lag_xx*lag_xx + lag_yy*lag_yy );
            feedback_lagphase = atan2( lag_yy, lag_xx );
        }
    }
    
    // Oscillators
    for (int k=0; k<nodecount; k++) {
        // Phase interaction function (PIF)
        xpred[k] = 0.0;
        for (int j=0; j<nodecount; j++) {
            if ( j != k )
                xpred[k] += sinT( x[j] - x[k] );        // Negative phases WILL occur here
        }
        // Connectivity strength (logarithm, +ve)
        if ( paramindex.connK != -1 )
            xpred[k] *= -exp(x[paramindex.connK])/nodecountf;
        
        // Advance by endogenous frequency
        xpred[k] += x[k+nodecount];
        // Peripheral feedback (scaled phase-interaction with feedback)
        if (( paramindex.feedbackStrength != -1 ) && ( t > feedbacklag )) {
            xpred[k] += -exp(feedbackStrength)*feedback_lagmag*sinT(feedback_lagphase - x[k]);
        }
        
        // Scale all (dx/dt) components by sampling rate
        xpred[k] *= dt;
        // Add current state
        xpred[k] += x[k];
        // Do NOT wrap as this messes up the sigma-vectors in transition
    }
    // All other nodes are propagated unaltered (random walk)
    for ( int k = nodecount; k < n_statevars; k++ )
        xpred[k] = x[k];
}
void KuramotoUKF::observationFunction( M1 x, M1 y ) {
    // Oscillator nodes
    y[0] = 0.0;
    for (int k=0; k<nodecount; k++) {
        y[0] += cosT( x[k] );
    }
    // Global scaling
    y[0] *= ascaling/nodecountf;
}
datatype KuramotoUKF::cosT( datatype x ) {
    return cos( x );
}
datatype KuramotoUKF::sinT( datatype x ) {
    return sin( x );
}
void KuramotoUKF::readMeasurementFile( const string filename, int n_obs ) {
    datatype num;
    
    UnscentedKalmanFilter::n_obs = n_obs;
    assert( n_obs == 1 );		// Only one measurement currently supported
    
    // Open text file
    cout << "Reading measurement file:";
    std::fstream myfile(filename, std::ios_base::in);
    
    // Get data length
    N = 0;
    while (myfile >> num)			// Count number of elements
        N++;
    assert( N > 0 );
    
    // Allocate memory
    y = allocMatrix( N, n_obs );
    
    // Reset file read position
    myfile.clear();					// Clear flags
    myfile.seekg(0);				//  then seek
    cout << " " << N << " samples...";
    
    // Read file
    int k=0;
    while (myfile >> y[k++][0]) ;
    cout << " done." << endl;
}
void KuramotoUKF::printState() {
    cout << "State at time step t = " << t << endl;
    printVector( state[t] );
    printMatrix( stateP[t] );
    printVector( y[t] );
}
void KuramotoUKF::printPrediction() {
    cout << "Predicted state at time step t = " << t << endl;
    printVector( statepred[t] );
    printMatrix( statePXXpred[t] );
    printVector( ypred[t] );
}
void KuramotoUKF::printSigmaX() {
    cout << "SigmaX matrix for time step t = " << t << endl;
    genSigmaX( state[t], stateP[t], n_statevars, sigmaX );
    for ( int k = 0; k < n_sigmavecs; k ++ ) {
        cout << k << ": " << endl;
        printVector( sigmaX[k] );
    }
    cout << "SigmaW scaling" << t << endl;
    genSigmaW();
    for ( int k = 0; k < n_sigmavecs; k ++ ) {
        cout << k << ": " << sigmaWm[k] << ", " << sigmaWc[k] << endl;
    }
}
void KuramotoUKF::printK( ) {
    cout << "Kalman gain K for time step t = " << t << endl;
    printMatrix(K);
}
void KuramotoUKF::setInitialConditionsFromPriors() {
    int n_priors;
    M1 priorvec = priorsToPriorVec( prior, &n_priors );
    M1 paramvec = priorVecToParamVec( priorvec, paramPriorList );
    stateConditions* statecond = unpackParamVec(paramvec);
    KuramotoUKF::setInitialConditions( statecond );
}
void KuramotoUKF::setInitialConditions( stateConditions* statecond ) {
    KuramotoUKF::setInitialConditions( statecond->x, statecond->P, statecond->stateNoise, statecond->obsNoise, statecond->ascaling, statecond->feedbackStrength, statecond->feedbacklag_secs );
}
void KuramotoUKF::setInitialConditions( M1 x0, M2 P0, M2 stateNoise, M2 obsNoise, datatype ascaling, datatype feedbackStrength, datatype feedbacklag_secs ) {
    
    assert( nodecount > 0 );
    assert( n_statevars > 0 );
    assert( n_obs > 0 );
    //assert( n_priors > 0 );
    
    t = 0;
    for ( int i = 0; i < n_statevars; i++ )
        state[t][i] = x0[i];
    KuramotoUKF::ascaling = ascaling;
    KuramotoUKF::feedbackStrength = feedbackStrength;
    KuramotoUKF::feedbacklag_secs = feedbacklag_secs;
    KuramotoUKF::feedbacklag = round(feedbacklag_secs/dt);
    
    // Scale noise
    datatype sqrtdt = sqrt(dt);
    UnscentedKalmanFilter::matmultbyscalar( P0, sqrtdt, stateP[t] );
    UnscentedKalmanFilter::matmultbyscalar( stateNoise,  sqrtdt, UnscentedKalmanFilter::stateNoise );
    UnscentedKalmanFilter::matmultbyscalar( obsNoise, sqrtdt*(KuramotoUKF::ascaling*KuramotoUKF::ascaling), UnscentedKalmanFilter::obsNoise );
}
void KuramotoUKF::applymask( M2 X, int dim1, int dim2, int** mask ) {
    for ( int i = 0; i < dim1; i++ ) {
        for ( int j = 0; j < dim2; j++ ) {
            if ( mask[i][j] == 0 )
                X[i][j] = 0.0;
        }
    }
}
void KuramotoUKF::mathCheck() {
    UnscentedKalmanFilter::mathCheck();
}
datatype KuramotoUKF::run() {
    datatype logLikeli = 0.0;
    
    // Set time to zero
    reset();
    
    // Model estimation
    for ( int n = 0; n < N-1 ; n++ ) {
        // Prediction step (priors)
        predict();
        // Update step (posteriors)
        update();
    }
    
    // Run the next section exactly as specified --- the call to persistant maths functions will preserve memory allocation so long as the previous call was of the same or larger dimension. If not, memory is deallocated and reallocated dynamically.
    
    // Evaluate log-likelihood of states given previous states: p( x(k) | x(k-1) )
    logLikeli += logLikeliMVNpersistent(state[0], n_statevars, state[0], stateP[0]);
    for ( int t = 1; t < N; t++ )
        logLikeli += logLikeliMVNpersistent(state[t], n_statevars, statepred[t], statePXXpred[t]);
    
    // Log-likelihood of observations given states: p( y(k) | x(k) )
    for ( int t = 0; t < N; t++ ) {
        observationFunction( state[t], yGivenX );
        if ( !isnan(y[t][0]) )
            logLikeli += logLikeliMVNpersistent(y[t], n_obs, yGivenX, obsNoise);
    }
    
    // Return log-likelihood of model
    return logLikeli;
}
void KuramotoUKF::testTrigonometricApproximations( const string filename ) {
    int K = 1e4;
    datatype dph = twopi/(K/3.0);
    datatype phaseoffset = -pi;
    M2 x = allocMatrix(K, 3);
    for ( int k = 0; k < K; k++ ) {
        // Col 0: Phase
        x[k][0] = phaseoffset + k*dph;
        // Col 1: Sine
        x[k][1] = sinT( x[k][0] );
        // Col 2: Cos
        x[k][2] = cosT( x[k][0] );
    }
    writeMatrixToFile(filename, x, K, 3);
}
datatype KuramotoUKF::negLogLikeliFcn( M1 x, Prior* prior, int n ) {
    return negLogLikeliFcn( x ) - logLikeliPriors( x, prior, n );
}
datatype KuramotoUKF::negLogLikeliFcn( M1 x ) {
    M1 paramvec = priorVecToParamVec(x, paramPriorList);
    stateConditions* initCond = unpackParamVec(paramvec);
    
    // Set conditions and return (negative) log-likelihood
    setInitialConditions( initCond );
    datatype logLikeli = run();
    
    return -logLikeli;
}
void KuramotoUKF::formPriors( KuramotoUKF::ModelParamsSimple params ) {
    
    /* Ordering:
            statevector  means
            statevector  Sigma0
            statevector  Sigma (transition)
            other params ( ascaling, sigmay )
     
       This ordering requires multiple definitions to account for parameter or variable definitions
     */
    
    // Clean-out old priors first
    if ( prior != nullptr ) delete[] prior;
    if ( paramPriorList != nullptr ) delete[] paramPriorList;
    
    // Determine some parameters and set in class
    nodecount = params.nodecount;
    n_statevars = 2*nodecount +
        ((int) isVariable(params.connK));
    n_params = 3*n_statevars + 2 + 2*((int) isParameter(params.feedback_strength));
    int n_obs = params.n_obs;
    feedbacklag_secs_min = params.feedbacklag_secs_min;
    
    // Initialise covariance masks
    statePmask = initCovarMask( n_statevars, n_statevars, params.covarMask );
    obsPmask = initCovarMask( n_obs, n_obs, params.covarMask );
    //crossPmask = initCovarMask( n_statevars, n_obs,
    //                (params.covarMask==Masktype::full) ? (Masktype::full) : (Masktype::zero) );
    crossPmask = initCovarMask( n_statevars, n_obs, Masktype::full );
    
    // Form new priors lists
    paramPriorList = new int[n_params];
    prior = new MatrixManip::Prior[n_params];
    int* blockindices = nullptr;
    
    
    
    // Starting phase and inst freqs
    int k;
    paramindex.phase0 = new int[nodecount];
    blockindices = new int[nodecount];
    for ( k = 0; k < nodecount; k++ ) {
        // Initial phase values (linear; fixed in state transition)
        prior[k].mu = 2.0*pi*k/nodecount; prior[k].sd = pi;
        prior[k].varclass = MatrixManip::Prior::VarClass::linear;
        blockindices[k] = k;
        paramindex.phase0[k] = k;
        paramPriorList[k] = k;
    }
    updateCovarMasks_StateVar( 0, params.theta.covarMask, blockindices, nodecount );
    k = nodecount;      // Next entry; prior list
    int j = k;          // Next entry; param list
    
    // Instantaneous frequencies (single estimate; linear; fixed in state transition)
    prior[k].mu = params.omega.value;
    prior[k].sd = prior[k].mu/10;
    prior[k].varclass = MatrixManip::Prior::VarClass::linear;
    paramindex.instfreq = new int[nodecount];
    for ( int i = 0; i < nodecount; i++ ) {
        //updateCovarMasks_StateVar( j, params.omega.covarMask );
        blockindices[i] = nodecount + i;
        paramindex.instfreq[i] = j;
        paramPriorList[j] = k;
        j++;
    }
    updateCovarMasks_StateVar( nodecount, params.omega.covarMask, blockindices, nodecount );
    delete[] blockindices;
    k++;
    
    // Coupling strength K (log-linear; fixed in state transition)
    if ( isVariable( params.connK ) ) {
        prior[k].mu = params.connK.value; prior[k].sd = prior[k].mu;
        prior[k].varclass = MatrixManip::Prior::VarClass::loglinear;
        updateCovarMasks_StateVar( j, params.connK.covarMask );
        paramindex.connK = j;
        paramPriorList[j++] = k++;
    } else
        paramindex.connK = -1;
    
    // Delayed feedback node
    if ( isVariable(params.feedback_strength) ) {
        prior[k].mu = params.feedback_strength.value; prior[k].sd = prior[k].mu/10;
        prior[k].varclass = MatrixManip::Prior::VarClass::loglinear;
        updateCovarMasks_StateVar( j, params.feedback_strength.covarMask );
        paramindex.feedbackStrength = j;
        paramPriorList[j++] = k++;
    } else
        paramindex.feedbackStrength = -1;
    
    
    
    /// Initial variance (invariant to sampling interval) --- some share common prior ///
    
    // Initial phase (variance)
    prior[k].mu = params.theta.P0; prior[k].sd = prior[k].mu;
    prior[k].varclass = MatrixManip::Prior::VarClass::variance;
    paramindex.phase0Sigma = new int[nodecount];
    for ( int i = 0; i < nodecount; i++ ) {
        //updateCovarMasks_StateVar( j, params.theta.covarMask );
        paramindex.phase0Sigma[i] = j;
        paramPriorList[j] = k;
        j++;
    }
    k++;
    // Inst freq (variance)
    prior[k].mu = params.omega.P0; prior[k].sd = prior[k].mu;
    prior[k].varclass = MatrixManip::Prior::VarClass::variance;
    paramindex.instfreq0Sigma = new int[nodecount];
    for ( int i = 0; i < nodecount; i++ ) {
        //updateCovarMasks_StateVar( j, params.omega.covarMask );
        paramindex.instfreq0Sigma[i] = j;
        paramPriorList[j] = k;
        j++;
    }
    k++;
    // Connection strength (variance)
    if ( isVariable( params.connK ) ) {
        prior[k].mu = params.connK.P0; prior[k].sd = prior[k].mu;
        prior[k].varclass = MatrixManip::Prior::VarClass::variance;
        //updateCovarMasks_StateVar( j, params.connK.covarMask );
        paramindex.connK0Sigma = j;
        paramPriorList[j++] = k++;
    }
    // Delayed-feedback strengh (variance)
    if ( isVariable( params.feedback_strength ) ) {
        prior[k].mu = params.feedback_strength.P0; prior[k].sd = prior[k].mu/10;
        prior[k].varclass = MatrixManip::Prior::VarClass::variance;
        //updateCovarMasks_StateVar( j, params.feedback_strength.covarMask );
        paramindex.feedbackStrength0Sigma = j;
        paramPriorList[j++] = k++;
    }
    
    
    
    /// State-transition variance (invariant to sampling interval) --- some share common prior ///
    
    // Phase (variance)
    prior[k].mu = params.theta.P; prior[k].sd = prior[k].mu/10;
    prior[k].varclass = MatrixManip::Prior::VarClass::variance;
    paramindex.phaseSigma = new int[nodecount];
    for ( int i = 0; i < nodecount; i++ ) {
        //updateCovarMasks_StateVar( j, params.theta.covarMask );
        paramindex.phaseSigma[i] = j;
        paramPriorList[j] = k;
        j++;
    }
    k++;
    // Inst freq (variance)
    prior[k].mu = params.omega.P; prior[k].sd = prior[k].mu/10;
    prior[k].varclass = MatrixManip::Prior::VarClass::variance;
    paramindex.instfreqSigma = new int[nodecount];
    for ( int i = 0; i < nodecount; i++ ) {
        //updateCovarMasks_StateVar( j, params.omega.covarMask );
        paramindex.instfreqSigma[i] = j;
        paramPriorList[j] = k;
        j++;
    }
    k++;
    // Connection strength (variance)
    if ( isVariable( params.connK ) ) {
        prior[k].mu = params.connK.P; prior[k].sd = prior[k].mu/10;
        prior[k].varclass = MatrixManip::Prior::VarClass::variance;
        //updateCovarMasks_StateVar( j, params.connK.covarMask );
        paramindex.connKSigma = j;
        paramPriorList[j++] = k++;
    }
    // Delayed-feedback strengh (variance)
    if ( isVariable( params.feedback_strength ) ) {
        prior[k].mu = params.feedback_strength.P; prior[k].sd = prior[k].mu/10;
        prior[k].varclass = MatrixManip::Prior::VarClass::variance;
        //updateCovarMasks_StateVar( j, params.feedback_strength.covarMask );
        paramindex.feedbackStrengthSigma = j;
        paramPriorList[j++] = k++;
    }
    
    
    
    // Parameters (not associated with the state vector)
    
    if ( isParameter( params.connK ) ) {
        prior[k].mu = params.connK.value; prior[k].sd = prior[k].mu;
        prior[k].varclass = MatrixManip::Prior::VarClass::loglinear;
        //updateCovarMasks_StateVar( j, params.connK.covarMask );
        paramindex.connK = j;
        paramPriorList[j++] = k++;
    } else
        paramindex.connK = -1;
    
    if ( isParameter( params.feedback_strength ) ) {
        prior[k].mu = params.feedback_strength.value; prior[k].sd = prior[k].mu;
        prior[k].varclass = MatrixManip::Prior::VarClass::loglinear;
        //updateCovarMasks_StateVar( j, params.feedback_strength.covarMask );
        paramindex.feedbackStrength = j;
        paramPriorList[j++] = k++;
    } else
        paramindex.feedbackStrength = -1;
    
    // Feedback lag is required to be a parameter (if present; for now)
    if ( params.feedback_strength.mode != ParamMode::none ) {
        prior[k].mu = params.feedbacklag_secs.value - params.feedbacklag_secs_min;
        prior[k].sd = 2.0*prior[k].mu;
        prior[k].varclass = MatrixManip::Prior::VarClass::loglinear;
        //updateCovarMasks_StateVar( j, params.feedback_lag_secs.covarMask );
        paramindex.feedbackLagSecs = j;
        paramPriorList[j++] = k++;
    } else
        paramindex.feedbackLagSecs = -1;
    
    // Scaling coefficient
    prior[k].mu = params.ascaling.value; prior[k].sd = prior[k].mu/10;
    prior[k].varclass = MatrixManip::Prior::VarClass::loglinear;
    //updateCovarMasks_StateVar( j, params.ascaling.covarMask );
    paramindex.ascaling = j;
    paramPriorList[j++] = k++;
    
    // Observation noise (invariant to sampling interval and ascaling)
    prior[k].mu = params.sigmay.value; prior[k].sd = prior[k].mu/10;
    prior[k].varclass = MatrixManip::Prior::VarClass::variance;
    //updateCovarMasks_StateVar( j, params.sigmay.covarMask );
    for ( int i = 0; i < n_obs; i++ ) {
        paramindex.ySigma = j;
        paramPriorList[j] = k;
        j++;
    }
    k++;
    
    // Check param count
    assert( n_params == j );
    
    // Set Priors count
    n_priors = k;
    
    // Perform parameter transformations
    for ( int i = 0; i < n_params; i++ ) {
        switch ( prior[i].varclass ) {
            case MatrixManip::Prior::VarClass::linear:
                // No transformation needed
                break;
            case MatrixManip::Prior::VarClass::loglinear:
                // Log-transform
                prior[i].sd /= prior[i].mu;
                prior[i].mu = log( prior[i].mu );
                break;
            case MatrixManip::Prior::VarClass::variance:
                if ( prior[i].mu == 0 ) {
                    prior[i].mu = NAN;
                    prior[i].sd = 0.0;
                } else {
                    // Variance to std dev...
                    prior[i].mu = sqrt(prior[i].mu);
                    //  ...then log-transform
                    prior[i].sd /= prior[i].mu;
                    prior[i].mu = log( prior[i].mu );
                }
                break;
            default:
                assert("KuramotoUKF::formPriors: Unexpected variable class.");
        }
    }
}
bool** KuramotoUKF::initCovarMask( int dim1, int dim2, Masktype masktype ) {
    // Initialise state covar mask
    bool** mask = new bool*[dim1];
    for(int i=0; i<dim1; i++) {
        mask[i] = new bool[dim2];
        for(int j=0; j<dim2; j++)
            mask[i][j] = false;
    }
    // Select mask type and populate matrix
    switch ( masktype ) {
        case Masktype::zero:
            // Initialise with all false mask
            for ( int i = 0; i < dim1; i++ )
                for ( int j = 0; j < dim2; j++ )
                    mask[i][j] = false;
            break;
        case Masktype::defaultmask:
        case Masktype::full:
            // Initialise with full mask
            for ( int i = 0; i < dim1; i++ )
                for ( int j = 0; j < dim2; j++ )
                    mask[i][j] = true;
            break;
        case Masktype::block:
        case Masktype::diagonal:
            // Initialise with diagonal mask
            for ( int i = 0; i < dim1; i++ )
                for ( int j = 0; j < dim2; j++ )
                    if ( i == j )
                        mask[i][j] = true;
                    else
                        mask[i][j] = false;
            break;
        default:
            assert( "KuramotoUKF::initCovarMask: Requested covariance mask-type not recognised." );
    }
    return mask;
}
void KuramotoUKF::updateCovarMasks_StateVar( int k, Masktype masktype ) {
    updateCovarMasks_StateVar( k, masktype, nullptr, 0 );
}
void KuramotoUKF::updateCovarMasks_StateVar( int k, Masktype masktype, int* blockindices, int blockcount ) {
    switch ( masktype ) {
        case Masktype::defaultmask:
            break;
        case Masktype::zero:
            for ( int i = 0; i < n_statevars; i++ ) {
                statePmask[i][k] = false;
                statePmask[k][i] = false;
            }
            for ( int i = 0; i < n_statevars; i++ )
                crossPmask[k][i] = false;
            break;
        case Masktype::diagonal:
            for ( int i = 0; i < n_statevars; i++ ) {
                if  ( i == k ) {
                    statePmask[i][i] = true;
                } else {
                    statePmask[i][k] = false;
                    statePmask[k][i] = false;
                }
            }
            for ( int i = 0; i < n_statevars; i++ )
                crossPmask[k][i] = false;
            break;
        case Masktype::full:
            for ( int i = 0; i < n_statevars; i++ ) {
                statePmask[i][k] = true;
                statePmask[k][i] = true;
            }
            for ( int i = 0; i < n_statevars; i++ )
                crossPmask[k][i] = true;
            break;
        case Masktype::block:
            statePmask[k][k] = true;
            for ( int i = 0; i < blockcount; i++ )
                for ( int j = 0; j < blockcount; j++ )
                    statePmask[ blockindices[i] ][ blockindices[j] ] = true;
            break;
        default:
            assert( "KuramotoUKF::updateCovarMask: Requested covariance mask-type not recognised." );
    }
}
bool KuramotoUKF::isMaskFull( bool** mask, int dim1, int dim2 ) {
    // Is covariance mask ::full ?
    for ( int i = 0; i < dim1; i++ ) {
        for ( int j = 0; j < dim2; j++ )
            if ( !mask[i][j] ) {
                return false;
            }
    }
    return true;
}
void KuramotoUKF::setPriors( int* paramPriorList, int n_priors, int n_params ) {
    // Not necessary for UKF operation, but is necessary for prior based optimisation
    KuramotoUKF::paramPriorList = paramPriorList;
    KuramotoUKF::n_priors = n_priors;
    KuramotoUKF::n_params = n_params;
}
KuramotoUKF::Prior* KuramotoUKF::getPriors() {
    return prior;
}
int* KuramotoUKF::getParamPriorList( ) {
    return paramPriorList;
}
void KuramotoUKF::printStateCond( stateConditions* statecond ) {
    cout << "Initial conditions report:" << endl << "x0:" << endl;
    printVector(statecond->x);
    cout << "P0:" << endl;
    printMatrix(statecond->P);
    cout << "stateNoise:" << endl;
    printMatrix(statecond->stateNoise);
    cout << "obsNoise:" << endl;
    printMatrix(statecond->obsNoise);
}
M1 KuramotoUKF::priorsToPriorVec( Prior* prior, int* n_priors ) {
    (*n_priors) = KuramotoUKF::n_priors;
    M1 x(*n_priors);
    for ( int k = 0; k < *n_priors; k++ )
        x[k] = prior[k].mu;
    return x;
}
KuramotoUKF::stateConditions* KuramotoUKF::unpackParamVec( M1 paramvec ) {
    
    // Initial state conditions structure
    stateConditions* statecond = new stateConditions;
    statecond->x          = MatrixManip::allocMatrix(n_statevars);
    statecond->P          = MatrixManip::allocMatrix(n_statevars,n_statevars);
    statecond->stateNoise = MatrixManip::allocMatrix(n_statevars,n_statevars);
    statecond->obsNoise   = MatrixManip::allocMatrix(n_obs,n_obs);
    
    // Initial states for generating data
    for ( int k = 0; k < n_statevars; k++ ) {
        statecond->x[k] = paramvec[k];
        statecond->P[k][k] = exp(paramvec[n_statevars+k]);
        statecond->P[k][k] = statecond->P[k][k]*statecond->P[k][k];
        if ( !isnan(paramvec[2*n_statevars+k]) ) {
            statecond->stateNoise[k][k] = exp(paramvec[2*n_statevars+k]);
            statecond->stateNoise[k][k] = statecond->stateNoise[k][k]*statecond->stateNoise[k][k];
        }
    }
    
    // Observation noise
    for ( int k = 0; k < n_obs; k++ ) {
        statecond->obsNoise[k][k] = exp(paramvec[paramindex.ySigma]);
        statecond->obsNoise[k][k] = statecond->obsNoise[k][k]*statecond->obsNoise[k][k];
    }
    
    // Sensory feedback
    if ( paramindex.feedbackStrength != -1 ) {
        statecond->feedbackStrength = exp(paramvec[paramindex.feedbackStrength]);
        statecond->feedbacklag_secs = feedbacklag_secs_min + exp(paramvec[paramindex.feedbackLagSecs]);
    } else {
        statecond->feedbackStrength = 0.0;
        statecond->feedbacklag_secs = feedbacklag_secs_min;
    }
    
    // Observation scaling
    statecond->ascaling = exp(paramvec[paramindex.ascaling]);
    
    return statecond;
}
M1 KuramotoUKF::priorVecToParamVec( M1 priorvec, int* paramPriorList ) {
    M1 paramvec = allocMatrix(n_params);
    for  (int k = 0; k < n_params; k++ )
        paramvec[k] = priorvec[paramPriorList[k]];
    return paramvec;
}
M1 KuramotoUKF::rmsy() {
    M1 rms = allocMatrix(n_obs);
    int M;
    for ( int k = 0; k < n_obs; k++ ) {
        M = 0;
        for ( int t = 0; t < N; t ++ ) {
            if ( !isnan(y[t][k]) ) {
                rms[k] += y[t][k]*y[t][k];
                M++;
            }
        }
        rms[k] /= M;
        rms[k] = sqrt(rms[k]);
    }
    return rms;
}
bool KuramotoUKF::isVariable( ModelVar var ) {
    return (var.mode == ParamMode::variable);
}
bool KuramotoUKF::isParameter( ModelVar var ) {
    return ( var.mode == ParamMode::parameter );
}
