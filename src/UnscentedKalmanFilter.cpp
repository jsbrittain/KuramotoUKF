//
//  UnscentedKalmanFilter.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 14/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#include "UnscentedKalmanFilter.hpp"

UnscentedKalmanFilter::UnscentedKalmanFilter() {
    for (int k=0; k<TRIGRES; k++) {
        costable[k] = cos(k*dph);
        sintable[k] = sin(k*dph);
    }
}
UnscentedKalmanFilter::~UnscentedKalmanFilter() {
}
void UnscentedKalmanFilter::reset( ) {
    // Initial state will still be defined at t=0, so just reset time
    t = 0;
}
void UnscentedKalmanFilter::genSigmaX( M1 x, M2 P, int dim, M2 sigmaX ) {
    
    // SigmaX are ordered: [L'][k]
    
    // Construct scaled sqrt-covariance
    cholesky( P, sigmasqrtP );
    matmultbyscalar( sigmasqrtP, sqrtPscaling, sigmasqrtP );
    
    // Traverse state variables
    for (int k = 0; k < dim; k++) {
        for ( int i = 0; i < n_sigmavecs; i ++ )
            sigmaX[i][k] = 0.0;
        sigmaX[0][k] = x[k];
        // Traverse sigma vectors (1..L)
        for (int i = 0; i < dim; i++) {
            sigmaX[i+1][k] = x[k] + sigmasqrtP[k][i];
            sigmaX[dim+i+1][k] = x[k] - sigmasqrtP[k][i];
        }
    }
}
void UnscentedKalmanFilter::genSigmaW( ) {
    sigmaWm[0] = clambda/( n_statevars + clambda );
    sigmaWc[0] = clambda/( n_statevars + clambda ) + ( 1.0-(calpha*calpha)+cbeta );
    for (int i=1; i<n_sigmavecs; i++) {
        sigmaWm[i] = 1.0/(2.0*( n_statevars + clambda ));
        sigmaWc[i] = sigmaWm[i];
    }
}
void UnscentedKalmanFilter::predict() {
    
    // Generate sigma vectors from current state
    genSigmaX( state[t], stateP[t], n_statevars, sigmaX );
    t++;
    
    // State transition on sigma vectors
    for ( int i = 0; i < n_sigmavecs; i++)
        stateTransitionFunction( sigmaX[i], sigmaXpred[i] );
    
    // Reform mean and covar from sigmaXpred
    for ( int k = 0; k < n_statevars; k++ ) {
        statepred[t][k] = 0;        // Reset to zero to permit consistent answer after time reset
        for ( int i = 0; i < n_sigmavecs; i++ ) {
            statepred[t][k] += sigmaWm[i]*sigmaXpred[i][k];
        }
    }
    // Mean subtract sigmaXpred
    for ( int i = 0; i < n_sigmavecs; i++ ) {
        for ( int k = 0; k<n_statevars; k++ ) {
            sigmaXpredMS[i][k] = sigmaXpred[i][k] - statepred[t][k];
        }
    }
    // Form covariance matrix
    for ( int i = 0; i < n_sigmavecs; i++ )
        outerproduct( sigmaXpredMS[i], sigmaPcovarX[i] );
    for ( int i = 0; i < n_statevars; i++ ) {
        for ( int j = 0; j < n_statevars; j++) {
            // Zero covar
            statePXXpred[t][i][j] = 0.0;
            // Only calculate required elements
            if ( statePmask[i][j] ) {
                // Weighted sum over sigma vectors
                for ( int k = 0; k < n_sigmavecs; k++ )
                    statePXXpred[t][i][j] += sigmaWc[k]*sigmaPcovarX[k][i][j];
                statePXXpred[t][i][j] += stateNoise[i][j];
            }
        }
    }
    
    // Measurement prediction
    for ( int k = 0; k < n_obs; k++ )
        ypred[t][k] = 0.0;
    for ( int i = 0; i < n_sigmavecs; i++ ) {
        observationFunction( sigmaXpred[i], sigmaY[i] );
        for ( int k = 0; k < n_obs; k++ )
            ypred[t][k] += sigmaWm[i]*sigmaY[i][k];
    }
}
void UnscentedKalmanFilter::update() {
    
    // Mean subtract sigmaYpred
    for ( int i = 0; i <  n_sigmavecs; i ++ )
        for ( int k = 0; k < n_obs; k++ )
            sigmaYpredMS[i][k] = sigmaY[i][k] - ypred[t][k];
    // Pyy
    for ( int i = 0; i < n_sigmavecs; i++ )
        outerproduct( sigmaYpredMS[i], sigmaPcovarY[i] );
    for ( int i = 0; i < n_obs; i++ ) {
        for ( int j = 0; j < n_obs; j++) {
            // Zero covar
            statePYYpred[t][i][j] = 0.0;
            // Only calculate required elements
            if ( obsPmask[i][j] ) {
                // Weighted sum over sigma vectors
                for ( int k = 0; k < n_sigmavecs; k++ )
                    statePYYpred[t][i][j] += sigmaWc[k]*sigmaPcovarY[k][i][j];
                statePYYpred[t][i][j] += obsNoise[i][j];
            }
        }
    }
    
    // Pxy
    for ( int i = 0; i <  n_sigmavecs; i ++ )
        outerproduct( sigmaXpredMS[i], sigmaYpredMS[i], sigmaPcovarXY[i] );
    for ( int i = 0; i < n_statevars; i++ ) {
        for ( int j = 0; j < n_obs; j++) {
            // Zero covar
            statePXYpred[t][i][j] = 0.0;
            // Only calculate required elements
            if ( crossPmask[i][j] ) {
                // Weighted sum over sigma vectors
                for ( int k = 0; k < n_sigmavecs; k++ )
                    statePXYpred[t][i][j] += sigmaWc[k]*sigmaPcovarXY[k][i][j];
            }
        }
    }
    
    // Gain matrix
    matinvPD( statePYYpred[t], PyyInv );
    matmult( statePXYpred[t], PyyInv, K );
    
    // Check for missing data
    if ( std::isnan(y[t][0]) ) {      // Need to check what to do if only one channel data missing
        // Update states without measurement data: x[k] = xp[k]
        for ( int k = 0; k < n_statevars; k++ ) {
            state[t][k] = statepred[t][k];
            // Propogate state covariance
            for ( int j = 0; j < n_statevars; j++ )
                stateP[t][k][j] = statePXXpred[t][k][j];
        }
        return;
    }
    
    // Update states: x[k] = xp[k] + K*(y-yp)
    for ( int k = 0; k < n_obs; k++ )
        ydiff[k] = y[t][k] - ypred[t][k];
    matmult( K, ydiff, state[t] );
    for ( int k = 0; k < n_statevars; k++ )
        state[t][k] += statepred[t][k];
    
    // Update state covariance
    mattranspose( K, KT );
    matmult( statePYYpred[t], KT, PyyKT );
    matmult( K, PyyKT, KPyyKT );
    matsub( statePXXpred[t], KPyyKT, stateP[t] );
}
void UnscentedKalmanFilter::initialise( ) {
    
    sigmasqrtP = allocMatrix( n_statevars, n_statevars );
    state = allocMatrix( N, n_statevars );
    statepred = allocMatrix( N, n_statevars );
    stateP = allocMatrix( N, n_statevars, n_statevars );
    statePXXpred = allocMatrix( N, n_statevars, n_statevars );
    statePYYpred = allocMatrix( N, n_obs, n_obs );
    statePXYpred = allocMatrix( N, n_statevars, n_obs );
    ypred = allocMatrix( N, n_obs );
    
    // Sigma vectors
    n_sigmavecs = 2*n_statevars+1;
    sigmaX = allocMatrix( n_sigmavecs, n_statevars );
    sigmaXpred = allocMatrix( n_sigmavecs, n_statevars );
    sigmaY = allocMatrix( n_sigmavecs, n_obs );
    sigmaYpredMS = allocMatrix( n_sigmavecs, n_obs );
    sigmaPcovarX = allocMatrix( n_sigmavecs, n_statevars, n_statevars );
    sigmaPcovarY = allocMatrix( n_sigmavecs, n_obs, n_obs );
    sigmaPcovarXY = allocMatrix( n_sigmavecs, n_statevars, n_obs );
    sigmaXpredMS = allocMatrix( n_sigmavecs, n_statevars );
    PyyInv = allocMatrix( n_obs, n_obs );
    ydiff = allocMatrix( n_obs );
    
    stateNoise = allocMatrix( n_statevars, n_statevars );
    obsNoise = allocMatrix( n_obs, n_obs );
    
    K = allocMatrix( n_statevars, n_obs );
    KT = allocMatrix( n_obs, n_statevars );
    PyyKT = allocMatrix( n_obs, n_statevars );
    KPyyKT = allocMatrix( n_statevars, n_statevars );
    L = allocMatrix( n_obs, n_obs );
    Linv = allocMatrix( n_obs, n_obs );
    LinvT = allocMatrix( n_obs, n_obs );
    
    t = 0;
    clambda = calpha*calpha*(n_statevars+ckappa)-n_statevars;
    sqrtPscaling = sqrt(n_statevars+clambda);
    sigmaWc = allocMatrix( n_sigmavecs );
    sigmaWm = allocMatrix( n_sigmavecs );
    
    genSigmaW();
}
void UnscentedKalmanFilter::checkCholesky() {
    datatype Xk[3][3] = { { 25, 15, -5 }, { 15, 18, 0 }, { -5, 0, 11 } };
    M2 X = allocMatrix(3,3);
    for ( int i = 0; i<3; i++) {
        for ( int j = 0; j<3; j++)
            X[i][j] = Xk[i][j];
    }
    M2 Y = allocMatrix(3, 3);
    datatype LXkcorrect[3][3] = { { 5, 0, 0}, {3, 3, 0}, {-1, 1, 3 } };
    M2 LXcorrect = allocMatrix(3,3);
    for ( int i = 0; i<3; i++) {
        for ( int j = 0; j<3; j++)
            LXcorrect[i][j] = LXkcorrect[i][j];
    }
    std::cout << "\nChecking cholesky decomposition using matrix X = " << std::endl;
    printMatrix( X );
    cholesky( X, Y );
    std::cout << "\nLower diagonal decomposition L = " << std::endl;
    printMatrix( Y );
    std::cout << "Should be:" << std::endl;
    printMatrix( LXcorrect );
    
    datatype X2k[4][4] = { { 18, 22, 54, 42 }, { 22, 70, 86, 62 }, { 54, 86, 174, 134 }, { 42, 62, 134, 106 } };
    M2 X2 = allocMatrix(4,4);
    for ( int i = 0; i<4; i++) for ( int j = 0; j<4; j++) X2[i][j] = X2k[i][j];
    
    M2 Y2 = allocMatrix( 4, 4 );
    datatype LX2kcorrect[4][4] = { { 4.24264, 0, 0, 0 }, { 5.18545, 6.56591, 0, 0 }, { 12.72792, 3.04604, 1.64974, 0 }, { 9.89949, 1.62455, 1.84971, 1.39262 } };
    M2 LX2correct = allocMatrix(4,4);
    for ( int i = 0; i<4; i++) for ( int j = 0; j<4; j++) LX2correct[i][j] = LX2kcorrect[i][j];
    
    std::cout << "\nChecking cholesky decomposition using matrix X = " << std::endl;
    printMatrix( X2 );
    cholesky( X2, Y2 );
    std::cout << "\nLower diagonal decomposition L = " << std::endl;
    printMatrix( Y2 );
    std::cout << "Should be (5 decimal places):" << std::endl;
    printMatrix( LX2correct );
    
    deallocMatrix(X2);
    deallocMatrix(LX2correct);
    deallocMatrix(Y2);
}
void UnscentedKalmanFilter::checkMatMultiply() {
    datatype Xk[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
    M2 X = allocMatrix(2,3);
    for ( int i = 0; i<2; i++) for ( int j = 0; j<3; j++) X[i][j] = Xk[i][j];
    
    datatype Yk[3][2] = { { 7, 8 }, { 9, 10 }, { 11, 12 } };
    M2 Y = allocMatrix(3,2);
    for ( int i = 0; i<3; i++) for ( int j = 0; j<2; j++) Y[i][j] = Yk[i][j];
    
    datatype XYkcorrect[2][2] = { { 58, 64 }, { 139, 154 } };
    M2 XYcorrect = allocMatrix(2,2);
    for ( int i = 0; i<2; i++) for ( int j = 0; j<2; j++) XYcorrect[i][j] = XYkcorrect[i][j];
    
    M2 XY = allocMatrix(2,2);
    matmult( X, Y, XY );
    
    std::cout << "\nChecking Matrix multipication..." << std::endl;
    std::cout << "X = " << std::endl;
    printMatrix( X );
    std::cout << "Y = " << std::endl;
    printMatrix( Y );
    std::cout << "Calculated XY = " << std::endl;
    printMatrix( XY );
    std::cout << "Correct XY = " << std::endl;
    printMatrix( XYcorrect );
}
void UnscentedKalmanFilter::checkMatInv() {
    datatype Xk[2][2] = { { 4, 7 }, { 2, 6 } };
    M2 X = allocMatrix(2,2);
    for ( int i = 0; i<2; i++) for ( int j = 0; j<2; j++) X[i][j] = Xk[i][j];
    
    datatype Ykcorrect[2][2] = { { 4, 7 }, { 2, 6 } };
    M2 Ycorrect = allocMatrix(2,2);
    for ( int i = 0; i<2; i++) for ( int j = 0; j<2; j++) Ycorrect[i][j] = Ykcorrect[i][j];
    
    M2 Y = allocMatrix( 2, 2 );
    
    matinvPD( X,  Y );
    std::cout << "\nChecking Matrix inverse [full]..." << std::endl << "X = " << std::endl;
    printMatrix( X );
    std::cout << "inv(X) = " << std::endl;
    printMatrix( Y );
    std::cout << "Shoud be:" << std::endl;
    printMatrix( Ycorrect );
    
    deallocMatrix( X );
    deallocMatrix( Y );
    deallocMatrix( Ycorrect );
    
    datatype X2k[3][3] = { { 2, 0, 0 }, { 0, 10, 0 }, { 0, 0, 0.5 } };
    M2 X2 = allocMatrix(3,3);
    for ( int i = 0; i<3; i++) for ( int j = 0; j<3; j++) X2[i][j] = X2k[i][j];
    
    datatype Y2kcorrect[3][3] = { { 0.5, 0, 0 }, { 0, 1.0/10.0, 0 }, { 0, 0, 2 } };
    M2 Y2correct = allocMatrix(3,3);
    for ( int i = 0; i<3; i++) for ( int j = 0; j<3; j++) Y2correct[i][j] = Y2kcorrect[i][j];
    
    M2 Y2 = allocMatrix( 3, 3 );
    
    matinvPD( X2, Y2 );
    std::cout << "\nChecking Matrix inverse [diagonal]..." << std::endl << "X = " << std::endl;
    printMatrix( X2 );
    std::cout << "inv(X) = " << std::endl;
    printMatrix( Y2 );
    std::cout << "Shoud be:" << std::endl;
    printMatrix( Y2correct );
}
void UnscentedKalmanFilter::mathCheck() {
    checkCholesky();
    checkMatMultiply();
    //checkMatInv();            // Routines optimised for Pyy inverse only (scalar!)
}
void UnscentedKalmanFilter::saveStates( const std::string filename ) {
    saveMatrixToTextFile(filename, state);
}
void UnscentedKalmanFilter::saveStatesCovar( const std::string filename ) {
    saveMatrixToTextFile(filename, stateP);
}
void UnscentedKalmanFilter::saveObsCovar( const std::string filename ) {
    saveMatrixToTextFile(filename, obsNoise);
}
void UnscentedKalmanFilter::generate( int steps ) {
    N = steps;
    M1 newstate = allocMatrix(n_statevars);
    M1 newy = allocMatrix(n_obs);
    for ( int t = 1; t < N; t++ ) {
        stateTransitionFunction(state[t-1], newstate);
        mvnrand( newstate, n_statevars, stateNoise, state[t] );
        
        observationFunction( state[t], newy );
        mvnrand( newy, n_obs, obsNoise, y[t] );
    }
}
void UnscentedKalmanFilter::saveObs( const std::string filename  ) {
    MatrixManip::saveMatrixToTextFile(filename, y);
}
void UnscentedKalmanFilter::saveObsPred( const std::string filename  ) {
    MatrixManip::saveMatrixToTextFile(filename, ypred);
}
