//
//  MatrixManip.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 14/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#ifndef MatrixManip_hpp
#define MatrixManip_hpp

#include <iostream>
#include <iomanip>
#include <fstream>

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <limits>

using namespace std;

using datatype=double;

using M1=std::vector<datatype>;
using M2=std::vector<M1>;
using M3=std::vector<M2>;
using M4=std::vector<M3>;

static const datatype pi = 3.14159265359;
static const datatype twopi = 2*pi;
static const datatype logtwopi = log(twopi);

template<class T>
T getmax(T a, T b) {
  return (a)>(b)?(a):(b);
}

template<class T>
T getmin(T a, T b) {
  return ((a)<(b)?(a):(b));
}

class MatrixManip {
private:
    datatype logtwopi;
    
    int persistentLogLikeliMVNdim = -1;
    M2 persistentLogLikeliMVNe;
    M2 persistentLogLikeliMVNeT;
    M2 persistentLogLikeliMVNiS;
    M2 persistentLogLikeliMVNeiS;
    M2 persistentLogLikeliMVNMahalanobisSqr;
public:
    
    struct Prior {
        datatype mu = 0.0;
        datatype sd = 0.0;
        enum VarClass { linear, loglinear, variance } varclass;
    };
    
    MatrixManip();
    ~MatrixManip();
    static void cholesky(M2 A, M2& L);
    static M2 mateye( int n );
    static M2 mateye( int n, datatype value );
    static void mattranspose( M2 A, M2 D );
    static void matmult( M2 A, M2 B, M2 D );
    static void matmult( M2 A, M1 B, M1 D );
    static void matmultbyscalar( M2 A, datatype scale, M2& D );
    static void matdivr( M2 A, M2 B, M2 D );
    static void matinvByHand2D( M2 A, M2 D );
    static void matinv( M2 A, int n, M2& D );
    static void matinvPD( M2 A, M2 D );
    static void matinvDiag( M2 A, M2 D );
    static void matinvL( M2 A, M2 D );
    static void matadd( M2 A, M2 B, M2& D );
    static void matsub( M2 A, M2 B, M2& D );
    static void printVector( M1 x );
    static void printMatrix( M1 X );
    static void printMatrix( M2 X );
    static void outerproduct( M1 x, M2& D );
    static void outerproduct( M2 x, M2& D );
    static void outerproduct( M1 x, M1 y, M2& D );
    static M1 allocVector( int dim );
    static M1 allocMatrix( int dim );
    static M2 allocMatrix( int dim1, int dim2 );
    static M3 allocMatrix( int dim1, int dim2, int dim3 );
    void uniformrand( int dim1, int dim2, M2& u );
    void boxmuller( M2 u, int dim, M1& z );
    void mvnrand( M1 mu, int dim, M2 Sigma, datatype& x );
    void mvnrand( M1 mu, int dim, M2 Sigma, M1& x );
    datatype randn( datatype mu, datatype sd );
    datatype randUniform( );
    datatype logRandUniform( );
    datatype logdetPD( M2 X, int n );
    datatype logLikeliMVN( M1 state, int statedim, M1 mu, M2 Sigma );
    datatype logLikeliMVNpersistent( M1 x, int n, M1 mu, M2 Sigma );
    void testRandomNumberGenerators( const std::string filestem );
    void writeMatrixToFile( const std::string filename, M2 x );
    static datatype normLikeli( datatype x, datatype mu, datatype sd );
    static datatype logLikeliPriors( M1 state, std::vector<Prior> prior );
    static void printPriors( std::vector<Prior> prior );
    static datatype matrms( M1 x );
    static void saveVectorToTextFile( std::string filename, M1 D );
    static void saveMatrixToTextFile( std::string filename, M1 D );
    static M1 loadVectorFromTextFile( std::string filename, int expected_dim );
    static void saveMatrixToTextFile( std::string filename, M2 D );
    static M2 loadMatrixFromTextFile( std::string filename, int expected_dim1, int expected_dim2 );
    static void saveMatrixToTextFile( std::string filename, M3 D );
};

#endif
