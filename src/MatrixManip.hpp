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

#define datatype double
#define getmax(a,b) ((a)>(b)?(a):(b))
#define getmin(a,b) ((a)<(b)?(a):(b))

static const datatype pi = 3.14159265359;
static const datatype twopi = 2*pi;
static const datatype logtwopi = log(2*pi);

class MatrixManip {
private:
    datatype logtwopi;
    
    int persistentLogLikeliMVNdim = -1;
    datatype** persistentLogLikeliMVNe = nullptr;
    datatype** persistentLogLikeliMVNeT = nullptr;
    datatype** persistentLogLikeliMVNiS = nullptr;
    datatype** persistentLogLikeliMVNeiS = nullptr;
    datatype** persistentLogLikeliMVNMahalanobisSqr = nullptr;
public:
    
    struct Prior {
        datatype mu = 0.0;
        datatype sd = 0.0;
        enum VarClass { linear, loglinear, variance } varclass;
    };
    
    MatrixManip();
    ~MatrixManip();
    template <class T=datatype> static void cholesky(T **A, int n, T** L);
    template <class T=datatype> static T** mateye( int n );
    template <class T=datatype> static T** mateye( int n, T value );
    template <class T=datatype> static void mattranspose( T** A, int n1, int n2, T** D );
    template <class T=datatype> static void matmult( T** A, int na1, int na2, T** B, int nb1, int nb2, T** D );
    template <class T=datatype> static void matmult( T** A, int na1, int na2, T* B, int nb1, T* D );
    template <class T=datatype> static void matmultbyscalar( T** A, int dim1, int dim2, T scale, T** D );
    template <class T=datatype> static void matdivr( T** A, int na1, int na2, T** B, int nb1, int nb2, T** D );
    template <class T=datatype> static void matinvByHand2D( T ** A, T** D );
    template <class T=datatype> static void matinv( T ** A, int n, T** D );
    template <class T=datatype> static void matinvPD( T** A, int n, T** D );
    template <class T=datatype> static void matinvDiag( T** A, int n, T** D );
    template <class T=datatype> static void matinvL( T** A, int n, T** D );
    template <class T=datatype> static void matadd( T** A, T** B, int dim1, int dim2, T** D );
    template <class T=datatype> static void matsub( T** A, T** B, int dim1, int dim2, T** D );
    template <class T=datatype> static void printVector( T* x, int n );
    template <class T=datatype> static void printVector( int* x, int n );
    template <class T=datatype> static void printMatrix( T* X, int n );
    template <class T=datatype> static void printMatrix( T** X, int n1, int n2 );
    template <class T=datatype> static void outerproduct( T* x, int n, T** D );
    template <class T=datatype> static void outerproduct( T** x, int n1, int n2, T** D );
    template <class T=datatype> static void outerproduct( T* x, int nx, T* y, int ny, T** D );
    template <class T=datatype> static T* allocVector( int dim );
    template <class T=datatype> static T* allocMatrix( int dim );
    template <class T=datatype> static T** allocMatrix( int dim1, int dim2 );
    template <class T=datatype> static T*** allocMatrix( int dim1, int dim2, int dim3 );
    template <class T=datatype> static void deallocMatrix( T** M, int dim );
    template <class T=datatype> static void deallocMatrix( T*** M, int dim1, int dim2 );
    template <class T=datatype> static void deallocMatrix( T**** M, int dim1, int dim2, int dim3 );
    template <class T=datatype> void uniformrand( int dim1, int dim2, T **u );
    template <class T=datatype> void boxmuller( T** u, int dim, T* z );
    template <class T=datatype> void mvnrand( T* mu, int dim, T** Sigma, T* x );
    template <class T=datatype> T randn( T mu, T sd );
    template <class T=datatype> T randUniform( );
    template <class T=datatype> T logRandUniform( );
    template <class T=datatype> T logdetPD( T** X, int n );
    template <class T=datatype> T logLikeliMVN( T* state, int statedim, T* mu, T** Sigma );
    template <class T=datatype> T logLikeliMVNpersistent( T* x, int n, T* mu, T** Sigma );
    template <class T=datatype> void testRandomNumberGenerators( const std::string filestem );
    template <class T=datatype> void writeMatrixToFile( const std::string filename, T** x, int dim1, int dim2 );
    template <class T=datatype> static T normLikeli( T x, T mu, T sd );
    template <class T=datatype> static T logLikeliPriors( T* state, Prior* prior, int n );
    template <class T=datatype> static void printPriors( Prior* prior, int n );
    template <class T=datatype> static T matrms( T* x, int n );
    template <class T=datatype> static void saveVectorToTextFile( std::string filename, T* D, int dim );
    template <class T=datatype> static void saveMatrixToTextFile( std::string filename, T* D, int dim );
    template <class T=datatype> static T* loadVectorFromTextFile( std::string filename, int expected_dim );
    template <class T=datatype> static void saveMatrixToTextFile( std::string filename, T** D, int dim1, int dim2 );
    template <class T=datatype> static T** loadMatrixFromTextFile( std::string filename, int expected_dim1, int expected_dim2 );
    template <class T=datatype> static void saveMatrixToTextFile( std::string filename, T*** D, int dim1, int dim2, int dim3 );
};

#include "MatrixManip.inl"

#endif /* MatrixManip_hpp */
