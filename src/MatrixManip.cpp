//
//  MatrixManip.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 14/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#include "MatrixManip.hpp"

MatrixManip::MatrixManip() {
    logtwopi = log(twopi);
}
MatrixManip::~MatrixManip() {
}

void MatrixManip::cholesky(M2 A, M2& L) {
    // (L)(L)^T = A
    // A = input matrix, n = (n x n) dimensionality, L = output matrix
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < (i+1); j++) {
            datatype s = 0.0;
            for (size_t k = 0; k < j; k++)
                s += L[i][k] * L[j][k];
            L[i][j] = (i == j) ? sqrt(A[i][i] - s) : (1.0 / L[j][j] * (A[i][j] - s));
        }
    }
}
M2 MatrixManip::mateye( int n ) {
    return mateye( n, static_cast<datatype>(1.0) );
}
M2 MatrixManip::mateye( int n, datatype value ) {
    M2 I = allocMatrix(n, n);
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
            if ( i == j )
                I[i][j] = value;
            else
                I[i][j] = static_cast<datatype>(0);
        }
    }
    return I;
}
void MatrixManip::mattranspose(M2 A, M2 D) {
    // Matrix transpose
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[i].size(); j++) {
            D[j][i] = A[i][j];
        }
    }
}
void MatrixManip::matmult(M2 A, M1 B, M1 D ) {
    // Matrix multipication --- slowest possible method!
    for (size_t i = 0; i < A.size(); i++) {
        D[i] = 0;
        for (size_t k = 0; k < A[i].size(); k++) {
            D[i] += A[i][k]*B[k];
        }
    }
}
void MatrixManip::matmult( M2 A, M2 B, M2 D ) {
    // Matrix multipication --- slowest possible method!
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < B[0].size(); j++) {
            D[i][j] = 0.0;
            for (size_t k = 0; k < A[i].size(); k++) {
                D[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}
void MatrixManip::matmultbyscalar(M2 A, datatype scale, M2& D) {
    for (size_t i = 0; i<A.size(); i++) {
        for (size_t j = 0; j<A[i].size(); j++)
            D[i][j] = scale*A[i][j];
    }
}
void MatrixManip::matdivr( M2 A, M2 B, M2 D ) {
    assert( 0 );
}
void MatrixManip::matinvByHand2D( M2 A, M2 D ) {
    datatype invAbsDet;
    invAbsDet = 1.0/fabs(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
    D[0][0] =  invAbsDet*A[1][1];
    D[1][1] =  invAbsDet*A[0][0];
    D[0][1] = -invAbsDet*A[0][1];
    D[1][0] = -invAbsDet*A[1][0];
}
void MatrixManip::matinv( M2 A, int n, M2& D ) {
    // Invert square matrices
    switch ( n ) {
        case 1:     // Scalar
            D[0][0] = 1/A[0][0];
            break;
        case 2:     // 2D
            matinvByHand2D( A, D );
            break;
        default:
            assert("General matrix inverse not yet supported. Use matinvPD for positive definite matrices.");
    }
    
}
void MatrixManip::matinvPD( M2 A, M2 D ) {
    // Inverse for Positive Definite Matrices
    int n = A.size();
    
    M2 L = allocMatrix(n, n);
    M2 Linv = allocMatrix(n, n);
    M2 LinvT = allocMatrix(n, n);
    
    cholesky( A, L );
    matinvL( L, Linv );
    mattranspose( Linv, LinvT );
    matmult( LinvT, Linv, D );
}
void MatrixManip::matinvDiag( M2 A, M2 D ) {
    // Inverse of diagonal matrix
    size_t n = A.size();
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++)
            D[i][j] = 0.0;
    }
    for (size_t k = 0; k < n; k++)
        D[k][k] = 1.0/A[k][k];
}
void MatrixManip::matinvL( M2 L, M2 D ) {
    // Inverse of lower-diagonal matrix
    size_t n = L.size();
    for (size_t i = 0; i < n; i++) {
        D[i][i] = 1.0/L[i][i];
        for (size_t j = (i+1); j < n; j++) {
            for (size_t k = i; k < j; k++)
                D[j][i] -= L[j][k]*D[k][i];
            D[j][i] /= L[j][j];
        }
    }
}
void MatrixManip::matadd( M2 A, M2 B, M2& D ) {
    // D = A + B
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[i].size(); j++)
            D[i][j] = A[i][j] + B[i][j];
    }
}
 void MatrixManip::matsub( M2 A, M2 B, M2& D ) {
    // D = A - B
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[i].size(); j++)
            D[i][j] = A[i][j] - B[i][j];
    }
}
void MatrixManip::printVector( M1 x ) {
    for (size_t k=0; k<x.size(); k++)
        std::cout << std::fixed << "    " << k << ": " << x[k] << std::endl;
}
void MatrixManip::printMatrix( M1 X ) {
    if ( X.size() == 0 ) {
        std::cout << "nullptr vector reference!" << std::endl;
        return;
    }
    printVector( X );
}
 void MatrixManip::printMatrix( M2 X ) {
    for (size_t i=0; i<X.size(); i++) {
        std::cout << "    [ ";
        for (size_t j=0; j<X[i].size(); j++)
            std::cout << " " << X[i][j] << " ";
        std::cout << " ]" << std::endl;
    }
}
void MatrixManip::outerproduct( M1 x, M2& D ) {
    for (size_t i=0; i<x.size(); i++) {
        for (size_t j=0; j<x.size(); j++) {
            D[i][j] = x[i]*x[j];
        }
    }
}
void MatrixManip::outerproduct( M2 x, M2& D ) {
    size_t n = x.size();
    assert( x[0].size() == 1 );
    for(size_t i=0; i<n; i++) {
        for(size_t j=0; j<n; j++) {
            D[i][j] = x[i][0]*x[j][0];
        }
    }
}
 void MatrixManip::outerproduct( M1 x, M1 y, M2& D ) {
    for (size_t i=0; i<x.size(); i++) {
        for (size_t j=0; j<y.size(); j++) {
            D[i][j] = x[i]*y[j];
        }
    }
}
M1 MatrixManip::allocVector( int dim ) {
    return allocMatrix( dim );
}
std::vector<datatype> MatrixManip::allocMatrix( int dim ) {
    M1 M(dim,0);
    return M;
}
M2 MatrixManip::allocMatrix( int dim1, int dim2 ) {
    // Allocate and initialise 2D matrix
    M2 M(dim1,M1(dim2,0));
    return M;
}
M3 MatrixManip::allocMatrix( int dim1, int dim2, int dim3 ) {
    // Allocate and initialise 2D matrix
    M3 M(dim1,M2(dim2,M1(dim3,0)));
    return M;
}
void MatrixManip::uniformrand( int dim1, int dim2, M2& u ) {
    for ( int i = 0; i < dim1; i++ ) {
        for ( int j = 0; j < dim2; j++ )
            u[i][j] = (static_cast<datatype>(rand()))/RAND_MAX;       // (0,1)
    }
}
void MatrixManip::boxmuller( M2 u, int dim, M1& z ) {
    for ( int k = 0; k < dim; k++ ) {
        z[k] = sqrt(-2*log(u[0][k])) * cos( twopi*u[1][k] );
    }
}
void MatrixManip::mvnrand( M1 mu, int dim, M2 Sigma, datatype& x ) {
  M1 Mx(1);
  mvnrand(mu, dim, Sigma, Mx);
  x = Mx[0];
}
void MatrixManip::mvnrand( M1 mu, int dim, M2 Sigma, M1& x ) {
    /*  1.  Find any real matrix A such that A AT = Σ.
     2.  Let z = (z1, …, zN)T be a vector whose components are N independent standard normal variates (which can be generated, for example, by using the Box–Muller transform).
     3.  Let x be μ + Az. This has the desired distribution due to the affine transformation property.
     https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
     */
    
    // Step 1: Find any real matrix A such that A AT = Σ
    M2 L = allocMatrix( dim, dim );
    cholesky( Sigma, L );
    
    // Step 2: Generate N i.i.d. standard normal variates
    M2 u = allocMatrix( 2, dim );
    uniformrand( 2, dim, u );
    M1 z = allocMatrix( dim );
    boxmuller( u, dim, z );
    
    // Step 3: x = μ + A.z
    matmult(L, z, x );
    for ( int k = 0; k < dim; k++ )
        x[k] += mu[k];
}
datatype MatrixManip::randn( datatype mu, datatype sd ) {
    M2 Sigma = allocMatrix(1,1);
    Sigma[0][0] = sd*sd;
    M1 x;
    M1 Mmu(1);
    Mmu[0] = mu;
    mvnrand( Mmu, 1, Sigma, x );
    datatype y = x[0];
    return y;
}
datatype MatrixManip::randUniform( ) {
    return static_cast<datatype>(rand() % RAND_MAX)/RAND_MAX;
}
datatype MatrixManip::logRandUniform( ) {
    return log( randUniform() );
}
datatype MatrixManip::logdetPD( M2 X, int n ) {
    // https://makarandtapaswi.wordpress.com/2011/07/08/cholesky-decomposition-for-matrix-inversion/
    
    // First, obtain Cholesky decomposition
    M2 L = allocMatrix(n, n);
    cholesky( X, L );
    
    // Determininant of L from Cholesky decomposition is product of diagonal entries
    datatype logdet = 0.0;
    for ( int i = 0; i < n; i++ )
        logdet += log( L[i][i] );
    
    // Determinant of X is then the square of the determinant of the lower diagonal L
    logdet *= 2.0;
    
    return logdet;
}
datatype MatrixManip::logLikeliMVN( M1 x, int n, M1 mu, M2 Sigma ) {
    // Returns (positive) log-likelihood for an observation given a multivariate normal distribution
    
    M2 e = allocMatrix( n, 1 );
    M2 eT = allocMatrix( 1, n );
    for ( int k = 0; k < n; k++ ) {
        e[k][0] = x[k] - mu[k];
        eT[0][k] = e[k][0];
    }
    M2 iS = allocMatrix( n, n );
    matinvPD(Sigma, iS);
    M2 eiS = allocMatrix(1, n);
    matmult(eT, iS, eiS);
    M2 MahalanobisSqr = allocMatrix(1, 1);
    matmult(eiS, e, MahalanobisSqr);
    
    datatype logLikeli = -0.5*( logdetPD(Sigma, n) + MahalanobisSqr[0][0] + n*logtwopi );
    
    return logLikeli;
}
datatype MatrixManip::logLikeliMVNpersistent( M1 x, int n, M1 mu, M2 Sigma ) {
    // Returns (positive) log-likelihood for an observation given a multivariate normal distribution
    assert( n > 0 );
    
    // If already allocated then leave alone, otherwise allocate
    if ( n > persistentLogLikeliMVNdim ) {
        persistentLogLikeliMVNe = allocMatrix( n, 1 );
        persistentLogLikeliMVNeT = allocMatrix( 1, n );
        persistentLogLikeliMVNiS = allocMatrix( n, n );
        persistentLogLikeliMVNeiS = allocMatrix(1, n);
        if ( persistentLogLikeliMVNMahalanobisSqr.size()==0 ) persistentLogLikeliMVNMahalanobisSqr = allocMatrix(1, 1);
        persistentLogLikeliMVNdim = n;
    }
    
    for(int k = 0; k < n; k++) {
        persistentLogLikeliMVNe[k][0] = x[k] - mu[k];
        persistentLogLikeliMVNeT[0][k] = persistentLogLikeliMVNe[k][0];
    }
    
    matinvPD(Sigma, persistentLogLikeliMVNiS);               // <-- Matrix inverse!!!
    matmult(persistentLogLikeliMVNeT, persistentLogLikeliMVNiS, persistentLogLikeliMVNeiS);
    matmult(persistentLogLikeliMVNeiS, persistentLogLikeliMVNe, persistentLogLikeliMVNMahalanobisSqr);
    return -0.5*( logdetPD(Sigma, n) + persistentLogLikeliMVNMahalanobisSqr[0][0] + n*logtwopi );
}
void MatrixManip::writeMatrixToFile( const std::string filename, M2 x ) {
    saveMatrixToTextFile(filename, x);
}
void MatrixManip::testRandomNumberGenerators( const std::string filestem ) {
    int N = 1e5;
    
    std::cout << "\nWriting " << N << " samples of random noise ( uniform / mvn{diag} / mvn{full} ) to disk..." << std::endl;
    
    // Output uniform random numbers to file
    M2 u = allocMatrix(N,1);
    uniformrand(N, 1, u);
    writeMatrixToFile( filestem + "/test_uniform.txt", u );
    
    // Output MVN (diagonal) to file
    int M = 3;
    M1 mu = allocMatrix(M);
    mu[0] = 0.0; mu[1] = 1.0; mu[2] = -3.0;
    M2 Sigma = mateye(M);
    M2 X = allocMatrix(N, M);
    for ( int n = 0; n < N; n++ )
        mvnrand( mu, M, Sigma, X[n] );
    writeMatrixToFile( filestem + "/test_mvnDiag.txt", X );
    
    // Output MVN (covar structure) to file
    M2 Sigma2 = allocMatrix(3,3);
    datatype Sigmak[3][3] = { {1, 0.2, 0.4}, {0.2, 0.5, -0.1}, {0.4, -0.1, 2} };
    for ( int i = 0; i < 3; i++ ) for ( int j = 0; j < 3; j++ ) Sigma2[i][j] = Sigmak[i][j];
    M2 X2 = allocMatrix(N, 3);
    for ( int n = 0; n < N; n++ )
        mvnrand( mu, M, Sigma2, X2[n] );
    writeMatrixToFile( filestem + "/test_mvnFull.txt", X2 );
}
datatype MatrixManip::normLikeli( datatype x, datatype mu, datatype sd ) {
    if ( sd == 0 )
        return 0;
    else {
        datatype sigma = sd*sd;
        return -0.5*( log(sigma) + (x-mu)*(x-mu)/sigma + log(twopi) );
    }
}
datatype MatrixManip::logLikeliPriors( M1 state, std::vector<Prior> prior ) {
    assert(state.size() <= prior.size());
    datatype logLikeli = 0.0;
    for (size_t k = 0; k<prior.size(); k++)
        logLikeli += normLikeli( state[k], prior[k].mu, prior[k].sd );
    return logLikeli;
}
void MatrixManip::printPriors( std::vector<Prior> prior ) {
    std::cout << "Prior structure:" << std::endl;
    int n = prior.size();
    for( int k = 0; k < n; k++) {
        std::cout << k << ": ( " << prior[k].mu << ", " << prior[k].sd << " ) [";
        switch ( prior[k].varclass ) {
            case MatrixManip::Prior::VarClass::linear:
                std::cout << "linear";
                break;
            case MatrixManip::Prior::VarClass::loglinear:
                std::cout << "loglinear";
                break;
            case MatrixManip::Prior::VarClass::variance:
                std::cout << "variance";
                break;
            default:
                assert("MatrixManip::printPriors: Unexpected variable class.");
        }
        std::cout << "]" << std::endl;
    }
}
datatype MatrixManip::matrms( M1 x ) {
    int n = x.size();
    datatype rms = 0;
    for ( int k = 0; k < n; k++ )
        rms += x[k]*x[k];
    rms = sqrt( rms/n );
    return rms;
}
void MatrixManip::saveVectorToTextFile( std::string filename, M1 D ) {
    saveMatrixToTextFile( filename, D );
}
void MatrixManip::saveMatrixToTextFile( std::string filename, M1 D ) {
    int dim = D.size();
    std::ofstream myfile(filename,std::ofstream::binary);
    myfile << std::scientific << std::setprecision(std::numeric_limits<datatype>::max_digits10);
    if (myfile.is_open()) {
        // Write params into file
        for ( int i = 0; i < dim; i++ )
            myfile << D[i] << std::endl;
        myfile.close();
    } else
        assert( "Cannot open file for writing." );
}
M1 MatrixManip::loadVectorFromTextFile( std::string filename, int expected_dim ) {
    M1 D = allocMatrix( expected_dim );
    std::ifstream myfile(filename,std::ofstream::binary);
    if (myfile.is_open()) {
        // Read from file
        for ( int i = 0; i < expected_dim; i++ )
            myfile >> D[i];
        myfile.close();
    } else
        assert( "Cannot open file for reading." );
    return D;
}
void MatrixManip::saveMatrixToTextFile( std::string filename, M2 D ) {
    int dim1 = D.size(), dim2 = D[0].size();
    // Save matrix to text file: 2-dimensions
    std::ofstream myfile(filename,std::ofstream::binary);
    myfile << std::scientific << std::setprecision(std::numeric_limits<datatype>::max_digits10);
    if (myfile.is_open()) {
        // Write params into file
        for ( int i = 0; i < dim1; i++ ) {
            for ( int j = 0; j < dim2; j++ )
                myfile << "\t" << D[i][j];
            myfile << std::endl;
        }
        myfile.close();
    } else
        assert( "Cannot open file for writing." );
}
void MatrixManip::saveMatrixToTextFile( std::string filename, M3 D ) {
    // Save matrix to text file: 3-dimensions
    int dim1=D.size(), dim2=D[0].size(), dim3=D[0][0].size();
    std::ofstream myfile(filename,std::ofstream::binary);
    myfile << std::scientific << std::setprecision(std::numeric_limits<datatype>::max_digits10);
    if (myfile.is_open()) {
        // Write params into file
        for ( int i = 0; i < dim1; i++ ) {
            for ( int j = 0; j < dim2; j++ )
                for ( int k = 0; k < dim3; k++ )
                    myfile << "\t" << D[i][j][k];
            myfile << std::endl;
        }
        myfile.close();
    } else
        assert( "Cannot open file for writing." );
}
M2 MatrixManip::loadMatrixFromTextFile( std::string filename, int expected_dim1, int expected_dim2 ) {
    M2 D = allocMatrix( expected_dim1, expected_dim2 );
    std::ifstream myfile(filename,std::ifstream::binary);
    if (myfile.is_open()) {
        // Read from file
        for ( int i = 0; i < expected_dim1; i++ ) {
            for ( int j = 0; j < expected_dim2; j++ )
                myfile >> D[i][j];
        }
        myfile.close();
    } else
        assert( "Cannot open file for reading." );
    return D;
}
