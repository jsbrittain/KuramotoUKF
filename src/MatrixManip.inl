//
//  MatrixManip.inl
//  mcmc
//
//  Created by John-Stuart Brittain on 14/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

// Inline (".inl") function implementations.
//
// Templates cannot be compiled as part of the cpp file and linked since the calling types
//  are not known without being instantiated as part of the broader project source. Instead
//  template implementations are declared in the header file. Inline files keep definition
//  and implementation details separate. This file is "included" at the bottom of its
//  corresponding header file.
//

template <class T> void MatrixManip::cholesky(T **A, int n, T** L) {
    // (L)(L)^T = A
    // A = input matrix, n = (n x n) dimensionality, tindex = time index in to A, L = output matrix
    int i,j,k;
    T s;
    for (i = 0; i < n; i++) {
        for (j = 0; j < (i+1); j++) {
            s = 0.0;
            for (k = 0; k < j; k++)
                s += L[i][k] * L[j][k];
            L[i][j] = (i == j) ? sqrt(A[i][i] - s) : (1.0 / L[j][j] * (A[i][j] - s));
        }
    }
    //return;
}
template <class T> T** MatrixManip::mateye( int n ) {
    return mateye<T>( n, (T) 1.0 );
}
template <class T> T** MatrixManip::mateye( int n, T value ) {
    T** I = allocMatrix<T>( n, n );
    for ( int i = 0; i<n; i++ ) {
        for ( int j = 0; j<n; j++) {
            if ( i == j )
                I[i][j] = value;
            else
                I[i][j] = (T) 0;
        }
    }
    return I;
}
template <class T> void MatrixManip::mattranspose( T** A, int n1, int n2, T** D ) {
    // Matrix transpose
    for ( int i = 0; i < n1; i++ ) {
        for ( int j = 0; j < n2; j++ ) {
            D[j][i] = A[i][j];
        }
    }
}
template <class T> void MatrixManip::matmult( T** A, int na1, int na2, T* B, int nb1, T* D ) {
    // Matrix multipication --- slowest possible method!
    for ( int i = 0; i < na1; i++) {
        D[i] = 0;
        for ( int k = 0; k < na2; k++) {
            D[i] += A[i][k]*B[k];
        }
    }
}
template <class T> void MatrixManip::matmult( T** A, int na1, int na2, T** B, int nb1, int nb2, T** D ) {
    // Matrix multipication --- slowest possible method!
    for ( int i = 0; i < na1; i++) {
        for ( int j = 0; j < nb2; j++) {
            D[i][j] = 0.0;
            for ( int k = 0; k < na2; k++) {
                // D[i,j] += A[i,k]*B[k,j];
                //std::cout << i << " " << j << " " << k << " " << A[i*na2+k] << " " << B[k*nb2+j] << " " << A[i*na2+k]*B[k*nb2+j] << std::endl;
                D[i][j] += A[i][k]*B[k][j];
            }
            //std::cout << D[i*nb2+j] << std::endl;
        }
    }
}
template <class T> void MatrixManip::matmultbyscalar( T** A, int dim1, int dim2, T scale, T** D )
{
    int i,j;
    for ( i = 0; i < dim1; i++ ) {
        for ( j = 0; j < dim2; j++ )
            D[i][j] = scale*A[i][j];
    }
}
template <class T> void MatrixManip::matdivr( T** A, int na1, int na2, T** B, int nb1, int nb2, T** D ) {
    assert( 0 );
}
template <class T> void MatrixManip::matinvByHand2D( T** A, T** D ) {
    T invAbsDet;
    invAbsDet = 1.0/fabs(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
    D[0][0] =  invAbsDet*A[1][1];
    D[1][1] =  invAbsDet*A[0][0];
    D[0][1] = -invAbsDet*A[0][1];
    D[1][0] = -invAbsDet*A[1][0];
}
template <class T> void MatrixManip::matinv( T** A, int n, T** D ) {
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
template <class T> void MatrixManip::matinvPD( T** A, int n, T** D ) {
    // Inverse for Positive Definite Matrices
    
    T** L = allocMatrix(n, n);
    T** Linv = allocMatrix(n, n);
    T** LinvT = allocMatrix(n, n);
    
    cholesky( A, n, L );
    matinvL( L, n, Linv );
    mattranspose( Linv, n, n, LinvT );
    matmult( LinvT, n, n, Linv, n, n, D );
    
    deallocMatrix(&L, n, n);
    deallocMatrix(&Linv, n, n);
    deallocMatrix(&LinvT, n, n);
}
template <class T> void MatrixManip::matinvDiag( T** A, int n, T** D ) {
    // Inverse of diagonal matrix
    for ( int i = 0; i < n; i++ ) {
        for ( int j = 0; j < n; j++ )
            D[i][j] = 0.0;
    }
    for ( int k = 0; k < n; k++ )
        D[k][k] = 1.0/A[k][k];
}
template <class T> void MatrixManip::matinvL( T** L, int n, T** D ) {
    // Inverse of lower-diagonal matrix
    int i,j,k;
    
    for ( i = 0; i < n; i++ ) {
        D[i][i] = 1.0/L[i][i];
        for ( j = (i+1); j < n; j++ ) {
            for ( k = i; k < j; k++ )
                D[j][i] -= L[j][k]*D[k][i];
            D[j][i] /= L[j][j];
        }
    }
}
template <class T> void MatrixManip::matadd( T** A, T** B, int dim1, int dim2, T** D ) {
    // D = A - B
    for ( int i = 0; i < dim1; i++ ) {
        for ( int j = 0; j < dim2; j++ )
            D[i][j] = A[i][j] + B[i][j];
    }
}
template <class T> void MatrixManip::matsub( T** A, T** B, int dim1, int dim2, T** D ) {
    // D = A - B
    for ( int i = 0; i < dim1; i++ ) {
        for ( int j = 0; j < dim2; j++ )
            D[i][j] = A[i][j] - B[i][j];
    }
}
template <class T> void MatrixManip::printVector( T* x, int n ) {
    int k;
    for ( k = 0; k<n; k++ )
        std::cout << std::fixed << "    " << k << ": " << x[k] << std::endl;
}
template <class T> void MatrixManip::printVector( int* x, int n ) {
    int k;
    for ( k = 0; k<n; k++ )
        std::cout << std::fixed << "    " << k << ": " << x[k] << std::endl;
}
template <class T> void MatrixManip::printMatrix( T* X, int n ) {
    if ( X == NULL ) {
        std::cout << "NULL vector reference!" << std::endl;
        return;
    }
    printVector( X, n );
}
template <class T> void MatrixManip::printMatrix( T** X, int n1, int n2 ) {
    for ( int i = 0; i < n1; i++ ) {
        std::cout  << "    [ ";
        for ( int j = 0; j < n2; j++)
            std::cout << " " << X[i][j] << " ";
        std::cout << " ]" << std::endl;
    }
}
template <class T> void MatrixManip::outerproduct( T* x, int n, T** D ) {
    for ( int i = 0; i<n; i++) {
        for ( int j = 0; j<n; j++ ) {
            D[i][j] = x[i]*x[j];
        }
    }
}
template <class T> void MatrixManip::outerproduct( T** x, int n1, int n2, T** D ) {
    assert( n2 == 1 );
    for ( int i = 0; i<n1; i++) {
        for ( int j = 0; j<n1; j++ ) {
            D[i][j] = x[i][0]*x[j][0];
        }
    }
}
template <class T> void MatrixManip::outerproduct( T* x, int nx, T* y, int ny, T** D ) {
    for ( int i = 0; i<nx; i++) {
        for ( int j = 0; j<ny; j++ ) {
            D[i][j] = x[i]*y[j];
        }
    }
}
template <class T> T* MatrixManip::allocVector( int dim ) {
    return allocMatrix<T>( dim );
}
template <class T> T* MatrixManip::allocMatrix( int dim ) {
    T* M = new T[dim];
    for ( int i = 0; i < dim; i++ )
        M[i] = 0.0;
    return M;
}
template <class T> T** MatrixManip::allocMatrix( int dim1, int dim2 ) {
    // Allocate and initialise 2D matrix
    T** M = new T*[dim1];
    for ( int i = 0; i<dim1; i++ ) {
        M[i] = allocMatrix<T>( dim2 );
    }
    return M;
}
template <class T> T*** MatrixManip::allocMatrix( int dim1, int dim2, int dim3 ) {
    // Allocate and initialise 2D matrix
    T*** M = new T**[dim1];
    for ( int i = 0; i<dim1; i++ ) {
        M[i] = allocMatrix<T>( dim2, dim3 );
    }
    return M;
}
template <class T> void MatrixManip::deallocMatrix( T** M, int dim ) {
    if ( *M == NULL ) return;
    delete[] *M;
    *M = NULL;
}
template <class T> void MatrixManip::deallocMatrix( T*** M, int dim1, int dim2 ) {
    if ( *M == NULL ) return;
    for ( int k = 0; k < dim1; k++ )
        delete[] (*M)[k];
    delete[] *M;
    *M = NULL;
}
template <class T> void MatrixManip::deallocMatrix( T**** M, int dim1, int dim2, int dim3 ) {
    if ( *M == NULL ) return;
    for ( int i = 0; i < dim1; i++ ) {
        for ( int j = 0; j < dim2; j++ )
            delete[] (*M)[i][j];
        delete[] (*M)[i];
    }
    delete[] *M;
    *M = NULL;
}
template <class T> void MatrixManip::uniformrand( int dim1, int dim2, T **u ) {
    for ( int i = 0; i < dim1; i++ ) {
        for ( int j = 0; j < dim2; j++ )
            u[i][j] = ((T) rand()) /RAND_MAX;       // (0,1)
    }
}
template <class T> void MatrixManip::boxmuller( T** u, int dim, T* z ) {
    for ( int k = 0; k < dim; k++ ) {
        z[k] = sqrt(-2*log(u[0][k])) * cos( twopi*u[1][k] );
    }
}
template <class T> void MatrixManip::mvnrand( T* mu, int dim, T** Sigma, T* x ) {
    /*  1.  Find any real matrix A such that A AT = Σ.
     2.  Let z = (z1, …, zN)T be a vector whose components are N independent standard normal variates (which can be generated, for example, by using the Box–Muller transform).
     3.  Let x be μ + Az. This has the desired distribution due to the affine transformation property.
     https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
     */
    
    // Step 1: Find any real matrix A such that A AT = Σ
    T** L = allocMatrix( dim, dim );
    cholesky( Sigma, dim, L );
    
    // Step 2: Generate N i.i.d. standard normal variates
    T** u = allocMatrix( 2, dim );
    uniformrand( 2, dim, u );
    T* z = allocMatrix( dim );
    boxmuller( u, dim, z );
    
    // Step 3: x = μ + A.z
    matmult(L, dim, dim, z, dim, x );
    for ( int k = 0; k < dim; k++ )
        x[k] += mu[k];
    
    // Clean-up
    deallocMatrix( &L, dim, dim );
    deallocMatrix( &u, 2, dim );
    deallocMatrix( &z, dim );
}
template <class T> T MatrixManip::randn( T mu, T sd ) {
    T** Sigma = allocMatrix(1,1);
    Sigma[0][0] = sd*sd;
    T* x = new T[1];
    mvnrand( &mu, 1, Sigma, x );
    deallocMatrix(&Sigma, 1, 1);
    T y = x[0];
    delete[] x;
    return y;
}
template <class T> T MatrixManip::randUniform( ) {
    return ((T) (rand() % RAND_MAX))/RAND_MAX;
}
template <class T> T MatrixManip::logRandUniform( ) {
    return log( randUniform() );
}
template <class T> T MatrixManip::logdetPD( T** X, int n ) {
    // https://makarandtapaswi.wordpress.com/2011/07/08/cholesky-decomposition-for-matrix-inversion/
    
    // First, obtain Cholesky decomposition
    T** L = allocMatrix(n, n);
    cholesky( X, n, L );
    
    // Determininant of L from Cholesky decomposition is product of diagonal entries
    T logdet = 0.0;
    for ( int i = 0; i < n; i++ )
        logdet += log( L[i][i] );
    
    // Determinant of X is then the square of the determinant of the lower diagonal L
    logdet *= 2.0;
    
    // Cleanup
    deallocMatrix(&L, n, n);
    
    return logdet;
}
template <class T> T MatrixManip::logLikeliMVN( T* x, int n, T* mu, T** Sigma ) {
    // Returns (positive) log-likelihood for an observation given a multivariate normal distribution
    
    T** e = allocMatrix( n, 1 );
    T** eT = allocMatrix( 1, n );
    for ( int k = 0; k < n; k++ ) {
        e[k][0] = x[k] - mu[k];
        eT[0][k] = e[k][0];
    }
    T** iS = allocMatrix( n, n );
    matinvPD(Sigma, n, iS);
    T** eiS = allocMatrix(1, n);
    matmult(eT, 1, n, iS, n, n, eiS);
    T** MahalanobisSqr = allocMatrix(1, 1);
    matmult(eiS, 1, n, e, n, 1, MahalanobisSqr);
    
    T logLikeli = -0.5*( logdetPD(Sigma, n) + MahalanobisSqr[0][0] + n*logtwopi );
    
    // Cleanup
    deallocMatrix( &e, n, 1 );
    deallocMatrix( &eT, 1, n );
    deallocMatrix( &iS, n, n );
    deallocMatrix( &eiS, 1, n );
    deallocMatrix( &MahalanobisSqr, 1, 1 );
    
    return logLikeli;
}
template <class T> T MatrixManip::logLikeliMVNpersistent( T* x, int n, T* mu, T** Sigma ) {
    // Returns (positive) log-likelihood for an observation given a multivariate normal distribution
    assert( n > 0 );
    
    // If already allocated then leave alone, otherwise allocate
    if ( n > persistentLogLikeliMVNdim ) {
        deallocMatrix(&persistentLogLikeliMVNe, persistentLogLikeliMVNdim, 1 );
        deallocMatrix(&persistentLogLikeliMVNeT, 1, persistentLogLikeliMVNdim);
        deallocMatrix(&persistentLogLikeliMVNiS, persistentLogLikeliMVNdim, persistentLogLikeliMVNdim);
        deallocMatrix(&persistentLogLikeliMVNeiS, 1, persistentLogLikeliMVNdim);
        persistentLogLikeliMVNe = allocMatrix( n, 1 );
        persistentLogLikeliMVNeT = allocMatrix( 1, n );
        persistentLogLikeliMVNiS = allocMatrix( n, n );
        persistentLogLikeliMVNeiS = allocMatrix(1, n);
        if ( persistentLogLikeliMVNMahalanobisSqr == NULL ) persistentLogLikeliMVNMahalanobisSqr = allocMatrix(1, 1);
        persistentLogLikeliMVNdim = n;
    }
    
    for ( int k = 0; k < n; k++ ) {
        persistentLogLikeliMVNe[k][0] = x[k] - mu[k];
        persistentLogLikeliMVNeT[0][k] = persistentLogLikeliMVNe[k][0];
    }
    
    matinvPD(Sigma, n, persistentLogLikeliMVNiS);               // <-- Matrix inverse!!!
    matmult(persistentLogLikeliMVNeT, 1, n, persistentLogLikeliMVNiS, n, n, persistentLogLikeliMVNeiS);
    matmult(persistentLogLikeliMVNeiS, 1, n, persistentLogLikeliMVNe, n, 1, persistentLogLikeliMVNMahalanobisSqr);
    return -0.5*( logdetPD(Sigma, n) + persistentLogLikeliMVNMahalanobisSqr[0][0] + n*logtwopi );
}
template <class T> void MatrixManip::writeMatrixToFile( const std::string filename, T** x, int dim1, int dim2 ) {
    /*std::ofstream myfile(filename,std::ofstream::binary);
     for ( int i = 0; i < dim1; i++ ) {
     for ( int j = 0; j < dim2; j++ )
     myfile << x[i][j] << '\t';
     myfile << std::std::endl;
     }*/
    saveMatrixToTextFile(filename, x, dim1, dim2);
}
template <class T> void MatrixManip::testRandomNumberGenerators( const std::string filestem ) {
    int N = 1e5;
    
    std::cout << "\nWriting " << N << " samples of random noise ( uniform / mvn{diag} / mvn{full} ) to disk..." << std::endl;
    
    // Output uniform random numbers to file
    T** u = allocMatrix(N,1);
    uniformrand(N, 1, u);
    writeMatrixToFile( filestem + "/test_uniform.txt", u, N, 1 );
    
    // Output MVN (diagonal) to file
    int M = 3;
    T* mu = allocMatrix(M);
    mu[0] = 0.0; mu[1] = 1.0; mu[2] = -3.0;
    T** Sigma = mateye(M);
    T** X = allocMatrix(N, M);
    for ( int n = 0; n < N; n++ )
        mvnrand( mu, M, Sigma, X[n] );
    writeMatrixToFile( filestem + "/test_mvnDiag.txt", X, N, M );
    
    // Output MVN (covar structure) to file
    T** Sigma2 = allocMatrix(3,3);
    T Sigmak[3][3] = { {1, 0.2, 0.4}, {0.2, 0.5, -0.1}, {0.4, -0.1, 2} };
    for ( int i = 0; i < 3; i++ ) for ( int j = 0; j < 3; j++ ) Sigma2[i][j] = Sigmak[i][j];
    T** X2 = allocMatrix(N, 3);
    for ( int n = 0; n < N; n++ )
        mvnrand( mu, M, Sigma2, X2[n] );
    writeMatrixToFile( filestem + "/test_mvnFull.txt", X2, N, 3 );
    
    deallocMatrix(&u,N,1);
    deallocMatrix(&mu,M);
    deallocMatrix(&Sigma, M, M);
    deallocMatrix(&X,N,M);
    deallocMatrix(&Sigma2, 3,3);
    deallocMatrix(&X2, N, 3);
}
template <class T> T MatrixManip::normLikeli( T x, T mu, T sd ) {
    if ( sd == 0 )
        return 0;
    else {
        T sigma = sd*sd;
        return -0.5*( log(sigma) + (x-mu)*(x-mu)/sigma + log(twopi) );
    }
}
template <class T> T MatrixManip::logLikeliPriors( T* state, Prior* prior, int n ) {
    T logLikeli = 0.0;
    for ( int k = 0; k < n; k++ )
        logLikeli += normLikeli( state[k], prior[k].mu, prior[k].sd );
    return logLikeli;
}
template <class T> void MatrixManip::printPriors( Prior* prior, int n ) {
    std::cout << "Prior structure:" << std::endl;
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
template <class T> T MatrixManip::matrms( T* x, int n ) {
    T rms = 0;
    for ( int k = 0; k < n; k++ )
        rms += x[k]*x[k];
    rms = sqrt( rms/n );
    return rms;
}
template <class T> void MatrixManip::saveVectorToTextFile( std::string filename, T* D, int dim ) {
    saveMatrixToTextFile( filename, D, dim );
}
template <class T> void MatrixManip::saveMatrixToTextFile( std::string filename, T* D, int dim ) {
    std::ofstream myfile(filename,std::ofstream::binary);
    myfile << std::scientific << std::setprecision(std::numeric_limits<T>::max_digits10);
    if (myfile.is_open()) {
        // Write params into file
        for ( int i = 0; i < dim; i++ )
            myfile << D[i] << std::endl;
        myfile.close();
    } else
        assert( "Cannot open file for writing." );
}
template <class T> T* MatrixManip::loadVectorFromTextFile( std::string filename, int expected_dim ) {
    T* D = allocMatrix( expected_dim );
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
template <class T> void MatrixManip::saveMatrixToTextFile( std::string filename, T** D, int dim1, int dim2 ) {
    // Save matrix to text file: 2-dimensions
    std::ofstream myfile(filename,std::ofstream::binary);
    myfile << std::scientific << std::setprecision(std::numeric_limits<T>::max_digits10);
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
template <class T> void MatrixManip::saveMatrixToTextFile( std::string filename, T*** D, int dim1, int dim2, int dim3 ) {
    // Save matrix to text file: 3-dimensions
    std::ofstream myfile(filename,std::ofstream::binary);
    myfile << std::scientific << std::setprecision(std::numeric_limits<T>::max_digits10);
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
template <class T> T** MatrixManip::loadMatrixFromTextFile( std::string filename, int expected_dim1, int expected_dim2 ) {
    T** D = allocMatrix( expected_dim1, expected_dim2 );
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
