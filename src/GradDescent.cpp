//
//  GradDescent.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 03/01/2017.
//  Copyright © 2017 John-Stuart Brittain. All rights reserved.
//

#include "GradDescent.hpp"

GradDescent::GradDescent( int paramcount, MatrixManip::Prior* prior ) : paramcount(paramcount), prior(prior) {
    //
}
GradDescent::GradDescent( int paramcount, MatrixManip::Prior* prior, Method method ) : paramcount(paramcount), prior(prior), method(method) {
}
void GradDescent::setStartingPosition( datatype* x0 ) {
    x = MatrixManip::allocMatrix(paramcount);
    for ( int k = 0; k < paramcount; k++ )
        x[k] = x0[k];
}
void GradDescent::usePriorsForStartingPosition( ) {
    x = MatrixManip::allocMatrix(paramcount);
    for ( int k = 0; k < paramcount; k++ )
        x[k] = prior[k].mu;
}
void GradDescent::run() {
    
    // Check that method is not nullptr (passthrough)
    if ( method == passthrough ) {
        cost = costFunction(x, prior, paramcount);
        return;
    }
    
    // Allocate memory
    datatype *newx = MatrixManip::allocMatrix(paramcount);
    datatype *xp = MatrixManip::allocMatrix(paramcount);
    datatype *xm = MatrixManip::allocMatrix(paramcount);
    datatype *dCdx = MatrixManip::allocMatrix(paramcount);
    datatype *dx = MatrixManip::allocMatrix(paramcount);
    for ( int i = 0; i < paramcount; i++ ) dx[i] = 1e-5;//sdscaling*prior[i].sd;
    datatype newcost, lastcost, propchange;
    datatype beta1t = 0, newbeta1t = 0, beta2t = 0, newbeta2t = 0, mhat, vhat;
    datatype *v = MatrixManip::allocMatrix(paramcount);
    datatype *newv = MatrixManip::allocMatrix(paramcount);
    datatype *m = MatrixManip::allocMatrix(paramcount);
    datatype *newm = MatrixManip::allocMatrix(paramcount);
    datatype *cache = MatrixManip::allocMatrix(paramcount);
    datatype *newcache = MatrixManip::allocMatrix(paramcount);

    // Initial conditions
    cost = costFunction(x, prior, paramcount);
    if ( verbose ) {
        std::cout << "Iter 0, cost " << cost << std::endl;
    }
    
    // Iterate gradient descent until convergence
    int k = 1;
    while ( true ) {
        // Local derivatives
        for ( int i = 0; i < paramcount; i++ ) {
            for ( int j = 0; j < paramcount; j++ ) {
                if ( method == nesterov ) {
                    // Take derivative about projected point
                    xp[j] = x[j] - momentumcoeff*v[j]; xm[j] = x[j] - momentumcoeff*v[j];
                } else {
                    xp[j] = x[j]; xm[j] = x[j];
                }
            }
            xp[i] += dx[i]; xm[i] -= dx[i];
            dCdx[i] = ( costFunction(xp, prior, paramcount) - costFunction(xm, prior, paramcount) )/(2.0*dx[i]);
        }
        
        // Update position and cost
        switch ( method ) {
            case vanilla:
                for ( int i = 0; i < paramcount; i++ )
                    newx[i] = x[i] - learningrate*dCdx[i];
                break;
            case momentum:
            case nesterov:
                for ( int i = 0; i < paramcount; i++ ) {
                    newv[i] = momentumcoeff*v[i] + learningrate*dCdx[i];
                    newx[i] = x[i] - newv[i];
                }
                break;
            case adagrad:
                for ( int i = 0; i < paramcount; i++ ) {
                    newcache[i] = cache[i] + dCdx[i]*dCdx[i];
                    newx[i] = x[i] - learningrate*dCdx[i]/(sqrt(newcache[i]) + eps);
                }
                break;
            case RMSprop:
                for ( int i = 0; i < paramcount; i++ ) {
                    // RMSProp
                    newcache[i] = decay_rate*cache[i] + (1-decay_rate)*dCdx[i]*dCdx[i];
                    newx[i] = x[i] - learningrate*dCdx[i]/(sqrt(newcache[i]) + eps);
                }
                break;
            case adam:      // Adaptive Moment Estimation (Adam)
                newbeta1t = beta1t*beta1;
                newbeta2t = beta2t*beta2;
                for ( int i = 0; i < paramcount; i++ ) {
                    newm[i] = beta1*m[i] + (1-beta1)*dCdx[i];
                    newv[i] = beta2*v[i] + (1-beta2)*dCdx[i]*dCdx[i];
                    mhat = newm[i]/(1-newbeta1t);
                    vhat = newv[i]/(1-newbeta2t);
                    newx[i] = x[i] - learningrate*mhat/(sqrt(vhat)+epsilon);
                }
                break;
            case passthrough:
                assert(method!=passthrough);
                break;
        }
        newcost = costFunction(newx, prior, paramcount);
        
        // Prevent divergence
        lastcost = cost;
        if (( !std::isnan(newcost) ) && (newcost < cost)) {
            for ( int i = 0; i < paramcount; i++ ) {
                x[i] = newx[i];
                v[i] = newv[i];
                m[i] = newm[i];
                beta1t = newbeta1t;
                beta2t = newbeta2t;
                cache[i] = newcache[i];
            }
            cost = newcost;
            propchange = (lastcost-cost)/std::abs(lastcost);
            if ( propchange < conv_th )
                learningrate /= 10;     // Continue to decrease step size to end algorithm
            else
                if ( method != adam )
                    learningrate *= 10; // Increase step when convergent
        } else {
            learningrate /= 10;         // Decrease step when divergent
            propchange = 0.0;
        }
        
        // Verbose
        if ( verbose ) {
            std::cout << "Iter " << k << ", alpha (log10) " << log10(learningrate) << "[" << log10(learningrate_th) << "], cost " << cost << " (" << newcost << "), change (log10) " << log10(propchange) << " [" << log10(conv_th) << "]" << std::endl;
        }
        
        // Check for convergence
        if ( (propchange < conv_th) && (learningrate < learningrate_th) )
            break;
        k++;
    }
    
    // Tidy up
    MatrixManip::deallocMatrix( &newx, paramcount );
    MatrixManip::deallocMatrix( &xp,   paramcount );
    MatrixManip::deallocMatrix( &xm,   paramcount );
    MatrixManip::deallocMatrix( &dCdx, paramcount );
    MatrixManip::deallocMatrix( &dx,   paramcount );
}
datatype GradDescent::getCost() {
    return cost;
}
datatype* GradDescent::getPos() {
    return x;
}
void GradDescent::setVerbose( bool value ) {
    verbose = value;
}



void GradDescent::calcHessian() {
    MatrixManip::deallocMatrix( &hess, paramcount, paramcount );
    hess = MatrixManip::allocMatrix( paramcount, paramcount );
    datatype* xpp = MatrixManip::allocMatrix( paramcount );
    datatype* xmm = MatrixManip::allocMatrix( paramcount );
    datatype* xmp = MatrixManip::allocMatrix( paramcount );
    datatype* xpm = MatrixManip::allocMatrix( paramcount );
    datatype *dx  = MatrixManip::allocMatrix( paramcount );
    for ( int i = 0; i < paramcount; i++ ) dx[i] = 1e-5;//sdscaling*prior[i].sd;
    int progress = 0, newprogress = 0, iter = 0, totaliter = paramcount + paramcount*(paramcount-1)/2;
    
    std::cout << "Computing Hessian..." << std::endl << "<" << std::flush;
    for ( int i = 0; i < paramcount; i++ ) {
        for ( int j = i; j < paramcount; j++ ) {
            
            // Display progress
            iter++;
            newprogress = 100*iter/(totaliter);
            if  ( newprogress > progress ) {
                //std::cout << " [" << i << ", " << j << "] (" << paramcount << ")" << std::endl;
                if (( progress / 10 ) > (newprogress / 10))
                    std::cout << newprogress/10 << std::flush;
                else
                    std::cout << '.' << std::flush;
                progress = newprogress;
            }
            
            // Determine interrogation points
            for ( int k = 0; k < paramcount; k++ ) {
                xpp[k] = x[k]; xpm[k] = x[k];
                xmp[k] = x[k]; xmm[k] = x[k];
            }
            xpp[i] += dx[i]; xpp[j] += dx[j];
            xpm[i] += dx[i]; xpm[j] -= dx[j];
            xmp[i] -= dx[i]; xmp[j] += dx[j];
            xmm[i] -= dx[i]; xmm[j] -= dx[j];
            
            // Evaluate Hessian
            hess[i][j] = ( costFunction(xpp, prior, paramcount) - costFunction(xpm, prior, paramcount) - costFunction(xmp, prior, paramcount) + costFunction(xmm, prior, paramcount) )/(4.0*dx[i]*dx[j]);
            if ( i != j )
                hess[j][i] = hess[i][j];
        }
    }
    std::cout << ">" << std::endl;
    
    if ( verbose ) {
        MatrixManip::printMatrix( hess, paramcount, paramcount );
    }
    
    // Tidy-up
    MatrixManip::deallocMatrix( &xpp, paramcount );
    MatrixManip::deallocMatrix( &xmm, paramcount );
    MatrixManip::deallocMatrix( &xmp, paramcount );
    MatrixManip::deallocMatrix( &xpm, paramcount );
    MatrixManip::deallocMatrix( &dx,  paramcount );
}
void GradDescent::saveToFile( datatype** D, std::string filename ) {
    MatrixManip::saveMatrixToTextFile(filename, D, paramcount, paramcount);
}
datatype** GradDescent::getHessian() {
    return hess;
}
datatype** GradDescent::getRegularisedHessian() {
    datatype** cholhess = MatrixManip::allocMatrix(paramcount,paramcount);
    MatrixManip::cholesky(hess, paramcount, cholhess);
    datatype** cholhessT = MatrixManip::allocMatrix(paramcount,paramcount);
    MatrixManip::mattranspose(cholhess, paramcount, paramcount, cholhessT);
    datatype** reghess = MatrixManip::allocMatrix(paramcount,paramcount);
    MatrixManip::matmult(cholhess, paramcount, paramcount, cholhessT, paramcount, paramcount, reghess);
    MatrixManip::deallocMatrix(&cholhess, paramcount, paramcount);
    MatrixManip::deallocMatrix(&cholhessT, paramcount, paramcount);
    return reghess;
}
datatype** GradDescent::getInverseHessian() {
    datatype** invhess = MatrixManip::allocMatrix(paramcount,paramcount);
    MatrixManip::matinvPD(hess, paramcount, invhess);
    return invhess;
}
