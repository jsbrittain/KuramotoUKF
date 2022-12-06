//
//  GradDescent.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 03/01/2017.
//  Copyright © 2017 John-Stuart Brittain. All rights reserved.
//

#include "GradDescent.hpp"

GradDescent::GradDescent( int paramcount, std::vector<MatrixManip::Prior> prior )
  : paramcount(paramcount), prior(prior) {
    //
}
GradDescent::GradDescent( int paramcount, std::vector<MatrixManip::Prior> prior, Method method )
  : paramcount(paramcount), prior(prior), method(method) {
}
void GradDescent::setStartingPosition( M1 x0 ) {
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
    M1 newx = MatrixManip::allocMatrix(paramcount);
    M1 xp = MatrixManip::allocMatrix(paramcount);
    M1 xm = MatrixManip::allocMatrix(paramcount);
    M1 dCdx = MatrixManip::allocMatrix(paramcount);
    M1 dx = MatrixManip::allocMatrix(paramcount);
    for ( int i = 0; i < paramcount; i++ ) dx[i] = 1e-5;//sdscaling*prior[i].sd;
    datatype newcost, lastcost, propchange;
    datatype beta1t = 0, newbeta1t = 0, beta2t = 0, newbeta2t = 0, mhat, vhat;
    M1 v = MatrixManip::allocMatrix(paramcount);
    M1 newv = MatrixManip::allocMatrix(paramcount);
    M1 m = MatrixManip::allocMatrix(paramcount);
    M1 newm = MatrixManip::allocMatrix(paramcount);
    M1 cache = MatrixManip::allocMatrix(paramcount);
    M1 newcache = MatrixManip::allocMatrix(paramcount);

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
}
datatype GradDescent::getCost() {
    return cost;
}
M1 GradDescent::getPos() {
    return x;
}
void GradDescent::setVerbose( bool value ) {
    verbose = value;
}
void GradDescent::calcHessian() {
    hess = MatrixManip::allocMatrix( paramcount, paramcount );
    M1 xpp = MatrixManip::allocMatrix( paramcount );
    M1 xmm = MatrixManip::allocMatrix( paramcount );
    M1 xmp = MatrixManip::allocMatrix( paramcount );
    M1 xpm = MatrixManip::allocMatrix( paramcount );
    M1 dx  = MatrixManip::allocMatrix( paramcount );
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
        MatrixManip::printMatrix( hess );
    }
}
void GradDescent::saveToFile( M2 D, std::string filename ) {
    MatrixManip::saveMatrixToTextFile(filename, D);
}
M2 GradDescent::getHessian() {
    return hess;
}
M2 GradDescent::getRegularisedHessian() {
    M2 cholhess = MatrixManip::allocMatrix(paramcount,paramcount);
    MatrixManip::cholesky(hess, cholhess);
    M2 cholhessT = MatrixManip::allocMatrix(paramcount,paramcount);
    MatrixManip::mattranspose(cholhess, cholhessT);
    M2 reghess = MatrixManip::allocMatrix(paramcount,paramcount);
    MatrixManip::matmult(cholhess, cholhessT, reghess);
    return reghess;
}
M2 GradDescent::getInverseHessian() {
    M2 invhess = MatrixManip::allocMatrix(paramcount,paramcount);
    MatrixManip::matinvPD(hess, invhess);
    return invhess;
}
