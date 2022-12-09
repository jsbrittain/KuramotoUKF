//
//  GradDescent.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 03/01/2017.
//  Copyright © 2017 John-Stuart Brittain. All rights reserved.
//

// Accelerated gradient descent follows:
//  sebastianruder.com/optimizing-gradient-descent/
//  cs231n.github.io/neural-networks-3/

#ifndef GradDescent_hpp
#define GradDescent_hpp

#include <iostream>
#include <cstdio>
#include "MatrixManip.hpp"

class GradDescent {
public:
    bool verbose = false;
    int paramcount;
    std::vector<MatrixManip::Prior> prior;
    M1 x;
    datatype cost = 0;
    datatype momentumcoeff = 0.9;
    datatype learningrate = 0.1;
    M2 hess;
    datatype sdscaling = 1e-4;
    datatype eps = 1e-6;
    datatype decay_rate = 0.99;
    datatype conv_th = 1e-4;
    datatype learningrate_th = 1e-20;
    enum Method { passthrough, vanilla, momentum, nesterov, adagrad, RMSprop, adam } method = adam;
    
    // Adam parameters
    datatype beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8;
    
    GradDescent( int paramcount, std::vector<MatrixManip::Prior> prior );
    GradDescent( int paramcount, std::vector<MatrixManip::Prior> prior, Method method );
    virtual datatype costFunction( M1 state, std::vector<MatrixManip::Prior> prior, int n ) = 0;
    void setStartingPosition( M1 x0 );
    void usePriorsForStartingPosition( );
    void run();
    datatype getCost();
    M1 getPos();
    void setVerbose( bool value );
    void calcHessian();
    void saveToFile( M2 D, std::string filename );
    M2 getHessian();
    M2 getRegularisedHessian();
    M2 getInverseHessian();
};

#endif /* GradDescent_hpp */
