//
//  KurGradDescent.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 03/01/2017.
//  Copyright © 2017 John-Stuart Brittain. All rights reserved.
//

#include "KurGradDescent.hpp"

KurGradDescent::KurGradDescent( KuramotoUKF* kurf, int paramcount, MatrixManip::Prior* prior ) : kurf(kurf), GradDescent( paramcount, prior ) {
    //
}

KurGradDescent::KurGradDescent( KuramotoUKF* kurf, int paramcount, MatrixManip::Prior* prior, Method method ) : kurf(kurf), GradDescent( paramcount, prior, method) {
    //
}

KurGradDescent::~KurGradDescent() {
    //
}

datatype KurGradDescent::costFunction( M1 state, MatrixManip::Prior* prior, int n ) {
    return kurf->negLogLikeliFcn( state, prior, n );
}
