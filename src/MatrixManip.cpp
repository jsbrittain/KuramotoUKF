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
    deallocMatrix(&persistentLogLikeliMVNe, persistentLogLikeliMVNdim, 1 );
    deallocMatrix(&persistentLogLikeliMVNeT, 1, persistentLogLikeliMVNdim);
    deallocMatrix(&persistentLogLikeliMVNiS, persistentLogLikeliMVNdim, persistentLogLikeliMVNdim);
    deallocMatrix(&persistentLogLikeliMVNeiS, 1, persistentLogLikeliMVNdim);
}
