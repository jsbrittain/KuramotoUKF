//
//  KurGradDescent.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 03/01/2017.
//  Copyright © 2017 John-Stuart Brittain. All rights reserved.
//

#ifndef KurGradDescent_hpp
#define KurGradDescent_hpp

#include <cstdio>

#include "GradDescent.hpp"
#include "KuramotoUKF.hpp"

class KurGradDescent : public GradDescent {
public:
    KuramotoUKF* kurf;
    
    KurGradDescent( KuramotoUKF* kurf, int paramcount, MatrixManip::Prior* prior );
    KurGradDescent( KuramotoUKF* kurf, int paramcount, MatrixManip::Prior* prior, Method method );
    ~KurGradDescent();
    datatype costFunction( datatype* state, MatrixManip::Prior* prior, int n ) override;
};

#endif /* KurGradDescent_hpp */
