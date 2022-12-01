//
//  KurPSOchains.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 22/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#ifndef KurPSOchains_hpp
#define KurPSOchains_hpp

#include <iostream>
#include <thread>

#include "MatrixManip.hpp"
#include "KuramotoUKF.hpp"
#include "KurPSO.hpp"

class KurPSOchains {
public:
    KurPSOchains();
    ~KurPSOchains();
    void runThread( KuramotoUKF::ModelParamsSimple modelparams, int paramcount, KurPSO::Prior* prior, int* paramPriorList, int n_priors, std::string loadfile, std::string savedir, datatype** map, datatype** cost, int threadno );
    datatype* run( KuramotoUKF::ModelParamsSimple modelparams, int paramcount, KurPSO::Prior* prior, int* paramPriorList, int n_priors, std::string loadfile, std::string savedir, int threadcount );
};

#endif /* KurPSOchains_hpp */
