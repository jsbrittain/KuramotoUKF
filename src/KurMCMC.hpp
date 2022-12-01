//
//  KurMCMC.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 22/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#ifndef KurMCMC_hpp
#define KurMCMC_hpp

#include <thread>

#include "KuramotoUKF.hpp"
#include "MetropolisChain.hpp"
#include "KuramotoMC.hpp"

class KurMCMC {
public:
    KurMCMC();
    ~KurMCMC();
    void runThread( KuramotoUKF::KurfParams kurfparams, KuramotoUKF::ModelParamsSimple, MetropolisChain::MetropolisChainParams mcmcparams, int paramcount, MatrixManip::Prior* prior, int* paramPriorList, int n_priors, std::string loadfile, std::string savedir, datatype** map, datatype** cost, int threadno );
    datatype* run( KuramotoUKF::KurfParams kurfparams, KuramotoUKF::ModelParamsSimple, MetropolisChain::MetropolisChainParams mcmcparams, int paramcount, MatrixManip::Prior* prior, int* paramPriorList, int n_priors, std::string loadfile, std::string savedir, int threadcount );
};

#endif /* KurMCMC_hpp */
