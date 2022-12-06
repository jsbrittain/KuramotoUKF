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
#include <vector>
#include <memory>

#include "KuramotoUKF.hpp"
#include "MetropolisChain.hpp"
#include "KuramotoMC.hpp"

class KurMCMC {
public:
    KurMCMC();
    ~KurMCMC();
    void runThread(
        KuramotoUKF::KurfParams kurfparams,
        KuramotoUKF::ModelParamsSimple modelparams,
        MetropolisChain::MetropolisChainParams mcmcparams,
        int paramcount,
        std::vector<MatrixManip::Prior> prior,
        std::vector<int> paramPriorList,
        int n_priors,
        std::string loadfile,
        std::string savedir,
        std::shared_ptr<M1> map,
        std::shared_ptr<datatype> cost,
        int threadno );
    M1 run(
        KuramotoUKF::KurfParams kurfparams,
        KuramotoUKF::ModelParamsSimple modelparams,
        MetropolisChain::MetropolisChainParams mcmcparams,
        int paramcount,
        std::vector<MatrixManip::Prior> prior,
        std::vector<int> paramPriorList,
        int n_priors,
        std::string loadfile,
        std::string savedir,
        int threadcount );
};

#endif /* KurMCMC_hpp */
