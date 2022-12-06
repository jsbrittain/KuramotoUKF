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
#include <memory>

#include "MatrixManip.hpp"
#include "KuramotoUKF.hpp"
#include "KurPSO.hpp"

class KurPSOchains {
public:
    KurPSOchains();
    ~KurPSOchains();
    void runThread(
        KuramotoUKF::ModelParamsSimple modelparams,
        int paramcount,
        std::vector<KurPSO::Prior> prior,
        std::vector<int> paramPriorList,
        int n_priors,
        std::string loadfile,
        std::string savedir,
        std::shared_ptr<M1> map,
        std::shared_ptr<datatype> cost,
        int threadno );
    M1 run(
        KuramotoUKF::ModelParamsSimple modelparams,
        int paramcount,
        std::vector<KurPSO::Prior> prior,
        std::vector<int> paramPriorList,
        int n_priors,
        std::string loadfile,
        std::string savedir,
        int threadcount );
};

#endif /* KurPSOchains_hpp */
