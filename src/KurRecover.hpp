//
//  KurRecover.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 10/01/2017.
//  Copyright © 2017 John-Stuart Brittain. All rights reserved.
//

#ifndef KurRecover_hpp
#define KurRecover_hpp

#include <cstdio>

#include "KuramotoMC.hpp"
#include "KurPSOchains.hpp"
#include "KurMCMC.hpp"
#include "KurGradDescent.hpp"

class KurRecover {
public:
    struct Options {
        std::string loadfile;
        std::string savedir;
        std::string initstatefile;
        bool useRmsSigmay = true;
        int  threadcount = 1;
        bool do_pso  = true;
        GradDescent::Method grad_method = GradDescent::adam;    // Usually passthrough or adam
        bool do_hess = true;
        bool do_mcmc = true;
        bool verbose = 1;
        struct {
            int burnin      = 200;
            int chainlength = 1e4;
            bool tuneProposal = true;
            std::string proposalfile = "";
        } mcmc;
    };
    
    void generateData( KuramotoUKF::ModelParamsSimple modelparams, std::string savedir );
    void parameterRecovery( KuramotoUKF::ModelParamsSimple modelparams, Options options );
    void saveParams( std::string filename, M1 x, int n );
    M1 loadinitstates( std::string filename, int expected_states );
    M2 loadproposal( std::string filename, int n );
    void saveproposal( std::string filename, M2 P, int n );
};

#endif /* KurRecover_hpp */
