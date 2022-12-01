//
//  main.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 13/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#include "KurRecover.hpp"

KuramotoUKF::ModelParamsSimple specifyModelParams( ) {
    
    // Model specification
    KuramotoUKF::ModelParamsSimple modelparams;
    
    // Kuramoto UKF parameters for initialisation
    modelparams.kurfparams.alpha       = 1e-3;
    modelparams.kurfparams.beta        = 2.0;
    modelparams.kurfparams.kappa       = 0.0;
    modelparams.kurfparams.dt          = 0.010;
    
    modelparams.nodecount              = 3;
    modelparams.n_obs                  = 1;
    modelparams.covarMask              = KuramotoUKF::Masktype::diagonal;
    modelparams.feedbacklag_secs_min   = 0.050;
    modelparams.feedback_obs_chan      = 0;
    
    modelparams.theta = {
        .mode=KuramotoUKF::ParamMode::variable,
        .covarMask=KuramotoUKF::Masktype::block,
        .value=NAN, .P=pow(0.1,2), .P0=pow(1.0,2)
    };
    modelparams.omega = {
        .mode=KuramotoUKF::ParamMode::variable,
        .covarMask=KuramotoUKF::Masktype::block,
        .value=2*pi*5.0, .P=pow(0.1,2), .P0=pow(2*pi*1.0,2)
    };
    modelparams.connK = {
        .mode=KuramotoUKF::ParamMode::variable,
        .covarMask=KuramotoUKF::Masktype::block,
        .value=0.00001, .P=pow(1,2), .P0=pow(1,2)
        //.value=5.0, .P=0.0, .P0=pow(1,2)
    };
    modelparams.feedback_strength = {
        .mode=KuramotoUKF::ParamMode::none,//parameter,
        .covarMask=KuramotoUKF::Masktype::defaultmask,
        .value=0.1, .P=pow(0.1,2), .P0=pow(1.0,2)
    };
    modelparams.feedbacklag_secs = {
        .mode=KuramotoUKF::ParamMode::none,//fixed,
        .covarMask=KuramotoUKF::Masktype::defaultmask,
        .value=0.051, .P=pow(0.1,2), .P0=pow(1.0,2)
    };
    
    modelparams.gen_from_priors_mode = KuramotoUKF::ModelParamsSimple::means;
    return modelparams;
}

int main(int argc, const char * argv[]) {
    
    /*
       Things to do:
            1) Smoothing
            2) Missing data ( partial observations on y )
            3) Multiple datasets; hierarchical model
            4) Working covar masks
    */
    
    // General settings
    std::string savedir = "/Users/jsb/mcmc/results/surr2";
    
    // Setup model specification
    KuramotoUKF::ModelParamsSimple modelparams = specifyModelParams();
    KurRecover kurrec;
    KurRecover::Options recopts;
    
    // Generate data to recover
    bool generate_surrogate_data = true;
    if ( generate_surrogate_data ) {
        kurrec.generateData( modelparams, savedir );
        recopts.loadfile = savedir + "/gen_y.txt";
        recopts.useRmsSigmay = false;
    } else {
        recopts.loadfile = "/Users/jsb/mcmc/y.txt";
        //recopts.loadfile = "/Users/jsb/mcmc/y_gap.txt";
        //recopts.loadfile = savedir + "/gen_y.txt";
    }
    
    // Specify recovery params
    recopts.savedir       = savedir;
    recopts.initstatefile = savedir + "/out_params.txt";
    recopts.threadcount   = 3;
    recopts.do_pso        = false;
    recopts.grad_method   = GradDescent::passthrough;    // Usually passthrough or adam
    recopts.do_hess       = false;
    recopts.do_mcmc       = false;
    recopts.verbose       = true;
    recopts.mcmc.tuneProposal = false;
    recopts.mcmc.proposalfile = "";//savedir + "/out_tuned_proposal.txt";
    recopts.mcmc.burnin   = 0;
    recopts.mcmc.chainlength = 1e4;
    
    // Attempt full parameter recovery
    kurrec.parameterRecovery( modelparams, recopts );
    
    // Return
    return 0;
}
