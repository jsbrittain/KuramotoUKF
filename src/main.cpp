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
    modelparams.covarMask              = KuramotoUKF::Masktype::full;
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

    // Check input arguments
    if ( argc < 2 ) {
        std::cout
          << "Data folder should be specified:" << std::endl
          << "  ./kukf data_folder" << std::endl << std::endl
          << "If data_folder/y.txt exists and contains a time-series vector" << std::endl
          << "  then this file is read. Otherwise, surrogate data in generated" << std::endl
          << "  and placed into the data_folder." << std::endl;
        return 1;
    }
    
    // General settings
    std::string savedir = argv[1];
    
    // Setup model specification
    KuramotoUKF::ModelParamsSimple modelparams = specifyModelParams();
    KurRecover kurrec;
    KurRecover::Options recopts;
    
    // Check if data folder contain a valid data file...
    bool generate_surrogate_data = false;
    recopts.loadfile = savedir + "/y.txt";
    if (FILE *file = fopen(recopts.loadfile.c_str(), "r")) {
        // File exists, read in later
        fclose(file);
    } else {
        // ...otherwise, generate data to recover...
        generate_surrogate_data = true;
        kurrec.generateData( modelparams, savedir );
        recopts.loadfile = savedir + "/gen_y.txt";
        recopts.useRmsSigmay = false;
    }

    // Specify recovery params
    recopts.savedir       = savedir;
    recopts.initstatefile = savedir + "/out_params.txt";
    recopts.threadcount   = 3;
    recopts.do_pso        = true;
    recopts.grad_method   = GradDescent::adam;//passthrough;
    recopts.do_hess       = true;
    recopts.do_mcmc       = true;
    recopts.verbose       = true;
    recopts.mcmc.tuneProposal = false;
    recopts.mcmc.proposalfile = "";//savedir + "/out_tuned_proposal.txt";
    recopts.mcmc.burnin   = 1e3;
    recopts.mcmc.chainlength = 1e4;
    
    // Attempt full parameter recovery
    kurrec.parameterRecovery( modelparams, recopts );
    
    // Return
    return 0;
}
