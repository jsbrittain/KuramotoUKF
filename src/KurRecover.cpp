//
//  KurRecover.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 10/01/2017.
//  Copyright © 2017 John-Stuart Brittain. All rights reserved.
//

#include "KurRecover.hpp"

void KurRecover::generateData( KuramotoUKF::ModelParamsSimple modelparams, std::string savedir ) {
    
    // Set initial conditions
    KuramotoUKF::KurfParams kurfparams = modelparams.kurfparams;
    cout << "Constructor" << endl;
    KuramotoUKF kuramoto( kurfparams );
    
    // Form priors from simplified parameters structure
    kuramoto.formPriors( modelparams );
    KuramotoUKF::Prior *prior = kuramoto.getPriors();
    int* paramPriorList = kuramoto.getParamPriorList();
    
    // Get (and print) initial states for generating data
    int n_priors;
    datatype* priorvec = kuramoto.priorsToPriorVec( prior, &n_priors );
    cout << "Number of priors: " << n_priors << endl;
    
    // Select parameters randomly from prior distributions
    srand((int)time(NULL));
    if ( modelparams.gen_from_priors_mode == KuramotoUKF::ModelParamsSimple::GenFromPriorsMode::noisy ) {
        for ( int k = 0; k < n_priors; k++ )
            priorvec[k] = kuramoto.randn( prior[k].mu, prior[k].sd );
    }
    
    //// Generate data ////
    
    cout << "Initialise" << endl;
    kuramoto.initialise(modelparams.n_obs,6000);
    datatype* paramvec = kuramoto.priorVecToParamVec( priorvec, paramPriorList );
    KuramotoUKF::stateConditions* statecond = kuramoto.unpackParamVec(paramvec);
    
    kuramoto.setInitialConditions( statecond );
    kuramoto.generate(6000);
    
    kuramoto.saveObs(    savedir + "/gen_y.txt" );
    kuramoto.saveStates( savedir + "/gen_x.txt" );
    kuramoto.saveVectorToTextFile( savedir + "/out_params.txt", priorvec, n_priors );
    
    kuramoto.reset();
    
    /*datatype logLikeli = kuramoto.run();
    std::cout << "Run (no priors): " << logLikeli << std::endl;
    kuramoto.saveStates(      savedir + "/out_x.txt"       );
    kuramoto.saveStatesCovar( savedir + "/out_x_covar.txt" );
    kuramoto.saveObsPred(     savedir + "/out_ypred.txt"   );*/
    
    cout << " Neg log likeli (without priors): " << kuramoto.negLogLikeliFcn(priorvec) << endl;
    cout << " Neg log likeli (with priors):    " << kuramoto.negLogLikeliFcn(priorvec, prior, n_priors) << endl;
    MatrixManip::printVector(priorvec, n_priors);
}

void KurRecover::parameterRecovery( KuramotoUKF::ModelParamsSimple modelparams, Options options ) {
    
    /// Initialise KuramotoUKF model to get generated parameters (such as prior counts)
    
    // Construct Kuromoto UKF
    datatype* stateMAP = NULL;
    KuramotoUKF kuramoto( modelparams.kurfparams );
    kuramoto.readMeasurementFile( options.loadfile, 1 );
    if ( options.useRmsSigmay ) {
        modelparams.ascaling.value = kuramoto.rmsy()[0]/0.7071;           // Approx sine amplitude of data
        modelparams.sigmay.value = pow( 0.01*modelparams.ascaling.value, 2 );   // Obs noise SD defaults to 1% of sine amplitude approximation
        if (options.verbose) std::cout << "Adjusting ascaling by RMS to " << modelparams.ascaling.value << std::endl;
    }
    kuramoto.formPriors(modelparams);
    kuramoto.initialise();
    //kuramoto.setInitialConditionsFromPriors();
    KuramotoUKF::Prior *prior = kuramoto.getPriors();
    int* paramPriorList = kuramoto.getParamPriorList();
    int n_priors;
    datatype* priorvec = kuramoto.priorsToPriorVec( prior, &n_priors );
    datatype* paramvec = kuramoto.priorVecToParamVec( priorvec, paramPriorList );
    KuramotoUKF::stateConditions* statecond = kuramoto.unpackParamVec(paramvec);
    // Check if initial state file is specified
    if ( options.initstatefile.compare("") != 0 ) {
        if ( options.verbose )
            std::cout << "Loading initial state from file: " << options.initstatefile << std::endl;
        stateMAP = loadinitstates(options.initstatefile, n_priors);
        MatrixManip::printMatrix(stateMAP, n_priors);
        if ( options.verbose ) {
            std::cout << "Initialised KuramotoUKF gives neg-log-likeli (without priors): " << kuramoto.negLogLikeliFcn(stateMAP) << std::endl;
            std::cout << "Initialised KuramotoUKF gives neg-log-likeli (with priors): " << kuramoto.negLogLikeliFcn(stateMAP, prior, n_priors) << std::endl;
        }
    } else
        if ( options.verbose )
            std::cout << "Initialised KuramotoUKF gives neg-log-likeli (with priors): " << kuramoto.negLogLikeliFcn(priorvec, prior, n_priors) << std::endl;
    kuramoto.setInitialConditions( statecond );
    
    /// Particle Swarm Optimisation
    
    if ( options.do_pso ) {
        KurPSOchains kurpsochains;
        if ( stateMAP == NULL )
            stateMAP = MatrixManip::allocMatrix(n_priors);
        stateMAP = kurpsochains.run( modelparams, kuramoto.n_priors, prior, paramPriorList, n_priors, options.loadfile, options.savedir, options.threadcount );
        // Backup MAP
        MatrixManip::saveMatrixToTextFile(options.savedir + "/out_params_pso.txt", stateMAP, n_priors);
        
    }
    
    // Gradient descent
    
    KurGradDescent kurgraddescent( &kuramoto, n_priors, prior, options.grad_method );
    if ( stateMAP != NULL )
        kurgraddescent.setStartingPosition( stateMAP );
    else {
        stateMAP = MatrixManip::allocMatrix(n_priors);
        kurgraddescent.usePriorsForStartingPosition( );
    }
    kurgraddescent.setVerbose( options.verbose );
    kurgraddescent.run();
    stateMAP = kurgraddescent.getPos();
    datatype costMAP = kurgraddescent.getCost();
    if ( options.verbose )
        cout << "Gradient Descent cost: " << costMAP << " at:" << endl;
    MatrixManip::printMatrix(stateMAP, n_priors);
    // Backup MAP
    if ( options.grad_method != KurGradDescent::Method::passthrough )
        MatrixManip::saveMatrixToTextFile(options.savedir + "/out_params_grad.txt", stateMAP, n_priors);
    
    // Hessian
    
    datatype** invhess = NULL;
    if ( options.do_hess ) {
        kurgraddescent.calcHessian();
        kurgraddescent.saveToFile( kurgraddescent.getHessian(), options.savedir + "/hess.txt" );
        datatype** hess = MatrixManip::allocMatrix(n_priors,n_priors);
        hess = kurgraddescent.getHessian();
        
        // Use Laplace approximation for MCMC proposals
        invhess = kurgraddescent.getInverseHessian( );
        kurgraddescent.saveToFile( invhess, options.savedir + "/invhess.txt" );
        MatrixManip::printMatrix(invhess, n_priors, n_priors);
        
        // Check that derived covar is +ve definite, otherwise reject
        bool valid_hess = true;
        for ( int i = 0; i < n_priors; i++ )
            for ( int j = 0; j < n_priors; j++ )
                if ( isnan(hess[i][j]) || isinf(hess[i][j]) )
                    valid_hess = false;
        options.do_hess = false;
    }
    
    // MCMC interrogation
    
    datatype** covarProposal = MatrixManip::allocMatrix( n_priors, n_priors );
    for ( int i = 0; i < n_priors; i++ ) {
        if ( options.do_hess ) {
            // Get from Laplace approximation
            covarProposal[i][i] = invhess[i][i];
        } else {
            // Global scaling on proposal distributions
            covarProposal[i][i] = (1e-3)*exp(prior[i].sd);     // log -> sd
            covarProposal[i][i] = pow(10.0,round(log10(covarProposal[i][i])));  // Round to nearest power of 10 (for convenience)
            covarProposal[i][i] *= covarProposal[i][i];         //  sd -> var
        }
    }
    
    // Load proposal distribution from file (if specified)
    if ( options.mcmc.proposalfile.compare("") != 0 ) {
        MatrixManip::deallocMatrix(&covarProposal, n_priors, n_priors);
        covarProposal = loadproposal( options.mcmc.proposalfile, n_priors );
        for ( int k = 0; k < n_priors; k++ ) {
            //covarProposal[k][k] /= 10;
            cout << k << ": " << log10(covarProposal[k][k]) << endl;
        }
        //options.mcmc.tuneProposal = false;
    }
    
    MetropolisChain::MetropolisChainParams mcmcparams;
    mcmcparams.state0        = stateMAP;
    mcmcparams.covarProposal = covarProposal;       // Will also return "tuned" proposal
    mcmcparams.tuneProposal  = options.mcmc.tuneProposal;
    mcmcparams.tuningiters   = 100;
    mcmcparams.prior         = prior;
    mcmcparams.statedim      = n_priors;
    mcmcparams.randseed      = (int) time(NULL);
    mcmcparams.verbose       = options.verbose;
    mcmcparams.burnin        = options.mcmc.burnin;
    mcmcparams.chainlength   = options.mcmc.chainlength;
    
    if ( options.do_mcmc ) {
        // Divide chain by permitted
        mcmcparams.chainlength /= options.threadcount;
        // Setup MCMC
        KurMCMC kurmc = KurMCMC();
        modelparams.kurfparams.nodecount = kuramoto.nodecount;
        modelparams.kurfparams.n_statevars = kuramoto.n_statevars;
        stateMAP = kurmc.run( modelparams.kurfparams, modelparams, mcmcparams, kuramoto.n_params, prior, paramPriorList, n_priors, options.loadfile, options.savedir, options.threadcount );
    }
    
    // Run simulation on MAP estimate for final state-vectors
    datatype negloglikeli = kuramoto.negLogLikeliFcn( stateMAP, prior, n_priors );
    if ( options.verbose ) std::cout << "Writing estimated states to file for log-neg-likeli = " << negloglikeli << std::endl;
    saveParams( options.savedir + "/out_params.txt", stateMAP, n_priors );
    kuramoto.saveStates(      options.savedir + "/out_x.txt"      );
    kuramoto.saveStatesCovar( options.savedir + "/out_xSigma.txt" );
    kuramoto.saveObs(         options.savedir + "/out_y.txt"      );
    kuramoto.saveObsPred(     options.savedir + "/out_ypred.txt"  );
    kuramoto.saveObsCovar(    options.savedir + "/out_ySigma.txt" );
    
    // Clean-up
    delete[] paramPriorList;
    delete[] priorvec;
    delete[] paramvec;
}

void KurRecover::saveParams( std::string filename, datatype* x, int n ) {
    MatrixManip::saveVectorToTextFile(filename, x, n);
}

datatype* KurRecover::loadinitstates( std::string filename, int expected_states ) {
    datatype* x = MatrixManip::allocMatrix(expected_states);
    std::ifstream myfile(filename,std::ifstream::binary);
    if (myfile.is_open()) {
        // Count number of params in file
        int n = 0; datatype num;
        while (myfile >> num)
            n++;
        myfile.close();
        myfile.open(filename,std::ifstream::binary);
        // Check that this matches the number of expected params
        assert( n == expected_states );
        // Read params into state vector
        for ( int k = 0; k < expected_states; k++ )
            myfile >> x[k];
        myfile.close();
    } else
        assert( "Cannot open initial states file." );
    return x;
}
datatype** KurRecover::loadproposal( std::string filename, int n ) {
    datatype** x = MatrixManip::allocMatrix(n,n);
    std::ifstream myfile(filename,std::ifstream::binary);
    if (myfile.is_open()) {
        // Count number of params in file
        int k = 0; datatype num;
        while (myfile >> num)
            k++;
        myfile.close();
        myfile.open(filename,std::ifstream::binary);
        // Check that this matches the number of expected params
        if ( k != (n*n) ) {
            std::cout << "Error load proposal definition file: " << filename << std::endl << " expected " << n << " parameters, found " << (int)sqrt(k) << std::endl;
            assert( k == (n*n) );
        }
        // Read params into state vector
        for ( int i = 0; i < n; i++ )
            for ( int j = 0; j < n; j++ )
                myfile >> x[i][j];
        myfile.close();
    } else
        assert( "Cannot open proposal file for reading." );
    return x;
}
