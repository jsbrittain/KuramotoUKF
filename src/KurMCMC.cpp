//
//  KurMCMC.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 22/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#include "KurMCMC.hpp"

KurMCMC::KurMCMC() {
}
KurMCMC::~KurMCMC() {
}
void KurMCMC::runThread( KuramotoUKF::KurfParams kurfparams, KuramotoUKF::ModelParamsSimple modelparams, MetropolisChain::MetropolisChainParams mcmcparams, int paramcount, MatrixManip::Prior* prior, int* paramPriorList, int n_priors, std::string loadfile, std::string savedir, datatype** map, datatype** cost, int threadno ) {
    
    // Construct Kuromoto UKF
    KuramotoUKF* kurf = new KuramotoUKF( kurfparams );
    kurf->readMeasurementFile( loadfile, 1 );
    kurf->formPriors(modelparams);
    kurf->initialise();
    kurf->setInitialConditionsFromPriors();
    
    // Run MCMC
    KuramotoMCChain* mc = new KuramotoMCChain( kurf, mcmcparams );
    std::string filename = savedir + "/chain" + std::to_string(threadno) + ".txt";
    if ( threadno == 0 ) {
        mc->saveMatrixToTextFile( savedir + "/out_proposal_thread0.txt", mc->covarProposal, mc->statedim, mc->statedim );
    }
    mc->run();
    mc->saveChain(filename);
    
    // Return MAP cost and associated state-vector
    datatype* cost0 = new datatype;
    (*cost0) = mc->getMAPcost();
    (*cost) = cost0;
    datatype* map0 = mc->getMAPstate();
    (*map) = map0;
    
    // Clean-up
    delete mc;
    delete kurf;
}

datatype* KurMCMC::run( KuramotoUKF::KurfParams kurfparams, KuramotoUKF::ModelParamsSimple modelparams, MetropolisChain::MetropolisChainParams mcmcparams, int paramcount, MatrixManip::Prior* prior, int* paramPriorList, int n_priors, std::string loadfile, std::string savedir, int threadcount ) {
    
    // Run burn-in on single thread and use adjusted covar proposal for final set
    if ( mcmcparams.tuneProposal ) {
        // Check chainlength sufficient for tuning iterations
        int chainlength = mcmcparams.chainlength;
        int burnin = mcmcparams.burnin;
        mcmcparams.chainlength = mcmcparams.tuningiters;
        mcmcparams.burnin = 0;
        // Construct Kuromoto UKF
        KuramotoUKF* kurf = new KuramotoUKF( kurfparams );
        kurf->readMeasurementFile( loadfile, 1 );
        kurf->formPriors( modelparams );
        kurf->initialise();
        kurf->setInitialConditionsFromPriors();
        // Run MCMC
        KuramotoMCChain* mc = new KuramotoMCChain( kurf, mcmcparams );
        mc->run();
        // Overwrite proposal distribution
        for ( int k = 0; k < mcmcparams.statedim; k++ )
            mcmcparams.covarProposal[k][k] = mc->covarProposal[k][k];
        // Save proposal to file
        mc->saveMatrixToTextFile( savedir + "/out_tuned_proposal.txt", mc->covarProposal, mc->statedim, mc->statedim );
        // Turn off proposal tuning for main run
        mcmcparams.tuneProposal = false;
        mcmcparams.chainlength = chainlength;
        mcmcparams.burnin = burnin;
        // Take MAP as new start point (after initial burnin)
        //delete[] mcmcparams.state0;
        //mcmcparams.state0 = mc->getMAPstate();
        // Clean-up
        delete mc;
        delete kurf;
    }
    
    // Form pointer array for return values from threads
    datatype** maps = new datatype*[threadcount];
    for ( int k = 0; k < threadcount; k++ )
        maps[k] = NULL;
    datatype** costs = new datatype*[threadcount];
    for ( int k = 0; k < threadcount; k++ )
        costs[k] = NULL;
    
    // Initialise threads list
    std::thread* mcmcThreads = new std::thread[threadcount];
    
    // Call threads
    for ( int threadno = 0; threadno < threadcount; threadno++ )
        mcmcThreads[threadno] = std::thread( &KurMCMC::runThread, this, kurfparams, modelparams, mcmcparams, paramcount, prior, paramPriorList, n_priors, loadfile, savedir, &maps[threadno], &costs[threadno], threadno );
    
    // Join threads
    for ( int threadno = 0; threadno < threadcount; threadno++ )
        mcmcThreads[threadno].join();
    
    // Evaluate best solution to return
    datatype bestCost = (*costs[0]);
    int bestIndex = 0;
    for ( int k = 1; k < threadcount; k++ )
        if ( (*costs[k]) < bestCost ) {
            bestCost = (*costs[k]);
            bestIndex = k;
        }
    cout << "Best solution has cost = " << bestCost << endl;
    
    // Allocate new memory location for MAP estimate
    datatype* bestmap = MatrixManip::allocMatrix(n_priors);
    for ( int k = 0; k < n_priors; k++ )
        bestmap[k] = maps[bestIndex][k];
    //MatrixManip::printVector( bestmap, n_priors );
    
    // Tidy up
    delete[] mcmcThreads;
    delete[] maps;
    delete[] costs;
    
    // Return MAP estimate of parameters
    return bestmap;
}
