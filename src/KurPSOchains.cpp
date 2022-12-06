//
//  KurPSOchains.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 22/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#include "KurPSOchains.hpp"

KurPSOchains::KurPSOchains() {
}
KurPSOchains::~KurPSOchains() {
}
void KurPSOchains::runThread(
    KuramotoUKF::ModelParamsSimple modelparams,
    int paramcount,
    std::vector<KurPSO::Prior> prior,
    std::vector<int> paramPriorList,
    int n_priors,
    std::string loadfile,
    std::string savedir,
    std::shared_ptr<M1> map,
    std::shared_ptr<datatype> cost,
    int threadno )
{ 
    // Construct Kuromoto UKF
    KuramotoUKF kurf( modelparams.kurfparams );
    kurf.readMeasurementFile(loadfile, 1);
    //kurf->setPriors( paramPriorList, n_priors, paramcount );
    kurf.formPriors( modelparams );
    kurf.initialise();
    kurf.setInitialConditionsFromPriors();
    
    // Run PSO
    KurPSO kurpso( kurf, paramcount, prior );
    kurpso.run();
    datatype* cost0 = new datatype;
    (*cost) = kurpso.getBestCost();
    (*map) = kurpso.getBestPos();
    
    cout << "MAP estimate, neg-log-likeli = " << *cost0 << endl;
    
    // Save MAP results
    //kurf->saveStates(  savedir + "/out_x_"     + std::to_string(threadno) + ".txt" );
    //kurf->saveObs(     savedir + "/out_y_"     + std::to_string(threadno) + ".txt" );
    //kurf->saveObsPred( savedir + "/out_ypred_" + std::to_string(threadno) + ".txt" );
}

M1 KurPSOchains::run(
    KuramotoUKF::ModelParamsSimple modelparams,
    int paramcount,
    std::vector<KurPSO::Prior> prior,
    std::vector<int> paramPriorList,
    int n_priors,
    std::string loadfile,
    std::string savedir,
    int threadcount )
{
    // Form pointer array for return values from threads
    std::vector<std::shared_ptr<M1>> maps(threadcount,make_shared<M1>());
    std::vector<std::shared_ptr<datatype>> costs(threadcount);

    // Initialise threads list
    std::vector<std::thread> psoThreads(threadcount);
    
    // Call threads
    for ( int threadno = 0; threadno < threadcount; threadno++ )
        psoThreads[threadno] = std::thread( &KurPSOchains::runThread, this, modelparams, paramcount, prior, paramPriorList, n_priors, loadfile, savedir, maps[threadno], costs[threadno], threadno );
    
    // Join threads
    for ( int threadno = 0; threadno < threadcount; threadno++ )
        psoThreads[threadno].join();
    
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
    M1 bestmap = MatrixManip::allocMatrix(n_priors);
    for ( int k = 0; k < n_priors; k++ )
        bestmap[k] = (*maps[bestIndex])[k];
    
    // Return MAP estimate of parameters
    return bestmap;
}
