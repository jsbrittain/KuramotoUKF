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
void KurPSOchains::runThread( KuramotoUKF::ModelParamsSimple modelparams, int paramcount, KurPSO::Prior* prior, int* paramPriorList, int n_priors, std::string loadfile, std::string savedir, datatype** map, datatype** cost, int threadno ) {
    
    // Construct Kuromoto UKF
    KuramotoUKF* kurf = new KuramotoUKF( modelparams.kurfparams );
    kurf->readMeasurementFile(loadfile, 1);
    //kurf->setPriors( paramPriorList, n_priors, paramcount );
    kurf->formPriors( modelparams );
    kurf->initialise();
    kurf->setInitialConditionsFromPriors();
    
    // Run PSO
    KurPSO* kurpso = new KurPSO( kurf, paramcount, prior );
    kurpso->run();
    datatype* cost0 = new datatype;
    (*cost0) = kurpso->getBestCost();
    (*cost) = cost0;
    M1 Mmap0 = kurpso->getBestPos();
    (*map) = &Mmap0[0];
    
    cout << "MAP estimate, neg-log-likeli = " << *cost0 << endl;
    //MatrixManip::printVector( map0, paramcount );
    
    // Save MAP results
    //kurf->saveStates(  savedir + "/out_x_"     + std::to_string(threadno) + ".txt" );
    //kurf->saveObs(     savedir + "/out_y_"     + std::to_string(threadno) + ".txt" );
    //kurf->saveObsPred( savedir + "/out_ypred_" + std::to_string(threadno) + ".txt" );
    
    // Clean-up
    delete kurpso;
    delete kurf;
}

M1 KurPSOchains::run( KuramotoUKF::ModelParamsSimple modelparams, int paramcount, KurPSO::Prior* prior, int* paramPriorList, int n_priors, std::string loadfile, std::string savedir, int threadcount ) {
    
    // Form pointer array for return values from threads
    datatype** maps = new datatype*[threadcount];
    for ( int k = 0; k < threadcount; k++ )
        maps[k] = nullptr;
    datatype** costs = new datatype*[threadcount];
    for ( int k = 0; k < threadcount; k++ )
        costs[k] = nullptr;
    
    // Initialise threads list
    std::thread* psoThreads = new std::thread[threadcount];
    
    // Call threads
    for ( int threadno = 0; threadno < threadcount; threadno++ )
        psoThreads[threadno] = std::thread( &KurPSOchains::runThread, this, modelparams, paramcount, prior, paramPriorList, n_priors, loadfile, savedir, &maps[threadno], &costs[threadno], threadno );
    
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
        bestmap[k] = maps[bestIndex][k];
    //MatrixManip::printVector( bestmap, n_priors );
    
    // Tidy up
    delete[] psoThreads;
    delete[] maps;
    delete[] costs;
    
    // Return MAP estimate of parameters
    return bestmap;
}
