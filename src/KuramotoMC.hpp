//
//  KuramotoMC.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 15/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#ifndef KuramotoMC_hpp
#define KuramotoMC_hpp

#include "KuramotoUKF.hpp"
#include "MetropolisChain.hpp"

class KuramotoMCChain : public MetropolisChain {
public:
    KuramotoUKF kuramotoukf;
    KuramotoMCChain( KuramotoUKF* kurf, MetropolisChainParams kurmcparams );
    datatype negLogLikeliFcn( M1 state, Prior* prior, int n ) override;
};

#endif /* KuramotoMC_hpp */
