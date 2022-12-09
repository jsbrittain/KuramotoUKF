//
//  KuramotoMC.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 15/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#include "KuramotoMC.hpp"

KuramotoMCChain::KuramotoMCChain(KuramotoUKF& kurf, MetropolisChainParams mcmcparams )
  : kuramotoukf(kurf), MetropolisChain( mcmcparams ) {
    //
};
datatype KuramotoMCChain::negLogLikeliFcn( M1 state, std::vector<Prior> prior ) {
    return kuramotoukf.negLogLikeliFcn( state, prior );
}
