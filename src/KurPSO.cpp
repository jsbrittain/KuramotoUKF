//
//  KurPSO.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 17/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#include "KurPSO.hpp"

KurPSO::KurPSO( KuramotoUKF kurf, int paramcount, std::vector<Prior> prior )
  : kurf(kurf), ParticleSwarmOptimiser( paramcount, prior ) {
    //
}
datatype KurPSO::costFunction( M1 state, std::vector<Prior> prior, int n ) {
    return kurf.negLogLikeliFcn( state, prior );
}
