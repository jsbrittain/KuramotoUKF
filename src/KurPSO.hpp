//
//  KurPSO.hpp
//  mcmc
//
//  Created by John-Stuart Brittain on 17/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#ifndef KurPSO_hpp
#define KurPSO_hpp

#include "ParticleSwarmOptimiser.hpp"
#include "KuramotoUKF.hpp"

class KurPSO : public ParticleSwarmOptimiser {
public:
    KuramotoUKF kurf;
    
    KurPSO( KuramotoUKF kurf, int paramcount, std::vector<Prior> prior );
    datatype costFunction( M1 state, std::vector<Prior> prior, int n ) override;
};

#endif /* KurPSO_hpp */
