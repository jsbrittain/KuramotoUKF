//
//  ParticleSwarmOptimiser.cpp
//  mcmc
//
//  Created by John-Stuart Brittain on 17/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#include "ParticleSwarmOptimiser.hpp"

ParticleSwarmOptimiser::ParticleSwarmOptimiser( int paramcount, std::vector<Prior> prior )
  : paramcount(paramcount), prior(prior) {
    //
}
ParticleSwarmOptimiser::~ParticleSwarmOptimiser() {
}
void ParticleSwarmOptimiser::initialise() {
    
    particlecount = 10*paramcount;
    
    // PSO parameters
    psoparams.phi1=2.05;
    psoparams.phi2=2.05;
    psoparams.phi=psoparams.phi1+psoparams.phi2;
    psoparams.chi=2/(psoparams.phi-2+std::sqrt(psoparams.phi*psoparams.phi-4*psoparams.phi));
    psoparams.w=psoparams.chi;                      // Inertia Weight
    psoparams.wdamp=1;                              // Inertia Weight Damping Ratio
    psoparams.c1=psoparams.chi*psoparams.phi1;      // Personal Learning Coefficient
    psoparams.c2=psoparams.chi*psoparams.phi2;
    psoparams.velocityRangeScale = 1.0;
    psoparams.velocityMin = allocMatrix(paramcount);
    psoparams.velocityMax = allocMatrix(paramcount);
    for ( int i = 0; i < paramcount; i++ ) {
        psoparams.velocityMin[i] = -psoparams.velocityRangeScale*prior[i].sd;
        psoparams.velocityMax[i] =  psoparams.velocityRangeScale*prior[i].sd;
    }
    
    // Convert priors to vectorised form
    M1 zeros = allocMatrix(paramcount);
    M1 priorMu = allocMatrix(paramcount);
    M2 priorVariance = allocMatrix(paramcount,paramcount);
    for ( int k = 0; k < paramcount; k++ ) {
        priorMu[k] = prior[k].mu;
        priorVariance[k][k] = prior[k].sd*prior[k].sd;
    }
    
    // Update display
    if ( psoparams.verbose )
        std::cout << "Initialising " << particlecount << " particles..." << std::endl;
    
    // Construct particles
    particle = std::vector<particleStruct>(particlecount);
    for ( int k = 0; k < particlecount; k++ ) {
        // Location
        particle[k].position = allocMatrix(paramcount);
        particle[k].velocity = allocMatrix(paramcount);
        particle[k].best.position = allocMatrix(paramcount);
        particle[k].best.velocity = allocMatrix(paramcount);
        particle[k].neighbourhoodbest.position = allocMatrix(paramcount);
        particle[k].neighbourhoodbest.velocity = allocMatrix(paramcount);
        
        // Repeat until particle has a valid cost
        particle[k].cost = NAN;
        int attempt_count = 0;
        while ( std::isnan(particle[k].cost) ) {
            // Try initialising each particle at most 15 times
            if ( (attempt_count++) > 15 )
                break;
            // Randomly draw starting location
            mvnrand( priorMu, paramcount, priorVariance, particle[k].position );
            mvnrand( zeros,   paramcount, priorVariance, particle[k].velocity );
            // Evaluate cost and assign
            particle[k].cost = costFunction( particle[k].position, prior, paramcount );
        }
        
        // Define neighbours
        particle[k].neighbour = std::vector<int>(neighbours);
        for( int i = 0; i < neighbours; i++ ) {
            int x = k+i-1-2*floor(neighbours/2)+neighbours;
            // Neighbours include a self-reference
            particle[k].neighbour[i] = ( x - particlecount*floor(x/particlecount) );
        }
        
        // Log particle best
        for ( int i = 0; i < paramcount; i++ ) {
            particle[k].best.position[i] = particle[k].position[i];
            particle[k].best.velocity[i] = particle[k].velocity[i];
        }
        particle[k].best.cost = particle[k].cost;
    }
    
    // Find global best
    int bestindex = 0;
    bestcost = particle[0].cost;
    for ( int k = 0; k < particlecount; k++ ) {
        if ( particle[k].cost < bestcost ) {
            bestindex = k;
            bestcost = particle[k].cost;
        }
    }
    globalbest.position = allocMatrix(paramcount);
    globalbest.velocity = allocMatrix(paramcount);
    for ( int i = 0; i < paramcount; i++ ) {
        globalbest.position[i] = particle[bestindex].position[i];
        globalbest.velocity[i] = particle[bestindex].velocity[i];
    }
    globalbest.cost = particle[bestindex].cost;
    
    // Update display
    if ( psoparams.verbose ) {
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Best starting cost = " << bestcost << std::endl;
        std::cout << "Starting optimisation..." << std::endl;
    }
}
void ParticleSwarmOptimiser::run() {
    // Initialise particles
    initialise();
    
    M2 unirndbest = allocMatrix(1,paramcount);
    M2 unirndneighbest = allocMatrix(1,paramcount);
    
    // Iterate until convergence
    bestcount = 0;
    for ( int it = 0; it < maxiters; it++ ) {
        
        // Update neighbourhood best
        for ( int k = 0; k < particlecount; k++ ) {
            // Overwrite neighbourhood best with particle best
            particle[k].neighbourhoodbest.cost = particle[k].cost;
            for ( int i = 0; i < paramcount; i++ ) {
                particle[k].neighbourhoodbest.position[i] = particle[k].position[i];
                particle[k].neighbourhoodbest.velocity[i] = particle[k].velocity[i];
            }
            // Traverse neighbours for best
            for ( int j = 0; j < neighbours; j++ ) {
                if ( particle[particle[k].neighbour[j]].best.cost < particle[k].neighbourhoodbest.cost ) {
                    for ( int i = 0; i < paramcount; i++ ) {
                        particle[k].neighbourhoodbest.position[i] = particle[particle[k].neighbour[j]].position[i];
                        particle[k].neighbourhoodbest.velocity[i] = particle[particle[k].neighbour[j]].velocity[i];
                    }
                    particle[k].neighbourhoodbest.cost = particle[particle[k].neighbour[j]].cost;
                }
            }
        }
        
        // Update particle trajectories
        for ( int k = 0; k < particlecount; k++ ) {
            
            // Update particle parameters
            uniformrand(1, paramcount, unirndbest);
            uniformrand(1, paramcount, unirndneighbest);
            for ( int i = 0; i < paramcount; i++ ) {
                
                // Update Velocity
                particle[k].velocity[i] = psoparams.w*particle[k].velocity[i] + psoparams.c1*unirndbest[0][i]*(particle[k].best.position[i]-particle[k].position[i]) + psoparams.c2*unirndneighbest[0][i]*(particle[k].neighbourhoodbest.position[i]-particle[k].position[i]);
                
                // Apply Velocity Limits
                particle[k].velocity[i] = getmax(particle[k].velocity[i],psoparams.velocityMin[i]);
                particle[k].velocity[i] = getmin(particle[k].velocity[i],psoparams.velocityMax[i]);
                
                // Update Position
                particle[k].position[i] = particle[k].position[i] + particle[k].velocity[i];
            }
            // Evaluation
            particle[k].cost = costFunction( particle[k].position, prior, paramcount );
            
            // If cost function is NaN, revert particle to previous best
            if ( std::isnan( particle[k].cost ) ) {
                for ( int i = 0; i < paramcount; i++ ) {
                    particle[k].position[i] = particle[k].best.position[i];
                    particle[k].velocity[i] = particle[k].best.velocity[i];
                }
                particle[k].cost = particle[k].best.cost;
            }
            
            // Update Personal Best
            if ( particle[k].cost < particle[k].best.cost ) {
                for ( int i = 0; i < paramcount; i++ ) {
                    particle[k].best.position[i] = particle[k].position[i];
                    particle[k].best.velocity[i] = particle[k].velocity[i];
                }
                particle[k].best.cost = particle[k].cost;
            }
        }
        
        // Identify global best
        for ( int k = 0; k < particlecount; k++ ) {
            if ( particle[k].best.cost < globalbest.cost ) {
                for ( int i = 0; i < paramcount; i++ ) {
                    globalbest.position[i] = particle[k].position[i];
                    globalbest.velocity[i] = particle[k].velocity[i];
                }
                globalbest.cost = particle[k].cost;
            }
        }
        lastbestcost = bestcost;
        bestcost = globalbest.cost;
        
        // Update dampening
        psoparams.w = psoparams.w*psoparams.wdamp;
        
        // Update display
        if ( psoparams.verbose )
            std::cout << "  Iteration " << it << ": Best Cost = " << bestcost << ", Change = " << std::abs( (bestcost-lastbestcost)/lastbestcost ) << " [Thr=" << psoparams.convthreshold << "]" << std::endl;
        
        // Check for convergence
        if ( std::abs( (bestcost-lastbestcost)/lastbestcost ) < psoparams.convthreshold ) {
            bestcount = bestcount + 1;
            if ( bestcount >= psoparams.convreps ) {
                // Cost has not changed for n steps --- end optimisation
                if ( psoparams.verbose )
                    std::cout << "   convergence condition met." << std::endl;
                break;
            }
        } else
            bestcount = 0;
    }
}
datatype ParticleSwarmOptimiser::getBestCost() {
    return globalbest.cost;
}
M1 ParticleSwarmOptimiser::getBestPos() {
    M1 map = allocMatrix(paramcount);
    for ( int i = 0; i < paramcount; i++ )
        map[i] = globalbest.position[i];
    return map;
}
