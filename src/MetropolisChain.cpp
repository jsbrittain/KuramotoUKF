//
//  MetropolisChain
//  mcmc
//
//  Created by John-Stuart Brittain on 14/12/2016.
//  Copyright © 2016 John-Stuart Brittain. All rights reserved.
//

#include "MetropolisChain.hpp"

MetropolisChain::MetropolisChain( MetropolisChainParams mcmcparams) : state0(mcmcparams.state0), statedim(mcmcparams.statedim), covarProposal(mcmcparams.covarProposal), tuneProposal(mcmcparams.tuneProposal), tuningiters(mcmcparams.tuningiters), prior(mcmcparams.prior), randseed(mcmcparams.randseed), burnin(mcmcparams.burnin), chainlength(mcmcparams.chainlength), verbose(mcmcparams.verbose) {
    initialise();
}
MetropolisChain::MetropolisChain( M1 state0, M2 covarProposal, bool tuneProposal, int tuningiters, std::vector<Prior>& prior, int statedim, int randseed, int burnin, int chainlength, int verbose ) : state0(state0), statedim(statedim), covarProposal(covarProposal), tuneProposal(tuneProposal), tuningiters(tuningiters), prior(prior), randseed(randseed), burnin(burnin), chainlength(chainlength), verbose(verbose) {
    initialise();
}
MetropolisChain::~MetropolisChain()
{
}
void MetropolisChain::initialise( ) {
    // Initialise at start of chain
    chaint = 0;
    // Set random seed
    srand( randseed );
    // Allocate memory for chain (current) and proposal state values
    assert( chainlength >= burnin );
    stateCurrent = allocMatrix(chainlength,statedim);
    stateNegLogLikeli = allocMatrix(chainlength);
    stateProposal = allocMatrix(statedim);
    // Update initial state vector
    for ( int k = 0; k < statedim; k++ )
        stateCurrent[chaint][k] = state0[k];
}
void MetropolisChain::run() {
    
    // Initial state
    chaint = 0;
    datatype acceptance_rate;
    clock_t timer = clock();
    logPcurrent = negLogLikeliFcn( stateCurrent[chaint], prior, statedim );
    stateNegLogLikeli[chaint] = logPcurrent;
    if ( verbose ) {
        timer = clock() - timer;
        float tsecs = ((float)timer)/CLOCKS_PER_SEC;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "First run diagnostics..." << std::endl;
        std::cout << "  Log-likelihood of model: " << logPcurrent << std::endl;
        std::cout << "  Singe sweep execution time: " << tsecs << " secs." << std::endl;
        std::cout << "  Estimated time for " << burnin << " + " << chainlength << " iterations: " << tsecs*(burnin+chainlength)/60 << " mins, " << tsecs*(burnin+chainlength) << " secs." << std::endl;
    }
    
    // Approximate covarProposal based on univariate searches
    if ( tuneProposal ) {
        // Keep best solution after burnin
        M1 beststate = getMAPstate();
        
        // Determine univariate acceptance rates
        datatype desired_acceptance_rate = 0.234;
        datatype univariate_acceptance_rate = 0.50;//pow( desired_acceptance_rate, 1.0/(datatype)statedim );
        datatype tolerance = 0.2;     // 20% tolerance (in atanh space)
        datatype lower_bound = tanh(atanh(univariate_acceptance_rate)*(1-tolerance));
        datatype upper_bound = tanh(atanh(univariate_acceptance_rate)*(1+tolerance));
        std::cout << "Desired acceptance rate: " << desired_acceptance_rate << std::endl;
        std::cout << "Target univariate acceptance rate: " << univariate_acceptance_rate << std::endl;
        
        // Backup full covarProposal
        M2 covarProposal0 = MatrixManip::allocMatrix(statedim,statedim);
        for ( int i = 0; i < statedim; i++ )
            covarProposal0[i][i] = covarProposal[i][i];
        
        // Tune variance based on acceptance rate per parameter
        datatype scale_adjust, prevadjust = 0, logPcurrent0 = logPcurrent;
        for ( int k = 0; k < statedim; k++ ) {
            // Reset some parameters
            scale_adjust = 1e1;
            prevadjust = 0;
            // Zero covar matrix except kth diagonal entry
            for ( int i = 0; i < statedim; i++ )
                covarProposal[i][i] = 0;
            covarProposal[k][k] = covarProposal0[k][k];
            while ( true ) {
                logPcurrent = logPcurrent0;
                acceptance_rate = steps( tuningiters );
                if ( verbose ) std::cout << "Acceptance rate for parameter [ " << k+1 << " / " << statedim << " ]: " << acceptance_rate << std::endl;
                // Adjust single parameter variance based on burn-in
                if ( acceptance_rate < lower_bound ) {
                    if ( prevadjust > 0 )
                        scale_adjust /= 2;
                    prevadjust = -1;
                    covarProposal[k][k] /= scale_adjust;
                    if ( verbose ) std::cout << "Decreasing proposal (log10) variance to " << log10(covarProposal[k][k]) << " (scale adjust " << scale_adjust << ")" << std::endl;
                    continue;
                }
                if ( acceptance_rate > upper_bound ) {
                    if ( prevadjust < 0 )
                        scale_adjust /= 2;
                    prevadjust = +1;
                    covarProposal[k][k] *= scale_adjust;
                    if ( log10(covarProposal[k][k]) > 20 ) {
                        if ( verbose )
                            std::cout << "Proposal distribution exceeds any reasonable range..." << std::endl << " This parameter appears to play no part in the log-likelihood evaluation." << std::endl;
                        break;
                    }
                    if ( verbose ) std::cout << "Increasing proposal (log10) variance to " << log10(covarProposal[k][k]) << " (scale adjust " << scale_adjust << ")" << std::endl;
                    continue;
                }
                if ( verbose ) std::cout << " parameter within tolerance." << std::endl;
                break;
            }
            // Save adjusted variance to proposal
            covarProposal0[k][k] = covarProposal[k][k];
            // Save univariate chain
            //saveChain("/Users/jsb/proposal_" + std::to_string(k) + ".txt");
        }
        // Restore full (adjusted) covariance matrix
        for ( int i = 0; i < statedim; i++ )
            covarProposal[i][i] = covarProposal0[i][i];
        // Report output
        if ( verbose) std::cout << "Proposal distribution tuned." << std::endl;
        
        // Make a chain of two elements; original and beststate after burn in
        for ( int k = 0; k < statedim; k++ )
            stateCurrent[1][k] = beststate[k];
        chaint = 1; chainlength = 2;
        
        //delete[] beststate;
        return;
    }
    
    // Burn-in
    if ( burnin > 0 ) {
        if ( verbose ) std::cout << "Burning-in..." << std::endl;
        // Run burn-in
        acceptance_rate = steps( burnin );
        // Report acceptance rate
        if ( verbose ) std::cout << "Burn-in acceptance rate " << acceptance_rate << std::endl;
    }
    // Set last state to initial condition and reset time
    if ( chaint > 0 ) {
        for ( int k = 0; k < statedim; k++ )
            stateCurrent[0][k] = stateCurrent[chaint-1][k];
        chaint = 0;
    }
    
    // Iterate through chain using last (t-1) sample as current
    if ( verbose ) std::cout << "Running chain..." << std::endl;
    acceptance_rate = steps( chainlength );
    
    // Report final acceptance rate
    if ( verbose ) std::cout << "Final acceptance rate " << acceptance_rate << std::endl;
}
datatype MetropolisChain::steps( int maxcount ) {
    // Iterate through Metropolis steps
    //
    // This extended function is mainly to provide feedback when in verbose mode
    
    int progress = 0, newprogress = 0, accepted = 0;
    if ( verbose ) std::cout << "<" << std::flush;
    for ( chaint = 1; chaint < maxcount; chaint++ ) {
        // Run Metropolis step
        if (step()) accepted++;
        
        // Report progress
        if ( verbose ) {
            newprogress = floor( 100*chaint/maxcount );
            if ( newprogress > progress ) {
                progress = newprogress;
                if ( ((progress % 10) == 0) && (progress!=100) )
                    std::cout << "|" << std::flush;
                else {
                    if ( progress < 100 )
                        std::cout << "." << std::flush;
                }
            }
        }
    }
    if ( verbose ) std::cout << ">" << std::endl;       // endl will flush the stream
    return ((datatype)accepted/((datatype)maxcount-1));
}
bool MetropolisChain::step() {
    // Metropolis step //
    
    // Generate proposal -- sample from proposal distribution
    if (!tuneProposal)
        mvnrand(stateCurrent[chaint-1], statedim, covarProposal, stateProposal);
    else {
        // Only one variables changes at a time when tuning
        M1 mu = allocMatrix(1);
        M2 prop = allocMatrix(1, 1);
        for ( int k = 0; k < statedim; k++ ) {
            if (covarProposal[k][k] != 0 ) {
                mu[0] = stateCurrent[chaint-1][k];
                prop[0][0] = covarProposal[k][k];
                mvnrand(mu, 1, prop, stateProposal[k]);
            } else
                stateProposal[k] = stateCurrent[chaint-1][k];
        }
    }
    
    // Determine proposal likelihood
    logPproposal = negLogLikeliFcn( stateProposal, prior, statedim );
    logPcrit = logPcurrent - logPproposal;
    
    // Acceptance criteria
    bool accepted = ( logRandUniform() < logPcrit );
    if ( accepted ) {
        // Accept proposal; update chain sample to proposal
        for ( int k = 0; k < statedim; k++ )
            stateCurrent[chaint][k] = stateProposal[k];
        logPcurrent = logPproposal;
    } else {
        // Reject proposal
        for ( int k = 0; k < statedim; k++ )
            stateCurrent[chaint][k] = stateCurrent[chaint-1][k];
    }
    stateNegLogLikeli[chaint] = logPcurrent;
    return accepted;
}
datatype MetropolisChain::negLogLikeliTestFcn( M1 state ) {
    M1 mu = allocMatrix(statedim);
    M2 Sigma = mateye( statedim );
    for ( int i = 0; i < statedim; i++ ) {
        mu[i] = i - static_cast<datatype>(statedim)/2.0;
        Sigma[i][i] = (i+1);
    }
    datatype logLikeli = logLikeliMVN( state, statedim, mu, Sigma );
    return -logLikeli;
}
void MetropolisChain::saveChain( const std::string filename ) {
    // First column is neg-log-likelihood
    std::ofstream myfile(filename,std::ofstream::binary);
    myfile << std::scientific << std::setprecision(std::numeric_limits<datatype>::max_digits10);
    for ( int t = 0; t < chainlength; t++ ) {
        myfile << stateNegLogLikeli[t] << '\t';
        for ( int k = 0; k < statedim; k++ )
            myfile << stateCurrent[t][k] << '\t';
        myfile << std::endl;
    }
}
datatype MetropolisChain::getMAPcost( ) {
    datatype bestcost = stateNegLogLikeli[0];
    for ( int k = 1; k < chainlength; k++ ) {
        if ( stateNegLogLikeli[k] < bestcost ) {
            bestcost = stateNegLogLikeli[k];
        }
    }
    return bestcost;
}
M1 MetropolisChain::getMAPstate( ) {
    int bestpos = 0;
    datatype bestcost = stateNegLogLikeli[0];
    for ( int k = 1; k < chainlength; k++ ) {
        if ( stateNegLogLikeli[k] < bestcost ) {
            bestpos = k;
            bestcost = stateNegLogLikeli[k];
        }
    }
    // Allocate new memory for return
    M1 MAPstate;
    for ( int k = 0; k < statedim; k++ )
        MAPstate[k] = stateCurrent[bestpos][k];
    return MAPstate;
}
