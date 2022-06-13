# -*- coding: utf-8 -*-
"""
Adaptive Sequential Monte Carlo inversion combined with Sequential Geostatistical Resampling for posterior inference and evidence computation.  This code controls the user-defined parameters,
 calls the main program asmc_sgr.py to performe the inversion and saves the results. 

This codes correspond to the article by Amaya et al. (2022). It is a Python 3.7 implementation of the Adaptive Sequential Monte Carlo (ASMC) method (Zhou et al., 2016; algorithm 4) to estimate 
the posterior probability density function (PDF) and the evidence (marginal likelihood) trough a Bayesian inversion. ASMC is a particle approach that relies on importance sampling over a sequence 
of intermediate distributions (power posteriors) that link the prior and the posterior PDF. Each power posterior is approximated by updating the particle importance weights and states using a small 
pre-defined number MCMC proposal steps. ASMC method adaptively tunes the sequence of power posteriors and performs resampling of particles when the variance of their importance weights 
becomes too large.

The test case is a synthetic groundwater transport problem from Laloy et al. (2016) with channelized categorical models representing the hydraulic conductivity spatial distribution. 

This ASMC implementation (referred to as ASMC-SGR, algorithm 1 in Amaya et al. (2022)) includes:

-an adaptation of the Sequential Geostatistical Resampling algorithm from Laloy et al. (2016) to generate new model proposals troughout the inversion using the direct sampling
approach (DS, Mariethoz et al., 2010). DS perfomrs multiple-point statistics conditioned siumlations based on a training image, and it is implemented using the
DeeSse algorithm (patented by the University of Neuchâtel, http://www.randlab.org/research/deesse/),

-an adaptation of the DREAMzs algorithm to perform MCMC steps approximating each power posterior (ter Braak and Vrugt, 2008; Vrugt, 2009; Laloy and Vrugt, 2012),

-the finite-volume open-source code MaFloT for transport simulations in porous media (Kunze & Lunati, 2012).

In case you have a question or if you find a bug, please write me an email to macarena.amaya@unil.ch. 


References:

Amaya, M., Linde, N., Laloy, E. (2022). Advances in Water Resources (under review). 

Amaya, M., Linde, N., & Laloy, E. (2021). Adaptive sequential Monte Carlo for posterior inference and model selection among complex geological priors. Geophysical Journal International, 226(2), 1220-1238.

Künze, R., & Lunati, I. (2012). An adaptive multiscale method for density-driven instabilities. Journal of Computational Physics, 231(17), 5557-5570.

Laloy, E., & Vrugt, J. A. (2012). High‐dimensional posterior exploration of hydrologic models using multiple‐try DREAM (ZS) and high‐performance computing. Water Resources Research, 48(1).

Laloy, E., Linde, N., Jacques, D., & Mariethoz, G. (2016). Merging parallel tempering with sequential geostatistical resampling for improved posterior exploration of high-dimensional subsurface categorical fields. Advances in water resources, 90, 57-69.

Mariethoz, G., Renard, P., & Straubhaar, J. (2010). The direct sampling method to perform multiple‐point geostatistical simulations. Water Resources Research, 46(11).

Ter Braak, C. J., & Vrugt, J. A. (2008). Differential evolution Markov chain with snooker updater and fewer chains. Statistics and Computing, 18(4), 435-446.

Vrugt, J. A., ter Braak, C., Diks, C., Robinson, B. A., Hyman, J. M., & Higdon, D. (2009). Accelerating Markov chain Monte Carlo simulation by differential evolution with self-adaptive randomized subspace sampling. International Journal of Nonlin ear Sciences and Numerical Simu- lation, 10(3), 273–290.

Zhou, Y., Johansen, A. M., & Aston, J. A., (2016). Toward automatic model comparison: an adaptive sequential Monte Carlo approach, Journal of Computational and Graphical Statistics,69925(3), 701–726.

                                                                                                                                                                                                               
"""


import os
import time
import numpy as np
import shutil

work_dir=os.getcwd()

import asmc_sgr

CaseStudy=1  
Restart=False
    

if  CaseStudy==1: 

    # User defined parameters:
    seq=72 # Number of particles (N)
    thin=36# Thinning parameter for saving the sampled likelihoods
    steps=36# Iterations per intermediate distribution (K)
    CESSf_div=0.9997 # targeted CESS (CESS_op) 
    ESSf_div=0.3 # ESS treshold (ESS*) 
    AR_max=35.0# Max acceptance rate before decreasing proposal scale
    AR_min=15.0 # Min acceptance rate before decreasing proposal scale
    
#    Init='True'
    
    ndraw=seq*500000# Set a high number of iterations to stop in case the beta sequence becomes too long due to a bad choice of CESS  
    tune_phi=seq*500000 # Iterations in which the proposal scale is tune it (for our case we tune it all along the run, so we define a high number of iterations that is not reached). 

    
    #Decide if to run forward solver in parallel
    
    DoParallel=True 
    parallel_jobs=seq
    MakeNewDir=False
    
            
     
#% Run the ASMC-SGR algorithm
            
if __name__ == '__main__':
    
    #start_time = time.time()

    q=asmc.Sampler(main_dir=work_dir,CaseStudy=CaseStudy,seq=seq,ndraw=ndraw,parallel_jobs=seq,steps=steps,
                   parallelUpdate = 1,pCR=False,thin=thin,nCR=3,DEpairs=1,pJumpRate_one=0.2,BoundHandling='Fold',
                   lik_sigma_est=False,DoParallel=DoParallel,CESSf_div=CESSf_div,ESSf_div=ESSf_div,AR_min=AR_min,AR_max=AR_max,tune_phi=tune_phi)
    
    print("Iterating")
    
    if Restart:
        tmpFilePath=work_dir+'/out_tmp.pkl'
    else:
        tmpFilePath=None 
    
#    Sequences, Z, OutDiag, fx, MCMCPar, MCMCVar, beta_run, jr_seq, weig_seq, CESS_ev, increment, ESS_ev, evid_cont, evid_ev, weights_unnorm, new_weight_ev, weig_cont, eve_seq = q.sample(RestartFilePath=tmpFilePath)
    Sequences, OutDiag, MCMCPar, MCMCVar, beta_run, phi_used_ev, phi_ev, weig_seq, CESS_ev, increment, ESS_ev, evid_cont, evid_ev, weights_unnorm, new_weight_ev, weig_cont, eve_seq, last_models,last_likelihoods,intermediate_models = q.sample(RestartFilePath=tmpFilePath)
    
    np.save('Sequences_states',Sequences) # Evolution of the states for every particle (latent parameters) and its likelihood
    np.save('AR',OutDiag.AR) # Acceptance Rate
    np.save('beta_seq',beta_run) # Sequence that defines the intermediate distributions (resulting from the adaptive procedure) 
    np.save('beta_increment',increment)
    np.save('CESS',CESS_ev)
    np.save('ESS',ESS_ev)    
    np.save('phi_used_ev',phi_used_ev)
    np.save('phi_ev',phi_ev)# Proposal scale evolution
    np.save('weig_ev',weig_seq) # Weights evolution	
    np.save('weig_cont',weig_cont) # Weights evolution	 	
    np.save('evid_ev',evid_ev) # Evidence evolution
    np.save('eve_ev',eve_seq) # Eve index evolution
    np.save('last_models.npy',last_models)
    np.save('last_likelihoods',last_likelihoods)
    np.save('intermediate_models ',intermediate_models )    