"""
Adaptive Sequential Monte Carlo inversion combined with Sequential Geostatistical Resampling for posterior inference and evidence computation.  This is the main code for the inversion, 
to modify the user-defined parameters and run the inversion please see run_asmc_sgr.py

This code corresponds to the article by Amaya et al. (2022). It is a Python 3.7 implementation of the Adaptive Sequential Monte Carlo (ASMC) method (Zhou et al., 2016; algorithm 4) to estimate 
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
from __future__ import print_function

import numpy as np
import numpy.matlib as matlib
try:
    import cPickle as pickle
except:
    import pickle

import time

from scipy.stats import triang

from asmc_sgr_func import* # This imports all ASMC and inverse problem-related functions

import sys

from attrdict import AttrDict

import scipy.io

import matlab.engine

MCMCPar=AttrDict()

MCMCVar=AttrDict()

Measurement=AttrDict()

OutDiag=AttrDict()

Extra=AttrDict()


from run_DS import runDS
from InitializeInput import InitializeInput
from SetInput import SetInput
from run_forward import run_forward
import matplotlib.pyplot as plt
import multiprocessing
from joblib import Parallel, delayed



class Sampler:

    

    
    def __init__(self, main_dir=None,CaseStudy=0,seq = 3,ndraw=10000,thin = 1,  nCR = 3, 
                 DEpairs = 3, parallelUpdate = 1.0, pCR=True,k=10,pJumpRate_one=0.2,
                 steps=100,savemodout=False, saveout=True,save_tmp_out=True,Prior='LHS',
                 DoParallel=True,eps=5e-2,BoundHandling='Reflect',
                 lik_sigma_est=False,parallel_jobs=4,rng_seed=123,it_b=10,jr_factor=0.2,CESSf_div=0.999993,ESSf_div=0.5,AR_min=15.0,AR_max=30.0,tune_phi=1000):
        
        
        
        
        self.CaseStudy=CaseStudy
        MCMCPar.seq = seq
        MCMCPar.ndraw=ndraw
        MCMCPar.thin=thin
        MCMCPar.nCR=nCR
        MCMCPar.DEpairs=DEpairs
        MCMCPar.parallelUpdate=parallelUpdate
        MCMCPar.Do_pCR=pCR
        MCMCPar.k=k
        MCMCPar.pJumpRate_one=pJumpRate_one
        MCMCPar.steps=steps
        MCMCPar.savemodout=savemodout
        MCMCPar.saveout=saveout  
        MCMCPar.save_tmp_out=save_tmp_out  
        MCMCPar.Prior=Prior
        MCMCPar.DoParallel=DoParallel
        MCMCPar.eps = eps
        MCMCPar.BoundHandling = BoundHandling
        MCMCPar.lik_sigma_est=lik_sigma_est
        Extra.n_jobs=parallel_jobs
        Extra.main_dir=main_dir
        np.random.seed(seed=None) 
        MCMCPar.it_b=it_b
        MCMCPar.jr_factor=jr_factor
        MCMCPar.AR_min=AR_min
        MCMCPar.AR_max=AR_max      
        MCMCPar.ESSf_div=ESSf_div
        MCMCPar.CESSf_div=CESSf_div
        MCMCPar.tune_phi=tune_phi



        if self.CaseStudy==1:   
            
### Load data        
            mat = scipy.io.loadmat('data_03.mat')
            an_array = np.array(list(mat.items()))
            meas_cc=an_array[4,1]
            Measurement.MeasData=meas_cc
            sigma=an_array[5,1]
            Measurement.Sigma=sigma[0,0]
            Measurement.N=np.size(Measurement.MeasData) 
            

### Model size, likelihood type, proposal scale and pumping wells location
            
            ModelName = 'runTransport'
            MCMCPar.lik=2
            MCMCPar.nx=75
            MCMCPar.ny=101
            MCMCPar.nz=1            
            MCMCPar.cond_type=2
            
            if MCMCPar.cond_type==2:
                MCMCPar.use_phi_dist='Triangular'
                MCMCPar.phi_cf='None'
                phi_min=np.ones((1,MCMCPar.seq))*5
#                print(phi_min)
                phi_fixed=np.ones((1,MCMCPar.seq))*50
                     
                MCMCPar.phi=np.ones(MCMCPar.seq)*50
                MCMCPar.phi_max=np.ones(MCMCPar.seq)*50
                MCMCPar.phi_min=np.ones(MCMCPar.seq)*5        
             
            MCMCPar.dc_set=np.array([[3,51,1,0],[10,51,1,1],[17,51,1,1],[24,51,1,0],[31,51,1,0],[38,51,1,1],[45,51,1,0],[52,51,1,0],[59,51,1,0],[66,51,1,0],[73,51,1,0]])
            MCMCPar.sf=1

                       
#           Define the beta list to choose from in the binary search, to crate an adaptative beta sequence.
                       
            gamma=np.linspace(-1,1,400)

            MCMCPar.gamma=gamma

                    
        self.MCMCPar=MCMCPar
        self.Measurement=Measurement
        self.Extra=Extra
        self.ModelName=ModelName

### Initialization       
        
    def _init_sampling(self):
        
        start_time = time.time()
        
        Iter=self.MCMCPar.seq
        iteration=2
        iloc=0
        T=0
        
        Init=True
        InitializeInput(MCMCPar.seq, Init)
        out=np.zeros((MCMCPar.seq,MCMCPar.nx*MCMCPar.ny))
        
        str0='unc_'
        
        m_current=np.zeros((MCMCPar.seq,MCMCPar.nx,MCMCPar.ny))
        
        for i in range(MCMCPar.seq):    
            out[i,:]=runDS(i+1,str0,0)
            m_current[i,:,:]=np.reshape(out[i,:],(MCMCPar.nx,MCMCPar.ny),order='F')
        

       #Run the forward solver Matflow

        u=np.linspace(1,MCMCPar.seq,MCMCPar.seq).astype(int)  
        njobs=MCMCPar.seq
        start_engine_paralel='No'  
        
       
        if start_engine_paralel=='No':
    
            eng=[] #empy engine as it is initialized at every iteration within the function run forward
            
            time1=time.time()
            
            d = run_forward_par(njobs,m_current,u,eng,start_engine_paralel) 

            time2=time.time()

            print('time run_forw',time2-time1)

           
        self.MCMCPar.n=self.MCMCPar.nx*self.MCMCPar.ny
        self.MCMCPar.CR=np.cumsum((1.0/self.MCMCPar.nCR)*np.ones((1,self.MCMCPar.nCR)))
        Nelem=np.floor(self.MCMCPar.ndraw/self.MCMCPar.seq)++self.MCMCPar.seq*2
        OutDiag.CR=np.zeros((np.int(np.floor(Nelem/self.MCMCPar.steps))+2,self.MCMCPar.nCR+1))
        OutDiag.AR=np.zeros((np.int(np.floor(Nelem/self.MCMCPar.steps))+2,2))
        OutDiag.AR[0,:] = np.array([self.MCMCPar.seq,-1])
        pCR = (1.0/self.MCMCPar.nCR) * np.ones((1,self.MCMCPar.nCR))
        
        # Calculate the actual CR values based on pCR
        CR,lCR = GenCR(self.MCMCPar,pCR)  
#        
        if self.MCMCPar.savemodout:
            self.fx = np.zeros((self.Measurement.N,np.int(np.floor(self.MCMCPar.ndraw/self.MCMCPar.thin))))
            MCMCVar.m_func = self.MCMCPar.seq     

        self.Sequences = np.empty((np.int(np.floor(Nelem/self.MCMCPar.thin))+1,1,self.MCMCPar.seq))
           
        self.MCMCPar.Table_JumpRate=np.zeros((self.MCMCPar.n,self.MCMCPar.DEpairs))
        for zz in range(0,self.MCMCPar.DEpairs):
            self.MCMCPar.Table_JumpRate[:,zz] = 2.38/np.sqrt(2 * (zz+1) * np.linspace(1,self.MCMCPar.n,self.MCMCPar.n).T)

        
### Compute likelihood from simulated data    

        of,log_p = CompLikelihood(m_current,d,self.MCMCPar,self.Measurement,self.Extra)
        
        X = np.concatenate((np.reshape(m_current,(self.MCMCPar.seq,MCMCPar.n)),of,log_p),axis=1)
        
        Xfx = d
    
        
        if self.MCMCPar.savemodout==True:
            self.fx=fx0
        else:
            self.fx=None


        aux_in=np.reshape(X.T,(1,self.MCMCPar.n+2,self.MCMCPar.seq))
        self.Sequences[0,0,:self.MCMCPar.seq] = aux_in[0,self.MCMCPar.n+1,:self.MCMCPar.seq]
        # Store N_CR
        OutDiag.CR[0,:MCMCPar.nCR+1] = np.concatenate((np.array([Iter]).reshape((1,1)),pCR),axis=1)
        delta_tot = np.zeros((1,self.MCMCPar.nCR))

     
        self.OutDiag=OutDiag    
        MCMCVar.Iter=Iter
        MCMCVar.iteration=iteration
        MCMCVar.iloc=iloc; MCMCVar.T=T; MCMCVar.X=X
        MCMCVar.Xfx=Xfx; MCMCVar.CR=CR; MCMCVar.pCR=pCR
        MCMCVar.lCR=lCR; MCMCVar.delta_tot=delta_tot
        self.MCMCVar=MCMCVar
        MCMCVar.m_current=m_current
        
        if self.MCMCPar.save_tmp_out==True:
            with open('out_tmp'+'.pkl','wb') as f:
                 pickle.dump({'Sequences':self.Sequences,
                 'OutDiag':self.OutDiag,'fx':self.fx,'MCMCPar':self.MCMCPar,
                 'MCMCVar':self.MCMCVar,'Measurement':self.Measurement,
                 'ModelName':self.ModelName,'Extra':self.Extra},f, protocol=pickle.HIGHEST_PROTOCOL)
    
        end_time = time.time()
        
    
        print("init_sampling took %5.4f seconds." % (end_time - start_time))
        

### Initialization for the main sampling loop      

    def sample(self,RestartFilePath=None):
        start_time1b = time.time()
        
        if not(RestartFilePath is None):
            print('This is a restart')
            with open(RestartFilePath, 'rb') as fin:
                tmp_obj = pickle.load(fin)
            self.Sequences=tmp_obj['Sequences']
            self.Z=tmp_obj['Z']
            self.OutDiag=tmp_obj['OutDiag']
            self.fx=tmp_obj['fx']
            self.MCMCPar=tmp_obj['MCMCPar']
            self.MCMCVar=tmp_obj['MCMCVar']
            self.Measurement=tmp_obj['Measurement']
            self.ModelName=tmp_obj['ModelName']
            self.Extra=tmp_obj['Extra']
            del tmp_obj
            
            self.ndim=self.MCMCPar.n
#                
            self.MCMCPar.ndraw = 2 * self.MCMCPar.ndraw
            
            # Reset rng
            np.random.seed(np.floor(time.time()).astype('int'))
            
### Extend Sequences, Z, OutDiag.AR,OutDiag.Rstat and OutDiag.CR
            self.Sequences=np.concatenate((self.Sequences,np.zeros((self.Sequences.shape))),axis=0)
            self.OutDiag.AR=np.concatenate((self.OutDiag.AR,np.zeros((self.OutDiag.AR.shape))),axis=0)
            self.OutDiag.CR=np.concatenate((self.OutDiag.CR,np.zeros((self.OutDiag.CR.shape))),axis=0)
      
            
        else:
            self._init_sampling()
        

### Initialize variables and arrays to store the results:

        prov_AR=30     # (no real meaning, just to not change the jr_scale on the first loop),   
        beta_run=[]
        increment=[]
        phi_used=np.zeros(self.MCMCPar.seq)
        phi_used_ev=np.zeros((1,self.MCMCPar.seq))
        phi_ev=np.zeros((1,self.MCMCPar.seq))
        
        likelihood_ev=[]
        weig_seq=[]
        weig_cont=[]
        weights_unnorm=[]
        weig_unn=np.ones(self.MCMCPar.seq)
        norm_weight=np.ones(self.MCMCPar.seq)/self.MCMCPar.seq
        new_weight_ev=[]
        weig_seq=np.append(weig_seq,norm_weight)
        
        CESSf=self.MCMCPar.CESSf_div*self.MCMCPar.seq
        CESS_ev=[]
        ESS_ev=[]
        beta=0.
        omega=-6
        beta_run=np.append(beta_run,beta)
        
        ind=0
        num_inc=self.MCMCPar.gamma.shape[0]
        
        evid_cont=[]
        evid_ev=[]
        evid_evolution=0
        
        eve_prev=np.arange(0,self.MCMCPar.seq,dtype=int)
        anc_prev=np.arange(0,self.MCMCPar.seq,dtype=int)
        
        position=0
        counter=0.
        intermediate_models=[]
        
        eve_seq=[]
        eve_seq=np.append(eve_seq,eve_prev)
             
        end_time1b = time.time()
        print("The initialization of SAMPLE (before the main loop) took %5.4f seconds." % (end_time1b - start_time1b))
        
        u=np.linspace(1,MCMCPar.seq,MCMCPar.seq).astype(int) 
        njobs=MCMCPar.seq
        start_engine_paralel='Yes'  
        
        if start_engine_paralel=='Yes':
            eng=start_engine(MCMCPar.seq,u) #initilize an engine at each core
            
        if start_engine_paralel=='No':
            eng=[] #empy engine as it is initialized at every iteration within the function run forward
                
            
### Main sampling loop  
       

        while self.MCMCVar.Iter < self.MCMCPar.ndraw:
            start_time1c = time.time()
            
            if (self.MCMCVar.Iter < self.MCMCPar.tune_phi):  
               
                if (prov_AR < self.MCMCPar.AR_min):
                                
                    for i in range(MCMCPar.phi.shape[0]):
                    
                        self.MCMCPar.phi[i] = MCMCPar.phi[i]*0.8  
                        
                        if self.MCMCPar.phi[i]>=self.MCMCPar.phi_max[i]:
                            self.MCMCPar.phi[i]=self.MCMCPar.phi_max[i]
                        elif self.MCMCPar.phi[i]<=self.MCMCPar.phi_min[i]:
                            self.MCMCPar.phi[i]=self.MCMCPar.phi_min[i]                        
                    
                if (prov_AR > self.MCMCPar.AR_max):
                                
                    for i in range(MCMCPar.phi.shape[0]):
                    
                        self.MCMCPar.phi[i] = MCMCPar.phi[i]*(1.2)
                    
                        if self.MCMCPar.phi[i]>=self.MCMCPar.phi_max[i]:
                            self.MCMCPar.phi[i]=self.MCMCPar.phi_max[i]
                        elif self.MCMCPar.phi[i]<=self.MCMCPar.phi_min[i]:
                            self.MCMCPar.phi[i]=self.MCMCPar.phi_min[i]
                        

            # Initialize totaccept
            totaccept = 0
         
            # Loop a number of K mcmc steps for each intermediate distribution

            for gen_number in range(0,self.MCMCPar.steps):
                
                time_complete_gen_start=time.time()
                # Update T
                self.MCMCVar.T = self.MCMCVar.T + 1
                
                for i in range(MCMCPar.phi.shape[0]):
                    phi_used[i]=MCMCPar.phi[i]      

                
                phi_ev=np.append(phi_ev,self.MCMCPar.phi.reshape(1,self.MCMCPar.seq),axis=0)
                
                phi_used_ev=np.append(phi_used_ev,phi_used.reshape(1,self.MCMCPar.seq),axis=0)                
                
                xold = np.array(self.MCMCVar.X[:self.MCMCPar.seq,:self.MCMCPar.n])
                log_p_xold = np.array(self.MCMCVar.X[:self.MCMCPar.seq,self.MCMCPar.n + 2-1])


                SetInput(MCMCPar.seq, MCMCPar.cond_type,phi_used,xold,'No',MCMCPar.dc_set,self.MCMCVar.T)
                

                if (np.random.rand(1) <= self.MCMCPar.parallelUpdate):
                    Update = 'Parallel_Direction_Update'
                else:
                    Update = 'Snooker_Update'

                              
                out=np.zeros((MCMCPar.seq,MCMCPar.nx*MCMCPar.ny))
                xnew=np.zeros((MCMCPar.seq,MCMCPar.nx,MCMCPar.ny))
                str0='con_' 
                
                u=np.linspace(1,MCMCPar.seq,MCMCPar.seq).astype(int)
                
                time10=time.time()

                njobs=MCMCPar.seq
                start_time_iteration= time.time()                
                result=ds_parallel(njobs,u,self.MCMCVar.T)
                


                for i in range(MCMCPar.seq):    
                    xnew[i,:,:]=np.reshape(result[i,:],(MCMCPar.nx,MCMCPar.ny),order='F')
                time20=time.time()
                print('time run DS',time20-time10)

       
                if start_engine_paralel=='No':
                
                    time1=time.time()
            
                    d = run_forward_par(njobs,xnew,u,eng,start_engine_paralel) 

                    time2=time.time()

                    print('time run_forw',time2-time1)


                if start_engine_paralel=='Yes':
                    
                    time3=time.time()
        
                    d=run_forward_par(njobs,xnew,u,eng,start_engine_paralel) 
        
                    time4=time.time()
        
                    print('time run_forw',time4-time3)

                end_time_iteration= time.time()
                print('time_DS+FS',end_time_iteration-start_time_iteration)

                
                # Compute the likelihood of each proposal in each chain
                of_xnew,log_p_xnew = CompLikelihood(xnew,d,self.MCMCPar,self.Measurement,self.Extra)

                # Calculate the Metropolis ratio
                accept = Metrop(self.MCMCPar,log_p_xnew,log_p_xold,Extra,beta)
                
                xnew=np.reshape(xnew,(self.MCMCPar.seq,self.MCMCPar.n))
                print('xnew-xold',np.sum(xnew-xold))
                # And update X and the model simulation
                idx_X= np.argwhere(accept==1);idx_X=idx_X[:,0]
  
                if not(idx_X.size==0):
                     
                    self.MCMCVar.X[idx_X,:] = np.concatenate((xnew[idx_X,:],of_xnew[idx_X,:],log_p_xnew[idx_X,:]),axis=1)
                    self.MCMCVar.Xfx[idx_X,:] = d[idx_X,:]
                
                                  
                # Check whether to add the current points to the chains or not?
                if self.MCMCVar.T == self.MCMCPar.thin:
                    # Store the current sample in Sequences
                    self.MCMCVar.iloc = self.MCMCVar.iloc + 1
                    aux=np.reshape(self.MCMCVar.X.T,(1,self.MCMCPar.n+2,self.MCMCPar.seq))
                    likelihood_ev=np.append(likelihood_ev,aux[0,self.MCMCPar.n+1,:])
                    self.Sequences[self.MCMCVar.iloc,0,:self.MCMCPar.seq] =aux[0,self.MCMCPar.n+1,:]
                   
                    # Check whether to store the simulation results of the function evaluations
                    if self.MCMCPar.savemodout==True:
                        self.fx=np.append(self.fx,self.MCMCVar.Xfx,axis=0)
                        # Update m_func
                        self.MCMCVar.m_func = self.MCMCVar.m_func + self.MCMCPar.seq
                    else:
                        self.MCMCVar.m_func=None
                        # And set the T to 0
                    self.MCMCVar.T = 0

                counter=counter+1.
                
                ### Save intermediate results
                
                if counter==2000.:
                    intermediate_models=np.append(intermediate_models,self.MCMCVar.X)
                    np.save('intermediate_weights.npy',weig_seq)
                    np.save('intermediate_models_run.npy',intermediate_models)                    
                    np.save('intermediate_beta.npy',beta_run)
                    np.save('intermediate_likelihood_ev',likelihood_ev)
                    np.save('intermediate_phi_ev',phi_used_ev)
                    np.save('intermediate_ESS',ESS_ev)                    
                    counter=0.
                    
                # Compute squared jumping distance for each CR value
                if (self.MCMCPar.Do_pCR==True and self.MCMCVar.Iter < 0.1 * self.MCMCPar.ndraw):
                   
                    # Calculate the standard deviation of each dimension of X
                    r = matlib.repmat(np.std(self.MCMCVar.X[:,:self.MCMCPar.n],axis=0),self.MCMCPar.seq,1)
                    # Compute the Euclidean distance between new X and old X
                    delta_normX = np.sum(np.power((xold[:,:self.MCMCPar.n] - self.MCMCVar.X[:,:self.MCMCPar.n])/r,2),axis=1)
                                        
                    # Use this information to update delta_tot which will be used to update the pCR values
                    self.MCMCVar.delta_tot = CalcDelta(self.MCMCPar.nCR,self.MCMCVar.delta_tot,delta_normX,self.MCMCVar.CR[:,gen_number])


                # Compute number of accepted moves
                totaccept = totaccept + np.sum(accept)

                # Update total number of MCMC iterations
                self.MCMCVar.Iter = self.MCMCVar.Iter + self.MCMCPar.seq
                
                time_complete_gen_end=time.time()
                print('time complete iteration',time_complete_gen_end-time_complete_gen_start)
             
            curr_log_lik=np.array(self.MCMCVar.X[:self.MCMCPar.seq,self.MCMCPar.n + 2-1])    

            # Store acceptance rate
            self.OutDiag.AR[self.MCMCVar.iteration-1,:] = np.concatenate((np.array([self.MCMCVar.Iter]).reshape((1,1)), np.array([100 * totaccept/(self.MCMCPar.steps * self.MCMCPar.seq)]).reshape((1,1))),axis=1)
            
            prov_AR=100 * totaccept/(self.MCMCPar.steps * self.MCMCPar.seq)
            
            # Store probability of individual crossover values
            self.OutDiag.CR[self.MCMCVar.iteration-1,:self.MCMCPar.nCR+1] = np.concatenate((np.array([self.MCMCVar.Iter]).reshape((1,1)), self.MCMCVar.pCR),axis=1)
            
            # Is pCR updating required?
            if (self.MCMCPar.Do_pCR==True): #and self.MCMCVar.Iter < 0.1 * self.MCMCPar.ndraw):

                # Update pCR values
                self.MCMCVar.pCR = AdaptpCR(self.MCMCPar.seq,self.MCMCVar.delta_tot,self.MCMCVar.lCR,self.MCMCVar.pCR)

            # Generate CR values from current pCR values
            self.MCMCVar.CR,lCRnew = GenCR(MCMCPar,self.MCMCVar.pCR); self.MCMCVar.lCR = self.MCMCVar.lCR + lCRnew

            # Calculate Gelman and Rubin Convergence Diagnostic
            start_idx = np.maximum(1,np.floor(0.5*self.MCMCVar.iloc)).astype('int64')-1; end_idx = self.MCMCVar.iloc
            
            # Update number of complete generation loops
            self.MCMCVar.iteration = self.MCMCVar.iteration + 1

            if self.MCMCPar.save_tmp_out==True:
                with open('out_tmp'+'.pkl','wb') as f:
                    pickle.dump({'Sequences':self.Sequences,
                    'OutDiag':self.OutDiag,'fx':self.fx,'MCMCPar':self.MCMCPar,
                    'MCMCVar':self.MCMCVar,'Measurement':self.Measurement,
                    'ModelName':self.ModelName,'Extra':self.Extra},f, protocol=pickle.HIGHEST_PROTOCOL)
            
            if beta>=1.:
                last_models=self.MCMCVar.X           
                last_likelihoods=curr_log_lik
                break
            

            next_beta,incr,CESS_found,omega_new = binary_search(curr_log_lik,CESSf,beta,norm_weight,self.MCMCPar.seq,self.MCMCPar.gamma,omega)            

            CESS_ev=np.append(CESS_ev,CESS_found)
                
### Calculate importance weights for current beta 
    
            contribution = np.exp((next_beta - beta) * curr_log_lik)
            
            weig_cont=np.append(weig_cont,contribution)
            
            new_weight = np.multiply(norm_weight,contribution)
            
            new_weight_ev=np.append(new_weight_ev,new_weight)
            
            ESS=(np.sum(norm_weight*contribution))**2 / np.sum(norm_weight**2*contribution**2)
            
            ESS_ev=np.append(ESS_ev,ESS)
            
            norm_weight = new_weight / np.sum(new_weight)
            
            weig_seq=np.append(weig_seq,norm_weight)
            
            weig_unn=np.multiply(weig_unn,contribution)
            
            weights_unnorm=np.append(weights_unnorm,weig_unn)
            
        
            evid=np.sum(new_weight)
            
            evid_log=np.log(evid)
                        
            evid_cont=np.append(evid_cont, evid_log)
            
            evid_evolution=evid_evolution+evid_log
            
            evid_ev=np.append(evid_ev, evid_evolution)
            
### Perform resampling if needed
            
            if (ESS/self.MCMCPar.seq < MCMCPar.ESSf_div):
                
                print('Resample')
                
                Xres, ind, eve = resampling(norm_weight,self.MCMCPar.seq,self.MCMCVar.X,self.MCMCPar.n, anc_prev,eve_prev)
                
                eve_seq=np.append(eve_seq,eve)
                
                anc_prev=ind
                
                eve_prev=eve
                
                self.MCMCVar.X = Xres
                
                norm_weight= (1 / self.MCMCPar.seq) * np.ones(self.MCMCPar.seq)
              
            else:
                
                eve_seq=np.append(eve_seq,eve_prev)
                          
            beta=next_beta
            omega=omega_new
            increment=np.append(increment,incr)
            beta_run=np.append(beta_run,beta)   
#                       
        
        start_time22 = time.time()
        # Remove zeros from pre-allocated variavbles if needed

        self.Sequences,self.OutDiag= Dreamzs_finalize(self.MCMCPar,self.Sequences,self.OutDiag,self.fx,self.MCMCVar.iteration,self.MCMCVar.iloc,self.MCMCVar.pCR)
       
        if self.MCMCPar.saveout==True:
            with open('dreamzs_out'+'.pkl','wb') as f:
                pickle.dump({'Sequences':self.Sequences,'OutDiag':self.OutDiag,'MCMCPar':self.MCMCPar,'Measurement':self.Measurement,'Extra':self.Extra},f
                , protocol=pickle.HIGHEST_PROTOCOL)
                
        end_time22 = time.time()
        print("This saving took %5.4f seconds." % (end_time22 - start_time22))

#        last_models=[]
#        last_likelihoods=[]
        self.Sequences=self.Sequences[1:,:,:]        
 #       last_models=[]
        eve_seq=eve_seq.reshape(int(eve_seq.shape[0]/self.MCMCPar.seq),self.MCMCPar.seq)
        
        return self.Sequences, self.OutDiag, self.MCMCPar, self.MCMCVar, beta_run, phi_used_ev, phi_ev, weig_seq, CESS_ev, increment, ESS_ev, evid_cont, evid_ev, weights_unnorm, new_weight_ev, weig_cont, eve_seq,last_models,last_likelihoods,intermediate_models
        
