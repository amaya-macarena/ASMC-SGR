# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 09:17:33 2021

@author: mamaya
"""

import numpy as np
import sys
import random
import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.pyplot as plt
from scipy.io import loadmat
#import os
import matlab.engine
#from oct2py import octave
import time
#import pymatlab

def run_forward(model,eng):
    
    T=model.shape[0]
#    K1=loadmat('trueK_2.mat')
#    K2=dict.items(K1)
#    K3=list(K2)
#    K4=np.array(K3)
#    K=K4[3,1]
    K=model.astype(float)
    K[K==1.]=1e-2
    K[K==0.]=1e-4
 #   model=matlab.double(K.tolist())

    
    

    
    wells=np.array([3,10,17,24,31,38,45,52,59,66,73])
    idat=11
    d=np.zeros((T,30,idat))
    
    for i in range (T):
#        cond=model[i,:,:]
#        cond=model
#        cond[cond==1.]=1e-2
#        cond[cond==0.]=1e-4
#        [~, conc]=MaFloT(COND)
        time1=time.time()
        A=matlab.double(K[i,:,:].tolist())
        time2=time.time()
        conc=eng.MaFloT(A)
        conc_np=np.array(conc)
        btc=np.zeros((30,idat))
        print('time2-time1',time2-time1)

        for j in range(idat):
            btc[:,j]=conc_np[wells[j],51,:]
        
        d[i,:,:]=btc
#        
    return d
        
        
        