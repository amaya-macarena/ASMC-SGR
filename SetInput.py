# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 08:59:12 2021

@author: mamaya
"""

import numpy as np
import shutil
import fileinput
import sys
import GenCond
from GenCond import GenCond


def SetInput(nc,cond_type,phi,m_current,Init,direct_cond,T):
    
        
    
#    seed=np.zeros(nc)
#    for i in range(nc):
#        seed[i]=1+int(np.random.uniform()*(1e6-100))    
# Now create the conditioning dataset
    for i in range (1,nc+1):
        
        a_file = open('test_con_'+str(i)+'.in')
        lines_to_read = [323]

        for position, line in enumerate(a_file):
            if position in lines_to_read:
                a=int(line)
#                print(line)
        
#        new_seed=1+int(np.random.uniform()*(1e6-100)) 
        new_seed=int(np.random.uniform(low=1.0,high=2.0)*(1e6-100))
        str_name='test_con_'+str(i)+'.in' 
        fin = open(str_name,'rt')
        data=fin.read()        
#        data = data.replace('test_simul', 'test_'+str(i)+'_simul_con')
#        data = data.replace('test_report_con.txt', 'test_report_'+str(i)+'.txt')        
#        data = data.replace('cond_.dat', 'cond_'+str(i)+'.dat') 
#        data = data.replace('cond_1.dat', 'cond_'+str(i)+'_T'+str(T)+'.dat')
        data = data.replace(''+str(a)+'',''+str(new_seed)+'') 
        fin = open(str_name,'wt')
        fin.write(data)
        fin.close()
        
        GenCond(i,cond_type,phi[i-1],m_current[i-1,:],direct_cond,T)

    return []