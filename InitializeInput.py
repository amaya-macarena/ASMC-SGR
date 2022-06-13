# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 12:58:56 2021

@author: mamaya
"""
import numpy as np
import shutil
import fileinput
import sys

def InitializeInput(nc,Init):
     
    dummy=1
    
    if (Init==True): #create nc copies of the original input file that handle conditioning data
        
        # Input file(s) used for conditional simulation
        str1_con='test_con.in'
        
        for i in range (1,nc+1):
        
            str2='test_con_'+str(i)+'.in'
#            if exist(str2,'file');
#                delete(str2);
            shutil.copy(str1_con,str2)
            
            
        # Input file(s) used for unconditional simulation
        str1_unc='test_unc.in'
        
        for i in range (1,nc+1):
        
            str2='test_unc_'+str(i)+'.in'
#            if exist(str2,'file');
#                delete(str2);
            shutil.copy(str1_unc,str2)       
        
# Loop over input file(s) used for conditional simulation  
    for i in range (1,nc+1):
        
        seed=1+int(np.random.uniform()*(1e6-100))
        str_name='test_con_'+str(i)+'.in' 
        fin = open(str_name,'rt')
        data=fin.read()        
        data = data.replace('test_simul', 'test_'+str(i)+'_simul_con')
        data = data.replace('test_report_con.txt', 'test_report_'+str(i)+'.txt')        
        data = data.replace('cond_.dat', 'cond_'+str(i)+'.dat') 
        data = data.replace('444',''+str(seed)+'') 
        fin = open(str_name,'wt')
        fin.write(data)
        fin.close()

# Loop over input file(s) used for unconditional simulation  
    for i in range (1,nc+1):
        
        seed=1000+int(np.random.uniform()*(1e7-1000))
        str_name='test_unc_'+str(i)+'.in' 
        fin = open(str_name,'rt')
        data=fin.read()        
        data = data.replace('test_simul', 'test_'+str(i)+'_simul_unc') 
        data = data.replace('444',''+str(seed)+'')
        fin = open(str_name,'wt')
        fin.write(data)
        fin.close()        
        
    return []
        

            
            