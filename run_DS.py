# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 09:03:59 2021

@author: mamaya
"""

import numpy as np
from subprocess import call
import matplotlib.pyplot as plt

def runDS(idx,str0,T):
    
    str1='./deesse'
    str1_b = 'test_'+str0+str(idx)+'.in'
    str2='test_'+str(idx)+'_simul_'+str0+'real00000.gslib'

    call([str1, str1_b])
    
    data = np.loadtxt(str2,skiprows=3)
    data2=np.reshape(data,(75,101),order='F')
    plt.imshow(data2, interpolation='none', aspect='equal', cmap='jet')
    plt.savefig('./after_resim_it'+str(T)+'_chain'+str(idx)+'.png', dpi=300)
    
#    np.save('data'+str(idx)+'T_'+str(T)+'.npy',data)
    return data