# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 09:24:30 2021

@author: mamaya
"""

import numpy as np
import sys
import random
import matplotlib.pyplot as plt

def GenCond(idx,cond_type,phi,D,direct_cond,T):
    
    
#    direct_cond=np.array([[3,51,1,0],[10,51,1,1],[17,51,1,1],[24,51,1,0],[31,51,1,0],[38,51,1,1],[45,51,1,0],[52,51,1,0],[59,51,1,0],[66,51,1,0],[73,51,1,0]])
#    phi=10
    nx=75
    ny=101
    nz=1
    cond_type=2
    
    D=np.reshape(D,(nx,ny))

#    print('D.shape',D.shape)
#    [nx,ny,nz]=D.shape
#    print('[nx,ny,nz]',[nx,ny,nz])
#    print('([nx,ny,nz])')
    x=np.linspace(1,nx,nx)
    y=np.linspace(1,ny,ny)
    z=np.linspace(1,nz,nz)
    [yy,xx,zz]=np.meshgrid(y,x,z)
    
    lim=np.empty(2)
    pos=np.empty(3)
    if cond_type==2:
    # Box type selection, taken from the SIPPI toolbox of Hansen et al. (2013) - works herein for 2D images only   
#        phi=5
        lim[0]=phi
        lim[1]=phi
         # phi defines half the side length of the box that is resimulated (the larger phi, the smaller the number of conditioning data, phi is actually half the box side length)   
        
        pos[0]=min(x)+random.uniform(0,1)*(max(x)-min(x))
        pos[1]=min(y)+random.uniform(0,1)*(max(y)-min(y))
        pos[2]=min(z)+random.uniform(0,1)*(max(z)-min(z))
        
#        print('pos',pos[0],pos[1],pos[2])
        
        used=np.ones(np.shape(xx))
        
        ind=np.zeros(xx.shape)
        ind[:,:,0]=xx[:,:,0]
        ind[abs(ind-pos[0])<lim[0]]=0
#        ind[abs(ind-pos[0])>=lim[0]]=1
        
        ind2=np.zeros(yy.shape)
        ind2[:,:,0]=yy[:,:,0]
        ind2[abs(ind2-pos[1])<lim[1]]=0
#        ind2[abs(ind2-pos[1])>=lim[1]]=1
        
        suma=ind+ind2
#        print('suma.shape',suma.shape)
        out_box=np.argwhere(suma!=0)
        
#        print('out_box',out_box)
#        print('out_box_shape',out_box.shape)
#        print(T)
        
        # Classical selection (appears to induce some model degradation due to the interaction 
        # of very small phi values with likelihood maximimzation)
        
        
#        plt.figure(idx+T);plt.imshow(ind[:,:,0], interpolation='none', aspect='equal', cmap='jet');plt.colorbar()

#        plt.figure(idx+T);plt.imshow(ind2[:,:,0], interpolation='none', aspect='equal', cmap='jet');plt.colorbar()
    
#        plt.figure(T);plt.imshow(suma[:,:,0], interpolation='none', aspect='equal', cmap='jet');plt.colorbar()
    
        if out_box.size==0: #Make sure there is a least one conditioning datum, otherwise the MPS code might complain
            print('empty conditining')
            a=int(np.random.uniform(0,nx-1))
            b=int(np.random.uniform(0,ny-1))
            out_box=np.array([[a,b,0]])
    
        d_cond=np.empty((out_box.shape[0],4))
        D_new=np.ones((75,101))*4
        
        for i in range (out_box.shape[0]):
#            print(i)
            a=out_box[i,0]
            b=out_box[i,1]
            c=out_box[i,2]
            d_cond[i,:]=[int(xx[a,b,c]),int(yy[a,b,c]),int(zz[a,b,c]), D[a,b]]
            D_new[a,b]=D[a,b]
#            plt.figure(7)
#            plt.imshow(D[a,b], interpolation='none', aspect='equal', cmap='jet')
        
 #       plt.figure(6);
 #       plt.plot(d_cond[:,0],d_cond[:,1],'or')
 #       plt.axis('equal')
        
 #       plt.figure(7);
 #       plt.imshow(D_new, interpolation='none', aspect='equal', cmap='jet')
 #       print('d_cond.shape',d_cond.shape)
 #       print('d_cond',d_cond)
 #       print('direct_cond.shape',direct_cond.shape)
  #      print('direct_cond',direct_cond)
     #   plt.savefig('./check_it'+str(T)+'_chain'+str(idx)+'.png', dpi=300)

#        print('direct_cond.shape',direct_cond.shape)
        # First find potential duplicates
        rows_to_del=np.array([0])
               
        for j in range(direct_cond.shape[0]):
            
            for u in range(d_cond.shape[0]):
                if (d_cond[u,0]==direct_cond[j,0] and d_cond[u,1]==direct_cond[j,1]):
#                    print('yes')
#                    print(int(u))
                    rows_to_del=np.append(rows_to_del,int(u))
        rows_to_del.astype(int)
        if rows_to_del.size==1:
            rows_to_del=[]
        else:
            rows_to_del=rows_to_del[1:]
#        print('rows_to_del_idx',rows_to_del,idx)
#        print('actual rows',d_cond[rows_to_del,:])
#        print('rows_to_del.size_idx',rows_to_del.size,idx)            
#        rows_to_del.astype(int)
        d_cond_del=np.delete(d_cond,rows_to_del,0)
 #       print('d_cond_del.shape',d_cond_del.shape)
 #       print('d_cond_del',d_cond_del)
        d_cond=np.append(d_cond_del,direct_cond,axis=0)
 #       print('d_cond_append.shape',d_cond.shape)
 #       print('d_cond_append',d_cond)
        
        test=np.ones((75,101))*4
#        for p in range(75):
#            for q in range(101):
#                ext=np.where(d_cond[:,0])==float(p) and d_cond[:,1]==float(q))
#                print('ext',ext)
#                test[p,q]=d_cond[ext,3]               
        for t in range(d_cond.shape[0]):
            a=int(d_cond[t,0])
            b=int(d_cond[t,1])
            test[a-1,b-1]=d_cond[t,3]
      
#        plt.imshow(test, interpolation='none', aspect='equal', cmap='jet')  
#        plt.savefig('./test_it'+str(T)+'_chain'+str(idx)+'.png', dpi=300)
        #Now replace the conditioning data in the files
#        tit='cond_'+str(idx)+'_T'+str(T)+'.dat'
        tit='cond_'+str(idx)+'.dat'
        nvar=1
        line01=str(d_cond.shape[0])
        line02=str(nvar+3)
        line03='X'
        line04='Y'
        line05='Z'
        line06='facies'
        fid = open(tit,'w')
        fid.write(line01)
        fid.write("\n")
        fid.write(line02)
        fid.write("\n")
        fid.write(line03)
        fid.write("\n")
        fid.write(line04)
        fid.write("\n")
        fid.write(line05)
        fid.write("\n")
        fid.write(line06)
        fid.write("\n")
        np.savetxt(fid,d_cond,fmt='%e')
        
#        np.save('d_cond_part'+str(idx)+'_T'+str(T)+'.npy',d_cond)
#        for i in range(d_cond.shape[0]):
#            
##            fid.write(str(d_cond[i,:]))
#            np.savetxt(fid,d_cond[i,:],fmt='%.2f')

        
        fid.close()
        
      
    
    return d_cond