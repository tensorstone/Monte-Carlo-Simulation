
# coding: utf-8

# In[4]:

import numpy as np
import matplotlib.pyplot as plt
import datetime
from copy import copy
import random
from multiprocessing.pool import Pool


# In[3]:

def mc_sim(seed):
    random.seed(seed)
    I =32
    J =32
    #IJ  are the shape of the simulated optical lattice
    A = np.empty((I,J))
    for i in range(I):
        for j in range(J):
            A[i,j]= 1

    JKBTMAX = 5 #range of temperature 
    Steplen = 0.025 # step length
    K = int(JKBTMAX/Steplen) +1 
    C = np.empty((I,J,K))
    for k in range(K):
        for i in range(I):
            for j in range(J):
                C[i,j,k]=0

    save = []
    qiguaile = []
    W =12560000
    for JKBT in range(K):
        Temp = K-JKBT 
        for w in range(W):
            i = int(np.round((I-1)*random.uniform(0,1)))
            j = int(np.round((J-1)*random.uniform(0,1)))
            A[i,j]=-A[i,j]
            if i==0:
                im=I-1
            else:
                im=i-1
            if i==I-1:
                ip=0
            else:
                ip=i+1
            if j==0:
                jm=J-1
            else:
                jm=j-1
            if j==J-1:
                jp=0
            else:
                jp=j+1


            dE =-2* (A[i,j]*A[ip,j]+A[i,j]*A[i,jp]+A[i,j]*A[im,j]+A[i,j]*A[i,jm])
            if dE<0:
                A[i,j] = A[i,j]
            else:
                if np.exp(-dE/(Temp*Steplen))<random.uniform(0,1):
                    A[i,j] = -A[i,j]

            if w>=W*0.9:
                C[:,:,JKBT] = C[:,:,JKBT]+A[:,:]
                if w%25600==0:
                    save.append(copy(A[:,:]))
        C[:,:,JKBT] = C[:,:,JKBT]/(W*0.1)            
        qiguaile.append(C[:,:,JKBT])

    # compute
    T = np.zeros((9849,))
    u=5.025
    for i in range(9849):
        if i%49==0:
            u = u- 0.025
        T[i]=u

    # save
    save = np.reshape(save,[-1,32*32])
    save = np.asarray(save)
    np.savetxt("IsingModel32_32simul_%d.csv"%seed, save, delimiter=',')
    np.savetxt("Temp32_32_%d.csv"%seed, T, delimiter=',')


# In[6]:

pool = Pool(processes=10)
params = range(10)
pool.map(mc_sim, params)


# In[ ]:



