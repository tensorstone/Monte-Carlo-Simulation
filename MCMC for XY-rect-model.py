
# coding: utf-8

# In[129]:

import numpy as np
import matplotlib.pyplot as plt
import datetime
from copy import copy

I =16
J =16
#IJ  are the shape of the simulated optical lattice
A = np.empty((I,J))
for i in range(I):
    for j in range(J):
        A[i,j]= 0
B = np.empty((I,J))
for i in range(I):
    for j in range(J):
        B[i,j]= 2*np.pi*random.random()-np.pi
        
        
JKBTMAX = 5 #range of temperature 
Steplen = 0.025 # step length
K = int(JKBTMAX/Steplen) +1 
C = np.empty((I,J,K))
for k in range(K):
    for i in range(I):
        for j in range(J):
            C[i,j,k]=0

import random
save = []
qiguaile = []
W =5120000
for JKBT in range(K):
    Temp = K-JKBT 
    for w in range(W):
        i = int(np.round((I-1)*random.uniform(0,1)))
        j = int(np.round((J-1)*random.uniform(0,1)))
        B[i,j]= 2*np.pi*random.random()-np.pi
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
        
        
        dE =-(np.cos(B[i,j]-A[ip,j])+np.cos(B[i,j]-A[i,jp])+np.cos(B[i,j]-A[im,j])+np.cos(B[i,j]-A[i,jm]))+(np.cos(A[i,j]-A[ip,j])+np.cos(A[i,j]-A[i,jp])+np.cos(A[i,j]-A[im,j])+np.cos(A[i,j]-A[i,jm]))
        if dE<0:
            A[i,j] = B[i,j]
        else:
            if np.exp(-dE/(Temp*Steplen))<random.uniform(0,1):
                A[i,j] = A[i,j]
            else:
                A[i,j] = B[i,j]
        if w>=W*0.9:
            C[:,:,JKBT] = C[:,:,JKBT]+A[:,:]
            if w%2560==0:
                save.append(copy(A[:,:]))
    C[:,:,JKBT] = C[:,:,JKBT]/(W*0.1)            
    qiguaile.append(C[:,:,JKBT])


# In[132]:

M=[]
for i in range(K):
    M.append(np.sum(C[:,:,i]))
M = np.asarray(M)


# In[133]:

plt.plot([5-0.025*i for i in range(K)],np.abs(M))
plt.show()


# In[134]:

np.shape(save)


# In[135]:

save_arr = np.asarray(save)


# In[136]:

for i in range(save_arr.shape[0]):
    for j in range(save_arr.shape[1]):
        for k in range(save_arr.shape[2]):
            save_arr[i,j,k] = np.cos(save_arr[i,j,k])


# In[144]:

np.shape(save)


# In[143]:

num_0 = 40198
plt.imshow(save_arr[num_0])
plt.show()
save_arr[num_0]


# In[145]:

T = np.zeros((40200,))
u=5.025
for i in range(40200):
    if i%200==0:
        u = u- 0.025
    T[i]=u
    
    
save = np.reshape(save,[-1,16*16])
save = np.asarray(save)
np.savetxt("rect-XYModel16_16simul.csv",save,delimiter=',')
np.savetxt("rect-Temp16_16.csv",T,delimiter=',')


# In[86]:

plt.imshow(save[20000])
plt.show()


# In[ ]:



