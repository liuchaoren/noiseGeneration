#!/et/cl221/usr/local/bin/python2.7
# colorgau.py # a colored gaussian random number, featuring a specified correlation function

import numpy as np
from numpy.random import multivariate_normal
from scipy.sparse import diags
#import matplotlib.pyplot as plt
import sys
sys.path.append('/dscrhome/cl221/lib/mypy')
#from scipy.fftpack import fft
from scipy.fftpack import ifft
#from cmath import sqrt as csqrt
from math import sqrt
from myfunc import fftshiftf
#import time


#start = time.time()

# define the constants
# correlation function: (epsilon/tau)*exp(-t/tau)
tau = 1*2.
epsilon = 1*2.
dt = 0.01
wlen = 2**17
dw = 1/(wlen*2*dt)

# correlation matrix, exponential decay
# nearest neighbor correlation strength
NNcorr = 0.7 # correlation strength between nearest neighbor sites
length=27 # number of sites
diagnals = []
for i in range(length-1,0,-1)+range(0,length):
    diagnals.append([NNcorr**i]*(length-i))
cov = diags(diagnals,range(-(length-1),length)).todense()
mean = np.array([0.]*length)

pi=np.pi
# rwlist is the fft of correlation function
rwlist = []
for i in range(-wlen,wlen+1):
    w = i*dw
    rwlist.append((1/dt)*2*epsilon/(1+(2*pi*tau*w)**2))
###########################################
a=np.zeros([length,2*wlen+1])
b=np.zeros([length,2*wlen+1])
a[:,wlen]=multivariate_normal(mean,cov)
b[:,wlen]=np.array([0.]*length)

for i in range(wlen):
    a[:,i]=sqrt(0.5)*multivariate_normal(mean,cov)
    a[:,2*wlen-i]=a[:,i]
    b[:,i]=sqrt(0.5)*multivariate_normal(mean,cov)
    b[:,2*wlen-i]=-b[:,i]
alpha = a+1j*b
###########################################
etafreq=np.sqrt(np.array(rwlist)*2*wlen)*alpha
#etafreq=[sqrt(i*2*wlen)*j for (i,j) in zip(rwlist,alpha)]
etafreq = fftshiftf(etafreq)

timeserials = np.zeros([length,2*wlen+1],dtype='complex128')
for i in range(length):
    timeserials[i] = ifft(etafreq[i])
# output
np.savetxt('randseq',np.real(timeserials))
