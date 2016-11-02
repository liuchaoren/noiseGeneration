#!/et/cl221/usr/local/bin/python2.7
# colorgau.py # a colored gaussian random number, featuring a specified correlation function

import numpy as np
#import matplotlib.pyplot as plt
import sys
sys.path.append('/home/dbchem/cl221/lib/mypy')
#from scipy.fftpack import fft
from scipy.fftpack import ifft
#from cmath import sqrt as csqrt
from math import sqrt
from myfunc import fftshiftf
#import time


#start = time.time()

# open files
timeserial=open('randseq','w')
# define the constants
# correlation function: (epsilon/tau)*exp(-t/tau)
tau = 1*2.
#epsilon = 0.16**2*tau
epsilon = 1*2.
dt = 0.01
wlen = 2**19
dw = 1/(wlen*2*dt)

pi=np.pi
# rwlist is the fft of correlation function
rwlist = []
for i in range(-wlen,wlen+1):
    w = i*dw
    rwlist.append((1/dt)*2*epsilon/(1+(2*pi*tau*w)**2))
###########################################
a=[0]*(2*wlen+1)
b=[0]*(2*wlen+1)
a[wlen]=np.random.randn()
b[wlen]=0

for i in range(wlen):
    a[i]=sqrt(0.5)*np.random.randn()
    a[2*wlen-i]=a[i]
    b[i]=sqrt(0.5)*np.random.randn()
    b[2*wlen-i]=-b[i]
alpha = [m+1j*n for (m,n) in zip(a,b)]
###########################################
etafreq=[sqrt(i*2*wlen)*j for (i,j) in zip(rwlist,alpha)]
etafreq = fftshiftf(etafreq)

rtlist = ifft(etafreq)
# output
for i in rtlist:
    timeserial.write('%s\n' % (i.real))

#print time.time()-start
#################################################
#rtlist = []
#for i in range(-wlen,wlen+1):
#    rt = etafreq[wlen]
#    for j in range(wlen+1,2*wlen+1):
#        rt = rt + 2*(etafreq[j]*np.exp(2j*np.pi*i*(j-wlen)/(2*wlen))).real
#    rtlist.append(rt)
#    
#for i in rtlist:
#    timeserial.write('%s\n' % (i.real))

#etatime=[i.real for i in etatime]
#for i in etatime:
#    timeserial.write('%s\n' %i)    

#plt.hist(etatime,100)
#plt.savefig('tmp.png')

