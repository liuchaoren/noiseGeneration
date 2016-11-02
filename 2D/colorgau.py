#!/bin/python
# colorgau.py
# a colored gaussian random number, featuring a specified correlation function

import numpy as np
#import matplotlib.pyplot as plt
import sys
sys.path.append('/home/dbchem/cl221/lib/mypy')
#from scipy.fftpack import fft
from scipy.fftpack import ifft
#from cmath import sqrt as csqrt
from math import sqrt
from myfunc import fftshiftf 
import time

start = time.time()
# open files
timeserial = []
for n in range(1,9):
    timeserial.append(open('randseq'+'%s' %n,'w'))
# define the constants
# temporal
# correlation function: (epsilon/tau)*exp(-t/tau)
epsilontime = 1*6.
tautime = 1*6.
dt = 0.01
wlen = 2**15
dw = 1/(wlen*2*dt)
# spatial
# correlation function: exp(-a*x^2)
#a=1
epsilonspace=1*1.
tauspace=1*1.
dx = 0.01
plen = 2**13
dp = 1/(plen*2*dx)
# pick up lines
xjump = 10
xlines = 8

pi=np.pi

# store the fft of time dim
rwlist = []
for i in range(-wlen,wlen+1):
    w = i*dw
    rwlist.append((1/dt)*2*epsilontime/(1+(2*pi*tautime*w)**2))
rwlist = fftshiftf(rwlist)

# store the fft of spatio dim
rplist = []
for i in range(-plen,plen+1):
    p = i*dp
#    rplist.append((1/dx)*sqrt(pi/a)*np.exp(-(p*pi)**2/a))
    rplist.append((1/dx)*2*epsilonspace/(1+(2*pi*tauspace*p)**2))
rplist = fftshiftf(rplist)

wseries = [[],[],[],[],[],[],[],[]]
 
# ifft of first column 
alpha = np.zeros(2*plen+1,dtype='complex')
alpha[0]=np.random.randn()
for i in range(1,plen+1):
    real = sqrt(0.5)*np.random.randn()
    imag = sqrt(0.5)*np.random.randn()
    alpha[i]=real+1j*imag
    alpha[2*plen+1-i]=real-1j*imag
wptmp = [sqrt(2*wlen*2*plen*m*rwlist[0]) for m in rplist]
wcolumn=[i*j for (i,j) in zip(wptmp,alpha)]
wcolumn = ifft(wcolumn).tolist()

for n in range(8):
    wseries[n].append(wcolumn[1000*n:1000*n+xlines*xjump:xjump])

# ifft of other columns
for i in range(1,wlen+1):
    alpha = np.zeros(2*plen+1,dtype='complex')
    wptmp = [sqrt(2*wlen*2*plen*m*rwlist[i]) for m in rplist]
    for j in range(2*plen+1):
        real = sqrt(0.5)*np.random.randn()
        imag = sqrt(0.5)*np.random.randn()
        alpha[j]=real+1j*imag
    wcolumn = [i*j for (i,j) in zip(alpha,wptmp)]
    wcolumn = ifft(wcolumn).tolist()
    for n in range(8):
        wseries[n].append(wcolumn[1000*n:1000*n+xlines*xjump:xjump])

for n in range(8):
    wseries[n]=(np.array(wseries[n])).transpose()
    wseries[n]=wseries[n].tolist()
#wseries = (np.array(wseries)).transpose()
#wseries = wseries.tolist()
# ifft of the rows
for n in range(8):
    for i in range(xlines):
        wseries[n][i]=wseries[n][i]+[m.conjugate() for m in ((wseries[n][i])[1:])[-1::-1]]
        wseries[n][i] = ifft(wseries[n][i]).tolist()
        for j in range(2*wlen+1):
            timeserial[n].write('%s  ' % (wseries[n][i])[j].real)
        timeserial[n].write('\n')

print time.time()-start
