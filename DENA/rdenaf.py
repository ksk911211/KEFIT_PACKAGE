#!/usr/local/anaconda3/bin/python3
import os,sys
import matplotlib.pyplot as plt
from matplotlib import cm
import pylab
from collections import OrderedDict
import numpy as np
import pickle

cmaps = pylab.get_cmap('Oranges')

xloc = np.linspace(0.3,0.91,8)
f = open(sys.argv[1],'rb')
a=pickle.load(f);
f.close()
time = a['times']; lent = len(time)
print('>>> TMIN > ',time[0],' TMAX > ',time[-1],' TDEL >',time[1]-time[0],'ms')
fitx = a['fit']['ne'][time[0]][0]
xind = []
try: tstart = int(sys.argv[2]); tend = int(sys.argv[3]); delt = int(sys.argv[4])
except: print('>>> Type [tstart] [tend] [delt] in ms'); exit()
ind1 = np.where(time>=tstart);
ind2 = np.where(time[ind1]<=tend);
for xx in xloc:
	xind.append(np.argmin(abs(fitx-xx)));
plt.figure('1D trace')
llegend = []
for ind in xind:
	yy = np.zeros(lent);
	for j in range(lent):
		yy[j] = a['fit']['ne'][time[j]][1][ind]
		yr    = a['infile']['tciv'][time[j]][-1]
		yf    = a['infile']['tcif'][time[j]][-1]
		ys    = a['infile']['tcis'][time[j]][-1]
		if ind==xind[-1]: print(time[j],yy[j],yr,yf,ys)
	plt.plot(time[ind1][ind2],yy[ind1][ind2]-0.*yy[ind1][ind2][0],'--x')
	llegend.append('%3.2f'%fitx[ind])

plt.legend(llegend)
plt.xlabel('time [ms]'); plt.ylabel('density')

plt.figure('2D trace')
llegend = []
lent2 = round((len(time[ind1][ind2])+4)*(time[2]-time[1])/delt)
count = 1; prev_t = -1; 
for ind in range(lent):
	if (time[ind]<tstart or time[ind]>tend): continue
	if time[ind] >= (prev_t+delt): prev_t = time[ind]
	else: continue
	color = cmaps(1.*count/lent2)
	count += 1
	plt.plot(fitx,a['fit']['ne'][time[ind]][1],c=color)
	llegend.append('%i ms'%time[ind])

plt.xlabel('\\psi_N'); plt.ylabel('density')
plt.legend(llegend)

plt.figure('RAW vs FIT')
llegend = []; slegend = [];
for i in range(5):
	yy = np.zeros(lent); yy2 = np.zeros(lent); yy3 = np.zeros(lent)
	for j in range(lent): 
		yy[j] = a['fit']['nel'][time[j]][i+2]
		yy2[j]=a['infile']['tciv'][time[j]][i+2]
		yy3[j]=a['infile']['tcis'][time[j]][i+2]
	line1, = plt.plot(time[ind1][ind2],yy[ind1][ind2])
	line2  = plt.errorbar(time[ind1][ind2],yy2[ind1][ind2],yerr=yy3[ind1][ind2])
	llegend.append('FIT%02i'%(i+1)); llegend.append('TCI%02i'%(i+1));
	slegend.append(line1); slegend.append(line2)
plt.legend(slegend,llegend)
plt.xlabel('time [ms]'); plt.ylabel('density')
plt.show()
