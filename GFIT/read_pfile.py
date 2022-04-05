#!/usr/local/anaconda3/bin/python3
import os, sys
import numpy
from scipy.interpolate import interp1d
import subprocess
from gefit_tool import *
import eqdsk
from shutil import copyfile,rmtree
from exec_dirs import efit_dir
import matplotlib.pyplot as plt

def readp(filename):

	f = open(filename,'r')
	ndat = int(f.readline())
	line = f.readline()
	xx = np.zeros(ndat)
	te = np.zeros(ndat)
	ne = np.zeros(ndat)
	ti = np.zeros(ndat)
	ni = np.zeros(ndat)
	vt = np.zeros(ndat)

	for i in range(ndat):
		line = f.readline().split()
		xx[i] = float(line[0])
		te[i] = float(line[1])
		ne[i] = float(line[2])
		ti[i] = float(line[3])
		ni[i] = float(line[4])
		if len(line) > 5:
			vt[i] = float(line[5])
			isvt = True
		else: isvt = False

	f.close()

	return xx,te,ne,ti,ni,vt,isvt


datlen = len(sys.argv)
cind = ['r','b','orange','magenta','g','gray']
fig,[ax1,ax2,ax3] = plt.subplots(1,3,figsize=(15.,5.))
slegend1 = []; slegend2 = []; slegend3 = [];
plegend1 = []; plegend2 = []; plegend3 = [];

for i in range(datlen-1):

	filename = sys.argv[i+1]
	print('>>> #%i -> %s'%(i+1,filename))
	xx,te,ne,ti,ni,vt,isvt = readp(filename)

	line1, = ax1.plot(xx,te,color=cind[i])
	line2, = ax1.plot(xx,ti,'--',color=cind[i])
	line3, = ax2.plot(xx,ne,color=cind[i])
	line4, = ax2.plot(xx,ni,'--',color=cind[i])

	if i==0: 
		slegend1.append(line1)
		slegend1.append(line2)
		plegend1.append('Te [keV]')
		plegend1.append('Ti [keV]')
		slegend2.append(line3)
		slegend2.append(line4)
		plegend2.append('Ne [$10^19 m^{-3}$]')
		plegend2.append('Ni [$10^19 m^{-3}$]')
	slegend1.append(line1)
	slegend2.append(line2)
	plegend1.append('#%i'%(i+1))
	plegend2.append('#%i'%(i+1))

	if isvt: 
		line5, = ax3.plot(xx,vt,color=cind[i])
		slegend3.append(line5)
		plegend3.append('#%i'%(i+1))

ax1.set_title('T [keV]')
ax2.set_title('ne [$10^19 m^{-3}$]')
ax3.set_title('$V_\phi$ [km/s]')

ax1.set_xlabel('$\psi_N$')
ax2.set_xlabel('$\psi_N$')
ax3.set_xlabel('$\psi_N$')

ax1.legend(slegend1,plegend1)
ax2.legend(slegend2,plegend2)
if isvt:
	ax3.legend(slegend3,plegend3)

fig.tight_layout()
plt.show(block=False)
input()



