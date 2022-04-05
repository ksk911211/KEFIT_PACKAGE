#!/usr/local/anaconda3/bin/python3
import os,sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import matplotlib as mpl
from exec_dirs import gzip_dir

nfile = sys.argv[1]; tmin = sys.argv[2]; tmax = sys.argv[3]; tfile = sys.argv[4]
if (not os.path.isfile(nfile) and not os.path.isfile(nfile+'.gz')): exit()
if not os.path.isfile(nfile):
	comm=gzip_dir+' -d '+nfile+'.gz'
	os.system(comm)
npzfile  = np.load(nfile, mmap_mode='r')
da = [];
for i in range(2): da.append(npzfile['arr_%i'%i])
comm=gzip_dir+' '+nfile
os.system(comm)

peaks = [];
if os.path.isfile(tfile):
	npzfile = np.load(tfile,mmap_mode='r')
	peaks = npzfile['arr_0']

ind1 = np.where(da[0]>=float(tmin)/1.e3);
ind2 = np.where(da[0][ind1]<=float(tmax)/1.e3);
plt.figure('ELM peaks')
plt.plot(da[0][ind1][ind2]*1.e3,da[1][ind1][ind2],linewidth=1)
for time in peaks:
	plt.axvline(x=time,linestyle='--',color='gold')
plt.xlabel('time [ms]')
plt.ylabel('[a.u.]')
plt.show()
