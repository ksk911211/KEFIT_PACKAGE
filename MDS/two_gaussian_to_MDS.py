import os,sys
import numpy as np
import pickle
from scipy.interpolate import interp1d, interp2d
import matplotlib.pyplot as plt

with open('twoCES@31189.txt','r') as f:

	for i in range(16):
		line = f.readline()
		if line.find('DimSize')>-1:
			[ntime,nradial] = np.array(line.split('=')[1].split(','),dtype='int')

		if line.find('ValNo')>-1:
			nval = int(line.split('=')[1])

		if line.find('ValName')>-1:

			dat  = {}
			dat['times'] = np.zeros(ntime)
			dat['radius'] = np.zeros(nradial)
			vals = []
			line2= line.split('=')[1].split("'")
			
			for i in range(nval):
				ind = 1+2*i
				dat[line2[ind]] = {}
				vals.append(line2[ind])

				for j in range(nradial):
					dat[line2[ind]][j] = np.zeros(ntime)

	line = f.readline()
	if line.find('[data]')>-1:
		for i in range(ntime):
			for j in range(nradial):
				line = f.readline()
				if not line: break
				line2= np.array(line.split(','),dtype='float')
				dat['times'][i] = line2[0]
				dat['radius'][j] = line2[1]

				for k in range(nval):
					dat[vals[k]][j][i] = line2[k+2]

def _get_save_zipf(nfile,data=[],ndata=4):
	gzip_dir = 'gzip'
	np.savez(nfile,data[0],data[1],data[2],data[3])
	comm=gzip_dir+' '+nfile
	os.system(comm)
	if os.path.isfile(nfile): os.system('rm -f %s'%nfile)
	return  

ces_t = {}
ces_v = {}
if not os.path.isdir('MDS'): os.mkdir('MDS')
for ch in range(dat['radius'].shape[0]):
	nfile_t = 'MDS/%s%i.npz'%('TI',ch+1)
	nfile_v = 'MDS/%s%i.npz'%('VT',ch+1)
	_get_save_zipf(nfile_t,[dat['times'],dat['Ti'][ch]*1.e+3,dat['Ti'][ch]*1.e+2,dat['radius'][ch]*1.e3])
	_get_save_zipf(nfile_v,[dat['times'],dat['Vc'][ch]*1.e+0,dat['Vc'][ch]*1.e-1,dat['radius'][ch]*1.e3])

	nfile_t = 'MDS/%s_size'%('TI')
	nfile_v = 'MDS/%s_size'%('VT')

with open(nfile_t,'w') as f: f.write('%i %i %i %i\n'%(dat['times'].shape[0],28000,dat['radius'].shape[0],dat['radius'].shape[0]))
with open(nfile_v,'w') as f: f.write('%i %i %i %i\n'%(dat['times'].shape[0],28000,dat['radius'].shape[0],dat['radius'].shape[0]))

len(self.mds['da'][0]),