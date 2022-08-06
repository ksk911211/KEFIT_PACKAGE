#!/usr/local/anaconda3/bin/python3
import os,sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
# Color names are: https://matplotlib.org/3.1.0/gallery/color/named_colors.html
# Directories to draw the plots
dir1  = '/Users/sk42/Works/2022_KSTAR_EXP/n12_RMP/PROFILES/6100/4/PROFILE/FITPROF/SPROFILES'
dir2  = '/Users/sk42/Works/2022_KSTAR_EXP/n12_RMP/PROFILES/8780/8/PROFILE/FITPROF/SPROFILES'
dir3  = '/Users/sk42/Works/2022_KSTAR_EXP/n12_RMP/PROFILES/10380/2/PROFILE/FITPROF/SPROFILES'
dirs  = [dir1,dir2,dir3]

labels= ['ELMy','n=1 RMP','n=2 RMP']  #
pcolor= ['black','orangered','dodgerblue','gold','darkgreen'] #colors for each plots
dodraw= [True,True,True]

# Plot channel setups ('te,''ne','ti','vt')
draw_ch = [True,True,True,True]     # draw or not
xmap    = [0,0,0,0]                 # 0(psi_norm),1(tor_rho_norm)
plot_ra = ['raw','raw','raw','raw'] # raw(scatter) or averaged
x_min   = [0.6,0.6,0.6,0.6]         # xmin
x_max   = [1.0,1.0,1.0,1.0]         # xmax
y_min   = [0.0,0.0,0.0,0.0]         # ymin
y_max   = [2.0,5.0,2.0,200.]        # ymax

# Plot diagnostics
te_diag = ['TS','TSE'] #ECE
ne_diag = ['TS','TSE']
ti_diag = ['CES']
vt_diag = ['CES']

def _read_files(tdir):
	profiles = dict()
	file_n = tdir + '/chease_kinprof_extended'
	profiles['fit'] = _read_chease_kin_ext_file(file_n)
	file_n = tdir + '/psi_rho.dat'
	profiles['map'] = _read_psi_rho_files(file_n)
	profiles['mapf']= interp1d(profiles['map'][:,0],profiles['map'][:,1])
	#Te
	chs = ['te','ne','ti','vt']
	for ch in chs: profiles[ch]  = dict();
	profiles['te']['flag'] = te_diag
	profiles['ne']['flag'] = ne_diag
	profiles['ti']['flag'] = ti_diag
	profiles['vt']['flag'] = vt_diag

	for ch in chs:
		for flag in profiles[ch]['flag']:
			profiles[ch][flag] = dict()
			file_n = tdir + '/%s_avg_%s.dat'%(ch.upper(),flag.upper())
			profiles[ch][flag]['avg'] = _read_avg_files(file_n)
			file_n = tdir + '/%s_raw_%s.dat'%(ch.upper(),flag.upper())
			profiles[ch][flag]['raw'] = _read_raw_files(file_n)
	return profiles

def _read_avg_files(file_n):
	# Psi_Norm  Rho_Norm, VAL, STD
	dat = [];
	with open(file_n,'r') as f:
		line = f.readline()
		while True:
			line = f.readline();
			if not line: break;
			dat.append(line.split())
	dat = np.array(dat,dtype='float')
	for i in range(len(dat[:,0])): 
		if dat[i,0]>1.2: dat[i,0] = 1.2;
	return dat

def _read_raw_files(file_n):
	# Psi_Norm, VAL, STD, R, Rho_Norm
	with open(file_n,'r') as f:
		lend = int(f.readline())
		dat  = np.zeros((lend,5))
		for i in range(lend):
			line = f.readline();
			dat[i,:] = np.array(line.split(),dtype='float')
			if dat[i,0]>=1.2: dat[i,0] = 1.2;
	return dat

def _read_psi_rho_files(file_n):
	# Psi_Norm, Rho_Norm
	with open(file_n,'r') as f:
		lend = int(f.readline())
		dat  = np.zeros((lend,2))
		line = f.readline()
		for i in range(lend):
			line = f.readline();
			dat[i,:] = np.array(line.split(),dtype='float')
	return dat

def _read_chease_kin_ext_file(file_n):
	# Psi_Norm, Te, Ne, Ti, Vtor, Rout, Rin
	with open(file_n,'r') as f:
		lend = int(f.readline())
		dat  = np.zeros((lend,7))
		line = f.readline();
		for i in range(lend):
			line = f.readline();
			dat[i,:] = np.array(line.split(),dtype='float')
			if dat[i,0]>=1.2: dat[i,0] = 1.2;
	return dat

# Hardcoded
chs   = ['te','ne','ti','vt']
xlabel= ['$\\psi_N$ [a.u.]','$\\rho_N$ [a.u.]']
ylabel= ['$T_e$[keV]','$n_e$[$10^{19}m^{-3}$]','$T_i$[keV]','$V_{\\phi}$[km/s]']

lend = len(dirs)
profiles = dict();
for i in range(lend):
	profiles[i] = _read_files(dirs[i])

for i in range(4):
	if not draw_ch[i]: continue;
	plt.figure(i+1)
	val = plot_ra[i];
	if val.lower() == 'raw': val_ind = 1;
	else: val_ind = 2;
	tch = chs[i]
	llegend = [];

	for ii in range(lend):
		mapf = profiles[ii]['mapf']
		if not dodraw[ii]: continue
		llegend.append(labels[ii])
		pp = profiles[ii]
		for flag in pp[tch]['flag']:
			if xmap[i] == 0: plt.errorbar(pp[tch][flag][val][:,0],pp[tch][flag][val][:,val_ind],pp[tch][flag][val][:,val_ind+1],fmt='o',markersize='3',c=pcolor[ii],ecolor=pcolor[ii],capthick=2)
			else: plt.errorbar(mapf(pp[tch][flag][val][:,0]),pp[tch][flag][val][:,val_ind],pp[tch][flag][val][:,val_ind+1],fmt='o',markersize='3',c=pcolor[ii],ecolor=pcolor[ii],capthick=2)
		if xmap[i] == 0: plt.plot(pp['fit'][:,0],pp['fit'][:,chs.index(tch)+1],c=pcolor[ii])
		else: plt.plot(mapf(pp['fit'][:,0]),pp['fit'][:,chs.index(tch)+1],c=pcolor[ii])

	plt.xlabel(xlabel[xmap[i]])
	plt.ylabel(ylabel[i])
	plt.legend(llegend)
	plt.axis([x_min[i],x_max[i],y_min[i],y_max[i]])
plt.show()
