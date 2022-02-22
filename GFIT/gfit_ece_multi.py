# ECE loading and smooting 

import os,sys
from MDS import mds 
import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time as time_check
from scipy.interpolate import interp1d
from multi_mds import multi_mds 
#import pickle

def gefit_ece(shot, t_sample = 0.0002, lfs_option = 'y'):
	# Only 2nd harmonic is loaded! (recommanded by Dr. K.D. Lee) 
	start = time_check.time()
	msg = 'ECE post-processing is started...'
	print(msg) 
	g = mds('kstar',shot)
#	f2 = mds('EFIT01', shot)
#	t, rs = getmds(f2, 'rsurf', 0., 0., 0)
#	t_start = min(t[0] /1.e3,2.)
#	t_end   = t[-1]/1.e3 - 0.5

	mu0 = 4. * np.pi * 1e-7
	ttn = 56 * 16           # number of coils and turns
	f = mds('pcs_kstar',shot)
	tmp, tf = getmds(f, 'PCITFMSRD', 0.1, 0.2, 0.01)
	i_tf = np.mean(tf)      # TF current

	## Boundary
#	print('Boundary from EFIT01')
#	t, rs = getmds(f2, 'rsurf', t_start*1000., t_end*1000., 0)
#	t, a  = getmds(f2, 'aminor',t_start*1000., t_end*1000., 0)
#	t = t/1000
#	rs, tmp = mv_avg(rs, 5)
#	a, tmp  = mv_avg(a, 5)
#	bd_lfs_tmp = rs + a 
#	bd_hfs_tmp = rs - a 
#	bd_lfs = bd_lfs_tmp.max() 
#	bd_hfs = bd_hfs_tmp.min() 
#	print(bd_lfs)
#	print(bd_hfs)
	bd_hfs = 1.27
	bd_lfs = 2.31

	
	## Resonance position
	freq = []
	for i in range(0,76):
		ch_name = 'ECE'+str(i+1).zfill(2)+':freq'
		tt, yy = getmds(g, ch_name)
		freq.append(yy[0])
	freq = np.array(freq)
	B_res = 2. * np.pi * (1e+9 * freq) * 9.10e-31 / 1.602e-19
	r1 = mu0 * ttn * i_tf / (2. * np.pi * B_res / 1.)       # R_res
	r2 = mu0 * ttn * i_tf / (2. * np.pi * B_res / 2.)
	r3 = mu0 * ttn * i_tf / (2. * np.pi * B_res / 3.)

	## Valid channel (finding the resonance location inside the plasma) 
	ind_ch1 = getbetween(r1, [bd_hfs, bd_lfs])
	ind_ch2 = getbetween(r2, [bd_hfs, bd_lfs])
	ind_ch3 = getbetween(r3, [bd_hfs, bd_lfs])
	valid_ch1 = np.zeros(76)
	valid_ch2 = np.zeros(76)
	valid_ch3 = np.zeros(76)
	valid_ch1[ind_ch1] = 1
	valid_ch2[ind_ch2] = 1
	valid_ch3[ind_ch3] = 1
	valid_r1 = np.zeros(76)
	valid_r2 = np.zeros(76)
	valid_r3 = np.zeros(76)
	for i in range(0,76):
		if valid_ch1[i] == 1:
			valid_r1[i] = r1[i]
		if valid_ch2[i] == 1:
			valid_r2[i] = r2[i]
		if valid_ch3[i] == 1:
			valid_r3[i] = r3[i]

	## Overlap channel removal (finding 2 or 3 harmonics are located inside the plasma for one channel) 
	R1, INX1 = ch_overlap(valid_r1,valid_r2,valid_r3,bd_lfs,bd_hfs)
	R2, INX2 = ch_overlap(valid_r2,valid_r1,valid_r3,bd_lfs,bd_hfs)
	R3, INX3 = ch_overlap(valid_r3,valid_r1,valid_r2,bd_lfs,bd_hfs)
	if lfs_option == 'y':
		print('Only LFS channel will be considered') 
		ind_lfs1 = np.where(R1 < 1.8)
		ind_lfs2 = np.where(R2 < 1.8)
		ind_lfs3 = np.where(R3 < 1.8) 
		R1 = np.delete(R1, ind_lfs1)
		INX1 = np.delete(INX1, ind_lfs1)
		R2 = np.delete(R2, ind_lfs2)
		INX2 = np.delete(INX2, ind_lfs2)
		R3 = np.delete(R3, ind_lfs3)
		INX3 = np.delete(INX3, ind_lfs3)

	## ECE data load 
	try: os.mkdir('DATASAVE')
	except: pass
	try: os.mkdir('DATASAVE/%i'%shot)
	except: pass

	print(INX2)

	## Parallel MDS loading
	CH = []
	for i in range(0,len(INX2)):
		CH.append(INX2[i]+1) 
		
	ECE_pll = multi_mds(shot,'ECE',CH,'*','*',t_sample)
	ECE_pll.multi_load() 
	
	count1 = 1
	multi_count = 0 
	for i in INX2:
		try:
#			ch_name = 'ECE' + str(i).zfill(2)
			print('ECE Ch. %02i/%02i is avail'%(i+1-INX2[0],INX2[-1]-INX2[0]+1))
#			tt, yy = getmds(g, ch_name,t_start,t_end,t_sample)
			tt = ECE_pll.t[multi_count]
			yy = ECE_pll.y[multi_count]
			filename= 'DATASAVE/%i/ECE%i.npz'%(shot,count1)
			np.savez(filename,tt,R2,yy)
			comm = 'gzip ' + filename
			os.system(comm)
			count1 = count1 + 1
		except: print('ECE Ch. %i/%i is unavail'%(i+1-INX2[0],INX2[-1]-INX2[0])) 
		multi_count = multi_count + 1
	f = open('DATASAVE/%i/ECE_size'%shot,'w')
	f.write('%i'%(count1-1))
	f.close
		
	endtime = time_check.time() - start
	timedone = 'Working time is ... ' + str(endtime) +' [s]'
	print(timedone)
	return


def load_ece(shot,time,dt,dirs,fig):

	# Time type modfiication
	if str(type(time)) == "<class 'str'>":
		time = np.array(time.strip('\n').split(','),dtype='float')
	f = open('DATASAVE/%i/ECE_size'%shot,'r')
	nch = int(float(f.readline()))
	f.close()
	tlen = len(time)
	profv = np.zeros((tlen,nch))
	profe = np.zeros((tlen,nch))
	pltime = time[0] - max(0.8*dt,0.15)
	putime = time[0] + max(0.8*dt,0.15)

	for k in range(nch):
		filename = 'DATASAVE/%i/ECE%i.npz'%(shot,k+1)
		comm = 'gzip -d ' + filename +'.gz'
		os.system(comm)
		npzfile = np.load(filename, mmap_mode='r')
		tt = npzfile['arr_0']
		rr = npzfile['arr_1']
		te = npzfile['arr_2']
		comm = 'gzip ' + filename
		os.system(comm)
		if (k == 0 and not fig==None):
			pltime = max(pltime,tt[0])
			putime = min(putime,tt[-1])
			pind1 = np.where(tt>=pltime)
			pind2 = np.where(tt[pind1]<=putime)
			ptlen = len(pind2[0])
			ptt = np.zeros(ptlen)
			ptt[:] = tt[pind1][pind2]
			prr = np.zeros(nch)
			pprofv = np.zeros((ptlen,nch))
		
		if not fig == None:
			prr[k] = rr[k]
			pprofv[:,k] = te[pind1][pind2]
	
		for i in range(tlen):
			ltime = time[i] - dt*0.5
			utime = time[i] + dt*0.5
			#print(tt)
			ltime = max(ltime,tt[0])
			utime = min(utime,tt[-1])
		
			ind1 = np.where(tt>=ltime)
			ind2 = np.where(tt[ind1]<=utime)

			profv[i,k] = np.mean(te[ind1][ind2])
			profe[i,k] = np.std(te[ind1][ind2])

		if np.mean(profv[:,k]) < 100.:
			for i in range(tlen): 
				profv[i,k] = 10.
				profe[i,k] = 10.

	for i in range(tlen):
		f = open(dirs+'te_%06i_%ims_ECE_2nd.dat'%(shot,time[i]*1.e3),'w')
		f.write('R[m],   Z[m],   Te[eV], Error(Std.)[eV]\n')
		ind = np.argsort(rr)
		for k in ind:
			f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(rr[k],0.,profv[i,k],profe[i,k]))
		f.close()

	if not fig==None: draw_ece(fig,time[0],ltime,utime,ptt,prr,pprofv)

	return	

def draw_ece(fig,ttime,ltime,utime,tt,rr,te):

	t = np.array(tt)
	r = np.array(rr)
	te = np.array(te)	
	ax = fig.axes
	ax[0].cla()
	A=ax[0].contourf(tt,rr,np.transpose(te)/1.e3)
	ax[0].set_title('ECE-RT [keV]')
	ax[0].set_xlabel('time [s]')
	ax[0].set_ylabel('R [m]')
	ax[0].axvline(x=ttime,linestyle='--',c='magenta')
	ax[0].axvline(x=utime,linestyle='--',c='gold')
	ax[0].axvline(x=ltime,linestyle='--',c='gold')
	divider = make_axes_locatable(ax[0])
	cax = divider.append_axes('right', size='5%', pad=0.01)
	fig.colorbar(A,cax=cax)
	plt.tight_layout()
	plt.draw()

	return

def mv_avg(x,n,sigma_opt = 'n'):
	# for miving average with moving standard deviation values
	if n%2 == 0 : win = int(n/2);
	else: win = int((n-1)/2);
	y = np.zeros(len(x))
	s = np.zeros(len(x))

	for i in range(0,len(x)):
		if i < 0 + win or i > len(x) - win -1 : 
			y[i] = x[i]
		else:
			y_tmp = 0
			for ii in range(0,n):
				y_tmp = y_tmp + x[i-win+ii]
			y[i] = y_tmp / float(n)

			if sigma_opt == 'y':
				s_tmp = 0 
				for ii in range(0,n):
					s_tmp = s_tmp + (x[i-win+ii] - y[i])**2
				s[i] = np.sqrt(s_tmp / float(n))
	return y, s


def getnearest(x,value):
	# for target time index
	idx = (np.abs(x-value)).argmin()
	return idx


def getbetween(x,value):
	# for b/w time index
	idx = []
	tmp0 = x - value[0]
	tmp1 = x - value[1]
	tmp2 = tmp0 * tmp1
	for i in range(0,len(tmp2)):
		if tmp2[i] < 0:
			idx.append(i)
	return idx


def ch_overlap(x,y1,y2,bd_lfs,bd_hfs):
	# for removing harmonic-ovelapped channel 
	yy = []
	ind_yy = [] 
	for i in range(0,len(x)):
		if (x[i] > bd_lfs) or (x[i] < bd_hfs): pass;
		elif (y1[i] < bd_lfs) and (y1[i] > bd_hfs): pass;
		elif (y2[i] < bd_lfs) and (y2[i] > bd_hfs): pass;
		else:
			yy.append(x[i])
			ind_yy.append(i) 
	return np.array(yy), np.array(ind_yy) 


def sort_ece(r, y, dy):
	# sorting ECE channel 
	if len(r) == 0:
		R = r; Y = y; DY =dy; 
	else:
		ly =len(y[0,:])
		size = np.shape(r)
		tmp = np.zeros([size[0],1+ly+ly])
		tmp[:,0] = r
		tmp[:,1+0:1+ly] = y 
		tmp[:,ly+1+0:ly+1+ly] = dy
		tmp2 = sorted(tmp, key = lambda x : x[0])
		tmp2 = np.array(tmp2) 
		R = tmp2[:,0]
		Y = tmp2[:,1+0:1+ly]
		DY = tmp2[:,ly+1+0:ly+1+ly]

	return R, Y, DY


def getmds(g, tag, t1=0, t2=0, sampling_rate = 0):
	# for MDS data aquisition (resampling) 
	if t1 == t2:
		t, y = g.get('\\' + tag)
	elif t2 > t1:
		if sampling_rate > 0 :
			res_name = 'resample(\\' + tag + ',' +str(t1) + ',' + str(t2) + ',' + str(sampling_rate) + ')'
			t, y = g.get(res_name)
		else:
			t, y = g.get('\\' + tag + '[' + str(t1) + ':' + str(t2) + ']')
			if len(t) <= int(t2-t1)+2:
				tt , yy = g.get('\\' + tag)
				# time range cutting
				ttt = []
				yyy = []
				for i in range(0,len(tt)):
					if tt[i]>= t1 and tt[i] <= t2:
						ttt.append(tt[i])
						yyy.append(yy[i])
				t = np.array(ttt)
				y = np.array(yyy)
	else:
		print('start time(t1) should be smaller than end time(t2)!!')
		exit()
	return t, y

#shot = 22705
#gefit_ece(shot,lfs_option='n')
#load_ece([1.3,1.7,1.9,2.1])



