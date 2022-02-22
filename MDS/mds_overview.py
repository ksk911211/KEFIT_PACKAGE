import os,sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from progress import *

ists1  = False; ists2  = False; isces1 = False; isces2 = False;
istci  = False; isece = False; isref  = False;
dirs = 'MDS/DATASAVE/%i/'%int(float(sys.argv[1]))
dts = [0.,0.,0.]
try: times = np.array(sys.argv[2].split(','),dtype='int')
except: times = []
try: plot_type = sys.argv[3]
except: plot_type = '1d'
try: dts = np.array(sys.argv[4].split(','),dtype='float')
except: dts = [0.,0.,0.]

if os.path.isfile(dirs+'TE_size'): ists1 = True
if os.path.isfile(dirs+'NE_size'): ists2 = True
if os.path.isfile(dirs+'TI_size'): isces1= True
if os.path.isfile(dirs+'VT_size'): isces2= True
if os.path.isfile(dirs+'TCI_size'):istci = True
if os.path.isfile(dirs+'ECE_size'):isece = True 
if os.path.isfile(dirs+'REF_size'):isref = True

if plot_type == '1d':
	gs1 = gridspec.GridSpec(6,2,hspace=0.,left=0.05,right=0.99,bottom=0.1,top=0.95)
	gs2 = gridspec.GridSpec(7,2,hspace=0.,left=0.01,right=0.99,bottom=0.1,top=0.95)
	axes1 = [0,0,0,0,0,0]
	axes2 = [0,0,0,0,0,0,0]
	fig = plt.figure("MDS Time Overview",figsize=(10,6))
if plot_type == '2d':
	gs1 = gridspec.GridSpec(2,3,hspace=0.3,wspace=0.36,left=0.06,right=0.96,bottom=0.08,top=0.95)
	axes1 = [0,0,0,0,0,0]
	fig = plt.figure("MDS 2D Time",figsize=(13,7))

filename = dirs+'TE2.npz'
scale = [1.e-3,1.e-19,1.e-3,1.e0,1.e-3,1.e-19]
maxv  = [6.e0, 8.e0,  1.e1, 6.e2,4.e0]
def draw_plot1(isch,ch_name,pind=0):
	
	filename = dirs+ch_name+'.npz'
	if pind == 0: axes1[pind] = fig.add_subplot(gs1[pind,0])
	else: axes1[pind] = fig.add_subplot(gs1[pind,0],sharex=axes1[0])	
	if isch:
		comm = 'gzip -d '+filename+'.gz'
		os.system(comm)
		npzfile = np.load(filename,mmap_mode='r')
		comm = 'gzip ' + filename
		os.system(comm)
		time = npzfile['arr_0']
		if not ch_name.lower() == 'ece5': val  = npzfile['arr_1']
		else: val = npzfile['arr_2']
		val = np.nan_to_num(val)
		axes1[pind].plot(time,val*scale[pind])
		print('>>> %s exists'%ch_name.upper())
		return min(time),max(time)
	else:   return 1.e30, 0.

def draw_plot2d_tsces(isch,ch_name,pind):
	if not isch: return 1.e30, 0.
	filename = dirs+ch_name+'_size'
	f = open(filename,'r')
	nch = int(f.readline().split()[-1])
	f.close()
	print('>>> %s ->'%ch_name.upper())
	for k in range(nch):
		filename = dirs+ch_name+'%i.npz'%(k+1)
		comm = 'gzip -d '+filename+'.gz'
		os.system(comm)
		npzfile = np.load(filename,mmap_mode='r')
		comm = 'gzip ' + filename
		os.system(comm)
		time = npzfile['arr_0']
		if not ch_name.lower() == 'ece':
			val  = npzfile['arr_1']
			rad  = npzfile['arr_3']
		else:
			val  = npzfile['arr_2']
			rad  = npzfile['arr_1'][k]*1.e3
		if k == 0:
			times = np.array(time,dtype='float')
			tlen = len(times)
			prof = np.zeros((tlen,nch))
			radius= np.zeros(nch)
		val = np.nan_to_num(val)
		radius[k] = rad/1.e3
		prof[:,k] = val*scale[pind]
		update_progress(float(((k+1)/nch)))	

	if not ch_name.lower() == 'ece':	
		ind1 = np.where(times>2.0)
		if len(ind1[0]) > 1.: maxv2 = max(prof[ind1[0],0])*1.1
		else: maxv2 = max(prof[:,0])
		
		for k in range(nch):
			for i in range(tlen):
				prof[i,k] = min(max(prof[i,k],-10.),maxv2)
	else:
		ind1 = np.where(times>1.5)
		prof = prof[ind1[0],:]
		times= times[ind1[0]]
		for k in range(nch):
			for i in range(len(times)):
				if radius[k] > 2.2: prof[i,k] = min(prof[i,k],1.)

	A = axes1[pind].contourf(times,radius,np.transpose(prof),vmin=0)
	dv = make_axes_locatable(axes1[pind])
	cax = dv.append_axes('right',size='5%',pad=0.01)
	cbar= fig.colorbar(A,cax=cax)

	return min(time),max(time)

def draw_plot1r(axes,fig,gs1):
	axes[5] = fig.add_subplot(gs1[5,0],sharex=axes[0])
	filename = dirs + '/REF.npz'
	if not os.path.isfile(filename+'.gz'): 
		axes[5].text(0.5,0.5,'No data',color='r')
		axes[5].set_xlim(0,1)
		return 1.e30,0.
	comm='gzip -d '+ filename+'.gz'
	os.system(comm)
	f=open(filename,'rb')
	prof = pickle.load(f)
	f.close()
	comm='gzip '+ filename
	os.system(comm)
	axes[5].set_xlabel('time [s]')
	axes[5].set_ylim(0,1)
	if (not prof['isref']):	
		axes[5].text(0.5,0.5,'No data',color='r')
		return 1.e30,0.
	axes[5].set_xlim(0,1)
	nomds = np.sum(prof['prof']['is']) == 0
	axes[5].plot([prof['time'][0]-0.1,prof['time'][0]],[0,0],color='gray')
	for k in range(8):
		if prof['prof']['is'][k] == 1:
			axes[5].plot([prof['time'][k],prof['time'][k],prof['time'][k]+0.02,prof['time'][k]+0.02,],[0,0.5,0.5,0.],color='lime')
		else: axes[5].plot([prof['time'][k],prof['time'][k],prof['time'][k]+0.02,prof['time'][k]+0.02,],[0,0.5,0.5,0.],color='gray')
		if k<7.: axes[5].plot([prof['time'][k],prof['time'][k+1]],[0,0],color='gray')
		axes[5].plot([prof['time'][7],prof['time'][7]+0.1],[0,0],color='gray')
		axes[5].set_xlim(prof['time'][0]-0.1,prof['time'][7]+0.1)
	return 1.e30,0.

if plot_type=='1d':	
	min1,max1 = draw_plot1(ists1, 'TE2' ,0)
	min2,max2 = draw_plot1(ists2, 'NE2' ,1)
	min3,max3 = draw_plot1(isces1,'TI2' ,2)
	min4,max4 = draw_plot1(isces2,'VT2' ,3)
	min5,max5 = draw_plot1(isece, 'ECE5',4)
	min6,max6 = draw_plot1r(axes1,fig,gs1)#isref, 'RE2',5)

	mintt = min(min1,min2,min3,min4,min5,min6)
	maxtt = max(max1,max2,max3,max4,max5,max6)
	logics = [ists1,ists2,isces1,isces2,isece,isref]
	names  = ['TS-TE','TS-NE','CES-TI','CES-VT','ECE-TE','REF-NE']
	tlen = len(times)
	for k in range(6):
		axes1[k].axes.get_xaxis().set_visible(False)
		axes1[k].set_xlim(mintt,maxtt)
		if not logics[k]: 
			axes1[k].set_ylim(0,1.0)
			if k<5: axes1[k].text(.5*(mintt+maxtt),0.5,'No Data',color='r')
		for t in range(tlen):
			axes1[k].axvline(x=times[t]/1.e3,linestyle='--',color='gold')
		axes1[k].legend([names[k]],loc='upper right')
	
	axes1[0].set_title('Profile Ch.')
	axes1[-1].axes.get_xaxis().set_visible(True)
	axes1[-1].set_xlabel('times [s]')

	mintt = 1.e30; maxtt = 0.
	if istci: print('>>> TCI ch. exists!')
	if istci:
		names = ['int1','int2','tci1','tci2','tci3','tci4','tci5']
		isch  = [0,0,0,0,0,0,0]
		for k in range(7):
			filename = dirs+'ne_%s.npz'%names[k]
			if k==0: axes2[k] = fig.add_subplot(gs2[k,1])
			else: axes2[k] = fig.add_subplot(gs2[k,1],sharex=axes2[0])
			if os.path.isfile(filename+'.gz'):
				comm = 'gzip -d ' + filename+'.gz'
				os.system(comm)
				npzfile = np.load(filename,mmap_mode='r')
				comm = 'gzip ' + filename
				os.system(comm)
				time = npzfile['arr_0']
				val  = npzfile['arr_1']
				lent = len(time); lenv = len(val); lens = min(lent,lenv)
				time = time[:lens]; val = val[:lens]
				axes2[k].plot(time,val)		
				mintt = min(mintt,min(time))
				maxtt = max(maxtt,max(time))
				isch[k] = 1.

		if maxtt == 0.: mintt = 0.; maxtt = 1.;
		for k in range(7):
			axes2[k].axes.get_xaxis().set_visible(False)	
			axes2[k].set_xlim(mintt,maxtt)
			if not isch[k] == 1:
				axes2[k].set_ylim(0,1.0)
				axes2[k].text(.5*(mintt+maxtt),0.5,'No Data',color='r')
			for t in range(tlen):
				axes2[k].axvline(x=times[t]/1.e3,linestyle='--',color='gold')
			axes2[k].legend([names[k].upper()],loc='upper right')
		axes2[0].set_title('INTERFEROMETRY')
		axes2[-1].axes.get_xaxis().set_visible(True)
		axes2[-1].set_xlabel('times [s]')

	plt.show(block=True)

if plot_type == '2d':
	axes1[0] = fig.add_subplot(gs1[0,0])
	axes1[1] = fig.add_subplot(gs1[1,0],sharex=axes1[0])
	axes1[2] = fig.add_subplot(gs1[0,1],sharex=axes1[0])
	axes1[3] = fig.add_subplot(gs1[1,1],sharex=axes1[0])
	axes1[4] = fig.add_subplot(gs1[0,2],sharex=axes1[0])
	min1,max1 = draw_plot2d_tsces(ists1,'TE',0)
	min2,max2 = draw_plot2d_tsces(ists2,'NE',1)
	min3,max3 = draw_plot2d_tsces(isces1,'TI',2)
	min4,max4 = draw_plot2d_tsces(isces2,'VT',3)
	min5,max5 = draw_plot2d_tsces(isece,'ECE',4)	
	
	mint = min(min1,min2,min3,min4,min5)
	maxt = max(max1,max2,max3,max4,max5)
	if not ists1: axes1[0].text(0.5*(mint+maxt),0.5,'No data',color='red')
	if not ists2: axes1[1].text(0.5*(mint+maxt),0.5,'No data',color='red')
	if not isces1: axes1[2].text(0.5*(mint+maxt),0.5,'No data',color='red')
	if not isces2: axes1[3].text(0.5*(mint+maxt),0.5,'No data',color='red')
	if not isece: axes1[4].text(0.5*(mint+maxt),0.5,'No data',color='red')
	
	if ists1: axes1[0].set_ylim(1.8)
	if ists2: axes1[1].set_ylim(1.8)


	for k in range(5): 
		axes1[k].set_xlim(mint,maxt)
		axes1[k].set_xlabel('time [s]')	
		axes1[k].set_ylabel('R [m]')
	axes1[0].set_title('TS-$T_e$ [keV]')
	axes1[1].set_title('TS-$n_e [10^{19}m^{-3}]$')
	axes1[2].set_title('CES-$T_i$ [keV]')
	axes1[3].set_title('CES-$V_\phi$ [km/s]')
	axes1[4].set_title('ECE-$T_e$ [keV]')

	ttlen = len(times)
	ddts = [dts[0],dts[0],dts[1],dts[1],dts[2]]
	for k in range(5):
		for j in range(ttlen):
			axes1[k].axvline(x=times[j]/1.e3,linestyle='--',color='magenta')
			if ddts[k] >0.:
				axes1[k].axvline(x=(times[j]+0.5*ddts[k])/1.e3,linestyle='--',color='gold')
				axes1[k].axvline(x=(times[j]-0.5*ddts[k])/1.e3,linestyle='--',color='gold')

	plt.show(block=True)
