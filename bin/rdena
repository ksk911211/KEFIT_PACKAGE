#!/usr/local/anaconda3/bin/python3
import os,sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
from exec_dirs import rdena_db_dir

global shot_list
global data_prof

data_prof = dict()

def _find_shot_list():
	global shot_list
	shot_list = []
	for i in range(10000,35000):
		if os.path.isfile(rdena_db_dir+'/dena_%05i.sav'%i): shot_list.append(i)
	shot_list = np.array(shot_list,dtype='int');

def _load_file(shotn):
	global data_prof, shot_list
	if not shotn in shot_list: print('>>> No such shot number, %i'%shotn); return
	f = open(rdena_db_dir+'/dena_%05i.sav'%shotn,'rb')
	data_prof[shotn] = pickle.load(f); 
	f.close()

def _check_shot(shot2):
	global data_prof, shot_list
	if not shot2 in shot_list: print('>>> No such shot number, %i'%shot2); return
	if not shot2 in data_prof.keys(): print('>>> Load files,%i'%shot2);_load_file(shot2)
	
def _draw_plot_1d(slist,rloc):
	global data_prof, shot_list
	fig = plt.figure(); legends = []
	for shot in slist:
		_check_shot(shot)
		time = data_prof[shot]['times']
		val1d= np.zeros(len(time))
		for i in range(len(time)):
			tt       = time[i]
			xx1d     = data_prof[shot]['fit']['ne'][tt][0]
			rind     = np.argmin(abs(xx1d-rloc));
			val1d[i] = data_prof[shot]['fit']['ne'][tt][1][rind]
		plt.plot(time,val1d,linestyle='--',marker='x')
		legends.append('%s_loc_%s'%(shot,rloc))
	fig.legend(legends)
	plt.xlabel('time [ms]')
	plt.ylabel('$n_e [10^{19}/m^3]$')
	fig.tight_layout()
	plt.show(block=False)

def _draw_plot_1d_time(loc_list,shot2):
	global data_prof, shot_list
	shot = int(shot2)
	_check_shot(shot)
	fig = plt.figure(); legends = []
	time = data_prof[shot]['times']; val1d= np.zeros(len(time))
	for rloc in loc_list:
		tt = time[0]
		xx1d     = data_prof[shot]['fit']['ne'][tt][0]
		rind     = np.argmin(abs(xx1d-rloc));
		for i in range(len(time)):
			tt = time[i]
			val1d[i] = data_prof[shot]['fit']['ne'][tt][1][rind]
		plt.plot(time,val1d,linestyle='--',marker='x')
		legends.append('%s_loc_%s'%(shot,rloc))
	fig.legend(legends)
	plt.xlabel('time [ms]')
	plt.ylabel('$n_e [10^{19}/m^3]$')
	fig.tight_layout()
	plt.show(block=False)

def _draw_plot_2d(slist,tloc):
	global data_prof, shot_list
	fig = plt.figure(); legends = []
	for shot in slist:
		_check_shot(shot)
		time = data_prof[shot]['times']
		tind = np.argmin(abs(time-tloc));
		tt   = time[tind]
		xx1d     = data_prof[shot]['fit']['ne'][tt][0]
		val2d    = data_prof[shot]['fit']['ne'][tt][1]
		plt.plot(xx1d,val2d,linestyle='--')
		legends.append('%s_%s_ms'%(shot,tt))
	fig.legend(legends)
	plt.xlabel('$\\psi_N$')
	plt.ylabel('$n_e [10^{19}/m^3]$')
	fig.tight_layout()
	plt.show(block=False)

def _draw_plot_2d_time(tlist,shot2):
	global data_prof, shot_list
	shot = int(shot2)
	_check_shot(shot)
	fig = plt.figure(); legends = []
	for tt in tlist:
		time = data_prof[shot]['times']
		tind = np.argmin(abs(time-tt));
		tt   = time[tind]
		xx1d = data_prof[shot]['fit']['ne'][tt][0]
		val2d= data_prof[shot]['fit']['ne'][tt][1]
		plt.plot(xx1d,val2d,linestyle='--')
		legends.append('%s_%s_ms'%(shot,tt))
	fig.legend(legends)
	plt.xlabel('$\\psi_N$')
	plt.ylabel('$n_e [10^{19}/m^3]$')
	fig.tight_layout()
	plt.show(block=False)

def _draw_plot_tci_check(shot):
	global data_prof, shot_list
	_check_shot(shot)
	fig = plt.figure();
	time = data_prof[shot]['times']
	xx   = np.zeros(len(time))
	xxv  = np.zeros(len(time))
	xxs  = np.zeros(len(time))
	for i in range(5):
		ax = fig.add_subplot(3,2,i+1)
		for j in range(len(time)): 
			xx[j] = data_prof[shot]['fit']['nel'][time[j]][i+2]
			xxv[j]= data_prof[shot]['infile']['tciv'][time[j]][i+2]
			xxs[j]= data_prof[shot]['infile']['tcis'][time[j]][i+2]
		ax.scatter(time,xx,linestyle='--',marker='x',color='r')
		ax.errorbar(time,xxv,xxs)
		ax.legend(['tci%02i'%(i+1)])
		ax.set_xlabel('time [ms]')
		ax.set_ylabel('$n_e [10^{19}/m^3]$')
	fig.tight_layout()
	plt.show(block=False)

def _write_profiles(shot2,locs,rtype):
	global data_prof, shot_list
	shot = int(shot2)
	_check_shot(shot)
	time = data_prof[shot]['times']; val1d= np.zeros(len(time))
	xx1d = data_prof[shot]['fit']['ne'][time[0]][0]
	if rtype=='times':
		f=open('profile_%05i_times.dat'%(shot),'w')
		tlist = []; line = '%13s'%('psi_norm')
		for t2 in locs:
			tind = np.argmin(abs(t2-time))
			tlist.append(time[tind])
			line = line+'\t%13s'%(str(time[tind])+' ms')
		f.write(line+'\n');
		for ii in range(len(xx1d)):
			line = '%13.7e'%xx1d[ii]
			for tt in tlist: line = line + '\t%13.7e'%data_prof[shot]['fit']['ne'][tt][1][ii]
			f.write(line+'\n');
		f.close()
	elif rtype=='locs':
		f=open('profile_%05i_locs.dat'%(shot),'w')
		llist = []; line = '%13s'%('time [ms]')
		for l2 in locs:
			lind = np.argmin(abs(l2-xx1d))
			llist.append(lind)
			line = line+'\t%13s'%('psi_%4.3f'%l2)
		f.write(line+'\n');
		for tt in time:
			line = '%13s'%str(tt)
			for ll in llist: line = line + '\t%13.7e'%data_prof[shot]['fit']['ne'][tt][1][ll]
			f.write(line+'\n');
		f.close()
			
def _print_shot_list():
	print('------------------------')
	count = 0; line = ''
	for shot in shot_list: 
		line = line + '%05i   '%shot
		count += 1
		if count == 15: 
			print(line); count = 0; line = ''
	if count > 0: print(line)
	print('------------------------')

def _get_interaction(strings,dtype):
	out = None
	while out == None:
		_find_shot_list()
		ans = input('>>> '+strings+', [e for exit] \n>>> ')
		if ans.find('e')>-1: return 'e'
		ans = ans.lower()
		if dtype=='logic':
			if ans=='y': out = True
			elif ans=='n': out = False
		elif dtype=='int':
			try: out = int(ans)
			except: pass
		elif dtype=='float':
			try: out = float(ans)
			except: pass
	return out

def _get_slist(msg,dtype='int'):
	global shot_list
	out = []
	while len(out)==0:	
		_find_shot_list()
		ans = input('>>> %s with spaces or , [e for exit]\n>>> '%msg);
		if ans.find('e')>-1: 
			if dtype=='logic': return False
			else: return 'e'
		try:
			if ans.find(',') > -1: ans = ans.split(',')
			elif (ans.find('-') > -1 and msg=='shot'): 
				ans = ans.split('-')
				ans = np.linspace(int(ans[0]),int(ans[1]),(int(ans[1])-int(ans[0]))+1)
			else:   ans = ans.split()
			if dtype=='int': out = np.array(ans,dtype='int')
			else: out = np.array(ans,dtype='float')
		except: pass
		if (len(out)>0 and msg=='shot'):
			for ss in out:
				if not ss in shot_list: out = []
	
	return out

run_type = -1
while not run_type == 8:
	_find_shot_list()
	run_type = _get_interaction('1) Loc comparison, 2) 2D comparison, 3) Loc times, 4) 2D times, 5) Check, 6) Shot list, 7) Get data','int')
	if run_type == 'e': exit()
	if run_type == 1:
		stay = True
		slist= _get_slist('shot')
		if 'e'==slist: continue
		while stay:
			rloc = _get_interaction('location [0-1]?','float')
			if 'e'== rloc: stay=False; continue
			_draw_plot_1d(slist,rloc)
			stay = _get_interaction('Stay[y/n]?','logic')
	elif run_type==2:
		stay = True
		slist= _get_slist('shot')
		if 'e'==slist: continue
		while stay:
			tloc = _get_interaction('time [ms]?','float')
			if 'e'== tloc: stay=False; continue
			_draw_plot_2d(slist,tloc)
			stay = _get_interaction('Stay[y/n]?','logic')
	elif run_type==3:
		stay = True
		shot =  _get_interaction('shot','int')
		if 'e'== shot: continue
		while stay:
			loc_list = _get_slist('location [0-1] ?','float')		
			if 'e'==loc_list: stay=False; continue
			_draw_plot_1d_time(loc_list,shot)
			stay = _get_interaction('Stay[y/n]?','logic')
	elif run_type==4:
		stay = True
		shot = _get_interaction('shot','int')
		if 'e'== shot: continue
		while stay:
			tlist = _get_slist('time [ms] ?','int')
			if 'e'==tlist: stay=False; continue
			_draw_plot_2d_time(tlist,shot)
			stay = _get_interaction('Stay[y/n]?','logic')
	elif run_type==5:
		stay = True
		shot = _get_interaction('shot','int')
		if 'e'== shot: stay=False; continue
		_draw_plot_tci_check(shot)
	elif run_type==6:
		_print_shot_list()

	elif run_type==7:
		stay = True
		shot = _get_interaction('shot','int')
		if 'e'== shot: stay=False; continue
		rtype= _get_interaction('1)times, 2)locations','int')
		if 'e'== shot: stay=False; continue
		if rtype==1:
			tlist = _get_slist('time [ms]?','int')
			if 'e'== tlist: stay=False; continue
			_write_profiles(shot,tlist,'times')
		elif rtype==2:
			rlist =  _get_slist('location [0-1] ?','float')
			if 'e'== rlist: stay=False; continue
			_write_profiles(shot,rlist,'locs')
