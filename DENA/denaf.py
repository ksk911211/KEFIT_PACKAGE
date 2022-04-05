#!/usr/local/anaconda3/bin/python3
import os,sys, copy, time, pickle
import matplotlib.pyplot as plt
from matplotlib import gridspec

from get_efit import *
from progress import *
import eqdsk
import fittool
import rtefit

import numpy as np
from scipy.interpolate import interp1d
from exec_dirs import dummy_dir
from shutil import move, copyfile, copytree, rmtree
from MDS import mds
import pickle

currdir = os.getcwd()

class kstar_density_tool:

	def _den_tool(self):

		self._make_directories()
		self._declare_variables()
		self._initialise_variables()
		self._read_input()
		self._load_diag()
		self._make_time_list()		
		self._fit()
		return	

	def _make_directories(self):
		for dirs in ['profiles','gfiles','output']:
			if not os.path.isdir(dirs.upper()): os.mkdir(dirs.upper())
		return

	def _make_time_list(self):
		self.profiles['times'] = np.array([],dtype='int')
		tnow = self.note_in['STIME'];
		while tnow<=self.note_in['ETIME']:
			self.profiles['times'] = np.append(self.profiles['times'],tnow)
			tnow += self.note_in['DTIME']		

		self.profiles['gtime'] = np.copy(self.profiles['times'])
		if self.note_in['GFLAG']>= 1: elist = self.efit_list['times'][self.note_in['GFLAG']];
		else: elist = self.out_rt1;
		for i in range(len(self.profiles['times'])):
			tt = self.profiles['times'][i]
			tind = np.argmin(abs(tt-elist))
			self.profiles['gtime'][i] = elist[tind]
			if (i>0 and self.note_in['FIXEQ'] == 1):
				self.profiles['gtime'][i] = self.profiles['gtime'][0]

	def _fit(self):

		tlen = len(self.profiles['times'])
		if tlen ==0: return
		self._set_fitopt()
		self._update_fit_param()

		self._get_gfiles()
		self._make_ts_raw_files()
		self._make_tci_raw_files()

		for k in range(tlen): 
			self._dofit(k)
			update_progress(float((k+1)/tlen))
		self.profiles['didfit'] = True
		if not self.brutal:
			self._draw_tci_fit()
		if self.note_in['SAVEF'] == '':
			self._dump_profiles('denaf_%05i.sav'%self.note_in['SHOTN'])
		else:
			self._dump_profiles('%s'%self.note_in['SAVEF'].split('\n')[0])

		return

	def _dump_profiles(self,filename):
		f = open(filename,'wb')
		pickle.dump(self.profiles,f)
		f.close()
		return

	def _load_diag(self):

		tmin = 0.; tmax = 500000;
		print('>>> Load TS CORE...')
		g = mds('kstar',self.note_in['SHOTN'])
		#Load TS_CORE
		radius = np.array([])
		self.ts['core']['rr'] = radius
		self.ts['edge']['rr'] = radius
		for ch in range(30):
			node_name1 = '\\TS_CORE%i.CORE%i_NE'%(ch+1,ch+1)
			node_name2 = '\\TS_CORE%i.CORE%i_NERRH'%(ch+1,ch+1)
			node_name3 = '\\TS_CORE%i.CORE%i_POS'%(ch+1,ch+1)
			self.ts['core'][ch] = g.get(node_name1)
			self.ts['core_err'][ch] = g.get(node_name2)
			if len(self.ts['core'][ch][0])==0: break			
			if len(self.ts['core_err'][ch][0])==0:self.ts['core_err'][ch] = [self.ts['core'][ch][0],self.ts['core'][ch][0]/10.]
			self.ts['core']['nch'] = ch + 1
			rr = g.get(node_name3)
			radius = np.append(radius,rr[1]/1000.)
		if self.ts['core']['nch'] > 0: 
			self.ts['core']['rr'] = np.copy(radius)
			tmin = max(tmin,self.ts['core'][0][0][0])
			tmax = min(tmax,self.ts['core'][0][0][-1])

		print('>>> Load TS EDGE...')
		radius = np.array([])
		for ch in range(30):
			node_name1 = '\\TS_EDGE%i.EDGE%i_NE'%(ch+1,ch+1)
			node_name2 = '\\TS_EDGE%i.EDGE%i_NERRH'%(ch+1,ch+1)
			node_name3 = '\\TS_EDGE%i.EDGE%i_POS'%(ch+1,ch+1)
			self.ts['edge'][ch] = g.get(node_name1)
			self.ts['edge_err'][ch] = g.get(node_name2)
			if len(self.ts['edge'][ch][0])==0: break			
			if len(self.ts['edge_err'][ch][0])==0:self.ts['edge_err'][ch] = [self.ts['edge'][ch][0],self.ts['edge'][ch][0]/10.]
			self.ts['edge']['nch'] = ch + 1
			rr = g.get(node_name3)
			radius = np.append(radius,rr[1]/1000.)
		if self.ts['edge']['nch'] > 0: 
			self.ts['edge']['rr'] = np.copy(radius)
			tmin = max(tmin,self.ts['edge'][0][0][0])
			tmax = min(tmax,self.ts['edge'][0][0][-1])
			self.ts['edge']['rr'][0] = self.ts['edge']['rr'][0] + 0.002

		if self.ts['edge']['nch'] == 0: print('>>> TS data unavail.');
		else: print('>>> TS CORE/EDGE CH. %i/%i [#] '%(self.ts['core']['nch'],self.ts['edge']['nch']))
		print('>>> Load INT...')
		line = '>>> INT CH '
		length = [1.9,2.75]
		for ch in range(self.int['nch']):
			node_name1 = '\\NE_INT%02i'%(ch+1)
			self.int[ch+1]['val'] = g.get(node_name1)
			if len(self.int[ch+1]['val'][0]) > 0: 
				self.int[ch+1]['is'] = True
				dt = (self.int[ch+1]['val'][0][10]-self.int[ch+1]['val'][0][9])*1000
				dscale = 2 ; #int(2./dt);
				self.int[ch+1]['val'] = [self.int[ch+1]['val'][0][0::dscale],self.int[ch+1]['val'][1][0::dscale]/length[ch]]
				line = line + '%i '%(ch+1)
				tmin = max(tmin,self.int[ch+1]['val'][0][0])
				tmax = min(tmax,self.int[ch+1]['val'][0][-1])

		if not line=='>>> INT CH ': print(line+ 'avail.')
		else: print(line+ 'unavail.')
		print('>>> Load TCI...')
		line = '>>> TCI CH '
		for ch in range(self.tci['nch']):
			node_name1 = '\\NE_TCI%02i'%(ch+1)
			self.tci[ch+1]['val'] = g.get(node_name1)
			if len(self.tci[ch+1]['val'][0]) > 0: 
				self.tci[ch+1]['is'] = True
				dt = (self.tci[ch+1]['val'][0][10]-self.tci[ch+1]['val'][0][9])*1000
				dscale = 2; #int(2./dt);
				self.tci[ch+1]['val'] = [self.tci[ch+1]['val'][0][0::dscale],self.tci[ch+1]['val'][1][0::dscale]]
				line = line + '%i '%(ch+1)
				tmin = max(tmin,self.tci[ch+1]['val'][0][0])
				tmax = min(tmax,self.tci[ch+1]['val'][0][-1])

		if not line=='>>> TCI CH.': print(line+ 'avail.')
		else: print(line+ 'unavail.')
		print('>>> Load Ip/Da...')
		node_name1 = '\\pcrc03'
		node_name2 = '\\tor_ha10'
		self.mds['ip'] = g.get(node_name1)
		self.mds['da'] = g.get(node_name2)

		self.ts_core_ch = []; self.ts_edge_ch = [];
		for ch in range(5): self.ts_core_ch.append('%02i'%(self.ts['core']['nch']-ch))
		for ch in range(5): self.ts_edge_ch.append('%02i'%(ch+1))		

		lena = len(self.mds['da'][0]); lenv = len(self.mds['da'][1]); lend = lena- lenv
		if not lend == 0:
			if lend > 0: 
				for i in range(lend): self.mds['da'][1].append(0)
			else:
				for i in range(lend): self.mds['da'][0].append(0)

		self.tmin = tmin*1000.; self.tmax = tmax*1000.;
		print('>>> Load EFIT...')

		efitok = False; self.out_rt1=None
		if self.note_in['GFLAG'] > 0:
			self.out_rt1 = None
			self.efit_list = get_efit_list2(self.note_in['SHOTN'])
			for i in range(1,5): 
				if (self.efit_list['isefit'][i] and not efitok): efitok = True

		else:
			currdir = os.getcwd()
			if not os.path.isdir('GFILES/RT1'): os.mkdir('GFILES/RT1')
			os.chdir('GFILES/RT1')
			try: rtefit.run(self.note_in['SHOTN'],0.,10000.); print('>>> RTEFIT avail.'); efitok=True
			except: print('>>> RTEFIT unavail.');
			self.out_rt1 = list()
			status, output = subprocess.getstatusoutput('ls')
			output = output.split()
			for g in output:
				if g.find('%s_'%self.note_in['SHOTN']) > -1:
					tt = g.split('.')[0].split('_')[-1]
					self.out_rt1.append(tt)
			os.chdir(currdir)
		if not efitok: exit()

		return
		
	def _make_ts_raw_files(self):

		scale = 1.e-18;
		print('>>> Generate TS inputs...')
		if self.note_in['docal'] == 1: self.profiles['scale']['docal'] = True
		else: self.profiles['scale']['docal'] = False

		if self.note_in['FCORE'] == 1: self.profiles['scale']['fcore'] = True
		else: self.profiles['scale']['fcore'] = False		

		tot_nch = self.ts['core']['nch']+self.ts['edge']['nch']
		profv = dict(); profe = dict(); profr = np.zeros(tot_nch)
		profr[0:self.ts['core']['nch']] = self.ts['core']['rr']; profr[self.ts['core']['nch']:tot_nch] = self.ts['edge']['rr'];
		for k in range(len(self.profiles['times'])):
			time1 = self.profiles['times'][k]
			if not self.ts['core']['nch'] > 0:
				self.profiles['infile']['ts'][time1] = dummy_dir
				continue
			
			ltime = time1 - float(self.note_in['ATS']) * 0.5
			utime = time1 + float(self.note_in['ATS']) * 0.5
			core_ind1 = np.where((self.ts['core'][0][0]*1000 >= ltime))
			core_ind2 = np.where((self.ts['core'][0][0][core_ind1]*1000 <= utime))
			edge_ind1 = np.where((self.ts['edge'][0][0]*1000 >= ltime))
			edge_ind2 = np.where((self.ts['edge'][0][0][edge_ind1]*1000 <= utime))
			
			for i in range(tot_nch):
				if i < self.ts['core']['nch']: label = 'core'; j = i; ind1 = core_ind1; ind2 = core_ind2
				else: label = 'edge'; j = i - self.ts['core']['nch']; ind2 = edge_ind1; ind2 = edge_ind2

				profv[i] = np.nan_to_num(self.ts[label][j][1][ind1][ind2]) * scale
				profe[i] = np.nan_to_num(self.ts[label+'_err'][j][1][ind1][ind2]) * scale

			core_ind = int(self.note_in['CALCHC'])-1; edge_ind = int(self.note_in['CALCHE'])-1; ratio = float(self.note_in['CALCH'])
			factor = np.mean(profv[core_ind]) / np.mean(profv[edge_ind+self.ts['core']['nch']])/ratio
			if not self.note_in['docal']: factor = 1.

			self.profiles['scale'][time1] = factor

			if self.note_in['FCORE'] == 1:
				for i in range(self.ts['core']['nch'],tot_nch): profv[i] = profv[i]*factor; profe[i] = profe[i]*factor;
			else:
				for i in range(self.ts['core']['nch']): profv[i] = profv[i]/factor; profe[i] = profe[i]/factor;

			self.profiles['infile']['ts'][time1] = 'PROFILES/TS_NE_%ims.dat'%time1
			f = open(self.profiles['infile']['ts'][time1],'w')
			f.write('R[m]\tZ[m]\tVAL\tERR\n')
			for i in range(len(profv[0])):
				for j in range(tot_nch):
					f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(profr[j],0.,profv[j][i],profe[j][i]))
			f.close()

		return

	def _make_tci_raw_files(self):

		L_weight=[4,4,7.23,5.51,4.60,3.81,2.48]
		L_weight=[1,1,1,1,1,1,1]
		print('>>> Generate TCI/INT inputs...')

		for k in range(len(self.profiles['times'])):

			time1 = self.profiles['times'][k]
			iltime = time1 - float(self.note_in['AINT']) * 0.5
			iutime = time1 + float(self.note_in['AINT']) * 0.5			
			tltime = time1 - float(self.note_in['ATCI']) * 0.5
			tutime = time1 + float(self.note_in['ATCI']) * 0.5

			tot_nch = self.int['nch']+self.tci['nch']
			self.profiles['infile']['tciv'][time1] = np.zeros(tot_nch)
			self.profiles['infile']['tcis'][time1] = np.zeros(tot_nch)
			
			for i in range(tot_nch):
				if i<self.int['nch']:
					j = i+1; ltime = iltime; utime = iutime
					if self.int[j]['is']: dat = self.int[j]['val']
					else: dat = None
				else:
					j = i + 1 - self.int['nch']; ltime = tltime; utime = tutime;
					if self.tci[j]['is']: dat = self.tci[j]['val']
					else: dat = None

				if not dat == None:
					lenx = len(dat[0]); lenv = len(dat[1]); lenm = min(lenx,lenv);
					dat[0] = dat[0][0:lenm]; dat[1] = dat[1][0:lenm]
					ind1 = np.where(dat[0]*1000>ltime)
					ind2 = np.where(dat[0][ind1]*1000<utime)

					tciavg = np.sum(dat[1][ind1][ind2])/len(dat[1][ind1][ind2])
					ind3 = np.where(abs(dat[1][ind1][ind2]-tciavg) < 0.5)
					if len(dat[1][ind1][ind2][ind3])==0: continue

					tciavg = np.sum(dat[1][ind1][ind2][ind3])/len(dat[1][ind1][ind2][ind3])
				#	tcisig = (np.max(dat[1][ind1][ind2][ind3]) - np.min(dat[1][ind1][ind2][ind3])) / 3.
					tcisig = np.std(dat[1][ind1][ind2][ind3])* 1.5 / L_weight[i] #* 4.

					if tciavg<=0.: tciavg = 0.; tcisig = 0.
					self.profiles['infile']['tciv'][time1][i] = tciavg
					self.profiles['infile']['tcis'][time1][i] = tcisig


			self.profiles['infile']['tcif'][time1] = 'PROFILES/TCI_%ims.dat'%time1
			f = open(self.profiles['infile']['tcif'][time1],'w')
			for k in range(tot_nch):
				if self.profiles['infile']['tciv'][time1][k] > 0: isdata = 1
				else: isdata = 0
				val = self.profiles['infile']['tciv'][time1][k]; sig = self.profiles['infile']['tcis'][time1][k];
				f.write('%i %4.2f %4.3f \n'%(isdata,max(0.,val),max(0,sig)))
			f.close()
		return

	def _get_gfiles(self):

		print('>>> Download GFILES')
		GFLAG = ['RT1','EFIT01','EFIT02','EFIT03','EFIT04','EFIT05']
		tlen = len(self.profiles['times'])
		for flag in GFLAG:
			if not os.path.isdir('GFILES/%s'%flag): os.mkdir('GFILES/%s'%flag);

		efit_dir = 'GFILES/%s'%GFLAG[self.note_in['GFLAG']]
		shot = self.note_in['SHOTN']

		if self.note_in['GFLAG'] >= 1: self.profiles['GFLAG'] = 'EFIT%02i'%self.note_in['GFLAG']
		else: self.profiles['GFLAG'] = 'EFITRT1'

		if self.note_in['GFLAG'] > 0:	
			for years in shotk.keys():
				if ((shot-shotk[years]['shot'][0])*(shot-shotk[years]['shot'][-1])) < 0:
					tyear = years; break;
			efitdir = self.efit_list['dirs'][self.note_in['GFLAG']]
			fileline = ''; prevt = -1;
			for i in range(tlen):
				ttime = self.profiles['times'][i]
				gtime = self.profiles['gtime'][i]
				self.profiles['infile']['gfile'][ttime] = '%s/g%06i.%06i'%(efit_dir,shot,gtime)
				if not os.path.isfile(self.profiles['infile']['gfile'][ttime]):
					if not ttime==prevt:
						fileline = fileline + ' %sg%06i.%06i'%(efitdir,shot,gtime)
						prevt = ttime;
	
			pwd = os.getcwd(); os.chdir(efit_dir); 
			if not fileline == '':
			        comm = 'cp %s .'%(fileline)
			        os.system(comm)
			os.chdir(pwd)
		else: 
			for i in range(tlen):
				ttime = self.profiles['times'][i]
				gtime = self.profiles['gtime'][i]
				self.profiles['infile']['gfile'][ttime] = '%s/kstar_EFITRT1_%i_%06i.geqdsk'%(efit_dir,shot,gtime)
		return

	def _set_fitopt(self):

		use_rho = False; excn = '';
		if self.note_in['RHOFIT'] == 1: use_rho = True
		self.fit.exp_lsq = True
		self.fit.fit_opt['func_type']['ne'] = 2
		for k in range(self.ts['core']['nch']):
			if self.note_in['ts_core%i'%(k+1)] == 0:
				if excn == '': excn = '%i'%k
				else: excn = excn + ',%i'%k
		for k in range(self.ts['edge']['nch']):
			if self.note_in['ts_edge%i'%(k+1)] == 0:
				if excn == '': excn = '%i'%(k+self.ts['core']['nch'])
				else: excn = excn + ',%i'%(k+self.ts['core']['nch'])

		print('>>> Exclude %s'%excn)

		self.fit.fit_opt['use_rho']['ne'] = use_rho
		self.fit.fit_opt['exclude']['ne']['ts'] = excn
		self.fit.fit_opt['weight']['ne']['ts'] = float(self.note_in['WTS'])
		if not self.ts['core']['nch'] > 0: self.fit.fit_opt['weight']['ne']['ts'] = 0.
		self.fit.fit_opt['weight']['ne']['refl'] = float(self.note_in['WREFL'])
		if not self.ref['is']: self.fit.fit_opt['weight']['ne']['refl'] = 0.

		self.fit.noprint = False

		self.fit.fit_opt['oli']['ne']['ts']['use'] = True
		self.fit.fit_opt['oli']['ne']['ts']['per'] = 0.9
		self.fit.fit_opt['psi_end']['ne']['ts'] = 1.01
		self.fit.fit_opt['avg']['ne']['ts'] = False

		self.fit.fit_opt['scale']['ne']['ts']['core'] = float(self.note_in['CM'])
		self.fit.fit_opt['scale']['ne']['ts']['edge'] = float(self.note_in['EM'])

		if self.note_in['docal'] == 1:
			if self.note_in['FCORE'] == 1: self.fit.fit_opt['scale']['ne']['ts']['edge'] = float(self.note_in['CM'])
			else: self.fit.fit_opt['scale']['ne']['ts']['core'] = float(self.note_in['EM'])

		self.profiles['scale']['core'] = self.fit.fit_opt['scale']['ne']['ts']['core']
		self.profiles['scale']['edge'] = self.fit.fit_opt['scale']['ne']['ts']['edge']

		return

	def _dofit(self,tind):

		time1 = self.profiles['times'][tind]
		self.fit.ne_prof['fit_old'] = [];
		self.fit.fit_opt['file']['gfile']     = None
		self.notci = True
		for k in range(len(self.fit.inter_list)):
			flag = self.fit.inter_list[k]
			self.fit.fit_opt[flag]['val'] =self.profiles['infile']['tciv'][time1][k]

			if k<2: insig = float(self.note_in['WINT%02i'%(k+1)]); usech = 0
			else: insig = float(self.note_in['WTCI%02i'%(k-1)]); usech = 1;
			if insig == 0.: insig = self.profiles['infile']['tcis'][time1][k]
			else: self.notci = False
			self.fit.fit_opt[flag]['sig'] = insig #*self.fit.fit_opt[flag]['val'];
			if usech == 0: self.fit.fit_opt[flag]['sig'] = 0.
			else: self.notci = False

		if self.notci: self.fit.fit_opt['weight']['ne']['ts'] = 1.

		self.fit.fit_opt['file']['ne']['kfile'] = None
		for j in self.fit.__dict__['ne_list']: self.fit.fit_opt['file']['ne'][j] = None		

		self.fit.fit_opt['file']['gfile'] = self.profiles['infile']['gfile'][time1]
		self.fit.fit_opt['file']['ne']['ts'] = self.profiles['infile']['ts'][time1]
		self.fit.post['didfit']['ne'] = False
		self.fit.first_run = True
		self.fit.fit_opt['mds']['time'] = time1
		self.fit.main_run()
		self.fit.post['time'] = time1
		nel = [self.fit.post['int01'],self.fit.post['int02'],self.fit.post['tci01'],self.fit.post['tci02'],self.fit.post['tci03'],self.fit.post['tci04'],self.fit.post['tci05']]
		self.profiles['fit']['nel'][time1] = np.copy(np.array(nel))
		fit_eq = self.fit.fit_eq; fit_ne = self.fit.ne_prof;
		if self.note_in['RHOFIT'] == 1: fit_ne['ts']['xxr'] = fit_eq['rho_to_psi'](fit_ne['ts']['xxr'])
		tsx = np.copy(fit_ne['ts']['xxr']); 
		for i in range(len(tsx)): tsx[i] = min(tsx[i],1.1)
		tsf = interp1d(fit_eq['psin2'],fit_ne['fit2p']); tsv = tsf(tsx);
		self.profiles['fit']['ne'][time1] = copy.deepcopy([fit_eq['psin2'],fit_ne['fit2p'],fit_ne['ts'],fit_eq['psi_to_rho'],fit_eq['psif'],tsv,self.fit.post['chi']['ne']])

		return

	def _declare_variables(self):
		self.note_in    		= dict()
		self.profiles 			= dict()
		self.ts 				= dict()
		self.tci 			    = dict()
		self.int  				= dict()
		self.ref 				= dict()
		self.mds                = dict()
		return

	def _initialise_variables(self):

		self.tci['nch'] = 5
		self.int['nch'] = 2

		self.profiles['shotn']		= 0
		self.profiles['didfit']     = False
		self.profiles['tmin']		= 0
		self.profiles['tmax']		= 0
		self.profiles['delt']		= 100
		self.profiles['times']		= np.array([])
		self.profiles['gtime']		= np.array([])
		self.profiles['infile']     = dict()
		self.profiles['scale']      = dict()
		self.profiles['infile']['ts']     = dict()
		self.profiles['infile']['tciv']   = dict()
		self.profiles['infile']['tcis']   = dict()
		self.profiles['infile']['tcif']   = dict()
		self.profiles['infile']['gfile']  = dict()
		self.profiles['fit']			  = dict()
		self.profiles['fit']['nel']       = dict()
		self.profiles['fit']['ne']		  = dict()
		self.profiles['GFLAG'] 		      = 'EFIT01'

		self.ts['core'] = dict()
		self.ts['edge'] = dict()
		self.ts['core_err'] = dict()
		self.ts['edge_err'] = dict()		

		self.ts['core']['nch'] = 0
		self.ts['edge']['nch'] = 0

		for i in range(1,self.tci['nch']+1): 
			self.tci[i]       = dict()
			self.tci[i]['is'] = False

		for i in range(1,self.int['nch']+1): 
			self.int[i]       = dict()
			self.int[i]['is'] = False			
		self.ref['is'] = False

		for i in range(1,30):
			self.note_in['ts_core%i'%i] = 0
			self.note_in['ts_edge%i'%i] = 0
			self.note_in['tci%i'%i] = 0
			self.note_in['int%i'%i] = 0

		self.note_in['SHOTN'] = 0
		self.note_in['STIME'] = 1000.
		self.note_in['ETIME'] = 5000.
		self.note_in['DTIME'] = 100.
		self.note_in['SAVEF'] = ''

		self.note_in['WTCI01']= 0.01
		self.note_in['WTCI02']= 0.03
		self.note_in['WTCI03']= 0.03
		self.note_in['WTCI04']= 0.03
		self.note_in['WTCI05']= 0.02

		self.note_in['ATCI']  = 100.
		self.note_in['HMODE'] = 1
		self.note_in['GFLAG'] = 0
		self.note_in['FCORE'] = 1
		self.note_in['RHOFIT']= 0
		self.note_in['FTRUNC']= 12
		self.note_in['DTRUNC']= 12
		self.note_in['ISINP'] = 0
		self.note_in['FIXEQ'] = 0
		self.note_in['EXTCOE']= 0

		#unused variables (for future?)
		self.note_in['ATS']   = 100.
		self.note_in['AINT']  = 100.
		self.note_in['AREFL'] = 100.
		self.note_in['CM']    = 1.0
		self.note_in['EM']    = 1.0
		self.note_in['WTS']   = 0
		self.note_in['WINT01']= 0
		self.note_in['WINT02']= 0
		self.note_in['WREFL'] = 0		

		self.note_in['CALCH'] = 1.0
		self.note_in['CALCHC']= 13
		self.note_in['CALCHE']= 2
		self.note_in['WIDTH'] = 0.05
		self.note_in['docal'] = 0

		self.note_in['param'] = dict()
		for i in range(1,9):
			self.note_in['param'][i] = dict()
			self.note_in['param'][i]['min'] = -np.inf
			self.note_in['param'][i]['max'] = +np.inf
			self.note_in['param'][i]['val'] = 0.
			self.note_in['param'][i]['vary']= True

		self.profiles['fitopt']		        = dict()
		self.profiles['fitopt']['weight']       = dict()
		self.profiles['fitopt']['weight']['ts'] = 1.
		self.profiles['fitopt']['weight']['ref']= 1.
		self.profiles['fitopt']['tci_sig']      = np.zeros(10)
		self.profiles['fitopt']['int_sig']      = np.zeros(10)

		self._input_variables()
		
		return

	def _update_fit_param(self):

		self.fit.fit_trunc  = self.note_in['FTRUNC']
		self.fit.line_trunc = self.note_in['DTRUNC']

		self.fit.fit_opt['sep_fix']['ne'] = True
		self.fit.fit_opt['sep_val']['ne'] = 0.1		

		self.fit.param['ne']['min'][0]  = 0.1
		self.fit.param['ne']['max'][0]  = 0.5
		self.fit.param['ne']['val'][0]  = 0.1
		self.fit.param['ne']['vary'][0] = False

		self.fit.param['ne']['min'][1]  = 0.1     #0.1
		self.fit.param['ne']['max'][1]  = 5.0     #3.0
		self.fit.param['ne']['val'][1]  = 1.5     #2.0
		self.fit.param['ne']['vary'][1] = True

		self.fit.param['ne']['min'][2]  = 0.02
		self.fit.param['ne']['max'][2]  = 0.2
		self.fit.param['ne']['val'][2]  = 0.05
		self.fit.param['ne']['vary'][2] = False

		self.fit.param['ne']['min'][3]  = 0.1
		self.fit.param['ne']['max'][3]  = 100.0
		self.fit.param['ne']['val'][3]  = 0.98 #1.0
		self.fit.param['ne']['vary'][3] = False #True

		self.fit.param['ne']['min'][4]  = 0.
		self.fit.param['ne']['max'][4]  = 2.
		self.fit.param['ne']['val'][4]  = 0.1
		self.fit.param['ne']['vary'][4] = True

		self.fit.param['ne']['min'][5]  = 0.
		self.fit.param['ne']['max'][5]  = 400.
		self.fit.param['ne']['val'][5]  = 3.0
		self.fit.param['ne']['vary'][5] = True

		self.fit.param['ne']['min'][6]  = 0.
		self.fit.param['ne']['max'][6]  = 0.5
		self.fit.param['ne']['val'][6]  = 0.1
		self.fit.param['ne']['vary'][6] = False

		self.fit.param['ne']['min'][7]  = 1.01
		self.fit.param['ne']['max'][7]  = 3.0
		self.fit.param['ne']['val'][7]  = 1.5
		self.fit.param['ne']['vary'][7] = True


		if not self.note_in['HMODE'] == 1:

			self.fit.param['ne']['min'][2]  = 0.02
			self.fit.param['ne']['max'][2]  = 0.2
			self.fit.param['ne']['val'][2]  = 0.15
			self.fit.param['ne']['vary'][2] = False			

			self.fit.param['ne']['min'][6]  = 0.
			self.fit.param['ne']['max'][6]  = 0.5
			self.fit.param['ne']['val'][6]  = 0.1
			self.fit.param['ne']['vary'][6] = True
	
			self.fit.param['ne']['min'][7]  = 1.01
			self.fit.param['ne']['max'][7]  = 3.0
			self.fit.param['ne']['val'][7]  = 1.5
			self.fit.param['ne']['vary'][7] = True

		if self.note_in['EXTCOE'] == 0: return

		for i in range(8):
			for j in ['min','max','val','vary']:
				self.fit.param['ne'][j][i] = self.note_in['param'][i+1][j]
		return

	def _input_variables(self):
		self.input_list = ['SHOTN','STIME','ETIME','DTIME','ATCI','HMODE','GFLAG','RHOFIT','FTRUNC','DTRUNC','FIXEQ','EXTCOE','SAVEF']
		for i in range(1,6): self.input_list.append('WTCI%02i'%i)
		return

	def _write_input(self):

		f = open('indenaf','w')
		f.write('!Fit parameters\n')
		for key in self.input_list:
			f.write('%s = %s \n'%(key,self.note_in[key]))
		f.write('\n!Fit coefficients\n')

		for i in range(8):
			for j in ['min','max','val','vary']:
				f.write('Param_%s_%i = %s \n'%(j,i+1,self.note_in['param'][i+1][j]))
		f.close()
		return

	def _read_input(self):

		if not os.path.isfile(self.infile):
			print('>> No infile, generate input...')
			self._write_input()
			exit()

		f = open(self.infile,'r')
		while True:
			line = f.readline()
			if not line: break
			line = line.replace(' ', '').split('!')
			if line[0] == '': continue
			try:
				key = line[0].split('=')[0].upper()
				val = line[0].split('=')[1]
			except: continue

			if not key.find('PARAM') > -1:
				if   type(self.note_in[key]) is str:   self.note_in[key] = val
				elif type(self.note_in[key]) is float: self.note_in[key] = float(val)
				elif type(self.note_in[key]) is int:   self.note_in[key] = int(val)
			else:
				key = key.split('_')
				ind = int(key[2])
				typ = key[1].lower()
				if typ=='vary':
					if val.lower().split('\n')[0] == 'true': self.note_in['param'][ind][typ] = True
					else: self.note_in['param'][ind][typ] = False
				else: self.note_in['param'][ind][typ] = float(val)

		f.close()
		
		self.profiles['shotn'] = self.note_in['SHOTN']
		self.note_in['ISINP'] = 1
		return

	def __init__(self,infile):
		self.brutal = True
		self.fit = fittool.fit_tool()
		self.infile = infile
		return

if __name__ == "__main__":

	import denaf as dcal

	try: os.mkdir('PROFILES')
	except: pass
	try: os.mkdir('MDS')
	except: pass

	now = time.gmtime(time.time())
	years = now.tm_year, now.tm_mon, now.tm_mday
	hours = now.tm_hour, now.tm_min, now.tm_sec
	line = '#-- %i/%02i/%02i -- %02ih %02im %02is --'%(years[0],years[1],years[2],hours[0],hours[1],hours[2])
	os.system("echo '%s' >> log.denaf"%line)
	print('>>> RUN DENAF ')
	try: infile = sys.argv[1]
	except: infile = 'indenaf'
	dcal = dcal.kstar_density_tool(infile)
	dcal._den_tool()
