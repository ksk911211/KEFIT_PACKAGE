#!/usr/local/miniconda3/bin/python3
import os,sys
import numpy as np
from MDS import mds
import time

import tkinter as tk
from tkinter import ttk

from shutil import copyfile
import copy
import eqdsk
import pickle

from scipy.interpolate import interp1d, interp2d

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.widgets import MultiCursor
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

from get_efit import *
from exec_dirs import gzip_dir, efit_rmp, efit_dir, version

from aeqdsk import _read_afile

class kstar_diagnostic_tool:

	def _mds_tool(self):
		self._declare_variables()
		self._initialise_variables()
		self._make_directories()

		self._load_diagnostics()
		if self.nogui==1: exit()
		self._load_opt()

		self.root.protocol('WM_DELETE_WINDOW')
		self._make_home_canvas()
		self._make_note_frame()

		for k in range(4): self.root.grid_columnconfigure(k,weight=1)

		self._make_input_frame()
		self._gui_preset()
						
		self.root.resizable(0,0)
		self.root.mainloop()		
		
		return

	def _make_directories(self):
		if not os.path.isdir('DATASAVE'): os.mkdir('DATASAVE')
		if not os.path.isdir('DATASAVE/%i'%self.shotn): os.mkdir('DATASAVE/%i'%self.shotn)
		return		

	def _declare_variables(self):

		self.figure    		   = dict()
		self.figure['name']    = dict()
		self.figure['size']    = dict()
		self.figure['canvas']  = dict()
		self.figure['widget']  = dict()
		self.figure['toolbar'] = dict()
		self.figure['legend']  = dict()
		self.figure['pegend']  = dict()
		self.figure['axes']    = dict()
		self.figure['gs']      = dict()

		self.home       	   = dict()
		self.note_in           = dict()
		self.opt               = dict()
	
		self.mds               = dict()
		self.ces               = dict()
		self.ts 			   = dict()
		self.int  		       = dict()
		self.tci 			   = dict()		
		self.rr                = dict()
		return

	def _initialise_variables(self):

		self.prev_xmap  = 1
		self.efit_no_g  = tk.IntVar()
		self.efit_no_k  = tk.IntVar()
		self.efit_no_g.set(1)
		self.efit_no_k.set(1)
		self.opt['g_no'] = 1
		self.opt['k_no'] = 1
		self.tci['nch'] = 5
		self.int['nch'] = 2	
		self.ces['nch'] = 0
		self.ts['core'] = dict()
		self.ts['edge'] = dict()
		self.ts['core']['nch'] = 0
		self.ts['edge']['nch'] = 0		

		self.tfig_list	= []
		self.pfig_list  = []
		self.check_list	= []

		for ch in ['ti','vt']:
			self.ces[ch] = dict()
			for flag in ['val','err']:
				self.ces[ch][flag] = dict()	
		for ch in ['te','ne']:
			self.ts[ch] = dict()
			for flag in ['core','edge','core_err','edge_err']:
				self.ts[ch][flag] = dict()

		for i in range(1,self.tci['nch']+1): 
			self.tci[i]       = dict()
			self.tci[i]['is'] = False

		for i in range(1,self.int['nch']+1): 
			self.int[i]       = dict()
			self.int[i]['is'] = False			

		self.opt['dacrit'] = 1.5
		self.opt['duty']   = 40
		self.opt['peaks']  = []

		self.opt['WTCI01'] = 0.02
		self.opt['WTCI02'] = 0.02
		self.opt['WTCI03'] = 0.02
		self.opt['WTCI04'] = 0.02
		self.opt['WTCI05'] = 0.02
		self.opt['WINT01'] = 0.02
		self.opt['WINT02'] = 0.02

		self.opt['TCES']  = self.time
		self.opt['TTS']   = self.time
		self.opt['TINT']  = self.time
		self.opt['TTCI']  = self.time

		self.opt['slist'] = dict()
		self.opt['slist']['ces'] = dict()
		self.opt['slist']['ts'] = dict()
		self.opt['slist']['ces']['ti'] = []
		self.opt['slist']['ces']['vt'] = []
		self.opt['slist']['ts']['te'] = []
		self.opt['slist']['ts']['ne'] = []

		self.opt['ACES']  = self.aces
		self.opt['ATS']   = self.ats
		self.opt['AINT']  = 50
		self.opt['ATCI']  = 50

		self.note_in['kflag'] = tk.IntVar()
		self.note_in['xmap']  = tk.IntVar()
		self.note_in['xmap'].set(1)
			
		self.opt['cpage'] = '';
		self.opt['c_ch']  = '';

		self.rr['ts'] = dict()
		self.rr['ts']['core'] = dict()
		self.rr['ts']['edge'] = dict()
		self.rr['ces'] = dict()


		self.opt['xmap']  = dict()	
		self.opt['xmap']['ces'] = dict()
		self.opt['xmap']['ts'] = dict()
		self.opt['xmap']['ces']['ti'] = 1
		self.opt['xmap']['ces']['vt'] = 1
		self.opt['xmap']['ts']['te'] = 1
		self.opt['xmap']['ts']['ne'] = 1

		self.opt['ready'] = dict()
		self.opt['ready']['ces'] = dict()
		self.opt['ready']['ts']  = dict()
		self.opt['ready']['ces']['ti'] = False
		self.opt['ready']['ces']['vt'] = False
		self.opt['ready']['ts']['te'] = False
		self.opt['ready']['ts']['ne'] = False

		self.opt['isdata']= dict()
		self.opt['isdata']['ces'] = False
		self.opt['isdata']['ts']  = False

		self.figure['size']['home'] = (11,7.4)
		self.figure['size']['plot'] = (3.8,0.9)

		self.opt['range'] = dict()
		for ch in ['ti','vt','te','ne']:
			self.opt['range'][ch] = dict()
			self.opt['range'][ch]['xmin'] = 1.7
			self.opt['range'][ch]['xmax'] = 2.3
			self.opt['range'][ch]['ymin'] = 0.0
			self.opt['range'][ch]['ymax'] = 5.0
		self.opt['range']['vt']['ymax']	= 300
		return		
		
	def _load_diagnostics(self):
		self.g = mds('kstar',self.shotn)
		self._get_year()
		self.savfile='./DATASAVE/'+str(self.shotn)+'/'
		self._load_ip_da()
		if not self.nogui: self._load_eq()
		self._load_ces()
		self._load_ts()
		self._load_interf()
		self.g.close()
		self._diagnostics_exceptions()
		if not self.nogui: self._get_xmap()
		return

	def _get_save_zipf(self,nfile,data=[],ndata=4):

		if len(data)==0:
			comm=gzip_dir+' -d '+nfile+'.gz'
			os.system(comm)
			npzfile  = np.load(nfile, mmap_mode='r')
			ans = [];
			for i in range(ndata): ans.append(npzfile['arr_%i'%i])
			comm=gzip_dir+' '+nfile
			os.system(comm)
			return ans
		else:
			if ndata ==4: np.savez(nfile,data[0],data[1],data[2],data[3])
			elif ndata ==2: np.savez(nfile,data[0],data[1])
			elif ndata ==3: np.savez(nfile,data[0],data[1],data[2])
			comm=gzip_dir+' '+nfile
			os.system(comm)
			return

	def _load_ces(self):

		print('>>> Load CES...')
		for flag in ['ti','vt']:
					
			if flag == 'ti': self.ces['rr'] = np.array([])	
			for ch in range(1,50):
				nfile      = self.savfile+'%s%i.npz'%(flag.upper(),ch)
				node_name1 = '\\CES_%s%02i'%(flag.upper(),ch)
				node_name2 = '\\CES_%s%02i:err_bar'%(flag.upper(),ch)
				node_name3 = '\\CES_RT%02i'%(ch)

				if not os.path.isfile(nfile+'.gz'):
					self.ces[flag]['val'][ch]    = self.g.get(node_name1)
					self.ces[flag]['err'][ch]    = self.g.get(node_name2)
					temp = self.g.get(node_name3);
					try: 
						if len(self.ces[flag]['val'][ch][0])==0: break
					except: 
						self.ces[flag]['val'][ch] = [[],[]]
						self.ces[flag]['err'][ch] = [[],[]]
						break
					rr = temp[1][0];			
					if len(self.ces[flag]['err'][ch][0])==0: self.ces[flag]['err'][ch] = [self.ces[flag]['val'][ch][0],self.ces[flag]['val'][ch][1]/10.]
					self._get_save_zipf(nfile,[self.ces[flag]['val'][ch][0],self.ces[flag]['val'][ch][1],self.ces[flag]['err'][ch][1],rr])
				else:
					data  = self._get_save_zipf(nfile)
					self.ces[flag]['val'][ch] = [data[0]]; self.ces[flag]['val'][ch].append(data[1])
					self.ces[flag]['err'][ch] = [data[0]]; self.ces[flag]['err'][ch].append(data[2])
					rr = data[3]
								
				if flag=='ti': self.ces['rr']=np.append(self.ces['rr'],rr/1.e3); self.ces['nch'] = ch

			if self.ces['nch'] == 0: print('>>> CES %s data unavail.'%(flag.upper()));
			else: print('>>> CES %s CH. %i [#] '%(flag.upper(),self.ces['nch']))

			nfile      = self.savfile+'%s_size'%(flag.upper())
			if not os.path.isfile(nfile):
				f = open(nfile,'w')
				f.write('%i %i %i %i\n'%(len(self.ces[flag]['val'][ch][0]),len(self.mds['da'][0]),self.ces['nch'],self.ces['nch']))
				f.close()	

		use_tgf = 'n'; tgf_dir = '';
		tgf_default = os.environ["TGF_PATH"]+'/twoCES@%i.txt'%self.shotn
		if self.ces['nch'] == 0: 
			if os.path.isfile(tgf_default):
				print('>>> No CES on MDSplus but found TGF file at TGF_PATH!')
				print('>>> TGF file -> %s'%tgf_default)
				tgf_dir = tgf_default
				use_tgf = 'y'
			else:
				print('>>> ----------------------------------------------------------------------------------------')
				print('>>> Do you have two gaussian fitting (TGF)? > Then, follow this procedure')
				print('>>> 1. Set env."TGF_PATH" where TGF file will be stored, ex)/home/%s/TGF'%(os.environ['USER']))
				print('>>> 2. Put TGF file to "TGF_PATH" with name two@xxxxx.txt, ex) two@30000 for shot No. 30000')
				print('>>> 3. Re-run MDS Load')
				print('>>> ----------------------------------------------------------------------------------------')

		if os.path.isfile(tgf_default) and self.ces['nch']>0:
			use_tgf=input('>>> Do you want to use TGF instead of MDS?[y/n]  ')
			if use_tgf.lower()=='y':
				tgf_dir = tgf_default
			else: use_tgf = 'n'

		if (use_tgf.lower()=='y' and os.path.isfile(tgf_dir)):
			self.ces['rr'] = np.array([])
			ces_dat = self._read_tgf(tgf_dir)
			self.ces['nch'] = ces_dat['radius'].shape[0]
			for ch in range(1,self.ces['nch']+1):
				nfile_t      = self.savfile+'%s%i.npz'%('TI',ch)
				nfile_v      = self.savfile+'%s%i.npz'%('VT',ch)
				self.ces['ti']['val'][ch]    = [ces_dat['times'],ces_dat['Ti'][ch-1]*1.e+3]
				self.ces['ti']['err'][ch]    = [ces_dat['times'],ces_dat['Ti'][ch-1]*1.e+2]
				self.ces['vt']['val'][ch]    = [ces_dat['times'],ces_dat['Vc'][ch-1]*1.e+0]
				self.ces['vt']['err'][ch]    = [ces_dat['times'],ces_dat['Vc'][ch-1]*1.e-1]

				rr = ces_dat['radius'][ch-1]*1.e3
				if os.path.isfile(nfile_t+'.gz'): os.remove(nfile_t+'.gz')
				if os.path.isfile(nfile_v+'.gz'): os.remove(nfile_v+'.gz')	
				self._get_save_zipf(nfile_t,[self.ces['ti']['val'][ch][0],self.ces['ti']['val'][ch][1],self.ces['ti']['err'][ch][1],rr])
				self._get_save_zipf(nfile_v,[self.ces['vt']['val'][ch][0],self.ces['vt']['val'][ch][1],self.ces['vt']['err'][ch][1],rr])
			
				self.ces['rr']=np.append(self.ces['rr'],rr/1.e3);
	
			nfile_t      = self.savfile+'%s_size'%('TI')
			nfile_v      = self.savfile+'%s_size'%('VT')
			for ch in range(1,self.ces['nch']+1):
				with open(nfile_t,'w') as f:
					f.write('%i %i %i %i\n'%(len(self.ces['ti']['val'][ch][0]),len(self.mds['da'][0]),self.ces['nch'],self.ces['nch']))
				with open(nfile_v,'w') as f:
					f.write('%i %i %i %i\n'%(len(self.ces['vt']['val'][ch][0]),len(self.mds['da'][0]),self.ces['nch'],self.ces['nch']))
			print('>>> CES TI/VT CH. %i [#] '%(self.ces['nch']))
		return

	def _read_tgf(self,nfile):

		with open(nfile,'r') as f:

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
		return dat

	def _load_ts(self):

		print('>>> Load TS CORE...')
		nmax = 30;
		sfile = self.savfile+'TE_size'
		if os.path.isfile(sfile):
			f = open(sfile,'r')
			line = f.readline(); f.close()
			nmax = int(line.split()[-1]) + 1
		for flag in ['te','ne']:
			if flag == 'te': factor = 1.e-3; self.ts['core']['rr'] = np.array([])
			else: factor = 1.e-19;
			for ch in range(1,nmax):
				nfile      = self.savfile+'%s%i.npz'%(flag.upper(),ch)
				node_name1 = '\\TS_CORE%i.CORE%i_%s'%(ch,ch,flag.upper())
				node_name2 = '\\TS_CORE%i.CORE%i_%sRRH'%(ch,ch,flag.upper())
				node_name3 = '\\TS_CORE%i.CORE%i_POS'%(ch,ch)

				if not os.path.isfile(nfile+'.gz'):
					self.ts[flag]['core'][ch]        = self.g.get(node_name1)
					self.ts[flag]['core_err'][ch]    = self.g.get(node_name2)
					temp = self.g.get(node_name3);
					if len(self.ts[flag]['core'][ch][0])==0: break			
					if len(self.ts[flag]['core_err'][ch][0])==0: self.ts[flag]['core_err'][ch] = [self.ts[flag]['core'][ch][0],self.ts[flag]['core'][ch][1]/10.]
					rr = temp[1][0];
					self._get_save_zipf(nfile,[self.ts[flag]['core'][ch][0],self.ts[flag]['core'][ch][1],self.ts[flag]['core_err'][ch][1],rr])
				else:
					data  = self._get_save_zipf(nfile)
					self.ts[flag]['core'][ch] = [data[0]]; self.ts[flag]['core'][ch].append(data[1])
					self.ts[flag]['core_err'][ch] = [data[0]]; self.ts[flag]['core_err'][ch].append(data[2])
					rr = data[3]
								
				if flag=='te': self.ts['core']['rr']=np.append(self.ts['core']['rr'],rr/1.e3); self.ts['core']['nch'] = ch	

		print('>>> Load TS EDGE...')
		for flag in ['te','ne']:			
			if flag == 'te': self.ts['edge']['rr'] = np.array([])	
			for ch in range(1,30):
				nfile      = self.savfile+'%s%i.npz'%(flag.upper(),self.ts['core']['nch']+ch)
				node_name1 = '\\TS_EDGE%i.EDGE%i_%s'%(ch,ch,flag.upper())
				node_name2 = '\\TS_EDGE%i.EDGE%i_%sRRH'%(ch,ch,flag.upper())
				node_name3 = '\\TS_EDGE%i.EDGE%i_POS'%(ch,ch)

				if not os.path.isfile(nfile+'.gz'):
					self.ts[flag]['edge'][ch]        = self.g.get(node_name1)
					self.ts[flag]['edge_err'][ch]    = self.g.get(node_name2)				
					temp = self.g.get(node_name3);
					if len(self.ts[flag]['edge'][ch][0])==0: break			
					if len(self.ts[flag]['edge_err'][ch][0])==0: self.ts[flag]['edge_err'][ch] = [self.ts[flag]['edge'][ch][0],self.ts[flag]['edge'][ch][1]/10.]
					rr = temp[1][0];
					self._get_save_zipf(nfile,[self.ts[flag]['edge'][ch][0],self.ts[flag]['edge'][ch][1],self.ts[flag]['edge_err'][ch][1],rr])
				else:
					data  = self._get_save_zipf(nfile)
					self.ts[flag]['edge'][ch] = [data[0]]; self.ts[flag]['edge'][ch].append(data[1])
					self.ts[flag]['edge_err'][ch] = [data[0]]; self.ts[flag]['edge_err'][ch].append(data[2])
					rr = data[3]

				if flag=='te': self.ts['edge']['rr']=np.append(self.ts['edge']['rr'],rr/1.e3); self.ts['edge']['nch'] = ch	

			if self.ts['edge']['nch'] == 0: print('>>> TS data unavail.');
			else: print('>>> TS %s CORE/EDGE CH. %i/%i [#] '%(flag.upper(),self.ts['core']['nch'],self.ts['edge']['nch']))		

			nfile      = self.savfile+'%s_size'%(flag.upper())
			if not os.path.isfile(nfile):
				f = open(nfile,'w')
				f.write('%i %i %i %i\n'%(len(self.ts[flag]['edge'][ch][0]),len(self.mds['da'][0]),self.ts['core']['nch']+self.ts['edge']['nch'],self.ts['core']['nch']))
				f.close()
		return

	def _load_interf(self):
		print('>>> Load INT...')
		line = '>>> INT CH '
		for ch in range(1,self.int['nch']+1):
			node_name1 = '\\NE_INT%02i'%(ch)
			nfile = self.savfile+'ne_int%i.npz'%ch
			self.int[ch]['val'] = self.g.get(node_name1)
			if not os.path.isfile(nfile+'.gz'):			
				if (len(self.int[ch]['val'][0]) > 0): 
					self.int[ch]['is'] = True
					dt = (self.int[ch]['val'][0][10]-self.int[ch]['val'][0][9])*1000
					dscale = int(2./dt);
					self.int[ch]['val'] = [self.int[ch]['val'][0][0::dscale],self.int[ch]['val'][1][0::dscale]]
					self._get_save_zipf(nfile,self.int[ch]['val'],2)
			else:
					data  = self._get_save_zipf(nfile,[],2)
					self.int[ch]['is'] = True
					self.int[ch]['val'] = [data[0]]; self.int[ch]['val'].append(data[1])
			if self.int[ch]['is']: line = line + '%i '%(ch)

		if not line=='>>> INT CH ': print(line+ 'avail.')
		else: print(line+ 'unavail.')
		print('>>> Load TCI...')
		line = '>>> TCI CH '
		for ch in range(1,self.tci['nch']+1):
			node_name1 = '\\NE_TCI%02i'%(ch)
			nfile = self.savfile+'ne_tci%i.npz'%ch
			self.tci[ch]['val'] = self.g.get(node_name1)
			if not os.path.isfile(nfile+'.gz'):					
				if (len(self.tci[ch]['val'][0]) > 0): 
					self.tci[ch]['is'] = True
					dt = (self.tci[ch]['val'][0][10]-self.tci[ch]['val'][0][9])*1000
					dscale = int(2./dt);
					self.tci[ch]['val'] = [self.tci[ch]['val'][0][0::dscale],self.tci[ch]['val'][1][0::dscale]]
					self._get_save_zipf(nfile,self.tci[ch]['val'],2)
			else:
					data  = self._get_save_zipf(nfile,[],2)
					self.tci[ch]['is'] = True
					self.tci[ch]['val'] = [data[0]]; self.tci[ch]['val'].append(data[1])
			if self.tci[ch]['is']: line = line + '%i '%(ch)

		if not line=='>>> TCI CH.': print(line+ 'avail.')
		else: print(line+ 'unavail.')

		nfile = self.savfile+'TCI_size'
		if not os.path.isfile(nfile):
			f = open(nfile,'w')
			for ch in range(1,self.int['nch']+1): 
				if self.int[ch]['is']: f.write('1 %i \n'%(len(self.int[ch]['val'][0])))
				else: f.write('0 0 \n')

			for ch in range(1,self.tci['nch']+1): 
				if self.tci[ch]['is']: f.write('1 %i \n'%(len(self.tci[ch]['val'][0])))
				else: f.write('0 0 \n')
			f.close()
		return

	def _load_ip_da(self):

		print('>>> Load Ip/Da...')
		node_name1 = '\\pcrc03'
		node_name2 = '\\tor_ha10'
		node_name3 = '\\wtot_dlm03'

		nfile1 = self.savfile+'IP.npz'
		nfile2 = self.savfile+'DA.npz'
		nfile3 = self.savfile+'WDIA.npz'

		nbi_pwr=['NB11_PNB','NB12_PNB','NB13_PNB','NB2A_PNB','NB2B_PNB','NB2C_PNB']
		nbi_vg =['NB11_VG1','NB12_VG1','NB13_VG1','NB2A_VG1','NB2B_VG1','NB2C_VG1']

		if not os.path.isfile(nfile1+'.gz'): 
			self.mds['ip'] = self.g.get(node_name1)
			self._get_save_zipf(nfile1,self.mds['ip'],2)
		else:
			data = self._get_save_zipf(nfile1,[],2)
			self.mds['ip'] = [data[0]]; self.mds['ip'].append(data[1])

		if not os.path.isfile(nfile2+'.gz'): 
			self.mds['da'] = self.g.get(node_name2)
			self._get_save_zipf(nfile2,self.mds['da'],2)
		else:
			data = self._get_save_zipf(nfile2,[],2)
			self.mds['da'] = [data[0]]; self.mds['da'].append(data[1])			

		if not os.path.isfile(nfile3+'.gz'): 
			self.mds['wdia'] = self.g.get(node_name3)
			self._get_save_zipf(nfile3,self.mds['wdia'],2)
		else:
			data = self._get_save_zipf(nfile3,[],2)
			self.mds['wdia'] = [data[0]]; self.mds['wdia'].append(data[1])				

		len_nbi = len(nbi_pwr)
		for i in range(len_nbi):
			nfile = self.savfile+'%s.npz'%nbi_pwr[i]
			if not os.path.isfile(nfile+'.gz'): 
				pnb = self.g.get('\\'+nbi_pwr[i])
				vg  = self.g.get('\\'+nbi_vg[i])
				self.mds['pnb%i'%(i+1)] = [pnb[0]];
				self.mds['pnb%i'%(i+1)].append(pnb[1]);
				self.mds['pnb%i'%(i+1)].append(vg[1]);
				self._get_save_zipf(nfile,self.mds['pnb%i'%(i+1)],3)
			else:
				data = self._get_save_zipf(nfile,[],3)
				self.mds['pnb%i'%(i+1)] = [data[0]]; self.mds['pnb%i'%(i+1)].append(data[1]);
				self.mds['pnb%i'%(i+1)].append(data[2])
		f = open(self.savfile+'/NBI.dat','w')
		f.write('NB11_PNB NB12_PNB NB13_PNB NB2A_PNB NB2B_PNB NB2C_PNB')
		f.close()

		lena = len(self.mds['da'][0]); lenv = len(self.mds['da'][1]); lend = lena- lenv
		if not lend == 0:
			if lend > 0: 
				for i in range(lend): self.mds['da'][1].append(0)
			else:
				for i in range(lend): self.mds['da'][0].append(0)
		return

	def _get_elm_peak(self):

		tavg = max(self.opt['ATS'],self.opt['ACES'])*3.

		tmin = (self.time-tavg)/1.e3
		tmax = (self.time+tavg)/1.e3
		ind1 = np.where(self.mds['da'][0]>tmin);
		ind2 = np.where(self.mds['da'][0][ind1]<tmax);
		ind3 = int(round(0.005/(self.mds['da'][0][4]-self.mds['da'][0][3]))) #0.003
		ind4 = int(round(0.00003/(self.mds['da'][0][4]-self.mds['da'][0][3])))

		ind3 = max(ind3,1); ind4 = max(ind4,1)

		tt = self.mds['da'][0][ind1][ind2]
		yy = self.mds['da'][1][ind1][ind2]
		lent = len(tt)

		self.opt['dacrit'] = float(self.note_in['e7'].get())
		self.opt['peaks']  = []

		mean_sig = np.mean(yy);
		if mean_sig>0: base = mean_sig*0.9; dbase = mean_sig*0.3;
		else: base = mean_sig*1.1; dbase = -mean_sig*0.3;
		dbase2 = np.max(yy)-base;
		dbase  = max(dbase,dbase2*0.2)

		ii=0
		while ii<lent:
			if (yy[ii]-base) > self.opt['dacrit']*dbase:
				self.opt['peaks'].append(tt[ii]*1.e3);
				ii += ind3
			else:
				ii += ind4

		self._draw_peaks()

		return

	def _get_year(self):
		self.ok,self.year = get_year(self.shotn)
		if not self.ok: print('>>> Unavail shot number'); exit()
		return

	def _load_eq(self):

		print('>>> Load EFIT...')
		self.efit_list = get_efit_list2(self.shotn)
		isefit = False
		for key in self.efit_list['isefit'].keys():
			if ((not isefit) and self.efit_list['isefit'][key]): isefit = True

		if not isefit: print('>>> No EFITs, exit!'); exit();

		self.gkfiles = dict()
		self.gkfiles['g'] = dict()
		self.gkfiles['k'] = dict()
		self.gkfiles['a'] = dict()
		self.gkfiles['eq']= dict()

		for efit_no in self.efit_list['isefit'].keys():
			self.gkfiles['g'][efit_no] = ''
			self.gkfiles['k'][efit_no] = ''
			if not self.efit_list['isefit'][efit_no]: continue
			ind1 = np.argmin(abs(self.efit_list['times'][efit_no]-self.time))
			t_time = self.efit_list['times'][efit_no][ind1]

			gfile_ori  = 'g%06i.%06i_%i_b'%(self.shotn,t_time,efit_no)
			gfile_sav  = 'g%06i.%06i_%i'%(self.shotn,t_time,efit_no)
			kfile_sav  = 'k%06i.%06i_%i'%(self.shotn,t_time,efit_no)
			afile_sav  = 'a%06i.%06i_%i'%(self.shotn,t_time,efit_no)
			gfile_dir  = self.efit_list['dirs'][efit_no] + '/g%06i.%06i'%(self.shotn,t_time)
			kfile_dir  = self.efit_list['dirs'][efit_no] + '/k%06i.%06i'%(self.shotn,t_time)
			afile_dir  = self.efit_list['dirs'][efit_no] + '/a%06i.%06i'%(self.shotn,t_time)

			if not os.path.isfile(gfile_sav): copyfile(gfile_dir,gfile_ori)
			if not os.path.isfile(kfile_sav): copyfile(kfile_dir,kfile_sav)
			if not os.path.isfile(afile_sav): copyfile(afile_dir,afile_sav)

			self.gkfiles['g'][efit_no] = gfile_sav;
			self.gkfiles['k'][efit_no] = kfile_sav;
			self.gkfiles['a'][efit_no] = _read_afile(afile_sav)

			if not os.path.isfile(gfile_sav):
				if efit_no < 3:
					os.system(efit_rmp+' '+kfile_sav)
					copyfile(kfile_sav,'kfile_run')
					efit_exec = efit_dir+'/efit65'
					os.system(efit_exec)
					gfile_sav2= 'g%06i.%06i'%(self.shotn,t_time)
					copyfile(gfile_sav2,gfile_sav)
				else:
					copyfile(gfile_ori,gfile_sav)

			eq = eqdsk.eqdsk(gfile_sav)
			self.gkfiles['eq'][efit_no] = eq
			self.gkfiles['eq'][efit_no].read_eqdsk(gfile_sav)
			self.gkfiles['eq'][efit_no].make_grid()
			self.gkfiles['eq'][efit_no].get_flux_contour()
			self.gkfiles['eq'][efit_no].make_rho_R_psin()

			self.gkfiles['eq'][efit_no].psif = interp2d(eq.R,eq.Z,(eq.psirz-eq.smag)/(eq.sbdy-eq.smag))
			rho_map   = eq.prhoR[:,1]
			psi_map   = eq.prhoR[:,0]
			Rlow = interp1d(eq.prhoR[:,0],eq.prhoR[:,2],'slinear')
			Rhigh= interp1d(eq.prhoR[:,0],eq.prhoR[:,3],'slinear')
			mu0  = 4.* np.pi * 1.e-7
			self.gkfiles['eq'][efit_no].jl = - Rlow(eq.psin) * eq.pp  - eq.ffp/mu0/Rlow(eq.psin)
			self.gkfiles['eq'][efit_no].jh = - Rhigh(eq.psin)* eq.pp  - eq.ffp/mu0/Rhigh(eq.psin)
			self.gkfiles['eq'][efit_no].ja = 0.5*(self.gkfiles['eq'][efit_no].jl+self.gkfiles['eq'][efit_no].jh)
			lenp      = len(psi_map)
			psi_map_ex= np.zeros(lenp+100)
			rho_map_ex= np.zeros(lenp+100)

			drdp = (rho_map[-1]-rho_map[-2])/(psi_map[-1]-psi_map[-2])

			for i in range(lenp+100):
				if i<lenp:
					psi_map_ex[i] = psi_map[i];
					rho_map_ex[i] = rho_map[i];
				else:
					psi_map_ex[i] = 1.+0.01*float(i-lenp+1);
					rho_map_ex[i] = 1.+drdp*(psi_map_ex[i]-1.);

			self.gkfiles['eq'][efit_no].psi_to_rho = interp1d(psi_map_ex,rho_map_ex)
		return

	def _get_xmap(self):

		if self.ts['core']['nch']>3:
			self.rr['ts']['core'][1] = self.ts['core']['rr']
			self.rr['ts']['core'][2] = self.gkfiles['eq'][self.efit_no_g.get()].psif(self.rr['ts']['core'][1],0.)
			self.rr['ts']['core'][3] = self.gkfiles['eq'][self.efit_no_g.get()].psi_to_rho(self.rr['ts']['core'][2])
	
			self.rr['ts']['edge'][1] = self.ts['edge']['rr']
			self.rr['ts']['edge'][2] = self.gkfiles['eq'][self.efit_no_g.get()].psif(self.rr['ts']['edge'][1],0.)
			self.rr['ts']['edge'][3] = self.gkfiles['eq'][self.efit_no_g.get()].psi_to_rho(self.rr['ts']['edge'][2])
		if self.ces['nch']>3:
			self.rr['ces'][1] = self.ces['rr']
			self.rr['ces'][2] = self.gkfiles['eq'][self.efit_no_g.get()].psif(self.rr['ces'][1],0.)
			self.rr['ces'][3] = self.gkfiles['eq'][self.efit_no_g.get()].psi_to_rho(self.rr['ces'][2])				

		return

	def _diagnostics_exceptions(self):

		if (self.shotn >= 21760 and self.shotn <= 24136):
			self.ts['core']['rr'] = [1.806, 1.826, 1.848, 1.871, 1.894, 1.917, 1.942, 1.966, 1.991, 2.016, 2.041, 2.068, 2.093, 2.120]
			self.ts['edge']['rr'] = [2.124, 2.137, 2.143, 2.149, 2.156, 2.162, 2.177, 2.191, 2.202, 2.216, 2.229, 2.242, 2.257, 2.271, 2.285, 2.297, 2.311]

		if self.year == 2011:
			self.ces['rr'] = [1.795,1.800,1.850,1.900,1.950,2.000,2.050,2.100,2.140,2.160,2.170,2.180,2.190,2.200,2.205,2.210,2.215,2.220,2.225,2.230,2.235,2.240,2.245,2.250,2.255,2.265,2.275,2.280,2.285,2.290,2.295,2.300] # 2011
		elif self.year == 2012:	
			self.ces['rr'] = [1.800,1.850,1.900,1.950,2.000,2.050,2.100,2.150,2.170,2.180,2.190,2.200,2.205,2.210,2.215,2.220,2.225,2.230,2.235,2.240,2.245,2.250,2.255,2.260,2.265,2.270,2.275,2.280,2.285,2.290,2.295,2.300] # 2012
		elif (self.year == 2013 or self.year == 2014):	
			self.ces['rr'] = [1.795,1.850,1.900,1.950,2.000,2.050,2.100,2.150,2.170,2.180,2.190,2.200,2.205,2.210,2.215,2.220,2.225,2.230,2.235,2.240,2.245,2.250,2.255,2.260,2.265,2.270,2.275,2.280,2.285,2.290,2.295,2.300] # 2013/2014
		elif self.year == 2015:	
			self.ces['rr'] = [1.801,1.822,1.843,1.874,1.895,1.945,1.995,2.016,2.047,2.078,2.099,2.125,2.150,2.171,2.192,2.203,2.213,2.223,2.228,2.233,2.238,2.243,2.248,2.253,2.259,2.264,2.269,2.273,2.280,2.286,2.291,2.296] # 2015
		return

	def _make_home_canvas(self):

		self.figure['name']['home'] = plt.figure(1,figsize=self.figure['size']['home'])
		self.figure['name']['home'].tight_layout()

		self.figure['canvas']['home'] = FigureCanvasTkAgg(self.figure['name']['home'],master=self.root)

		self.figure['widget']['home'] = self.figure['canvas']['home'].get_tk_widget()
		self.figure['widget']['home'].grid(row=1,column=4,rowspan=2,sticky='w')

		toolbar_frame = tk.Frame(self.root)
		toolbar_frame.grid(row=0,column=4)
		
		self.figure['toolbar']['home'] = NavigationToolbar2Tk(self.figure['canvas']['home'],toolbar_frame)

		return

	def _make_note_frame(self):

		titles = ['MDS']
		self.titles = titles
		self.home['nb'] = ttk.Notebook(self.root,width=400,height=720)
		for i in range(1,len(titles)+1):
			self.home['page%i'%i] = ttk.Frame(self.home['nb'])
			self.home['nb'].add(self.home['page%i'%i],text=titles[i-1])
			if i > 1: self.home['nb'].tab(i-1,state='disabled')
			for k in range(10): self.home['page%i'%i].grid_columnconfigure(k,weight=1)

		self.home['nb'].grid(row=1,column=0,columnspan=4,sticky='n')

		return

	def _make_input_frame(self):

		irow = 0
		self.l1 = tk.Label(self.home['page1'], text="===================== DIAGNOSTICS ======================",justify='center')
		self.l1.grid(row=irow,column=0,columnspan=10,pady=5)
		self.note_in['frame1'] = tk.Frame(self.home['page1'])
		self.note_in['frame1'].grid(row=irow+1,column=0,columnspan=10,padx=20)
		self.note_in['b5'] = tk.Button(self.note_in['frame1'],  text="EFIT",   width = 2,command=lambda: self._efit_page('efit'))
		self.note_in['b1'] = tk.Button(self.note_in['frame1'],  text="CES-TI", width = 4,command=lambda: self._ces_page('ti'))
		self.note_in['b2'] = tk.Button(self.note_in['frame1'],  text="CES-VT", width = 4,command=lambda: self._ces_page('vt'))
		self.note_in['b3'] = tk.Button(self.note_in['frame1'],  text="TS-TE",  width = 4,command=lambda: self._ts_page('te'))
		self.note_in['b4'] = tk.Button(self.note_in['frame1'],  text="TS-NE",  width = 4,command=lambda: self._ts_page('ne'))
		self.note_in['b5'].pack(side='left',padx=5); self.note_in['b1'].pack(side='left',padx=5); self.note_in['b2'].pack(side='left',padx=5);
		self.note_in['b3'].pack(side='left',padx=5); self.note_in['b4'].pack(side='left',padx=5);

		irow += 2
		self.l2 = tk.Label(self.home['page1'], text="====================== TIME INFO. =======================",justify='center')
		self.l2.grid(row=irow,column=0,columnspan=10,pady=5)
		self.note_in['frame2'] = tk.Frame(self.home['page1'])
		self.note_in['frame2'].grid(row=irow+1,column=0,columnspan=10,padx=20)
		self.l3 = tk.Label(self.note_in['frame2'], text="Time [ms]",justify='center')
		self.l3.pack(side='left')
		self.note_in['e1'] = tk.Entry(self.note_in['frame2'],width=8,justify='center')
		self.note_in['e1'].insert(10,'%i'%self.time)
		self.note_in['e1'].pack(side='left')
		self.l4 = tk.Label(self.note_in['frame2'], text="Avg. [ms]",justify='center')
		self.l4.pack(side='left',padx=10)
		self.note_in['e2'] = tk.Entry(self.note_in['frame2'],width=8,justify='center')
		self.note_in['e2'].insert(10,'100')
		self.note_in['e2'].pack(side='left')
		b5 = tk.Button(self.note_in['frame2'],  text="SET", width = 5,command=lambda: self._set_time())
		b5.pack(side='right',padx=5);

		irow += 5
		self.l2 = tk.Label(self.home['page1'], text="====================== TIME LIST. =======================",justify='center')
		self.l2.grid(row=irow,column=0,columnspan=10,pady=5)		
		self.l5 = tk.Label(self.home['page1'], text="== AVAILABLE [ms] ==",justify='center')
		self.l5.grid(row=irow+1,column=0,columnspan=5,pady=5,padx=45,sticky='w')
		self.l6 = tk.Label(self.home['page1'], text="== SELECTED [ms] ==",justify='center')
		self.l6.grid(row=irow+1,column=5,columnspan=5,pady=5,padx=45,sticky='e')

		irow += 2
		self.note_in['frame3'] = tk.Frame(self.home['page1'])
		self.note_in['frame3'].grid(row=irow,column=0,columnspan=5,padx=20,sticky='w')
		self.note_in['s1'] = tk.Scrollbar(self.note_in['frame3'])
		self.note_in['s1'].pack(side='right',fill='y')
		self.note_in['l1'] = tk.Listbox(self.note_in['frame3'],yscrollcommand = self.note_in['s1'].set,height=7)
		self.note_in['l1'].pack(side='left',fill='x')
		self.note_in['s1']["command"] = self.note_in['l1'].yview
		self.note_in['l1'].bind('<ButtonRelease-3>',lambda x: self._select_time(list_no=1))
		self.note_in['l1'].bind('<Double-1>',lambda x: self._click_list1())
			
		self.note_in['frame4'] = tk.Frame(self.home['page1'])
		self.note_in['frame4'].grid(row=irow,column=5,columnspan=5,padx=20,sticky='e')
		self.note_in['s2'] = tk.Scrollbar(self.note_in['frame4'])
		self.note_in['s2'].pack(side='right',fill='y')
		self.note_in['l2'] = tk.Listbox(self.note_in['frame4'],yscrollcommand = self.note_in['s2'].set,height=7)
		self.note_in['l2'].pack(side='left',fill='x')
		self.note_in['s2']["command"] = self.note_in['l2'].yview
		self.note_in['l2'].bind('<ButtonRelease-3>',lambda x: self._click_list2())
		self.note_in['l2'].bind('<Double-1>',lambda x: self._select_time(list_no=2))

		self.note_in['frame5'] = tk.Frame(self.home['page1'])
		self.note_in['frame5'].grid(row=irow+1,column=3,columnspan=7,padx=20,sticky='e')
		b6 = tk.Button(self.note_in['frame5'],  text="SELECT ALL", width = 8,command=lambda: self._select_all())
		b7 = tk.Button(self.note_in['frame5'],  text="DELETE ALL", width = 8,command=lambda: self._remove_all())
		b6.pack(side='left',pady=5); b7.pack(side='left',pady=5);

		irow += 2
		self.l11 = tk.Label(self.home['page1'], text="===================== ELM SETUP ======================",justify='center')
		self.l11.grid(row=irow,column=0,columnspan=10,pady=5)
		self.note_in['frame9'] = tk.Frame(self.home['page1'])
		self.note_in['frame9'].grid(row=irow+1,column=0,columnspan=10,padx=20,sticky='w')

		self.l12 = tk.Label(self.note_in['frame9'], text="Crit",justify='center')
		self.l12.pack(side='left')
		self.note_in['e7'] = tk.Entry(self.note_in['frame9'],width=6,justify='center')
		self.note_in['e7'].insert(10,'%3.1f'%self.opt['dacrit'])
		self.note_in['e7'].pack(side='left',padx=5)

		self.l13 = tk.Label(self.note_in['frame9'], text="Duty [%]",justify='center')
		self.l13.pack(side='left')
		self.note_in['e8'] = tk.Entry(self.note_in['frame9'],width=6,justify='center')
		self.note_in['e8'].insert(10,'%i'%self.opt['duty'])
		self.note_in['e8'].pack(side='left',padx=5)

		b11 = tk.Button(self.note_in['frame9'],  text="FIND", width = 4,command=lambda: self._get_elm_peak())
		b11.pack(side='left',padx=5);
		b12 = tk.Button(self.note_in['frame9'],  text="PRESELECT", width = 9,command=lambda: self._preselect())
		b12.pack(side='left',padx=5);

		irow += 2
		self.l7 = tk.Label(self.home['page1'], text="====================== PLOT SETUP =======================",justify='center')
		self.l7.grid(row=irow,column=0,columnspan=10,pady=5)
		self.note_in['frame6'] = tk.Frame(self.home['page1'])
		self.note_in['frame6'].grid(row=irow+1,column=0,columnspan=10,padx=20,sticky='w')

		self.l8 = tk.Label(self.note_in['frame6'], text="XMAP",justify='center')
		self.l8.pack(side='left',padx=5)
		self.note_in['r1']=tk.Radiobutton(self.note_in['frame6'], text="R[m]",value=1, variable=self.note_in['xmap'], command=lambda: self._change_xmap())
		self.note_in['r2']=tk.Radiobutton(self.note_in['frame6'], text="psi", value=2, variable=self.note_in['xmap'], command=lambda: self._change_xmap())
		self.note_in['r3']=tk.Radiobutton(self.note_in['frame6'], text="rho", value=3, variable=self.note_in['xmap'], command=lambda: self._change_xmap())
		self.note_in['r1'].pack(side='left'); self.note_in['r2'].pack(side='left'); self.note_in['r3'].pack(side='left');

		b8 = tk.Button(self.note_in['frame6'],  text="SET", width = 4,command=lambda: self._set_range())
		b9 = tk.Button(self.note_in['frame6'],  text="GENERATE", width = 8,command=lambda: self._generate())
		b8.pack(side='left',padx=5);b9.pack(side='left',padx=5);
	

		irow += 2
		self.note_in['frame7'] = tk.Frame(self.home['page1'])
		self.note_in['frame7'].grid(row=irow,column=0,columnspan=10,padx=20,sticky='w')

		self.l9 = tk.Label(self.note_in['frame7'], text="XLIM[min/max]",justify='center')
		self.l9.pack(side='left')
		self.note_in['e3'] = tk.Entry(self.note_in['frame7'],width=6,justify='center')
		self.note_in['e3'].insert(10,'1.5')
		self.note_in['e3'].pack(side='left')
		self.note_in['e4'] = tk.Entry(self.note_in['frame7'],width=6,justify='center')
		self.note_in['e4'].insert(10,'2.4')
		self.note_in['e4'].pack(side='left')

		self.l10 = tk.Label(self.note_in['frame7'], text="YLIM[min/max]",justify='center')
		self.l10.pack(side='left')
		self.note_in['e5'] = tk.Entry(self.note_in['frame7'],width=6,justify='center')
		self.note_in['e5'].insert(10,'0.0')
		self.note_in['e5'].pack(side='left')
		self.note_in['e6'] = tk.Entry(self.note_in['frame7'],width=6,justify='center')
		self.note_in['e6'].insert(10,'5.0')
		self.note_in['e6'].pack(side='left')

		self.note_in['frame8'] = tk.Frame(self.home['page1'])
		self.note_in['frame8'].grid(row=irow+1,column=0,columnspan=10,pady=20)		
		b10 = tk.Button(self.note_in['frame8'],  text="DONE", width = 10,command=lambda: self._done())
		b10.pack(side='left');		

		irow = 4
		self.l14 = tk.Label(self.home['page1'], text="===================== EQU SETUP ======================",justify='center')
		self.l14.grid(row=irow,column=0,columnspan=10,pady=5)
		self.note_in['frame10'] = tk.Frame(self.home['page1'])
		self.note_in['frame10'].grid(row=irow+1,column=0,columnspan=10,padx=20,sticky='w')
	
		self.l15 = tk.Label(self.note_in['frame10'], text='[EFITg]')
		self.l15.pack(side='left',padx=6)
		self.note_in['r1'] = tk.Radiobutton(self.note_in['frame10'], text="01[MAG]", value=1, variable=self.efit_no_g)
		self.note_in['r2'] = tk.Radiobutton(self.note_in['frame10'], text="02[MSE]", value=2, variable=self.efit_no_g)
		self.note_in['r3'] = tk.Radiobutton(self.note_in['frame10'], text="03[MSE]", value=3, variable=self.efit_no_g)
		self.note_in['r4'] = tk.Radiobutton(self.note_in['frame10'], text="04[MAG+]", value=4, variable=self.efit_no_g)
		self.note_in['r5'] = tk.Radiobutton(self.note_in['frame10'], text="05[MSE+]", value=5, variable=self.efit_no_g)
	
		self.note_in['r1'].pack(side='left')
		self.note_in['r2'].pack(side='left')
#		self.note_in['r3'].pack(side='left')
		self.note_in['r4'].pack(side='left')
		self.note_in['r5'].pack(side='left')

		self.note_in['frame11'] = tk.Frame(self.home['page1'])
		self.note_in['frame11'].grid(row=irow+2,column=0,columnspan=10,padx=20,sticky='w')
		
		self.l16 = tk.Label(self.note_in['frame11'], text='[EFITk]')
		self.l16.pack(side='left',padx=6)
		self.note_in['r6'] = tk.Radiobutton(self.note_in['frame11'], text="01[MAG]", value=1, variable=self.efit_no_k)
		self.note_in['r7'] = tk.Radiobutton(self.note_in['frame11'], text="02[MSE]", value=2, variable=self.efit_no_k)
		self.note_in['r8'] = tk.Radiobutton(self.note_in['frame11'], text="03", value=3, variable=self.efit_no_k)
		self.note_in['r9'] = tk.Radiobutton(self.note_in['frame11'], text="04[MAG+]", value=4, variable=self.efit_no_k)
		self.note_in['r10']= tk.Radiobutton(self.note_in['frame11'], text="05[MSE+]", value=5, variable=self.efit_no_k)

		self.note_in['r6'].pack(side='left')
		self.note_in['r7'].pack(side='left')
#		self.note_in['r8'].pack(side='left')
		self.note_in['r9'].pack(side='left')
		self.note_in['r10'].pack(side='left')
		for i in range(1,11): self.note_in['r%i'%i].config(state='disabled')
		return

	def _set_time(self):

		self.opt['T%s'%self.opt['cpage'].upper()] = int(self.note_in['e1'].get())
		self.opt['A%s'%self.opt['cpage'].upper()] = int(self.note_in['e2'].get())

		tmin = self.opt['T%s'%self.opt['cpage'].upper()] - self.opt['A%s'%self.opt['cpage'].upper()]*0.8
		tmax = self.opt['T%s'%self.opt['cpage'].upper()] + self.opt['A%s'%self.opt['cpage'].upper()]*0.8

		tmin = tmin*1.e-3; tmax = tmax*1.e-3

		if not (self.opt['cpage']=='' or self.opt['cpage']=='tci' or self.opt['cpage']=='efit'):
			self.figure['name']['home'].canvas.draw_idle()
			self.figure['axes']['home'][0].set_xlim([tmin,tmax])

		diag = self.opt['cpage']
		if not (diag == 'ces' or diag == 'ts'): return
		self._get_list(diag)
		self._get_elm_peak()
		return

	def _get_list(self,diag):

		tmin = self.opt['T%s'%self.opt['cpage'].upper()] - self.opt['A%s'%self.opt['cpage'].upper()]*0.5
		tmax = self.opt['T%s'%self.opt['cpage'].upper()] + self.opt['A%s'%self.opt['cpage'].upper()]*0.5	
		self.note_in['l1'].delete(0,'end')
		if diag == 'ces': 
			ind1 = np.where(self.ces['ti']['val'][1][0]>=tmin/1.e3); ind2 = np.where(self.ces['ti']['val'][1][0][ind1]<=tmax/1.e3)
			tt   = self.ces['ti']['val'][1][0][ind1][ind2]
		elif diag== 'ts': 
			ind1 = np.where(self.ts['te']['core'][1][0]>=tmin/1.e3); ind2 = np.where(self.ts['te']['core'][1][0][ind1]<=tmax/1.e3)
			tt   = self.ts['te']['core'][1][0][ind1][ind2]
		for time in tt: self.note_in['l1'].insert('end','%i'%(time*1.e3))
		self.opt['isdata'][diag] = len(tt)>0
		return

	def _sync_list(self,save=True):

		if (self.opt['cpage']=='' or self.opt['cpage']=='tci' or self.opt['cpage']=='efit'): return
		diag = self.opt['cpage']; ch = self.opt['c_ch'];
		if save: 
			self.opt['slist'][diag][ch] = copy.deepcopy(list(self.note_in['l2'].get(0,'end')))
		else:
			self.note_in['l2'].delete(0,'end')
			for item in self.opt['slist'][diag][ch]:
				self.note_in['l2'].insert('end',item)

	def _click_list1(self):

		self._clear_select()
		selection = self.note_in['l1'].curselection()
		if len(selection) == 0: return
		line = self.note_in['l1'].get(selection[0])

		lists_old = list(self.note_in['l2'].get(0,'end'));
		lists     = list(self.note_in['l2'].get(0,'end'));
		self.note_in['l2'].delete(0,'end');
		if not line in lists: lists.append(line)
		lists.sort()
		self.figure['name']['home'].canvas.draw_idle()
		newlist = False
		for item in lists:
			self.note_in['l2'].insert('end',item)
			if not item in lists_old:
				time = int(item)
				self._draw_time(time)
			newlist = True;
		if newlist: self._reset_gui()
		return		

	def _click_list2(self):

		self._clear_select()
		selection = self.note_in['l2'].curselection()
		if len(selection) == 0: return
		item = self.note_in['l2'].get(selection[0])
		time = int(item)
		self.note_in['l2'].delete(selection[0]);
		self._delete_time(time)
		self._reset_gui()
		return

	def _select_time(self,list_no=2):

		self._clear_select()
		if list_no==1:
			selection = self.note_in['l1'].curselection()
			if len(selection) == 0: return
			line = self.note_in['l1'].get(selection[0])
			color = 'g'
		else:
			selection = self.note_in['l2'].curselection()
			if len(selection) == 0: return
			line = self.note_in['l2'].get(selection[0])
			color = 'r'			

		time = int(line)
		if self.opt['cpage']=='ces': 
			data = self.ces[self.opt['c_ch']]
			xx = self.rr['ces'][self.note_in['xmap'].get()]
			nch = self.ces['nch']
			self._select_time_ces(time,xx,data,nch,color)
			

		if self.opt['cpage']=='ts': 			
			data = self.ts[self.opt['c_ch']]
			xx1  = self.rr['ts']['core'][self.note_in['xmap'].get()]
			nch1 = self.ts['core']['nch']
			xx2  = self.rr['ts']['edge'][self.note_in['xmap'].get()]
			nch2 = self.ts['edge']['nch']			
			self._select_time_ts(time,xx1,nch1,xx2,nch2,data,color)			
		self._set_range()
		return		

	def _select_time_ces(self,time,xx,data,nch,color):

		self.figure['name']['home'].canvas.draw_idle()

		if self.opt['c_ch'] == 'ti': factor = 1.e-3; 
		else: factor = 1.e0;

		ind = np.argmin(abs(data['val'][1][0]-time/1.e3))
		yy = []; yerr= [];
		for i in range(1,nch+1):
			yy.append(data['val'][i][1][ind]);
			yerr.append(data['err'][i][1][ind]);
		np.nan_to_num(yy, copy=False); np.nan_to_num(yerr, copy=False)
		yy = np.array(yy)*factor; yerr = np.array(yerr)*factor
		line1 = self.figure['axes']['home'][1].errorbar(xx,yy,yerr,fmt='o',c=color)
		self.figure['pegend']['stime'].append(line1)
		self.figure['legend']['stime'].append(time)

		line2 = self.figure['axes']['home'][0].axvline(x=time/1.e3,c=color,linestyle='--')
		self.figure['pegend']['stime'].append(line2)
		self.figure['legend']['stime'].append(time)

		return

	def _select_time_ts(self,time,xx1,nch1,xx2,nch2,data,color):

		self.figure['name']['home'].canvas.draw_idle()

		if self.opt['c_ch'] == 'te': factor = 1.e-3; 
		else: factor = 1.e-19;

		ind = np.argmin(abs(data['core'][1][0]-time/1.e3))
		xx = []; yy = []; yerr= [];
		for i in range(1,nch1+1):
			xx.append(xx1[i-1])
			yy.append(data['core'][i][1][ind]);
			yerr.append(data['core_err'][i][1][ind]);
		for i in range(1,nch2+1):
			xx.append(xx2[i-1])
			yy.append(data['edge'][i][1][ind]);
			yerr.append(data['edge_err'][i][1][ind]);			
		np.nan_to_num(yy, copy=False); np.nan_to_num(yerr, copy=False)
		yy = np.array(yy)*factor; yerr = np.array(yerr)*factor
		line1 = self.figure['axes']['home'][1].errorbar(xx,yy,yerr,fmt='o',c=color)
		self.figure['pegend']['stime'].append(line1)
		self.figure['legend']['stime'].append(time)

		line2 = self.figure['axes']['home'][0].axvline(x=time/1.e3,c=color,linestyle='--')
		self.figure['pegend']['stime'].append(line2)
		self.figure['legend']['stime'].append(time)
		return				

	def _clear_select(self):

		lenp = len(self.figure['legend']['stime'])
		if lenp>0:
			self.figure['pegend']['stime'][0].remove()
			self.figure['pegend']['stime'][1].remove()
			self.figure['pegend']['stime'] = []
			self.figure['legend']['stime'] = []
		return

	def _select_all(self):

		self._clear_select()
		list1 = list(self.note_in['l1'].get(0,'end'));
		list2 = list(self.note_in['l2'].get(0,'end'));
		list2_old = list(self.note_in['l2'].get(0,'end'));
		self.note_in['l2'].delete(0,'end');

		for item in list1:
			if not item in list2: list2.append(item)
		list2.sort()
		self.figure['name']['home'].canvas.draw_idle()
		for item in list2:
			self.note_in['l2'].insert('end',item)
			if not item in list2_old:
				time = int(item)
				self._draw_time(time)

		self._reset_gui()
		return		

	def _remove_all(self):

		self._clear_select()
		self.note_in['l2'].delete(0,'end')
		times = copy.deepcopy(self.figure['legend']['time'])
		for time in times:
			self._delete_time(time)

		self._reset_gui()
		return

	def _get_duty(self,time):

		lenp = len(self.opt['peaks'])
		if  lenp< 3: return 1000.

		for i in range(lenp-1):
			left = time - self.opt['peaks'][i]
			right = time - self.opt['peaks'][i+1]
			if left*right <=0.: break

		if i==(lenp-2): return 2000.

		duty = -100.*(right)/(self.opt['peaks'][i+1]-self.opt['peaks'][i]);

		if duty <5.: return 3000;
		else: return duty

	def _preselect(self):

		self._remove_all()
		list1 = list(self.note_in['l1'].get(0,'end'));
		list2 = list(self.note_in['l2'].get(0,'end'));
		list2_old = list(self.note_in['l2'].get(0,'end'));
		self.note_in['l2'].delete(0,'end');

		self.opt['duty'] = float(self.note_in['e8'].get())
		for item in list1:
			duty = self._get_duty(float(item))
			if ((not item in list2) and duty<=self.opt['duty']): list2.append(item)

		list2.sort()
		self.figure['name']['home'].canvas.draw_idle()
		for item in list2:
			self.note_in['l2'].insert('end',item)
			if not item in list2_old:
				time = int(item)
				self._draw_time(time)

		self._reset_gui()

		return

	def _draw_efit(self):

		llegend = []
		for i in range(1,6):
			if not self.efit_list['isefit'][i]: continue
			eq = self.gkfiles['eq'][i]
			self.figure['axes']['home'][0].plot(eq.rzbdy[:,0],eq.rzbdy[:,1],linestyle='--')
			self.figure['axes']['home'][0].set_xlabel('R [m]')
			self.figure['axes']['home'][0].set_ylabel('Z [m]')
			
			self.figure['axes']['home'][1].plot(eq.psin,eq.pres/1.e3,linestyle='-')
			self.figure['axes']['home'][1].set_xlabel('$\psi_N$')
			self.figure['axes']['home'][1].set_ylabel('P [kPa]')

			self.figure['axes']['home'][2].plot(eq.psin[:-2],eq.q[:-2],linestyle='-')
			self.figure['axes']['home'][2].set_xlabel('$\psi_N$')
			self.figure['axes']['home'][2].set_ylabel('q')

			self.figure['axes']['home'][3].plot(eq.psin,eq.ja/1.e6,linestyle='-')
			self.figure['axes']['home'][3].set_xlabel('$\psi_N$')
			self.figure['axes']['home'][3].set_ylabel('$j_{\phi} [MA/m^2]$')
			llegend.append('EFIT%02i'%i)
	
		self.figure['axes']['home'][0].legend(llegend)
		self.figure['axes']['home'][1].legend(llegend)
		self.figure['axes']['home'][2].legend(llegend)
		self.figure['axes']['home'][3].legend(llegend)

		glist = [1,2,4,5]
		line = '{:8s}'.format('')
		for i in glist: line = line + '{:7s} '.format(' EFIT%02i'%i)
		line = line + '\n{:8s}'.format('chisq')
		for i in glist: 
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['tsaisq'])
			else: line = line + '{:7s} '.format('  -')
		line = line + '\n{:8s}'.format('betap')
		for i in glist: 
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['betap'])
			else: line = line + '{:7s} '.format('  -')
		line = line + '\n{:8s}'.format('betan')
		for i in glist:
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['betan'])
			else: line = line + '{:7s} '.format('  -')
		line = line + '\n{:8s}'.format('wmhd')
		for i in glist:
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['wplasm']/1.e3)
			else: line = line + '{:7s} '.format('  -')
		line = line + '\n{:8s}'.format('q95')
		for i in glist:
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['qpsib'])
			else: line = line + '{:7s} '.format('  -')
		line = line + '\n{:8s}'.format('q0')
		for i in glist:
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['qqmagx'])
			else: line = line + '{:7s} '.format('  -')
		line = line + '\n{:8s}'.format('Raxis')
		for i in glist:
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['rmagx']/1.e2)
			else: line = line + '{:7s} '.format('  -')
		line = line + '\n{:8s}'.format('drsep')
		for i in glist:
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['ssep'])
			else: line = line + '{:7s} '.format('  -')
		line = line + '\n{:8s}'.format('ip')
		for i in glist:
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['cpasma']/1.e3)
			else: line = line + '{:7s} '.format('  -')
		line = line + '\n{:8s}'.format('li')
		for i in glist: 
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['ali'])
			else: line = line + '{:7s} '.format('  -')
		line = line + '\n{:8s}'.format('bcentr')
		for i in glist: 
			if self.efit_list['isefit'][i]:line = line + '{:7.2f} '.format(self.gkfiles['a'][i]['bcentr'])
			else: line = line + '{:7s} '.format('  -')

	
		self.figure['axes']['home'][4].text(-0.1,0.3,line,horizontalalignment='left',transform = self.figure['axes']['home'][4].transAxes,linespacing=1.2,family='monospace')
		self.figure['axes']['home'][4].set_xlim([0, 1])
		self.figure['axes']['home'][4].set_ylim([0, 1])
		self.figure['axes']['home'][4].axis('off')

		for i in range(1,6):
			if not self.efit_list['isefit'][i]: continue
			eq = self.gkfiles['eq'][i]
			self.figure['axes']['home'][0].scatter(eq.rmag,eq.zmag,marker='x')

		self.figure['name']['home'].tight_layout()
		return

	def _draw_mds(self):

		self.figure['pegend']['mds']   = []; self.figure['legend']['mds']   = [];
		self.figure['pegend']['peak']  = []; self.figure['legend']['peak']  = [];
		self.figure['pegend']['time']  = []; self.figure['legend']['time']  = [];
		self.figure['pegend']['time2'] = []; self.figure['legend']['time2'] = [];
		self.figure['pegend']['stime'] = []; self.figure['legend']['stime'] = [];

		plot1, = self.figure['axes']['home'][0].plot(self.mds['ip'][0],  self.mds['ip'][1]/-1.e5)

		t0 = float(self.note_in['e1'].get())/1.e3
		try: dt = float(self.note_in['e2'].get())/1.e3
		except: dt = 0.3;
		
		ip_max = max(abs(self.mds['ip'][1]))/1.e5

		tind   = np.where(abs(self.mds['da'][0]-t0) < dt)
		da_max = max(self.mds['da'][1][tind]);
		da_off = np.mean(self.mds['da'][1][tind])
		da_sca = (da_max-da_off)*2./ip_max
		if da_off > 0: da_off = 0.8*da_off
		else: da_off = 1.2*da_off
		
		plot2, = self.figure['axes']['home'][0].plot(self.mds['da'][0],  (self.mds['da'][1]-da_off)/da_sca)
		plot3, = self.figure['axes']['home'][0].plot(self.mds['pnb1'][0],self.mds['pnb1'][1]/1.e0)
		plot4, = self.figure['axes']['home'][0].plot(self.mds['pnb2'][0],self.mds['pnb2'][1]/1.e0)
		plot5, = self.figure['axes']['home'][0].plot(self.mds['pnb3'][0],self.mds['pnb3'][1]/1.e0)

		self.figure['pegend']['mds'].append(plot1); self.figure['legend']['mds'].append('IP[100kA]')
		self.figure['pegend']['mds'].append(plot2); self.figure['legend']['mds'].append('DA')
		self.figure['pegend']['mds'].append(plot3); self.figure['legend']['mds'].append('NB1A [MW]')
		self.figure['pegend']['mds'].append(plot4); self.figure['legend']['mds'].append('NB1B [MW]')
		self.figure['pegend']['mds'].append(plot5); self.figure['legend']['mds'].append('NB1C [MW]')
		self.figure['axes']['home'][0].legend(self.figure['pegend']['mds'],self.figure['legend']['mds'],loc='upper right')

		self.figure['axes']['home'][0].set_ylim(0.,ip_max*1.1)
		self.figure['axes']['home'][0].set_xlabel('time [s]')

		if    self.note_in['xmap'].get() ==1: self.figure['axes']['home'][1].set_xlabel('R[m]')
		elif  self.note_in['xmap'].get() ==2: self.figure['axes']['home'][1].set_xlabel('$\\psi_N$')
		else:                             self.figure['axes']['home'][1].set_xlabel('$\\rho_N$')
		
		self.figure['name']['home'].tight_layout()
		return

	def _draw_peaks(self):

		self.figure['name']['home'].canvas.draw_idle()
		for item in self.figure['pegend']['peak']: item.remove();
		self.figure['pegend']['peak'] = []; self.figure['legend']['peak'] = [];

		for time in self.opt['peaks']:
			line2 = self.figure['axes']['home'][0].axvline(x=time/1.e3,c='gold',linestyle='--')
			self.figure['pegend']['peak'].append(line2)
			self.figure['legend']['peak'].append(time)
		return

	def _draw_time(self,time):

		if self.opt['cpage']=='ces': 
			data = self.ces[self.opt['c_ch']]
			xx = self.rr['ces'][self.note_in['xmap'].get()]
			nch = self.ces['nch']
			self._draw_time_ces(time,xx,data,nch)
			
		if self.opt['cpage']=='ts': 			
			data = self.ts[self.opt['c_ch']]
			xx1  = self.rr['ts']['core'][self.note_in['xmap'].get()]
			nch1 = self.ts['core']['nch']
			xx2  = self.rr['ts']['edge'][self.note_in['xmap'].get()]
			nch2 = self.ts['edge']['nch']			
			self._draw_time_ts(time,xx1,nch1,xx2,nch2,data)
		self._set_range()
		return

	def _draw_time_ces(self,time,xx,data,nch):

		self.figure['name']['home'].canvas.draw_idle()

		if self.opt['c_ch'] == 'ti': factor = 1.e-3; 
		else: factor = 1.e0;

		ind = np.argmin(abs(data['val'][1][0]-time/1.e3))
		yy = []; yerr= [];
		for i in range(1,nch+1):
			yy.append(data['val'][i][1][ind]);
			yerr.append(data['err'][i][1][ind]);
		np.nan_to_num(yy, copy=False); np.nan_to_num(yerr, copy=False)
		yy = np.array(yy)*factor; yerr = np.array(yerr)*factor
		line1 = self.figure['axes']['home'][1].errorbar(xx,yy,yerr,fmt='o',c='b')
		self.figure['pegend']['time'].append(line1)
		self.figure['legend']['time'].append(time)

		line2 = self.figure['axes']['home'][0].axvline(x=time/1.e3,c='b',linestyle='--')
		self.figure['pegend']['time2'].append(line2)
		self.figure['legend']['time2'].append(time)

		return

	def _draw_time_ts(self,time,xx1,nch1,xx2,nch2,data):

		self.figure['name']['home'].canvas.draw_idle()

		if self.opt['c_ch'] == 'te': factor = 1.e-3; 
		else: factor = 1.e-19;

		ind = np.argmin(abs(data['core'][1][0]-time/1.e3))
		xx = []; yy = []; yerr= [];
		for i in range(1,nch1+1):
			xx.append(xx1[i-1])
			yy.append(data['core'][i][1][ind]);
			yerr.append(data['core_err'][i][1][ind]);
		for i in range(1,nch2+1):
			xx.append(xx2[i-1])
			yy.append(data['edge'][i][1][ind]);
			yerr.append(data['edge_err'][i][1][ind]);			
		np.nan_to_num(yy, copy=False); np.nan_to_num(yerr, copy=False)
		yy = np.array(yy)*factor; yerr = np.array(yerr)*factor
		line1 = self.figure['axes']['home'][1].errorbar(xx,yy,yerr,fmt='o',c='b')
		self.figure['pegend']['time'].append(line1)
		self.figure['legend']['time'].append(time)

		line2 = self.figure['axes']['home'][0].axvline(x=time/1.e3,c='b',linestyle='--')
		self.figure['pegend']['time2'].append(line2)
		self.figure['legend']['time2'].append(time)
		return		

	def _set_range(self):

		self.figure['name']['home'].canvas.draw_idle()
		xmin = float(self.note_in['e3'].get())
		xmax = float(self.note_in['e4'].get())
		ymin = float(self.note_in['e5'].get())
		ymax = float(self.note_in['e6'].get())

		ch = self.opt['c_ch']
		self.opt['range'][ch]['xmin'] = xmin;
		self.opt['range'][ch]['xmax'] = xmax;
		self.opt['range'][ch]['ymin'] = ymin;
		self.opt['range'][ch]['ymax'] = ymax;

		self.figure['axes']['home'][1].set_xlim([xmin,xmax])
		self.figure['axes']['home'][1].set_ylim([ymin,ymax])
		return

	def _change_xmap(self):
		if self.note_in['xmap'].get() == self.prev_xmap: return
		self.prev_xmap = self.note_in['xmap'].get()
		diag = self.opt['cpage']; ch = self.opt['c_ch'];
		self.opt['xmap'][diag][ch] = self.prev_xmap;
		times = copy.deepcopy(self.figure['legend']['time'])
		self._clear_home_canvas()

		if self.prev_xmap == 1: xmin = 1.7; xmax = 2.3
		elif self.prev_xmap == 1: xmin = 0.; xmax = 1.4
		else: xmin = 0.; xmax = 1.3

		self.note_in['e3'].delete(0,'end'); self.note_in['e3'].insert(10,'%2.1f'%xmin);
		self.note_in['e4'].delete(0,'end'); self.note_in['e4'].insert(10,'%2.1f'%xmax);

		if    self.note_in['xmap'].get() ==1: self.figure['axes']['home'][1].set_xlabel('R[m]')
		elif  self.note_in['xmap'].get() ==2: self.figure['axes']['home'][1].set_xlabel('$\\psi_N$')
		else:                             self.figure['axes']['home'][1].set_xlabel('$\\rho_N$')		

		for time in times: self._draw_time(int(time))
		return

	def _delete_time(self,time):

		self.figure['name']['home'].canvas.draw_idle()
		ind = self.figure['legend']['time'].index(time)
		self.figure['pegend']['time'][ind].remove()
		self.figure['pegend']['time'].remove(self.figure['pegend']['time'][ind])
		self.figure['legend']['time'].remove(time)

		ind = self.figure['legend']['time2'].index(time)
		self.figure['pegend']['time2'][ind].remove()
		self.figure['pegend']['time2'].remove(self.figure['pegend']['time2'][ind])
		self.figure['legend']['time2'].remove(time)
		self._set_range()
		return

	def _efit_page(self,ch='efit'):

		if self.opt['c_ch'] == ch: return

		for i in range(1,6):
			if self.efit_list['isefit'][i]: 
				self.note_in['r%i'%i].config(state='normal')
				self.note_in['r%i'%(i+5)].config(state='normal')
			else:
				self.note_in['r%i'%i].config(state='disabled')
				self.note_in['r%i'%(i+5)].config(state='disabled')

		self._sync_list(True)		
		self.opt['cpage'] = 'efit'; self.opt['c_ch'] = ch; 
		self._clear_home_canvas(True)

		return

	def _ces_page(self,ch='ti'):

		if self.opt['c_ch'] == ch: return
		self._sync_list(True)
		self.note_in['xmap'].set(self.opt['xmap']['ces'][ch]);
		self._clear_home_canvas()
		self.opt['cpage'] = 'ces'; self.opt['c_ch'] = ch;
		self.prev_xmap = self.opt['xmap']['ces'][ch]
		self.note_in['e1'].delete(0,'end'); self.note_in['e1'].insert(10,'%i'%self.opt['TCES']);
		self.note_in['e2'].delete(0,'end'); self.note_in['e2'].insert(10,'%i'%self.opt['ACES']);
		self.note_in['l1'].delete(0,'end'); self.note_in['l2'].delete(0,'end');
		self.note_in['e3'].delete(0,'end'); self.note_in['e3'].insert(10,'%2.1f'%self.opt['range'][ch]['xmin']);
		self.note_in['e4'].delete(0,'end'); self.note_in['e4'].insert(10,'%2.1f'%self.opt['range'][ch]['xmax']);
		self.note_in['e5'].delete(0,'end'); self.note_in['e5'].insert(10,'%2.1f'%self.opt['range'][ch]['ymin']);
		self.note_in['e6'].delete(0,'end'); self.note_in['e6'].insert(10,'%2.1f'%self.opt['range'][ch]['ymax']);				
		self._sync_list(False)
		self._set_time()
		self._get_list('ces')
		tlist = list(self.note_in['l2'].get(0,'end'))
		for item in tlist: self._draw_time(int(item))
		return

	def _ts_page(self,ch='te'):

		if self.opt['c_ch'] == ch: return
		self._sync_list(True)
		self.note_in['xmap'].set(self.opt['xmap']['ts'][ch]);
		self._clear_home_canvas()
		self.opt['cpage'] = 'ts'; self.opt['c_ch'] = ch; 
		self.prev_xmap = self.opt['xmap']['ts'][ch]
		self.note_in['e1'].delete(0,'end'); self.note_in['e1'].insert(10,'%i'%self.opt['TTS']);
		self.note_in['e2'].delete(0,'end'); self.note_in['e2'].insert(10,'%i'%self.opt['ATS']);
		self.note_in['l1'].delete(0,'end'); self.note_in['l2'].delete(0,'end');
		self.note_in['e3'].delete(0,'end'); self.note_in['e3'].insert(10,'%2.1f'%self.opt['range'][ch]['xmin']);
		self.note_in['e4'].delete(0,'end'); self.note_in['e4'].insert(10,'%2.1f'%self.opt['range'][ch]['xmax']);
		self.note_in['e5'].delete(0,'end'); self.note_in['e5'].insert(10,'%2.1f'%self.opt['range'][ch]['ymin']);
		self.note_in['e6'].delete(0,'end'); self.note_in['e6'].insert(10,'%2.1f'%self.opt['range'][ch]['ymax']);				
		self._sync_list(False)
		self._set_time()
		self._get_list('ts')
		tlist = list(self.note_in['l2'].get(0,'end'))
		for item in tlist: self._draw_time(int(item))
		return

	def _generate(self):

		diag = self.opt['cpage']; ch = self.opt['c_ch']
		slist = list(self.note_in['l2'].get(0,'end'))
		if not len(slist) > 0: return
		if diag=='ces':
			nfile = os.getcwd()+'/'+'%s_%i_%i.dat'%(ch,self.shotn,self.time)
			f = open(nfile,'w')
			if ch=='ti': f.write('R[m]     Z[m]     TI[eV]   Error[eV]  \n')
			else:        f.write('R[m]     Z[m]     VT[km/s] Error[km/s]\n')

			xx = self.rr['ces'][1]; 
			tt = self.ces[ch]['val'][1][0]; yy = self.ces[ch]['val']; yerr = self.ces[ch]['err']
			for time in slist:
				tind = np.argmin(abs(tt-float(time)/1.e3))
				for i in range(1,self.ces['nch']+1):
					f.write('%9.7f  %9.7f  %9.4f  %9.4f\n'%(xx[i-1],0.,np.nan_to_num(yy[i][1][tind]),np.nan_to_num(yerr[i][1][tind])))
			f.close()

		if diag=='ts':
			nfile1 = os.getcwd()+'/'+'%s_%i_%i.dat'%(ch,self.shotn,self.time)
			nfile2 = os.getcwd()+'/'+'%s_%i_%i_edge.dat'%(ch,self.shotn,self.time)
			f1 = open(nfile1,'w'); f2 = open(nfile2,'w')
			if ch=='te': 
				f1.write('R[m]     Z[m]     TE[eV]   Error[eV]  \n')
				f2.write('R[m]     Z[m]     TE[eV]   Error[eV]  \n')
				factor = 1.
			else:        
				f1.write('R[m]     Z[m]     NE[18/m3] Error[18/m3]\n')
				f2.write('R[m]     Z[m]     NE[18/m3] Error[18/m3]\n')
				factor = 1.e-18;

			tt = self.ts[ch]['core'][1][0];
			xx1 = self.rr['ts']['core'][1]; yy1 = self.ts[ch]['core']; yerr1 = self.ts[ch]['core_err'];
			xx2 = self.rr['ts']['edge'][1]; yy2 = self.ts[ch]['edge']; yerr2 = self.ts[ch]['edge_err'];
			 
			for time in slist:
				tind = np.argmin(abs(tt-float(time)/1.e3))
				for i in range(1,self.ts['core']['nch']+1):
					f1.write('%9.7f  %9.7f  %9.4f  %9.4f\n'%(xx1[i-1],0.,np.nan_to_num(yy1[i][1][tind])*factor,np.nan_to_num(yerr1[i][1][tind])*factor))
				for i in range(1,self.ts['edge']['nch']+1):
					f2.write('%9.7f  %9.7f  %9.4f  %9.4f\n'%(xx2[i-1],0.,np.nan_to_num(yy2[i][1][tind])*factor,np.nan_to_num(yerr2[i][1][tind])*factor))				
			f1.close(); f2.close();

		self.opt['ready'][diag][ch] = True
		self._gui_preset()

		return

	def _gui_preset(self):

		if self.ces['nch'] < 3: 
			self.note_in['b1'].config(state='disabled');
			self.note_in['b2'].config(state='disabled');
		if self.ts['core']['nch'] < 3: 
			self.note_in['b3'].config(state='disabled');
			self.note_in['b4'].config(state='disabled');

		if self.opt['ready']['ces']['ti']: self.note_in['b1'].config(bg='green');
		else: self.note_in['b1'].config(bg='red');
		if self.opt['ready']['ces']['vt']: self.note_in['b2'].config(bg='green');
		else: self.note_in['b2'].config(bg='red');
		if self.opt['ready']['ts']['te']: self.note_in['b3'].config(bg='green');
		else: self.note_in['b3'].config(bg='red');
		if self.opt['ready']['ts']['ne']: self.note_in['b4'].config(bg='green');
		else: self.note_in['b4'].config(bg='red');						
		return

	def _reset_gui(self):

		diag = self.opt['cpage']; ch = self.opt['c_ch'];
		if not self.opt['ready'][diag][ch]: return
		self.opt['ready'][diag][ch] = False
		self._gui_preset()
		return

	def _get_value(self,tt,yy,time,factor=1.):

		if len(tt)>0: 
			ind = np.argmin(abs(tt-time/1.e3))
			val = yy[ind]*factor;
			if val>10.: return val;
			else: return ''
		else: return ''

	def _done(self):
		# isdata variable is included to check available data for given time slices
		if (not self.opt['ready']['ces']['ti'] and self.opt['isdata']['ces']): return
		if (not self.opt['ready']['ces']['vt'] and self.opt['isdata']['ces']): return
		if (not self.opt['ready']['ts']['te']  and self.opt['isdata']['ts']):  return
		if (not self.opt['ready']['ts']['ne']  and self.opt['isdata']['ts']):  return

		currdir = os.getcwd()
		if not self.opt['isdata']['ts']:
			self.te_file = ''
			self.ne_file = ''
			self.te_file_edge = ''
			self.ne_file_edge = ''
		else:
			self.te_file = currdir+'/te_%i_%i.dat'%(self.shotn,self.time)
			self.ne_file = currdir+'/ne_%i_%i.dat'%(self.shotn,self.time)
			self.te_file_edge = currdir+'/te_%i_%i_edge.dat'%(self.shotn,self.time)
			self.ne_file_edge = currdir+'/ne_%i_%i_edge.dat'%(self.shotn,self.time)

		if not self.opt['isdata']['ces']:
			self.ti_file = ''
			self.vt_file = ''
		else:
			self.ti_file = currdir+'/ti_%i_%i.dat'%(self.shotn,self.time)
			self.vt_file = currdir+'/vt_%i_%i.dat'%(self.shotn,self.time)		

		f = open('result.dat','w')
		f.write('g_file         %s\n'%(self.gkfiles['g'][self.efit_no_g.get()]))
		f.write('k_file         %s\n'%(self.gkfiles['k'][self.efit_no_k.get()]))
		f.write('te_edge_file   %s\n'%(self.te_file_edge))
		f.write('ne_edge_file   %s\n'%(self.ne_file_edge))
		f.write('te_file        %s\n'%(self.te_file))
		f.write('ne_file        %s\n'%(self.ne_file))
		f.write('ti_file        %s\n'%(self.ti_file))
		f.write('vt_file        %s\n'%(self.vt_file))

		wdia = self._get_value(self.mds['wdia'][0],self.mds['wdia'][1],self.time,1.e3)
		pnb1 = self._get_value(self.mds['pnb1'][0],self.mds['pnb1'][1],self.time,1.e6)
		pnb2 = self._get_value(self.mds['pnb2'][0],self.mds['pnb2'][1],self.time,1.e6)
		pnb3 = self._get_value(self.mds['pnb3'][0],self.mds['pnb3'][1],self.time,1.e6)
		pnb4 = self._get_value(self.mds['pnb4'][0],self.mds['pnb4'][1],self.time,1.e4)
		pnb5 = self._get_value(self.mds['pnb5'][0],self.mds['pnb5'][1],self.time,1.e4)
		pnb6 = self._get_value(self.mds['pnb6'][0],self.mds['pnb6'][1],self.time,1.e4)

		vg1 = self._get_value(self.mds['pnb1'][0],self.mds['pnb1'][2],self.time)
		vg2 = self._get_value(self.mds['pnb2'][0],self.mds['pnb2'][2],self.time)
		vg3 = self._get_value(self.mds['pnb3'][0],self.mds['pnb3'][2],self.time)
		vg4 = self._get_value(self.mds['pnb4'][0],self.mds['pnb4'][2],self.time)
		vg5 = self._get_value(self.mds['pnb5'][0],self.mds['pnb5'][2],self.time)
		vg6 = self._get_value(self.mds['pnb6'][0],self.mds['pnb6'][2],self.time)

		f.write('wdia    %s\n'%wdia)
		f.write('NBI1AP  %s\n'%pnb1)
		f.write('NBI1AE  %s\n'%vg1)
		f.write('NBI1BP  %s\n'%pnb2)
		f.write('NBI1BE  %s\n'%vg2)
		f.write('NBI1CP  %s\n'%pnb3)
		f.write('NBI1cE  %s\n'%vg3)
		f.write('NBI2AP  %s\n'%pnb4)
		f.write('NBI2AE  %s\n'%vg4)
		f.write('NBI2BP  %s\n'%pnb5)
		f.write('NBI2BE  %s\n'%vg5)
		f.write('NBI2CP  %s\n'%pnb6)
		f.write('NBI2CE  %s\n'%vg6)
		f.write('ECHPW   %s\n'%'')
		f.write('SHOTN   %i\n'%self.shotn)
		f.write('TIME    %i\n'%self.time)
		f.close()
		self._save_opt()
		print('>>> Finished!')
		exit()
		return

	def _save_opt(self):

		self._sync_list(True)
		self.opt['g_no'] = self.efit_no_g.get()
		self.opt['k_no'] = self.efit_no_k.get()

		f = open('mds.save','wb')
		pickle.dump(self.opt,f)
		f.close()
		return

	def _load_opt(self):

		if not os.path.isfile('mds.save'): return
		print('>>> Load previous opt.')
		f = open('mds.save','rb')
		temp = pickle.load(f)
		f.close()

		for key in temp.keys():
			self.opt[key] = copy.deepcopy(temp[key])

		self.efit_no_g.set(self.opt['g_no'])
		self.efit_no_k.set(self.opt['k_no'])

		self.opt['cpage'] = ''
		self.opt['c_ch']  = ''

		return

	def _make_canvas_frame(self,ptype=1):

		if ptype == 1:
			gs = gridspec.GridSpec(7,1)
			self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(gs[0:2]))
			self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(gs[2:7]))
		elif ptype == 2:
			gs = gridspec.GridSpec(2,3)
			self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(gs[:,0]))
			self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(gs[0,1]))
			self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(gs[1,1]))
			self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(gs[0,2]))
			self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(gs[1,2]))
		return		

	def _clear_home_canvas(self,efit_page=False):

		if efit_page:
			self.figure['name']['home'].canvas.draw_idle()
			self.figure['name']['home'].clf()
			self.figure['pegend'] = dict(); self.figure['legend'] = dict();
			self.figure['axes']['home'] = [];
			self._make_canvas_frame(ptype=2)			
			self._draw_efit()
			return

		for i in range(1,11): self.note_in['r%i'%i].config(state='disabled')

		if (self.opt['cpage']=='' or self.opt['cpage']=='tci' or self.opt['cpage']=='efit'):
			self.figure['name']['home'].canvas.draw_idle()
			self.figure['name']['home'].clf()
			self.figure['pegend'] = dict(); self.figure['legend'] = dict();
			self.figure['axes']['home'] = [];
			self._make_canvas_frame()
			self._draw_mds()
		else:
			self.figure['name']['home'].canvas.draw_idle()
			self.figure['axes']['home'][1].cla()
			for item in self.figure['pegend']['time']: item.remove()
			for item in self.figure['pegend']['time2']: item.remove()
			for item in self.figure['pegend']['stime']: item.remove()
			self.figure['pegend']['time'] = []; self.figure['legend']['time'] = [];
			self.figure['pegend']['time2'] = []; self.figure['legend']['time2'] = [];
			self.figure['pegend']['stime'] = []; self.figure['legend']['stime'] = [];
			self.figure['axes']['home'][0].set_xlabel('time [s]')

			if    self.note_in['xmap'].get() ==1: self.figure['axes']['home'][1].set_xlabel('R[m]')
			elif  self.note_in['xmap'].get() ==2: self.figure['axes']['home'][1].set_xlabel('$\\psi_N$')
			else:                             self.figure['axes']['home'][1].set_xlabel('$\\rho_N$')			
		return

	def __init__(self,shotn,time,aces,ats,nogui):
		self.shotn = shotn
		self.time  = time
		self.aces  = aces
		self.ats   = ats
		self.nogui = nogui
		return

if __name__ == "__main__":

	import gui_mds as gmds

	now = time.gmtime(time.time())
	years = now.tm_year, now.tm_mon, now.tm_mday
	hours = now.tm_hour, now.tm_min, now.tm_sec
	line = '#-- %i/%02i/%02i -- %02ih %02im %02is --'%(years[0],years[1],years[2],hours[0],hours[1],hours[2])
	os.system("echo '%s' >> log.mds"%line)

	try: shotn = int(sys.argv[1]); time = int(sys.argv[2])
	except: exit() 
	try: aces = int(sys.argv[3]); ats = int(sys.argv[4])
	except: aces = 150; ats = 150;
	try:    nogui = int(sys.argv[5])
	except: nogui = 0

	print(' -------------------------------------------------------------')
	print('||            KSTAR MDS post-processing tool Ver %s         ||'%version['mds'])
	print('||                 Profile & Channel comparison              ||')
	print('||            Developed by S.K.Kim (PU) & PLARE (SNU)        ||')
	print('--------------------------------------------------------------')


	gmds = gmds.kstar_diagnostic_tool(shotn,time,aces,ats,nogui)

	gmds.root = tk.Tk()
	gmds.root.title('KSTAR MDS tool #%i'%shotn)
	gmds._mds_tool()

