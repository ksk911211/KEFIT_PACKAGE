#!/usr/local/anaconda3/bin/python3
import os, sys
import numpy as np
import tkinter as tk
import eqdsk
from scipy.interpolate import interp1d
from shutil import move, copyfile, copytree, rmtree
from tkinter.filedialog import askopenfilename,asksaveasfilename, askdirectory
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import matplotlib.pyplot as plt
from gefit_tool import call_mse_data, read_mse_data
from exec_dirs import gfit_dir,gped_dir,pedscan_dir,version,author,dena_dir,shotk,years
from get_efit import get_year
import pickle
import time

class gui_toolp:

	def gfit_pack(self):

		self.root.resizable(0,0)

		self.l1 = tk.Label(self.root, text="==== INPUT =====  PROFILE FITTING PKG  ===== STATUS ====",justify='center')
		self.l1.grid(row=0, column=0,columnspan=3)
		self.b1 = tk.Button(self.root, text="RAW PROFILE FIT", bg = "lightgray",command=lambda: self.button_func1(),height = 2,width = 20)
		self.b1.grid(row=1, column=1, columnspan=1)

		self.b1 = tk.Button(self.root, text="EPED PREDICTOR", bg = "lightgray",command=lambda: self.button_func2(),height = 2,width = 20)
		self.b1.grid(row=2, column=1, columnspan=1)
		self.b1 = tk.Button(self.root, text="PED SCANNER", bg = "lightgray",command=lambda: self.button_func3(),height = 2,width = 20)
		self.b1.grid(row=3, column=1, columnspan=1)
		self.b1 = tk.Button(self.root, text="EXPORT PROFILE", bg = "lightgray",command=lambda: self.button_func4(),height = 2,width = 20)
		self.b1.grid(row=4, column=1, columnspan=1)
		self.b1 = tk.Button(self.root, text="DONE", bg = "lightgray",command=lambda: self.button_func5(),height = 1,width = 5)
		self.b1.grid(row=5, column=1, columnspan=1)	

		self.l100 = tk.Label(self.root, text="N/A",fg='red',anchor='w',width=8,justify='left')
		self.l100.grid(row=1, column=2,columnspan=1)	

		self.l101 = tk.Label(self.root, text="N/A",fg='red',anchor='w',width=8,justify='left')
		self.l101.grid(row=2, column=0,columnspan=1,sticky='e')

		self.l102 = tk.Label(self.root, text="N/A",fg='red',anchor='w',width=8,justify='left')
		self.l102.grid(row=2, column=2,columnspan=1)

		self.l103 = tk.Label(self.root, text="N/A",fg='red',anchor='w',width=8,justify='left')
		self.l103.grid(row=3, column=0,columnspan=1,sticky='e')	

		self.l104 = tk.Label(self.root, text="N/A",fg='red',anchor='w',width=8,justify='left')
		self.l104.grid(row=3, column=2,columnspan=1)

		self.l105 = tk.Label(self.root, text="N/A",fg='red',anchor='w',width=8,justify='left')
		self.l105.grid(row=4, column=2,columnspan=1)					

		self.l106 = tk.Label(self.root, text="Ver %s"%version['gfitp'],anchor='e',width=8,justify='left')
		self.l106.grid(row=5,column=2,sticky='se')


		ok,year = get_year(self.shot)
		if ok==0: exit()
		#self.b1 = tk.Button(self.root, text="DENA", bg = "lightgray",command=lambda: self.button_func6(),height = 1,width = 5)
		#self.b1.grid(row=5, column=0, columnspan=1)
		#if float(year) < 2020: self.b1.configure(state='disabled')

		self.update_status()

		self.root.mainloop()
		return

	def initialise_vars(self):

		self.t1_close = True

		self.eq_file = ''
		self.te_file = ''
		self.te_edge_file= ''
		self.ne_file = ''
		self.ne_edge_file= ''
		self.ti_file = ''
		self.vt_file = ''
		self.refl_file = None
		self.ece_file = None
		self.bnd_psi = 0.995
		self.lineden = 0.0

		self.zeff = 2.0
		self.amain = 2.0
		self.aimp = 12.0
		self.zimp = 6.0

		self.israw_fit = False
		self.iseped_prof = False
		self.iseped_tiprof = False

		self.iseped_run = False
		self.iseped_runprof = False

		self.ispedscan_run = False
		self.ispedscan_runprof = False

		self.bs_model = 'csauter'
		self.bsmulti = 1.0

		self.isvt = False
		self.xtype = 1

		self.xi_te0 = 0.
		self.xi_ne0 = 0.
		self.xi_ti0 = 0.
		self.xi_vt0 = 0.

		self.xi_te = 0.
		self.xi_ne = 0.
		self.xi_ti = 0.
		self.xi_vt = 0.		

		self.shot = 0
		self.time = 0
		self.mds_dir = './MDS'
		return

	def detect_close(self,win_num):

		if (win_num == 1):
			self.t1.destroy()
			self.t1_close = True			

		return			

	def write_history(self):

		f = open('gfitp_history.dat','w')

		if self.israw_fit:	f.write('raw_fit\n')
		if self.iseped_prof:	f.write('epedprof\n')
		if self.iseped_tiprof:	f.write('epedtiprof\n')

		if self.iseped_run:	f.write('eped_run\n')
		if self.iseped_runprof:	f.write('eped_runprof\n')
		if self.ispedscan_run:	f.write('pedscan_run\n')
		if self.ispedscan_runprof:	f.write('pedscan_runprof\n')

		f.write('zeff %f\n'%self.zeff)
		f.write('zimp %f\n'%self.zimp)
		f.write('amain %f\n'%self.amain)
		f.write('aimp %f\n'%self.aimp)
		if not self.eq_file == '':		f.write('eq_file %s\n'%self.eq_file)
		if not self.te_file == '':		f.write('te_file %s\n'%self.te_file)
		if not self.te_edge_file == '': f.write('te_edge_file %s\n'%self.te_edge_file)	
		if not self.ne_file == '':		f.write('ne_file %s\n'%self.ne_file)
		if not self.ne_edge_file == '': f.write('ne_edge_file %s\n'%self.ne_edge_file)
		if not self.ti_file == '':		f.write('ti_file %s\n'%self.ti_file)
		if not self.vt_file == '':		f.write('vt_file %s\n'%self.vt_file)		
		if not self.refl_file == None:	f.write('refl_file %s\n'%self.refl_file)
		if not self.ece_file  == None:	f.write('ece_file %s\n'%self.ece_file)				
		f.write('run_name EPED\n')
		f.write('run_name_ped ped_scan\n')
		f.write('bnd_psi %f\n'%self.bnd_psi)
		f.write('lineden %f\n'%self.lineden)
		f.write('bs_model %s\n'%self.bs_model)
		f.write('bsmulti %f\n'%self.bsmulti)
		f.write('shot %i\n'%self.shot)
		f.write('time %i\n'%self.time)		
		f.write('kin_file PROFILES/chease_eped_mod\n')
		f.write('ti_file2 PROFILES/chease_kinprof_fit\n')
		f.write('mds_dir %s\n'%self.mds_dir)

		f.close()

		return

	def write_history_file(self,dirs):

		if dirs[0] == '/': return dirs
		elif dirs[:2] == '..': return dirs[3:]
		else: return 'FITPROF/%s'%dirs

	def read_history(self):
		try:	f=open('gfitp_history.dat','r')
		except:	return

		while True:
			line = f.readline()
			if not line: break
			line = line.split()
			if line[0].lower() == 'raw_fit':		self.israw_fit = True
			if line[0].lower() == 'epedprof':		self.iseped_prof = True
			if line[0].lower() == 'epedtiprof':		self.iseped_tiprof = True

			if line[0].lower() == 'eped_run':		self.iseped_run = True
			if line[0].lower() == 'eped_runprof':	self.iseped_runprof = True

			if line[0].lower() == 'pedscan_run':	self.ispedscan_run = True
			if line[0].lower() == 'pedscan_runprof':self.ispedscan_runprof = True

			if line[0].lower() == 'zeff':			self.zeff = float(line[1])
			if line[0].lower() == 'amain':			self.amain = float(line[1])
			if line[0].lower() == 'zimp':			self.zimp = float(line[1])
			if line[0].lower() == 'aimp':			self.aimp = float(line[1])
			if line[0].lower() == 'eq_file':		self.eq_file = line[1]
			if line[0].lower() == 'te_file':		self.te_file = line[1]
			if line[0].lower() == 'ne_file':		self.ne_file = line[1]
			if line[0].lower() == 'te_edge_file':   self.te_edge_file = line[1]
			if line[0].lower() == 'ne_edge_file':   self.ne_edge_file = line[1]
			if line[0].lower() == 'ti_file':		self.ti_file = line[1]
			if line[0].lower() == 'vt_file':		self.vt_file = line[1]	
			if line[0].lower() == 'refl_file':		self.refl_file = line[1]
			if line[0].lower() == 'ece_file':		self.ece_file = line[1]				
			if line[0].lower() == 'bnd_psi':		self.bnd_psi = float(line[1])
			if line[0].lower() == 'lineden':		self.lineden = float(line[1])
			if line[0].lower() == 'bs_model':		self.bs_model = line[1].lower()
			if line[0].lower() == 'bsmulti':		self.bsmulti = float(line[1])
			if line[0].lower() == 'shot':			self.shot = int(float(line[1]))
			if line[0].lower() == 'time':			self.time = int(float(line[1]))
			if line[0].lower() == 'mds_dir':                self.mds_dir = line[1]

		f.close()

		return		

	def read_fitopt(self,filename):

		f = open(filename,'rb')
		fit_opt = pickle.load(f)
		f.close()

		self.iseped_prof = fit_opt['ped_scan_fit'] 
		self.iseped_tiprof = fit_opt['use_ti_eped']
		if fit_opt['use_density_scale']: self.lineden = fit_opt['target_density']

		self.eq_file = self.write_history_file(fit_opt['file']['gfile'])
		if not fit_opt['file']['te']['ece'] == None:  self.ece_file = self.write_history_file(fit_opt['file']['te']['ece'])
		if not fit_opt['file']['ne']['refl'] == None: self.refl_file = self.write_history_file(fit_opt['file']['ne']['refl'])

		if not fit_opt['file']['ne']['ts'] == None:   self.ne_file = self.write_history_file(fit_opt['file']['ne']['ts'])
		if not fit_opt['file']['ne']['tse'] == None:  self.ne_edge_file = self.write_history_file(fit_opt['file']['ne']['tse'])
		if not fit_opt['file']['te']['ts'] == None:   self.te_file = self.write_history_file(fit_opt['file']['te']['ts'])
		if not fit_opt['file']['te']['tse'] == None:  self.te_edge_file = self.write_history_file(fit_opt['file']['te']['tse'])		

		if not fit_opt['file']['ti']['ces'] == None:   self.ti_file = self.write_history_file(fit_opt['file']['ti']['ces'])
		if not fit_opt['file']['vt']['ces'] == None:   self.vt_file = self.write_history_file(fit_opt['file']['vt']['ces'])

		self.zeff   = fit_opt['zeff']
		self.aimp   = fit_opt['aimp']
		self.amain  = fit_opt['amain']
		self.zimp   = fit_opt['zimp']

		return

	def read_scanopt(self,filename):
		f = open(filename,'r')
		while True:
			line = f.readline()
			if not line: break
			line = line.split('=')
			
			if(len(line[0].split())==0):	linen = 'none'
			else:	linen = line[0].split()[0]

			if linen.lower() == 'use_neo':
				log = line[1].split()[0]
				if log.lower() == 'true':
					self.use_neo = True
				else:
					self.use_neo = False

			if linen.lower() == 'use_hager':
				log = line[1].split()[0]
				if log.lower() == 'true':
					self.use_hager = True
				else:
					self.use_hager = False

			if linen.lower() == 'use_chang':
				log = line[1].split()[0]
				if log.lower() == 'true':
					self.use_chang = True
				else:
					self.use_chang = False					

			if linen.lower() == 'bsmulti':
				self.bsmulti = float(line[1].split()[0])

		if self.use_neo:
			self.bs_model = 'neo'
		elif self.use_hager:
			self.bs_model = 'hager'
		elif self.use_chang:
			self.bs_model = 'csauter'
		else:
			self.bs_model = 'sauter'

		f.close()
		return			

	def read_epedopt(self,filename):

		f = open(filename,'r')
		while True:
			line = f.readline()
			if not line: break
			line = line.split('=')
			
			if(len(line[0].split())==0):	linen = 'none'
			else:	linen = line[0].split()[0]

			if linen.lower() == 'use_neo':
				log = line[1].split()[0]
				if log.lower() == 'true':
					self.use_neo = True
				else:
					self.use_neo = False

			if linen.lower() == 'use_hager':
				log = line[1].split()[0]
				if log.lower() == 'true':
					self.use_hager = True
				else:
					self.use_hager = False

			if linen.lower() == 'use_chang':
				log = line[1].split()[0]
				if log.lower() == 'true':
					self.use_chang = True
				else:
					self.use_chang = False					

			if linen.lower() == 'bsmulti':
				self.bsmulti = float(line[1].split()[0])

		if self.use_neo:
			self.bs_model = 'neo'
		elif self.use_hager:
			self.bs_model = 'hager'
		elif self.use_chang:
			self.bs_model = 'csauter'
		else:
			self.bs_model = 'sauter'

		f.close()
		return	

	def button_func1(self):
		try:	os.mkdir('FITPROF')
		except:	pass
		self.write_history()
		copyfile('gfitp_history.dat','FITPROF/gfitp_history.dat')
		os.chdir('FITPROF')
		if self.efit_mode: os.system(gfit_dir+' -efit')
		else: os.system(gfit_dir+' -fitp')
		os.chdir(self.currdir)
		isfile = os.path.isfile('FITPROF/fit_opt.save')
		iskfile = os.path.isfile('FITPROF/PROFILES/chease_eped_mod')

		try:	
			self.read_fitopt('FITPROF/fit_opt.save')
			#os.remove('FITPROF/fit_opt.save')
		except:
			print('>>> No saved gfit opt')
			return

		try:
			copyfile('FITPROF/PROFILES/chease_kinprof_extended','PROFILES/chease_kinprof_extended')
		except:
			print('>>> No extended kinectic profile')

		try:
			copyfile('FITPROF/PROFILES/VT_fit.dat','PROFILES/VT_fit.dat')
		except:
			print('>>> No rotation profile')			
				
		if (isfile):	self.israw_fit = True
		if (iskfile and self.iseped_prof):	self.iseped_prof = True
		if (self.iseped_prof and self.iseped_tiprof):	self.iseped_tiprof = True
		self.write_history()
		self.update_status()
		return

	def button_func1a(self):
		try:	os.mkdir(self.mds_dir)
		except:	pass
		dirs = os.getcwd() + '/'+self.mds_dir
		rdirs = dirs + '/result.dat'
		call_mse_data(dirs,None,100)
		dat = read_mse_data(rdirs)

		self.eq_file = 	    self.mds_dir+'/'+dat[0].split('/')[-1]
		self.te_edge_file = self.mds_dir+'/'+dat[2].split('/')[-1]		
		self.ne_edge_file = self.mds_dir+'/'+dat[3].split('/')[-1]
		self.te_file =	    self.mds_dir+'/'+dat[4].split('/')[-1]
		self.ne_file = 	    self.mds_dir+'/'+dat[5].split('/')[-1]
		self.ti_file = 	    self.mds_dir+'/'+dat[6].split('/')[-1]
		self.vt_file = 	    self.mds_dir+'/'+dat[7].split('/')[-1]
		self.write_history()

		print('>>> MDS data is loaded..')

		return

	def button_func2(self):
		self.ask_delete_file('eped')
		if self.iseped_prof:	
			copyfile('FITPROF/PROFILES/chease_eped_mod','PROFILES/chease_eped_mod')
			try:	copyfile('FITPROF/PROFILES/TI_fit.dat','PROFILES/TI_fit.dat')
			except:	pass
			try:	copyfile('FITPROF/PROFILES/NE_raw.dat','PROFILES/NE_raw.dat')
			except:	pass
			try:	copyfile('FITPROF/PROFILES/TE_raw.dat','PROFILES/TE_raw.dat')
			except:	pass
			try:	copyfile('FITPROF/PROFILES/TI_raw.dat','PROFILES/TI_raw.dat')
			except:	pass			
 
		else:	
			print('>>> No eped profile !')
			return

		if not self.israw_fit:	
			print('>>> Profile fitting is not ready !')
			return
		self.check_eped_run()
		self.write_history()
		os.system(gped_dir)
		self.check_eped_run()
		self.update_status()
		try:	self.read_epedopt('EPED/eped_opt')
		except:	pass
		self.write_history()
		return

	def button_func3(self):
		self.ask_delete_file('pedscan')
		if self.iseped_prof:	
			copyfile('FITPROF/PROFILES/chease_eped_mod','PROFILES/chease_eped_mod')
			copyfile('FITPROF/PROFILES/chease_kinprof_fit','PROFILES/chease_kinprof_fit')
			try:	copyfile('FITPROF/PROFILES/TI_fit.dat','PROFILES/TI_fit.dat')
			except:	pass
			try:	copyfile('FITPROF/PROFILES/NE_raw.dat','PROFILES/NE_raw.dat')
			except:	pass					
			try:	copyfile('FITPROF/PROFILES/TE_raw.dat','PROFILES/TE_raw.dat')
			except:	pass
			try:	copyfile('FITPROF/PROFILES/TI_raw.dat','PROFILES/TI_raw.dat')
			except:	pass			

		else:	
			print('>>> No eped profile !')
			return

		if not self.israw_fit:	
			print('>>> Profile fitting is not ready !')
			return	

		self.check_pedscan_run()
		self.write_history()
		os.system(pedscan_dir)
		self.check_pedscan_run()
		self.update_status()
		try:	self.read_scandopt('ped_scan/scan_opt')
		except:	pass		
		self.write_history()

		return	

	def button_func4(self):

		if not self.israw_fit:
			print('>>> No fitted profiles...')
			return
		else:
			copyfile('FITPROF/PROFILES/chease_kinprof_fit','PROFILES/chease_kinprof_fit')
		if self.iseped_prof:	
			copyfile('FITPROF/PROFILES/chease_eped_mod','PROFILES/chease_eped_mod')
			copyfile('FITPROF/PROFILES/chease_kinprof_fit','PROFILES/chease_kinprof_fit')			

		try:	copyfile('FITPROF/PROFILES/NE_raw_TS.dat','PROFILES/NE_raw_TS.dat')
		except:	pass
		try:    copyfile('FITPROF/PROFILES/NE_raw_TSE.dat','PROFILES/NE_raw_TSE.dat')
		except: pass
		try:	copyfile('FITPROF/PROFILES/NE_raw_REFL.dat','PROFILES/NE_raw_REFL.dat')
		except:	pass		
		try:	copyfile('FITPROF/PROFILES/TE_raw_TS.dat','PROFILES/TE_raw_TS.dat')
		except:	pass
		try:    copyfile('FITPROF/PROFILES/TE_raw_TSE.dat','PROFILES/TE_raw_TSE.dat')
		except: pass
		try:	copyfile('FITPROF/PROFILES/TE_raw_ECE.dat','PROFILES/TE_raw_ECE.dat')
		except:	pass		
		try:	copyfile('FITPROF/PROFILES/TI_raw_CES.dat','PROFILES/TI_raw_CES.dat')
		except:	pass
		try:	copyfile('FITPROF/PROFILES/VT_raw_CES.dat','PROFILES/VT_raw_CES.dat')
		except:	pass
		try:	copyfile('FITPROF/PROFILES/VT_fit.dat','PROFILES/VT_fit.dat')
		except:	pass

		self.CheckVar1.set(1)
		eq = eqdsk.eqdsk(self.eq_file)
		eq.read_eqdsk_file()

		self.Vf = interp1d(eq.avolp[:,0],eq.avolp[:,2],'cubic')
		self.wmhd = eq.wmhd;

		psin2 = np.linspace(0,2.0,301)
		rho2 = np.copy(psin2);	R2 =np.copy(psin2);  R3 =np.copy(psin2);
		num = len(eq.prhoR[:,0])
		drhodp = (eq.prhoR[-1,1]-eq.prhoR[num-2,1])/(eq.prhoR[-1,0]-eq.prhoR[num-2,0])
		dRdp = (eq.prhoR[-1,2]-eq.prhoR[num-2,2])/(eq.prhoR[-1,0]-eq.prhoR[num-2,0])
		dRdp2 = (eq.prhoR[-1,3]-eq.prhoR[num-2,3])/(eq.prhoR[-1,0]-eq.prhoR[num-2,0])

		rhof = interp1d(eq.prhoR[:,0],eq.prhoR[:,1],'cubic')
		Rf = interp1d(eq.prhoR[:,0],eq.prhoR[:,2],'cubic')
		Rf2 = interp1d(eq.prhoR[:,0],eq.prhoR[:,3],'cubic')

		for i in range(301):
			if psin2[i] < 1.0:
				rho2[i] = rhof(psin2[i])
				R2[i] = Rf(psin2[i])
				R3[i] = Rf2(psin2[i])
			else:
				rho2[i] = drhodp*(psin2[i]-1.) + eq.prhoR[-1,1]
				R2[i] = dRdp*(psin2[i]-1.) + eq.prhoR[-1,2]
				R3[i] = dRdp2*(psin2[i]-1.) + eq.prhoR[-1,3]

		self.psinf = interp1d(psin2,psin2)
		self.rhof = interp1d(psin2,rho2,'cubic')
		self.Rf = interp1d(psin2,R2,'cubic')
		self.Rf2 = interp1d(psin2,R3 ,'cubic')

		if not self.t1_close:
			print('>>> Spectrum Window is already opened...')
			return

		self.t1 = tk.Toplevel(self.root)
		self.t1.wm_title("Profile generator")
		self.t1_close = False		
		self.t1.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(1))


		self.t1.resizable(0,0)
		self.fig1, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2,2,figsize=(11,8))

		ax1.set_title('$T_e$ [keV]')
		ax2.set_title('$n_e$ [$10^{19}$/m3]')
		ax3.set_title('$T_i$ [keV]')
		ax4.set_title('$V_\phi$ [km/s]')
		if self.xtype == 1:
			ax1.set_xlabel('$\psi_n$ [a.u]')
			ax2.set_xlabel('$\psi_n$ [a.u]')
			ax3.set_xlabel('$\psi_n$ [a.u]')
			ax4.set_xlabel('$\psi_n$ [a.u]')
		elif self.xtype == 2:
			ax1.set_xlabel('$\\rho_t$ [a.u]')
			ax2.set_xlabel('$\\rho_t$ [a.u]')
			ax3.set_xlabel('$\\rho_t$ [a.u]')
			ax4.set_xlabel('$\\rho_t$ [a.u]')
		else:	
			ax1.set_xlabel('R [m]')
			ax2.set_xlabel('R [m]')
			ax3.set_xlabel('R [m]')
			ax4.set_xlabel('R [m]')

		ax1.set_ylabel('$T_e$ [keV]')
		ax2.set_ylabel('$n_e$ [$10^{19}$/m3]')
		ax3.set_ylabel('$T_i$ [keV]')
		ax4.set_ylabel('$V_\phi$ [km/s]')		

		self.canvas = FigureCanvasTkAgg(self.fig1,master=self.t1)
		self.plot_widget = self.canvas.get_tk_widget()
		self.plot_widget.grid(rowspan=30,row=1,column=10,columnspan=10)

		toolbar_frame = tk.Frame(self.t1)
		toolbar_frame.grid(column=10,row=0)
	
		self.toolbar = NavigationToolbar2Tk(self.canvas,toolbar_frame)
		self.fig1.canvas.draw_idle()

		self.l1 = tk.Label(self.t1, text="====== Profile  Generator ======",justify='center')
		self.l1.grid(row=0, column=0,columnspan=4)

		self.l1 = tk.Label(self.t1, text="======== Fit - Option ========",justify='center')
		self.l1.grid(row=1, column=0,columnspan=4)

		self.l2 = tk.Label(self.t1, text="RAW ",anchor='e')	
		self.l2.grid(row=2, column=1,columnspan=1)
		self.c1= tk.Checkbutton(self.t1,variable=self.CheckVar1,command=lambda: self.check_button_func1(1))
		self.c1.grid(row=2, column=2)
		
		count = 3

		if self.iseped_runprof:

			self.l2 = tk.Label(self.t1, text="EPED ",anchor='e')	
			self.l2.grid(row=count, column=1,columnspan=1)
			self.c1= tk.Checkbutton(self.t1,variable=self.CheckVar2,command=lambda: self.check_button_func1(2))
			self.c1.grid(row=count, column=2)
			count = count + 1

		if self.ispedscan_runprof:

			self.l2 = tk.Label(self.t1, text="PEDSCAN ",anchor='e')	
			self.l2.grid(row=count, column=1,columnspan=1)
			self.c1= tk.Checkbutton(self.t1,variable=self.CheckVar3,command=lambda: self.check_button_func1(3))
			self.c1.grid(row=count, column=2)

		self.l1 = tk.Label(self.t1, text="========= R - Option =========",justify='center')
		self.l1.grid(row=5, column=0,columnspan=4)

		self.l2 = tk.Label(self.t1, text="PSIN ",anchor='e')	
		self.l2.grid(row=6, column=1,columnspan=1)
		self.c1= tk.Checkbutton(self.t1,variable=self.CheckVar4,command=lambda: self.check_button_func2(1))
		self.c1.grid(row=6, column=2)

		self.l2 = tk.Label(self.t1, text="RHO ",anchor='e')	
		self.l2.grid(row=7, column=1,columnspan=1)
		self.c1= tk.Checkbutton(self.t1,variable=self.CheckVar5,command=lambda: self.check_button_func2(2))
		self.c1.grid(row=7, column=2)	

		self.l2 = tk.Label(self.t1, text="R [m] ",anchor='e')	
		self.l2.grid(row=8, column=1,columnspan=1)
		self.c1= tk.Checkbutton(self.t1,variable=self.CheckVar6,command=lambda: self.check_button_func2(3))
		self.c1.grid(row=8, column=2)	

		self.l2 = tk.Label(self.t1, text="R+ [m] ",anchor='e')	
		self.l2.grid(row=9, column=1,columnspan=1)
		self.c1= tk.Checkbutton(self.t1,variable=self.CheckVar7,command=lambda: self.check_button_func2(4))
		self.c1.grid(row=9, column=2)

		self.l1 = tk.Label(self.t1, text="======== 0-D Params ========",justify='center')
		self.l1.grid(row=10, column=0,columnspan=4)

		self.l2 = tk.Label(self.t1, text="NL [10(19)/m3] ",anchor='e')	
		self.l2.grid(row=11, column=0,columnspan=2,sticky='e')
		self.e1= tk.Entry(self.t1,width=6,justify='center')
		self.e1.grid(row=11, column=1,columnspan=2,sticky='e')
		self.e1.insert(0,'0')
		self.e1.configure(state='readonly')

		self.l2 = tk.Label(self.t1, text="WMHD [kJ] ",anchor='e')	
		self.l2.grid(row=12, column=0,columnspan=2,sticky='e')		
		self.e2 = tk.Entry(self.t1,width=6,justify='center')
		self.e2.grid(row=12, column=1,columnspan=2,sticky='e')
		self.e2.insert(0,'0')
		self.e2.configure(state='readonly')

		self.l2 = tk.Label(self.t1, text="WTHERM [kJ] ",anchor='e')	
		self.l2.grid(row=13, column=0,columnspan=2,sticky='e')
		self.e3= tk.Entry(self.t1,width=6,justify='center')
		self.e3.grid(row=13, column=1,columnspan=2,sticky='e')
		self.e3.insert(0,'0')
		self.e3.configure(state='readonly')

		self.l1 = tk.Label(self.t1, text="======== Square error ========",justify='center')
		self.l1.grid(row=14, column=0,columnspan=4)

		self.l2 = tk.Label(self.t1, text="RAW  ",anchor='e')	
		self.l2.grid(row=15, column=0,columnspan=2,sticky='e')
		self.l2 = tk.Label(self.t1, text="SELECT",anchor='e')	
		self.l2.grid(row=15, column=1,columnspan=2,sticky='e')		

		self.l2 = tk.Label(self.t1, text="TE ",anchor='e')	
		self.l2.grid(row=16, column=0,columnspan=1,sticky='e')
		self.e4= tk.Entry(self.t1,width=6,justify='center')
		self.e4.grid(row=16, column=0,columnspan=2,sticky='e')
		self.e4.insert(0,'0')
		self.e4.configure(state='readonly')
		self.e8= tk.Entry(self.t1,width=6,justify='center')
		self.e8.grid(row=16, column=1,columnspan=2,sticky='e')
		self.e8.insert(0,'0')
		self.e8.configure(state='readonly')		

		self.l2 = tk.Label(self.t1, text="NE ",anchor='e')	
		self.l2.grid(row=17, column=0,columnspan=1,sticky='e')		
		self.e5 = tk.Entry(self.t1,width=6,justify='center')
		self.e5.grid(row=17, column=0,columnspan=2,sticky='e')
		self.e5.insert(0,'0')
		self.e5.configure(state='readonly')
		self.e9 = tk.Entry(self.t1,width=6,justify='center')
		self.e9.grid(row=17, column=1,columnspan=2,sticky='e')
		self.e9.insert(0,'0')
		self.e9.configure(state='readonly')		

		self.l2 = tk.Label(self.t1, text="TI ",anchor='e')	
		self.l2.grid(row=18, column=0,columnspan=1,sticky='e')
		self.e6= tk.Entry(self.t1,width=6,justify='center')
		self.e6.grid(row=18, column=0,columnspan=2,sticky='e')
		self.e6.insert(0,'0')
		self.e6.configure(state='readonly')
		self.e10= tk.Entry(self.t1,width=6,justify='center')
		self.e10.grid(row=18, column=1,columnspan=2,sticky='e')
		self.e10.insert(0,'0')
		self.e10.configure(state='readonly')			

		self.l2 = tk.Label(self.t1, text="VT ",anchor='e')	
		self.l2.grid(row=19, column=0,columnspan=1,sticky='e')
		self.e7= tk.Entry(self.t1,width=6,justify='center')
		self.e7.grid(row=19, column=0,columnspan=2,sticky='e')
		self.e7.insert(0,'0')
		self.e7.configure(state='readonly')	
		self.e11= tk.Entry(self.t1,width=6,justify='center')
		self.e11.grid(row=19, column=1,columnspan=2,sticky='e')
		self.e11.insert(0,'0')
		self.e11.configure(state='readonly')							

		b1 = tk.Button(self.t1, text="SAVE", bg = "lightgray",command=lambda: self.button_func4a(self.xtype),height = 1,width = 4)
		b1.grid(row=21, column=1,columnspan=1,sticky='w')
		b2= tk.Button(self.t1, text="EXIT", bg = "lightgray",command=lambda: self.button_func4b(),height = 1,width = 4)
		b2.grid(row=21, column=2,columnspan=1,sticky='e')						
		

		self.check_button_func1(1)

		self.xi_te0 = self.xi_te
		self.xi_ne0 = self.xi_ne
		self.xi_ti0 = self.xi_ti
		self.xi_vt0 = self.xi_vt

		return		

	def read_profile(self,type='raw'):

		if os.path.isfile('PROFILES/chease_kinprof_extended'):	self.isvt = True

		if (self.isvt):
			f = open('PROFILES/chease_kinprof_extended','r')
			num = int(float(f.readline()))
			line = f.readline().split()
			self.psin2 = np.zeros(num);	self.te2 = np.zeros(num);  self.ne2 = np.zeros(num);  self.ti2 = np.zeros(num);	self.vt2 = np.zeros(num)
			for i in range(num):
				line = f.readline().split()
				self.psin2[i] = float(line[0])
				self.te2[i] = float(line[1])
				self.ne2[i] = float(line[2])
				self.ti2[i] = float(line[3])
				self.vt2[i] = float(line[4])
			f.close()

			vtf = interp1d(self.psin2,self.vt2,'cubic')

	
		if (type == 'eped' and self.iseped_runprof):	f = open('EPED/result/chease_kinprof_eped','r')
		elif (type == 'pedscan' and self.ispedscan_runprof):	f = open('ped_scan/result/chease_kinprof_pedscan','r')
		elif type == 'raw':	f = open('PROFILES/chease_kinprof_fit','r')
		num = int(float(f.readline()))
		line = f.readline().split()
		self.zeff = float(line[0]);	self.zimp = float(line[1]);	self.amain = float(line[2]); self.aimp = float(line[3]);
		self.psin = np.zeros(num);	self.te = np.zeros(num);  self.ne = np.zeros(num);  self.ti = np.zeros(num);
		for i in range(num):
			line = f.readline().split()
			self.psin[i] = float(line[0])
			self.te[i] = float(line[1])
			self.ne[i] = float(line[2])
			self.ti[i] = float(line[3])

		f.close()

		if self.isvt:	self.vt = vtf(self.psin)

		return
		
	def draw_plot(self,fig_ex,type,xtype):

		self.xi_te = 0;	self.xi_ne = 0;	self.xi_ti = 0;	self.xi_vt = 0;

		fig = fig_ex
		[ax1,ax2, ax3, ax4] = fig.axes
		ax1.cla()
		ax2.cla()
		ax3.cla()
		ax4.cla()

		if xtype == 1:	xf = self.psinf
		elif xtype == 2:	xf = self.rhof
		else:	xf = self.Rf


		ax1.set_title('$T_e$ [keV]')
		ax2.set_title('$n_e$ [$10^{19}$/m3]')
		ax3.set_title('$T_i$ [keV]')
		ax4.set_title('$V_\phi$ [km/s]')
		if xtype == 1:
			ax1.set_xlabel('$\psi_n$ [a.u]')
			ax2.set_xlabel('$\psi_n$ [a.u]')
			ax3.set_xlabel('$\psi_n$ [a.u]')
			ax4.set_xlabel('$\psi_n$ [a.u]')
			ax1.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
			ax2.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
			ax3.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
			ax4.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
		elif xtype == 2:
			ax1.set_xlabel('$\\rho_t$ [a.u]')
			ax2.set_xlabel('$\\rho_t$ [a.u]')
			ax3.set_xlabel('$\\rho_t$ [a.u]')
			ax4.set_xlabel('$\\rho_t$ [a.u]')
			ax1.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
			ax2.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
			ax3.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
			ax4.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
		else:	
			ax1.set_xlabel('R [m]')
			ax2.set_xlabel('R [m]')
			ax3.set_xlabel('R [m]')
			ax4.set_xlabel('R [m]')
			ax1.axvline(x=xf(1.0),color='goldenrod',linestyle='--',linewidth=1.0)
			ax2.axvline(x=xf(1.0),color='goldenrod',linestyle='--',linewidth=1.0)
			ax3.axvline(x=xf(1.0),color='goldenrod',linestyle='--',linewidth=1.0)
			ax4.axvline(x=xf(1.0),color='goldenrod',linestyle='--',linewidth=1.0)
			
		ax1.set_ylabel('$T_e$ [keV]')
		ax2.set_ylabel('$n_e$ [$10^{19}$/m3]')
		ax3.set_ylabel('$T_i$ [keV]')
		ax4.set_ylabel('$V_\phi$ [km/s]')	

		legend1 = [];	legend2 = [];	legend3 = [];	legend4 = [];
		plegend1 = []; plegend2 = []; plegend3 = [];	plegend4 = [];				

		self.read_profile('raw')
		if self.isvt:
			line1, = ax1.plot(xf(self.psin2),self.te2,'b')
			line2, = ax2.plot(xf(self.psin2),self.ne2,'b')
			line3, = ax3.plot(xf(self.psin2),self.ti2,'b')
		else:
			line1, = ax1.plot(xf(self.psin),self.te,'b')
			line2, = ax2.plot(xf(self.psin),self.ne,'b')
			line3, = ax3.plot(xf(self.psin),self.ti,'b')			
		max1 = max(self.te); max2 = max(self.ne); max3 = max(self.ti);	max4 = 1.
		plegend1.append(line1);	plegend2.append(line2);	plegend3.append(line3)
		legend1.append('Raw fit');	legend2.append('Raw fit');	legend3.append('Raw fit');	
		if self.isvt:	
			line4, = ax4.plot(xf(self.psin2),self.vt2,'b')
			max4 = max(self.vt)
			plegend4.append(line4)
			legend4.append('Raw fit')

		if xtype==4:
			if self.isvt:
				ax1.plot(self.Rf2(self.psin),self.te,'b')
				ax2.plot(self.Rf2(self.psin),self.ne,'b')
				ax3.plot(self.Rf2(self.psin),self.ti,'b')
			else:
				ax1.plot(self.Rf2(self.psin2),self.te2,'b')
				ax2.plot(self.Rf2(self.psin2),self.ne2,'b')
				ax3.plot(self.Rf2(self.psin2),self.ti2,'b')				
			if self.isvt:	ax4.plot(self.Rf2(self.psin2),self.vt2,'b')		
			ax1.axvline(x=self.Rf2(1.0),color='goldenrod',linestyle='--',linewidth=1.0)
			ax2.axvline(x=self.Rf2(1.0),color='goldenrod',linestyle='--',linewidth=1.0)
			ax3.axvline(x=self.Rf2(1.0),color='goldenrod',linestyle='--',linewidth=1.0)
			ax4.axvline(x=self.Rf2(1.0),color='goldenrod',linestyle='--',linewidth=1.0)	

		if type == 'eped':
			legend1.append('EPED fit');	legend2.append('EPED fit');	legend3.append('EPED fit');	
		elif type == 'pedscan':
			legend1.append('PEDSCAN fit');	legend2.append('PEDSCAN fit');	legend3.append('PEDSCAN fit');	

		if (type == 'eped' or type == 'pedscan'):
			self.read_profile(type)
			line1, = ax1.plot(xf(self.psin),self.te,'r')
			line2, = ax2.plot(xf(self.psin),self.ne,'r')
			line3, = ax3.plot(xf(self.psin),self.ti,'r')
			plegend1.append(line1);	plegend2.append(line2);	plegend3.append(line3)
			if self.isvt:	
				line4, = ax4.plot(xf(self.psin),self.vt,'r')
				plegend4.append(line4)
				if type == 'eped':	legend4.append('EPED fit')
				else:	legend4.append('PEDSCAN fit')

			if xtype==4:
				ax1.plot(self.Rf2(self.psin),self.te,'r')
				ax2.plot(self.Rf2(self.psin),self.ne,'r')
				ax3.plot(self.Rf2(self.psin),self.ti,'r')
				if self.isvt:	ax4.plot(self.Rf2(self.psin),self.vt,'r')

		try:
			f = open('PROFILES/TE_raw_TS.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,4))
			for i in range(num):
				line = f.readline().split()
				if xtype < 3:	dat[i,0] = xf(float(line[0]))
				else:	dat[i,0] = float(line[3])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
				dat[i,3] = float(line[0])
			f.close()
			line2 = ax1.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='lime',ecolor='lime',capthick=2)
			legend1.append('Raw data TS');	plegend1.append(line2)
			self.xi_te = self.calculate_xi2(self.psin,self.te,dat[:,3],dat[:,1],dat[:,2])
		except:
			pass

		try:
			f = open('PROFILES/TE_raw_TSE.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,4))
			for i in range(num):
				line = f.readline().split()
				if xtype < 3:	dat[i,0] = xf(float(line[0]))
				else:	dat[i,0] = float(line[3])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
				dat[i,3] = float(line[0])
			f.close()
			line2 = ax1.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='magenta',ecolor='magenta',capthick=2)
			legend1.append('Raw data TSE');	plegend1.append(line2)
			self.xi_te = self.calculate_xi2(self.psin,self.te,dat[:,3],dat[:,1],dat[:,2])
		except:
			pass
		
		try:
			f = open('PROFILES/TE_raw_ECE.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,4))
			for i in range(num):
				line = f.readline().split()
				if xtype < 3:	dat[i,0] = xf(float(line[0]))
				else:	dat[i,0] = float(line[3])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
				dat[i,3] = float(line[0])
			f.close()
			line2 = ax1.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='orange',ecolor='orange',capthick=2)
			legend1.append('Raw data ECE');	plegend1.append(line2)
			self.xi_te = self.calculate_xi2(self.psin,self.te,dat[:,3],dat[:,1],dat[:,2])
		except:
			pass			

		try:
			f = open('PROFILES/NE_raw_TS.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,4))
			for i in range(num):
				line = f.readline().split()
				if xtype < 3:	dat[i,0] = xf(float(line[0]))
				else:	dat[i,0] = float(line[3])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
				dat[i,3] = float(line[0])
			f.close()
			line1 = ax2.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='lime',ecolor='lime',capthick=2)
			legend2.append('Raw data TS');	plegend2.append(line1)
			self.xi_ne = self.calculate_xi2(self.psin,self.ne,dat[:,3],dat[:,1],dat[:,2])
		except:
			pass

		try:
			f = open('PROFILES/NE_raw_TSE.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,4))
			for i in range(num):
				line = f.readline().split()
				if xtype < 3:	dat[i,0] = xf(float(line[0]))
				else:	dat[i,0] = float(line[3])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
				dat[i,3] = float(line[0])
			f.close()
			line1 = ax2.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='magenta',ecolor='magenta',capthick=2)
			legend2.append('Raw data TSE');	plegend2.append(line1)
			self.xi_ne = self.calculate_xi2(self.psin,self.ne,dat[:,3],dat[:,1],dat[:,2])
		except:
			pass			

		try:
			f = open('PROFILES/NE_raw_REFL.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,4))
			for i in range(num):
				line = f.readline().split()
				if xtype < 3:	dat[i,0] = xf(float(line[0]))
				else:	dat[i,0] = float(line[3])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
				dat[i,3] = float(line[0])
			f.close()
			line1 = ax2.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='orange',ecolor='orange',capthick=2)
			legend2.append('Raw data REFL');	plegend2.append(line1)
			self.xi_ne = self.calculate_xi2(self.psin,self.ne,dat[:,3],dat[:,1],dat[:,2])
		except:
			pass

		try:
			f = open('PROFILES/TI_raw_CES.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,4))
			for i in range(num):
				line = f.readline().split()
				if xtype < 3:	dat[i,0] = xf(float(line[0]))
				else:	dat[i,0] = float(line[3])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
				dat[i,3] = float(line[0])
			f.close()
			line3 = ax3.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='lime',ecolor='lime',capthick=2)
			legend3.append('Raw data CES');	plegend3.append(line3)
			self.xi_ti = self.calculate_xi2(self.psin,self.ti,dat[:,3],dat[:,1],dat[:,2])
		except:
			pass

		try:
			f = open('PROFILES/VT_raw_CES.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,4))
			for i in range(num):
				line = f.readline().split()
				if xtype < 3:	dat[i,0] = xf(float(line[0]))
				else:	dat[i,0] = float(line[3])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
				dat[i,3] = float(line[0])
			f.close()
			line4 = ax4.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='lime',ecolor='lime',capthick=2)
			legend4.append('Raw data CES');	plegend4.append(line4)
			self.xi_vt = self.calculate_xi2(self.psin,self.vt,dat[:,3],dat[:,1],dat[:,2])
		except:
			pass				

		ax1.legend(plegend1,legend1)
		ax2.legend(plegend2,legend2)
		ax3.legend(plegend3,legend3)
		ax4.legend(plegend4,legend4)

		ax1.set_xlim((-0.05,1.15))
		ax2.set_xlim((-0.05,1.15))
		ax3.set_xlim((-0.05,1.15))
		ax4.set_xlim((-0.05,1.15))

		self.update_xi2()

		if (xtype >=3):
			if xtype == 4:	xmin = self.Rf2(1.15)
			else:	xmin = xf(0.0)-0.01
			xmax = xf(1.15)
			ax1.set_xlim((xmin,xmax))
			ax2.set_xlim((xmin,xmax))
			ax3.set_xlim((xmin,xmax))
			ax4.set_xlim((xmin,xmax))

		ax1.set_ylim((-0.1,1.4*max1))					
		ax2.set_ylim((-0.1,1.4*max2))					
		ax3.set_ylim((-0.1,1.4*max3))					
		ax4.set_ylim((-10,1.2*max4))						
		fig.tight_layout()
		

		return

	def button_func5(self):

		self.root.destroy()
		exit()

		return
	def button_func6(self):
		
		if os.path.isfile('run_state'): os.remove('run_state')
		os.system('%s %i r &'%(dena_dir,self.shot))
		time.sleep(5)
		return

	def check_button_func1(self,xtype):

		if xtype == 1:
			self.CheckVar1.set(1)
			if self.CheckVar2.get() == 1:	self.CheckVar2.set(0)
			if self.CheckVar3.get() == 1:	self.CheckVar3.set(0)
		elif xtype == 2:
			self.CheckVar2.set(1)
			if self.CheckVar1.get() == 1:	self.CheckVar1.set(0)
			if self.CheckVar3.get() == 1:	self.CheckVar3.set(0)						
		elif xtype == 3:
			self.CheckVar3.set(1)
			if self.CheckVar1.get() == 1:	self.CheckVar1.set(0)
			if self.CheckVar2.get() == 1:	self.CheckVar2.set(0)

		if (self.CheckVar1.get()+self.CheckVar2.get()+self.CheckVar3.get()) == 0:	self.CheckVar1.set(1)

		if self.CheckVar1.get() == 1:
			ftype = 'raw'
		elif self.CheckVar2.get() == 1:
			ftype = 'eped'
		elif self.CheckVar3.get() == 1:	
			ftype = 'pedscan'

		if self.CheckVar4.get() == 1:
			xxtype = 1
		elif self.CheckVar5.get() == 1:
			xxtype = 2
		elif self.CheckVar6.get() == 1:	
			xxtype = 3
		elif self.CheckVar7.get() == 1:	
			xxtype = 4			

		self.fig1.canvas.draw_idle()
		self.draw_plot(self.fig1,ftype,xxtype)
		self.post_process()

		self.xtype = xxtype

		return

	def check_button_func2(self,xtype):

		if xtype == 1:
			self.CheckVar4.set(1)
			if self.CheckVar5.get() == 1:	self.CheckVar5.set(0)
			if self.CheckVar6.get() == 1:	self.CheckVar6.set(0)
			if self.CheckVar7.get() == 1:	self.CheckVar7.set(0)
		elif xtype == 2:
			self.CheckVar5.set(1) 
			if self.CheckVar4.get() == 1:	self.CheckVar4.set(0)
			if self.CheckVar6.get() == 1:	self.CheckVar6.set(0)
			if self.CheckVar7.get() == 1:	self.CheckVar7.set(0)				
		elif xtype == 3:
			self.CheckVar6.set(1)
			if self.CheckVar4.get() == 1:	self.CheckVar4.set(0)
			if self.CheckVar5.get() == 1:	self.CheckVar5.set(0)
			if self.CheckVar7.get() == 1:	self.CheckVar7.set(0)
		elif xtype == 4:
			self.CheckVar7.set(1)
			if self.CheckVar4.get() == 1:	self.CheckVar4.set(0)
			if self.CheckVar5.get() == 1:	self.CheckVar5.set(0)
			if self.CheckVar6.get() == 1:	self.CheckVar6.set(0)

		if (self.CheckVar4.get()+self.CheckVar5.get()+self.CheckVar6.get()+self.CheckVar7.get()) == 0:	self.CheckVar4.set(1)

		if self.CheckVar1.get() == 1:
			ftype = 'raw'
		elif self.CheckVar2.get() == 1:
			ftype = 'eped'
		elif self.CheckVar3.get() == 1:	
			ftype = 'pedscan'

		if self.CheckVar4.get() == 1:
			xxtype = 1
		elif self.CheckVar5.get() == 1:
			xxtype = 2
		elif self.CheckVar6.get() == 1:	
			xxtype = 3
		elif self.CheckVar7.get() == 1:	
			xxtype = 4

		self.fig1.canvas.draw_idle()
		self.draw_plot(self.fig1,ftype,xxtype)

		self.post_process()
		self.xtype = xxtype

		return		

	def declare_var(self):

		self.CheckVar1 = tk.IntVar()
		self.CheckVar2 = tk.IntVar()
		self.CheckVar3 = tk.IntVar()
		self.CheckVar4 = tk.IntVar()
		self.CheckVar5 = tk.IntVar()
		self.CheckVar6 = tk.IntVar()		
		self.CheckVar7 = tk.IntVar()

		self.StrVar1 = tk.StringVar()

		self.CheckVar1.set(1)
		self.CheckVar2.set(0)
		self.CheckVar3.set(0)
		self.CheckVar4.set(1)
		self.CheckVar5.set(0)
		self.CheckVar6.set(0)

		return		

	def check_eped_run(self):

		if os.path.isfile('EPED/result.dat'):	self.iseped_run = True
		else:	self.iseped_run = False
		if os.path.isfile('EPED/result/chease_kinprof_eped'):	self.iseped_runprof = True
		else:	self.iseped_runprof = False
		return

	def check_pedscan_run(self):

		if os.path.isfile('ped_scan/result.dat'):	self.ispedscan_run = True
		else:	self.ispedscan_run = False
		if os.path.isfile('ped_scan/result/chease_kinprof_pedscan'):	self.ispedscan_runprof = True
		else:	self.ispedscan_runprof = False

		return		

	def update_status(self):

		if self.israw_fit:	self.l100.config(text='DONE',fg='blue')
		else:	self.l100.config(text='N/A',fg='red')

		if self.iseped_prof:	self.l101.config(text='READY',fg='green')
		else:	self.l101.config(text='N/A',fg='red')

		if self.iseped_prof:	self.l103.config(text='READY',fg='green')
		else:	self.l103.config(text='N/A',fg='red')		

		if self.iseped_runprof:	self.l102.config(text='DONE',fg='blue')
		elif self.iseped_run:	self.l102.config(text='READY',fg='green')
		else:	self.l102.config(text='N/A',fg='red')

		if self.ispedscan_runprof:	self.l104.config(text='DONE',fg='blue')
		elif self.ispedscan_run:	self.l104.config(text='READY',fg='green')
		else:	self.l104.config(text='N/A',fg='red')

		if self.israw_fit:	self.l105.config(text='READY',fg='green')
		else:	self.l105.config(text='N/A',fg='red')

		return

	def post_process(self):

		self.linen = 0.;	self.wkin = 0.
		pres = self.ne*(self.te + (1.0 - (self.zeff-1.)/self.zimp)*self.ti)
		pres = pres * 1.602*1.e-19*1.e19*1.e3
		vol = self.Vf(self.psin)
		nef = interp1d(self.psin,self.ne,'cubic')
		RR = np.zeros(201)
		NN = np.zeros(201)
		for i in range(101):
			psin = float(i)/100.
			RR[100-i] = self.Rf2(psin)
			RR[100+i] = self.Rf(psin)
			NN[100-i] = nef(psin)
			NN[100+i] = NN[100-i]

		self.linen = np.trapz(NN,x=RR) / 1.9*2.0
		self.wkin = np.trapz(pres,x=vol) * 1.5

		self.e1.configure(state='normal')
		self.e1.delete(0, 'end')
		self.e1.insert(0,str(round(self.linen,2)))
		self.e1.configure(state='readonly')

		self.e2.configure(state='normal')
		self.e2.delete(0, 'end')
		self.e2.insert(0,str(round(self.wmhd/1.e3,1)))
		self.e2.configure(state='readonly')

		self.e3.configure(state='normal')
		self.e3.delete(0, 'end')
		self.e3.insert(0,str(round(self.wkin/1.e3,1)))
		self.e3.configure(state='readonly')		

		return

	def ask_delete_file(self,flag='eped'):

		if flag.lower() == 'eped':
			if not (os.path.isdir('EPED')):	return
		if flag.lower() == 'pedscan':
			if not (os.path.isdir('ped_scan')):	return
		self.t2 = tk.Toplevel(self.root)
		self.t2.wm_title('DELETE FILE')
		self.t2.resizable(0,0)


		self.l2 = tk.Label(self.t2, text = '     Do you want to remove previous run ?    ',justify='center')
		self.l2.grid(row=0,column=0,columnspan=4,pady=15)

		self.b1 = tk.Button(self.t2, text="YES", bg = "lightgray",command=lambda: self.button_func1b(True,flag),height = 1,width = 5)
		self.b1.grid(row=1, column=1, columnspan=1)	
		self.b1 = tk.Button(self.t2, text="NO", bg = "lightgray",command=lambda: self.button_func1b(False,flag),height = 1,width = 5)
		self.b1.grid(row=1, column=2, columnspan=1)	

		self.t2.mainloop()

		return

	def button_func1b(self,delete=False,flag='eped'):

		if delete:
			if flag.lower() == 'eped':
				try:	copyfile('EPED/eped_opt','eped_opt')
				except:	print('>>> error1')
				try:	rmtree('EPED_OLD')
				except:	print('>>> error2')
				try:	move('EPED','EPED_OLD')
				except:	print('>>> error3')
			if flag.lower() == 'pedscan':
				try:	copyfile('ped_scan/scan_opt','scan_opt')
				except:	print('>>> error4')			
				try:	rmtree('ped_scan_old')
				except:	print('>>> error5')			
				try:	move('ped_scan','ped_scan_old')
				except:	print('>>> error6')
		self.t2.destroy()
		self.t2.quit()
		return

	def button_func4a(self,xtype):

		#input2 = asksaveasfilename()
		#if len(input2) == 0:
		#	return
		input2 = self.currdir + '/profile.save'
		self.write_kinprof(input2,xtype)
		self.button_func4b()
		return

	def button_func4b(self):

		self.t1_close = True
		self.t1.destroy()

		return

	def write_kinprof(self,filename,xtype):

		f = open(filename,'w')
		f2 = open('PROFILES/chease_kinprof.out','w')
		num = len(self.psin);	num2= num;
		f2.write('%i\n'%num2);
		if xtype == 4:	
			num2 =2*num - 1
		f.write('%i\n'%num2);

		if xtype == 1:	xf = self.psinf
		elif xtype == 2:	xf = self.rhof
		else:	xf	= self.Rf


		if xtype == 1: xname = 'PSI_NORM'
		elif xtype == 2: xname = 'RHO_NORM'
		else: xname = 'R [m]'
		if self.isvt:
			line ='%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\n'%(xname,'TE[keV]','NE[19/m3]','TI[keV]','NI[19/m3]','VT[km/s]','WT[krad/s]')
		else:
			line ='%9s\t%9s\t%9s\t%9s\t%9s\n'%(xname,'TE[keV]','TE[keV]','NE[19/m3]','TI[keV]','NI[19/m3]')

		f.write(line)
		if xtype == 4:
			if self.isvt:
				for i in range(num-1):
					ind = num-1-i
					xx = self.Rf2(self.psin[ind])
					scale = 1.0 - (self.zeff-1.)/self.zimp
					f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(xx,self.te[ind],self.ne[ind],self.ti[ind],scale*self.ne[ind],self.vt[ind]))
			else:
				for i in range(num-1):
					xx = self.Rf2(self.psin[num-1-i])
					scale = 1.0 - (self.zeff-1.)/self.zimp
					f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(xx,self.te[i],self.ne[i],self.ti[i],scale*self.ne[i]))

		if self.isvt:
			for i in range(num):
				xx = xf(self.psin[i])
				R = self.Rf(self.psin[i])
				scale = 1.0 - (self.zeff-1.)/self.zimp
				f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(xx,self.te[i],self.ne[i],self.ti[i],scale*self.ne[i],self.vt[i],self.vt[i]/R))
		else:
			for i in range(num):
				xx = xf(self.psin[i])
				scale = 1.0 - (self.zeff-1.)/self.zimp
				f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(xx,self.te[i],self.ne[i],self.ti[i],scale*self.ne[i]))

		f.close()

		f2.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.zeff,self.zimp,self.amain,self.aimp))
		for i in range(num):
			if not self.isvt: f2.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.psin[i],self.te[i],self.ne[i],self.ti[i],scale*self.ne[i]))
			else: f2.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.psin[i],self.te[i],self.ne[i],self.ti[i],scale*self.ne[i],self.vt[i]))

		f2.close()

		print('>>> Profile is saved to %s'%(filename))
		print('>>> Chease_kinprof is saved to PROFILES/chease_kinprof.out')

		if not self.isvt:	return

		f1 = open('PROFILES/TE_pre.dat','w')
		f2 = open('PROFILES/NE_pre.dat','w')
		f3 = open('PROFILES/TI_pre.dat','w')
		f4 = open('PROFILES/VT_pre.dat','w')
		line = '%9s\t%9s\t%9s\t%9s\n'%('R[m]','Z[m]','Value','Error')
		f1.write(line);f2.write(line);f3.write(line);f4.write(line);

		i = 0;
		while i < (num-1):
			ind = num-1-i
			xx = self.psin[ind]
			RR = self.Rf2(xx)
			if RR > 1.8:
				f1.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(RR,0.,self.te[ind]*1.e3,self.te[ind]*1.e2))
				f2.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(RR,0.,self.ne[ind]*1.e1,self.ne[ind]*1.e0))
				f3.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(RR,0.,self.ti[ind]*1.e3,self.ti[ind]*1.e2))
				f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(RR,0.,self.vt[ind],self.vt[ind]*0.1))
			if xx < 0.8:	i = i + 10
			elif (xx > 0.8 and xx < 0.9):	i = i + 8
			else:	i = i + 2

		i = 0;
		while i < (num):
			ind = i
			xx = self.psin[ind]
			RR = self.Rf(xx)
			f1.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(RR,0.,self.te[ind]*1.e3,self.te[ind]*1.e2))
			f2.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(RR,0.,self.ne[ind]*1.e1,self.ne[ind]*1.e0))
			f3.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(RR,0.,self.ti[ind]*1.e3,self.ti[ind]*1.e2))
			f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(RR,0.,self.vt[ind],self.vt[ind]*0.1))		
			if xx < 0.8:	i = i + 12
			elif (xx > 0.8 and xx < 0.9):	i = i + 8
			else:	i = i + 1	

		f1.close()
		f2.close()
		f3.close()
		f4.close()

		return

	def calculate_xi2(self,x,y,xx,yy,ss):

		len2 = len(xx);	sum1 = 0;	count = 0;
		prof = interp1d(x,y)

		for i in range(len2):
			if xx[i] < 1.0:
				count = count + 1
				yf = prof(xx[i])
				sum1 = sum1 + (yf-yy[i])**2/(ss[i]+1.e-2)**2

		sum1 = sum1 / count

		return sum1

	def update_xi2(self):

		self.e4.configure(state='normal')
		self.e4.delete(0, 'end')
		self.e4.insert(0,str(round(self.xi_te0,2)))
		self.e4.configure(state='readonly')

		self.e5.configure(state='normal')
		self.e5.delete(0, 'end')
		self.e5.insert(0,str(round(self.xi_ne0,2)))
		self.e5.configure(state='readonly')		

		self.e6.configure(state='normal')
		self.e6.delete(0, 'end')
		self.e6.insert(0,str(round(self.xi_ti0,2)))
		self.e6.configure(state='readonly')

		self.e7.configure(state='normal')
		self.e7.delete(0, 'end')
		self.e7.insert(0,str(round(self.xi_vt0,2)))
		self.e7.configure(state='readonly')

		self.e8.configure(state='normal')
		self.e8.delete(0, 'end')
		self.e8.insert(0,str(round(self.xi_te,2)))
		self.e8.configure(state='readonly')

		self.e9.configure(state='normal')
		self.e9.delete(0, 'end')
		self.e9.insert(0,str(round(self.xi_ne,2)))
		self.e9.configure(state='readonly')		

		self.e10.configure(state='normal')
		self.e10.delete(0, 'end')
		self.e10.insert(0,str(round(self.xi_ti,2)))
		self.e10.configure(state='readonly')

		self.e11.configure(state='normal')
		self.e11.delete(0, 'end')
		self.e11.insert(0,str(round(self.xi_vt,2)))
		self.e11.configure(state='readonly')			

		return		

	def __init__(self):
		try:	os.mkdir('PROFILES')
		except:	pass
		self.currdir = os.getcwd()
		self.initialise_vars()
		self.read_history()
		self.check_eped_run()
		self.check_pedscan_run()

		return

if __name__ == "__main__":

	import gfitp

	gfit = gfitp.gui_toolp()	
	gfit.isplare = True
	gfit.efit_mode = False
	try:
		if (sys.argv[1].lower() == '-efit'): 
			gfit.efit_mode = True
			print('>>> EFIT MODE')
	except: pass
	gfit.root = tk.Tk()
	gfit.declare_var()
	gfit.root.title('Integrated Profile Fitting PKG')
	gfit.gfit_pack()
