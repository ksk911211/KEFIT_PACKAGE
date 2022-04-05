#!/usr/local/anaconda3/bin/python3
import os, sys
import numpy as np
import tkinter as tk
import tkinter.font
from shutil import move, copyfile, copytree
from tkinter.filedialog import askopenfilename,asksaveasfilename, askdirectory
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import matplotlib.pyplot as plt
import eqdsk
import eped
from scipy.optimize import curve_fit
from read_namelist import read_namelist_str
from batch_run import *
from exec_dirs import gfit2_dir, eped_dir, node_force, node_default, node_init,version,author,comment

class gped:

	def declare_vars(self):
		#Vars
		for i in range(1,30):
			self.__dict__['CheckVar%d'%i] = tk.IntVar()

		for i in range(1,71):
			self.__dict__['StrVar%d'%i] = tk.StringVar()

		for i in range(1,10):
			self.__dict__['DoubleVar%d'%i] = tk.DoubleVar()

		for i in range(1,30):
			self.__dict__['CheckVar%d'%i].set(0)			

		return

	def draw_psirz(self,iseqdsk):

		if not iseqdsk:
			return
		self.fig.canvas.draw_idle()
		[ax1] = self.fig.axes
		ax1.cla()
		ax1.contourf(self.eq.RR,self.eq.ZZ,self.eq.psirz,50)
		ax1.plot(self.eq.rzbdy[:,0],self.eq.rzbdy[:,1],'r')
		ax1.plot(self.eq.rzlim[:,0],self.eq.rzlim[:,1],'black')
		ax1.axis('scaled')
		ax1.set_title('Equilibrium')
		ax1.set_xlabel('R [m]')
		ax1.set_ylabel('Z [m]')
		ax1.set_xlim([min(self.eq.RR[0,:]),max(self.eq.RR[0,:])])
		ax1.set_ylim([min(self.eq.ZZ[:,0]),max(self.eq.ZZ[:,0])])
		#self.fig.tight_layout()		
		return

	def update_psirz(self):

		self.draw_psirz(self.iseqdsk)

		if not self.iseqdsk:
			self.fig.canvas.draw_idle()
		[ax1] = self.fig.axes
		ax1.plot(self.rzbdy[:,0],self.rzbdy[:,1],'--',color='pink')
		ax1.axis('scaled')


		R_max = (max(self.rzbdy[:,0])+min(self.rzbdy[:,0]))*0.5
		Z_max = 0.
		R_min = (max(self.rzbdy[:,0])+min(self.rzbdy[:,0]))*0.5
		Z_min = 0.
		ind = np.zeros(4)

		for i in range(len(self.rzbdy[:,0])):

			if (self.rzbdy[i,0] > R_max):
				R_max = self.rzbdy[i,0]
				ind[0] = i
			if (self.rzbdy[i,0] < R_min):
				R_min = self.rzbdy[i,0]
				ind[1] = i
			if (self.rzbdy[i,1] > Z_max):
				Z_max = self.rzbdy[i,1]
				ind[2] = i
			if (self.rzbdy[i,1] < Z_min):
				Z_min = self.rzbdy[i,1]
				ind[3] = i

		r0 = (R_min + R_max) * 0.5
		a0 = (R_max - R_min) * 0.5
		z0 = (Z_min + Z_max) * 0.5
	
		tri_u  = (r0 - self.rzbdy[int(ind[2]),0])/a0
		tri_l  = (r0 - self.rzbdy[int(ind[3]),0])/a0
		elong  = (Z_max - Z_min) / a0 / 2.

		inform = 'tri_u = %4.3f\n'%tri_u
		inform = inform + 'tri_l = %4.3f\n'%tri_l
		inform = inform + 'elong = %4.3f\n'%elong

		if self.iseqdsk:
			ax1.text(self.eq.rmag,self.eq.zmag,inform,color='lime')
			self.rmag = self.eq.rmag
		else:
			ax1.text(r0,z0,inform,color='lime')		

		plt.draw()	

		return

	def change_checkvar(self):

		if self.CheckVar1.get() == 1:
			self.CheckVar3.set(1)
		else:
			self.CheckVar3.set(0)

		if self.tfc_init:
			self.CheckVar3.set(1)	

		return

	def change_checkvar2(self):

		if self.CheckVar3.get() == 1:
			self.tfc_init = True
		else:
			self.tfc_init = False

		return

	def detect_close(self,win_num):

		#if tk.tkMessageBox.askokcancel("Quit", "Do you really wish to quit?"):
		if (win_num == 3):
			os.chdir(self.currdir)
			self.t3.destroy()
			self.t3_close = True
		elif (win_num == 2):
			self.t2.destroy()
			self.t2_close = True
		elif (win_num == 1):
			self.t1.destroy()
			self.t1_close = True
		elif (win_num == 4):
			self.t4.destroy()
			self.t4_close = True	
		elif (win_num == 5):
			self.t5.destroy()
			self.t5_close = True					
		return		

	def scroll_func1(self):

		if (self.CheckVar2.get() == 1):
			self.use_param_shape = True
		else:
			self.use_param_shape = False

		self.target_bnd	= self.DoubleVar1.get()

		if self.use_param_shape:
			return

		self.use_ext_bnd = False

		if (self.iseqdsk):

			self.eq.target_psin = self.target_bnd
			self.eq.get_target_bnd()
		try:
			self.rzbdy = np.copy(self.eq.rzbdyt)
		except:
			print('>>> Have to load equilibrium first')
		self.update_psirz()

		return

	def renew_eped_coef(self):
		fit_file = 'PROFILES/chease_eped_mod'
		if os.path.isfile(fit_file):

			f = open(fit_file,'r')
			for i in range(3):
				line = f.readline().split()
				for j in range(8):
						self.eped_prof[i,j] = round(float(line[j]),3)
			f.close()

			self.alpt1 = self.eped_prof[1,6]
			self.alpt2 = self.eped_prof[1,7]
			self.alpn1 = self.eped_prof[0,6]
			self.alpn2 = self.eped_prof[0,7]
			self.StrVar9.set(self.trans_vars(self.alpt1,2))
			self.StrVar10.set(self.trans_vars(self.alpt2,2))
			self.StrVar11.set(self.trans_vars(self.alpn1,2))
			self.StrVar12.set(self.trans_vars(self.alpn2,2))

			self.at1 = self.eped_prof[1,2]
			self.an1 = self.eped_prof[0,2]
			self.StrVar13.set(self.trans_vars(self.at1,2))
			self.StrVar14.set(self.trans_vars(self.an1,2))		

			self.tsep = self.eped_prof[1,1]
			self.nsep = round(self.eped_prof[0,1]/(self.eped_prof[0,1]+self.eped_prof[0,0]*2*np.tanh(1)),3)
			self.neped = round(self.eped_prof[0,1] + 2.0 * np.tanh(1) * self.eped_prof[0,0],3)
			self.teped = round(self.eped_prof[1,1] + 2.0 * np.tanh(1) * self.eped_prof[1,0],3)
			self.tewidth = 0.5 * (self.eped_prof[0,4]+self.eped_prof[1,4])
			self.StrVar15.set(self.trans_vars(self.tsep,2))
			self.StrVar16.set(self.trans_vars(self.nsep,2))
			self.StrVar34.set(self.trans_vars(self.neped,2))

			self.e16.delete(0, 'end' )
			self.e16.insert(10,self.StrVar15.get())
			self.e14.delete(0, 'end' )
			self.e14.insert(10,self.StrVar13.get())
			self.e10.delete(0, 'end' )
			self.e10.insert(10,self.StrVar9.get())
			self.e11.delete(0, 'end' )
			self.e11.insert(10,self.StrVar10.get())					
			self.e17.delete(0, 'end' )
			self.e17.insert(10,self.StrVar16.get())
			self.e15.delete(0, 'end' )
			self.e15.insert(10,self.StrVar14.get())
			self.e12.delete(0, 'end' )
			self.e12.insert(10,self.StrVar11.get())		
			self.e13.delete(0, 'end' )
			self.e13.insert(10,self.StrVar12.get())
			self.e35.delete(0, 'end' )
			self.e35.insert(10,self.StrVar34.get())

			return
	
		return

	def make_kbm_coef(self):

		if (self.eped_prof[0,0] == 0):
			return

		if not self.iseqdsk:
			return	

		width = 0.5 * (self.eped_prof[0,4]+self.eped_prof[1,4])
		#neped = self.eped_prof[0,1] + self.eped_prof[0,0] * 2.*np.tanh(1)
		#teped = self.eped_prof[1,1] + self.eped_prof[1,0] * 2.*np.tanh(1)

		neped = self.eped_prof[0,1] + self.eped_prof[0,0] * (np.tanh(1)+np.tanh(2.0/self.eped_prof[0,4]*(width-0.5*self.eped_prof[0,4])))
		teped = self.eped_prof[1,1] + self.eped_prof[1,0] * (np.tanh(1)+np.tanh(2.0/self.eped_prof[1,4]*(width-0.5*self.eped_prof[1,4])))


		peped = (2.0 - (self.zeff-1.)/self.zimp) * neped * 1.e19 * teped * 1.e3 * 1.602 * 1.e-19

		bped = 4.*np.pi * 1.e-7 * abs(self.ip) * 1.e6 / self.perim

		bpped = peped / (bped**2) * 2 * (4.*np.pi*1.e-7)

		self.kbmc1 = width / bpped**float(self.e23.get()) #self.kbmc2
		print('>>> New fitting KBM Coef1 = %f (0.076 default)'%self.kbmc1)

		return

	def make_bnd_file(self,filename):

		f = open(filename,'w')
		f.write('%i\n'%len(self.rzbdy))
		for i in range(len(self.rzbdy)):
			f.write('%f\t%f\n'%(self.rzbdy[i,0],self.rzbdy[i,1]))

		f.close()

		return

	def write_eped_input(self):

		f = open('eped_opt','w')
		if self.iseqdsk:
			f.write('!--EQDSK -> '+ self.e2.get() +'\n')
			f.write('!--TARGET_BND_PSIN -> %f\n'%self.target_bnd)
		f.write('!-- Plasma geometry.\n')
		f.write('EQDSK = None \n')
		if (self.CheckVar2.get() == 1):
			f.write('USE_PARAM_SHAPE = True \n')
			f.write('ELONG  = %s\n'%self.e5.get())
			f.write('TRIA   = %s\n'%self.e4.get())
			f.write('SQUARE = %s\n'%self.e6.get())
			f.write('AMINOR = %s\n'%self.e9.get())
			f.write('RMAJOR = %s\n'%self.e7.get())
			self.make_bnd_file('eped_bnd')
		else:
			f.write('USE_PARAM_SHAPE = False \n')
			f.write('BND_FILE = eped_bnd \n')
			self.make_bnd_file('eped_bnd')

		f.write('\n')
		f.write('!-- Plasma 0-d parameters \n')
		f.write('IP= %f \n'%(1.e6*float(self.e20.get())))
		f.write('BVAC= %s \n'%self.e21.get())
		f.write('ZIMP  = %s \n'%self.e25.get())
		f.write('ZEFF  = %s \n'%self.e24.get())
		f.write('AMAIN = %s \n'%self.e26.get())
		f.write('AIMP  = %s \n'%self.e27.get())
		if (self.StrVar2.get().lower() == 'bpol'):
			f.write('Beta_val = %s \n'%self.e18.get())
		elif (self.StrVar2.get().lower() == 'betan'):
			f.write('Beta_val = %s \n'%self.e8.get())
		else:
			f.write('Beta_val = %f \n'%self.eq.rmag)
		f.write('LITARGET = %s \n'%self.e19.get())	
		f.write(' \n')
		f.write('!-- Fast ion profile option \n')
		f.write('APF = %s \n'%self.StrVar35.get())
		f.write('BPF = %s \n'%self.StrVar36.get())
		f.write('CPF = %s \n'%self.StrVar37.get())
		f.write('DPF = %s \n'%self.StrVar38.get())
		f.write(' \n')
		f.write('!-- Current profile option \n')
		if (self.StrVar29.get().lower() == 'sauter'):
			f.write('USE_NEO = False \n')
			f.write('USE_HAGER = False \n')
			f.write('USE_CHANG = False \n')
		elif (self.StrVar29.get().lower() == 'csauter'):
			f.write('USE_NEO = False \n')
			f.write('USE_HAGER = False \n')
			f.write('USE_CHANG = True \n')
		elif (self.StrVar29.get().lower() == 'hager'):
			f.write('USE_NEO = False \n')
			f.write('USE_HAGER = True \n')
			f.write('USE_CHANG = True \n')
		else:
			f.write('USE_NEO = True \n')
			f.write('USE_HAGER = True \n')
			f.write('USE_CHANG = True \n')

		if (self.CheckVar5.get() == 1):
			f.write('HAG_CORE_MOD=True \n')
		else:
			f.write('HAG_CORE_MOD=False \n')
		f.write('HAG_CORE_MOD_PSIN = %s\n'%self.StrVar43.get())
		f.write('Core_neo = %s \n'%self.StrVar44.get())
		f.write('BSMULTI = %s \n'%self.StrVar60.get())
		f.write('AJF = %s \n'%self.StrVar39.get())
		f.write('BJF = %s \n'%self.StrVar40.get())
		f.write('CJF = %s \n'%self.StrVar41.get())
		f.write('DJF = %s \n'%self.StrVar42.get())
		f.write('Current_ITERN = %s \n'%self.StrVar45.get())
		f.write('RELAX = %s \n'%self.StrVar46.get())
		f.write(' \n')
		f.write('!-- Eped_opt \n')
		f.write('KBM_COEF = %s \n'%self.e22.get())
		f.write('KBM_POWER = %s \n'%self.e23.get())
		f.write('NEPED = %s \n'%self.e35.get())
		f.write('NESEP = %s  \n'%self.e17.get())
		f.write('TESEP = %s \n'%self.e16.get())
		f.write('W_MIN = %s \n'%self.e31.get())
		f.write('W_MAX = %s \n'%self.e32.get())
		f.write('W_STEP = %s \n'%self.e33.get())
		f.write('AT1 = %s \n'%self.e14.get())
		f.write('AN1 = %s \n'%self.e15.get())
		f.write('ALPT1 = %s \n'%self.e10.get())
		f.write('ALPT2 = %s \n'%self.e11.get())
		f.write('ALPN1 = %s \n'%self.e12.get())
		f.write('ALPN2 = %s \n'%self.e13.get())

		f.write('!-- Stab opt \n')
		f.write('MODEN = %s \n'%self.e29.get())
		f.write('RUN_STAB = %s \n'%self.StrVar27.get().lower())
		f.write('NS = %s \n'%self.StrVar47.get())
		f.write('NT = %s \n'%self.StrVar48.get())
		f.write('MAP_NS = %s \n'%self.StrVar49.get())
		f.write('MAP_NT = %s \n'%self.StrVar50.get())
		f.write('MIS_GRIDN = %s \n'%self.StrVar51.get())
		f.write('MIS_PSISTART = %s \n'%self.StrVar52.get())
		f.write('MIS_XR1 = %s \n'%self.StrVar53.get())
		f.write('MIS_SIG1 = %s \n'%self.StrVar54.get())
		f.write('MIS_XR2 = %s \n'%self.StrVar55.get())
		f.write('MIS_SIG2 = %s \n'%self.StrVar56.get())
		if (self.CheckVar6.get() == 1):
			f.write('ELI_COMP = True \n')
		else:
			f.write('ELI_COMP = False \n')
		f.write('ELI_GRIDN = %s \n'%self.StrVar57.get())
		f.write('ELI_PSISTART = %s \n'%self.StrVar59.get())
		f.write('ELI_NDIST = %s \n'%self.StrVar58.get())
		f.write('NQA = %s \n'%self.e34.get())
		f.write('EXCLUDE = \n')
		f.write(' \n')
		f.write('!-- run option (EXT_VLOOP = 0.0 then automatic run) \n')
		f.write('NODE_LIST = %s \n'%self.nodelist)
		f.write(' \n')
		f.write('!-- convergence option\n')
		f.write('EPSILON = %s \n'%self.StrVar64.get())
		f.write('Beta_crit = %s \n'%self.StrVar65.get())
		f.write('Li_crit = %s \n'%self.StrVar66.get())
		f.write('Ip_crit = %s \n'%self.StrVar67.get())
		f.write('Bs_crit = %s \n'%self.StrVar68.get())
		f.write('Li_ITERN = %s \n'%self.StrVar69.get())


		f.write('!-- default option (Do not adjust it!) \n')
		line = 'RUN_MODE = EPED \n'
		if (self.CheckVar3.get() == 1):
			line = 'RUN_MODE = NORMAL \n'
		if (self.CheckVar1.get() == 1):
			if not(self.CheckVar7.get() == 1):
				line = 'RUN_MODE = EPED2 \n'
			else:
				line = 'RUN_MODE = EPED3 \n'

		f.write(line)

		if (self.StrVar2.get().lower() == 'bpol'):
			f.write('Beta_criterion_type = 1 \n')
		elif (self.StrVar2.get().lower() == 'betan'):
			f.write('Beta_criterion_type = 3 \n')
		else:
			f.write('Beta_criterion_type = 2 \n')

		f.write('ADJUST_PROF = False \n')
		f.write('kinetic_profile_type = 3 \n')
		f.write('chease_kinetic_file = chease_eped \n')
		f.write('NIDEAL = 8\n')
		f.write('TEPED = %f\n'%self.teped)
		f.write('TEWIDTH = %f\n'%self.tewidth)
		f.write('NCUTOFF = %i\n'%self.ncrit)
		f.write('GRCRIT1 = 0.03 \n')
		f.write('GRCRIT2 = 0.25 \n')
		f.write('QDELFIX = %s \n'%self.e75.get())
		if self.use_raw_prof:
			f.write('USERAW = True\n')
		else:
			f.write('USERAW = False\n')

		f.write('TIFILE = PROFILES/TI_fit.dat\n')

		f.close()

		return

	def open_result_list(self):

		#if self.outdir == '-search-':
		run_list =['-search-']
#		else:
#			run_list =[self.outdir,'-search-']

		filename = os.getcwd()+'/run.log'
		if os.path.isfile(filename):
			f = open(filename,'r')
			while True:
				line = f.readline()
				if not line: break
				line = line.split('\n')
				run_list.append(line[0])
			f.close()

		return run_list

	def modify_eped_opt(self):

		f = open('eped_opt','r')
		f2 = open('eped_opt_temp','w')

		while True:
			line = f.readline()
			if not line: break

			if (line.find('NQA') > -1):
				line = 'NQA = %s \n'%self.e65.get()
			elif (line.find('NCUTOFF') > -1):
				line = 'NCUTOFF = %s \n'%self.e64.get()
			elif (line.find('GRCRIT1') > -1):
				line = 'GRCRIT1 = %s \n'%self.e66.get()
			elif (line.find('GRCRIT2') > -1):
				line = 'GRCRIT2 = %s \n'%self.e67.get()	
			elif (line.find('EXCLUDE') > -1):
				line = 'EXCLUDE = %s \n'%self.e74.get()
		
			f2.write(line)

		f.close()
		f2.close()
		os.remove('eped_opt')
		move('eped_opt_temp','eped_opt')

		return

	def button_func1(self,skip=False,file=None):
		if not skip:
			self.input = askopenfilename()
			if (len(self.input) == 0):
				return
		else:
			self.input = file
		self.eq = eqdsk.eqdsk(self.input)
		self.eq.read_eqdsk_file()
		self.iseqdsk = True

		self.draw_psirz(self.iseqdsk)
		self.scroll_func1()

		self.update_psirz()

		#filename
		self.e2.delete(0, 'end' )
		self.e2.insert(10,self.input)

		R_max = (max(self.eq.rzbdy[:,0])+min(self.eq.rzbdy[:,0]))*0.5
		Z_max = 0.
		R_min = (max(self.eq.rzbdy[:,0])+min(self.eq.rzbdy[:,0]))*0.5
		Z_min = 0.
		ind = np.zeros(4)

		for i in range(len(self.eq.rzbdy[:,0])):

			if (self.eq.rzbdy[i,0] > R_max):
				R_max = self.eq.rzbdy[i,0]
				ind[0] = i
			if (self.eq.rzbdy[i,0] < R_min):
				R_min = self.eq.rzbdy[i,0]
				ind[1] = i
			if (self.eq.rzbdy[i,1] > Z_max):
				Z_max = self.eq.rzbdy[i,1]
				ind[2] = i
			if (self.eq.rzbdy[i,1] < Z_min):
				Z_min = self.eq.rzbdy[i,1]
				ind[3] = i								

		self.z0 = self.eq.rzbdy[int(ind[0]),1]
		self.r0 = (R_max + R_min) * 0.5
		self.a0 = (R_max - R_min) * 0.5
		self.elong = (Z_max - Z_min) / self.a0 / 2
		self.tri = (self.r0 - self.eq.rzbdy[int(ind[2]),0])/self.a0
		self.square = 0.

		self.z0 = round(self.z0,3)
		self.r0 = round(self.r0,3)
		self.a0 = round(self.a0,3)
		self.elong = round(self.elong,3)
		self.tri = round(self.tri,3)
		self.square = round(self.square,3)

		self.StrVar3.set(self.trans_vars(self.tri,2))
		self.StrVar4.set(self.trans_vars(self.elong,2))
		self.StrVar5.set(self.trans_vars(self.square,2))
		self.StrVar6.set(self.trans_vars(self.r0,2))
		self.StrVar8.set(self.trans_vars(self.a0,2))
		self.e4.delete(0, 'end' )
		self.e5.delete(0, 'end' )
		self.e6.delete(0, 'end' )
		self.e7.delete(0, 'end' )
		self.e9.delete(0, 'end' )

		self.e4.insert(10,self.StrVar3.get())
		self.e5.insert(10,self.StrVar4.get())
		self.e6.insert(10,self.StrVar5.get())
		self.e7.insert(10,self.StrVar6.get())
		self.e9.insert(10,self.StrVar8.get())

		self.eq.epslon = 4
		self.eq.make_chease_input()
		self.eq.run_chease()
		os.remove('chease_namelist')
		os.remove('EXPEQ')
		self.bp = round(self.eq.bp,3)
		self.li = round(self.eq.li,3)
		self.ip = round(abs(self.eq.ip) / 1.e6,3)
		self.bt = round(abs(self.eq.bcentr),3)
		self.betan = round(self.eq.betan,3)
		self.StrVar7.set(self.trans_vars(self.betan,2))
		self.StrVar17.set(self.trans_vars(self.bp,2))
		self.StrVar18.set(self.trans_vars(self.li,2))
		self.StrVar19.set(self.trans_vars(self.ip,2))
		self.StrVar20.set(self.trans_vars(self.bt,2))

		self.e8.delete(0, 'end' )
		self.e8.insert(10,self.StrVar7.get())
		self.e18.delete(0, 'end' )
		self.e18.insert(10,self.StrVar17.get())
		self.e19.delete(0, 'end' )
		self.e19.insert(10,self.StrVar18.get())
		self.e20.delete(0, 'end' )
		self.e20.insert(10,self.StrVar19.get())		
		self.e21.delete(0, 'end' )
		self.e21.insert(10,self.StrVar20.get())

		self.perim = self.eq.perim

		return		

	def button_func2(self):

		filename = askopenfilename()
		if (len(filename) == 0):
			return

		self.use_ext_bnd = True

		f = open(filename,'r')
		line = f.readline()
		dat_num = int(line.split()[0])

		self.rzbdy = np.zeros(shape=(dat_num,2))
		for i in range(dat_num):
			line = f.readline().split()
			self.rzbdy[i,0] = float(line[0])
			self.rzbdy[i,1] = float(line[1])

		self.update_psirz()
		return

	def button_func3(self):

		if (self.CheckVar2.get() == 1):
			self.use_param_shape = True
		else:
			self.use_param_shape = False		

		if not (self.use_param_shape):
			if not self.use_ext_bnd:
				self.scroll_func1()
			return

		self.chi = np.linspace(0,2*np.pi,256)
		self.rzbdy = np.zeros(shape=(256,2))

		self.StrVar3.set(self.e4.get())
		self.StrVar4.set(self.e5.get())
		self.StrVar5.set(self.e6.get())
		self.StrVar6.set(self.e7.get())
		self.StrVar8.set(self.e9.get())

		self.tri = float(self.StrVar3.get())
		self.elong = float(self.StrVar4.get())
		self.square = float(self.StrVar5.get())
		self.r0 = float(self.StrVar6.get())
		self.a0 = float(self.StrVar8.get())


		for i in range(256):

			theta = self.chi[i]
			
			self.rzbdy[i,0] = self.r0 + self.a0 * np.cos(theta + self.tri * np.sin(theta) + self.square*np.sin(2.0*theta))
			self.rzbdy[i,1] = self.elong * self.a0 * np.sin(theta + 0.0*self.square*np.cos(theta)) + self.z0

		self.update_psirz()		
		return			

	def button_func4(self):

		os.system(gfit2_dir)
		
		self.renew_eped_coef()

		self.use_raw_prof = True

		return

	def button_func5(self):


		self.zeff = float(self.e24.get())
		self.zimp = float(self.e25.get())
		self.amain = float(self.e26.get())
		self.aimp = float(self.e27.get())

		self.make_kbm_coef()

		self.kbmc1 = round(self.kbmc1,3)

		if (self.iseqdsk):
			print('>>> Tped = %f [keV], nped = %f [10(19)/m3], width = %f'%(self.teped,self.neped,self.tewidth))
		#	print('>>> New fitting KBM Coef1 = %f (0.076 default)'%self.kbmc1)
		else:
			print('>>> Have to load eqdsk first...')

		self.StrVar21.set(self.trans_vars(self.kbmc1,2))

		self.e22.delete(0, 'end' )
		self.e22.insert(10,self.StrVar21.get())

		return

	def button_func6(self):

		if not self.t1_close:
			print('>>> NUM Param Window is already opened...')
			return

		self.t1 = tk.Toplevel(self.root)
		self.t1.wm_title("NUM Params")
		self.t1_close = False		
		self.t1.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(1))


		self.l2 = tk.Label(self.t1, text = '====== Beam profile ======',justify='center')
		self.l2.grid(row=1,column=0,columnspan=4)

		self.l2 = tk.Label(self.t1, text = 'APF ',justify='center')
		self.l2.grid(row=2,column=0,columnspan=2)
		self.e36 = tk.Entry(self.t1,width=5,justify='center')
		self.e36.insert(0,self.StrVar35.get())
		self.e36.grid(row=2, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'BPF [0-1]',justify='center')
		self.l2.grid(row=3,column=0,columnspan=2)
		self.e37 =tk.Entry(self.t1,width=5,justify='center')
		self.e37.insert(0,self.StrVar36.get())
		self.e37.grid(row=3, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'CPF ',justify='center')
		self.l2.grid(row=4,column=0,columnspan=2)
		self.e38 = tk.Entry(self.t1,width=5,justify='center')
		self.e38.insert(0,self.StrVar37.get())
		self.e38.grid(row=4, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'DPF ',justify='center')
		self.l2.grid(row=5,column=0,columnspan=2)
		self.e39 = tk.Entry(self.t1,width=5,justify='center')
		self.e39.insert(0,self.StrVar38.get())
		self.e39.grid(row=5, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t1, text = 'AJF ',justify='center')
		self.l2.grid(row=6,column=0,columnspan=2)
		self.e40 = tk.Entry(self.t1,width=5,justify='center')
		self.e40.insert(0,self.StrVar39.get())
		self.e40.grid(row=6, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t1, text = 'BJF [0-1]',justify='center')
		self.l2.grid(row=7,column=0,columnspan=2)
		self.e41 = tk.Entry(self.t1,width=5,justify='center')
		self.e41.insert(0,self.StrVar40.get())
		self.e41.grid(row=7, column=2, columnspan = 2)						

		self.l2 = tk.Label(self.t1, text = 'CJF ',justify='center')
		self.l2.grid(row=8,column=0,columnspan=2)
		self.e42 = tk.Entry(self.t1,width=5,justify='center')
		self.e42.insert(0,self.StrVar41.get())
		self.e42.grid(row=8, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'DJF ',justify='center')
		self.l2.grid(row=9,column=0,columnspan=2)
		self.e43 = tk.Entry(self.t1,width=5,justify='center')
		self.e43.insert(0,self.StrVar42.get())
		self.e43.grid(row=9, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = '====== Current option ======',justify='center')
		self.l2.grid(row=10,column=0,columnspan=4)		

		self.l2 = tk.Label(self.t1, text = 'HAGMOD ',justify='center')
		self.l2.grid(row=11,column=0,columnspan=2)
		self.c5 = tk.Checkbutton(self.t1,variable=self.CheckVar5)
		self.c5.grid(row=11, column=2, columnspan = 2)
		self.l2 = tk.Label(self.t1, text = 'LI_ITER2 ',justify='center')
		self.l2.grid(row=12,column=0,columnspan=2)
		self.c7 = tk.Checkbutton(self.t1,variable=self.CheckVar7)
		self.c7.grid(row=12, column=2, columnspan = 2)		
		self.l2 = tk.Label(self.t1, text = 'MODPSIN [0-1]',justify='center')
		self.l2.grid(row=13,column=0,columnspan=2)
		self.e44 = tk.Entry(self.t1,width=5,justify='center')
		self.e44.insert(0,self.StrVar43.get())
		self.e44.grid(row=13, column=2, columnspan = 2)
		self.l2 = tk.Label(self.t1, text = 'CORENEO [0-1]',justify='center')
		self.l2.grid(row=14,column=0,columnspan=2)
		self.e45 = tk.Entry(self.t1,width=5,justify='center')
		self.e45.insert(0,self.StrVar44.get())
		self.e45.grid(row=14, column=2, columnspan = 2)	
		self.l2 = tk.Label(self.t1, text = 'BSMULTI [0-1]',justify='center')
		self.l2.grid(row=15,column=0,columnspan=2)
		self.e61 = tk.Entry(self.t1,width=5,justify='center')
		self.e61.insert(0,self.StrVar60.get())
		self.e61.grid(row=15, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t1, text = '====== Numerical option ======',justify='center')
		self.l2.grid(row=16,column=0,columnspan=4)		

		self.l2 = tk.Label(self.t1, text = 'ITERN [#]',justify='center')
		self.l2.grid(row=17,column=0,columnspan=2)
		self.e46 = tk.Entry(self.t1,width=5,justify='center')
		self.e46.insert(0,self.StrVar45.get())
		self.e46.grid(row=17, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t1, text = 'RELAX [0-1]',justify='center')
		self.l2.grid(row=18,column=0,columnspan=2)
		self.e47 = tk.Entry(self.t1,width=5,justify='center')
		self.e47.insert(0,self.StrVar46.get())
		self.e47.grid(row=18, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'NS [#]',justify='center')
		self.l2.grid(row=19,column=0,columnspan=2)
		self.e48 = tk.Entry(self.t1,width=5,justify='center')
		self.e48.insert(0,self.StrVar47.get())
		self.e48.grid(row=19, column=2, columnspan = 2)				

		self.l2 = tk.Label(self.t1, text = 'NT [#]',justify='center')
		self.l2.grid(row=20,column=0,columnspan=2)
		self.e49 = tk.Entry(self.t1,width=5,justify='center')
		self.e49.insert(0,self.StrVar48.get())
		self.e49.grid(row=20, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t1, text = 'NS_MAP [#]',justify='center')
		self.l2.grid(row=21,column=0,columnspan=2)
		self.e50 = tk.Entry(self.t1,width=5,justify='center')
		self.e50.insert(0,self.StrVar49.get())
		self.e50.grid(row=21, column=2, columnspan = 2)				

		self.l2 = tk.Label(self.t1, text = 'NT_MAP [#]',justify='center')
		self.l2.grid(row=22,column=0,columnspan=2)
		self.e51 = tk.Entry(self.t1,width=5,justify='center')
		self.e51.insert(0,self.StrVar50.get())
		self.e51.grid(row=22, column=2, columnspan = 2)	


		self.l2 = tk.Label(self.t1, text = '====== MISHKA opt. ======',justify='center')
		self.l2.grid(row=1,column=4,columnspan=4)

		self.l2 = tk.Label(self.t1, text = 'GRIDN ',justify='center')
		self.l2.grid(row=2,column=4,columnspan=2)
		self.e52 = tk.Entry(self.t1,width=5,justify='center')
		self.e52.insert(0,self.StrVar51.get())
		self.e52.grid(row=2, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'PSISTART ',justify='center')
		self.l2.grid(row=3,column=4,columnspan=2)
		self.e53 = tk.Entry(self.t1,width=5,justify='center')
		self.e53.insert(0,self.StrVar52.get())
		self.e53.grid(row=3, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'XR1 ',justify='center')
		self.l2.grid(row=4,column=4,columnspan=2)
		self.e54 = tk.Entry(self.t1,width=5,justify='center')
		self.e54.insert(0,self.StrVar53.get())
		self.e54.grid(row=4, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'SIG1 ',justify='center')
		self.l2.grid(row=5,column=4,columnspan=2)
		self.e55 = tk.Entry(self.t1,width=5,justify='center')
		self.e55.insert(0,self.StrVar54.get())
		self.e55.grid(row=5, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'XR2 ',justify='center')
		self.l2.grid(row=6,column=4,columnspan=2)
		self.e56 = tk.Entry(self.t1,width=5,justify='center')
		self.e56.insert(0,self.StrVar55.get())
		self.e56.grid(row=6, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'SIG2 ',justify='center')
		self.l2.grid(row=7,column=4,columnspan=2)
		self.e57 = tk.Entry(self.t1,width=5,justify='center')
		self.e57.insert(0,self.StrVar56.get())
		self.e57.grid(row=7, column=6, columnspan = 2)						


		self.l2 = tk.Label(self.t1, text = '====== ELITE opt. ======',justify='center')
		self.l2.grid(row=10,column=4,columnspan=4)

		self.l2 = tk.Label(self.t1, text = 'COMPRESS ',justify='center')
		self.l2.grid(row=11,column=4,columnspan=2)
		self.c6 = tk.Checkbutton(self.t1,variable=self.CheckVar6)
		self.c6.grid(row=11, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t1, text = 'GRIDN ',justify='center')
		self.l2.grid(row=12,column=4,columnspan=2)
		self.e58 = tk.Entry(self.t1,width=5,justify='center')
		self.e58.insert(0,self.StrVar57.get())
		self.e58.grid(row=12, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t1, text = 'NDIST ',justify='center')
		self.l2.grid(row=13,column=4,columnspan=2)
		self.e59 = tk.Entry(self.t1,width=5,justify='center')
		self.e59.insert(0,self.StrVar58.get())
		self.e59.grid(row=13, column=6, columnspan = 2)			

		self.l2 = tk.Label(self.t1, text = 'PSISTART ',justify='center')
		self.l2.grid(row=14,column=4,columnspan=2)
		self.e60 = tk.Entry(self.t1,width=5,justify='center')
		self.e60.insert(0,self.StrVar59.get())
		self.e60.grid(row=14, column=6, columnspan = 2)			

		self.l2 = tk.Label(self.t1, text = 'PSISTART ',justify='center')
		self.l2.grid(row=14,column=4,columnspan=2)
		self.e60 = tk.Entry(self.t1,width=5,justify='center')
		self.e60.insert(0,self.StrVar59.get())
		self.e60.grid(row=14, column=6, columnspan = 2)		

		self.l2 = tk.Label(self.t1, text = '====== Converge opt. ======',justify='center')
		self.l2.grid(row=16,column=4,columnspan=4)

		self.l2 = tk.Label(self.t1, text = 'EPSLON ',justify='center')
		self.l2.grid(row=17,column=4,columnspan=2)
		self.e68 = tk.Entry(self.t1,width=5,justify='center')
		self.e68.insert(0,self.StrVar64.get())
		self.e68.grid(row=17, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t1, text = 'Beta crit ',justify='center')
		self.l2.grid(row=18,column=4,columnspan=2)
		self.e69 = tk.Entry(self.t1,width=5,justify='center')
		self.e69.insert(0,self.StrVar65.get())
		self.e69.grid(row=18, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t1, text = 'Li crit ',justify='center')
		self.l2.grid(row=19,column=4,columnspan=2)
		self.e70 = tk.Entry(self.t1,width=5,justify='center')
		self.e70.insert(0,self.StrVar66.get())
		self.e70.grid(row=19, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t1, text = 'Ip crit ',justify='center')
		self.l2.grid(row=20,column=4,columnspan=2)
		self.e71 = tk.Entry(self.t1,width=5,justify='center')
		self.e71.insert(0,self.StrVar67.get())
		self.e71.grid(row=20, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t1, text = 'BS crit ',justify='center')
		self.l2.grid(row=21,column=4,columnspan=2)
		self.e72 = tk.Entry(self.t1,width=5,justify='center')
		self.e72.insert(0,self.StrVar68.get())
		self.e72.grid(row=21, column=6, columnspan = 2)						

		self.l2 = tk.Label(self.t1, text = 'ITERL ',justify='center')
		self.l2.grid(row=22,column=4,columnspan=2)
		self.e73 = tk.Entry(self.t1,width=5,justify='center')
		self.e73.insert(0,self.StrVar69.get())
		self.e73.grid(row=22, column=6, columnspan = 2)	

		b1 = tk.Button(self.t1, text="SAVE", bg = "lightgray",command=lambda: self.button_func8(),height = 1,width = 4)
		b1.grid(row=23, column=0,columnspan=8,pady=15)


		return

	def button_func7(self):

		if node_force:	
			self.nodelist = node_default
			return

		load2,used2,name2,nastat2 = read_node_status()

		if not self.t2_close:
			print('>>> NOTE LIST Window is already opened...')
			return

		self.t2 = tk.Toplevel(self.root)
		self.t2.wm_title("NODE LIST")
		self.t2_close = False		
		self.t2.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(2))

		self.l2 = tk.Label(self.t2, text = '====== NODE LIST ======',justify='center')
		self.l2.grid(row=0,column=0,columnspan=8)		

		self.l2 = tk.Label(self.t2, text = '--Node--',justify='center')
		self.l2.grid(row=1,column=0,columnspan=2,padx=5)
		self.l2 = tk.Label(self.t2, text = '--Load--',justify='center')
		self.l2.grid(row=1,column=2,columnspan=2,padx=5)
		self.l2 = tk.Label(self.t2, text = '--Used--',justify='center')
		self.l2.grid(row=1,column=4,columnspan=2,padx=5)				

		for i in range(len(load2)):

			self.l2 = tk.Label(self.t2, text = name2[i] ,justify='center')
			self.l2.grid(row=(i+2),column=0,columnspan=2)

			self.l2 = tk.Label(self.t2, text = load2[i] ,justify='center')
			self.l2.grid(row=(i+2),column=2,columnspan=2)	

			if not nastat2[i] == 'na':
				line = '%i/%i'%(used2[i,1],used2[i,2])
			else:
				line = 'N/A'
			self.l2 = tk.Label(self.t2, text = line ,justify='center')
			self.l2.grid(row=(i+2),column=4,columnspan=2)


		self.noden = len(load2)

		for i in range(10,10+len(load2)):

			self.__dict__['c%d'%i] = tk.Checkbutton(self.t2,variable=self.__dict__['CheckVar%d'%i])
			self.__dict__['c%d'%i].grid(row=i-8,column=6,columnspan=2)

		b1 = tk.Button(self.t2, text="SAVE", bg = "lightgray",command=lambda: self.button_func9(),height = 1,width = 4)
		b1.grid(row=self.noden+4, column=0,columnspan=8,pady=15)			

		return

	def button_func8(self):

		for i in range(35,61):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%(i+1)].get())

		for i in range(64,70):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%(i+4)].get())			

		self.t1.destroy()
		self.t1_close = True

		return	

	def button_func9(self):

		nodelist = ''
		count = 0
		for i in range(10,10+self.noden):

			if (self.__dict__['CheckVar%d'%i].get() == 1):

				if count == 0:
					nodelist = 'node%02i'%(i-9)
				else:
					nodelist = nodelist + ', node%02i'%(i-9)

				count = count + 1
		if count == 0:
			print(">>> Have to choose more than one...")
			return
		else:
			print(">>> Selected nodes ->",nodelist)
			self.nodelist = nodelist

		self.t2.destroy()
		self.t2_close = True

		return								

	def button_func10(self):

		self.renew_eped_coef()
		self.use_raw_prof = True

		return

	def button_func11(self):

		if not self.rinput:
			print('>>> Input is not ready...')
			return

		self.write_eped_input()	

		run_dir = os.getcwd() + '/' + self.e62.get()

		try:
			os.mkdir(run_dir)
		except:
			print('>>> There is dir already!')
			return

		move('eped_bnd',run_dir + '/eped_bnd')
		move('eped_opt',run_dir + '/eped_opt')
		try:
			move('fit_opt',run_dir + '/fit_opt')
		except:
			pass
		try:
			move('fit_opt_param',run_dir + '/fit_opt_param')
		except:
			pass
		if (self.use_raw_prof):
			try:
				copytree('PROFILES',run_dir + '/PROFILES')
			except:
				pass

		if self.iseqdsk:
			copyfile(self.e2.get(),run_dir + '/geqdsk')
	
		filename = run_dir+'/eped_bat'
		log_e = run_dir+'/eped_batch.e'
		log_o = run_dir+'/eped_batch.o'
		nodelist = node_default
		command = 'cd '+ run_dir + '\n'
		command = command + eped_dir + ' eped_opt \n'
		make_batch_script(filename,None,log_e,log_o,command,'EPED')
		
		
		currdir = os.getcwd()
		#self.root.quit()
		self.root.destroy()
		os.chdir(run_dir)
		os.system('chmod 777 eped_bat')
		if not (self.CheckVar4.get() == 1):		
			os.system('./eped_bat ')
		else:
			runid = submit_batch_script('eped_bat')
		os.chdir(currdir)
		return

	def button_func12(self):

		try:
			result_dir = self.StrVar62.get().split()[0]
		except:
			return

		result_dir2 = os.getcwd() + '/'+result_dir

		if result_dir == '-search-':
			input = askdirectory()
			if len(input) == 0:
				return
			result_dir = input
			result_dir2 = input

			ise = False
			isfile = True
			input = input.split('/')
			try:
				f = open('run.log','r')
			except:
				f = open('run.log','w')
				f.close()
				isfile = False
			while isfile:
				line = f.readline()
				if not line: break
				if (line.find(input[-1])>-1):
					ise = True
			if isfile: f.close()
			if not ise:
				os.system('echo ' + input[-1] + ' >> run.log')
			self.run_list = self.open_result_list()
			menu = self.e63["menu"]
			menu.delete(0, "end")
			for string in self.run_list:
				menu.add_command(label=string, command=lambda value=string: self.StrVar62.set(value))

		print('>>> Result dir >',result_dir2)
		if not os.path.isdir(result_dir2):
			print('No available result')
			return

		currdir = os.getcwd()
		os.chdir(result_dir2)

		if not (os.path.isfile('log.batch_stab') or os.path.isfile('result.dat')) :
			print('>>> No available result')
			return			

		if not self.t3_close:
			print('>>> NOTE LIST Window is already opened...')
			return

		self.sim = eped.eped('eped_opt',False)
		self.sim.make_equ = False
		self.sim.plot_only = False
		self.sim.make_ja = False

		print('>>> Collect_data...')
		self.sim.collect_data()

		self.t3 = tk.Toplevel(self.root)
		self.t3.wm_title("OPEN RESULT")
		self.t3_close = False		
		self.t3.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(3))


		self.t3.resizable(0,0)
		self.fig2, [ax1, ax2, ax3] = plt.subplots(1,3,figsize=(14,4.5))

		self.fig2.tight_layout()

		self.canvas2 = FigureCanvasTkAgg(self.fig2,master=self.t3)
		self.plot_widget2 = self.canvas2.get_tk_widget()
		self.plot_widget2.grid(rowspan=32,row=1,column=10,columnspan=10)

		toolbar_frame = tk.Frame(self.t3)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar2 = NavigationToolbar2Tk(self.canvas2,toolbar_frame)
		self.fig2.canvas.draw_idle()	

		self.l2 = tk.Label(self.t3, text = '======== POST PROCESSING ========',justify='center')
		self.l2.grid(row=0,column=0,columnspan=8)
		b8 = tk.Button(self.t3,  text="Collect data", bg = "light gray",command=lambda: self.button_func16(),height = 1,width = 15)
		b8.grid(row=1, column=0, columnspan=8)	
		b9 = tk.Button(self.t3,  text="Clear log", bg = "light gray",command=lambda: self.sim.clear_files(),height = 1,width = 15)
		b9.grid(row=2, column=0, columnspan=8)	
		b10 = tk.Button(self.t3, text="Draw result", bg = "light gray",command=lambda: self.button_func14(),height = 1,width = 15)
		b10.grid(row=3, column=0,columnspan=8)			
		b11 = tk.Button(self.t3, text="Make Equ", bg = "light gray",command=lambda: self.button_func15(),height = 1,width = 15)
		b11.grid(row=4, column=0,columnspan=8)
		b12 = tk.Button(self.t3, text="Draw bnd", bg = "light gray",command=lambda: self.button_func17(),height = 1,width = 15)
		b12.grid(row=6, column=0,columnspan=8)
		b13 = tk.Button(self.t3, text="Make Prof", bg = "light gray",command=lambda: self.button_func18(),height = 1,width = 15)
		b13.grid(row=5, column=0,columnspan=8)
		b14 = tk.Button(self.t3, text="EXIT", bg = "light gray",command=lambda: self.button_func20(self.t3,3),height = 1,width = 15)
		b14.grid(row=7, column=0,columnspan=8)			

		self.l4 = tk.Label(self.t3, text="NCUT",anchor='e')	
		self.l4.grid(row=8, column=0,columnspan=2)
		self.e64 = tk.Entry(self.t3,width=5,justify='center')
		self.e64.insert(10,str(self.sim.ncrit))
		self.e64.grid(row=8, column=2,columnspan=2)

		self.l4 = tk.Label(self.t3, text="NQA",anchor='e')	
		self.l4.grid(row=8, column=4,columnspan=2)
		self.e65 = tk.Entry(self.t3,width=5,justify='center')
		self.e65.insert(10,str(self.sim.nqcrit))
		self.e65.grid(row=8, column=6,columnspan=2)

		self.l4 = tk.Label(self.t3, text="CRIT1",anchor='e')	
		self.l4.grid(row=9, column=0,columnspan=2)
		self.e66 = tk.Entry(self.t3,width=5,justify='center')
		self.e66.insert(10,str(self.sim.grcrit1))
		self.e66.grid(row=9, column=2,columnspan=2)

		self.l4 = tk.Label(self.t3, text="CRIT2",anchor='e')	
		self.l4.grid(row=9, column=4,columnspan=2)
		self.e67 = tk.Entry(self.t3,width=5,justify='center')
		self.e67.insert(10,str(self.sim.grcrit2))
		self.e67.grid(row=9, column=6,columnspan=2)

		self.l4 = tk.Label(self.t3, text='EXC',anchor='e')
		self.l4.grid(row=10, column=0,columnspan=2)
		self.e74 = tk.Entry(self.t3,width=15,justify='center')
	
		line = ''
		for i in range(len(self.sim.exclude)):
			line = line + '%i'%(self.sim.exclude[i])
			if not i == (len(self.sim.exclude)-1):
				line = line + ','

		self.e74.insert(10,line)
		self.e74.grid(row=10, column=2,columnspan=6)

		frame = tk.Frame(self.t3)
		frame.grid(row=11,column=0,columnspan=8,rowspan=14)
		scrollbar = tk.Scrollbar(frame)
		scrollbar.pack(side="right", fill="y")
		scrollbar2 = tk.Scrollbar(frame,orient='horizontal')
		scrollbar2.pack(side="bottom", fill="x")
		listbox = tk.Listbox(frame,yscrollcommand = scrollbar.set,xscrollcommand = scrollbar2.set,width=30)
	
		os.chdir(result_dir2)
		f = open('eped_opt','r')
		count = 1
		while True:
			line = f.readline()
			if not line: break
			listbox.insert(count,line.split('\n')[0])
			count = count + 1
		f.close()
		listbox.pack(side="left", fill="x")
		scrollbar["command"]=listbox.yview
		scrollbar2['command'] = listbox.xview
		#frame.pack()

		return

	def button_func13(self):

		self.rinput = True

		if (len(self.rzbdy)==1):
			self.rinput = False
			print('>>> Plasma boundary is not generated!')

		if(self.nodelist == 'none'):
			self.rinput = False	
			print('>>> Nodelist is not ready!')

		if (os.path.isdir(self.e62.get())):
			self.rinput = False
			print('>>> Run_dir already exits -> Change the name')

		try:	
			if (len(self.e62.get().split()) == 0):
				self.rinput = False
				print('>>> Run name is not chosed')
		except:
			if (self.e62.get() == ''):
				self.rinput = False
				print('>>> Run name is not chosed')		


		if self.rinput:
			print('>>> Input is ready!')
			self.l100.config(text='READY',fg='lime')
		else:
			self.l100.config(text='NOT READY',fg='red')

		return

	def button_func14(self):
		self.fig2.canvas.draw_idle()	
		self.sim.plot_result(self.sim.gra,self.sim.grd,self.sim.width1,self.sim.width2,self.sim.fit1,self.sim.fit2,self.fig2)
		return

	def button_func15(self):

		if not self.t5_close:
			print('>>> Equil Window is already opened...')
			return

		self.t5 = tk.Toplevel(self.root)
		self.t5.wm_title("Equil selection")
		self.t5_close = False
		self.t5.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(5))

		self.t5.resizable(0,0)

		self.l2 = tk.Label(self.t5, text = '======= Equil select =======',justify='center')
		self.l2.grid(row=0,column=0,columnspan=3)		

		line = 1

		if self.sim.width1 >0.:
			b2 = tk.Button(self.t5, text="USE CRIT1", bg = "lightgray",command=lambda: self.button_func152(self.sim.width1),height = 1,width = 7)
			b2.grid(row=line, column=0,columnspan=3)
			line = line + 1
		if self.sim.width2 >0.:
			b3 = tk.Button(self.t5, text="USE CRIT2", bg = "lightgray",command=lambda: self.button_func152(self.sim.width1),height = 1,width = 7)
			b3.grid(row=line, column=0,columnspan=3)
			line = line + 1		

		b3 = tk.Button(self.t5, text="EXIT", bg = "lightgray",command=lambda: self.button_func19() ,height = 1,width = 7)
		b3.grid(row=line, column=0,columnspan=3)	
	

		return

	def button_func152(self,width,type=1):
		self.t5_close = True
		self.t5.destroy()
		self.sim.make_equ = False
		self.sim.collect_data()
		self.sim.make_result_equ(width,0)
		move('result/geqdsk1','result/geqdsk_eped')

		print('>>> Equilibrium is saved in %s/result/geqdsk_eped'%os.getcwd())
		self.button_func19()
		return		

	def button_func16(self):

		self.sim.nqcrit = float(self.e65.get())
		self.sim.ncrit = float(self.e64.get())
		self.sim.grcrit1 = float(self.e66.get())
		self.sim.grcrit2 = float(self.e67.get())
		exclude = self.e74.get().replace(" ","")
		self.sim.exclude = np.array([])
		if not exclude == '':
			exclude = exclude.split(',')
			for i in range(len(exclude)):
				self.sim.exclude = np.append(self.sim.exclude,int(exclude[i]))

		self.modify_eped_opt()
		self.sim.collect_data()
		return

	def button_func17(self):

		self.t4 = tk.Toplevel(self.root)
		self.t4.wm_title("Plasma Boundary")
		self.t4_close = False		
		self.t4.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(4))

		self.t4.resizable(0,0)
		self.fig3, ax1 = plt.subplots(1,1,figsize=(5,9))
		ax1.set_title('Equilibrium')
		ax1.set_xlabel('R [m]')
		ax1.set_ylabel('Z [m]')

		self.fig3.tight_layout()
		self.t4.resizable(0,0)		

		self.canvas3 = FigureCanvasTkAgg(self.fig3,master=self.t4)
		self.plot_widget3 = self.canvas3.get_tk_widget()
		self.plot_widget3.grid(rowspan=20,row=1,column=10,columnspan=10)

		toolbar_frame = tk.Frame(self.t4)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar3 = NavigationToolbar2Tk(self.canvas3,toolbar_frame)


		f = open('eped_bnd','r')
		line = f.readline()
		dat_num = int(line)
		bnd = np.zeros(shape=(dat_num,2))
		for i in range(dat_num):
			line = f.readline().split()
			bnd[i,0] = float(line[0])
			bnd[i,1] = float(line[1])
		f.close()

		ax1.plot(bnd[:,0],bnd[:,1],'--',color='pink')
		ax1.axis('scaled')


		R_max = (max(bnd[:,0])+min(bnd[:,0]))*0.5
		Z_max = 0.
		R_min = (max(bnd[:,0])+min(bnd[:,0]))*0.5
		Z_min = 0.
		ind = np.zeros(4)

		for i in range(dat_num):

			if (bnd[i,0] > R_max):
				R_max = bnd[i,0]
				ind[0] = i
			if (bnd[i,0] < R_min):
				R_min = bnd[i,0]
				ind[1] = i
			if (bnd[i,1] > Z_max):
				Z_max = bnd[i,1]
				ind[2] = i
			if (bnd[i,1] < Z_min):
				Z_min = bnd[i,1]
				ind[3] = i

		r0 = (R_min + R_max) * 0.5
		a0 = (R_max - R_min) * 0.5
	
		tri_u  = (r0 - bnd[int(ind[2]),0])/a0
		tri_l  = (r0 - bnd[int(ind[3]),0])/a0
		elong  = (Z_max - Z_min) / a0 / 2.

		inform = 'tri_u = %4.3f\n'%tri_u
		inform = inform + 'tri_l = %4.3f\n'%tri_l
		inform = inform + 'elong = %4.3f\n'%elong
		inform = inform + 'r0    = %4.3f\n'%r0
		inform = inform + 'a0    = %4.3f\n'%a0

		ax1.text(r0,0.0,inform,color='r')
		plt.draw()	

		return

	def button_func18(self):

		self.istifile = False
		self.tiprof = np.zeros(shape=(101,2))
		if not self.t5_close:
			print('>>> Kin_gen Window is already opened...')
			return

		if (self.sim.width1 == 0. and self.sim.width2 ==0.):
			print('>>> EPED result is unavailable...')
			return

		if self.sim.use_raw_prof:
				f = open('PROFILES/chease_eped_mod','r')
				try:
					f2 = open(self.sim.ti_file,'r')
					self.istifile = True
				except:
					pass
		else:
				f = open('equil/scan_0/chease_eped','r')
		eped_coef1 = np.zeros((3,8))

		for i in range(3):
			line = f.readline().split()
			for j in range(8):
				eped_coef1[i,j] = float(line[j])

		try:
			line = f.readline().split()
			self.zeff = float(line[0])
			self.zimp = float(line[1])
			self.amain = float(line[2])
			self.aimp = float(line[3])

			self.e24.delete(0, 'end')
			self.e24.insert(0, self.zeff)
			self.e25.delete(0, 'end')
			self.e25.insert(0, self.zimp)
			self.e26.delete(0, 'end')
			self.e26.insert(0, self.amain)
			self.e27.delete(0, 'end')
			self.e27.insert(0, self.aimp)									
		except:
			pass	
			
		f.close()		
	
		if self.istifile:
			line = f2.readline()
			line = f2.readline()
			for i in range(101):
				line = f2.readline().split()
				self.tiprof[i,0] = float(line[1])
				self.tiprof[i,1] = float(line[2])/1.e3
				
			f2.close()	

		eped_coef2 = np.copy(eped_coef1)
		psin = np.linspace(0,1,401)

		self.eped_coef = np.copy(eped_coef1)

		eped_coef1[2,1] = eped_coef1[1,1]
		eped_coef2[2,1] = eped_coef2[1,1]		

		for i in range(3):
			eped_coef1[i,3] = 1.0 - 0.5*self.sim.width1
			eped_coef1[i,4] = self.sim.width1
			eped_coef1[i,5] = 1.0 - self.sim.width1
			eped_coef2[i,3] = 1.0 - 0.5*self.sim.width2
			eped_coef2[i,4] = self.sim.width2
			eped_coef2[i,5] = 1.0 - self.sim.width2	

		for i in range(1,3):
			eped_coef1[i,0] = (self.sim.tped1-eped_coef1[i,1])/2./np.tanh(1)
			eped_coef2[i,0] = (self.sim.tped2-eped_coef2[i,1])/2./np.tanh(1)

		for i in range(1,3):
			if (i == 2 and self.istifile):
				eped_coef1[i,2] = self.tiprof[0,1] - eped_coef1[i,0]*(1.+np.tanh(1)) - eped_coef1[i,1]
				eped_coef2[i,2] = self.tiprof[0,1] - eped_coef2[i,0]*(1.+np.tanh(1)) - eped_coef2[i,1]
			else:
				eped_coef1[i,2] = self.eped_coef[i,2] + self.eped_coef[i,0]*(1.+np.tanh(1)) + self.eped_coef[i,1] - eped_coef1[i,0]*(1.+np.tanh(1)) - eped_coef1[i,1]
				eped_coef2[i,2] = self.eped_coef[i,2] + self.eped_coef[i,0]*(1.+np.tanh(1)) + self.eped_coef[i,1] - eped_coef2[i,0]*(1.+np.tanh(1)) - eped_coef2[i,1]
				
		if self.sim.use_raw_prof:	
			for i in range(3):
				if self.sim.width1 > 0.:
					if i==2:	eped_coef1[i,6], eped_coef1[i,7] = self.eped_fit(self.eped_coef[i,:],eped_coef1[i,:],'ti')
					else:		eped_coef1[i,6], eped_coef1[i,7] = self.eped_fit(self.eped_coef[i,:],eped_coef1[i,:])
					
				if self.sim.width2 > 0.:
					if i==2:	eped_coef2[i,6], eped_coef2[i,7] = self.eped_fit(self.eped_coef[i,:],eped_coef2[i,:],'ti')
					else:		eped_coef2[i,6], eped_coef2[i,7] = self.eped_fit(self.eped_coef[i,:],eped_coef2[i,:])
				

		if self.sim.width1 > 0.:
			self.ne1 = self.sim.eped_fun(psin,eped_coef1[0,:])
			self.te1 = self.sim.eped_fun(psin,eped_coef1[1,:])
			self.ti1 = self.sim.eped_fun(psin,eped_coef1[2,:])

		if self.sim.width2 > 0.:
			self.ne2 = self.sim.eped_fun(psin,eped_coef2[0,:])
			self.te2 = self.sim.eped_fun(psin,eped_coef2[1,:])
			self.ti2 = self.sim.eped_fun(psin,eped_coef2[2,:])

		if self.sim.use_raw_prof:
			self.ne3 = self.sim.eped_fun(psin,self.eped_coef[0,:])
			self.te3 = self.sim.eped_fun(psin,self.eped_coef[1,:])
			self.ti3 = self.sim.eped_fun(psin,self.eped_coef[2,:])

		self.t5 = tk.Toplevel(self.root)
		self.t5.wm_title("Kinetic prof")
		self.t5_close = False
		self.t5.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(5))

		self.t5.resizable(0,0)
		self.fig5, [ax1, ax2, ax3] = plt.subplots(1,3,figsize=(13,5))

		self.canvas5 = FigureCanvasTkAgg(self.fig5,master=self.t5)
		self.plot_widget5 = self.canvas5.get_tk_widget()
		self.plot_widget5.grid(rowspan=40,row=1,column=10,columnspan=40)

		toolbar_frame = tk.Frame(self.t5)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar5 = NavigationToolbar2Tk(self.canvas5,toolbar_frame)
		self.fig5.canvas.draw_idle()

		legend1 = []

		try:
			f = open('PROFILES/NE_raw.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,3))
			for i in range(num):
				line = f.readline().split()
				dat[i,0] = float(line[0])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
			f.close()
			ax1.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='gold',ecolor='gold',capthick=2)
		except:
			pass

		try:
			f = open('PROFILES/TE_raw.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,3))
			for i in range(num):
				line = f.readline().split()
				dat[i,0] = float(line[0])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
			f.close()
			ax2.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='gold',ecolor='gold',capthick=2)
		except:
			pass

		try:
			f = open('PROFILES/TI_raw.dat','r')
			num = int(float(f.readline()))
			dat = np.zeros(shape=(num,3))
			for i in range(num):
				line = f.readline().split()
				dat[i,0] = float(line[0])
				dat[i,1] = float(line[1])
				dat[i,2] = float(line[2])
			f.close()
			ax3.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='gold',ecolor='gold',capthick=2)
		except:
			pass

		max1 = 0.; max2 = 0.; max3 = 0.

		if self.sim.width1 > 0.:

			p1= ax1.plot(psin,self.ne1)
			p2= ax2.plot(psin,self.te1)
			p3= ax3.plot(psin,self.ti1)
			max1 = max(self.ne1)*1.2
			max2 = max(self.te1)*1.2
			max3 = max(self.ti1)*1.2
			legend1.append('Crit1')

		if self.sim.width2 > 0.:

			p1= ax1.plot(psin,self.ne2)
			p2= ax2.plot(psin,self.te2)
			p3= ax3.plot(psin,self.ti2)
			max1 = max(self.ne2)*1.2
			max2 = max(self.te2)*1.2
			max3 = max(self.ti2)*1.2						
			legend1.append('Crit2')

		if self.sim.use_raw_prof:

			p1= ax1.plot(psin,self.ne3)
			p2= ax2.plot(psin,self.te3)
			if self.istifile:	p3= ax3.plot(self.tiprof[:,0],self.tiprof[:,1])
			else:	p3= ax3.plot(psin,self.ti3)
			legend1.append('RAW FIT')

		ax1.set_title("NE [10(19)/m3]")
		ax1.set_xlabel('Normalized radius ($\psi_N$)')
		ax1.set_ylabel('NE [10(19)/m3]')
		ax1.legend(legend1)
		ax1.set_xlim((-0.05,1.05))			
		ax1.set_ylim((-0.1,max1))

		ax2.set_title("TE [keV]")
		ax2.set_xlabel('Normalized radius ($\psi_N$)')
		ax2.set_ylabel('TE [keV]')
		ax2.legend(legend1)
		ax2.set_xlim((-0.05,1.05))
		ax2.set_ylim((-0.1,max2))

		ax3.set_title("TI [keV]")
		ax3.set_xlabel('Normalized radius ($\psi_N$)')
		ax3.set_ylabel('TI [keV]')
		ax3.legend(legend1)
		ax3.set_xlim((-0.05,1.05))
		ax3.set_ylim((-0.1,max3))

		self.fig5.tight_layout()

		line = 2

		if self.sim.use_raw_prof:
			b1 = tk.Button(self.t5, text="USE RAW", bg = "lightgray",command=lambda: self.write_kinprof(psin,3),height = 1,width = 7)
			b1.grid(row=line, column=1,columnspan=3)
			line = line + 1
		if self.sim.width1 >0.:
			b2 = tk.Button(self.t5, text="USE CRIT1", bg = "lightgray",command=lambda: self.write_kinprof(psin,1),height = 1,width = 7)
			b2.grid(row=line, column=1,columnspan=3)
			line = line + 1
		if self.sim.width2 >0.:
			b3 = tk.Button(self.t5, text="USE CRIT2", bg = "lightgray",command=lambda: self.write_kinprof(psin,2),height = 1,width = 7)
			b3.grid(row=line, column=1,columnspan=3)
			line = line + 1

		b3 = tk.Button(self.t5, text="EXIT", bg = "lightgray",command=lambda: self.button_func19() ,height = 1,width = 7)
		b3.grid(row=line, column=1,columnspan=3)			

		return

	def button_func19(self):

		self.t5_close = True
		self.t5.destroy()
		return	

	def button_func20(self,obj,values):

		if values == 8:	self.t8_close = True
		elif values == 7:	self.t7_close = True
		elif values == 6:	self.t6_close = True
		elif values == 5:	self.t5_close = True
		elif values == 4:	self.t4_close = True
		elif values == 3:	
			self.t3_close = True
			os.chdir(self.currdir)
		elif values == 2:	self.t2_close = True
		elif values == 1:	self.t1_close = True
		obj.destroy()
		return

	def eped_fun(self,x,a7,a8):
		
		val = self.temp[1]
		val = val + self.temp[0]*(np.tanh(2.0/self.temp[4]*(1.0-self.temp[3]))-np.tanh(2.0/self.temp[4]*(x-self.temp[3]))) 
		val = val + self.temp[2] * ((1.0 - (x/self.temp[5])**(1.01+abs(a7)))**(1.01+abs(a8))) * 0.5 * (1.0+np.sign(self.temp[5]-x))
		return val
	
	def eped_fun2(self,x,a7,a8):
		
		val = self.temp[1]
		val = val + self.temp[0]*(np.tanh(2.0/self.temp[4]*(1.0-self.temp[3]))-np.tanh(2.0/self.temp[4]*(x-self.temp[3])))
		val = val + self.temp[2] * ((1.0 - (x/self.temp[5])**(abs(a7)))**(abs(a8))) * 0.5 * (1.0+np.sign(self.temp[5]-x))
		return val

	def eped_fit(self,eped_prof1,eped_prof2,flag='ne'):
		
		if (flag.lower() == 'ti' and self.istifile):
			num = int(max(0.5,eped_prof2[5]-0.30)*100.)
			psin = self.tiprof[0:num,0]
		else:	
			num = 60
			psin = np.linspace(0.0,max(0.5,eped_prof2[5]-0.30),num)

		prof = np.copy(psin)
		self.temp = np.copy(eped_prof1)

		if (flag.lower() == 'ti' and self.istifile):
			prof = self.tiprof[0:num,1]
		else:
			for i in range(num):
				prof[i] = self.eped_fun(psin[i],eped_prof1[6]-1.01,eped_prof1[7]-1.01)

		self.temp = np.copy(eped_prof2)
		count = 0
		if (flag.lower() == 'ti' and self.istifile):	pp = [1.1,2.1]
		else:	pp = [0.,1.0]
		while count < 8:
			if (flag.lower() == 'ti' and self.istifile):	
				popt, pcov = curve_fit(self.eped_fun2,psin,prof,p0=pp,maxfev=300000)
			else:	
				popt, pcov = curve_fit(self.eped_fun,psin,prof,p0=pp,maxfev=300000)	
				
			if (flag.lower() == 'ti' and self.istifile):
				if popt[1] == pp[1]:
					count = count + 1
					if count == 1:	pp = [0.,0.]
					elif count == 2:	pp = [2.1,1.1]
					elif count == 3:	pp = [2.1,2.1]
					elif count == 4:	pp = [1.6,1.6]
					elif count == 5:	pp = [eped_prof1[6],eped_prof1[7]]
					elif count == 6:	pp = [1.1,1.6]
					elif count == 7:	pp = [1.6,1.1]
				else:
					count = 8
			else:
				if popt[1] == pp[1]:
					count = count + 1
					if count == 1:	pp = [0.,0.]
					elif count == 2:	pp = [1.,0.]
					elif count == 3:	pp = [1.,1.]
					elif count == 4:	pp = [.5,.5]
					elif count == 5:	pp = [eped_prof1[6]-1.01,eped_prof1[7]-1.01]
					elif count == 6:	pp = [0.,.5]
					elif count == 7:	pp = [.5,0.]
				else:
					count = 8


		if popt[1] == pp[1]:
			popt[0] = abs(eped_prof1[6]) - 1.01
			popt[1] = abs(eped_prof1[7]) - 1.01
			print('Fitting failed')
			print(popt[0],popt[1])
		if (flag.lower() == 'ti' and self.istifile):	return (abs(popt[0]),abs(popt[1]))	
		else:	return (abs(popt[0])+1.01,abs(popt[1])+1.01)

	def write_kinprof(self,psin,type=1):

		if type==1:
			y1 = self.ne1
			y2 = self.ti1
			y3 = self.te1
		elif type==2:
			y1 = self.ne2
			y2 = self.ti2
			y3 = self.te2
		else:
			y1 = self.ne3
			y2 = self.ti3
			y3 = self.te3

		try:
			os.mkdir('result')
		except:
			pass

		f = open('result/chease_kinprof_eped','w')
		f.write('401\n')
		f.write('%f\t%f\t%f\t%f\n'%(self.sim.ch.zeff,self.sim.ch.zimp,self.sim.ch.amain,self.sim.ch.aimp))
		denw = 1.0 - (self.sim.ch.zeff-1.)/self.sim.ch.zimp
		for i in range(401):
			f.write('%f\t%f\t%f\t%f\t%f\n'%(psin[i],y3[i],y1[i],y2[i],y1[i]*denw))
		f.close()
		f = open('result/logk.out','w')
		if type==1:	f.write('type1\n')
		elif type ==2:	f.write('type2\n')
		else:		f.write('type0\n')
		f.close()

		print('>>> Profile is saved in %s/result/chease_kinprof_eped'%os.getcwd())
		return

	def gui_eped(self):

		self.fig, ax1 = plt.subplots(1,1,figsize=(5,9))
		ax1.set_title('Equilibrium')
		ax1.set_xlabel('R [m]')
		ax1.set_ylabel('Z [m]')

		self.fig.tight_layout()
		self.root.resizable(0,0)
		self.canvas = FigureCanvasTkAgg(self.fig,master=self.root)
		self.plot_widget = self.canvas.get_tk_widget()
		self.plot_widget.grid(rowspan=40,row=1,column=10,columnspan=10)

		toolbar_frame = tk.Frame(self.root)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar = NavigationToolbar2Tk(self.canvas,toolbar_frame)

		self.l1 = tk.Label(self.root, text="============== EQDSK option ==============",justify='center')
		self.l1.grid(row=0, column=0,columnspan=9)

		self.l2 = tk.Label(self.root, text="DIR ",anchor='e')
		self.e2 = tk.Entry(self.root,width=30,justify='center')
		self.e2.insert(10,self.StrVar1.get())
		self.l2.grid(row=1, column=0)
		self.e2.grid(row=1, column=1,columnspan=7)
		b1 = tk.Button(self.root, text="Load Eq", bg = "lightgray",command=lambda: self.button_func1(),height = 1,width = 4)
		b1.grid(row=1, column=8)

		self.l3 = tk.Label(self.root, text="=============== Constraint ===============",justify='center')
		self.l3.grid(row=2, column=0,columnspan=9)
		self.l4 = tk.Label(self.root, text="Beta const ",anchor='e')	
		self.l4.grid(row=3, column=0,columnspan=2)
		self.e3 = tk.OptionMenu(self.root,self.StrVar2,'BPOL','BETAN','RMAG')
		self.e3.config(width=5)
		self.e3.grid(row=3,column=2,columnspan=3)
		self.l5 = tk.Label(self.root, text="li const ",anchor='e')
		self.l5.grid(row=3,column=6,columnspan=2)
		self.c1 = tk.Checkbutton(self.root,variable=self.CheckVar1, command=lambda: self.change_checkvar())
		self.c1.grid(row=3, column=8)

		self.l6 = tk.Label(self.root, text="============== Equilibrium ==============",justify='center')
		self.l6.grid(row=4, column=0,columnspan=9)

		self.l8 = tk.Label(self.root, text="PARAM_BND ",anchor='e')	
		self.l8.grid(row=5, column=0,columnspan=4)
		self.c2 = tk.Checkbutton(self.root,variable=self.CheckVar2)
		self.c2.grid(row=5, column=4)

		b2 = tk.Button(self.root, text="OPEN", bg = "lightgray",command=lambda: self.button_func2(),height = 1,width = 4)
		b2.grid(row=5, column=6,columnspan=2)		

		b3 = tk.Button(self.root, text="SAVE", bg = "lightgray",command=lambda: self.button_func3(),height = 1,width = 4)
		b3.grid(row=5, column=8,columnspan=2)

		self.l7 = tk.Label(self.root, text="BND_PSI ",anchor='e')	
		self.l7.grid(row=6, column=0,columnspan=3,sticky='se')
		self.s1 = tk.Scale(self.root, variable=self.DoubleVar1, command=lambda x: self.scroll_func1(), orient='horizontal', showvalue=True, from_=0.980,to=0.999,resolution=0.001,length=150)
		self.s1.grid(row=6,column=3,columnspan=7,sticky='s')

		self.l9 = tk.Label(self.root, text="TRI",anchor='e')
		self.l9.grid(row=7, column=0,columnspan=2)
		self.e4 = tk.Entry(self.root,width=5,justify='center')
		self.e4.insert(10,self.StrVar3.get())
		self.e4.grid(row=8, column=0,columnspan=2)

		self.l10 = tk.Label(self.root, text="ELON",anchor='e')
		self.l10.grid(row=7, column=2,columnspan=2)		
		self.e5 = tk.Entry(self.root,width=5,justify='center')
		self.e5.insert(10,self.StrVar4.get())
		self.e5.grid(row=8, column=2,columnspan=2)

		self.l11 = tk.Label(self.root, text="SQUA",anchor='e')
		self.l11.grid(row=7, column=4,columnspan=2)	
		self.e6 = tk.Entry(self.root,width=5,justify='center')
		self.e6.insert(10,self.StrVar5.get())
		self.e6.grid(row=8, column=4,columnspan=2)

		self.l12 = tk.Label(self.root, text="R0 [m]",anchor='e')
		self.l12.grid(row=7, column=6,columnspan=2)	
		self.e7 = tk.Entry(self.root,width=5,justify='center')
		self.e7.insert(10,self.StrVar6.get())
		self.e7.grid(row=8, column=6,columnspan=2)

		self.l13 = tk.Label(self.root, text="Betan",anchor='e')
		self.l13.grid(row=9, column=0,columnspan=2)	
		self.e8 = tk.Entry(self.root,width=5,justify='center')
		self.e8.insert(10,self.StrVar7.get())
		self.e8.grid(row=10, column=0,columnspan=2)

		self.l14 = tk.Label(self.root, text="a0 [m]",anchor='e')
		self.l14.grid(row=7, column=8,columnspan=2)							
		self.e9 = tk.Entry(self.root,width=5,justify='center')
		self.e9.insert(10,self.StrVar8.get())	
		self.e9.grid(row=8, column=8,columnspan=2)	

		self.l14 = tk.Label(self.root, text="Bpol",anchor='e')
		self.l14.grid(row=9, column=2,columnspan=2)							
		self.e18 = tk.Entry(self.root,width=5,justify='center')
		self.e18.insert(10,self.StrVar17.get())	
		self.e18.grid(row=10, column=2,columnspan=2)	

		self.l14 = tk.Label(self.root, text="li",anchor='e')
		self.l14.grid(row=9, column=4,columnspan=2)							
		self.e19 = tk.Entry(self.root,width=5,justify='center')
		self.e19.insert(10,self.StrVar18.get())	
		self.e19.grid(row=10, column=4,columnspan=2)					

		self.l14 = tk.Label(self.root, text="Ip [MA]",anchor='e')
		self.l14.grid(row=9, column=6,columnspan=2)							
		self.e20 = tk.Entry(self.root,width=5,justify='center')
		self.e20.insert(10,self.StrVar19.get())	
		self.e20.grid(row=10, column=6,columnspan=2)	

		self.l14 = tk.Label(self.root, text="B0 [T]",anchor='e')
		self.l14.grid(row=9, column=8,columnspan=2)							
		self.e21 = tk.Entry(self.root,width=5,justify='center')
		self.e21.insert(10,self.StrVar20.get())	
		self.e21.grid(row=10, column=8,columnspan=2)	

		self.l15 = tk.Label(self.root, text="================ Profile ================",justify='center')
		self.l15.grid(row=11, column=0,columnspan=9)
		self.l16 = tk.Label(self.root, text="Fix core ",anchor='e')	
		self.l16.grid(row=12, column=0,columnspan=3)
		self.c3= tk.Checkbutton(self.root,variable=self.CheckVar3,command=lambda: self.change_checkvar2())
		self.c3.grid(row=12, column=3)
		b4 = tk.Button(self.root, text="USE RAW", bg = "lightgray",command=lambda: self.button_func4(),height = 1,width = 5)
		b4.grid(row=12, column=8,columnspan=2)
		b8 = tk.Button(self.root, text="Load Prev", bg = "lightgray",command=lambda: self.button_func10(),height = 1,width = 7)
		b8.grid(row=12, column=5,columnspan=3)		

		self.l17 = tk.Label(self.root, text="SEP",anchor='e')
		self.l17.grid(row=13, column=2,columnspan=2)
		self.l18 = tk.Label(self.root, text="CORE",anchor='e')
		self.l18.grid(row=13, column=4,columnspan=2)	
		self.l19 = tk.Label(self.root, text="ALP1",anchor='e')
		self.l19.grid(row=13, column=6,columnspan=2)	
		self.l20 = tk.Label(self.root, text="ALP2",anchor='e')
		self.l20.grid(row=13, column=8,columnspan=2)	
										
		self.l21 = tk.Label(self.root, text="TE",anchor='e')
		self.l21.grid(row=14, column=0,columnspan=2)
		self.l22 = tk.Label(self.root, text="NE",anchor='e')
		self.l22.grid(row=15, column=0,columnspan=2)

		self.e16 = tk.Entry(self.root,width=5,justify='center')
		self.e16.insert(10,self.StrVar15.get())
		self.e16.grid(row=14, column=2,columnspan=2)

		self.e14 = tk.Entry(self.root,width=5,justify='center')
		self.e14.insert(10,self.StrVar13.get())
		self.e14.grid(row=14, column=4,columnspan=2)			

		self.e10 = tk.Entry(self.root,width=5,justify='center')
		self.e10.insert(10,self.StrVar9.get())
		self.e10.grid(row=14, column=6,columnspan=2)	

		self.e11 = tk.Entry(self.root,width=5,justify='center')
		self.e11.insert(10,self.StrVar10.get())
		self.e11.grid(row=14, column=8,columnspan=2)	

		self.e17 = tk.Entry(self.root,width=5,justify='center')
		self.e17.insert(10,self.StrVar16.get())
		self.e17.grid(row=15, column=2,columnspan=2)	

		self.e15 = tk.Entry(self.root,width=5,justify='center')
		self.e15.insert(10,self.StrVar14.get())
		self.e15.grid(row=15, column=4,columnspan=2)	

		self.e12 = tk.Entry(self.root,width=5,justify='center')
		self.e12.insert(10,self.StrVar11.get())
		self.e12.grid(row=15, column=6,columnspan=2)	

		self.e13 = tk.Entry(self.root,width=5,justify='center')
		self.e13.insert(10,self.StrVar12.get())
		self.e13.grid(row=15, column=8,columnspan=2)

		self.l25 = tk.Label(self.root, text="NPED",anchor='e')
		self.l25.grid(row=16, column=0,columnspan=2)
		self.l25 = tk.Label(self.root, text="ZEFF",anchor='e')
		self.l25.grid(row=16, column=2,columnspan=2)
		self.l26 = tk.Label(self.root, text="ZIMP",anchor='e')
		self.l26.grid(row=16, column=4,columnspan=2)	
		self.l27 = tk.Label(self.root, text="AMAIN",anchor='e')
		self.l27.grid(row=16, column=6,columnspan=2)	
		self.l28 = tk.Label(self.root, text="AIMP",anchor='e')
		self.l28.grid(row=16, column=8,columnspan=2)	

		self.e35 = tk.Entry(self.root,width=5,justify='center')
		self.e35.insert(10,self.StrVar34.get())
		self.e35.grid(row=17, column=0,columnspan=2)

		self.e24 = tk.Entry(self.root,width=5,justify='center')
		self.e24.insert(10,self.StrVar23.get())
		self.e24.grid(row=17, column=2,columnspan=2)

		self.e25 = tk.Entry(self.root,width=5,justify='center')
		self.e25.insert(10,self.StrVar24.get())
		self.e25.grid(row=17, column=4,columnspan=2)

		self.e26 = tk.Entry(self.root,width=5,justify='center')
		self.e26.insert(10,self.StrVar25.get())
		self.e26.grid(row=17, column=6,columnspan=2)

		self.e27 = tk.Entry(self.root,width=5,justify='center')
		self.e27.insert(10,self.StrVar26.get())
		self.e27.grid(row=17, column=8,columnspan=2)

		self.l23 = tk.Label(self.root, text="============= KBM constraint =============",justify='center')
		self.l23.grid(row=18, column=0,columnspan=9)

		self.l24 = tk.Label(self.root, text="Coef1",anchor='e')
		self.l24.grid(row=19, column=0,columnspan=2)
		self.e22 = tk.Entry(self.root,width=5,justify='center')
		self.e22.insert(10,self.StrVar21.get())
		self.e22.grid(row=19, column=2,columnspan=2)	

		self.l25 = tk.Label(self.root, text="Coef2",anchor='e')
		self.l25.grid(row=19, column=4,columnspan=2)
		self.e23 = tk.Entry(self.root,width=5,justify='center')
		self.e23.insert(10,self.StrVar22.get())
		self.e23.grid(row=19, column=6,columnspan=2 )
		b5 = tk.Button(self.root, text="USE RAW", bg = "lightgray",command=lambda: self.button_func5(),height = 1,width = 5)
		b5.grid(row=19, column=8,columnspan=2)		


		self.l26 = tk.Label(self.root, text="================ Stability ================",justify='center')
		self.l26.grid(row=20, column=0,columnspan=9)
		self.e28 = tk.OptionMenu(self.root,self.StrVar27,'Mishka','Elite')
		self.e28.config(width=5)
		self.e28.grid(row=21, column=0,columnspan=4)

		self.l27 = tk.Label(self.root, text="n",anchor='e')
		self.l27.grid(row=21, column=4,columnspan=1)		
		self.e29 = tk.Entry(self.root,width=16,justify='center')
		self.e29.insert(10,self.StrVar28.get())
		self.e29.grid(row=21, column=4,columnspan=6)


		self.l28 = tk.Label(self.root, text="================ Options ================",justify='center')
		self.l28.grid(row=22, column=0,columnspan=9)
		self.l29 = tk.Label(self.root, text="BS model",anchor='e')
		self.l29.grid(row=23, column=0,columnspan=2)	
		self.e30 = tk.OptionMenu(self.root,self.StrVar29,'Sauter','CSauter','Hager','Neo')
		self.e30.config(width=5)
		self.e30.grid(row=23, column=2,columnspan=3)

		self.l30 = tk.Label(self.root, text="Wmin",anchor='e')
		self.l30.grid(row=24, column=0,columnspan=2)
		self.l31 = tk.Label(self.root, text="Wmax",anchor='e')
		self.l31.grid(row=25, column=0,columnspan=2)
		self.l32 = tk.Label(self.root, text="Wstep",anchor='e')
		self.l32.grid(row=26, column=0,columnspan=2)		
		self.e31 = tk.Entry(self.root,width=7,justify='center')
		self.e31.insert(10,self.StrVar30.get())
		self.e31.grid(row=24, column=2,columnspan=3)
		self.e32 = tk.Entry(self.root,width=7,justify='center')
		self.e32.insert(10,self.StrVar31.get())
		self.e32.grid(row=25, column=2,columnspan=3)
		self.e33 = tk.Entry(self.root,width=7,justify='center')
		self.e33.insert(10,self.StrVar32.get())
		self.e33.grid(row=26, column=2,columnspan=3)		

		self.l31 = tk.Label(self.root, text="Use batch",anchor='e')
		self.l31.grid(row=30, column=5,columnspan=3)
		self.c4= tk.Checkbutton(self.root,variable=self.CheckVar4)
		self.c4.grid(row=30, column=8)

		self.l32 = tk.Label(self.root, text="Nqa bilinear",anchor='e')
		self.l32.grid(row=25, column=5,columnspan=3)
		self.e34 = tk.Entry(self.root,width=5,justify='center')
		self.e34.insert(10,self.StrVar33.get())
		self.e34.grid(row=25, column=8,columnspan=2)		

		self.l32 = tk.Label(self.root, text="qdel",anchor='e')
		self.l32.grid(row=26, column=5,columnspan=3)
		self.e75 = tk.Entry(self.root,width=5,justify='center')
		self.e75.insert(10,self.StrVar70.get())
		self.e75.grid(row=26, column=8,columnspan=2)

		b6 = tk.Button(self.root, text="NUM_PARAM", bg = "lightgray",command=lambda: self.button_func6(),height = 1,width = 9)
		b6.grid(row=23, column=5,columnspan=5)	
		b7 = tk.Button(self.root, text="NODE_LIST", bg = "lightgray",command=lambda: self.button_func7(),height = 1,width = 9)
		b7.grid(row=24, column=5,columnspan=5)				

		self.l33 = tk.Label(self.root, text="================ Commands ================",justify='center')
		self.l33.grid(row=27, column=0,columnspan=9)

		self.l32 = tk.Label(self.root, text="RUN_NAME",anchor='e')
		self.l32.grid(row=28, column=0,columnspan=3)
		self.e62 = tk.Entry(self.root,width=20,justify='center')
		self.e62.insert(10,self.StrVar61.get())
		self.e62.grid(row=28, column=3,columnspan=5)		

		b9 = tk.Button(self.root, text="RUN EPED", bg = "lightgray",command=lambda: self.button_func11(),height = 1,width = 10)
		b9.grid(row=30, column=0, columnspan=5)

		b10 = tk.Button(self.root, text="CHECK INP", bg = "lightgray",command=lambda: self.button_func13(),height = 1,width = 10)
		b10.grid(row=29, column=0, columnspan=5)

		b11 = tk.Button(self.root, text="OPEN RESULT", bg = "lightgray",command=lambda: self.button_func12(),height = 1,width = 10)
		b11.grid(row=31, column=0, columnspan=5)

		b11 = tk.Button(self.root, text="EXIT", bg = "lightgray",command=lambda: self.button_func20(self.root,1),height = 1,width = 10)
		b11.grid(row=33, column=0, columnspan=5)			

		self.e63 = tk.OptionMenu(self.root,self.StrVar62,*self.run_list)
		#menu = self.e63.children["menu"]

		#for value in self.run_list:
		self.l100 = tk.Label(self.root, text="NOT READY",fg='red',anchor='center',bg='white',width=15,justify='left')
		self.l100.grid(row=29, column=5,columnspan=4)		

#		self.l32 = tk.Label(self.root, text="-OUTPUT LIST-",anchor='e')
#		self.l32.grid(row=31, column=0,columnspan=5)

		self.e63.config(width=10)
		self.e63.grid(row=31,column=5,columnspan=4)


		font = tkinter.font.Font(weight='bold',size=10)
		self.l34 = tk.Label(self.root, text="-Ver %s "%(version['eped']),anchor='w',width=20,justify='left',fg='dodgerblue',font=font)
		self.l34.grid(row=37, column=0,columnspan=10,sticky='sw')

		self.l35 = tk.Label(self.root, text="-%s"%(author['eped2']),anchor='w',width=40,justify='left',fg='dodgerblue',font=font)
		self.l35.grid(row=38, column=0,columnspan=10,sticky='sw')

		self.l36 = tk.Label(self.root, text="-%s"%(comment['eped']),anchor='w',width=40,justify='left',fg='magenta',font=font)
		self.l36.grid(row=39, column=0,columnspan=10,sticky='sw')

		self.read_preset(2)

		self.root.mainloop()
		return

	def initialise_vars(self):

		self.iseqdsk = False
		self.eped_prof = np.zeros(shape=(3,8))
		self.use_li = False						#check1
		self.use_ext_bnd = False				
		self.constraint = 'BPOL'				#str2
		self.fix_core = False					#check3
		self.t1_close = True
		self.t2_close = True
		self.t3_close = True
		self.t4_close = True
		self.t5_close = True

		self.epslon = 1.e-8
		self.beta_crit = 1.e-2
		self.li_crit = 2.e-2
		self.ip_crit = 1.e-5
		self.bs_crit = 1.e-5
		self.niterl = 10

		self.use_li_iter2 = True

		self.run_name = ''

		self.ncrit = 1

		self.eqdsk_name = None					#str1
		self.target_bnd = 0.995					#double1
		self.use_param_shape = False			#check2
		self.tri = 0.3							#str3
		self.elong = 1.4						#str4
		self.square = 0.0						#str5
		self.r0 = 1.8							#str6
		self.z0 = 0.0							#str7
		self.a0 = 0.5							#str8
		self.rzbdy = np.zeros(shape=(1,2))
		self.rinput = False

		self.betan = 1.0							#str7
		self.bp = 1.0								#str 17
		self.li = 0.8								#str 18
		self.ip = 0.6								#str 19
		self.bt = 2.0								#str 20
		self.zeff = 2.0								#str 23
		self.zimp = 6.0								#str 24
		self.amain = 2.0							#str 25
		self.aimp = 12.0							#str 26

		self.apf = 0.2								#str 35
		self.bpf = 0.8								#str 36		
		self.cpf = 2.5								#str 37		
		self.dpf = 3.0								#str 38		
		self.ajf = 0.2								#str 39
		self.bjf = 1.0								#str 40		
		self.cjf = 2.0								#str 41		
		self.djf = 2.0								#str 42

		self.hag_core_mod = True					#check5
		self.hag_core_mod_psin = 0.3				#str 43
		self.core_neo = 0.01						#str 44

		self.currenti = 35							#str 45
		self.relax = 0.8							#str 46
		self.ns = 150								#str 47	
		self.nt = 200								#str 48
		self.map_ns = 300							#str 49
		self.map_nt = 512							#str 50

		self.mis_gridn = 301						#str 51
		self.mis_psistart = 0.75					#str 52
		self.xr1 = 1.0								#str 53
		self.sig1 = 0.07							#str 54
		self.xr2 = 0.9								#str 55
		self.sig2 = 0.1								#str 56

		self.eli_comp = False						#str 57
		self.eli_gridn = 2000						#str 58
		self.eli_ndist = 50							#str 59
		self.eli_psistart = 0.5						#str 60

		#EPED opt
		self.kbmc1 = 0.076		#str 21
		self.kbmc2 = 0.5		#str 22
		self.alpt1 = 1.2		#str 9
		self.alpt2 = 1.4		#str 10
		self.alpn1 = 1.1		#str 11
		self.alpn2 = 1.2		#str 12
		self.at1 = 2.0			#str 13
		self.an1 = 1.0			#str 14
		self.wmin = 0.02		#str 30
		self.wmax = 0.07		#str 31
		self.wstep = 15			#str 32
		self.neped = 1.6		#str 34
		self.teped = 0.0
		self.tewidth = 0.0
		self.tsep  = 0.1		#str 15
		self.nsep  = 0.25	        #str 16
		self.bsmulti = 1.0		#str 61

		#STAB opt
		self.stab_code = 'Mishka'		#str 27
		self.mode_n = '5,7,10,15,20'	#str 28
		self.bsmodel = 'csauter'			#str 29
		self.nqa = 27.7					#str 33
		self.batch_run = True			#check 4
		self.nodelist = node_init
		self.qdelfix = 0.3
		self.use_raw_prof = False

		self.outdir = '-search-'
		self.ti_file = None

		try:	self.read_eped_opt('EPED/eped_opt')
		except:	print('>>> No prescribed setting 1')

		try:	self.read_eped_opt('eped_opt')
		except:	print('>>> No prescribed setting 2') 			

		self.read_preset()

		self.DoubleVar1.set(self.target_bnd)
		
		self.StrVar1.set(self.trans_vars(self.eqdsk_name,3))	
		self.StrVar2.set(self.trans_vars(self.constraint,1))			
		self.StrVar3.set(self.trans_vars(self.tri,2))
		self.StrVar4.set(self.trans_vars(self.elong,2))
		self.StrVar5.set(self.trans_vars(self.square,2))
		self.StrVar6.set(self.trans_vars(self.r0,2))	
		self.StrVar7.set(self.trans_vars(self.betan,2))
		self.StrVar8.set(self.trans_vars(self.a0,2))

		self.CheckVar1.set(self.trans_vars(self.use_li,4))
		self.CheckVar2.set(self.trans_vars(self.use_param_shape,4))
		self.CheckVar3.set(self.trans_vars(self.fix_core,4))
		self.CheckVar4.set(self.trans_vars(self.batch_run,4))
		self.CheckVar5.set(self.trans_vars(self.hag_core_mod,4))
		self.CheckVar6.set(self.trans_vars(self.eli_comp,4))
		self.CheckVar7.set(self.trans_vars(self.use_li_iter2,4))

		self.StrVar9.set(self.trans_vars(self.alpt1,2))
		self.StrVar10.set(self.trans_vars(self.alpt2,2))
		self.StrVar11.set(self.trans_vars(self.alpn1,2))
		self.StrVar12.set(self.trans_vars(self.alpn2,2))

		self.StrVar13.set(self.trans_vars(self.at1,2))
		self.StrVar14.set(self.trans_vars(self.an1,2))		

		self.StrVar15.set(self.trans_vars(self.tsep,2))
		self.StrVar16.set(self.trans_vars(self.nsep,2))

		self.StrVar17.set(self.trans_vars(self.bp,2))
		self.StrVar18.set(self.trans_vars(self.li,2))
		self.StrVar19.set(self.trans_vars(self.ip,2))
		self.StrVar20.set(self.trans_vars(self.bt,2))

		self.StrVar21.set(self.trans_vars(self.kbmc1,2))
		self.StrVar22.set(self.trans_vars(self.kbmc2,2))	

		self.StrVar23.set(self.trans_vars(self.zeff,2))
		self.StrVar24.set(self.trans_vars(self.zimp,2))	
		self.StrVar25.set(self.trans_vars(self.amain,2))
		self.StrVar26.set(self.trans_vars(self.aimp,2))
		self.StrVar27.set(self.trans_vars(self.stab_code,1))
		self.StrVar28.set(self.trans_vars(self.mode_n,1))	
		self.StrVar29.set(self.trans_vars(self.bsmodel,1))
		self.StrVar30.set(self.trans_vars(self.wmin,2))
		self.StrVar31.set(self.trans_vars(self.wmax,2))
		self.StrVar32.set(self.trans_vars(self.wstep,2))
		self.StrVar33.set(self.trans_vars(self.nqa,2))	
		self.StrVar34.set(self.trans_vars(self.neped,2))

		self.StrVar35.set(self.trans_vars(self.apf,2))
		self.StrVar36.set(self.trans_vars(self.bpf,2))	
		self.StrVar37.set(self.trans_vars(self.cpf,2))
		self.StrVar38.set(self.trans_vars(self.dpf,2))
		self.StrVar39.set(self.trans_vars(self.ajf,2))
		self.StrVar40.set(self.trans_vars(self.bjf,2))	
		self.StrVar41.set(self.trans_vars(self.cjf,2))
		self.StrVar42.set(self.trans_vars(self.djf,2))
		self.StrVar43.set(self.trans_vars(self.hag_core_mod_psin,2))
		self.StrVar44.set(self.trans_vars(self.core_neo,2))	
		self.StrVar45.set(self.trans_vars(self.currenti,2))
		self.StrVar46.set(self.trans_vars(self.relax,2))
		self.StrVar47.set(self.trans_vars(self.ns,2))	
		self.StrVar48.set(self.trans_vars(self.nt,2))
		self.StrVar49.set(self.trans_vars(self.map_ns,2))
		self.StrVar50.set(self.trans_vars(self.map_nt,2))

		self.StrVar51.set(self.trans_vars(self.mis_gridn,2))
		self.StrVar52.set(self.trans_vars(self.mis_psistart,2))
		self.StrVar53.set(self.trans_vars(self.xr1,2))	
		self.StrVar54.set(self.trans_vars(self.sig1,2))
		self.StrVar55.set(self.trans_vars(self.xr2,2))
		self.StrVar56.set(self.trans_vars(self.sig2,2))	
		self.StrVar57.set(self.trans_vars(self.eli_gridn,2))
		self.StrVar58.set(self.trans_vars(self.eli_ndist,2))		
		self.StrVar59.set(self.trans_vars(self.eli_psistart,2))	
		self.StrVar60.set(self.trans_vars(self.bsmulti,1))
		self.StrVar61.set(self.trans_vars(self.run_name,1))
		self.StrVar62.set(self.outdir)
		self.StrVar63.set(self.trans_vars(self.ncrit,2))

		self.StrVar64.set(self.trans_vars(self.epslon,2))
		self.StrVar65.set(self.trans_vars(self.beta_crit,2))
		self.StrVar66.set(self.trans_vars(self.li_crit,2))
		self.StrVar67.set(self.trans_vars(self.ip_crit,2))
		self.StrVar68.set(self.trans_vars(self.bs_crit,2))
		self.StrVar69.set(self.trans_vars(self.niterl,2))

		self.StrVar70.set(self.trans_vars(self.qdelfix,2))

		return

	def trans_vars(self,var1,type):

		if (type == 4):
			if (var1):
				return 1
			else:
				return 0
		elif (type == 2):
			return str(var1)
		
		elif (type == 3):
			if (var1 == None):
				return ''
			else:
				return var1
		elif (type == 1):
			return var1
		elif (type == 5):		
			return self.var_list[var1-1]
		elif (type == 6):
			ans = ''
			if not (len(var1) == 0):
				for i in range(len(var1)-1):
					ans = ans +str(var1[i])+','
				ans = ans + str(var1[-1])
			return ans

	def read_preset(self,type=1):

		try:	f=open('gfitp_history.dat','r')
		except:	return

		if (type==1):
			while True:
				line = f.readline()
				if not line: break
				line = line.split()
				if line[0].lower() == 'eq_file':
					self.eqdsk_name = line[1]
					self.use_param_shape = False
					self.fix_core = True
				elif line[0].lower() == 'bnd_psi':
					self.target_bnd = float(line[1])
				elif line[0].lower() == 'run_name':
					self.run_name = line[1]
					self.outdir = line[1]					
				elif line[0].lower() == 'zeff':
					self.zeff = float(line[1])
				elif line[0].lower() == 'zimp':
					self.zimp = float(line[1])
				elif line[0].lower() == 'aimp':
					self.aimp = float(line[1])
				elif line[0].lower() == 'amain':
					self.amain = float(line[1])
				elif line[0].lower() == 'bs_model':
					self.bsmodel = line[1].lower()
				elif line[0].lower() == 'bsmulti':
					self.bsmulti = float(line[1])					

			f.close()
			return
		else:
			while True:
				line = f.readline()
				if not line: break
				line = line.split()
				if line[0].lower() == 'eq_file':
					self.e2.delete(0, 'end')
					self.e2.insert(0, line[1])
					self.button_func1(True,line[1])
					self.use_param_shape = False
				elif line[0].lower() == 'bnd_psi':
					self.DoubleVar1.set(float(line[1]))
				elif line[0].lower() == 'run_name':
					self.e62.delete(0, 'end')
					self.e62.insert(0, line[1])				

			f.close()
		return		
		
	def read_eped_opt(self,filename):

		f = open(filename,'r')

		while True:
			line = f.readline()
			if not line:	break

			beta_val = 1.;	run_mode = 'normal';	use_neo = False;	use_hager = False;	use_chang = True;
			beta_criterion = 1;

			self.eqdsk_name = read_namelist_str(line,'EQDSK',self.eqdsk_name,3)
			self.use_param_shape = read_namelist_str(line,'USE_PARAM_SHAPE',self.use_param_shape,4)
			self.ip = read_namelist_str(line,'IP',self.ip,2)
			self.bt = read_namelist_str(line,'BVAC',self.bt,2)
			self.zimp = read_namelist_str(line,'ZIMP',self.zimp,2)
			self.zeff = read_namelist_str(line,'ZEFF',self.zeff,2)
			self.amain = read_namelist_str(line,'AMAIN',self.amain,2)
			self.aimp = read_namelist_str(line,'AIMP',self.aimp,2)
			beta_val = read_namelist_str(line,'Beta_val',beta_val,2)
			self.li = read_namelist_str(line,'LITARGET',self.li,2)		 
			self.apf = read_namelist_str(line,'APF',self.apf,2)
			self.bpf = read_namelist_str(line,'BPF',self.bpf,2)
			self.cpf = read_namelist_str(line,'CPF',self.cpf,2)
			self.dpf = read_namelist_str(line,'DPF',self.dpf,2)		 
			use_neo = read_namelist_str(line,'USE_NEO',use_neo,4)
			use_hager = read_namelist_str(line,'USE_HAGER',use_hager,4)
			use_chang = read_namelist_str(line,'USE_CHANG',use_chang,4)
			self.hag_core_mod = read_namelist_str(line,'HAG_CORE_MOD',self.hag_core_mod,4)
			self.hag_core_mod_psin = read_namelist_str(line,'HAG_CORE_MOD_PSIN',self.hag_core_mod_psin,2)
			self.core_neo = read_namelist_str(line,'Core_neo',self.core_neo,2)
			self.bsmulti = read_namelist_str(line,'BSMULTI',self.bsmulti,2)
			self.ajf = read_namelist_str(line,'AJF',self.ajf,2)
			self.bjf = read_namelist_str(line,'BJF',self.bjf,2)
			self.cjf = read_namelist_str(line,'CJF',self.cjf,2)
			self.djf = read_namelist_str(line,'DJF',self.djf,2)
			self.currenti = read_namelist_str(line,'Current_ITERN',self.currenti,1)
			self.relax = read_namelist_str(line,'RELAX',self.relax,2)
			 
			self.kbmc1 = read_namelist_str(line,'KBM_COEF',self.kbmc1,2)
			self.kbmc2 = read_namelist_str(line,'KBM_POWER',self.kbmc2,2)
			self.neped = read_namelist_str(line,'NEPED',self.neped,2)
			self.nsep = read_namelist_str(line,'NESEP',self.nsep,2)
			self.tsep = read_namelist_str(line,'TESEP',self.tsep,2)
			self.wmin = read_namelist_str(line,'W_MIN',self.wmin,2)
			self.wmax = read_namelist_str(line,'W_MAX',self.wmax,2)
			self.wstep = read_namelist_str(line,'W_STEP',self.wstep,1)
			self.at1 = read_namelist_str(line,'AT1',self.at1,2)
			self.an1 = read_namelist_str(line,'AN1',self.an1,2)
			self.alpt1 = read_namelist_str(line,'ALPT1',self.alpt1,2)
			self.alpt2 = read_namelist_str(line,'ALPT2',self.alpt2,2)
			self.alpn1 = read_namelist_str(line,'ALPN1',self.alpn1,2)
			self.alpn2 = read_namelist_str(line,'ALPN2',self.alpn2,2)

			self.mode_n = read_namelist_str(line,'MODEN',self.mode_n,3)
			self.stab_code = read_namelist_str(line,'RUN_STAB',self.stab_code,3)
			self.ns = read_namelist_str(line,'NS',self.ns,1)
			self.nt = read_namelist_str(line,'NT',self.nt,1)
			self.map_ns = read_namelist_str(line,'MAP_NS',self.map_ns,1)
			self.map_nt = read_namelist_str(line,'MAP_NT',self.map_nt,1)
			self.mis_gridn = read_namelist_str(line,'MIS_GRIDN',self.mis_gridn,1)
			self.mis_psistart = read_namelist_str(line,'MIS_PSISTART',self.mis_psistart,2)
			self.xr1 = read_namelist_str(line,'MIS_XR1',self.xr1,2)
			self.sig1 = read_namelist_str(line,'MIS_SIG1',self.sig1,2)
			self.xr2 = read_namelist_str(line,'MIS_XR2',self.xr2,2)
			self.sig2 = read_namelist_str(line,'MIS_SIG2',self.sig2,2)
			self.eli_comp = read_namelist_str(line,'ELI_COMP',self.eli_comp,4)
			self.eli_gridn = read_namelist_str(line,'ELI_GRIDN',self.eli_gridn,1)
			self.eli_psistart = read_namelist_str(line,'ELI_PSISTART',self.eli_psistart,2)
			self.eli_ndist = read_namelist_str(line,'ELI_NDIST',self.eli_ndist,1)
			self.nqa = read_namelist_str(line,'NQA',self.nqa,2)

			self.epslon = read_namelist_str(line,'EPSILON',self.epslon,2)
			self.beta_crit = read_namelist_str(line,'Beta_crit',self.beta_crit,2)
			self.li_crit = read_namelist_str(line,'Li_crit',self.li_crit,2)
			self.ip_crit = read_namelist_str(line,'Ip_crit',self.ip_crit,2)
			self.bs_crit = read_namelist_str(line,'Bs_crit',self.bs_crit,2)
			self.niterl = read_namelist_str(line,'Li_ITERN',self.niterl,1)

			run_mode = read_namelist_str(line,'RUN_MODE',run_mode,3)
			beta_criterion = read_namelist_str(line,'Beta_criterion_type',beta_criterion,1)
			self.teped = read_namelist_str(line,'TEPED',self.teped,2)
			self.tewidth = read_namelist_str(line,'TEWIDTH',self.tewidth,2)
			self.qdelfix = read_namelist_str(line,'QDELFIX',self.qdelfix,2)
			self.use_raw_prof = read_namelist_str(line,'USERAW',self.use_raw_prof,4)

		f.close()

		if beta_criterion == 1:	
			self.constraint = 'bpol'
			self.bp = beta_val;
		elif beta_criterion == 2:	
			self.constraint = 'rmag'
			self.rmag = beta_val;
		else:	
			self.constraint = 'betan'
			self.betan = beta_val

		if use_neo:	self.bsmodel = 'neo'
		elif use_hager:	self.bsmodel = 'hager'
		elif use_chang:	self.bsmodel = 'csauter'
		else:	self.bsmodel = 'sauter'

		run_mode = run_mode.lower()

		if run_mode == 'normal':
			self.use_li = False
			self.use_li_iter2 = False
			self.fix_core = True
		elif run_mode == 'eped':
			self.use_li = False
			self.use_li_iter2 = False
			self.fix_core = False
		elif run_mode == 'eped2':
			self.use_li = True
			self.use_li_iter2 = False
			self.fix_core = True
		else:
			self.use_li = True
			self.use_li_iter2 = True
			self.fix_core = True
		return

	def __init__(self):

		self.currdir = os.getcwd()
		self.tfc_init = False

		return

if __name__ == "__main__":

	import gui_eped

	print(' ---------------------------------------------------------------')
	print('||              Python based EPED tool Ver  %s                ||'%version['eped'])
	print('||                  Edge pedestal predictor                    ||')
	print('%s'%author['eped'])
	print(' ---------------------------------------------------------------\n')

	gped = gui_eped.gped()	
	gped.root = tk.Tk()
	gped.root.title('G-EPED')
	gped.declare_vars()
	gped.initialise_vars()
	gped.run_list = []
	gped.run_list = gped.open_result_list()	
	gped.gui_eped()
