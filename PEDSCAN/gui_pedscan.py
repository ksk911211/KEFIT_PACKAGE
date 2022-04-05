#!/usr/local/anaconda3/bin/python3
import os, sys
import tkinter as tk
import numpy as np
from shutil import move, copyfile, copytree
from scipy.interpolate import interp1d
from tkinter.filedialog import askopenfilename,asksaveasfilename, askdirectory
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import matplotlib.pyplot as plt
import eqdsk
import subprocess
from nodelist import make_nodelist
import pedscan
from read_namelist import read_namelist_str
from batch_run import *
from exec_dirs import gfit_dir,pedscane_dir,node_force,node_init,node_default,node_machine,version,author

currdir0 = os.getcwd()
class gpedscan:

	def declare_vars(self):
		#Vars
		for i in range(1,40):
			self.__dict__['CheckVar%d'%i] = tk.IntVar()

		for i in range(1,80):
			self.__dict__['StrVar%d'%i] = tk.StringVar()

		for i in range(1,10):
			self.__dict__['DoubleVar%d'%i] = tk.DoubleVar()

		for i in range(1,40):
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

		ax1.text(self.eq.rmag,self.eq.zmag,inform,color='lime')

		plt.draw()	

		return

	def checkvar_func1(self):

		if self.CheckVar2.get() == 1:	self.CheckVar1.set(1)


		if self.CheckVar1.get() == 1:
			self.fixed_ti = True
		else:
			self.fixed_ti = False

		return

	def checkvar_func2(self):

		if not self.iskinti:	self.CheckVar2.set(0)

		if self.CheckVar2.get() == 1:	self.CheckVar1.set(1)
		else:	
			if not self.fixed_ti:	self.CheckVar1.set(0)
		return

	def detect_close(self,win_num):

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
			os.chdir(self.currdir)
			self.t5.destroy()
			self.t5_close = True

		elif (win_num == 6):
			self.t6.destroy()
			self.t6_close = True						
		elif (win_num == 7):
			self.t7.destroy()
			self.t7_close = True
		elif (win_num == 8):
			self.t8.destroy()
			self.t8_close = True							

		return		

	def scroll_func1(self):

		self.target_bnd	= self.DoubleVar1.get()

		if (self.iseqdsk):

			self.eq.target_psin = self.target_bnd
			self.eq.get_target_bnd()
		try:
			self.rzbdy = np.copy(self.eq.rzbdyt)
		except:
			print('>>> Have to load equilibrium first')
		self.update_psirz()

		return

	def scroll_func2(self,sim):

		sim.target_i2 = int(self.DoubleVar4.get())
		sim.target_j2 = int(self.DoubleVar5.get())

		self.fig4.canvas.draw_idle()
		if self.CheckVar14.get() == 1:
			sim.draw_plot3(self.fig4,self.DoubleVar4.get(),self.DoubleVar5.get(),True)
		else:
			sim.draw_plot3(self.fig4,self.DoubleVar4.get(),self.DoubleVar5.get(),False)

		self.button_func16(sim)	

		return

	def renew_kin_file(self):

		dirs = 'PROFILES/chease_eped_mod'
		if (os.path.isfile(dirs)):
			self.e2.delete(0,'end')
			self.e2.insert(10,dirs)
		self.read_kin_file(dirs)
		dirs = 'PROFILES/chease_kinprof_fit'
		if (os.path.isfile(dirs)):
			self.e3.delete(0,'end')
			self.e3.insert(10,dirs)		
		return

	def read_kin_file(self,dirs):

		if (os.path.isfile(dirs)):
			try:
				f = open(dirs,'r')
				for i in range(3):
					line = f.readline().split()
					for j in range(8):
						self.eped_prof[i,j] = float(line[j])
				line = f.readline().split()
				self.zeff  = float(line[0])
				self.zimp  = float(line[1])
				self.amain = float(line[2])
				self.aimp  = float(line[3])
				self.e18.delete(0,'end')
				self.e18.insert(10,str(round(self.zeff,3)))
				self.e19.delete(0,'end')
				self.e19.insert(10,str(round(self.zimp,3)))
				self.e20.delete(0,'end')
				self.e20.insert(10,str(round(self.amain,3)))
				self.e21.delete(0,'end')
				self.e21.insert(10,str(round(self.aimp,3)))									
				self.iskin = True
			except:
				self.iskin = False
		else:
			self.iskin = False
		return

	def make_input_files(self):

		rundir = self.currdir + '/' + self.e65.get()
		try:
			os.mkdir(rundir+'/input')
		except:
			pass
		currdir = os.getcwd()
		inputdir = rundir + '/input'	

		if self.iseqdsk:
			name = self.e1.get()
			name2 = name.split('/')[-1]
			copyfile(name,inputdir+'/'+name2.split()[0])

		if self.iskin:
			name = self.e2.get()
			name2 = name.split('/')[-1]
			copyfile(name,inputdir+'/'+name2.split()[0])

		if self.iskinti:
			name = self.e3.get()
			name2 = name.split('/')[-1]
			copyfile(name,inputdir+'/'+name2.split()[0])			

		if self.isrot:
			name = self.e55.get()
			name2 = name.split('/')[-1]
			copyfile(name,inputdir+'/'+name2.split()[0])

		return	

	def write_pedscan_input(self,filename='ja_namelist'):

		f = open(filename,'w')
		f.write('--- Run option ---\n')
		f.write('run_name = %s \n'%self.e65.get())
		if self.CheckVar3.get() == 1:	f.write('USE_li = True\n')
		else:	f.write('USE_li = False \n')
		if self.CheckVar6.get() == 1:	f.write('USE_li2 = True\n')
		else:	f.write('USE_li2 = False \n')
		if self.CheckVar4.get() == 1:	f.write('Adjust_prof = True\n')
		else:	f.write('Adjust_prof = False \n')		
		if self.CheckVar5.get() == 1:	f.write('HAG_CORE_MOD = True\n')
		else:	f.write('HAG_CORE_MOD = False \n')
		if self.CheckVar7.get() == 1:	f.write('USE_COMP = True\n')
		else:	f.write('USE_COMP = False \n')
		if self.CheckVar8.get() == 1:	f.write('USE_ROT = True\n')
		else:	f.write('USE_ROT = False \n')		
		if self.CheckVar1.get() == 1:	f.write('FIX_TI = True\n')
		else:	f.write('FIX_TI = False \n')		
		if self.CheckVar2.get() == 1:	f.write('EXT_TI = True\n')
		else:	f.write('EXT_TI = False \n')


		f.write('!-- bootstrap option\n')
		if (self.StrVar13.get().lower() == 'nhager' or self.StrVar13.get().lower() == 'neo'):
			f.write('use_neo = True \n')
		else:
			f.write('use_neo = False \n')
		if self.StrVar13.get().lower() == 'hager':	
			f.write('use_hager = True \n')
		else:
			f.write('use_hager = False \n')
		if self.StrVar13.get().lower() == 'csauter':
			f.write('use_chang = True \n')
		else:
			f.write('use_chang = False \n')
		f.write('\n')
		f.write('!-- equil option\n')
		f.write('BND_PSIN = %f\n'%self.target_bnd)
		f.write('EQDSK = input/%s \n'%self.e1.get().split('/')[-1])
		f.write('EPED_PROF_FILE = input/%s \n'%self.e2.get().split('/')[-1])
		f.write('TI_FILE = input/%s\n'%self.e3.get().split('/')[-1])
		f.write('scan_nw = %s\n'%self.e4.get())
		f.write('scan_nh = %s\n'%self.e7.get())
		f.write('wmin = %s\n'%self.e5.get())
		f.write('wmax = %s\n'%self.e6.get())
		f.write('hmin = %s\n'%self.e8.get())
		f.write('hmax = %s\n'%self.e9.get())
		f.write('\n')
		f.write('stab_code = %s\n'%self.StrVar11.get())
		f.write('moden = %s\n'%self.e12.get())
		f.write('\n')
		if (self.StrVar14.get().lower() == 'bpol'):
			f.write('Beta_criterion_type = 1\n')
			f.write('Beta_val = %s\n'%self.e15.get())
		else:
			f.write('Beta_criterion_type = 2\n')
			f.write('Beta_val = %s\n'%self.e16.get())
		f.write('BPOL = %s\n'%self.e15.get())
		f.write('RMAG = %s\n'%self.e16.get())
		f.write('LI = %s\n'%self.e17.get())
		f.write('ZEFF = %s\n'%self.e18.get())
		f.write('ZIMP = %s\n'%self.e19.get())
		f.write('AMAIN = %s\n'%self.e20.get())
		f.write('AIMP = %s\n'%self.e21.get())
		f.write('LINEDEN = %s\n'%self.e22.get())
		f.write('BSMULTI = %s\n'%self.e23.get())
		f.write('\n')

		f.write('APF = %s\n'%self.StrVar25.get())
		f.write('BPF = %s\n'%self.StrVar26.get())
		f.write('CPF = %s\n'%self.StrVar27.get())
		f.write('DPF = %s\n'%self.StrVar28.get())
		f.write('AJF = %s\n'%self.StrVar29.get())
		f.write('BJF = %s\n'%self.StrVar30.get())
		f.write('CJF = %s\n'%self.StrVar31.get())
		f.write('DJF = %s\n'%self.StrVar32.get())
		f.write('hag_core_mod_psin = %s\n'%self.StrVar33.get())
		f.write('Current_ITERN = %s\n'%self.StrVar34.get())
		f.write('RELAX = %s\n'%self.StrVar35.get())
		f.write('\n')

		f.write('NS = %s\n'%self.StrVar36.get())
		f.write('NT = %s\n'%self.StrVar37.get())
		f.write('MAP_NS = %s\n'%self.StrVar38.get())
		f.write('MAP_NT = %s\n'%self.StrVar39.get())
		f.write('NSE = %s\n'%self.StrVar40.get())
		f.write('NTE = %s\n'%self.StrVar41.get())
		f.write('MAP_NSE = %s\n'%self.StrVar42.get())
		f.write('MAP_NTE = %s\n'%self.StrVar43.get())
		f.write('\n')

		f.write('MIS_GRIDN = %s\n'%self.StrVar44.get())
		f.write('MIS_PSISTART = %s\n'%self.StrVar45.get())
		f.write('XR1 = %s\n'%self.StrVar46.get())
		f.write('SIG1 = %s\n'%self.StrVar47.get())
		f.write('XR2 = %s\n'%self.StrVar48.get())
		f.write('SIG2 = %s\n'%self.StrVar49.get())
		f.write('\n')

		f.write('delQ = %s\n'%self.e10.get())
		f.write('ELI_GRIDN = %s\n'%self.StrVar50.get())
		f.write('ELI_NDIST = %s \n'%self.StrVar51.get())
		f.write('ELI_PSISTART = %s\n'%self.StrVar52.get())
		f.write('ROT_FILE = %s\n'%self.e64.get())

		f.write('EPSLON = %s\n'%self.StrVar53.get())
		f.write('Beta_crit = %s\n'%self.StrVar54.get())
		f.write('li_crit = %s\n'%self.StrVar55.get())
		f.write('ip_crit = %s\n'%self.StrVar56.get())
		f.write('bs_crit = %s\n'%self.StrVar57.get())
		f.write('NITERL = %s\n'%self.StrVar58.get())
		f.write('\n')

		f.write('NQA = %s\n'%self.StrVar59.get())
		f.write('GRCRIT1 = %s\n'%self.StrVar61.get())
		f.write('GRCRIT2 = %s\n'%self.StrVar60.get())
		f.write('NCUT = %s\n'%self.StrVar62.get())
		f.write('GRMAX = %s\n'%self.StrVar63.get())

		if self.CheckVar11.get() ==1:	f.write('USE_BILINEAR = True\n')
		else:	f.write('USE_BILINEAR = False\n')
		if self.CheckVar12.get() ==1:	f.write('USE_DIA = True\n')
		else:	f.write('USE_DIA = False\n')
		f.write('TARGETI = %f\n'%self.DoubleVar2.get())
		f.write('TARGETJ = %f\n'%self.DoubleVar3.get())
		if self.CheckVar9.get() ==1:	f.write('HIGHQ = True\n')
		else:	f.write('HIGHQ = False\n')		
			
		f.write('TI_FILE_TYPE = 2\n')
		if node_machine.lower() == 'fusma':	self.node_list = make_nodelist(self.nodelist)
		if node_force:	self.node_list = node_default
		f.write('node_list = %s\n'%self.node_list)
		f.write('TARGETN = 10\n')
		#### NODELIST
		f.write('\n')
		f.close()

		return

	def modify_pedscan_input(self,sim,filename='scan_opt'):

		f = open(filename,'r')
		f2 = open(filename+'_temp','w')
		nqa = 0; gr1 = 0; gr2 = 0; ncut = 0; grm = 0; bili = 0; dia = 0; ti = 0; tj = 0; tn = 0; 
		while True:
			line = f.readline()
			if not line:	break
			if line.lower().find('nqa') > -1:
				line2 = 'NQA = %f\n'%sim.nqa; nqa = 1;
			elif line.lower().find('grcrit1') > -1:
				line2 = 'GRCRIT1 = %f\n'%sim.gr_crit1; gr1 = 1;
			elif line.lower().find('grcrit2') > -1:
				line2 = 'GRCRIT2 = %f\n'%sim.gr_crit2; gr2 = 1;
			elif line.lower().find('ncut') > -1:
				line2 = 'NCUT = %i\n'%sim.ncut; ncut = 1;
			elif line.lower().find('grmax') > -1:
				line2 = 'GRMAX = %f\n'%sim.grmax; grm = 1;
			elif line.lower().find('bilinear') > -1:
				line2 = 'USE_BILINEAR = %s\n'%sim.use_bilinear; bili = 1;
			elif line.lower().find('use_dia') > -1:
				line2 = 'USE_DIA = %s\n'%sim.use_dia;	dia = 1;
			elif line.lower().find('targeti') > -1:
				line2 = 'TARGETI = %f\n'%sim.target_i;	ti = 1;
			elif line.lower().find('targetj') > -1:
				line2 = 'TARGETJ = %f\n'%sim.target_j;	tj = 1;
			elif line.lower().find('targetn') > -1:
				line2 = 'TARGETN = %i\n'%sim.targetn;	tn = 1;
			else:
				line2 = line
			f2.write(line2)
		if nqa == 0:
			line2 = 'NQA = %f\n'%sim.nqa;
			f2.write(line2)
		if gr1 == 0:
			line2 = 'GRCRIT1 = %f\n'%sim.gr_crit1;
			f2.write(line2)
		if gr2 == 0:
			line2 = 'GRCRIT2 = %f\n'%sim.gr_crit2; 
			f2.write(line2)
		if ncut == 0:
			line2 = 'NCUT = %i\n'%sim.ncut;
			f2.write(line2)
		if grm == 0:
			line2 = 'GRMAX = %f\n'%sim.grmax;
			f2.write(line2)
		if bili == 0:
			line2 = 'USE_BILINEAR = %s\n'%sim.use_bilinear;
			f2.write(line2)
		if dia == 0:
			line2 = 'USE_DIA = %s\n'%sim.use_dia;
			f2.write(line2)
		if ti == 0:
			line2 = 'TARGETI = %f\n'%sim.target_i;
			f2.write(line2)
		if tj == 0:
			line2 = 'TARGETJ = %f\n'%sim.target_j;
			f2.write(line2)
		if tn == 0:
			line2 = 'TARGETN = %i\n'%sim.targetn;
			f2.write(line2)

		f2.close()

		os.remove(filename)
		move(filename+'_temp',filename)

		return

	def open_result_list(self):

		run_list =['-search-']
		filename = os.getcwd()+'/pedscan.log'
		if os.path.isfile(filename):
			f = open(filename,'r')
			while True:
				line = f.readline()
				if not line: break
				line = line.split('\n')
				run_list.append(line[0])
			f.close()

		return run_list

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
		self.e1.delete(0, 'end')
		self.e1.insert(10,self.input)

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

		self.eq.epslon = 4
		self.eq.make_chease_input()
		self.eq.run_chease()
		os.remove('chease_namelist')
		os.remove('EXPEQ')
		self.bp = round(self.eq.bp,3)
		self.li = round(self.eq.li,3)
		self.ip = round(abs(self.eq.ip) / 1.e6,3)
		self.bt = round(abs(self.eq.bcentr),3)

		self.e15.delete(0,'end')
		self.e15.insert(10,str(self.bp))
		self.e16.delete(0,'end')
		self.e16.insert(10,str(round(self.eq.rmag,3)))
		self.e17.delete(0,'end')
		self.e17.insert(10,str(self.li))

		return		

	def button_func2(self):

		if not self.iseqdsk:
			print('>> EQDSK must be selected first...')
			return

		if not self.t1_close:
			print('>>> NUM Param Window is already opened...')
			return			

		self.t1 = tk.Toplevel(self.root)
		self.t1.wm_title("EFIT Inform")
		self.t1_close = False		
		self.t1.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(1))

		self.t1.resizable(0,0)
		self.fig2, [ax1, ax2] = plt.subplots(1,2,figsize=(12,5))

		#self.fig2.tight_layout()

		self.canvas2 = FigureCanvasTkAgg(self.fig2,master=self.t1)
		self.plot_widget2 = self.canvas2.get_tk_widget()
		self.plot_widget2.grid(rowspan=40,row=2,column=0,columnspan=40)

		toolbar_frame = tk.Frame(self.t1)
		toolbar_frame.grid(column=0,row=0)
		
		self.toolbar2 = NavigationToolbar2Tk(self.canvas2,toolbar_frame)
		self.fig2.canvas.draw_idle()

		ax1.plot(self.eq.psin,-self.eq.pp/max(abs(self.eq.pp)),self.eq.psin,self.eq.ffp/max(abs(self.eq.ffp)))
		ax1.set_title("P'($\psi$), FF'($\psi$)")
		ax1.set_xlabel('Normalized radius ($\psi_N$)')
		ax1.set_ylabel('Normalized value [a.u]')
		ax1.legend(["P'($\psi$)","FF'($\psi$)"])
		ax1.set_xlim((-0.05,1.05))

		jplowf = np.copy(self.eq.psin)
		jphighf = np.copy(self.eq.psin)

		Rlow = interp1d(self.eq.prhoR[:,0],self.eq.prhoR[:,2],'slinear')
		Rhigh  = interp1d(self.eq.prhoR[:,0],self.eq.prhoR[:,3],'slinear')

		mu0 = np.pi * 4.0 * 1.e-7

		for i in range(len(self.eq.psin)):
			jplowf[i] = - Rlow(self.eq.psin[i]) * self.eq.pp[i]  - self.eq.ffp[i]/mu0/Rlow(self.eq.psin[i])
			jphighf[i] = -Rhigh(self.eq.psin[i]) * self.eq.pp[i] - self.eq.ffp[i]/mu0/Rhigh(self.eq.psin[i])

		ax2.plot(self.eq.psin,jplowf/1.e6,self.eq.psin,jphighf/1.e6,self.eq.javg[:,0],self.eq.javg[:,1]/self.eq.javg[0,1]*jphighf[0]/1.e6)
		ax2.set_title("$J_{\phi}$($\psi$)")
		ax2.set_xlabel('Normalized radius ($\psi_N$)')
		ax2.set_ylabel('Current density [$MA/m^2$]')
		ax2.legend(["LFS","HFS","AVERAGE"])
		ax2.set_xlim((-0.05,1.05))

		return

	def button_func3(self,skip=False,file=None):

		if not skip:
			self.input2 = askopenfilename()
			if (len(self.input2) == 0):
				return
		else:
			self.input2 = file

		self.e2.delete(0,'end')
		self.e2.insert(10,self.input2)
		#self.CheckVar1.set(1)
		self.read_kin_file(self.input2)

		return

	def button_func4(self):

		os.system(gfit_dir)
		self.renew_kin_file()

		return

	def button_func5(self,skip=False,file=None):

		if not skip:
			self.input2 = askopenfilename()
			if (len(self.input2) == 0):
				return
		else:
			self.input2 = file
		self.e3.delete(0,'end')
		self.e3.insert(10,self.input2)
		self.ti_file = self.input2
		#self.CheckVar1.set(1)
		self.iskinti = True
		self.CheckVar2.set(1)

		return			

	def button_func6(self):

		if not self.t7_close:
			print('>>> Kinetic profile Window is already opened...')
			return

		psin = np.linspace(0,1,201)
		prof1 = np.zeros((201,3))
		
		if not (self.iskin or self.iskinti):	
			print('>>> No available kinetic profile...')
			return
		if self.iskinti:
			try:
				f = open(self.ti_file,'r')
				line = f.readline()
				dat_num = int(float(line))
				prof2 = np.zeros((dat_num,2))
				line = f.readline()

				for i in range(dat_num):	
					line = f.readline().split()
					prof2[i,0] = float(line[0])
					prof2[i,1] = float(line[3])
				f.close()
			except:
				print('>>> Kinetic profile is invalid')
				self.iskinti = False

		if self.iskin:
			for i in range(201):
				psi = psin[i]
				for j in range(3):
					tanh1 = np.tanh(2.*(self.eped_prof[j,3]-self.eped_prof[j,5])/self.eped_prof[j,4])
					tanh2 = np.tanh(2.*(self.eped_prof[j,3]-psi)/self.eped_prof[j,4])
					core1 = 1. - (psi/self.eped_prof[j,5])**self.eped_prof[j,6]
					if psi > self.eped_prof[j,5]: core1 = 0.
					
					core2 = core1 ** self.eped_prof[j,7]
					val = self.eped_prof[j,1]
					val = val + self.eped_prof[j,0] * (tanh1+tanh2)
					val = val + self.eped_prof[j,2] * core2

					prof1[i,j] = val

		if not (self.iskin or self.iskinti):	
			print('>>> No available kinetic profile...')
			return

		self.t7 = tk.Toplevel(self.root)
		self.t7.wm_title("Kinetic profile")
		self.t7_close = False		
		self.t7.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(7))

		self.t7.resizable(0,0)
		self.fig5, [ax1, ax2, ax3] = plt.subplots(1,3,figsize=(13,5))

		self.canvas5 = FigureCanvasTkAgg(self.fig5,master=self.t7)
		self.plot_widget5 = self.canvas5.get_tk_widget()
		self.plot_widget5.grid(rowspan=40,row=2,column=0,columnspan=40)

		toolbar_frame = tk.Frame(self.t7)
		toolbar_frame.grid(column=0,row=0)
		
		self.toolbar5 = NavigationToolbar2Tk(self.canvas5,toolbar_frame)
		self.fig5.canvas.draw_idle()


		
		ax1.plot(psin,prof1[:,1])
		ax1.set_title("TE[keV]")
		ax1.set_xlabel('Normalized radius ($\psi_N$)')
		ax1.set_ylabel('TE[keV]')
			
		ax1.set_xlim((-0.05,1.05))

		ax2.plot(psin,prof1[:,0])
		ax2.set_title("NE[10(19)/m3]")
		ax2.set_xlabel('Normalized radius ($\psi_N$)')
		ax2.set_ylabel('TE[10(19)/m3]')
		ax2.set_xlim((-0.05,1.05))

		llegend = [];	plegend = [];
		if self.iskin:
			line, = ax3.plot(psin,prof1[:,2])
			plegend.append('EPED_PROF_TI')
			llegend.append(line)
			line, = ax3.plot(psin,prof1[:,1],linestyle='--')
			plegend.append('EPED_PROF_TE')
			llegend.append(line)

		if self.iskinti: 
			line, = ax3.plot(prof2[:,0],prof2[:,1])
			plegend.append('EXT_PROF_TI')
			llegend.append(line)

		ax3.set_title("TI[keV]")
		ax3.set_xlabel('Normalized radius ($\psi_N$)')
		ax3.set_ylabel('TI[keV]')
		ax3.set_xlim((-0.05,1.05))
		ax3.legend(llegend,plegend)			

		self.fig5.tight_layout()

		return

	def button_func7(self):

		if not self.t1_close:
			print('>>> NUM Param Window is already opened...')
			return

		self.t2 = tk.Toplevel(self.root)
		self.t2.wm_title("NUM Params")
		self.t2_close = False		
		self.t2.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(2))


		self.l2 = tk.Label(self.t2, text = '====== Beam profile ======',justify='center')
		self.l2.grid(row=1,column=0,columnspan=4)

		self.l2 = tk.Label(self.t2, text = 'APF ',justify='center')
		self.l2.grid(row=2,column=0,columnspan=2)
		self.e25 = tk.Entry(self.t2,width=5,justify='center')
		self.e25.insert(0,self.StrVar25.get())
		self.e25.grid(row=2, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'BPF [0-1]',justify='center')
		self.l2.grid(row=3,column=0,columnspan=2)
		self.e26 =tk.Entry(self.t2,width=5,justify='center')
		self.e26.insert(0,self.StrVar26.get())
		self.e26.grid(row=3, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'CPF ',justify='center')
		self.l2.grid(row=4,column=0,columnspan=2)
		self.e27 = tk.Entry(self.t2,width=5,justify='center')
		self.e27.insert(0,self.StrVar27.get())
		self.e27.grid(row=4, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'DPF ',justify='center')
		self.l2.grid(row=5,column=0,columnspan=2)
		self.e28 = tk.Entry(self.t2,width=5,justify='center')
		self.e28.insert(0,self.StrVar28.get())
		self.e28.grid(row=5, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t2, text = 'AJF ',justify='center')
		self.l2.grid(row=6,column=0,columnspan=2)
		self.e29 = tk.Entry(self.t2,width=5,justify='center')
		self.e29.insert(0,self.StrVar29.get())
		self.e29.grid(row=6, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t2, text = 'BJF [0-1]',justify='center')
		self.l2.grid(row=7,column=0,columnspan=2)
		self.e30 = tk.Entry(self.t2,width=5,justify='center')
		self.e30.insert(0,self.StrVar30.get())
		self.e30.grid(row=7, column=2, columnspan = 2)						

		self.l2 = tk.Label(self.t2, text = 'CJF ',justify='center')
		self.l2.grid(row=8,column=0,columnspan=2)
		self.e31 = tk.Entry(self.t2,width=5,justify='center')
		self.e31.insert(0,self.StrVar31.get())
		self.e31.grid(row=8, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'DJF ',justify='center')
		self.l2.grid(row=9,column=0,columnspan=2)
		self.e32 = tk.Entry(self.t2,width=5,justify='center')
		self.e32.insert(0,self.StrVar32.get())
		self.e32.grid(row=9, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = '====== Current option ======',justify='center')
		self.l2.grid(row=10,column=0,columnspan=4)		

		self.l2 = tk.Label(self.t2, text = 'HAGMOD ',justify='center')
		self.l2.grid(row=11,column=0,columnspan=2)
		self.c6 = tk.Checkbutton(self.t2,variable=self.CheckVar5)
		self.c6.grid(row=11, column=2, columnspan = 2)
		self.l2 = tk.Label(self.t2, text = 'MODPSIN [0-1]',justify='center')
		self.l2.grid(row=12,column=0,columnspan=2)
		self.e33 = tk.Entry(self.t2,width=5,justify='center')
		self.e33.insert(0,self.StrVar33.get())
		self.e33.grid(row=12, column=2, columnspan = 2)
		self.l2 = tk.Label(self.t2, text = 'LI_ITER2 ',justify='center')
		self.l2.grid(row=13,column=0,columnspan=2)
		self.c6 = tk.Checkbutton(self.t2,variable=self.CheckVar6)
		self.c6.grid(row=13, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t2, text = '====== Numerical option ======',justify='center')
		self.l2.grid(row=16,column=0,columnspan=4)		

		self.l2 = tk.Label(self.t2, text = 'ITERN [#]',justify='center')
		self.l2.grid(row=17,column=0,columnspan=2)
		self.e34 = tk.Entry(self.t2,width=5,justify='center')
		self.e34.insert(0,self.StrVar34.get())
		self.e34.grid(row=17, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t2, text = 'RELAX [0-1]',justify='center')
		self.l2.grid(row=18,column=0,columnspan=2)
		self.e35 = tk.Entry(self.t2,width=5,justify='center')
		self.e35.insert(0,self.StrVar35.get())
		self.e35.grid(row=18, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'NS [#]',justify='center')
		self.l2.grid(row=19,column=0,columnspan=2)
		self.e36 = tk.Entry(self.t2,width=5,justify='center')
		self.e36.insert(0,self.StrVar36.get())
		self.e36.grid(row=19, column=2, columnspan = 2)				

		self.l2 = tk.Label(self.t2, text = 'NT [#]',justify='center')
		self.l2.grid(row=20,column=0,columnspan=2)
		self.e37 = tk.Entry(self.t2,width=5,justify='center')
		self.e37.insert(0,self.StrVar37.get())
		self.e37.grid(row=20, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'NS_MAP [#]',justify='center')
		self.l2.grid(row=21,column=0,columnspan=2)
		self.e38 = tk.Entry(self.t2,width=5,justify='center')
		self.e38.insert(0,self.StrVar38.get())
		self.e38.grid(row=21, column=2, columnspan = 2)				

		self.l2 = tk.Label(self.t2, text = 'NT_MAP [#]',justify='center')
		self.l2.grid(row=22,column=0,columnspan=2)
		self.e39 = tk.Entry(self.t2,width=5,justify='center')
		self.e39.insert(0,self.StrVar39.get())
		self.e39.grid(row=22, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'NSE [#]',justify='center')
		self.l2.grid(row=23,column=0,columnspan=2)
		self.e40 = tk.Entry(self.t2,width=5,justify='center')
		self.e40.insert(0,self.StrVar40.get())
		self.e40.grid(row=23, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'NTE [#]',justify='center')
		self.l2.grid(row=24,column=0,columnspan=2)
		self.e41 = tk.Entry(self.t2,width=5,justify='center')
		self.e41.insert(0,self.StrVar41.get())
		self.e41.grid(row=24, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'NS_MAPE [#]',justify='center')
		self.l2.grid(row=25,column=0,columnspan=2)
		self.e42 = tk.Entry(self.t2,width=5,justify='center')
		self.e42.insert(0,self.StrVar42.get())
		self.e42.grid(row=25, column=2, columnspan = 2)				

		self.l2 = tk.Label(self.t2, text = 'NT_MAPE [#]',justify='center')
		self.l2.grid(row=26,column=0,columnspan=2)
		self.e43 = tk.Entry(self.t2,width=5,justify='center')
		self.e43.insert(0,self.StrVar43.get())
		self.e43.grid(row=26, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = '====== MISHKA opt. ======',justify='center')
		self.l2.grid(row=1,column=4,columnspan=4)

		self.l2 = tk.Label(self.t2, text = 'GRIDN ',justify='center')
		self.l2.grid(row=2,column=4,columnspan=2)
		self.e44 = tk.Entry(self.t2,width=5,justify='center')
		self.e44.insert(0,self.StrVar44.get())
		self.e44.grid(row=2, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'PSISTART ',justify='center')
		self.l2.grid(row=3,column=4,columnspan=2)
		self.e45 = tk.Entry(self.t2,width=5,justify='center')
		self.e45.insert(0,self.StrVar45.get())
		self.e45.grid(row=3, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'XR1 ',justify='center')
		self.l2.grid(row=4,column=4,columnspan=2)
		self.e46 = tk.Entry(self.t2,width=5,justify='center')
		self.e46.insert(0,self.StrVar46.get())
		self.e46.grid(row=4, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'SIG1 ',justify='center')
		self.l2.grid(row=5,column=4,columnspan=2)
		self.e47 = tk.Entry(self.t2,width=5,justify='center')
		self.e47.insert(0,self.StrVar47.get())
		self.e47.grid(row=5, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'XR2 ',justify='center')
		self.l2.grid(row=6,column=4,columnspan=2)
		self.e48 = tk.Entry(self.t2,width=5,justify='center')
		self.e48.insert(0,self.StrVar48.get())
		self.e48.grid(row=6, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'SIG2 ',justify='center')
		self.l2.grid(row=7,column=4,columnspan=2)
		self.e49 = tk.Entry(self.t2,width=5,justify='center')
		self.e49.insert(0,self.StrVar49.get())
		self.e49.grid(row=7, column=6, columnspan = 2)	
		self.l2 = tk.Label(self.t2, text = 'HIGHQ ',justify='center')
		self.l2.grid(row=8,column=4,columnspan=2)		
		self.c9 = tk.Checkbutton(self.t2,variable=self.CheckVar9)
		self.c9.grid(row=8, column=6, columnspan = 2)								

		self.l2 = tk.Label(self.t2, text = '====== ELITE opt. ======',justify='center')
		self.l2.grid(row=10,column=4,columnspan=4)

		self.l2 = tk.Label(self.t2, text = 'COMPRESS ',justify='center')
		self.l2.grid(row=11,column=4,columnspan=2)
		self.c7 = tk.Checkbutton(self.t2,variable=self.CheckVar7)
		self.c7.grid(row=11, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'ROTATION ',justify='center')
		self.l2.grid(row=12,column=4,columnspan=2)
		self.c8 = tk.Checkbutton(self.t2,variable=self.CheckVar8)
		self.c8.grid(row=12, column=6, columnspan = 2)		

		self.l2 = tk.Label(self.t2, text = 'GRIDN ',justify='center')
		self.l2.grid(row=13,column=4,columnspan=2)
		self.e50 = tk.Entry(self.t2,width=5,justify='center')
		self.e50.insert(0,self.StrVar50.get())
		self.e50.grid(row=13, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'NDIST ',justify='center')
		self.l2.grid(row=14,column=4,columnspan=2)
		self.e51 = tk.Entry(self.t2,width=5,justify='center')
		self.e51.insert(0,self.StrVar51.get())
		self.e51.grid(row=14, column=6, columnspan = 2)			

		self.l2 = tk.Label(self.t2, text = 'PSISTART ',justify='center')
		self.l2.grid(row=15,column=4,columnspan=2)
		self.e52 = tk.Entry(self.t2,width=5,justify='center')
		self.e52.insert(0,self.StrVar52.get())
		self.e52.grid(row=15, column=6, columnspan = 2)	


		self.l2 = tk.Label(self.t2, text = '====== Converge opt. ======',justify='center')
		self.l2.grid(row=16,column=4,columnspan=4)

		self.l2 = tk.Label(self.t2, text = 'EPSILON ',justify='center')
		self.l2.grid(row=17,column=4,columnspan=2)
		self.e53 = tk.Entry(self.t2,width=5,justify='center')
		self.e53.insert(0,self.StrVar53.get())
		self.e53.grid(row=17, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'BETA CRIT ',justify='center')
		self.l2.grid(row=18,column=4,columnspan=2)
		self.e54 = tk.Entry(self.t2,width=5,justify='center')
		self.e54.insert(0,self.StrVar54.get())
		self.e54.grid(row=18, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'li CRIT ',justify='center')
		self.l2.grid(row=19,column=4,columnspan=2)
		self.e55 = tk.Entry(self.t2,width=5,justify='center')
		self.e55.insert(0,self.StrVar55.get())
		self.e55.grid(row=19, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'IP CRIT ',justify='center')
		self.l2.grid(row=20,column=4,columnspan=2)
		self.e56 = tk.Entry(self.t2,width=5,justify='center')
		self.e56.insert(0,self.StrVar56.get())
		self.e56.grid(row=20, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'BS CRIT ',justify='center')
		self.l2.grid(row=21,column=4,columnspan=2)
		self.e57 = tk.Entry(self.t2,width=5,justify='center')
		self.e57.insert(0,self.StrVar57.get())
		self.e57.grid(row=21, column=6, columnspan = 2)						

		self.l2 = tk.Label(self.t2, text = 'ITERL [#] ',justify='center')
		self.l2.grid(row=22,column=4,columnspan=2)
		self.e58 = tk.Entry(self.t2,width=5,justify='center')
		self.e58.insert(0,self.StrVar58.get())
		self.e58.grid(row=22, column=6, columnspan = 2)	

		b1 = tk.Button(self.t2, text="SAVE", bg = "lightgray",command=lambda: self.button_func9(),height = 1,width = 4)
		b1.grid(row=27, column=0,columnspan=8,pady=15)


		return

	def button_func8(self):
		
		if not node_machine.lower()=='fusma':	
			self.nodelist = node_default
			return

		load2,used2,name2,nastat2 = read_node_status()

		if not self.t2_close:
			print('>>> NOTE LIST Window is already opened...')
			return

		self.t3 = tk.Toplevel(self.root)
		self.t3.wm_title("NODE LIST")
		self.t3_close = False		
		self.t3.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(3))

		self.l2 = tk.Label(self.t3, text = '====== NODE LIST ======',justify='center')
		self.l2.grid(row=0,column=0,columnspan=8)		

		self.l2 = tk.Label(self.t3, text = '--Node--',justify='center')
		self.l2.grid(row=1,column=0,columnspan=2,padx=5)
		self.l2 = tk.Label(self.t3, text = '--Load--',justify='center')
		self.l2.grid(row=1,column=2,columnspan=2,padx=5)
		self.l2 = tk.Label(self.t3, text = '--Used--',justify='center')
		self.l2.grid(row=1,column=4,columnspan=2,padx=5)				

		for i in range(len(load2)):

			self.l2 = tk.Label(self.t3, text = name2[i] ,justify='center')
			self.l2.grid(row=(i+2),column=0,columnspan=2)

			self.l2 = tk.Label(self.t3, text = load2[i] ,justify='center')
			self.l2.grid(row=(i+2),column=2,columnspan=2)	

			if not nastat2[i] == 'na':
				line = '%i/%i'%(used2[i,1],used2[i,2])
			else:
				line = 'N/A'
			self.l2 = tk.Label(self.t3, text = line ,justify='center')
			self.l2.grid(row=(i+2),column=4,columnspan=2)


		self.noden = len(load2)

		for i in range(20,20+len(load2)):

			self.__dict__['c%d'%i] = tk.Checkbutton(self.t3,variable=self.__dict__['CheckVar%d'%i])
			self.__dict__['c%d'%i].grid(row=i-18,column=6,columnspan=2)

		b1 = tk.Button(self.t3, text="SAVE", bg = "lightgray",command=lambda: self.button_func10(),height = 1,width = 4)
		b1.grid(row=self.noden+4, column=0,columnspan=8,pady=15)			

		return

	def button_func9(self):

		for i in range(25,59):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%(i)].get())					

		self.t2.destroy()
		self.t2_close = True

		return	

	def button_func10(self):

		nodelist = ''
		count = 0
		for i in range(20,20+self.noden):

			if (self.__dict__['CheckVar%d'%i].get() == 1):

				if count == 0:
					nodelist = 'node%02i'%(i-19)
				else:
					nodelist = nodelist + ',node%02i'%(i-19)

				count = count + 1
		if count == 0:
			print(">>> Have to choose more than one...")
			return
		else:
			print(">>> Selected nodes ->",nodelist)
			self.nodelist = nodelist

		self.t3.destroy()
		self.t3_close = True

		return								

	def button_func11(self):

		if not self.t4_close:
			print('>>> NUM Param Window is already opened...')
			return

		self.t4 = tk.Toplevel(self.root)
		self.t4.wm_title("Stability Params")
		self.t4_close = False		
		self.t4.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(4))


		self.l2 = tk.Label(self.t4, text = '=========== Stability opt ===========',justify='center')
		self.l2.grid(row=1,column=0,columnspan=4)

		self.l2 = tk.Label(self.t4, text = 'NQA ',justify='center')
		self.l2.grid(row=2,column=0,columnspan=2)
		self.e59 = tk.Entry(self.t4,width=5,justify='center')
		self.e59.insert(0,self.StrVar59.get())
		self.e59.grid(row=2, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t4, text = 'CRIT_DIA',justify='center')
		self.l2.grid(row=3,column=0,columnspan=2)
		self.e60 =tk.Entry(self.t4,width=5,justify='center')
		self.e60.insert(0,self.StrVar60.get())
		self.e60.grid(row=3, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t4, text = 'CRIT_ALF',justify='center')
		self.l2.grid(row=4,column=0,columnspan=2)
		self.e61 =tk.Entry(self.t4,width=5,justify='center')
		self.e61.insert(0,self.StrVar61.get())
		self.e61.grid(row=4, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t4, text = 'CUT_N',justify='center')
		self.l2.grid(row=5,column=0,columnspan=2)
		self.e62 =tk.Entry(self.t4,width=5,justify='center')
		self.e62.insert(0,self.StrVar62.get())
		self.e62.grid(row=5, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t4, text = 'GRMAX',justify='center')
		self.l2.grid(row=6,column=0,columnspan=2)
		self.e63 =tk.Entry(self.t4,width=5,justify='center')
		self.e63.insert(0,self.StrVar63.get())
		self.e63.grid(row=6, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t4, text = 'USE_BILINEAR',justify='center')
		self.l2.grid(row=7,column=0,columnspan=2)
		self.c12 = tk.Checkbutton(self.t4,variable=self.CheckVar11)
		self.c12.grid(row=7, column=3)


		b1 = tk.Button(self.t4, text="SAVE", bg = "lightgray",command=lambda: self.button_func12(),height = 1,width = 4)
		b1.grid(row=11, column=0,columnspan=8,pady=15)

		return

	def button_func12(self):

		for i in range(59,64):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%(i)].get())
		
		self.t4.destroy()
		self.t4_close = True

	def button_func102(self):

		self.input3 = askopenfilename()
		if (len(self.input3) == 0):
			return
		self.e64.delete(0,'end')
		self.e64.insert(10,self.input3)
		self.CheckVar8.set(1)
		self.isrot = True
		return

	def button_func13(self):

		os.remove('scan_opt_temp')
		if not self.rinput:
			print('>>> Input is not ready...')
			return
		inputname = 'scan_opt'
		self.write_pedscan_input(inputname)

		run_dir = self.currdir + '/' + self.e65.get()
		os.mkdir(run_dir)
		self.make_input_files()
		try:	copytree('PROFILES',run_dir + '/PROFILES')
		except:	pass

		move(inputname,run_dir+'/'+inputname)

		filename = run_dir+'/pedscan_bat'
		log_e = run_dir+'/pedscan_batch.e'
		log_o = run_dir+'/pedscan_batch.o'
		nodelist = node_default

		command = 'cd '+ self.currdir + '\n'
		command = command + pedscane_dir + ' run ' + run_dir+'/'+inputname +' \n'
		make_batch_script(filename,nodelist,log_e,log_o,command,'PEDSCAN')

		os.system('chmod 777 ' + filename)
		if not (self.CheckVar10.get() == 1):
			os.system(filename)
		else:
			runid = submit_batch_script(filename)
		
		self.root.destroy()

		return

	def button_func14(self):

		try:
			result_dir = self.StrVar66.get().split()[0]
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
				f = open('pedscan.log','r')
			except:
				f = open('pedscan.log','w')
				f.close()
				isfile = False
			while isfile:
				line = f.readline()
				if not line: break
				if (line.find(input[-1])>-1):
					ise = True
			if isfile: f.close()
			if not ise:
				os.system('echo ' + input[-1] + ' >> pedscan.log')
			self.run_list = self.open_result_list()
			menu = self.e66["menu"]
			menu.delete(0, "end")
			for string in self.run_list:
				menu.add_command(label=string, command=lambda value=string: self.StrVar66.set(value))

		print('>>> Result dir >',result_dir2)
		if not os.path.isdir(result_dir2):
			print('No available result')
			return

		currdir = os.getcwd()

		import pedscan
		os.chdir(result_dir2)
		ped = pedscan.pedscan('scan_opt',False,True,True)
		ped.work_dir = os.getcwd()
		if not os.path.isfile(ped.work_dir+'/result.dat'):
			print('>>> No available data -> Start collect')
			ped.collect_data()
	
		try:	
			ped.analyse_data('result.dat')
		except:	
			print('>>> No available data -> Start collect')
			ped.collect_data()
			ped.analyse_data('result.dat')
			print('>>> No available result!')
			return

		if not self.t5_close:
			print('>>> NOTE LIST Window is already opened...')
			return

		self.StrVar59.set(str(ped.nqa))
		self.StrVar60.set(str(ped.gr_crit2))
		self.StrVar61.set(str(ped.gr_crit1))
		self.StrVar62.set(str(ped.ncut))
		self.StrVar63.set(str(ped.grmax))
		self.StrVar67.set(str(ped.targetn))
		self.DoubleVar2.set(ped.target_i)
		self.DoubleVar3.set(ped.target_j)

		if ped.use_bilinear:
			self.CheckVar11.set(1)
		else:
			self.CheckVar11.set(0)

		self.use_dia = ped.use_dia
		if ped.use_dia:
			self.CheckVar14.set(1)
		else:
			self.CheckVar14.set(0)

		self.t5 = tk.Toplevel(self.root)
		self.t5.wm_title("OPEN RESULT")
		self.t5_close = False		
		self.t5.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(5))


		self.t5.resizable(0,0)
		self.fig3, ax1 = plt.subplots(1,1,figsize=(7,6))

		#self.fig3.tight_layout()

		self.canvas3 = FigureCanvasTkAgg(self.fig3,master=self.t5)
		self.plot_widget3 = self.canvas3.get_tk_widget()
		self.plot_widget3.grid(rowspan=30,row=1,column=10,columnspan=10)

		toolbar_frame = tk.Frame(self.t5)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar3 = NavigationToolbar2Tk(self.canvas3,toolbar_frame)
		self.fig3.canvas.draw_idle()	

		self.l2 = tk.Label(self.t5, text = '========= POST PROCESSING =========',justify='center')
		self.l2.grid(row=0,column=0,columnspan=8)
		b8 = tk.Button(self.t5, text="Renew data", bg = "light gray",command=lambda: self.button_func18(ped),height = 1,width = 15)
		b8.grid(row=1, column=0,columnspan=8)	
		b9 = tk.Button(self.t5, text="Clear log", bg = "light gray",command=lambda: self.button_func17(ped),height = 1,width = 15)
		b9.grid(row=2, column=0,columnspan=8)
		b9 = tk.Button(self.t5, text="Draw result", bg = "light gray",command=lambda: self.button_func16(ped),height = 1,width = 15)
		b9.grid(row=3, column=0,columnspan=8)			
		b10 = tk.Button(self.t5, text="Make profile", bg = "light gray",command=lambda: self.button_func20(ped),height = 1,width = 15)
		b10.grid(row=4, column=0,columnspan=8)
		b10 = tk.Button(self.t5, text="Plot specturm", bg = "light gray",command=lambda: self.button_func19(ped),height = 1,width = 15)
		b10.grid(row=5, column=0,columnspan=8)							
		b11 = tk.Button(self.t5, text="Change opt", bg = "light gray",command=lambda: self.button_func11(),height = 1,width = 15)
		b11.grid(row=6, column=0,columnspan=8)
		b11 = tk.Button(self.t5, text="EXIT", bg = "light gray",command=lambda: self.button_func23(self.t5,5),height = 1,width = 15)
		b11.grid(row=7, column=0,columnspan=8)		


		self.l7 = tk.Label(self.t5, text="TARGETN ")	
		self.l7.grid(row=8, column=1,columnspan=3)
		self.e67 = tk.Entry(self.t5,width=5,justify='center')
		self.e67.insert(0,self.StrVar67.get())
		self.e67.grid(row=8, column=3, columnspan = 2)

		self.l6 = tk.Label(self.t5, text="Use_dia_effect ",anchor='e')
		self.l6.grid(row=9,column=1,columnspan=3)
		self.c14 = tk.Checkbutton(self.t5,variable=self.CheckVar14)
		self.c14.grid(row=9, column=4)

		self.l6 = tk.Label(self.t5, text=" --INPUT-- ",anchor='e')
		self.l6.grid(row=10,column=0,columnspan=8)

		frame = tk.Frame(self.t5)
		frame.grid(row=11,column=0,columnspan=8,rowspan=14)
		scrollbar = tk.Scrollbar(frame)
		scrollbar.pack(side="right", fill="y")
		scrollbar2 = tk.Scrollbar(frame,orient='horizontal')
		scrollbar2.pack(side="bottom", fill="x")
		listbox = tk.Listbox(frame,yscrollcommand = scrollbar.set,xscrollcommand = scrollbar2.set,width=30)
	
		os.chdir(result_dir2)
		f = open('scan_opt','r')
		count = 1
		while True:
			line = f.readline()
			if not line: break
			line = line.split('\t')[0]
			listbox.insert(count,line.split('\n')[0])
			count = count + 1
		f.close()
		listbox.pack(side="left", fill="x")
		scrollbar["command"]=listbox.yview
		scrollbar2['command'] = listbox.xview

		return

	def button_func15(self):

		self.rinput = True

		if not(self.iseqdsk==1):
			self.rinput = False
			print('>>> EFIT is not selected!')

		if(self.nodelist == 'none'):
			self.rinput = False	
			print('>>> Nodelist is not ready!')

		if (os.path.isdir(self.e65.get())):
			self.rinput = False
			print('>>> Run_dir already exits -> Change the name')

		try:	
			if (len(self.e65.get().split()) == 0):
				self.rinput = False
				print('>>> Run name is not chosed')
		except:
			if (self.e65.get() == ''):
				self.rinput = False
				print('>>> Run name is not chosed')

		if not self.iskin:
			print('>>> Kinetic profile is not ready!')	
		if not (os.path.isfile(self.e2.get())):
			print('>>> Kinetic profile is not valid')
			self.rinput = False				

		if (self.CheckVar2.get() == 1):
			if not self.iskinti:
				print('>>> TI Kinetic profile is not ready')
				self.rinput = False
			if not (os.path.isfile(self.e3.get())):
				print('>>> TI Kinetic profile is not valid')
				self.rinput = False

		if (self.CheckVar8.get() == 1):
			if self.e11.get().lower() == 'mishka':
				print('>>> MISHKA do not support rotation effect')
				self.rinput = False
			if not self.isrot:
				print('>>> Rotation profile is not ready')
				self.rinput = False
			if not (os.path.isfile(self.e64.get())):
				print('>>> Rotation profile is not valid')
				self.rinput = False

		if self.rinput:
			print('>>> Input is ready!')
			self.l100.config(text='READY',fg='lime')
		else:
			self.l100.config(text='NOT READY',fg='red')

		if self.rinput:
			self.write_pedscan_input('scan_opt_temp')
			#self.make_input_files()

		return

	def button_func16(self,sim):

		self.fig3.canvas.draw_idle()
		sim.targetn = int(float(self.e67.get()))
		if self.CheckVar14.get() == 1:
			sim.draw_plot2(self.fig3,True)
		else:
			sim.draw_plot2(self.fig3,False)

		if not self.t6_close:
			self.fig4.canvas.draw_idle()
			if self.CheckVar14.get() == 1:
				sim.draw_plot3(self.fig4,self.DoubleVar4.get(),self.DoubleVar5.get(),True)
			else:
				sim.draw_plot3(self.fig4,self.DoubleVar4.get(),self.DoubleVar5.get(),False)
			
		return

	def button_func17(self,sim):

		sim.clear_data()
		return

	def button_func18(self,sim):

		sim.nqa = float(self.StrVar59.get())
		sim.gr_crit1 = float(self.StrVar61.get())
		sim.gr_crit2 = float(self.StrVar60.get())
		sim.ncut  = int(float(self.StrVar62.get()))
		sim.grmax = float(self.StrVar63.get())
		sim.targetn =int(float(self.e67.get()))
		sim.target_i = self.DoubleVar2.get()
		sim.target_j = self.DoubleVar3.get()


		if self.CheckVar11.get() == 1:
			sim.use_bilinear = True
		else:
			sim.use_bilinear = False		

		if self.CheckVar14.get() == 1:
			sim.use_dia = True
		else:
			sim.use_dia = False		

		self.modify_pedscan_input(sim)

		sim.read_namelist('scan_opt')
		sim.analyse_data('result.dat')

		return

	def button_func19(self,sim):

		if not self.t6_close:
			print('>>> Spectrum Window is already opened...')
			return

		self.t6 = tk.Toplevel(self.root)
		self.t6.wm_title("Growth spectrum")
		self.t6_close = False		
		self.t6.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(6))


		self.t6.resizable(0,0)
		self.fig4, ax1 = plt.subplots(1,1,figsize=(6,6))

		#self.fig3.tight_layout()

		self.canvas4 = FigureCanvasTkAgg(self.fig4,master=self.t6)
		self.plot_widget4 = self.canvas4.get_tk_widget()
		self.plot_widget4.grid(rowspan=30,row=1,column=10,columnspan=10)

		toolbar_frame = tk.Frame(self.t6)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar4 = NavigationToolbar2Tk(self.canvas4,toolbar_frame)
		self.fig4.canvas.draw_idle()			

		self.l7 = tk.Label(self.t6, text='==== Options ====')	
		self.l7.grid(row=0, column=0,columnspan=1)

		self.l7 = tk.Label(self.t6, text="HEIGHT [0-%i] "%(sim.nh-1))	
		self.l7.grid(row=1, column=0)
		self.s2 = tk.Scale(self.t6, variable=self.DoubleVar4, command=lambda x: self.scroll_func2(sim), orient='horizontal', showvalue=True, from_=0,to=sim.nh-1,resolution=1,length=150)
		self.s2.grid(row=2,column=0)

		self.l7 = tk.Label(self.t6, text="WIDTH [0-%i] "%(sim.nw-1))	
		self.l7.grid(row=3, column=0)
		self.s3 = tk.Scale(self.t6, variable=self.DoubleVar5, command=lambda x: self.scroll_func2(sim), orient='horizontal', showvalue=True, from_=0,to=sim.nw-1,resolution=1,length=150)
		self.s3.grid(row=4,column=0)		

		b11 = tk.Button(self.t6, text="EXIT", bg = "light gray",command=lambda: self.button_func23(self.t6,6),height = 1,width = 15)
		b11.grid(row=6, column=0,columnspan=1)		

		if self.CheckVar14.get() == 1:
			sim.draw_plot3(self.fig4,self.DoubleVar4.get(),self.DoubleVar5.get(),True)
		else:
			sim.draw_plot3(self.fig4,self.DoubleVar4.get(),self.DoubleVar5.get(),False)

		return

	def button_func20(self,sim):

		if not self.t8_close:
			print('>>> Profile Window is already opened...')
			return

		self.t8 = tk.Toplevel(self.root)
		self.t8.wm_title("Fitted profile")
		self.t8_close = False		
		self.t8.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(8))

		self.t8.resizable(0,0)
		self.fig6, [ax1, ax2, ax3] = plt.subplots(1,3,figsize=(14,5))

		self.canvas6 = FigureCanvasTkAgg(self.fig6,master=self.t8)
		self.plot_widget6 = self.canvas6.get_tk_widget()
		self.plot_widget6.grid(rowspan=30,row=1,column=10,columnspan=10)

		toolbar_frame = tk.Frame(self.t8)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar6 = NavigationToolbar2Tk(self.canvas6,toolbar_frame)
		self.fig6.canvas.draw_idle()			

		self.l7 = tk.Label(self.t8, text='==== Options ====')	
		self.l7.grid(row=0, column=0,columnspan=4)

		self.l7 = tk.Label(self.t8, text=" HEIGHT [%2.1f - %2.1f] "%(sim.hmin,sim.hmax))	
		self.l7.grid(row=1, column=0, columnspan=2)
		self.e70 = tk.Entry(self.t8,width=5,justify='center')
		self.e70.insert(0,str(self.DoubleVar2.get()))
		self.e70.grid(row=1, column=2, columnspan = 2,padx=10)

		self.l7 = tk.Label(self.t8, text=" WIDTH [%2.1f - %2.1f] "%(sim.wmin,sim.wmax))		
		self.l7.grid(row=2, column=0, columnspan=2)
		self.e71 = tk.Entry(self.t8,width=5,justify='center')
		self.e71.insert(0,str(self.DoubleVar3.get()))
		self.e71.grid(row=2, column=2, columnspan = 2,padx=10)

		b8 = tk.Button(self.t8, text="Draw", bg = "light gray",command=lambda: self.button_func21(sim),height = 1,width = 5)
		b8.grid(row=4, column=1,columnspan=1,sticky='w')
		b8 = tk.Button(self.t8, text="Save", bg = "light gray",command=lambda: self.button_func22(sim),height = 1,width = 5)
		b8.grid(row=4, column=2,columnspan=1,sticky='w',padx=20)

		b11 = tk.Button(self.t8, text="EXIT", bg = "light gray",command=lambda: self.button_func23(self.t8,8),height = 1,width = 15)
		b11.grid(row=5, column=0,columnspan=4)	

		return

	def button_func21(self,sim):

		self.DoubleVar2.set(float(self.e70.get()))
		self.DoubleVar3.set(float(self.e71.get()))

		sim.target_i = self.DoubleVar2.get()
		sim.target_j = self.DoubleVar3.get()
		self.fig6.canvas.draw_idle()
		sim.draw_plot4(self.fig6)		

		self.fig3.canvas.draw_idle()
		if self.CheckVar14.get() == 1:
			sim.draw_plot2(self.fig3,True)
		else:
			sim.draw_plot2(self.fig3,False)

		self.isdraw = True	
		return

	def button_func22(self,sim):

		if not self.isdraw:	
			print('>>> You have to draw it first')
			return

		self.DoubleVar2.set(float(self.e70.get()))
		self.DoubleVar3.set(float(self.e71.get()))

		sim.target_i = self.DoubleVar2.get()
		sim.target_j = self.DoubleVar3.get()

		self.button_func18(sim)
		try:	os.mkdir(sim.work_dir+'/result')
		except:	pass

		nef = interp1d(sim.psin,sim.ne,'cubic')
		tef = interp1d(sim.psin,sim.te,'cubic')
		tif = interp1d(sim.psin,sim.ti,'cubic')

		scale = 1.0 - (sim.zeff-1.)/sim.zimp

		f = open(sim.work_dir+'/result/chease_kinprof_pedscan','w')
		f.write('401\n')
		f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(sim.zeff,sim.zimp,sim.amain,sim.aimp))
		for i in range(401):
			x = float(i)/400.
			f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(x,tef(x),nef(x),tif(x),nef(x)*scale))
		f.close()
		print('>>> Profile is saved to result/chease_kinprof_pedscan')

		self.t8_close = True
		self.t8.destroy()

		return

	def button_func23(self,obj,values):

		if values == 8:	self.t8_close = True
		elif values == 7:	self.t7_close = True
		elif values == 6:	self.t6_close = True
		elif values == 5:	
			self.t5_close = True
			os.chdir(self.currdir)
		elif values == 4:	self.t4_close = True
		elif values == 3:	self.t3_close = True
		elif values == 2:	self.t2_close = True
		elif values == 1:	self.t1_close = True
		obj.destroy()
		return

	def gui_pedscan(self):

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

		#self.initialise_vars()	

		self.l1 = tk.Label(self.root, text="================= INPUT files =================",justify='center')
		self.l1.grid(row=0, column=0,columnspan=9)

		self.l2 = tk.Label(self.root, text="EQU ",anchor='e')
		self.e1 = tk.Entry(self.root,width=25,justify='center')
		self.e1.insert(10,self.StrVar1.get())
		self.l2.grid(row=1, column=0)
		self.e1.grid(row=1, column=1,columnspan=5)
		b1 = tk.Button(self.root, text="Load Eq", bg = "lightgray",command=lambda: self.button_func1(),height = 1,width = 6)
		b1.grid(row=1, column=7)
		b2 = tk.Button(self.root, text="Check Eq", bg = "lightgray",command=lambda: self.button_func2(),height = 1,width = 6)
		b2.grid(row=1, column=8)

		self.l3 = tk.Label(self.root, text="KIN ",anchor='e')
		self.e2 = tk.Entry(self.root,width=25,justify='center')
		self.e2.insert(10,self.StrVar2.get())
		self.l3.grid(row=2, column=0)
		self.e2.grid(row=2, column=1,columnspan=5)
		b3 = tk.Button(self.root, text="Load Kin", bg = "lightgray",command=lambda: self.button_func3(),height = 1,width = 6)
		b3.grid(row=2, column=7)		
		b4 = tk.Button(self.root, text="Use Raw", bg = "lightgray",command=lambda: self.button_func4(),height = 1,width = 6)
		b4.grid(row=2, column=8)

		self.l4 = tk.Label(self.root, text="TIK ",anchor='e')
		self.e3 = tk.Entry(self.root,width=25,justify='center')
		self.e3.insert(10,self.StrVar2.get())
		self.l4.grid(row=3, column=0)
		self.e3.grid(row=3, column=1,columnspan=5)
		b5 = tk.Button(self.root, text="Load Kin", bg = "lightgray",command=lambda: self.button_func5(),height = 1,width = 6)
		b5.grid(row=3, column=7)		
		b6 = tk.Button(self.root, text="Check K", bg = "lightgray",command=lambda: self.button_func6(),height = 1,width = 6)
		b6.grid(row=3, column=8)	

		self.l5 = tk.Label(self.root, text="=============== PEDSCAN option ===============",justify='center')
		self.l5.grid(row=4, column=0,columnspan=9)

		self.l6 = tk.Label(self.root, text="BND_PSI ",anchor='s')	
		self.l6.grid(row=5, column=0,columnspan=3,sticky="se")
		self.s1 = tk.Scale(self.root, variable=self.DoubleVar1, command=lambda x: self.scroll_func1(), orient='horizontal', showvalue=True, from_=0.980,to=0.999,resolution=0.001,length=200)
		self.s1.grid(row=5,column=3,columnspan=7,sticky='s')

		self.l7 = tk.Label(self.root, text="NW",anchor='e')
		self.l7.grid(row=7, column=0,columnspan=2)
		self.e4 = tk.Entry(self.root,width=7,justify='center')
		self.e4.insert(10,self.StrVar4.get())
		self.e4.grid(row=8, column=0,columnspan=2)

		self.l9 = tk.Label(self.root, text="WMIN",anchor='e')
		self.l9.grid(row=7, column=2,columnspan=2)		
		self.e5 = tk.Entry(self.root,width=7,justify='center')
		self.e5.insert(10,self.StrVar5.get())
		self.e5.grid(row=8, column=2,columnspan=2)

		self.l10 = tk.Label(self.root, text="WMAX",anchor='e')
		self.l10.grid(row=7, column=4,columnspan=2)	
		self.e6 = tk.Entry(self.root,width=7,justify='center')
		self.e6.insert(10,self.StrVar6.get())
		self.e6.grid(row=8, column=4,columnspan=2)

		self.l15 = tk.Label(self.root, text="FIXED TI",anchor='e')
		self.l15.grid(row=7, column=6,columnspan=2)							
		self.c1= tk.Checkbutton(self.root,variable=self.CheckVar1, command=lambda: self.checkvar_func1())
		self.c1.grid(row=8, column=6,columnspan=2)

		self.l15 = tk.Label(self.root, text="EXT TI",anchor='e')
		self.l15.grid(row=7, column=8,columnspan=2)							
		self.c2= tk.Checkbutton(self.root,variable=self.CheckVar2, command=lambda: self.checkvar_func2())
		self.c2.grid(row=8, column=8,columnspan=2)		

		self.l12 = tk.Label(self.root, text="NH",anchor='e')
		self.l12.grid(row=9, column=0,columnspan=2)							
		self.e7 = tk.Entry(self.root,width=7,justify='center')
		self.e7.insert(10,self.StrVar7.get())	
		self.e7.grid(row=10, column=0,columnspan=2)	

		self.l13 = tk.Label(self.root, text="HMIN",anchor='e')
		self.l13.grid(row=9, column=2,columnspan=2)							
		self.e8 = tk.Entry(self.root,width=7,justify='center')
		self.e8.insert(10,self.StrVar8.get())	
		self.e8.grid(row=10, column=2,columnspan=2)	

		self.l14 = tk.Label(self.root, text='HMAX',anchor='e')
		self.l14.grid(row=9, column=4,columnspan=2)							
		self.e9 = tk.Entry(self.root,width=7,justify='center')
		self.e9.insert(10,self.StrVar9.get())	
		self.e9.grid(row=10, column=4,columnspan=2)					

		self.l15 = tk.Label(self.root, text="QDEL ",anchor='w')	
		self.l15.grid(row=9, column=6,columnspan=2)
		self.e10 = tk.Entry(self.root,width=7,justify='center')
		self.e10.insert(10,self.StrVar10.get())
		self.e10.grid(row=10, column=6,columnspan=2)		

		self.l16 = tk.Label(self.root, text="================== Stability ==================",justify='center')
		self.l16.grid(row=11, column=0,columnspan=9)
		self.e11 = tk.OptionMenu(self.root,self.StrVar11,'Mishka','Elite')
		self.e11.config(width=5)
		self.e11.grid(row=12, column=0,columnspan=4)

		self.l17 = tk.Label(self.root, text="moden  ",anchor='center')
		self.l17.grid(row=12, column=4,columnspan=2,sticky='e')		
		self.e12 = tk.Entry(self.root,width=20,justify='center')
		self.e12.insert(10,self.StrVar12.get())
		self.e12.grid(row=12, column=4,columnspan=6,sticky='e')


		self.l18 = tk.Label(self.root, text="================ Equil option  ================",justify='center')
		self.l18.grid(row=13, column=0,columnspan=9)
		self.l18 = tk.Label(self.root, text="BS model",anchor='e')
		self.l18.grid(row=14, column=0,columnspan=2)	
		self.e13 = tk.OptionMenu(self.root,self.StrVar13,'Sauter','CSauter','Hager','Neo')
		self.e13.config(width=5)
		self.e13.grid(row=14, column=2,columnspan=3)

		self.l19 = tk.Label(self.root, text="Use_li_crit ",anchor='w')	
		self.l19.grid(row=14, column=5,columnspan=3)
		self.c3= tk.Checkbutton(self.root,variable=self.CheckVar3)
		self.c3.grid(row=14, column=8)		

		self.l20 = tk.Label(self.root, text="Beta const",anchor='e')
		self.l20.grid(row=15, column=0,columnspan=2)	
		self.e14 = tk.OptionMenu(self.root,self.StrVar14,'Bpol','Rmag')
		self.e14.config(width=5)
		self.e14.grid(row=15, column=2,columnspan=3)

		self.l21 = tk.Label(self.root, text="Adjust_prof ",anchor='w')	
		self.l21.grid(row=15, column=5,columnspan=3)
		self.c4= tk.Checkbutton(self.root,variable=self.CheckVar4)
		self.c4.grid(row=15, column=8)		

		self.l22 = tk.Label(self.root, text="BPOL",anchor='e')
		self.l22.grid(row=16, column=0,columnspan=2)
		self.e15 = tk.Entry(self.root,width=5,justify='center')
		self.e15.insert(10,self.StrVar15.get())
		self.e15.grid(row=17, column=0,columnspan=2)		

		self.l23 = tk.Label(self.root, text="RMAG",anchor='e')
		self.l23.grid(row=16, column=2,columnspan=2)
		self.e16 = tk.Entry(self.root,width=5,justify='center')
		self.e16.insert(10,self.StrVar16.get())
		self.e16.grid(row=17, column=2,columnspan=2)

		self.l24 = tk.Label(self.root, text="li",anchor='e')
		self.l24.grid(row=16, column=4,columnspan=2)
		self.e17 = tk.Entry(self.root,width=5,justify='center')
		self.e17.insert(10,self.StrVar17.get())
		self.e17.grid(row=17, column=4,columnspan=2)		

		self.l25 = tk.Label(self.root, text="ZEFF",anchor='e')
		self.l25.grid(row=16, column=6,columnspan=2)
		self.e18 = tk.Entry(self.root,width=5,justify='center')
		self.e18.insert(10,self.StrVar18.get())
		self.e18.grid(row=17, column=6,columnspan=2)

		self.l26 = tk.Label(self.root, text="ZIMP",anchor='e')
		self.l26.grid(row=16, column=8,columnspan=2)
		self.e19 = tk.Entry(self.root,width=5,justify='center')
		self.e19.insert(10,self.StrVar19.get())
		self.e19.grid(row=17, column=8,columnspan=2)

		self.l27 = tk.Label(self.root, text="AMAIN",anchor='e')
		self.l27.grid(row=18, column=0,columnspan=2)
		self.e20 = tk.Entry(self.root,width=5,justify='center')
		self.e20.insert(10,self.StrVar20.get())
		self.e20.grid(row=19, column=0,columnspan=2)		

		self.l28 = tk.Label(self.root, text="AIMP",anchor='e')
		self.l28.grid(row=18, column=2,columnspan=2)
		self.e21 = tk.Entry(self.root,width=5,justify='center')
		self.e21.insert(10,self.StrVar21.get())
		self.e21.grid(row=19, column=2,columnspan=2)

		self.l29 = tk.Label(self.root, text="LINEN",anchor='e')
		self.l29.grid(row=18, column=4,columnspan=2)
		self.e22 = tk.Entry(self.root,width=5,justify='center')
		self.e22.insert(10,self.StrVar22.get())
		self.e22.grid(row=19, column=4,columnspan=2)		

		self.l30 = tk.Label(self.root, text="BSMULTI",anchor='e')
		self.l30.grid(row=18, column=6,columnspan=2)
		self.e23 = tk.Entry(self.root,width=5,justify='center')
		self.e23.insert(10,self.StrVar23.get())
		self.e23.grid(row=19, column=6,columnspan=2)

		self.l30 = tk.Label(self.root, text="CNEO",anchor='e')
		self.l30.grid(row=18, column=8,columnspan=2)
		self.e24 = tk.Entry(self.root,width=5,justify='center')
		self.e24.insert(10,self.StrVar24.get())
		self.e24.grid(row=19, column=8,columnspan=2)

		self.l31 = tk.Label(self.root, text="=============== Other options ===============",justify='center')
		self.l31.grid(row=22, column=0,columnspan=9)
		b5 = tk.Button(self.root, text="NUM opt", bg = "lightgray",command=lambda: self.button_func7(),height = 1,width = 6)
		b5.grid(row=23, column=1,columnspan=3,sticky='e')
		b6 = tk.Button(self.root, text="NODE opt", bg = "lightgray",command=lambda: self.button_func8(),height = 1,width = 6)
		b6.grid(row=23, column=4,columnspan=3)
		b7 = tk.Button(self.root, text="PLOT opt", bg = "lightgray",command=lambda: self.button_func11(),height = 1,width = 6)
		b7.grid(row=23, column=7,columnspan=3,sticky='w')

		self.l32 = tk.Label(self.root, text="ROT_FILE",anchor='e')
		self.l32.grid(row=24, column=0,columnspan=3)
		self.e64 = tk.Entry(self.root,width=20,justify='center')
		self.e64.insert(10,self.StrVar64.get())
		self.e64.grid(row=24, column=3,columnspan=5)	

		b4 = tk.Button(self.root, text="OPEN", bg = "lightgray",command=lambda: self.button_func102(),height = 1,width = 4)
		b4.grid(row=24, column=8)		

		self.l32 = tk.Label(self.root, text="================ Commands ================",justify='center')
		self.l32.grid(row=25, column=0,columnspan=9)

		self.l32 = tk.Label(self.root, text="RUN_NAME",anchor='e')
		self.l32.grid(row=26, column=0,columnspan=3)
		self.e65 = tk.Entry(self.root,width=20,justify='center')
		self.e65.insert(10,self.StrVar65.get())
		self.e65.grid(row=26, column=3,columnspan=5)

		b9 = tk.Button(self.root, text="RUN PEDSCAN", bg = "lightgray",command=lambda: self.button_func13(),height = 1,width = 10)
		b9.grid(row=30, column=0, columnspan=5)

		b10 = tk.Button(self.root, text="CHECK INP", bg = "lightgray",command=lambda: self.button_func15(),height = 1,width = 10)
		b10.grid(row=29, column=0, columnspan=5)

		b11 = tk.Button(self.root, text="OPEN RESULT", bg = "lightgray",command=lambda: self.button_func14(),height = 1,width = 10)
		b11.grid(row=31, column=0, columnspan=5)

		b12 = tk.Button(self.root, text="EXIT", bg = "light gray",command=lambda: self.button_func23(self.root,1),height = 1,width = 10)
		b12.grid(row=32, column=0,columnspan=5)			

		self.l33 = tk.Label(self.root, text="Use batch",anchor='e')
		self.l33.grid(row=30, column=5,columnspan=3)
		self.c9= tk.Checkbutton(self.root,variable=self.CheckVar10)
		self.c9.grid(row=30, column=8)


		self.e66  = tk.OptionMenu(self.root,self.StrVar66,*self.run_list)
		self.e66.config(width=10)
		self.e66.grid(row=31,column=5,columnspan=4)		
		#menu = self.e63.children["menu"]

		#for value in self.run_list:
		self.l100 = tk.Label(self.root, text="NOT READY",fg='red',anchor='center',bg='white',width=15,justify='left')
		self.l100.grid(row=29, column=5,columnspan=4)		

		self.read_preset(2)

		self.root.mainloop()
		return

	def initialise_vars(self):

		self.iseqdsk = False
		self.rinput = False
		self.iskin = False
		self.iskinti = False
		self.isrot = False
		self.isdraw = False

		self.rzbdy = np.zeros(shape=(1,2))
		self.eped_prof = np.zeros(shape=(3,8))
		self.nodelist = node_init

		self.t1_close = True
		self.t2_close = True
		self.t3_close = True
		self.t4_close = True
		self.t5_close = True
		self.t6_close = True
		self.t7_close = True
		self.t8_close = True

		self.fixed_ti = True					#check1
		self.ext_ti = False						#check2		
		self.use_li = False						#check3
		self.adjust_prof = True					#check4
		self.hag_core_mod = True				#check5
		self.use_li2 = False					#check6
		self.use_comp = False					#check7
		self.use_rot = False					#check8
		self.highq = True						#check9
		self.batch_run = True					#check10
		self.use_bilinear = True				#check11
		self.use_dia = True						#check12

		self.target_bnd = 0.995					#double1
		self.target_i = 1						#double2
		self.target_j = 1						#double3
		self.target_i2 = 1						#double4
		self.target_j2 = 1						#double5

		self.equ_name = None					#str1
		self.kin_name = None 					#str2
		self.ti_file = None 					#str3
		self.nw = 8								#str4
		self.wmin = +0.7						#str5
		self.wmax = +1.3						#str6
		self.nh = 8								#str7
		self.hmin = +0.7						#str8
		self.hmax = +1.3						#str9
		self.qdel = 0.3							#str10		
		self.stab_code = 'mishka'				#str11
		self.mode_n = '5,7,10,15,20'			#str12
		self.bsmodel = 'csauter'					#str13
		self.beta_type = 'Bpol'					#str14

		self.bp = 1.0							#str15
		self.rmag = 1.83						#str16
		self.li = 1.0							#str17
		self.zeff = 2.0							#str18
		self.zimp = 6.0							#str19
		self.amain = 2.0						#str20
		self.aimp = 12.0						#str21
		self.linden = 0.0						#str22
		self.bsmulti = 1.0						#str23
		self.core_neo = 0.01					#str24

		self.apf = 1.0							#str 25
		self.bpf = 1.0							#str 26		
		self.cpf = 1.1							#str 27		
		self.dpf = 1.5							#str 28		
		self.ajf = 0.2							#str 29
		self.bjf = 1.0							#str 30		
		self.cjf = 2.0							#str 31		
		self.djf = 2.0							#str 32			
		self.hag_core_mod_psin = 0.3			#str 33

		self.currenti = 35						#str 34
		self.relax = 0.8						#str 35
		self.ns = 80							#str 36	
		self.nt = 80							#str 37
		self.map_ns = 200						#str 38
		self.map_nt = 200						#str 39

		self.nse = 150							#str 40
		self.nte = 200							#str 41
		self.map_nse = 300						#str 42
		self.map_nte = 512						#str 43		

		self.mis_gridn = 301					#str 44
		self.mis_psistart = 0.75				#str 45
		self.xr1 = 1.0							#str 46
		self.sig1 = 0.07						#str 47
		self.xr2 = 0.9							#str 48
		self.sig2 = 0.1							#str 49

		self.eli_gridn = 2000					#str 50
		self.eli_ndist = 50						#str 51
		self.eli_psistart = 0.5					#str 52

		self.epslon = 1.e-8						#str 53
		self.beta_crit = 2.e-2					#str 54
		self.li_crit = 2.e-2					#str 55
		self.ip_crit = 1.e-5					#str 56
		self.bs_crit = 1.e-5					#str 57
		self.niterl = 10						#str 58

		#japlot opt
		
		self.nq_cut_off = 27.7					#str 59
		self.crit_dia = 0.25					#str 60
		self.crit_alf = 0.03					#str 61
		self.cutoff_n = 1						#str 62
		self.grmax = 100						#str 63

		self.eli_rot_file = ''					#str 64	
		self.node_list = node_default	#str 
		self.run_name = ''						#str 65
		self.outdir = '-search-'

		try:	self.read_scanopt('ped_scan/scan_opt')
		except:	print('>>> No prescribed setting 1')

		try:	self.read_scanopt('scan_opt')
		except:	print('>>> No prescribed setting 2') 		

		self.read_preset()
				
		self.DoubleVar1.set(self.target_bnd)
		self.DoubleVar2.set(self.target_i)
		self.DoubleVar3.set(self.target_j)		
		self.DoubleVar4.set(self.target_i2)
		self.DoubleVar5.set(self.target_j2)			

		self.CheckVar1.set(self.trans_vars(self.fixed_ti,4))
		self.CheckVar2.set(self.trans_vars(self.ext_ti,4))
		self.CheckVar3.set(self.trans_vars(self.use_li,4))
		self.CheckVar4.set(self.trans_vars(self.adjust_prof,4))
		self.CheckVar5.set(self.trans_vars(self.hag_core_mod,4))
		self.CheckVar6.set(self.trans_vars(self.use_li2,4))
		self.CheckVar7.set(self.trans_vars(self.use_comp,4))
		self.CheckVar8.set(self.trans_vars(self.use_rot,4))
		self.CheckVar9.set(self.trans_vars(self.highq,4))
		self.CheckVar10.set(self.trans_vars(self.batch_run,4))
		self.CheckVar11.set(self.trans_vars(self.use_bilinear,4))
		self.CheckVar12.set(self.trans_vars(self.use_dia,4))
		
		self.StrVar1.set(self.trans_vars(self.equ_name,3))	
		self.StrVar2.set(self.trans_vars(self.kin_name,3))
		self.StrVar3.set(self.trans_vars(self.ti_file,3))
		self.StrVar4.set(self.trans_vars(self.nw,2))
		self.StrVar5.set(self.trans_vars(self.wmin,2))
		self.StrVar6.set(self.trans_vars(self.wmax,2))	
		self.StrVar7.set(self.trans_vars(self.nh,2))
		self.StrVar8.set(self.trans_vars(self.hmin,2))
		self.StrVar9.set(self.trans_vars(self.hmax,2))
		self.StrVar10.set(self.trans_vars(self.qdel,2))

		self.StrVar11.set(self.trans_vars(self.stab_code,2))
		self.StrVar12.set(self.trans_vars(self.mode_n,1))
		self.StrVar13.set(self.trans_vars(self.bsmodel,1))
		self.StrVar14.set(self.trans_vars(self.beta_type,2))

		self.StrVar15.set(self.trans_vars(self.bp,2))
		self.StrVar16.set(self.trans_vars(self.rmag,2))
		self.StrVar17.set(self.trans_vars(self.li,2))
		self.StrVar18.set(self.trans_vars(self.zeff,2))
		self.StrVar19.set(self.trans_vars(self.zimp,2))
		self.StrVar20.set(self.trans_vars(self.amain,2))
		self.StrVar21.set(self.trans_vars(self.aimp,2))		
		self.StrVar22.set(self.trans_vars(self.linden,2))
		self.StrVar23.set(self.trans_vars(self.bsmulti,2))	
		self.StrVar24.set(self.trans_vars(self.core_neo,2))

		self.StrVar25.set(self.trans_vars(self.apf,2))
		self.StrVar26.set(self.trans_vars(self.bpf,2))	
		self.StrVar27.set(self.trans_vars(self.cpf,2))
		self.StrVar28.set(self.trans_vars(self.dpf,2))
		self.StrVar29.set(self.trans_vars(self.ajf,2))
		self.StrVar30.set(self.trans_vars(self.bjf,2))	
		self.StrVar31.set(self.trans_vars(self.cjf,2))
		self.StrVar32.set(self.trans_vars(self.djf,2))
		self.StrVar33.set(self.trans_vars(self.hag_core_mod_psin,2))
		self.StrVar34.set(self.trans_vars(self.currenti,2))
		self.StrVar35.set(self.trans_vars(self.relax,2))
		self.StrVar36.set(self.trans_vars(self.ns,2))	
		self.StrVar37.set(self.trans_vars(self.nt,2))
		self.StrVar38.set(self.trans_vars(self.map_ns,2))
		self.StrVar39.set(self.trans_vars(self.map_nt,2))
		self.StrVar40.set(self.trans_vars(self.nse,2))	
		self.StrVar41.set(self.trans_vars(self.nte,2))
		self.StrVar42.set(self.trans_vars(self.map_nse,2))
		self.StrVar43.set(self.trans_vars(self.map_nte,2))

		self.StrVar44.set(self.trans_vars(self.mis_gridn,2))
		self.StrVar45.set(self.trans_vars(self.mis_psistart,2))
		self.StrVar46.set(self.trans_vars(self.xr1,2))	
		self.StrVar47.set(self.trans_vars(self.sig1,2))
		self.StrVar48.set(self.trans_vars(self.xr2,2))
		self.StrVar49.set(self.trans_vars(self.sig2,2))

		self.StrVar50.set(self.trans_vars(self.eli_gridn,2))
		self.StrVar51.set(self.trans_vars(self.eli_ndist,2))		
		self.StrVar52.set(self.trans_vars(self.eli_psistart,2))


		self.StrVar53.set(self.trans_vars(self.epslon,2))
		self.StrVar54.set(self.trans_vars(self.beta_crit,2))
		self.StrVar55.set(self.trans_vars(self.li_crit,2))
		self.StrVar56.set(self.trans_vars(self.ip_crit,2))
		self.StrVar57.set(self.trans_vars(self.bs_crit,2))
		self.StrVar58.set(self.trans_vars(self.niterl,2))		

		self.StrVar59.set(self.trans_vars(self.nq_cut_off,2))
		self.StrVar60.set(self.trans_vars(self.crit_dia,2))
		self.StrVar61.set(self.trans_vars(self.crit_alf,2))
		self.StrVar62.set(self.trans_vars(self.cutoff_n,2))
		self.StrVar63.set(self.trans_vars(self.grmax,2))

		self.StrVar64.set(self.trans_vars(self.eli_rot_file,3))		
		self.StrVar65.set(self.trans_vars(self.run_name,3))
		self.StrVar66.set(self.outdir)

		return

	def read_scanopt(self,filename):

		use_neo = False
		use_hager = False
		use_chang = True
		beta_criterion = 1

		f = open(filename,'r')
		while True:
			line = f.readline()
			if not line: break
			self.run_name = read_namelist_str(line,'run_name',self.run_name,3)
			self.use_li = read_namelist_str(line,'use_li',self.use_li,4)
			self.use_li2 = read_namelist_str(line,'use_li2',self.use_li2,4)			
			self.adjust_prof = read_namelist_str(line,'adjust_prof',self.adjust_prof,4)
			self.hag_core_mod = read_namelist_str(line,'hag_core_mod',self.hag_core_mod,4)
			self.use_comp = read_namelist_str(line,'USE_COMP',self.use_comp,4)
			self.use_rot = read_namelist_str(line,'USE_ROT',self.use_rot,4)
			self.fixed_ti = read_namelist_str(line,'FIX_TI',self.fixed_ti,4)
			self.ext_ti = read_namelist_str(line,'EXT_TI',self.ext_ti,4)
			use_neo = read_namelist_str(line,'use_neo',use_neo,4)
			use_hager = read_namelist_str(line,'use_hager',use_hager,4)
			use_chang = read_namelist_str(line,'use_chang',use_chang,4)
			self.target_bnd = read_namelist_str(line,'BND_PSIN',self.target_bnd,2)
			self.equ_name = read_namelist_str(line,'EQDSK',self.equ_name,3)
			self.kin_name = read_namelist_str(line,'EPED_PROF_FILE',self.kin_name,3)
			self.ti_file = read_namelist_str(line,'TI_FILE',self.ti_file,3)
			self.nw  = read_namelist_str(line,'scan_nw',self.nw,1)
			self.nh  = read_namelist_str(line,'scan_nh',self.nh,1)
			self.wmin = read_namelist_str(line,'wmin',self.wmin,2)
			self.wmax = read_namelist_str(line,'wmax',self.wmax,2)
			self.hmin = read_namelist_str(line,'hmin',self.hmin,2)
			self.hmax = read_namelist_str(line,'hmax',self.hmax,2)

			self.stab_code = read_namelist_str(line,'stab_code',self.stab_code,3)
			self.mode_n = read_namelist_str(line,'moden',self.mode_n,3)

			beta_criterion = read_namelist_str(line,'Beta_criterion_type',beta_criterion,1)
			self.bp = read_namelist_str(line,'BPOL',self.bp,2)
			self.rmag = read_namelist_str(line,'RMAG',self.rmag,2)
			self.li = read_namelist_str(line,'LI',self.li,2)
			self.zeff = read_namelist_str(line,'ZEFF',self.zeff,2)
			self.zimp = read_namelist_str(line,'ZIMP',self.zimp,2)
			self.amain = read_namelist_str(line,'AMAIN',self.amain,2)
			self.aimp = read_namelist_str(line,'AIMP',self.aimp,2)
			self.linden = read_namelist_str(line,'LINEDEN',self.linden,2)
			self.bsmulti = read_namelist_str(line,'BSMULTI',self.bsmulti,2)

			self.apf = read_namelist_str(line,'APF',self.apf,2)
			self.bpf = read_namelist_str(line,'BPF',self.bpf,2)
			self.cpf = read_namelist_str(line,'CPF',self.cpf,2)
			self.dpf = read_namelist_str(line,'DPF',self.dpf,2)
			self.ajf = read_namelist_str(line,'AJF',self.ajf,2)
			self.bjf = read_namelist_str(line,'BJF',self.bjf,2)
			self.cjf = read_namelist_str(line,'CJF',self.cjf,2)
			self.djf = read_namelist_str(line,'DJF',self.djf,2)
			self.hag_core_mod_psin = read_namelist_str(line,'hag_core_mod_psin',self.hag_core_mod_psin,2)
			self.currenti = read_namelist_str(line,'Current_ITERN',self.currenti,1)
			self.relax = read_namelist_str(line,'RELAX',self.relax,2)

			self.ns = read_namelist_str(line,'NS',self.ns,1)
			self.nt = read_namelist_str(line,'NT',self.nt,1)
			self.map_ns = read_namelist_str(line,'MAP_NS',self.map_ns,1)
			self.map_nt = read_namelist_str(line,'MAP_NT',self.map_nt,1)
			self.nse = read_namelist_str(line,'NSE',self.nse,1)
			self.nte = read_namelist_str(line,'NTE',self.nte,1)
			self.map_nse = read_namelist_str(line,'MAP_NSE',self.map_nse,1)
			self.map_nte = read_namelist_str(line,'MAP_NTE',self.map_nte,1)

			self.mis_gridn = read_namelist_str(line,'MIS_GRIDN',self.mis_gridn,1)
			self.mis_psistart = read_namelist_str(line,'MIS_PSISTART',self.mis_psistart,2)
			self.xr1 = read_namelist_str(line,'XR1',self.xr1,2)
			self.sig1 = read_namelist_str(line,'SIG1',self.sig1,2)
			self.xr2 = read_namelist_str(line,'XR2',self.xr2,2)
			self.sig2 = read_namelist_str(line,'SIG2',self.sig2,2)

			self.qdel = read_namelist_str(line,'delQ',self.qdel,2)
			self.eli_gridn = read_namelist_str(line,'ELI_GRIDN',self.eli_gridn,1)
			self.eli_ndist = read_namelist_str(line,'ELI_NDIST',self.eli_ndist,1)
			self.eli_psistart = read_namelist_str(line,'ELI_PSISTART',self.eli_psistart,2)
			self.eli_rot_file = read_namelist_str(line,'ROT_FILE',self.eli_rot_file,3)
			self.epslon = read_namelist_str(line,'EPSLON',self.epslon,2)
			self.beta_crit = read_namelist_str(line,'Beta_crit',self.beta_crit,2)
			self.li_crit = read_namelist_str(line,'li_crit',self.li_crit,2)
			self.ip_crit = read_namelist_str(line,'ip_crit',self.ip_crit,2)
			self.bs_crit = read_namelist_str(line,'bs_crit',self.bs_crit,2)
			self.niterl = read_namelist_str(line,'NITERL',self.niterl,1)

			self.nq_cut_off = read_namelist_str(line,'NQA',self.nq_cut_off,2)
			self.crit_alf = read_namelist_str(line,'GRCRIT1',self.crit_alf,2)
			self.crit_dia = read_namelist_str(line,'GRCRIT2',self.crit_dia,2)
			self.cuttoff_n = read_namelist_str(line,'NCUT',self.cutoff_n,1)
			self.grmax = read_namelist_str(line,'GRMAX',self.grmax,2)
			self.use_bilinear = read_namelist_str(line,'USE_BILINEAR',self.use_bilinear,4)
			self.use_dia = read_namelist_str(line,'USE_DIA',self.use_dia,4)
			self.target_i = read_namelist_str(line,'TARGETI',self.target_i,2)
			self.target_j = read_namelist_str(line,'TARGETJ',self.target_j,2)
			self.highq = read_namelist_str(line,'HIGHQ',self.highq,4)

		f.close()
		if use_neo:	self.bsmodel = 'neo'
		elif use_hager:	self.bsmodel = 'hager'
		elif use_chang:	self.bsmodel = 'csauter'
		else:	self.bsmodel = 'sauter'

		if beta_criterion == 1:	self.beta_type = 'bpol'
		else:	self.beta_type = 'rmag'

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
					self.equ_name = line[1]
					self.beta_type = 'rmag'
				elif line[0].lower() == 'kin_file':
					self.kin_name = line[1]
				elif line[0].lower() == 'ti_file2':
					self.ti_file = line[1]
				elif line[0].lower() == 'bnd_psi':
					self.target_bnd = float(line[1])
				elif line[0].lower() == 'run_name_ped':
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
				elif line[0].lower() == 'epedtiprof':
					self.fixed_ti = False

			f.close()
			return
		else:
			while True:
				line = f.readline()
				if not line: break
				line = line.split()
				if line[0].lower() == 'eq_file':
					self.e1.delete(0, 'end')
					self.e1.insert(0, line[1])
					self.button_func1(True,line[1])
				if line[0].lower() == 'kin_file':
					self.e2.delete(0, 'end')
					self.e2.insert(0, line[1])
					self.button_func3(True,line[1])
				if line[0].lower() == 'ti_file2':
					self.e3.delete(0, 'end')
					self.e3.insert(0, line[1])
					self.button_func5(True,line[1])										
				elif line[0].lower() == 'bnd_psi':
					self.DoubleVar1.set(float(line[1]))
				elif line[0].lower() == 'run_name_ped':
					self.e65.delete(0, 'end')
					self.e65.insert(0, line[1])
				elif line[0].lower() == 'epedtiprof':
					self.fixed_ti = False
	

			f.close()
		return					

	def __init__(self):

		self.currdir = os.getcwd()
		return

if __name__ == "__main__":

	import gui_pedscan

	print(' ---------------------------------------------------------------')
	print('||              Python based PedScanner Ver %s                ||'%version['pedscanner'])
	print('||             Linear PBM based pedestal modifier              ||')
	print('%s'%author['pedscanner'])
	print(' ---------------------------------------------------------------\n')

	gps = gui_pedscan.gpedscan()	
	gps.root = tk.Tk()
	gps.root.title('G-PEDSCAN')
	gps.declare_vars()
	gps.initialise_vars()
	gps.run_list = []
	gps.run_list = gps.open_result_list()	
	gps.gui_pedscan()

