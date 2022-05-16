#!/usr/local/anaconda3/bin/python3
import os, sys
import numpy as np
import tkinter as tk
import tkinter.font
from shutil import move, copyfile, copytree
from scipy.interpolate import interp1d
from tkinter.filedialog import askopenfilename,asksaveasfilename, askdirectory
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import matplotlib.pyplot as plt
import eqdsk
import subprocess
from batch_run import *
from exec_dirs import gfit_dir,jastab_dir,helena_exec,chease_exec,mis_exec,elite_dir,python3_exec,node_force,node_init,node_default,version,author,comment

currdir0 = os.getcwd()
class gjastab:

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

	def change_checkvar1(self):

		if self.CheckVar1.get() == 1:
			self.use_kin = True
		else:
			self.use_kin = False

		return

	def change_checkvar2(self):

		if self.CheckVar14.get() == 1:
			self.use_dia = True
		else:
			self.use_dia = False

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
			os.chdir(self.currdir)
			self.t5.destroy()
			self.t5_close = True
#			os.chdir(currdir0)
		elif (win_num == 6):
			self.t6.destroy()
			self.t6_close = True						
		elif (win_num == 7):
			self.t7.destroy()
			self.t7_close = True			

		return		

	def scroll_func1(self):


		self.target_bnd	= self.DoubleVar1.get()

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

	def renew_kin_file(self):
		dirs = os.getcwd() +'/PROFILES/chease_kinprof_fit'

		if (os.path.isfile(dirs)):
			self.e2.delete(0,'end')
			self.e2.insert(10,dirs)
		self.read_kin_file(dirs)
		return

	def read_kin_file(self,dirs):

		if (os.path.isfile(dirs)):

			f = open(dirs,'r')
			line = f.readline()
			line = f.readline().split()

			self.zeff  = float(line[0])
			self.zimp  = float(line[1])
			self.amain = float(line[2])
			self.aimp  = float(line[3])
			self.e19.delete(0,'end')
			self.e19.insert(10,str(round(self.zeff,3)))
			self.e20.delete(0,'end')
			self.e20.insert(10,str(round(self.zimp,3)))
			self.e21.delete(0,'end')
			self.e21.insert(10,str(round(self.amain,3)))
			self.e22.delete(0,'end')
			self.e22.insert(10,str(round(self.aimp,3)))									
			self.iskin = True
		else:
			self.iskin = False
		return

	def make_input_files(self):

		try:
			os.mkdir('input')
		except:
			pass
		currdir = os.getcwd()
		inputdir = currdir + '/input'	

		if self.iseqdsk:
			name = self.e1.get()
			name2 = name.split('/')[-1]
			copyfile(name,inputdir+'/'+name2.split()[0])

		if self.iskin:
			name = self.e2.get()
			name2 = name.split('/')[-1]
			copyfile(name,inputdir+'/'+name2.split()[0])

		if self.isrot:
			name = self.e55.get()
			name2 = name.split('/')[-1]
			copyfile(name,inputdir+'/'+name2.split()[0])

		return	

	def write_jastab_input(self,filename='ja_namelist'):

		f = open(filename,'w')
		f.write('***** JA diagram namelist ***** (by S.Kim ver.1)\n')
		f.write('--- Run option ---\n')
		if self.StrVar3.get().lower() == 'kefit':
			f.write('equ_type = 0\n')
		else:
			f.write('equ_type = 1\n')
		f.write('equ_name = %s \n'%self.e1.get().split('/')[-1])
		f.write('target_bnd = %f\n'%self.target_bnd)
		f.write('run_id = %s \n'%self.e68.get())
		if self.CheckVar9.get() == 1:
			f.write('batch_run = True \n')
		else:
			f.write('batch_run = False\n')
		f.write('	\n')
		f.write('--- Kinetic profile option\n')
		if self.CheckVar1.get() == 1:
			f.write('use_kin_prof = True \n')
		else:
			f.write('use_kin_prof = False \n')
		f.write('kin_prof_name = %s	\n'%self.e2.get().split('/')[-1])
		f.write('\n')
		f.write('--- Scan options ---\n')
		f.write('pedestal_center = %s \n'%self.e7.get())
		if self.CheckVar2.get() == 1:
			f.write('use_ped_finder = True \n')
		else:
			f.write('use_ped_finder = False \n')
		f.write('grid_nw = %s	\n'%self.e4.get())
		f.write('grid_nh = %s	\n'%self.e8.get())
		f.write('pmin = %s		\n'%self.e5.get())
		f.write('pmax = %s		\n'%self.e6.get())
		f.write('jmin = %s		\n'%self.e9.get())
		f.write('jmax = %s		\n'%self.e10.get())
		f.write('\n')
		f.write('--- Equil solver options ---\n')
		f.write('run_equ = %s \n'%self.StrVar11.get().lower())
		if self.CheckVar3.get() == 1:
			f.write('use_prev = True \n')
		else:
			f.write('use_prev = False \n')
		f.write('\n')
		f.write('--- Equilibrium options (equ_type = 1 only)\n')
		f.write('!--constraint\n')
		if self.StrVar15.get().lower() == 'bpol':
			f.write('beta_type = 1\n')
			f.write('beta_target = %s\n'%self.e16.get())
		else:
			f.write('beta_type = 2\n')
			f.write('beta_target = %s\n'%self.e17.get())
		f.write('li_target = %s\n'%self.e18.get())
		if self.CheckVar4.get() == 1:
			f.write('use_li = True    \n')
		else:
			f.write('use_li = False    \n')
		if self.CheckVar15.get() == 1:
			f.write('use_li2 = True    \n')
		else:
			f.write('use_li2 = False    \n')			
		f.write('\n')
		f.write('!-- bootstrap option\n')
		if self.StrVar14.get().lower() == 'nhager':
			f.write('use_neo = True \n')
		else:
			f.write('use_neo = False \n')
		if self.StrVar14.get().lower() == 'hager':	
			f.write('use_hager = True \n')
		else:
			f.write('use_hager = False \n')
		if self.StrVar14.get().lower() == 'csauter':
			f.write('use_chang = True \n')
		else:
			f.write('use_chang = False \n')
		if self.CheckVar6.get() == 1:
			f.write('hag_core_mod=True \n')
		else:
			f.write('hag_core_mod=False \n')
		f.write('hag_core_mod_psin = %s\n'%self.StrVar34.get())
		f.write('core_neo = %s \n'%self.e25.get())
		f.write('bsmulti =  %s \n'%self.e24.get())
		f.write('\n')
		f.write('!-- fast current\n')
		f.write('ajf = %s \n'%self.StrVar30.get())
		f.write('bjf = %s \n'%self.StrVar31.get())
		f.write('cjf = %s \n'%self.StrVar32.get())
		f.write('djf = %s \n'%self.StrVar33.get())
		f.write('\n')
		f.write('!-- fast beam\n')
		f.write('APF = %s \n'%self.StrVar26.get())
		f.write('BPF = %s \n'%self.StrVar27.get())
		f.write('CPF = %s \n'%self.StrVar28.get())
		f.write('DPF = %s \n'%self.StrVar29.get())
		f.write('\n')
		f.write('!-- prof option\n')
		if self.CheckVar5.get() == 1:
			f.write('adjust_prof = True \n')
		else:
			f.write('adjust_prof = False \n')
		f.write('line_den = %s \n'%self.e23.get())
		f.write('\n')
		f.write('--- Stability options ---\n')
		f.write('run_stab =  %s\n'%self.StrVar12.get().lower())
		f.write('highq = False \n')
		f.write('qdelfix =  %s \n'%self.e51.get())
		f.write('qdel_crit = 0.003 \n')
		f.write('\n')
		f.write('--- Stability option (ELITE only) ---\n')
		if self.CheckVar7.get() == 1:
			f.write('use_compression = True		\n')
		else:
			f.write('use_compression = False		\n')
		if self.CheckVar8.get() == 1:
			f.write('use_rotation = True\n')
		else:
			f.write('use_rotation = False\n')
		if not len(self.e55.get().split()) == 0:
			f.write('rot_prof_name = %s\n'%self.e55.get())
		else:
			f.write('rot_prof_name = chease_rot\n')
		f.write('\n')
		f.write('--- CHEASE NUM options ---\n')
		f.write('CNS =    %s\n'%self.StrVar37.get())
		f.write('CNT =    %s\n'%self.StrVar38.get())
		f.write('CNPSI =  %s\n'%self.StrVar39.get())
		f.write('CNCHI =  %s\n'%self.StrVar40.get())
		f.write('\n')
		f.write('CNSE =   %s\n'%self.StrVar41.get())
		f.write('CNTE =   %s\n'%self.StrVar42.get())
		f.write('CNPSIE = %s\n'%self.StrVar43.get())
		f.write('CNCHIE = %s\n'%self.StrVar44.get())
		f.write('\n')
		f.write('EPSILON =   %s\n'%self.StrVar56.get())
		f.write('Beta_crit = %s\n'%self.StrVar57.get())
		f.write('Li_crit =   %s\n'%self.StrVar58.get())
		f.write('Ip_crit =   %s\n'%self.StrVar59.get())
		f.write('Bs_crit =   %s\n'%self.StrVar60.get())
		f.write('Current_ITERN = %s\n'%self.StrVar35.get())
		f.write('Li_ITERN =  %s\n'%self.StrVar61.get())
		f.write('RELAX =     %s\n'%self.StrVar36.get())
		f.write('\n')
		f.write('--- MISHKA NUM options ---\n')
		f.write('moden = %s	\n'%self.e13.get())
		f.write('gridn = %s \n'%self.StrVar45.get())
		f.write('psis =  %s \n'%self.StrVar46.get())
		f.write('xr1 =   %s \n'%self.StrVar47.get())
		f.write('sig1 =  %s \n'%self.StrVar48.get())
		f.write('xr2 =   %s \n'%self.StrVar49.get())
		f.write('sig2 =  %s \n'%self.StrVar50.get())
		f.write('\n')
		f.write('--- ELITE NUM options ---\n')
		f.write('modene = %s \n'%self.e13.get())
		f.write('ngride = %s \n'%self.StrVar52.get())
		f.write('psise =  %s \n'%self.StrVar54.get())
		f.write('ndist =  %s \n'%self.StrVar53.get())
		f.write('\n')
		f.write('--- HELENA NUM options  ---\n')
		f.write('NPT = 301                       \n')
		f.write('b = 0.5                         \n')
		f.write('use_prescribed_helena = False   \n')
		f.write('\n')
		f.write('--- JAPLOT options ---\n')
		f.write('nqa =       %s\n'%self.StrVar63.get())
		f.write('ncutoff =   %s\n'%self.StrVar66.get())
		f.write('grcrit1 =   %s\n'%self.StrVar64.get())
		f.write('grcrit2 =   %s\n'%self.StrVar65.get())
		if self.CheckVar11.get() == 1:
			f.write('use_adja = True\n')
		else:
			f.write('use_adja = False\n')
		if self.CheckVar12.get() == 1:
			f.write('use_bilin = True\n')
		else:
			f.write('use_bilin = False\n')
		f.write('\n')
		f.write('--- Machine options ---\n')
		f.write('nodelist = %s\n'%self.nodelist)
		f.write('\n')
		f.write('--- Code directories ---\n')
		f.write('helena_dir = %s\n'%helena_exec)
		f.write('chease_dir = %s\n'%chease_exec)
		f.write('mis_dir = %s\n'%mis_exec)
		f.write('elite_dir = %s\n'%elite_dir)
		f.write('python_dir = %s\n'%python3_exec)
		f.write('\n')
		f.close()

		return

	def write_jaopt_input(self,filename='ja_opt'):

		f = open(filename,'w')
		if self.CheckVar10.get() == 1:
			f.write('use_elite_a = True \n')
		else:
			f.write('use_elite_a = False \n')
		if self.CheckVar11.get() == 1:	
			f.write('use_adj_a = True \n')
		else:
			f.write('use_adj_a = False \n')
		if self.CheckVar12.get() == 1:	
			f.write('use_bilinear = True \n')
		else:
			f.write('use_bilinear = False \n')
		if self.CheckVar13.get() == 1:	
			f.write('fill_up = True \n')
		else:
			f.write('fill_up = False \n')
		if self.CheckVar16.get() == 1:
			f.write('use_j2 = True \n')
		else:
			f.write('use_j2 = False \n')
		f.write('nq_cut_off = %s\n'%self.StrVar63.get())
		f.write('crit_dia = %s\n'%self.StrVar64.get())
		f.write('crit_alf = %s\n'%self.StrVar65.get())
		f.write('cutoff_n = %s\n'%self.StrVar66.get())
		f.write('grmax = %s\n'%self.StrVar67.get())
		f.close()

		return

	def open_result_list(self):

		run_list =['-search-']
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

	def button_func1(self):

		self.input = askopenfilename()
		if (len(self.input) == 0):
			return
		self.eq = eqdsk.eqdsk(self.input)
		self.eq.read_eqdsk_file()
		self.iseqdsk = True

		self.draw_psirz(self.iseqdsk)
		self.scroll_func1()

		self.update_psirz()

		#filename
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

		self.e16.delete(0,'end')
		self.e16.insert(10,str(self.bp))
		self.e17.delete(0,'end')
		self.e17.insert(10,str(round(self.eq.rmag,3)))
		self.e18.delete(0,'end')
		self.e18.insert(10,str(self.li))

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

	def button_func3(self):

		self.input2 = askopenfilename()
		if (len(self.input2) == 0):
			return
		self.e2.delete(0,'end')
		self.e2.insert(10,self.input2)
		self.CheckVar1.set(1)
		self.read_kin_file(self.input2)
		self.button_func32(self.input2)

		return	

	def button_func32(self,dir):

		if not self.t7_close:
			print('>>> Kinetic profile Window is already opened...')
			return

		try:

			f = open(dir,'r')
			dat_num = int(float(f.readline()))

			result = np.zeros((dat_num,4))
			line = f.readline()
			for i in range(dat_num):
				line = f.readline().split()
				result[i,0] = line[0]
				result[i,1] = line[1]
				result[i,2] = line[2]
				result[i,3] = line[3]

			f.close()
		except:
			print('>>> Kinetic profile is invalid')
			self.iskin = False
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


		ax1.plot(result[:,0],result[:,1])
		ax1.set_title("TE[keV]")
		ax1.set_xlabel('Normalized radius ($\psi_N$)')
		ax1.set_ylabel('TE[keV]')
		ax1.set_xlim((-0.05,1.05))

		ax2.plot(result[:,0],result[:,2])
		ax2.set_title("NE[10(19)/m3]")
		ax2.set_xlabel('Normalized radius ($\psi_N$)')
		ax2.set_ylabel('TE[10(19)/m3]')
		ax2.set_xlim((-0.05,1.05))

		ax3.plot(result[:,0],result[:,3])
		ax3.set_title("TI[keV]")
		ax3.set_xlabel('Normalized radius ($\psi_N$)')
		ax3.set_ylabel('TI[keV]')
		ax3.set_xlim((-0.05,1.05))				

		self.fig5.tight_layout()

		return

	def button_func4(self):

		os.system(gfit_dir)
		self.renew_kin_file()

		return

	def button_func5(self):

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
		self.e26 = tk.Entry(self.t2,width=5,justify='center')
		self.e26.insert(0,self.StrVar26.get())
		self.e26.grid(row=2, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'BPF [0-1]',justify='center')
		self.l2.grid(row=3,column=0,columnspan=2)
		self.e27 =tk.Entry(self.t2,width=5,justify='center')
		self.e27.insert(0,self.StrVar27.get())
		self.e27.grid(row=3, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'CPF ',justify='center')
		self.l2.grid(row=4,column=0,columnspan=2)
		self.e28 = tk.Entry(self.t2,width=5,justify='center')
		self.e28.insert(0,self.StrVar28.get())
		self.e28.grid(row=4, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'DPF ',justify='center')
		self.l2.grid(row=5,column=0,columnspan=2)
		self.e29 = tk.Entry(self.t2,width=5,justify='center')
		self.e29.insert(0,self.StrVar29.get())
		self.e29.grid(row=5, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t2, text = 'AJF ',justify='center')
		self.l2.grid(row=6,column=0,columnspan=2)
		self.e30 = tk.Entry(self.t2,width=5,justify='center')
		self.e30.insert(0,self.StrVar30.get())
		self.e30.grid(row=6, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t2, text = 'BJF [0-1]',justify='center')
		self.l2.grid(row=7,column=0,columnspan=2)
		self.e31 = tk.Entry(self.t2,width=5,justify='center')
		self.e31.insert(0,self.StrVar31.get())
		self.e31.grid(row=7, column=2, columnspan = 2)						

		self.l2 = tk.Label(self.t2, text = 'CJF ',justify='center')
		self.l2.grid(row=8,column=0,columnspan=2)
		self.e32 = tk.Entry(self.t2,width=5,justify='center')
		self.e32.insert(0,self.StrVar32.get())
		self.e32.grid(row=8, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'DJF ',justify='center')
		self.l2.grid(row=9,column=0,columnspan=2)
		self.e33 = tk.Entry(self.t2,width=5,justify='center')
		self.e33.insert(0,self.StrVar33.get())
		self.e33.grid(row=9, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = '====== Current option ======',justify='center')
		self.l2.grid(row=10,column=0,columnspan=4)		

		self.l2 = tk.Label(self.t2, text = 'HAGMOD ',justify='center')
		self.l2.grid(row=11,column=0,columnspan=2)
		self.c6 = tk.Checkbutton(self.t2,variable=self.CheckVar6)
		self.c6.grid(row=11, column=2, columnspan = 2)
		self.l2 = tk.Label(self.t2, text = 'MODPSIN [0-1]',justify='center')
		self.l2.grid(row=12,column=0,columnspan=2)
		self.e34 = tk.Entry(self.t2,width=5,justify='center')
		self.e34.insert(0,self.StrVar34.get())
		self.e34.grid(row=12, column=2, columnspan = 2)
		self.l2 = tk.Label(self.t2, text = 'LI_ITER2 ',justify='center')
		self.l2.grid(row=13,column=0,columnspan=2)
		self.c6 = tk.Checkbutton(self.t2,variable=self.CheckVar15)
		self.c6.grid(row=13, column=2, columnspan = 2)


		self.l2 = tk.Label(self.t2, text = '====== Numerical option ======',justify='center')
		self.l2.grid(row=16,column=0,columnspan=4)		

		self.l2 = tk.Label(self.t2, text = 'ITERN [#]',justify='center')
		self.l2.grid(row=17,column=0,columnspan=2)
		self.e35 = tk.Entry(self.t2,width=5,justify='center')
		self.e35.insert(0,self.StrVar35.get())
		self.e35.grid(row=17, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t2, text = 'RELAX [0-1]',justify='center')
		self.l2.grid(row=18,column=0,columnspan=2)
		self.e36 = tk.Entry(self.t2,width=5,justify='center')
		self.e36.insert(0,self.StrVar36.get())
		self.e36.grid(row=18, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'NS [#]',justify='center')
		self.l2.grid(row=19,column=0,columnspan=2)
		self.e37 = tk.Entry(self.t2,width=5,justify='center')
		self.e37.insert(0,self.StrVar37.get())
		self.e37.grid(row=19, column=2, columnspan = 2)				

		self.l2 = tk.Label(self.t2, text = 'NT [#]',justify='center')
		self.l2.grid(row=20,column=0,columnspan=2)
		self.e38 = tk.Entry(self.t2,width=5,justify='center')
		self.e38.insert(0,self.StrVar38.get())
		self.e38.grid(row=20, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'NS_MAP [#]',justify='center')
		self.l2.grid(row=21,column=0,columnspan=2)
		self.e39 = tk.Entry(self.t2,width=5,justify='center')
		self.e39.insert(0,self.StrVar39.get())
		self.e39.grid(row=21, column=2, columnspan = 2)				

		self.l2 = tk.Label(self.t2, text = 'NT_MAP [#]',justify='center')
		self.l2.grid(row=22,column=0,columnspan=2)
		self.e40 = tk.Entry(self.t2,width=5,justify='center')
		self.e40.insert(0,self.StrVar40.get())
		self.e40.grid(row=22, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'NSE [#]',justify='center')
		self.l2.grid(row=23,column=0,columnspan=2)
		self.e41 = tk.Entry(self.t2,width=5,justify='center')
		self.e41.insert(0,self.StrVar41.get())
		self.e41.grid(row=23, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'NTE [#]',justify='center')
		self.l2.grid(row=24,column=0,columnspan=2)
		self.e42 = tk.Entry(self.t2,width=5,justify='center')
		self.e42.insert(0,self.StrVar42.get())
		self.e42.grid(row=24, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'NS_MAPE [#]',justify='center')
		self.l2.grid(row=25,column=0,columnspan=2)
		self.e43 = tk.Entry(self.t2,width=5,justify='center')
		self.e43.insert(0,self.StrVar43.get())
		self.e43.grid(row=25, column=2, columnspan = 2)				

		self.l2 = tk.Label(self.t2, text = 'NT_MAPE [#]',justify='center')
		self.l2.grid(row=26,column=0,columnspan=2)
		self.e44 = tk.Entry(self.t2,width=5,justify='center')
		self.e44.insert(0,self.StrVar44.get())
		self.e44.grid(row=26, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = '====== MISHKA opt. ======',justify='center')
		self.l2.grid(row=1,column=4,columnspan=4)

		self.l2 = tk.Label(self.t2, text = 'GRIDN ',justify='center')
		self.l2.grid(row=2,column=4,columnspan=2)
		self.e45 = tk.Entry(self.t2,width=5,justify='center')
		self.e45.insert(0,self.StrVar45.get())
		self.e45.grid(row=2, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'PSISTART ',justify='center')
		self.l2.grid(row=3,column=4,columnspan=2)
		self.e46 = tk.Entry(self.t2,width=5,justify='center')
		self.e46.insert(0,self.StrVar46.get())
		self.e46.grid(row=3, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'XR1 ',justify='center')
		self.l2.grid(row=4,column=4,columnspan=2)
		self.e47 = tk.Entry(self.t2,width=5,justify='center')
		self.e47.insert(0,self.StrVar47.get())
		self.e47.grid(row=4, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'SIG1 ',justify='center')
		self.l2.grid(row=5,column=4,columnspan=2)
		self.e48 = tk.Entry(self.t2,width=5,justify='center')
		self.e48.insert(0,self.StrVar48.get())
		self.e48.grid(row=5, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'XR2 ',justify='center')
		self.l2.grid(row=6,column=4,columnspan=2)
		self.e49 = tk.Entry(self.t2,width=5,justify='center')
		self.e49.insert(0,self.StrVar49.get())
		self.e49.grid(row=6, column=6, columnspan = 2)

		self.l2 = tk.Label(self.t2, text = 'SIG2 ',justify='center')
		self.l2.grid(row=7,column=4,columnspan=2)
		self.e50 = tk.Entry(self.t2,width=5,justify='center')
		self.e50.insert(0,self.StrVar50.get())
		self.e50.grid(row=7, column=6, columnspan = 2)						

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
		self.e52 = tk.Entry(self.t2,width=5,justify='center')
		self.e52.insert(0,self.StrVar52.get())
		self.e52.grid(row=13, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'NDIST ',justify='center')
		self.l2.grid(row=14,column=4,columnspan=2)
		self.e53 = tk.Entry(self.t2,width=5,justify='center')
		self.e53.insert(0,self.StrVar53.get())
		self.e53.grid(row=14, column=6, columnspan = 2)			

		self.l2 = tk.Label(self.t2, text = 'PSISTART ',justify='center')
		self.l2.grid(row=15,column=4,columnspan=2)
		self.e54 = tk.Entry(self.t2,width=5,justify='center')
		self.e54.insert(0,self.StrVar54.get())
		self.e54.grid(row=15, column=6, columnspan = 2)	


		self.l2 = tk.Label(self.t2, text = '====== Converge opt. ======',justify='center')
		self.l2.grid(row=16,column=4,columnspan=4)

		self.l2 = tk.Label(self.t2, text = 'EPSILON ',justify='center')
		self.l2.grid(row=17,column=4,columnspan=2)
		self.e56 = tk.Entry(self.t2,width=5,justify='center')
		self.e56.insert(0,self.StrVar56.get())
		self.e56.grid(row=17, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'BETA CRIT ',justify='center')
		self.l2.grid(row=18,column=4,columnspan=2)
		self.e57 = tk.Entry(self.t2,width=5,justify='center')
		self.e57.insert(0,self.StrVar57.get())
		self.e57.grid(row=18, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'li CRIT ',justify='center')
		self.l2.grid(row=19,column=4,columnspan=2)
		self.e58 = tk.Entry(self.t2,width=5,justify='center')
		self.e58.insert(0,self.StrVar58.get())
		self.e58.grid(row=19, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'IP CRIT ',justify='center')
		self.l2.grid(row=20,column=4,columnspan=2)
		self.e59 = tk.Entry(self.t2,width=5,justify='center')
		self.e59.insert(0,self.StrVar59.get())
		self.e59.grid(row=20, column=6, columnspan = 2)	

		self.l2 = tk.Label(self.t2, text = 'BS CRIT ',justify='center')
		self.l2.grid(row=21,column=4,columnspan=2)
		self.e60 = tk.Entry(self.t2,width=5,justify='center')
		self.e60.insert(0,self.StrVar60.get())
		self.e60.grid(row=21, column=6, columnspan = 2)						

		self.l2 = tk.Label(self.t2, text = 'ITERL [#] ',justify='center')
		self.l2.grid(row=22,column=4,columnspan=2)
		self.e61 = tk.Entry(self.t2,width=5,justify='center')
		self.e61.insert(0,self.StrVar61.get())
		self.e61.grid(row=22, column=6, columnspan = 2)	

		b1 = tk.Button(self.t2, text="SAVE", bg = "lightgray",command=lambda: self.button_func7(),height = 1,width = 4)
		b1.grid(row=27, column=0,columnspan=8,pady=15)


		return

	def button_func6(self):

		if node_force:
			self.nodelist = 'compute'
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

		b1 = tk.Button(self.t3, text="SAVE", bg = "lightgray",command=lambda: self.button_func8(),height = 1,width = 4)
		b1.grid(row=self.noden+4, column=0,columnspan=8,pady=15)			

		return

	def button_func7(self):

		for i in range(26,51):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%(i)].get())

		for i in range(52,55):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%(i)].get())	

		for i in range(56,62):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%(i)].get())						

		self.t2.destroy()
		self.t2_close = True

		return	

	def button_func8(self):

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

	def button_func9(self):

		if not self.t4_close:
			print('>>> NUM Param Window is already opened...')
			return

		self.t4 = tk.Toplevel(self.root)
		self.t4.wm_title("JAPLOT Params")
		self.t4_close = False		
		self.t4.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(4))


		self.l2 = tk.Label(self.t4, text = '=========== JAPLOT ===========',justify='center')
		self.l2.grid(row=1,column=0,columnspan=4)

		self.l2 = tk.Label(self.t4, text = 'NQA ',justify='center')
		self.l2.grid(row=2,column=0,columnspan=2)
		self.e63 = tk.Entry(self.t4,width=5,justify='center')
		self.e63.insert(0,self.StrVar63.get())
		self.e63.grid(row=2, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t4, text = 'CRIT_DIA',justify='center')
		self.l2.grid(row=3,column=0,columnspan=2)
		self.e64 =tk.Entry(self.t4,width=5,justify='center')
		self.e64.insert(0,self.StrVar64.get())
		self.e64.grid(row=3, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t4, text = 'CRIT_ALF',justify='center')
		self.l2.grid(row=4,column=0,columnspan=2)
		self.e65 =tk.Entry(self.t4,width=5,justify='center')
		self.e65.insert(0,self.StrVar65.get())
		self.e65.grid(row=4, column=2, columnspan = 2)		

		self.l2 = tk.Label(self.t4, text = 'CUT_N',justify='center')
		self.l2.grid(row=5,column=0,columnspan=2)
		self.e66 =tk.Entry(self.t4,width=5,justify='center')
		self.e66.insert(0,self.StrVar66.get())
		self.e66.grid(row=5, column=2, columnspan = 2)

		self.l2 = tk.Label(self.t4, text = 'GRMAX',justify='center')
		self.l2.grid(row=6,column=0,columnspan=2)
		self.e67 =tk.Entry(self.t4,width=5,justify='center')
		self.e67.insert(0,self.StrVar67.get())
		self.e67.grid(row=6, column=2, columnspan = 2)	

		self.l2 = tk.Label(self.t4, text = 'USE_ELITE_A',justify='center')
		self.l2.grid(row=7,column=0,columnspan=2)
		self.c10 = tk.Checkbutton(self.t4,variable=self.CheckVar10)
		self.c10.grid(row=7, column=3)

		self.l2 = tk.Label(self.t4, text = 'USE_ADJ_A',justify='center')
		self.l2.grid(row=8,column=0,columnspan=2)
		self.c11 = tk.Checkbutton(self.t4,variable=self.CheckVar11)
		self.c11.grid(row=8, column=3)

		self.l2 = tk.Label(self.t4, text = 'USE_BILINEAR',justify='center')
		self.l2.grid(row=9,column=0,columnspan=2)
		self.c12 = tk.Checkbutton(self.t4,variable=self.CheckVar12)
		self.c12.grid(row=9, column=3)

		self.l2 = tk.Label(self.t4, text = 'USE_FILLUP',justify='center')
		self.l2.grid(row=10,column=0,columnspan=2)
		self.c13 = tk.Checkbutton(self.t4,variable=self.CheckVar13)
		self.c13.grid(row=10, column=3)

		self.l2 = tk.Label(self.t4, text = 'USE_JEDGE',justify='center')
		self.l2.grid(row=11,column=0,columnspan=2)
		self.c14 = tk.Checkbutton(self.t4,variable=self.CheckVar16)
		self.c14.grid(row=11, column=3)

		b1 = tk.Button(self.t4, text="SAVE", bg = "lightgray",command=lambda: self.button_func10(),height = 1,width = 4)
		b1.grid(row=12, column=0,columnspan=8,pady=15)

		return

	def button_func10(self):

		for i in range(63,68):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%(i)].get())
		
		self.t4.destroy()
		self.t4_close = True

	def button_func102(self):

		self.input3 = askopenfilename()
		if (len(self.input3) == 0):
			return
		self.e55.delete(0,'end')
		self.e55.insert(10,self.input3)
		self.CheckVar8.set(1)
		self.isrot = True
		return

	def button_func11(self):

		if not self.rinput:
			print('>>> Input is not ready...')
			return
		inputname = 'ja_namelist_'+self.e68.get()
		os.remove('ja_namelist_temp')
		self.write_jastab_input(inputname)
		self.root.destroy()
		os.system(jastab_dir + ' ' + inputname)
		return

	def button_func12(self):

		try:
			result_dir = self.StrVar69.get().split()[0]
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
			menu = self.e69["menu"]
			menu.delete(0, "end")
			for string in self.run_list:
				menu.add_command(label=string, command=lambda value=string: self.StrVar69.set(value))

		print('>>> Result dir >',result_dir2)
		if not os.path.isdir(result_dir2):
			print('No available result')
			return

		currdir = os.getcwd()
		os.chdir(result_dir2)

		stat,out = subprocess.getstatusoutput('ls ja_namelist*')
		out = out.split('\n')

		import ja_tool
		import gjaplot

		input_name = out[0]
		sim = ja_tool.simulation(input_name)
		sim.run_dir = os.getcwd()
		sim.batch_run = False

		if (sim.run_stab == 'mishka'):
			resname = 'ja_result_' + sim.RUN_ID+'_mis'
		elif (sim.run_stab == 'elite'):
			resname = 'ja_result_' + sim.RUN_ID+'_eli'

		if os.path.isfile(resname):
			f = open(resname,'r')
			for i in range(2): line = f.readline();
			f.close()
			if not len(line.split())>3: 
				print('>>> OLD result file, delete...')
				os.remove(resname)	
		
		if not (os.path.isfile(resname)) :
			print('>>> No available result -> Do collect')
			try:
				if (sim.run_equ == 'helena'):
					sim.generate_output_helena()
				elif (sim.run_equ == 'chease'):
					sim.generate_output_chease()
			except:
				pass

			if not (os.path.isfile(resname)) :	
				print('>>> No available result -> Problem in runs')
				return			

		if not self.t5_close:
			print('>>> NOTE LIST Window is already opened...')
			return

		self.gja = gjaplot.japlot(resname)
		self.StrVar63.set(str(self.gja.nq_cut_off))
		self.StrVar64.set(str(self.gja.crit_dia))
		self.StrVar65.set(str(self.gja.crit_alf))
		self.StrVar66.set(str(self.gja.cutoff_n))
		self.StrVar67.set(str(self.gja.grmax))
		self.CheckVar10.set(self.trans_vars(self.gja.use_elite_a,4))
		self.CheckVar11.set(self.trans_vars(self.gja.use_adj_a,4))
		self.CheckVar12.set(self.trans_vars(self.gja.use_bilinear,4))
		self.CheckVar13.set(self.trans_vars(self.gja.fill_up,4))
		self.CheckVar16.set(self.trans_vars(self.gja.use_j2,4))
		self.DoubleVar2.set(self.gja.target_i)
		self.DoubleVar3.set(self.gja.target_j)
		self.gja.post_processing()

		self.t5 = tk.Toplevel(self.root)
		self.t5.wm_title("OPEN RESULT")
		self.t5_close = False		
		self.t5.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(5))


		self.t5.resizable(0,0)
		self.fig3, ax1 = plt.subplots(1,1,figsize=(6,6))

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
		b8 = tk.Button(self.t5, text="Renew data", bg = "light gray",command=lambda: self.button_func16(),height = 1,width = 15)
		b8.grid(row=1, column=0,columnspan=8)	
		b9 = tk.Button(self.t5, text="Clear log", bg = "light gray",command=lambda: self.button_func15(sim),height = 1,width = 15)
		b9.grid(row=2, column=0,columnspan=8)
		b9 = tk.Button(self.t5, text="Clear all log", bg = "light gray",command=lambda: self.button_func15(sim,True),height = 1,width = 15)
		b9.grid(row=3, column=0,columnspan=8)			
		b10 = tk.Button(self.t5, text="Draw result", bg = "light gray",command=lambda: self.button_func14(),height = 1,width = 15)
		b10.grid(row=4, column=0,columnspan=8)
		b10 = tk.Button(self.t5, text="Plot specturm", bg = "light gray",command=lambda: self.button_func17(),height = 1,width = 15)
		b10.grid(row=5, column=0,columnspan=8)							
		b11 = tk.Button(self.t5, text="Change opt", bg = "light gray",command=lambda: self.button_func9(),height = 1,width = 15)
		b11.grid(row=6, column=0,columnspan=8)
		b12 = tk.Button(self.t5, text="EXIT", bg = "light gray",command=lambda: self.button_func18(self.t5,5),height = 1,width = 15)
		b12.grid(row=7, column=0,columnspan=8)		


		self.l7 = tk.Label(self.t5, text="NWI ",anchor='s')	
		self.l7.grid(row=8, column=0,columnspan=3,sticky="se")
		self.s2 = tk.Scale(self.t5, variable=self.DoubleVar2, command=lambda x: self.scroll_func2(), orient='horizontal', showvalue=True, from_=0,to=self.gja.gridn1-1,resolution=1,length=200)
		self.s2.grid(row=8,column=3,columnspan=7,sticky='s')

		self.l7 = tk.Label(self.t5, text="NHI ",anchor='s')	
		self.l7.grid(row=9, column=0,columnspan=3,sticky="se")
		self.s3 = tk.Scale(self.t5, variable=self.DoubleVar3, command=lambda x: self.scroll_func2(), orient='horizontal', showvalue=True, from_=0,to=self.gja.gridn2-1,resolution=1,length=200)
		self.s3.grid(row=9,column=3,columnspan=7,sticky='s')		

		self.l6 = tk.Label(self.t5, text="Use_dia_effect ",anchor='e')
		self.l6.grid(row=10,column=3,columnspan=3)
		self.c14 = tk.Checkbutton(self.t5,variable=self.CheckVar14, command=lambda: self.change_checkvar2())
		self.c14.grid(row=10, column=6)

		self.l6 = tk.Label(self.t5, text="--INPUT-- ",anchor='e')
		self.l6.grid(row=11,column=3,columnspan=4)

		frame = tk.Frame(self.t5)
		frame.grid(row=12,column=0,columnspan=8,rowspan=14)
		scrollbar = tk.Scrollbar(frame)
		scrollbar.pack(side="right", fill="y")
		scrollbar2 = tk.Scrollbar(frame,orient='horizontal')
		scrollbar2.pack(side="bottom", fill="x")
		listbox = tk.Listbox(frame,yscrollcommand = scrollbar.set,xscrollcommand = scrollbar2.set,width=30)
	
		os.chdir(result_dir2)
		f = open(input_name,'r')
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

	def button_func13(self):

		self.rinput = True

		if not(self.iseqdsk==1):
			self.rinput = False
			print('>>> EFIT is not selected!')

		if(self.nodelist == 'none'):
			self.rinput = False	
			print('>>> Nodelist is not ready!')

		if (os.path.isdir(self.e68.get())):
			self.rinput = False
			print('>>> Run_dir already exits -> Change the name')

		if (self.StrVar3.get().lower()=='mefit'):
			if not self.iskin:
				print('>>> Kinetic profile must be given for MEFIT type')

		try:	
			if (len(self.e68.get().split()) == 0):
				self.rinput = False
				print('>>> Run name is not chosed')
		except:
			if (self.e68.get() == ''):
				self.rinput = False
				print('>>> Run name is not chosed')

		if (self.CheckVar1.get() == 1):
			if not self.iskin:
				print('>>> Kinetic profile is not ready')
				self.rinput = False
			if not (os.path.isfile(self.e2.get())):
				print('>>> Kinetic profile is not valid')
				self.rinput = False

		if (self.CheckVar8.get() == 1):
			if self.StrVar12.get().lower() == 'mishka':
				print('>>> MISHKA do not support rotation effect')
				self.rinput = False
			if not self.isrot:
				print('>>> Rotation profile is not ready')
				self.rinput = False
			if not (os.path.isfile(self.e55.get())):
				print('>>> Rotation profile is not valid')
				self.rinput = False

		if self.rinput:
			print('>>> Input is ready!')
			self.l100.config(text='READY',fg='lime')
		else:
			self.l100.config(text='NOT READY',fg='red')

		if self.rinput:
			self.write_jastab_input('ja_namelist_temp')
			self.make_input_files()

		return

	def button_func14(self,use_dia=True):

		self.fig3.canvas.draw_idle()
		if self.use_dia:
			self.gja.draw_plot(2,self.fig3)
		else:
			self.gja.draw_plot(1,self.fig3)

		if not self.t6_close:
			self.fig4.canvas.draw_idle()
			if self.use_dia:
				self.gja.draw_plot(6,self.fig4)
			else:
				self.gja.draw_plot(5,self.fig4)			
			
		return

	def button_func15(self,sim,call=False):

		sim.clear_equ_dat(call)
		return

	def button_func16(self):

		self.gja.nq_cut_off = float(self.StrVar63.get())
		self.gja.crit_dia = float(self.StrVar64.get())
		self.gja.crit_alf = float(self.StrVar65.get())
		self.gja.cutoff_n = int(float(self.StrVar66.get()))
		self.gja.grmax = float(self.StrVar67.get())

		if self.CheckVar10.get() == 1:
			self.gja.use_elite_a = True
		else:
			self.gja.use_elite_a = False

		if self.CheckVar11.get() == 1:
			self.gja.use_adj_a = True
		else:
			self.gja.use_adj_a = False			

		if self.CheckVar12.get() == 1:
			self.gja.use_bilinear = True
		else:
			self.gja.use_bilinear = False		

		if self.CheckVar13.get() == 1:
			self.gja.fill_up = True
		else:
			self.gja.fill_up = False

		if self.CheckVar16.get() == 1:
			self.gja.use_j2 = True
		else:
			self.gja.use_j2 = False	

		self.write_jaopt_input()
		self.gja.post_processing()

		return

	def button_func17(self,use_dia=True):

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

		self.gja.target_i = int(self.DoubleVar2.get())
		self.gja.target_j = int(self.DoubleVar3.get())
		self.fig4.canvas.draw_idle()
		if self.use_dia:
			self.gja.draw_plot(6,self.fig4)
		else:
			self.gja.draw_plot(5,self.fig4)
		return

	def button_func18(self,obj,values):

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

	def scroll_func2(self,use_dia=True):

		self.gja.target_i = int(self.DoubleVar2.get())
		self.gja.target_j = int(self.DoubleVar3.get())	

		if not self.t6_close:
			self.fig4.canvas.draw_idle()
			if self.use_dia:
				self.gja.draw_plot(6,self.fig4)
			else:
				self.gja.draw_plot(5,self.fig4)
		
		self.fig3.canvas.draw_idle()

		if self.use_dia:
			self.gja.draw_plot(2,self.fig3)
		else:
			self.gja.draw_plot(1,self.fig3)
		return		


		return

	def gui_jastab(self):

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

		self.initialise_vars()	

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

		self.l4 = tk.Label(self.root, text="=============== JASTAB option ===============",justify='center')
		self.l4.grid(row=3, column=0,columnspan=9)
		self.l5 = tk.Label(self.root, text="EFIT Type ",anchor='e')	
		self.l5.grid(row=4, column=0,columnspan=2)
		self.e3 = tk.OptionMenu(self.root,self.StrVar3,'KEFIT','MEFIT')
		self.e3.config(width=5)
		self.e3.grid(row=4,column=2,columnspan=3)
		self.l6 = tk.Label(self.root, text="Use_Kin ",anchor='e')
		self.l6.grid(row=4,column=6,columnspan=2)
		self.c1 = tk.Checkbutton(self.root,variable=self.CheckVar1, command=lambda: self.change_checkvar1())
		self.c1.grid(row=4, column=8)

		self.l7 = tk.Label(self.root, text="BND_PSI ",anchor='s')	
		self.l7.grid(row=5, column=0,columnspan=3,sticky="se")
		self.s1 = tk.Scale(self.root, variable=self.DoubleVar1, command=lambda x: self.scroll_func1(), orient='horizontal', showvalue=True, from_=0.980,to=0.999,resolution=0.001,length=200)
		self.s1.grid(row=5,column=3,columnspan=7,sticky='s')

		self.l8 = tk.Label(self.root, text="NW",anchor='e')
		self.l8.grid(row=7, column=0,columnspan=2)
		self.e4 = tk.Entry(self.root,width=7,justify='center')
		self.e4.insert(10,self.StrVar4.get())
		self.e4.grid(row=8, column=0,columnspan=2)

		self.l9 = tk.Label(self.root, text="PMIN",anchor='e')
		self.l9.grid(row=7, column=2,columnspan=2)		
		self.e5 = tk.Entry(self.root,width=7,justify='center')
		self.e5.insert(10,self.StrVar5.get())
		self.e5.grid(row=8, column=2,columnspan=2)

		self.l10 = tk.Label(self.root, text="PMAX",anchor='e')
		self.l10.grid(row=7, column=4,columnspan=2)	
		self.e6 = tk.Entry(self.root,width=7,justify='center')
		self.e6.insert(10,self.StrVar6.get())
		self.e6.grid(row=8, column=4,columnspan=2)

		self.l11 = tk.Label(self.root, text="PEDCENTER",anchor='e')
		self.l11.grid(row=7, column=6,columnspan=4)	
		self.e7 = tk.Entry(self.root,width=10,justify='center')
		self.e7.insert(10,self.StrVar7.get())
		self.e7.grid(row=8, column=6,columnspan=4)

		self.l12 = tk.Label(self.root, text="NH",anchor='e')
		self.l12.grid(row=9, column=0,columnspan=2)							
		self.e8 = tk.Entry(self.root,width=7,justify='center')
		self.e8.insert(10,self.StrVar8.get())	
		self.e8.grid(row=10, column=0,columnspan=2)	

		self.l13 = tk.Label(self.root, text="JMIN",anchor='e')
		self.l13.grid(row=9, column=2,columnspan=2)							
		self.e9 = tk.Entry(self.root,width=7,justify='center')
		self.e9.insert(10,self.StrVar9.get())	
		self.e9.grid(row=10, column=2,columnspan=2)	

		self.l14 = tk.Label(self.root, text='JMAX',anchor='e')
		self.l14.grid(row=9, column=4,columnspan=2)							
		self.e10 = tk.Entry(self.root,width=7,justify='center')
		self.e10.insert(10,self.StrVar10.get())	
		self.e10.grid(row=10, column=4,columnspan=2)					

		self.l15 = tk.Label(self.root, text="AutoPedFinder",anchor='e')
		self.l15.grid(row=9, column=6,columnspan=4)							
		self.c2= tk.Checkbutton(self.root,variable=self.CheckVar2)
		self.c2.grid(row=10, column=6,columnspan=4)

		self.l15 = tk.Label(self.root, text="================ EQU Solver ================",justify='center')
		self.l15.grid(row=11, column=0,columnspan=9)
		self.e11 = tk.OptionMenu(self.root,self.StrVar11,'Chease','Helena')
		self.e11.config(width=5)
		self.e11.grid(row=12, column=0,columnspan=4)			
		self.l17 = tk.Label(self.root, text="      Use_prev ",anchor='e')	
		self.l17.grid(row=12, column=7,columnspan=2,sticky='w')
		self.c3= tk.Checkbutton(self.root,variable=self.CheckVar3)
		self.c3.grid(row=12, column=8)

		self.l17 = tk.Label(self.root, text="QDEL ",anchor='w')	
		self.l17.grid(row=12, column=4,columnspan=2,sticky='w')
		self.e51 = tk.Entry(self.root,width=4,justify='center')
		self.e51.insert(10,self.StrVar51.get())
		self.e51.grid(row=12, column=5,columnspan=2,sticky='w')		

		self.l18 = tk.Label(self.root, text="================== Stability ==================",justify='center')
		self.l18.grid(row=13, column=0,columnspan=9)
		self.e12 = tk.OptionMenu(self.root,self.StrVar12,'Mishka','Elite')
		self.e12.config(width=5)
		self.e12.grid(row=14, column=0,columnspan=4)

		self.l19 = tk.Label(self.root, text="moden  ",anchor='center')
		self.l19.grid(row=14, column=4,columnspan=2,sticky='e')		
		self.e13 = tk.Entry(self.root,width=20,justify='center')
		self.e13.insert(10,self.StrVar13.get())
		self.e13.grid(row=14, column=4,columnspan=6,sticky='e')


		self.l20 = tk.Label(self.root, text="=========== Equil option (MEFIT only) ===========",justify='center')
		self.l20.grid(row=15, column=0,columnspan=9)
		self.l21 = tk.Label(self.root, text="BS model",anchor='e')
		self.l21.grid(row=16, column=0,columnspan=2)	
		self.e14 = tk.OptionMenu(self.root,self.StrVar14,'Sauter','CSauter','Hager','NHager')
		self.e14.config(width=5)
		self.e14.grid(row=16, column=2,columnspan=3)

		self.l22 = tk.Label(self.root, text="Use_li_crit ",anchor='w')	
		self.l22.grid(row=16, column=5,columnspan=3)
		self.c4= tk.Checkbutton(self.root,variable=self.CheckVar4)
		self.c4.grid(row=16, column=8)		

		self.l23 = tk.Label(self.root, text="Beta const",anchor='e')
		self.l23.grid(row=17, column=0,columnspan=2)	
		self.e15 = tk.OptionMenu(self.root,self.StrVar15,'Bpol','Rmag')
		self.e15.config(width=5)
		self.e15.grid(row=17, column=2,columnspan=3)

		self.l23 = tk.Label(self.root, text="Adjust_prof ",anchor='w')	
		self.l23.grid(row=17, column=5,columnspan=3)
		self.c5= tk.Checkbutton(self.root,variable=self.CheckVar5)
		self.c5.grid(row=17, column=8)		

		self.l24 = tk.Label(self.root, text="BPOL",anchor='e')
		self.l24.grid(row=18, column=0,columnspan=2)
		self.e16 = tk.Entry(self.root,width=5,justify='center')
		self.e16.insert(10,self.StrVar16.get())
		self.e16.grid(row=19, column=0,columnspan=2)		

		self.l24 = tk.Label(self.root, text="RMAG",anchor='e')
		self.l24.grid(row=18, column=2,columnspan=2)
		self.e17 = tk.Entry(self.root,width=5,justify='center')
		self.e17.insert(10,self.StrVar17.get())
		self.e17.grid(row=19, column=2,columnspan=2)

		self.l24 = tk.Label(self.root, text="li",anchor='e')
		self.l24.grid(row=18, column=4,columnspan=2)
		self.e18 = tk.Entry(self.root,width=5,justify='center')
		self.e18.insert(10,self.StrVar18.get())
		self.e18.grid(row=19, column=4,columnspan=2)		

		self.l25 = tk.Label(self.root, text="ZEFF",anchor='e')
		self.l25.grid(row=18, column=6,columnspan=2)
		self.e19 = tk.Entry(self.root,width=5,justify='center')
		self.e19.insert(10,self.StrVar19.get())
		self.e19.grid(row=19, column=6,columnspan=2)

		self.l26 = tk.Label(self.root, text="ZIMP",anchor='e')
		self.l26.grid(row=18, column=8,columnspan=2)
		self.e20 = tk.Entry(self.root,width=5,justify='center')
		self.e20.insert(10,self.StrVar20.get())
		self.e20.grid(row=19, column=8,columnspan=2)

		self.l27 = tk.Label(self.root, text="AMAIN",anchor='e')
		self.l27.grid(row=20, column=0,columnspan=2)
		self.e21 = tk.Entry(self.root,width=5,justify='center')
		self.e21.insert(10,self.StrVar21.get())
		self.e21.grid(row=21, column=0,columnspan=2)		

		self.l28 = tk.Label(self.root, text="AIMP",anchor='e')
		self.l28.grid(row=20, column=2,columnspan=2)
		self.e22 = tk.Entry(self.root,width=5,justify='center')
		self.e22.insert(10,self.StrVar22.get())
		self.e22.grid(row=21, column=2,columnspan=2)

		self.l29 = tk.Label(self.root, text="LINEN",anchor='e')
		self.l29.grid(row=20, column=4,columnspan=2)
		self.e23 = tk.Entry(self.root,width=5,justify='center')
		self.e23.insert(10,self.StrVar23.get())
		self.e23.grid(row=21, column=4,columnspan=2)		

		self.l30 = tk.Label(self.root, text="BSMULTI",anchor='e')
		self.l30.grid(row=20, column=6,columnspan=2)
		self.e24 = tk.Entry(self.root,width=5,justify='center')
		self.e24.insert(10,self.StrVar24.get())
		self.e24.grid(row=21, column=6,columnspan=2)

		self.l30 = tk.Label(self.root, text="CNEO",anchor='e')
		self.l30.grid(row=20, column=8,columnspan=2)
		self.e25 = tk.Entry(self.root,width=5,justify='center')
		self.e25.insert(10,self.StrVar25.get())
		self.e25.grid(row=21, column=8,columnspan=2)

		self.l31 = tk.Label(self.root, text="=============== Other options ===============",justify='center')
		self.l31.grid(row=22, column=0,columnspan=9)
		b5 = tk.Button(self.root, text="NUM opt", bg = "lightgray",command=lambda: self.button_func5(),height = 1,width = 6)
		b5.grid(row=23, column=1,columnspan=3,sticky='e')
		b6 = tk.Button(self.root, text="NODE opt", bg = "lightgray",command=lambda: self.button_func6(),height = 1,width = 6)
		b6.grid(row=23, column=4,columnspan=3)
		b7 = tk.Button(self.root, text="PLOT opt", bg = "lightgray",command=lambda: self.button_func9(),height = 1,width = 6)
		b7.grid(row=23, column=7,columnspan=3,sticky='w')

		self.l32 = tk.Label(self.root, text="ROT_FILE",anchor='e')
		self.l32.grid(row=24, column=0,columnspan=3)
		self.e55 = tk.Entry(self.root,width=20,justify='center')
		self.e55.insert(10,self.StrVar55.get())
		self.e55.grid(row=24, column=3,columnspan=5)	

		b4 = tk.Button(self.root, text="OPEN", bg = "lightgray",command=lambda: self.button_func102(),height = 1,width = 4)
		b4.grid(row=24, column=8)		

		self.l32 = tk.Label(self.root, text="================ Commands ================",justify='center')
		self.l32.grid(row=25, column=0,columnspan=9)

		self.l32 = tk.Label(self.root, text="RUN_NAME",anchor='e')
		self.l32.grid(row=26, column=0,columnspan=3)
		self.e68 = tk.Entry(self.root,width=20,justify='center')
		self.e68.insert(10,self.StrVar68.get())
		self.e68.grid(row=26, column=3,columnspan=5)

		b9 = tk.Button(self.root, text="RUN JASTAB", bg = "lightgray",command=lambda: self.button_func11(),height = 1,width = 10)
		b9.grid(row=30, column=0, columnspan=5)

		b10 = tk.Button(self.root, text="CHECK INP", bg = "lightgray",command=lambda: self.button_func13(),height = 1,width = 10)
		b10.grid(row=29, column=0, columnspan=5)

		b11 = tk.Button(self.root, text="OPEN RESULT", bg = "lightgray",command=lambda: self.button_func12(),height = 1,width = 10)
		b11.grid(row=31, column=0, columnspan=5)

		b11 = tk.Button(self.root, text="EXIT", bg = "lightgray",command=lambda: self.button_func18(self.root,1),height = 1,width = 10)
		b11.grid(row=32, column=0, columnspan=5)

		self.l33 = tk.Label(self.root, text="Use batch",anchor='e')
		self.l33.grid(row=30, column=5,columnspan=3)
		self.c9= tk.Checkbutton(self.root,variable=self.CheckVar9)
		self.c9.grid(row=30, column=8)

		font = tkinter.font.Font(weight='bold',size=10)
		self.l34 = tk.Label(self.root, text="-Ver %s "%(version['jatool']),anchor='w',width=20,justify='left',fg='dodgerblue',font=font)
		self.l34.grid(row=37, column=0,columnspan=10,sticky='sw')

		self.l35 = tk.Label(self.root, text="-%s"%(author['jatool2']),anchor='w',width=40,justify='left',fg='dodgerblue',font=font)
		self.l35.grid(row=38, column=0,columnspan=10,sticky='sw')

		self.l36 = tk.Label(self.root, text="-%s"%(comment['jatool']),anchor='w',width=40,justify='left',fg='magenta',font=font)
		self.l36.grid(row=39, column=0,columnspan=10,sticky='sw')

		self.e69  = tk.OptionMenu(self.root,self.StrVar69,*self.run_list)
		self.e69.config(width=10)
		self.e69.grid(row=31,column=5,columnspan=4)		
		#menu = self.e63.children["menu"]

		#for value in self.run_list:
		self.l100 = tk.Label(self.root, text="NOT READY",fg='red',anchor='center',bg='white',width=15,justify='left')
		self.l100.grid(row=29, column=5,columnspan=4)		

		self.root.mainloop()
		return

	def initialise_vars(self):

		self.iseqdsk = False
		self.rzbdy = np.zeros(shape=(1,2))
		self.rinput = False
		self.iskin = False
		self.isrot = False

		self.t1_close = True
		self.t2_close = True
		self.t3_close = True
		self.t4_close = True
		self.t5_close = True
		self.t6_close = True
		self.t7_close = True

		self.run_name = ''						#str 68
		self.nodelist = node_init					#str 69

		self.use_kin  = True					#check1
		self.use_ped_finder = True				#check2
		self.use_prev = False					#check3
		self.use_li = False						#check4
		self.adjust_prof = True					#check5
		self.hag_core_mod = True				#check6
		self.eli_comp = False					#check7
		self.eli_rot = False					#check8
		self.batch_run = True					#check9
		self.use_elite_a = True					#check10
		self.use_adj_a = True					#check11
		self.use_bilinear = True				#check12
		self.fill_up = False					#check13
		self.use_dia = True						#check14
		self.use_li_iter2 = True				#check15
		self.use_j2 = False                                             #check16
		self.target_bnd = 0.995					#double1
		self.target_i = 4						#double2
		self.target_j = 4						#double3			

		self.equ_name = None					#str1
		self.kin_name = None 					#str2
		self.equ_type = 'KEFIT'					#str3
		self.nw = 8								#str4
		self.pmin = -0.5						#str5
		self.pmax = +0.5						#str6
		self.pedc = 0.98						#str7
		self.nh = 8								#str8
		self.jmin = -0.5						#str9
		self.jmax = +0.5						#str10
		self.equ_code = 'Chease'				#str11
		self.stab_code = 'Mishka'				#str12
		self.mode_n = '5,7,10,15,20'			#str13
		self.bsmodel = 'Hager'					#str14
		self.beta_type = 'Bpol'					#str15

		self.bp = 1.0							#str16
		self.rmag = 1.83						#str17
		self.li = 1.0							#str18
		self.zeff = 2.0							#str19
		self.zimp = 6.0							#str20
		self.amain = 2.0						#str21
		self.aimp = 12.0						#str22
		self.linden = 0.0						#str23
		self.bsmulti = 1.0						#str24
		self.core_neo = 0.01					#str25

		self.apf = 1.0							#str 26
		self.bpf = 0.9							#str 27		
		self.cpf = 1.1							#str 28		
		self.dpf = 1.5							#str 29		
		self.ajf = 0.2							#str 30
		self.bjf = 1.0							#str 31		
		self.cjf = 2.0							#str 32		
		self.djf = 2.0							#str 33			
		self.hag_core_mod_psin = 0.3			#str 34

		self.currenti = 35						#str 35
		self.relax = 0.8						#str 36
		self.ns = 80							#str 37	
		self.nt = 80							#str 38
		self.map_ns = 200						#str 39
		self.map_nt = 200						#str 40

		self.nse = 150							#str 41
		self.nte = 200							#str 42
		self.map_nse = 300						#str 43
		self.map_nte = 512						#str 44		

		self.mis_gridn = 301					#str 45
		self.mis_psistart = 0.75				#str 46
		self.xr1 = 1.0							#str 47
		self.sig1 = 0.07						#str 48
		self.xr2 = 0.9							#str 49
		self.sig2 = 0.1							#str 50

		self.qdel = 0.3							#str 51
		self.eli_gridn = 2000					#str 52
		self.eli_ndist = 50						#str 53
		self.eli_psistart = 0.5					#str 54
		self.eli_rot_file = ''					#str 55

		self.epslon = 1.e-8						#str 56
		self.beta_crit = 1.e-2					#str 57
		self.li_crit = 2.e-2					#str 58
		self.ip_crit = 1.e-5					#str 59
		self.bs_crit = 1.e-5					#str 60
		self.niterl = 10						#str 61

		#japlot opt
		
		self.nqa = 27.7							#str 62
		self.nq_cut_off = 27.7					#str 63
		self.crit_dia = 0.25					#str 64
		self.crit_alf = 0.03					#str 65
		self.cutoff_n = 1						#str 66
		self.grmax = 100						#str 67
				
		self.DoubleVar1.set(self.target_bnd)
		self.DoubleVar2.set(self.target_i)
		self.DoubleVar3.set(self.target_j)		

		self.CheckVar1.set(self.trans_vars(self.use_kin,4))
		self.CheckVar2.set(self.trans_vars(self.use_ped_finder,4))
		self.CheckVar3.set(self.trans_vars(self.use_prev,4))
		self.CheckVar4.set(self.trans_vars(self.use_li,4))
		self.CheckVar5.set(self.trans_vars(self.adjust_prof,4))
		self.CheckVar6.set(self.trans_vars(self.hag_core_mod,4))
		self.CheckVar7.set(self.trans_vars(self.eli_comp,4))
		self.CheckVar8.set(self.trans_vars(self.eli_rot,4))
		self.CheckVar9.set(self.trans_vars(self.batch_run,4))
		self.CheckVar10.set(self.trans_vars(self.use_elite_a,4))
		self.CheckVar11.set(self.trans_vars(self.use_adj_a,4))
		self.CheckVar12.set(self.trans_vars(self.use_bilinear,4))
		self.CheckVar13.set(self.trans_vars(self.fill_up,4))			
		self.CheckVar14.set(self.trans_vars(self.use_dia,4))
		self.CheckVar15.set(self.trans_vars(self.use_li_iter2,4))
		self.CheckVar16.set(self.trans_vars(self.use_j2,4))
		
		self.StrVar1.set(self.trans_vars(self.equ_name,3))	
		self.StrVar2.set(self.trans_vars(self.kin_name,3))			
		self.StrVar3.set(self.trans_vars(self.equ_type,1))
		self.StrVar4.set(self.trans_vars(self.nw,2))
		self.StrVar5.set(self.trans_vars(self.pmin,2))
		self.StrVar6.set(self.trans_vars(self.pmax,2))	
		self.StrVar7.set(self.trans_vars(self.pedc,2))
		self.StrVar8.set(self.trans_vars(self.nh,2))
		self.StrVar9.set(self.trans_vars(self.jmin,2))
		self.StrVar10.set(self.trans_vars(self.jmax,2))

		self.StrVar11.set(self.trans_vars(self.equ_code,1))
		self.StrVar12.set(self.trans_vars(self.stab_code,2))
		self.StrVar13.set(self.trans_vars(self.mode_n,1))
		self.StrVar14.set(self.trans_vars(self.bsmodel,1))
		self.StrVar15.set(self.trans_vars(self.beta_type,2))

		self.StrVar16.set(self.trans_vars(self.bp,2))
		self.StrVar17.set(self.trans_vars(self.rmag,2))
		self.StrVar18.set(self.trans_vars(self.li,2))
		self.StrVar19.set(self.trans_vars(self.zeff,2))
		self.StrVar20.set(self.trans_vars(self.zimp,2))
		self.StrVar21.set(self.trans_vars(self.amain,2))
		self.StrVar22.set(self.trans_vars(self.aimp,2))		
		self.StrVar23.set(self.trans_vars(self.linden,2))
		self.StrVar24.set(self.trans_vars(self.bsmulti,2))	
		self.StrVar25.set(self.trans_vars(self.core_neo,2))

		self.StrVar26.set(self.trans_vars(self.apf,2))
		self.StrVar27.set(self.trans_vars(self.bpf,2))	
		self.StrVar28.set(self.trans_vars(self.cpf,2))
		self.StrVar29.set(self.trans_vars(self.dpf,2))
		self.StrVar30.set(self.trans_vars(self.ajf,2))
		self.StrVar31.set(self.trans_vars(self.bjf,2))	
		self.StrVar32.set(self.trans_vars(self.cjf,2))
		self.StrVar33.set(self.trans_vars(self.djf,2))
		self.StrVar34.set(self.trans_vars(self.hag_core_mod_psin,2))
		self.StrVar35.set(self.trans_vars(self.currenti,2))
		self.StrVar36.set(self.trans_vars(self.relax,2))
		self.StrVar37.set(self.trans_vars(self.ns,2))	
		self.StrVar38.set(self.trans_vars(self.nt,2))
		self.StrVar39.set(self.trans_vars(self.map_ns,2))
		self.StrVar40.set(self.trans_vars(self.map_nt,2))
		self.StrVar41.set(self.trans_vars(self.nse,2))	
		self.StrVar42.set(self.trans_vars(self.nte,2))
		self.StrVar43.set(self.trans_vars(self.map_nse,2))
		self.StrVar44.set(self.trans_vars(self.map_nte,2))

		self.StrVar45.set(self.trans_vars(self.mis_gridn,2))
		self.StrVar46.set(self.trans_vars(self.mis_psistart,2))
		self.StrVar47.set(self.trans_vars(self.xr1,2))	
		self.StrVar48.set(self.trans_vars(self.sig1,2))
		self.StrVar49.set(self.trans_vars(self.xr2,2))
		self.StrVar50.set(self.trans_vars(self.sig2,2))

		self.StrVar51.set(self.trans_vars(self.qdel,2))

		self.StrVar52.set(self.trans_vars(self.eli_gridn,2))
		self.StrVar53.set(self.trans_vars(self.eli_ndist,2))		
		self.StrVar54.set(self.trans_vars(self.eli_psistart,2))
		self.StrVar55.set(self.trans_vars(self.eli_rot_file,1))

		self.StrVar56.set(self.trans_vars(self.epslon,2))
		self.StrVar57.set(self.trans_vars(self.beta_crit,2))
		self.StrVar58.set(self.trans_vars(self.li_crit,2))
		self.StrVar59.set(self.trans_vars(self.ip_crit,2))
		self.StrVar60.set(self.trans_vars(self.bs_crit,2))
		self.StrVar61.set(self.trans_vars(self.niterl,2))		

		self.StrVar62.set(self.trans_vars(self.nqa,2))
		self.StrVar63.set(self.trans_vars(self.nq_cut_off,2))
		self.StrVar64.set(self.trans_vars(self.crit_dia,2))
		self.StrVar65.set(self.trans_vars(self.crit_alf,2))
		self.StrVar66.set(self.trans_vars(self.cutoff_n,2))
		self.StrVar67.set(self.trans_vars(self.grmax,2))
	
		self.StrVar68.set(self.trans_vars(self.run_name,1))
		self.StrVar69.set('-search-')

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

	def __init__(self):

		self.run_list = []
		self.run_list = self.open_result_list()
		self.currdir = os.getcwd()
		return

if __name__ == "__main__":

	import gui_jastab

	print(' ---------------------------------------------------------------')
	print('||                Python based JA tool Ver %3s                 ||'%version['jatool'])
	print('||             Linear PBM stability analysis tool              ||')
	print('%s'%author['jatool'])
	print(' ---------------------------------------------------------------\n')

	gja = gui_jastab.gjastab()	
	gja.root = tk.Tk()
	gja.root.title('G-JATOOL')
	gja.declare_vars()
	gja.gui_jastab()
