#!/usr/local/anaconda3/bin/python3

import os,sys
import numpy as np
from shutil import copyfile
from scipy.optimize import curve_fit
import ch_tool 
import eqdsk
from scipy.interpolate import interp1d
import time
import matplotlib.pyplot as plt
from progress import update_progress
from batch_run import *
from exec_dirs import pedstab_dir,python3_exec,node_default,node_init

currdir = os.getcwd()
class pedscan:

	def make_scan_script(self,filename,nw,nh):

		w = np.linspace(self.wmin,self.wmax,self.nw) * (self.eped_prof[0,4] + self.eped_prof[1,4]) * 0.5
	
		log_e = self.work_dir+'/'+str(nw)+'_'+str(nh)+'/run.e'
		log_o = self.work_dir+'/'+str(nw)+'_'+str(nh)+'/run.o'
		command = 'cd '+self.work_dir+'/'+str(nw)+'_'+str(nh)+'\n'
		command = command + python3_exec + ' ' + pedstab_dir + ' ' + str(nw) + ' ' + str(nh) + ' ' + str(w[nw]) + ' \n'

		make_batch_script(filename,self.node_list,log_e,log_o,command,'PEDSCAN')

		return
		
	def submit_script(self):
	
		batch_num = []
	
		script_dir = self.work_dir + '/script'
		try:
			os.mkdir(script_dir)
		except:
			pass
	
		for i in range(self.nw):
		
			for j in range(self.nh):
			
				filename = script_dir + '/script_'+str(i)+'_'+str(j)
				
				self.make_scan_script(filename,i,j)
				runid = submit_batch_script(filename)
				batch_num.append(str(runid))
		
		print('Scripts are submitted')
		print('Collect stability data')		
		batch_chk = 1

		bat_count0 = 1
		while ( bat_count0 > 0):
			output = batch_script_status()
			bat_count = 0
			for i in range(self.nw):
				for j in range(self.nh):
					index1 = (self.nh)*i + j
					batch_chk = output.find(batch_num[index1])
					if (batch_chk > -1):
						bat_count = bat_count + 1;

			if not (bat_count == bat_count0):
				bat_count0 = bat_count
				update_progress(1.0 - bat_count0/self.nw/self.nh)
			if (bat_count0 > 0):
				time.sleep(30)

		print('Job scripts are finished')
		print('Finalizing process --> Collect data')
		
		return
		
	def check_inp(self):
		os.chdir(self.work_dir)
		if len(self.eped_file.split('/')) == 0.:	self.eped_file = currdir + '/' + self.eped_file
		if self.ti_file == None:	
			self.ext_ti = Fasle
		else:	
			if len(self.ti_file.split('/')) == 0.:	self.ti_file = currdir + '/' + self.ti_file
		if len(self.eqdsk_file.split('/')) == 0 :	self.eqdsk_file = currdir + '/' + self.eqdsk_file
		if self.rot_file == None:	
			self.use_rot = False
		else:	
			if len(self.rot_file.split('/')) == 0:	self.rot_file = currdir + '/' + self.rot_file

		if not (os.path.isfile(self.eped_file)):
			print('No EPED prof')
			exit()
		print('EPED PROF exits')
		
		if (self.ext_ti):
			if not (os.path.isfile(self.ti_file)):
				print('No Ti prof')
				exit()
			print('TI PROF exits')
		
		if not (os.path.isfile(self.eqdsk_file)):	
			print('No EQDSK file')
			exit()

		print('EQDSK file exits')
		
		if (self.use_rot):
			if not (os.path.isfile(self.rot_file)):
				print('No rot prof exits')
				exit()
			print('ROT PROF exits')			
		os.chdir(currdir)
		return
	
	def read_eped_prof(self,filename):
	
		self.ch.read_kinprof_eped_fun(filename)
		
		self.eped_prof = np.copy(self.ch.eped_prof)
	
		self.ncore = self.eped_prof[0,1] + self.eped_prof[0,0]*(np.tanh(1)+1.0) + self.eped_prof[0,2]
		self.nsep = self.eped_prof[0,1]		

		return
	
	def read_kin_prof(self,filename):

		if self.ti_file_type == 1:
			self.tik = self.ch.read_kinprof_kprofile_fun(filename)/1.e3
			self.tif = interp1d(self.psin,self.tik,'cubic')
		else:
			self.ch.read_kinprof_chease(self.ti_file)
			self.tif = interp1d(self.ch.psink,self.ch.tik,'cubic')
		
		return

	def eped_fun(self,x,a7,a8):
		
		val = self.temp[1]
		val = val + self.temp[0]*(np.tanh(2.0/self.temp[4]*(1.0-self.temp[3]))-np.tanh(2.0/self.temp[4]*(x-self.temp[3]))) 
		val = val + self.temp[2] * ((1.0 - (x/self.temp[5])**(1.01+abs(a7)))**(1.01+abs(a8))) * 0.5 * (1.0+np.sign(self.temp[5]-x))

		return val
		
	def eped_fit(self,eped_prof1,eped_prof2):
		
		num = 60
		psin = np.linspace(0.0,max(0.5,eped_prof2[5]-0.20),num)
		prof = np.copy(psin)
		self.temp = np.copy(eped_prof1)
		
		for i in range(num):
			prof[i] = self.eped_fun(psin[i],eped_prof1[6]-1.01,eped_prof1[7]-1.01)

		self.temp = np.copy(eped_prof2)
		count = 0
		pp = [0.,1.0]
		while count < 8:
			popt, pcov = curve_fit(self.eped_fun,psin,prof,p0=pp,maxfev=300000)
		#	print(popt,count,psin[-1])
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
			popt[0] = eped_prof1[6] - 1.01
			popt[1] = eped_prof1[7] - 1.01
			print('Fitting failed')
			print(popt[0],popt[1])
	
		return (abs(popt[0])+1.01,abs(popt[1])+1.01)

	def make_eped_prof(self,eped_prof,psin):
	
		num = len(psin)
		
		te = np.zeros(num)
		ti = np.zeros(num)
		ne = np.zeros(num)
		
		for i in range(num):
		
			psint = psin[i]
			
			ne[i] = eped_prof[0,1] + eped_prof[0,0]*(np.tanh(2.0/eped_prof[0,4]*(1.0-eped_prof[0,3]))-np.tanh(2.0/eped_prof[0,4]*(psint-eped_prof[0,3])))
			if (psint < eped_prof[0,5]):
				ne[i] = ne[i] + eped_prof[0,2] * ((1.0 - (psint/eped_prof[0,5])**eped_prof[0,6])**eped_prof[0,7])
				
			te[i] = eped_prof[1,1] + eped_prof[1,0]*(np.tanh(2.0/eped_prof[1,4]*(1.0-eped_prof[1,3]))-np.tanh(2.0/eped_prof[1,4]*(psint-eped_prof[1,3])))
			if (psint < eped_prof[1,5]):
				te[i] = te[i] + eped_prof[1,2] * ((1.0 - (psint/eped_prof[1,5])**eped_prof[1,6])**eped_prof[1,7])
			
			ti[i] = eped_prof[2,1] + eped_prof[2,0]*(np.tanh(2.0/eped_prof[2,4]*(1.0-eped_prof[2,3]))-np.tanh(2.0/eped_prof[2,4]*(psint-eped_prof[2,3])))
			if (psint < eped_prof[2,5]):
				ti[i] = ti[i] + eped_prof[2,2] * ((1.0 - (psint/eped_prof[2,5])**eped_prof[2,6])**eped_prof[2,7])
				
		return (ne,te,ti)
		
	def make_prof(self,nw,nh,final=False,onlyprof=False):

		eped_prof = np.copy(self.eped_prof)
		
		
		ncore = eped_prof[0,1] + eped_prof[0,0]*(np.tanh(1)+1.0) + eped_prof[0,2]
		tcoree = eped_prof[1,1] + eped_prof[1,0]*(np.tanh(1)+1.0) + eped_prof[1,2]
		tcorei = eped_prof[2,1] + eped_prof[2,0]*(np.tanh(1)+1.0) + eped_prof[2,2]

	
		width = (self.wmin + (self.wmax-self.wmin)*float(nw)/float(self.nw-1))
		height = (self.hmin + (self.hmax-self.hmin)*float(nh)/float(self.nh-1))
		height_n = height
		width_n  = width
		if os.getenv("PEDSCAN_TYPE") == 'DENFIX': height_n = 1.
		elif os.getenv("PEDSCAN_TYPE") == 'TWOWID':
			width_n = height
			height  = 1.
			height_n= 1.
		
		eped_prof[0,0] = eped_prof[0,0] * height_n
		eped_prof[1,0] = eped_prof[1,0] * height
		if not self.fixed_ti:	eped_prof[2,0] = eped_prof[2,0] * height
		
		eped_prof[0,2] = ncore - eped_prof[0,1] - eped_prof[0,0]*(np.tanh(1)+1.0)
		eped_prof[1,2] = tcoree - eped_prof[1,1] - eped_prof[1,0]*(np.tanh(1)+1.0)
		if not self.fixed_ti:	eped_prof[2,2] = tcorei - eped_prof[2,1] - eped_prof[2,0]*(np.tanh(1)+1.0)
		
		eped_prof[0,4] = eped_prof[0,4] * width_n
		eped_prof[1,4] = eped_prof[1,4] * width
		if not self.fixed_ti:	eped_prof[2,4] = eped_prof[2,4] * width
		
		eped_prof[0,3] = 1.0 - 0.5*eped_prof[0,4]
		eped_prof[1,3] = 1.0 - 0.5*eped_prof[1,4]
		if not self.fixed_ti:	eped_prof[2,3] = 1.0 - 0.5*eped_prof[2,4]
		
		eped_prof[0,5] = 1.0 - eped_prof[0,4]
		eped_prof[1,5] = 1.0 - eped_prof[1,4]
		if not self.fixed_ti:	eped_prof[2,5] = 1.0 - eped_prof[2,4]
		
		eped_prof[0,6], eped_prof[0,7] = self.eped_fit(self.eped_prof[0,:],eped_prof[0,:])
		eped_prof[1,6], eped_prof[1,7] = self.eped_fit(self.eped_prof[1,:],eped_prof[1,:])
		if not self.fixed_ti:	eped_prof[2,6], eped_prof[2,7] = self.eped_fit(self.eped_prof[2,:],eped_prof[2,:])
		
		self.psin2 = np.linspace(0,1.2,301)
	
		ne2,te2,ti2 = self.make_eped_prof(eped_prof,self.psin2)
		ni2 = np.copy(ne2 * (1.0 - (self.zeff-1.0)/self.zimp))
		
		nef = interp1d(self.psin2,ne2,'cubic')
		tef = interp1d(self.psin2,te2,'cubic')
		nif = interp1d(self.psin2,ni2,'cubic')
		tif = interp1d(self.psin2,ti2,'cubic')

		ne = nef(self.psin)
		te = tef(self.psin)
		ni = nif(self.psin)
		ti = tif(self.psin)

		self.ne = ne
		self.te = te

		self.ne2 = ne2
		self.te2 = te2

		if not (self.ext_ti):
			self.ti = ti
		else:
			self.ti = self.tif(self.psin)

		if onlyprof:	return

		rho = np.linspace(0,1.0,101)
		self.ch.read_rho_psi_R(self.currdir+'/CHEASE/RHO_PSI_R')		
		psinf = interp1d(self.ch.rho_map,self.ch.psi_map,'cubic')

		psin = psinf(rho)

		ne2 = nef(psin)
		te2 = tef(psin)

		if (final):
			try:
				os.mkdir('PROFILES')
			except:
				pass

		if (final):	
			f4 = open('PROFILES/NE.dat_mod','w')
			f5 = open('PROFILES/TE.dat_mod','w')
			f7 = open('PROFILES/chease_kinprof_mod','w')
		else:
			f4 = open('NE.dat','w')
			f5 = open('TE.dat','w')
			f7 = open('chease_kinprof','w')

		f7.write('%i\n'%401)
		f7.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.zeff,self.zimp,self.amain,self.aimp))

		for i in range(401):

			if not (self.ext_ti):
				f7.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.psin[i],te[i],ne[i],ti[i],ni[i]))
			else:
				f7.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.psin[i],te[i],ne[i],self.tif(self.psin[i]),ni[i]))
		f7.close()

		if not (self.ext_ti):
			if (final):
				f6 = open('PROFILES/TI.dat_mod','w')
			else:
				f6 = open('TI.dat','w')

			f6.write('R(m)       Z(m)       TI(eV)\n')
		
		f4.write('R(m)       Z(m)       NE(1E18 m-3)\n')
		f5.write('R(m)       Z(m)       TE(eV)\n')
		
		for i in range(101):

			f4.write('%9.6f\t%9.6f\t%9.6f\n'%(rho[i],psin[i],ne2[i]*1.e1))
			f5.write('%9.6f\t%9.6f\t%9.6f\n'%(rho[i],psin[i],te2[i]*1.e3))
		
			if not (self.ext_ti):
				f6.write('%9.6f\t%9.6f\t%9.6f\n'%(rho[i],psin[i],te2[i]*1.e3))
		if not (self.ext_ti):
			f6.close()
		
		f4.close()
		f5.close()
		
		return
		
	def make_stab_input(self,filename='stab_input'):
	
		f4 = open(filename,'w')
		f4.write('RUNSTAB = %s \n'%self.run_stab)
		f4.write('NS = %i \n'%self.ns)
		f4.write('NT = %i \n'%self.nt)
		f4.write('MAP_NS = %i \n'%self.map_nse)
		f4.write('MAP_NT = %i \n'%self.map_nte)
		f4.write('NSE = %i \n'%self.ns)
		f4.write('NTE = %i \n'%self.nt)
		f4.write('MAP_NSE = %i \n'%self.map_nse)
		f4.write('MAP_NTE = %i \n'%self.map_nte)
		f4.write('\n')	

		f4.write('MIS_GRIDN = %i \n'%self.mis_gridn)
		f4.write('MIS_PSISTART = %f \n'%self.mis_psistart)
		f4.write('XR1 = %f \n'%self.xr1)
		f4.write('SIG1 = %f \n'%self.sig1)
		f4.write('XR2 = %f \n'%self.xr2)
		f4.write('SIG2 = %f \n'%self.sig2)
		f4.write('ELI_GRIDN = %i \n'%self.eli_gridn)
		f4.write('ELI_PSISTART = %f \n'%self.eli_psistart)
		f4.write('ELI_NDIST = %i \n'%self.eli_ndist)
		f4.write('\n')	

		f4.write('USE_COMP = %s\n'%self.use_comp)
		f4.write('USE_ROT = %s\n'%self.use_rot)

		f4.write('QDEL = %f \n'%self.qdelfix)
		f4.write('EPSLON = %e \n'%self.epslon)
		f4.write('moden = %s \n'%self.mode_line)
		
		f4.close()
		
		return

	def make_chease_opt(self,filename='chease_opt',pedw=0.05):
	
		f = open(filename,'w')

		f.write('!-- Plasma 0-d.\n')
		f.write('EQDSK = geqdsk \n')
		f.write('BND_PSIN = %f'%self.target_psin)
		f.write('\n')

		f.write('ZIMP  = %f \n'%self.zimp)
		f.write('ZEFF  = %f \n'%self.zeff)
		f.write('AMAIN = %f \n'%self.amain)
		f.write('AIMP  = %f \n'%self.aimp)
		if (self.beta_criterion == 1):	self.beta_critval = self.bp
		else:	self.beta_critval = self.rmag
		f.write('Beta_val = %f \n'%self.beta_critval)
		f.write('LITARGET = %f \n'%self.li)	
		f.write(' \n')
		f.write('!-- Fast ion profile option \n')
		f.write('APF = %f \n'%self.apf)
		f.write('BPF = %f \n'%self.bpf)
		f.write('CPF = %f \n'%self.cpf)
		f.write('DPF = %f \n'%self.dpf)
		f.write(' \n')
		f.write('!-- Current profile option \n')

		f.write('USE_NEO = %s \n'%self.use_neo)
		f.write('USE_HAGER = %s \n'%self.use_hager)
		f.write('USE_CHANG = %s \n'%self.use_chang)
		f.write('HAG_CORE_MOD = %s \n'%self.hag_core_mod)
		f.write('HAG_CORE_MOD_PSIN = %f \n'%self.hag_core_mod_psin)
		f.write('Core_neo = %f \n'%self.core_neo)
		f.write('BSMULTI = %f \n'%self.bsmulti)
		f.write('AJF = %f \n'%self.ajf)
		f.write('BJF = %f \n'%self.bjf)
		f.write('CJF = %f \n'%self.cjf)
		f.write('DJF = %f \n'%self.djf)
		f.write('Current_ITERN = %f \n'%self.niterc)
		f.write('RELAX = %f \n'%self.relax)
		f.write(' \n')

		f.write('!-- convergence option\n')
		f.write('EPSILON = %e \n'%self.epslon)
		f.write('Beta_crit = %f \n'%self.beta_crit)
		f.write('Li_crit = %f \n'%self.li_crit)
		f.write('Ip_crit = %f \n'%self.ip_crit)
		f.write('Bs_crit = %f \n'%self.bs_crit)
		f.write('Li_ITERN = %f \n'%self.niterl)
		f.write(' \n')

		f.write('!-- grid option\n')
		f.write('NS = %i \n'%self.ns)
		f.write('NT = %i \n'%self.nt)
		f.write('MAP_NS = %i \n'%self.map_ns)
		f.write('MAP_NT = %i \n'%self.map_nt)
		f.write(' \n')

		f.write('!-- default option (Do not adjust it!) \n')
		line = 'RUN_MODE = NORMAL \n'
		if (self.use_li and self.use_li2):
			line = 'RUN_MODE = EPED3 \n'
		elif (self.use_li):
			line = 'RUN_MODE = EPED2 \n'			
		f.write(line)
		f.write('Beta_criterion_type = %i \n'%self.beta_criterion)
		f.write('ADJUST_PROF = %s \n'%self.adjust_prof)
		f.write('kinetic_profile_type = 1 \n')
		f.write('chease_kinetic_file = chease_kinprof \n')
		f.write('NIDEAL = 8\n')
		f.write('PED_WIDTH = %f \n'%pedw)
		f.close()

		return
		
	def make_directories(self):
	
		currdir2 = os.getcwd()
	
		try:
			os.mkdir(self.work_dir)
		except:
			print('WORK DIR EXISTS')
			
		os.chdir(self.work_dir)
		
		for i in range(self.nw):
		
			for j in range(self.nh):
			
				try:
					os.mkdir(str(i)+'_'+str(j))
					
				except:
					pass
		
				os.chdir(str(i)+'_'+str(j))
				copyfile(self.eqdsk_file,'geqdsk')
				self.make_prof(i,j)

				ww = (self.wmin + (self.wmax-self.wmin)*float(i)/float(self.nw-1))
				w1 = self.eped_prof[0,4]*ww
				w2 = self.eped_prof[1,4]*ww

				self.make_chease_opt('chease_opt',0.5*(w1+w2))
				if (self.use_rot):
					copyfile(self.rot_file,'chease_vtor')
				self.make_stab_input()
				os.chdir('../')

		os.chdir(currdir2)
		
		return
	
	def initialise_vars(self):
	
		self.ch = ch_tool.chease()
		self.currdir = os.getcwd()	
		self.psin = np.linspace(0,1.0,401)

		self.run_name = 'ped_scan'
		self.node_list = node_default					#str 69

		self.use_li = False						#check4
		self.adjust_prof = True					#check5
		self.hag_core_mod = True				#check6
		self.use_comp = False					#check7
		self.use_rot = False					#check8
		self.batch_run = False					#check9
		self.use_bilinear = True				#check12
		self.use_dia = True						#check14
		self.use_li2 = False					#check15
		self.fixed_ti = True
		self.ext_ti = True
		self.highq = False

		self.use_hager = True
		self.use_chang = True
		self.use_neo = True

		self.target_psin = 0.995				#double1
		self.target_i = 1						#double2
		self.target_j = 1						#double3			
		self.target_i2 = 1						#double2
		self.target_j2 = 1						#double3		

		self.eqdsk_file = None					#str1
		self.eped_file  = None 					#str2
		self.ti_file    = None					
		self.ti_file_type = 1
		self.nw = 8								#str4
		self.wmin = -0.5						#str5
		self.wmax = +0.5						#str6
		self.nh = 8								#str8
		self.hmin = -0.5						#str9
		self.hmax = +0.5						#str10

		self.run_stab = 'mishka'				#str12
		self.mode_line = '5,7,10,15,20'			#str13

		self.beta_criterion = 2	
		self.beta_critval = 0.0

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

		self.niterc = 35						#str 35
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

		self.qdelfix = 0.3						#str 51
		self.eli_gridn = 2000					#str 52
		self.eli_ndist = 50						#str 53
		self.eli_psistart = 0.5					#str 54
		self.rot_file 	= None					#str 55

		self.epslon = 1.e-8						#str 56
		self.beta_crit = 1.e-2					#str 57
		self.li_crit = 2.e-2					#str 58
		self.ip_crit = 1.e-5					#str 59
		self.bs_crit = 1.e-5					#str 60
		self.niterl = 10						#str 61

		#japlot opt
		
		self.nqa = 27.7							#str 62
		self.gr_crit1 = 0.03					#str 64
		self.gr_crit2 = 0.25					#str 65
		self.ncut = 1							#str 66
		self.grmax = 100						#str 67
		self.grmulti = 1.0
		self.targetn = 12

		self.use_kin_prof = True
#		self.work_dir = currdir+'/' + self.run_name
	
		return

	def read_namelist(self,filename):
	
		f4 = open(filename,'r')
		while True:
		
			line = f4.readline()
			if not line: break

			if (line.lower().find('node_list') > -1):
				self.node_list = line.split('=')[1]	

			self.run_name = self.ch.read_namelist_str(line,'RUN_NAME',self.run_name,3)

			self.use_li = self.ch.read_namelist_str(line,'USE_li',self.use_li,4)
			self.use_li2 = self.ch.read_namelist_str(line,'USE_li2',self.use_li2,4)
			self.adjust_prof = self.ch.read_namelist_str(line,'Adjust_prof',self.adjust_prof,4)
			self.hag_core_mod = self.ch.read_namelist_str(line,'HAG_CORE_MOD',self.hag_core_mod,4)
			self.use_comp = self.ch.read_namelist_str(line,'USE_COMP',self.use_comp,4)
			self.use_rot = self.ch.read_namelist_str(line,'USE_ROT',self.use_rot,4)
			self.fixed_ti = self.ch.read_namelist_str(line,'FIX_TI',self.fixed_ti,4)
			self.ext_ti = self.ch.read_namelist_str(line,'EXT_TI',self.ext_ti,4)			

			self.use_hager = self.ch.read_namelist_str(line,'USE_HAGER',self.use_hager,4)
			self.use_neo   = self.ch.read_namelist_str(line,'USE_NEO',self.use_neo,4)
			self.use_chang   = self.ch.read_namelist_str(line,'USE_CHANG',self.use_neo,4)			

			self.target_psin = self.ch.read_namelist_str(line,'BND_PSIN',self.target_psin,2)
	
			self.eqdsk_file = self.ch.read_namelist_str(line,'EQDSK',self.eqdsk_file,3)
			self.eped_file = self.ch.read_namelist_str(line,'eped_prof_file',self.eped_file,3)
			self.ti_file = self.ch.read_namelist_str(line,'Ti_file',self.ti_file,3)
			self.ti_file_type = self.ch.read_namelist_str(line,'TI_file_type',self.ti_file_type,1)
			self.nw = self.ch.read_namelist_str(line,'scan_nw',self.nw,1)
			self.nh = self.ch.read_namelist_str(line,'scan_nh',self.nh,1)
			self.wmin = self.ch.read_namelist_str(line,'wmin',self.wmin,2)
			self.wmax = self.ch.read_namelist_str(line,'wmax',self.wmax,2)
			self.hmin = self.ch.read_namelist_str(line,'hmin',self.hmin,2)
			self.hmax = self.ch.read_namelist_str(line,'hmax',self.hmax,2)			

			self.run_stab = self.ch.read_namelist_str(line,'Stab_code',self.run_stab,3)			
			if (line.lower().find('moden') > -1):
				self.mode_line = line.split('=')[1]

			self.beta_criterion = self.ch.read_namelist_str(line,'Beta_criterion_type',self.beta_criterion,1)			
			self.beta_critval = self.ch.read_namelist_str(line,'Beta_val',self.beta_critval,2)
			
			self.bp = self.ch.read_namelist_str(line,'BPOL',self.bp,2)
			self.rmag = self.ch.read_namelist_str(line,'RMAG',self.rmag,2)
			self.li = self.ch.read_namelist_str(line,'LI',self.li,2)
			self.zeff = self.ch.read_namelist_str(line,'ZEFF',self.zeff,2)
			self.zimp = self.ch.read_namelist_str(line,'ZIMP',self.zimp,2)
			self.amain = self.ch.read_namelist_str(line,'AMAIN',self.amain,2)
			self.aimp = self.ch.read_namelist_str(line,'AIMP',self.aimp,2)
			self.linden = self.ch.read_namelist_str(line,'LINEDEN',self.linden,2)
			self.bsmulti = self.ch.read_namelist_str(line,'BSMULTI',self.bsmulti,2)			
			self.core_neo = self.ch.read_namelist_str(line,'CORENEO',self.core_neo,2)			
			
			self.apf = self.ch.read_namelist_str(line,'APF',self.apf,2)
			self.bpf = self.ch.read_namelist_str(line,'BPF',self.bpf,2)
			self.cpf = self.ch.read_namelist_str(line,'CPF',self.cpf,2)
			self.dpf = self.ch.read_namelist_str(line,'DPF',self.dpf,2)

			self.ajf = self.ch.read_namelist_str(line,'AJF',self.ajf,2)
			self.bjf = self.ch.read_namelist_str(line,'BJF',self.bjf,2)
			self.cjf = self.ch.read_namelist_str(line,'CJF',self.cjf,2)
			self.djf = self.ch.read_namelist_str(line,'DJF',self.djf,2)
			self.hag_core_mod_psin = self.ch.read_namelist_str(line,'HAG_CORE_MOD_PSIN',self.hag_core_mod_psin,4)

			self.niterc = self.ch.read_namelist_str(line,'Current_ITERN',self.niterc,1)				
			self.relax = self.ch.read_namelist_str(line,'RELAX',self.relax,2)				
			self.ns = self.ch.read_namelist_str(line,'NS',self.ns,1)
			self.nt = self.ch.read_namelist_str(line,'NT',self.nt,1)
			self.map_ns = self.ch.read_namelist_str(line,'MAP_NS',self.map_ns,1)
			self.map_nt = self.ch.read_namelist_str(line,'MAP_NT',self.map_nt,1)

			self.nse = self.ch.read_namelist_str(line,'NSE',self.nse,1)
			self.nte = self.ch.read_namelist_str(line,'NTE',self.nte,1)
			self.map_nse = self.ch.read_namelist_str(line,'MAP_NSE',self.map_nse,1)
			self.map_nte = self.ch.read_namelist_str(line,'MAP_NTE',self.map_nte,1)

			self.mis_gridn = self.ch.read_namelist_str(line,'MIS_GRIDN',self.mis_gridn,1)
			self.mis_psistart = self.ch.read_namelist_str(line,'MIS_PSISTART',self.mis_psistart,2)
			self.xr1 = self.ch.read_namelist_str(line,'XR1',self.xr1,2)
			self.sig1 = self.ch.read_namelist_str(line,'SIG1',self.sig1,2)
			self.xr2 = self.ch.read_namelist_str(line,'XR2',self.xr2,2)
			self.sig2 = self.ch.read_namelist_str(line,'SIG2',self.sig2,2)
						
			self.qdelfix = self.ch.read_namelist_str(line,'delQ',self.qdelfix,2)
			self.eli_gridn = self.ch.read_namelist_str(line,'ELI_GRIDN',self.eli_gridn,2)
			self.eli_ndist = self.ch.read_namelist_str(line,'ELI_NDIST',self.eli_ndist,2)
			self.eli_psistart = self.ch.read_namelist_str(line,'ELI_PSISTART',self.eli_psistart,2)
			self.rot_file = self.ch.read_namelist_str(line,'rot_file',self.rot_file,3)

			self.epslon = self.ch.read_namelist_str(line,'EPSLON',self.epslon,2)	
			self.beta_crit = self.ch.read_namelist_str(line,'Beta_crit',self.beta_crit,2)
			self.li_crit = self.ch.read_namelist_str(line,'li_crit',self.li_crit,2)
			self.ip_crit = self.ch.read_namelist_str(line,'ip_crit',self.ip_crit,2)
			self.bs_crit = self.ch.read_namelist_str(line,'bs_crit',self.bs_crit,2)
			self.niterl = self.ch.read_namelist_str(line,'NITERL',self.niterl,1)
			
			self.nqa = self.ch.read_namelist_str(line,'NQA',self.nqa,1)	
			self.gr_crit1 = self.ch.read_namelist_str(line,'GRCRIT1',self.gr_crit1,2)
			self.gr_crit2 = self.ch.read_namelist_str(line,'GRCRIT2',self.gr_crit2,2)			
			self.ncut = self.ch.read_namelist_str(line,'NCUT',self.ncut,1)
			self.grmax = self.ch.read_namelist_str(line,'GRMAX',self.grmax,2)
			self.grmulti = self.ch.read_namelist_str(line,'GRMULTI',self.grmulti,2)
			self.targetn = self.ch.read_namelist_str(line,'targetn',self.targetn,1)
				
			self.use_bilinear = self.ch.read_namelist_str(line,'USE_BILINEAR',self.use_bilinear,4)
			self.use_dia = self.ch.read_namelist_str(line,'USE_dia',self.use_dia,4)
			self.target_i = self.ch.read_namelist_str(line,'TARGETI',self.target_i,2)
			self.target_j = self.ch.read_namelist_str(line,'TARGETJ',self.target_j,2)

			self.highq = self.ch.read_namelist_str(line,'HIGHQ',self.highq,4)
			
	
		if (self.run_stab.lower() == 'mishka'):
			self.use_rot = False
			self.use_comp = False

		if (self.ext_ti):
			self.fixed_ti = True

		self.work_dir = currdir+'/' + self.run_name			
		if self.ingfitp: self.work_dir = os.getcwd()
		return
		
	def collect_data(self):
		os.chdir(self.work_dir)
		f = open(self.work_dir+'/result.dat','w')
		
		for i in range(self.nw):
		
			for j in range(self.nh):
			
				os.chdir(str(i)+'_'+str(j))
					
				try:
					f2 = open('growth_list','r')
					line2 = ''
					while True:
						line = f2.readline()
						if not line: break
						f.write(line)
						line2 = line2 + line
					f2.close()
				except:
					f.write(line2)
				os.chdir('../')
				
		f.close()
		
		return

	def clear_data(self):
		
		os.chdir(self.work_dir)
		print('Remove temporary data...')		
		for i in range(self.nw):
			for j in range(self.nh):

				update_progress((self.nh*i+j+1)/self.nw/self.nh)
				os.chdir(str(i)+'_'+str(j))
				try:	os.remove('NMISHKA')
				except:	pass
				try:	os.remove('geqdsk')
				except:	pass
				try:	os.remove('CHEASE/EQDSK_COCOS_02.OUT')
				except:	pass
				try:	os.remove('CHEASE/EQDSK_COCOS_02_POS.OUT')
				except: pass
				try:	os.remove('CHEASE/NELITE')
				except:	pass
				try:	os.remove('CHEASE/NELITE2')
				except:	pass
				try:	os.remove('CHEASE/NDES')
				except:	pass
				try:	os.remove('CHEASE/NOUT')
				except:	pass
				try:	os.remove('MISHKA/CASPLOT')
				except:	pass
				try:	os.remove('MISHKA/eig_str')
				except:	pass
				try:	os.remove('MISHKA/fort.12')
				except:	pass
				try:	os.remove('MISHKA/fort.22')
				except:	pass
	
				os.chdir('../')

		return
	
	def analyse_data(self,filename='result.dat'):
	
		moden = self.mode_line.split(',')
		
		nn = len(moden)
		self.moden = np.zeros(nn)
		for i in range(nn):
			self.moden[i] = int(moden[i])

		self.nn = nn
		w = np.linspace(self.wmin,self.wmax,self.nw)
		h = np.linspace(self.hmin,self.hmax,self.nh)
		
		self.ww, self.hh = np.meshgrid(w,h)
		
		self.gr1 = np.zeros(shape=(self.nw,self.nh,nn))
		self.gr2 = np.zeros(shape=(self.nw,self.nh,nn))
		self.qa  = np.zeros(shape=(self.nw,self.nh,nn))
		
		self.gr3 = np.zeros(shape=(self.nw,self.nh))
		self.gr4 = np.zeros(shape=(self.nw,self.nh))
		
		self.n1 = np.zeros(shape=(self.nw,self.nh))
		self.n2 = np.zeros(shape=(self.nw,self.nh))
		
		currdir = os.getcwd()
		os.chdir(self.work_dir)
		try:	f = open(self.work_dir+'/'+filename,'r')
		except:	f = open(filename,'r')
	
		for i in range(self.nw):
			for j in range(self.nh):
				for nnc in range(nn):
					line = f.readline().split()
					try:
						if not (int(float(line[-1])) == -1):
							self.gr1[i,j,nnc] = float(line[9])
							self.gr2[i,j,nnc] = float(line[10])*self.grmulti
						else:
							self.gr1[i,j,nnc] = np.nan
							sefl.gr2[i,j,nnc] = np.nan
						
					except:
						pass

					if (np.isnan(self.gr1[i,j,nnc]) and int(float(line[-1])) > -1):
						self.gr1[i,j,nnc] = 0.0
						self.gr2[i,j,nnc] = 0.0
						
					self.qa[i,j,nnc] = float(line[12])
				
		f.close()	
		
		for i in range(self.nw):
			for j in range(self.nh):
				for k in range(nn):
					n_crit = self.nqa / self.qa[i,j,k]
					n = int(moden[k])
					if self.use_bilinear:	self.gr2[i,j,k] = self.gr2[i,j,k] * n / min(n_crit,n)
					if (n < self.ncut):
						self.gr2[i,j,k] = 0.0
					
		for i in range(self.nw):
			for j in range(self.nh):
			
				nmax1 = np.argmax(self.gr1[i,j,:])
				nmax2 = np.argmax(self.gr2[i,j,:])
				
				self.gr3[i,j] = self.gr1[i,j,nmax1]
				if not(np.isnan(self.gr3[i,j])):
					if self.gr3[i,j] == 0.:
						self.n1[i,j] = 100.
					else:	
						self.n1[i,j] = int(moden[nmax1])
				else:
					self.n1[i,j] = np.nan
				self.gr4[i,j] = self.gr2[i,j,nmax2]
				if not(np.isnan(self.gr4[i,j])):
					if self.gr4[i,j] == 0.:
						self.n2[i,j] = 100.
					else:
						self.n2[i,j] = int(moden[nmax2])
				else:
					self.n2[i,j] = np.nan

			os.chdir(currdir)
				
		return
		
	def draw_plot(self,isblock=False):

		self.read_eped_prof(self.eped_file)
		grdia_ratio = np.sqrt(self.ncore/self.nsep)
			
		plt.figure('PBM scan alf')		
		cs=plt.contourf(self.hh,self.ww,np.transpose(self.gr3))
		plt.colorbar()

		for i in range(self.nw):
			for j in range(self.nh):
				if (self.gr4[i,j] > 0.0):
					plt.text(self.hh[j,i],self.ww[j,i],str(int(self.n1[i,j])),color='r')
		
		plt.scatter(1,1,s=150,c='orange',marker='*')
		plt.contour(cs,levels=[self.gr_crit1],colors='pink')
		plt.title('Pedestal scan ($\gamma$/$w_A$)')
		if not os.getenv("PEDSCAN_TYPE") == 'TWOWID': 
			plt.xlabel('Pedestal height weight')
			plt.ylabel('Pedestal width weight')
		else:
			plt.xlabel('n Pedestal width weight')
			plt.ylabel('T Pedestal width weight')

		plt.xlim(self.hmin-0.05,self.hmax+0.05)
		plt.ylim(self.wmin-0.05,self.wmax+0.05)

		plt.figure('PBM scan dia')	
		cs = plt.contourf(self.hh,self.ww,np.transpose(self.gr4))
		plt.colorbar()
		
		for i in range(self.nw):
			for j in range(self.nh):
				if (self.gr4[i,j] > 0.0):
					plt.text(self.hh[j,i],self.ww[j,i],str(int(self.n2[i,j])),color='r')
					
		plt.scatter(1,1,s=150,c='orange',marker='*')
		plt.scatter(self.hh[self.target_i2,self.target_j2],self.ww[self.target_i2,self.target_j2],s=150,c='gray',marker='x')
		plt.contour(cs,levels=[self.gr_crit2],colors='pink')
		plt.contour(cs,levels=[self.gr_crit2*grdia_ratio],colors='orange')
		plt.title('Pedestal scan ($\gamma$/$w_*$)')
		if not os.getenv("PEDSCAN_TYPE") == 'TWOWID':
			plt.xlabel('Pedestal height weight')
			plt.ylabel('Pedestal width weight')
		else:
			plt.xlabel('n Pedestal width weight')
			plt.ylabel('T Pedestal width weight')
		plt.xlim(self.hmin-0.05,self.hmax+0.05)
		plt.ylim(self.wmin-0.05,self.wmax+0.05)

		plt.show(block=False)
		if isblock:	input()
		
		return

	def draw_plot2(self,fig_ex=None,use_dia=True):

		self.read_eped_prof(self.eped_file)
		grdia_ratio = np.sqrt(self.ncore/self.nsep)

		if fig_ex == None:	return
		else:
			fig = fig_ex
			len2 = len(fig.axes)
			if len2 == 1:	[ax1] = fig.axes
			else:
				[ax1, ax2] = fig.axes
				ax2.cla()
			ax1.cla()

		if not os.getenv("PEDSCAN_TYPE") == 'TWOWID':
			ax1.set_xlabel('Pedestal height [$h/h_{ref}$]')
			ax1.set_ylabel('Pedestal width [$w/w_{ref}$]')
		else:
			ax1.set_xlabel('n Pedestal width [$h/h_{ref}$]')
			ax1.set_ylabel('T Pedestal width [$h/h_{ref}$]')

		plegend = [];	slegend = [];

		if use_dia:
			ax1.set_title('$\gamma$ / $\omega_{*i}$')
			cs = ax1.contourf(self.hh,self.ww,np.transpose(self.gr4))
			cs1 = ax1.contour(cs,levels=[self.gr_crit2],colors='magenta')
			ax1.contour(cs,levels=[self.gr_crit2*grdia_ratio],colors='orange')
			cs2 = ax1.contour(self.hh,self.ww,np.transpose(self.n2),levels=[self.targetn],colors='silver')
			cs.set_clim(0,1)
			h1,_ = cs1.legend_elements()
			h2,_ = cs2.legend_elements()
			plegend.append(h1[0])
			plegend.append(h2[0])
			slegend.append('Marginal stability')
			slegend.append('Target n = %i'%self.targetn)
		else:
			ax1.set_title('$\gamma$ / $\omega_{A}$')
			cs = ax1.contourf(self.hh,self.ww,np.transpose(self.gr3))
			cs1 = ax1.contour(cs,levels=[self.gr_crit1],colors='magenta')
			cs2 = ax1.contour(self.hh,self.ww,np.transpose(self.n1),levels=[self.targetn],colors='silver')
			cs.set_clim(0,0.1)
			h1, = cs1.legend_elements()
			h2, = cs2.legend_elements()
			plegend.append(h1[0])
			plegend.append(h2[0])
			slegend.append('Marginal stability')
			slegend.append('Target n = %i'%self.targetn)			
		if len2 == 1:
			fig.colorbar(cs,orientation='vertical')
		else:
			fig.colorbar(cs,cax=ax2,orientation='vertical')
	
		for i in range(self.nw):
			for j in range(self.nh):
				if use_dia:				
					if self.gr4[i,j]>0.0:
						ax1.text(self.hh[j,i],self.ww[j,i],str(int(self.n2[i,j])),color='r')
				else:
					if self.gr3[i,j]>0.0:
						ax1.text(self.hh[j,i],self.ww[j,i],str(int(self.n1[i,j])),color='r')

		line = ax1.scatter(self.target_i,self.target_j,s=150,c='magenta',marker='*')
		plegend.append(line);	slegend.append('Selected point (profile)');
		#ax1.errorbar(self.target_i,self.target_j, yerr = [-0.05, 0.05], xerr =[-0.05, 0.05], fmt='o',markersize='10',c='orange',ecolor='orange',capthick=2)
		line = ax1.scatter(self.hh[self.target_i2,self.target_j2],self.ww[self.target_i2,self.target_j2],s=150,c='gray',marker='x')
		plegend.append(line);	slegend.append('Selected point (spectrum)');

		fig.tight_layout()
		ax1.legend(plegend,slegend)
		ax1.axis('equal')
		ax1.set_xlim(self.hmin-0.1,self.hmax+0.1)
		ax1.set_ylim(self.wmin-0.1,self.wmax+0.1)		

		return

	def draw_plot3(self,fig_ex,ind_i=0,ind_j=0,use_dia=True):

		ind_i = int(ind_i)
		ind_j = int(ind_j)

		if fig_ex == None:	return
		else:
			fig = fig_ex
			len2 = len(fig.axes)
			if len2 == 1:	[ax1] = fig.axes
			else:
				[ax1, ax2] = fig.axes
				ax2.cla()
			ax1.cla()

		if not use_dia:	ax1.plot(self.moden,self.gr1[ind_j,ind_i,:],'--',c='gold') 
		else:	ax1.plot(self.moden,self.gr2[ind_j,ind_i,:],'--',c='gold') 
		for k in range(self.nn):
			if not use_dia:
				if self.gr1[ind_j,ind_i,k]<self.gr_crit1:
					ax1.scatter(self.moden[k],self.gr1[ind_j,ind_i,k],c='b')
				else:
					ax1.scatter(self.moden[k],self.gr1[ind_j,ind_i,k],c='r')
			else:   
				if self.gr2[ind_j,ind_i,k]<self.gr_crit2:
					ax1.scatter(self.moden[k],self.gr2[ind_j,ind_i,k],c='b')
				else:
					ax1.scatter(self.moden[k],self.gr2[ind_j,ind_i,k],c='r')


		titlen='growth rate  [$h/h_{ref}$] = %3.1f,  [$w/w_{ref}$] = %3.1f ' %(self.hh[ind_i,ind_j],self.ww[ind_i,ind_j])
		ax1.set_title('%s' %titlen)
		ax1.set_xlabel('mode n')

		if not use_dia:	ax1.set_ylabel('$\gamma$ / $\omega_{A}$')
		else:	ax1.set_ylabel('$\gamma$ / $\omega_{i*}$')
		ax1.set_xlim([min(self.moden)-1,max(self.moden)+1])
		minx = min(self.moden)-1
		maxx = max(self.moden)+1            
		delx = maxx-minx
		delx = int(delx/6)

		return

	def draw_plot4(self,fig_ex):

		self.read_eped_prof(self.eped_file)
		if (self.ext_ti):
			self.read_kin_prof(self.ti_file)		

		if fig_ex == None:	return
		else:
			fig = fig_ex
			len2 = len(fig.axes)
			if len2 == 1:	[ax1] = fig.axes
			elif len2 == 2:	
				[ax1, ax2] = fig.axes
				ax2.cla()
			else:
				[ax1, ax2, ax3] = fig.axes
				ax2.cla()
				ax3.cla()
			ax1.cla()

		legend1 = [];	legend2 = [];	legend3 = [];
		plegend1 = [];	plegend2 = [];	plegend3 = [];
	
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
			line = ax1.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='lime',ecolor='lime',capthick=2)
			legend1.append('Raw_data')
			plegend1.append(line)
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
			line = ax2.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='lime',ecolor='lime',capthick=2)
			legend2.append('Raw_data')
			plegend2.append(line)
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
			line = ax3.errorbar(dat[:,0],dat[:,1], yerr = dat[:,2],fmt='x',markersize='5',c='lime',ecolor='lime',capthick=2)
			legend3.append('Raw_data')
			plegend3.append(line)
		except:
			pass

		ind_j = (1.-self.wmin)/(self.wmax-self.wmin)*(self.nw-1)
		ind_i = (1.-self.hmin)/(self.hmax-self.hmin)*(self.nh-1)

		self.make_prof(ind_j,ind_i,False,True)

		p1, = ax1.plot(self.psin,self.ne)
		p2, = ax2.plot(self.psin,self.te)
		p3, = ax3.plot(self.psin,self.ti)
		legend1.append('Reference');	plegend1.append(p1)
		legend2.append('Reference');	plegend2.append(p2)
		legend3.append('Reference');	plegend3.append(p3)

		ind_j = (self.target_j-self.wmin)/(self.wmax-self.wmin)*(self.nw-1)
		ind_i = (self.target_i-self.hmin)/(self.hmax-self.hmin)*(self.nh-1)

		self.make_prof(ind_j,ind_i,False,True)
		p1, = ax1.plot(self.psin,self.ne)
		p2, = ax2.plot(self.psin,self.te)
		p3, = ax3.plot(self.psin,self.ti)
		legend1.append('Target');	plegend1.append(p1)
		legend2.append('Target');	plegend2.append(p2)
		legend3.append('Target');	plegend3.append(p3)

		max1 = max(self.ne)*1.2
		max2 = max(self.te)*1.2
		max3 = max(self.ti)*1.2

		ax1.set_title("NE [10(19)/m3]")
		ax1.set_xlabel('Normalized radius ($\psi_N$)')
		ax1.set_ylabel('NE [10(19)/m3]')
		ax1.legend(plegend1,legend1)
		ax1.set_xlim((-0.05,1.05))			
		ax1.set_ylim((-0.1,max1))

		ax2.set_title("TE [keV]")
		ax2.set_xlabel('Normalized radius ($\psi_N$)')
		ax2.set_ylabel('TE [keV]')
		ax2.legend(plegend2,legend2)
		ax2.set_xlim((-0.05,1.05))
		ax2.set_ylim((-0.1,max2))

		ax3.set_title("TI [keV]")
		ax3.set_xlabel('Normalized radius ($\psi_N$)')
		ax3.set_ylabel('TI [keV]')
		ax3.legend(plegend3,legend3)
		ax3.set_xlim((-0.05,1.05))
		ax3.set_ylim((-0.1,max3))
		fig.tight_layout()

		return

	def make_target_prof(self,h_weight,w_weight):
	
		nw = (w_weight - self.wmin)/(self.wmax - self.wmin) * (self.nw-1)
		nh = (h_weight - self.hmin)/(self.hmax - self.hmin) * (self.nh-1)
	
		self.read_eped_prof(self.eped_file)
		if (self.ext_ti):
			self.read_kin_prof(self.ti_file)
	
		print('generate pedestal with height weight = %5.3f and width weight = %5.3f'%(h_weight,w_weight))	
		self.make_prof(nw,nh,True)
		
		return		
		
	
	def __init__(self,filename='scan_opt',clearf = False,nocheaserun=False,ingfitp=False):
	
		self.ingfitp = ingfitp
		self.initialise_vars()
		self.read_namelist(filename)

		if not ingfitp:
			try: os.mkdir(self.work_dir+'/input')
			except:	pass
			indir = self.work_dir+'/input/'
			try:	copyfile(self.eqdsk_file,indir+self.eqdsk_file.split('/')[-1])
			except:	pass
			try:	copyfile(self.ti_file,indir+self.ti_file.split('/')[-1])
			except:	pass
			try:	copyfile(self.eped_file,indir+self.eped_file.split('/')[-1])
			except:	pass

			self.eqdsk_file = self.work_dir+'/input/'+self.eqdsk_file.split('/')[-1]
			self.ti_file = self.work_dir+'/input/'+self.ti_file.split('/')[-1]
			self.eped_file = self.work_dir+'/input/'+self.eped_file.split('/')[-1]
			#print(os.getcwd(),'a1')
		else:
			self.eqdsk_file = 'input/'+self.eqdsk_file.split('/')[-1]
			self.ti_file = 'input/'+self.ti_file.split('/')[-1]
			self.eped_file = 'input/'+self.eped_file.split('/')[-1]
			#print(os.getcwd(),'a2')
			self.work_dir = os.getcwd()

		self.ch.eqdsk_name = self.eqdsk_file
		if clearf:	return
		currdir = os.getcwd()
		if not nocheaserun:	self.ch.load_eqdsk()
		if(self.beta_critval == 0.0):	
			eq = eqdsk.eqdsk(self.eqdsk_file)
			eq.read_eqdsk_file()
			if self.beta_criterion == 2:
				self.beta_critval = eq.rmag
				print('EQ axis',eq.rmag,'[m]')
			elif self.beta_criterion == 1:
				self.beta_critval = eq.bp
				print('Bp',eq.bp)
		os.chdir(currdir)
		return
	
	if __name__ == "__main__":

		import pedscan
		inputf = 'scan_opt'
		try:	inputf = sys.argv[1]
		except:	pass

		ped = pedscan.pedscan(inputf)
			
		ped.check_inp()
		ped.read_eped_prof(ped.eped_file)
		if (ped.ext_ti):
			ped.read_kin_prof(ped.ti_file)
		
		ped.make_directories()
		ped.submit_script()

		ped.collect_data()
		ped.clear_data()

