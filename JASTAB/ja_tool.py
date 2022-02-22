#tool.py
##!!!!!!!!!!!!!!!!!!!!!!! alf_over_omega ne[0] -> ne[-1]
import os,stat
import sys
import scipy
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import subprocess
from shutil import copyfile
from progress import update_progress
from exec_dirs import python3_exec,helena_exec,chease_exec,mis_exec,elite_dir,node_default
from batch_run import *

currdir = os.getcwd()
e0 = 1.601* 1.e-19;
mu0 = 4.0 * np.pi *1.e-7;

class simulation:

	def str_to_logic(self,s):
		s = s.lower()
		if (s == 'true'):
			return True
		if (s == 'false'):
			return False
				
	def read_namelist(self,namelist):
	
		try:
			f = open(namelist,'r')
		except:
			print('>>> no namelist file')
			exit()

		yesline = True
		line = f.readline()
		while(yesline):
			line = f.readline()
			if not line: break
			try:
				var_name = line.split('=')[0].split()[0].lower()
				vars = line.split('=')[1]
			except:
				var_name = 'NONE'
			
			if (var_name == 'equ_name'):
				self.init_HELinp = vars.split()[0]
			elif(var_name == 'chease_name'):
				self.init_CHEinp = vars.split()[0]
			elif (var_name == 'equ_type'):
				self.equ_type = int(vars.split()[0])
			elif (var_name == 'run_equ'):
				self.run_equ = vars.split()[0]	
			elif (var_name == 'run_id'):
				self.RUN_ID = vars.split()[0]
			elif (var_name == 'grid_nw'):
				self.grid_n1 = int(vars.split()[0])
			elif (var_name == 'grid_nh'):
				self.grid_n2 = int(vars.split()[0])		
			elif (var_name == 'run_dir'):
				self.run_dir = vars.split()[0]
			elif (var_name == 'pedestal_center'):
				self.pedc = float(vars.split()[0])
			elif (var_name == 'use_ped_finder'):
				self.use_ped_find = vars.split()[0]
			elif (var_name == 'use_ped_plot'):
				self.use_ped_plot = vars.split()[0]
			elif (var_name == 'npt'):
				self.NPT = int(vars.split()[0])
			elif (var_name == 'pmin'):
				self.plow = float(vars.split()[0])
			elif (var_name == 'pmax'):
				self.phigh = float(vars.split()[0])
			elif (var_name == 'jmin'):
				self.jlow = float(vars.split()[0])
			elif (var_name == 'jmax'):
				self.jhigh = float(vars.split()[0])
			elif (var_name == 'b'):
				self.b = float(vars.split()[0])
			elif (var_name == 'coreneo'):
				self.coreneo = float(vars.split()[0])
			elif (var_name == 'use_prev'):
				self.use_prev = vars.split()[0]
			elif (var_name == 'use_hag'):
				self.use_hag = vars.split()[0]
			elif (var_name == 'use_neo'):
				self.use_neo = vars.split()[0]
			elif (var_name == 'use_chang'):
				self.use_chang = vars.split()[0]				
			elif (var_name == 'qsub_dir'):
				self.qsub_dir = vars.split()[0]
			elif (var_name == 'qstat_dir'):
				self.qstat_dir = vars.split()[0]
			elif (var_name == 'nodelist'):
				self.nodelist = vars.split()[0]
			elif (var_name == 'python_dir'):
				self.python_dir = vars.split()[0]
			elif (var_name == 'batch_run'):
				self.batch_run = vars.split()[0]
			elif (var_name == 'helena_dir'):
				self.helena_dir = vars.split()[0]
			elif (var_name == 'chease_dir'):
				self.chease_dir = vars.split()[0]
			elif (var_name == 'mis_dir'):
				self.mis_dir = vars.split()[0]
			elif (var_name == 'elite_dir'):
				self.elite_dir = vars.split()[0]
			elif (var_name == 'qdelfix'):
				self.qdelfix = float(vars.split()[0])
			elif (var_name == 'qdel_crit'):
				self.qdel_crit = float(vars.split()[0])
			elif (var_name == 'highq'):
				self.highq = vars.split()[0]
			elif (var_name == 'run_stab'):
				self.run_stab = vars.split()[0]
			elif (var_name == 'ias'):
				self.ias = vars.split()[0]
			elif (var_name == 'psis'):
				self.psis = float(vars.split()[0])
			elif (var_name == 'xr2'):
				self.xr2 = float(vars.split()[0])
			elif (var_name == 'sig2'):
				self.sig2 = float(vars.split()[0])
			elif (var_name == 'xr1'):
				self.xr1 = float(vars.split()[0])
			elif (var_name == 'sig1'):
				self.sig1 = float(vars.split()[0])
			elif (var_name == 'ngride'):
				self.ngride = int(vars.split()[0])
			elif (var_name == 'psise'):
				self.psise = float(vars.split()[0])
			elif (var_name == 'ndist'):
				self.ndist = int(vars.split()[0])
			elif (var_name == 'moden'):
				vars = vars.split('\n')[0].split('\t')[0]
				self.modens = vars
			elif (var_name == 'modene'):
				vars = vars.split('\n')[0].split('\t')[0]
				self.modenes = vars
			elif (var_name == 'run_stab'):
				self.run_stab = vars.split()[0]	
			elif (var_name == 'gridn'):
				self.gridn = int(vars.split()[0])
			elif (var_name == 'cns'):
				self.CNS = int(vars.split()[0])
			elif (var_name == 'cnt'):
				self.CNT = int(vars.split()[0])
			elif (var_name == 'cnpsi'):
				self.CNPSI = int(vars.split()[0])
			elif (var_name == 'cnchi'):
				self.CNCHI = int(vars.split()[0])		
			elif (var_name == 'cnse'):
				self.CNSE = int(vars.split()[0])
			elif (var_name == 'cnte'):
				self.CNTE = int(vars.split()[0])
			elif (var_name == 'cnpsie'):
				self.CNPSIE = int(vars.split()[0])
			elif (var_name == 'cnchie'):
				self.CNCHIE = int(vars.split()[0])				
			elif (var_name == 'use_kin_prof'):
				self.use_kin_prof = vars.split()[0]
			elif (var_name == 'kin_prof_name'):
				self.kin_prof_name = vars.split()[0]
			elif (var_name == 'rot_prof_name'):
				self.rot_prof_name = vars.split()[0]
			elif (var_name == 'use_prescribed_helena'):
				self.use_prescribed_helena = vars.split()[0]
			elif (var_name == 'use_rotation'):
				self.use_rot = vars.split()[0]
			elif (var_name == 'use_compression'):
				self.use_comp = vars.split()[0]
			elif (var_name == 'del_data'):
				self.del_data = vars.split()[0]
			elif (var_name == 'beta_type'):
				self.beta_criterion = float(vars.split()[0])
			elif (var_name == 'adjust_prof'):
				self.adjust_prof = vars.split()[0]
			elif (var_name == 'line_den'):
				self.line_den = float(vars.split()[0])
			elif (var_name == 'target_bnd'):
				self.target_bnd = float(vars.split()[0])
			elif (var_name == 'epsilon'):
				self.epsilon = float(vars.split()[0])
			elif (var_name == 'beta_crit'):
				self.beta_crit = float(vars.split()[0])
			elif (var_name == 'li_crit'):
				self.li_crit = float(vars.split()[0])
			elif (var_name == 'ip_crit'):
				self.ip_crit = float(vars.split()[0])
			elif (var_name == 'bs_crit'):
				self.bs_crit = float(vars.split()[0])
			elif (var_name == 'current_itern'):
				self.iternc = float(vars.split()[0])			
			elif (var_name == 'li_intern'):
				self.iternl = float(vars.split()[0])
			elif (var_name == 'relax'):
				self.relax = float(vars.split()[0])
			elif (var_name == 'li_target'):
				self.li_target = float(vars.split()[0])
			elif (var_name == 'beta_target'):
				self.beta_target = float(vars.split()[0])				
			elif (var_name == 'use_li'):
				self.use_li = vars.split()[0]
			elif (var_name == 'use_li2'):
				self.use_li2 = vars.split()[0]				
			elif (var_name == 'use_core_mod'):
				self.use_core_mod = vars.split()[0]
			elif (var_name == 'use_core_mod_psin'):
				self.use_core_mod_psin = float(vars.split()[0])
			elif (var_name == 'bsmulti'):
				self.bsmulti = float(vars.split()[0])
			elif (var_name == 'apf'):
				self.apf = float(vars.split()[0])
			elif (var_name == 'bpf'):
				self.bpf = float(vars.split()[0])								
			elif (var_name == 'cpf'):
				self.cpf = float(vars.split()[0])
			elif (var_name == 'dpf'):
				self.dpf = float(vars.split()[0])	
			elif (var_name == 'ajf'):
				self.ajf = float(vars.split()[0])
			elif (var_name == 'bjf'):
				self.bjf = float(vars.split()[0])	
			elif (var_name == 'cjf'):
				self.cjf = float(vars.split()[0])
			elif (var_name == 'djf'):
				self.djf = float(vars.split()[0])
			elif (var_name == 'nqa'):
				self.nqa = float(vars.split()[0])
			elif (var_name == 'ncutoff'):
				self.ncut = float(vars.split()[0])
			elif (var_name == 'grcrit1'):
				self.grcrit1 = float(vars.split()[0])
			elif (var_name == 'grcrit2'):
				self.grcrit2 = float(vars.split()[0])
			elif (var_name == 'use_adja'):
				self.use_adja = vars.split()[0]
			elif (var_name == 'use_bilin'):
				self.use_bilin = vars.split()[0]
			elif (var_name == 'rmag'):	
				self.rmag = float(vars.split()[0])
			elif (var_name == 'epsilon'):
				self.epslon = float(vars.split()[0])

		self.run_stab = self.run_stab.lower()
		self.run_equ = self.run_equ.lower()
		
		self.use_prev = self.str_to_logic(self.use_prev)
		self.use_hag = self.str_to_logic(self.use_hag)
		self.use_neo = self.str_to_logic(self.use_neo)
		self.use_chang = self.str_to_logic(self.use_chang)
		self.batch_run = self.str_to_logic(self.batch_run)
		self.highq = self.str_to_logic(self.highq)
		self.ias = self.str_to_logic(self.ias)
		self.use_kin_prof = self.str_to_logic(self.use_kin_prof)
		self.use_prescribed_helena = self.str_to_logic(self.use_prescribed_helena)
		self.del_data = self.str_to_logic(self.del_data)
		self.adjust_prof = self.str_to_logic(self.adjust_prof)
		self.use_core_mod = self.str_to_logic(self.use_core_mod)

		self.use_ped_find = self.str_to_logic(self.use_ped_find)
		self.use_ped_plot = self.str_to_logic(self.use_ped_plot)
		self.use_rot = self.str_to_logic(self.use_rot)
		self.use_comp = self.str_to_logic(self.use_comp)
		self.use_li = self.str_to_logic(self.use_li)
		self.use_li2 = self.str_to_logic(self.use_li2)
		self.use_adja = self.str_to_logic(self.use_adja)
		self.use_bilin = self.str_to_logic(self.use_bilin)
	
		self.init_HELinp = currdir +'/input/' + self.init_HELinp

		if (self.run_equ == 'chease'):
			self.init_CHEinp = currdir+'/input/' + self.init_CHEinp
		if (self.use_kin_prof):
			self.kin_prof_name = currdir+ '/input/' + self.kin_prof_name
		
		if not(self.use_kin_prof):
			self.use_prescribed_helena = True

		if (self.use_rot):
			self.use_comp = True
			self.rot_prof_name = currdir+'/input/' + self.rot_prof_name 
		
		self.run_dir = currdir + '/ja_' + self.RUN_ID
		f.close()
		return
	
	def make_directory(self,directory,dir_name):
		try:
			os.mkdir(directory)
		except:
			print('>>> ' + dir_name+' directory exists!')
		return()
		
	def read_nlines_and_split(self,file,n):
		b2=[]
		for i in range(n):
			ifl=file.readline()
			a=ifl.split()
			for arg in a:
				b2.append(float(arg))
		return np.array(b2)

	def read_2d_array(self,file,n1,n2):
		b3=[]
		for i in range(n2):
			a=self.read_nlines_and_split(file,n1)
			b3.append(a)
		return np.array(b3)
		
	def read_elite(self,filename):
	
		file = open(filename,"r")
			
		line = file.readline()
		line = file.readline() # number of radial and poloidal points
		sp=line.split()

		npsi=int(sp[0])
		nchi=int(sp[1])

		mod1=npsi%5

		nlines1=int(npsi/5)+mod1
		nlines2=nchi
		line = file.readline() # psi mesh
		psi=self.read_nlines_and_split(file,nlines1)
		line = file.readline() # dp/dpsi
		dpdpsi=self.read_nlines_and_split(file,nlines1)
		line = file.readline() # d2p/dpsi2
		dp2dpsi=self.read_nlines_and_split(file,nlines1)
		line = file.readline() # fpol
		fpol=self.read_nlines_and_split(file,nlines1)
		line = file.readline() # fdf/dpsi
		ffp=self.read_nlines_and_split(file,nlines1)
		line = file.readline() # dfdf/dpsi2
		dffp=self.read_nlines_and_split(file,nlines1)
		line = file.readline() # q
		q =self.read_nlines_and_split(file,nlines1)
		line = file.readline() # R
		R=self.read_2d_array(file,nlines1,nlines2)
		line = file.readline() # Z
		Z=self.read_2d_array(file,nlines1,nlines2)
		line = file.readline() # Bp
		Bp=self.read_2d_array(file,nlines1,nlines2)
		line = file.readline() # ne
		ne =self.read_nlines_and_split(file,nlines1)
		line = file.readline() # dne/dpsi
		dne =self.read_nlines_and_split(file,nlines1)
		line = file.readline() # te
		te =self.read_nlines_and_split(file,nlines1)
		line = file.readline() # dte/dpsi
		dte =self.read_nlines_and_split(file,nlines1)
		line = file.readline() # ti
		ti =self.read_nlines_and_split(file,nlines1)
		line = file.readline() # dti/dpsi
		dti =self.read_nlines_and_split(file,nlines1)
		line = file.readline() # n main ion
		nmainion =self.read_nlines_and_split(file,nlines1)
		line = file.readline() # nZ
		nZ =self.read_nlines_and_split(file,nlines1)
		line = file.readline() # Zeff
		zeff = float(file.readline())
		line = file.readline() # Zimp
		zimp=float(file.readline())
		line = file.readline() # Aimp
		Aimp=float(file.readline())
		
		massnumber=2.0

		# calculates diamagnetic and Alfven frequencies

		rav = (min(R[:,npsi-1])+R[1,npsi-1])/2

		B=fpol[0]/rav
		mu0=4e-7*np.pi
		mp= 1.6726231e-27 # proton mass
		md=mp*massnumber
		rho=ne[-1]*md   ##############################ne[-1] - > ne[0]
#		rho=nmainion[0]*md
		omegan=2.9979e10/ne/4.8032e-10*dne*1.6022e-12*ti
		omegaspi=(omegan+2.9979e10/ne/4.8032e-10*dti*1.6022e-12*ne)/1e8
		omega_alf=np.sqrt(B*B/rho/rav/rav/mu0)#/q[0]/q[0])
		alf_over_maxomega=-omega_alf/np.min(omegaspi)#[0:200])

		file.close()
		return(npsi,psi,dpdpsi,ffp,te,ne,ti,zeff,zimp,alf_over_maxomega)
				
	def get_equ_prof(self,filename,linen):
		
		run_command1 = 'cat '+filename+' | grep -ri ALPHA -A ' + str(linen +1);
		run_command2 = 'cat '+filename+' | grep -ri AVERAGE -A ' + str(linen);

		ALPt = np.zeros(shape=(linen-1,1))
		ALPt2 = np.zeros(shape=(linen-1,1))
		JPHt = np.zeros(shape=(linen-2,1))
		PSIN = np.zeros(shape=(linen,1))

		stat, out = subprocess.getstatusoutput(run_command1)
		for i in range(linen-1):
			linet = out.split('\n')[i+3]
			ALPt[i] = float(linet.split()[7])
			ALPt2[i] = float(linet.split()[6])
			PSIN[i+1] = float(linet.split()[1])

		PSIN[0] = 0.0;
		
		count = 0
		target_psi = 0.8
		
		for i in range(linen-1):
			if ((PSIN[i]-target_psi)*(PSIN[i+1]-target_psi) < 0):
				count = i;
		
		stat, out = subprocess.getstatusoutput(run_command2)
		for i in range(linen-2):
			linet = out.split('\n')[i+2]
			JPHt[i] = float(linet.split()[1])

		linet = out.split('\n')[linen]
		bsfrac = float(linet.split()[1])
		alpm = max(ALPt)
		alpm2 = max(ALPt2)
		jphim = max(JPHt[count:linen])
		return (alpm, alpm2, jphim, PSIN, bsfrac)
		
	def get_corep(self,filename):
		
		run_command1 = 'cat '+filename+' | grep -ri REAL -A15'

		stat, out = subprocess.getstatusoutput(run_command1)

		linet = out.split('\n')[14]
		linet2 = out.split('\n')[15]
		
		psin1 = float(linet.split()[0])
		pmax1 = float(linet.split()[1])

		psin2 = float(linet2.split()[0])
		pmax2 = float(linet2.split()[1])	
		
		corep = (pmax2-pmax1)/(psin2-psin1)*(0.-psin1) + pmax1
			
		return (corep)	
		
	def get_prof_opt(self,prof_data):

		for i in range(np.size(prof_data)):
			line1 = prof_data[i]

			if(line1.split()[0] == 'NPTS'):
				linet = line1
				npts = int(linet.split('=')[1].split(',')[0])
				
		return(npts)

	def get_num_opt(self,num_data):
		
		for i in range(np.size(num_data)):
			line1 = num_data[i]
			
			linen = line1.split('=')[0].split()[0]
			if(linen == 'NRMAP' or linen == 'nrmap'):
				mapn = int(line1.split('=')[1].split(',')[0])
				
			if(linen == 'NPMAP' or linen == 'npmap'):
				mapnp = int(line1.split('=')[1].split(',')[0])

		return(mapn,mapnp)
	
	def get_HELENA(self,filename):

		f1 = open(filename,'r')
		count = 0
		shape_data = []
		prof_data = []
		phys_data = []
		num_data = []	
		while True:
			line1 = f1.readline()
			if not line1: break
			print(line1)		
			if(line1.split()[0] == '&SHAPE'): 
				count = 1
			if(line1.split()[0] == '&PROFILE' or line1.split()[0] == '&profile'): 
				count = 2
			if(line1.split()[0] == '&PHYS'): 
				count = 3
			if(line1.split()[0] == '&NUM'): 
				count = 4
			if(line1.split()[0] == '&END'): 
				count = 0
			
			if(count == 1):
				shape_data.append(line1)
			if(count == 2):
				prof_data.append(line1)
			if(count == 3):
				phys_data.append(line1)
			if(count == 4):
				num_data.append(line1)
			if(line1.split()[0] == 'XIAB'):
				xiab = float(line1.split('=')[1])
			if(line1.split()[0] == 'BVAC'):
				bvac = float(line1.split('=')[1])
			if(line1.split()[0] == 'RVAC'):
				rvac = float(line1.split('=')[1])
			if(line1.split()[0] == 'EPS'):
				eps = float(line1.split('=')[1])
		
		f1.close()			
		return (shape_data, prof_data, phys_data, num_data, xiab, bvac, rvac, eps)

	def read_zjz_data1(self,filename,npts):

		f3 = open(filename,'r')
		zjz_read = []
		for i in range(npts):
			
			line1 = f3.readline()
			linet = line1.split('=')[1].split(',')[0]

			zjz_read.append(float(linet))
	
		f3.close()
		return (zjz_read)

	def read_zjz_data2(self,filename,npts):

		f3 = open(filename,'r')
		zjz_read = np.zeros(shape=(npts-2,2))
		while True:
			line = f3.readline()
			if (line.find('AVERAGE') > -1):
				line = f3.readline()
				for i in range(npts-2):	
					line = f3.readline()
					zjz_read[i,0] = float(line.split()[0]) ** 2
					zjz_read[i,1] = float(line.split()[1])
				break
		f3.close()
		return zjz_read

	def write_prof_data(self,zjz,dpr,df2,ne,te,ti,npts):
		
		prof_data_temp = []

		linet = '&PROFILE\n'
		prof_data_temp.append(linet)
		
		linet = "NPTS = %i, \n"%(npts)
		prof_data_temp.append(linet)
		
		linet = "  IGAM=   7, AGA =-1.000, \n"
		prof_data_temp.append(linet)
		linet = "  IPAI=  7, EPI = 0.95, FPI = 1.0, \n"
		prof_data_temp.append(linet)
		linet = "  ICUR=  2, ECUR=0.95,  FCUR = 1.0,\n"
		prof_data_temp.append(linet)
		
		for i in range(npts):
			linet = "zjz(%i) = %f, dpr(%i) = %f, df2(%i) = %f \n"%(i+1,zjz[i]/zjz[0],i+1,dpr[i],i+1,df2[i]/df2[0])
			prof_data_temp.append(linet)

		for i in range(npts):
			if (self.use_prescribed_helena):
				linet = "neprof(%i) = %6.4f, teprof(%i) = %6.4f, tiprof(%i) = %6.4f, \n"%(i+1,ne[i]/1.e19,i+1,te[i]/1.e3,i+1,ti[i]/1.e3)
			elif(self.use_kin_prof):
				psintt = float(i)/float(npts-1)
				nef = interp1d(self.PSIX,self.NE,'cubic')
				tef = interp1d(self.PSIX,self.TE,'cubic')
				tif = interp1d(self.PSIX,self.TI,'cubic')
				nif = interp1d(self.PSIX,self.NI,'cubic')
				
				linet = "neprof(%i) = %6.4f, teprof(%i) = %6.4f, niprof(%i) = %6.4f, tiprof(%i) = %6.4f, \n"%(i+1,nef(psintt),i+1,tef(psintt),i+1,nif(psintt),i+1,tif(psintt))
				
			prof_data_temp.append(linet)
		return(prof_data_temp)

	def write_phys_data(self,xiab,bvac,rvac,eps,zeff,zimp,corep,coreneo,b):
		
		phys_data_temp = []
		
		linet = '&PHYS\n'
		phys_data_temp.append(linet)
		linet = 'IDETE = 10 \n'
		phys_data_temp.append(linet)
		linet = 'NEOSIG = .True. \n'
		phys_data_temp.append(linet)
		linet = 'B = %f \n'%(b)
		phys_data_temp.append(linet)
		linet = 'XIAB = %f \n'%xiab
		phys_data_temp.append(linet)
		linet = 'BVAC = %f \n'%(bvac)
		phys_data_temp.append(linet)
		linet = 'RVAC = %f \n'%(rvac)
		phys_data_temp.append(linet)
		linet = 'EPS = %f \n'%(eps)
		phys_data_temp.append(linet)
		linet = 'SCALEMIX = 0.700000 \n'
		phys_data_temp.append(linet)
		linet = 'ZJZSMOOTH = 0.300000 \n'
		phys_data_temp.append(linet)
		linet = 'BSMULTIP = 1.000000 \n'
		phys_data_temp.append(linet)
		if not(self.use_kin_prof):
			linet = 'ZEFF = %f \n'%(zeff)
		elif (self.use_prescribed_helena):
			linet = 'ZEFF = %f \n'%(self.zeff)
		phys_data_temp.append(linet)
		if not(self.use_kin_prof):
			linet = 'ZIMP = %f \n'%(zimp)
		elif (self.use_prescribed_helena):
			linet = 'ZIMP = %f \n'%(self.zimp)
		phys_data_temp.append(linet)
		linet = 'COREP = %f \n'%(corep)
		phys_data_temp.append(linet)
		linet = 'BETAP = -1 \n'
		phys_data_temp.append(linet)
		linet = 'BKEEP = .False. \n'
		phys_data_temp.append(linet)
		linet = 'CHANGFIX = .False. \n'
		phys_data_temp.append(linet)
		if (self.use_neo):
			linet = 'NEOFIX = .True. \n'
			phys_data_temp.append(linet)
		if (self.use_hag):
			linet = 'HAGFIX = .True. \n'
			phys_data_temp.append(linet)
		linet = 'IFASTSCA = .true. \n'
		phys_data_temp.append(linet)
		linet = 'SAFEMODE = .true. \n'
		phys_data_temp.append(linet)
		linet = 'CORENEO = %f \n'%(coreneo)
		phys_data_temp.append(linet)

		return(phys_data_temp)
			
	def RUN_HEL1(self, workdir,shape_data, prof_data, phys_data, num_data, b, run_save):
		
		RUN_DIR = workdir
		
		os.system('rm -r ' + RUN_DIR)
		os.mkdir(RUN_DIR)
		os.chdir(RUN_DIR)
		
		stat1 = 1
		stat2 = 1
		lastb = b	
		count = 1
		while(stat1 > 0 or stat2 > 0 or stat3>0):
		
			if (count > 1):
				wegt = float(int((count)/2)*np.power(-1,count))
				weg = np.power(1.5,wegt)
				lastb = b * weg
				bwarn = '>>> HEL_CRASHED -> B Changed: ' + str(b)+' -> '+str(lastb)
				print(bwarn)	
			if (count > 40):
				print('>>> HELENA CRAHSED FOR ALL B')
				exit()

			f3 = open('fort.10','w')

			for i in range(np.size(shape_data)):
				f3.write(shape_data[i])
			f3.write('&END\n')
			for i in range(np.size(prof_data)):
				f3.write(prof_data[i])
			f3.write('&END\n')
			for i in range(np.size(phys_data)):
				f3.write(phys_data[i])
			f3.write('B = %f \n'%(lastb))
			
			f3.write('&END\n')
			for i in range(np.size(num_data)):
				f3.write(num_data[i])
			f3.write('&END\n')
			f3.write('&PRI\n NPR1 = 1 \n &END \n &PLOT \n NPL1 = 1 \n &END \n &BALL \n &END \n')
			f3.close()
		
			stat,out = subprocess.getstatusoutput(self.helena_dir)
		
			stat1 = out.find('CRASHED')
			stat2 = out.find('severe')
			stat3 = out.find('SPLINE')
			count = count + 1
		os.chdir(workdir)
		
		if (run_save == 1):
		
			os.chdir(RUN_DIR)
			os.system('cp fort.10 ' + workdir + '/fort_res_inp')
			os.system('cp fort.20 ' + workdir + '/fort_res')
			os.chdir(workdir)

		return(RUN_DIR,lastb)

	def RUN_HEL2(self,filename,xiab2,btor2,initB2):
        
		status, output = subprocess.getstatusoutput(self.helena_dir)

		stat1 = 1
		stat2 = 1
		lastb = initB2
		count = 1
		isrun = True
		while(stat1 > 0 or stat2 > 0 or stat3>0):
			if (count > 1):
				wegt = float(int((count)/2)*np.power(-1,count))
				weg = np.power(1.3,wegt)
				lastb = initB2 * weg
				bwarn = '>>> HEL_CRASHED -> B Changed: ' + str(initB2)+' -> '+str(lastb)
				print(bwarn)
				
			if (count > 30):
				print('>>> HELENA CRAHSED FOR ALL B')
				isrun = False
				exit()		

			t1, t2, t3 = self.new_helinp('tempin',xiab2,btor2,lastb)
			stat,out = subprocess.getstatusoutput(self.helena_dir)

			stat1 = out.find('CRASHED')
			stat2 = out.find('severe')
			stat3 = out.find('SPLINE')
			count = count + 1

		return(lastb,isrun)
		
	def Make_HEL_inp(self, workdir,shape_data, prof_data, phys_data, num_data, lastb, filename):

		f3 = open(filename,'w')

		for i in range(np.size(shape_data)):
			f3.write(shape_data[i])
		f3.write('&END\n')
		for i in range(np.size(prof_data)):
			f3.write(prof_data[i])
		f3.write('&END\n')
		for i in range(np.size(phys_data)):
			f3.write(phys_data[i])
		f3.write('B = %f \n'%(lastb))

		f3.write('&END\n')
		for i in range(np.size(num_data)):
			f3.write(num_data[i])
		f3.write('&END\n')
		f3.write('&PRI\n NPR1 = 1 \n &END \n &PLOT \n NPL1 = 1 \n &END \n &BALL \n &END \n')
		f3.close()
		
		os.system('mv ' + filename + ' ' + workdir + '/.')
		
		return()
		
	def get_intep_data(self, npt, npt1, ps, data):

		data1 = interp1d(ps, data, kind='slinear');
		
		out = np.zeros(shape=(npt,1))
		
		for i in range(npt):
			
			psint = float(i)/float(npt-1)*ps[npt1-1]
			
			out[i] = data1(psint)

		return(out)

	def mult_prof(self, npt, prof, weightc, weightw, weight):

		out = np.zeros(shape=(npt,1))
		
		for i in range(npt):
		
			psint = float(i)/float(npt-1)
			weig = 1. + weight * scipy.exp(-1 * np.power(((psint-weightc)/weightw),2));
			
			out[i] = prof[i] * weig

		return(out)	
		
	def mult_prof2(self,psint, weightc, weightw, weight):

		weig = 1. + weight * scipy.exp(-1 * np.power(((psint-weightc)/weightw),2));
			
		return(weig)		
		
	def make_scan_equilibrium_helena(self):

		# --- main
		b = self.b
		
		# --- run_dir
		work_dir = self.run_dir + '/EQU_HELENA'
		equ_save_dir = self.run_dir + '/equ_inp'

		self.make_directory(work_dir,'work')
		self.make_directory(equ_save_dir,'equ save')

		if (self.use_kin_prof):
			if not (self.equ_type == 1):
				copyfile(self.kin_prof_name,work_dir+'/chease_kinprof')
				copyfile(self.kin_prof_name,equ_save_dir+'/chease_kinprof')
			elif (self.adjust_prof):
				copyfile(self.kin_prof_name+'_new',work_dir+'/chease_kinprof')
				copyfile(self.kin_prof_name+'_new',equ_save_dir+'/chease_kinprof')
			else:
				copyfile(self.kin_prof_name,work_dir+'/chease_kinprof')
				copyfile(self.kin_prof_name,equ_save_dir+'/chease_kinprof')
				
	
		if (self.use_rot):
			copyfile(self.rot_prof_name,equ_save_dir+'/chease_vtor')

		if (self.use_kin_prof):
			if not (self.equ_type == 1):
				copyfile(self.kin_prof_name,'chease_kinprof')
			elif (self.adjust_prof):
				copyfile(self.kin_prof_name+'_new','chease_kinprof')
			else:
				copyfile(self.kin_prof_name,'chease_kinprof')

		if (self.use_rot):
			copyfile(self.rot_prof_name,'chease_vtor')

		print('>>> Calculate given equilibrium')

		hel_shape,hel_prof,hel_phys,hel_num,xiab,bvac,rvac,eps = self.get_HELENA(self.init_HELinp)
		npts = self.get_prof_opt(hel_prof)
		mappn,mappnp = self.get_num_opt(hel_num)

		HELout = work_dir + '/fort.20'
		ELIout = work_dir + '/eliteinp'
		ZJZout = work_dir + '/final_zjz'

		try:
			alpm, alpm2, jphim, PSIN, bsfrac = self.get_equ_prof(HELout,mappn)
			yes_prev = True
		except:
			yes_prev = False

		if not(self.use_prev and yes_prev):
			out,b = self.RUN_HEL1(work_dir,hel_shape,hel_prof,hel_phys,hel_num,b,0)
		else:
			out = work_dir
			print('>>> use previous fort.20')

		ZJZint = self.read_zjz_data2(HELout,mappn)
		ZJZin = np.zeros(shape=(mappn,2))
		for i in range(mappn-2):
			ZJZin[i+1,0] = ZJZint[i,0]
			ZJZin[i+1,1] = ZJZint[i,1]
		ZJZin[0,0] = 0.0
		ZJZin[0,1] = (ZJZin[2,1] - ZJZin[1,1]) / (ZJZin[2,0] - ZJZin[1,0]) * (0.0 - ZJZin[1,0]) + ZJZin[1,1]
		ZJZin[mappn-1,0] = 1.0
		ZJZin[mappn-1,1] = (ZJZin[mappn-2,1] - ZJZin[mappn-3,1]) / (ZJZin[mappn-2,0] - ZJZin[mappn-3,0]) * (1.0 - ZJZin[mappn-3,0]) + ZJZin[mappn-3,1]

		corep = self.get_corep(HELout)
		alpm, alpm2, jphim, PSIN, bsfrac = self.get_equ_prof(HELout,mappn)
		npsi,psi,dpdpsi,fdf,te,ne,ti,zeff,zimp,alf_over_maxomega = self.read_elite(ELIout)

		if (self.use_ped_find):
			psi2 = (psi - psi[0])/(psi[-1] - psi[0])
			dpdpf = interp1d(psi2,abs(dpdpsi),'cubic')
			psi_temp = np.linspace(0.75,1.0,1001)
			dpdp_temp = dpdpf(psi_temp)
			ind = np.argmax(dpdp_temp)
			print('>>> pedestal_center %f -> %f'%(self.pedc,psi_temp[ind]))
			self.pedc = psi_temp[ind]

		pedw = 1.01 - self.pedc
		print('>>> pedestal weight option...')
		print('half width',round(pedw,5),'pedestal center',self.pedc)		

		try:
			zjz_prof = self.read_zjz_data(ZJZout,npts)
		except:
			zjzfit = interp1d(ZJZin[:,0],ZJZin[:,1],'cubic')
			zjz_prof = np.zeros(npts)
			for i in range(npts):
				psi_temp = float(i) / float(npts-1)
				zjz_prof[i] = zjzfit(psi_temp)	

		psitt = []
		for i in range(npts):
			psitt.append(float(i)/float(npts-1))
			
		dpdpf = self.get_intep_data(self.NPT,npsi,psi,dpdpsi)
		dpdpf[0] = dpdpf[1]*2 - dpdpf[2]

		fdff = self.get_intep_data(self.NPT,npsi,psi,fdf)
		tef = self.get_intep_data(self.NPT,npsi,psi,te)
		nef = self.get_intep_data(self.NPT,npsi,psi,ne)
		tif = self.get_intep_data(self.NPT,npsi,psi,ti)
		zjzf = self.get_intep_data(self.NPT,npts,psitt,zjz_prof)


		print('>>> Generate equilibrium files for scan')

		for i in range(self.grid_n1):

			pweight  =((self.phigh-self.plow)*float(i)/float(self.grid_n1-1) + self.plow)
			dpdp1 = self.mult_prof(self.NPT, dpdpf, self.pedc, pedw, pweight)
			
			
			for j in range(self.grid_n2):
				jweight  =((self.jhigh-self.jlow)*float(j)/float(self.grid_n2-1) + self.jlow)
				zjz1 = self.mult_prof(self.NPT, zjzf, self.pedc, pedw, jweight)		
				prof_data1 = self.write_prof_data(zjz1,dpdp1,fdff,nef,tef,tif,self.NPT)
				b = 0.15
				phys_data1 = self.write_phys_data(xiab,bvac,rvac,eps,zeff,zimp,corep,self.coreneo,b)
				
				fname = 'equ_scan_'+str(i+1)+'_'+str(j+1)		
				self.Make_HEL_inp(equ_save_dir,hel_shape,prof_data1,phys_data1,hel_num, b, fname)

		print('>>> Equilibrium run is finished')
		
		return()
		
	def read_hel1(self,filename,npts):
		f1 = open(filename,'r')
		finda = -1
		while (finda < 0):
			line = f1.readline()
			if not line: break
			finda = line.find('ALPHA')

		line = f1.readline()
		line = f1.readline()

		alpha = np.zeros(shape=(npts-1,1))
		alpha2 = np.zeros(shape=(npts-1,1))
		alpha3 = np.zeros(shape=(npts-1,1))
		PSIN = np.zeros(shape=(npts-1,1))
		jphi = np.zeros(shape=(npts-2,1))
		qprof = np.zeros(shape=(npts-1,1))
		for i in range(npts-1):
			line = f1.readline()
			alpha[i] = line.split()[6]
			alpha2[i] = line.split()[7]
			try:
				alpha3[i] = line.split()[10]
			except:
				alpha3[i] = 100.0
			if (alpha3[i] == 0.0):
				alpha3[i] = 100
			qprof[i] = line.split()[3]
			PSIN[i] = float(line.split()[1])

		finda = -1
		while (finda < 0):
			line = f1.readline()
			if not line: break
			finda = line.find('AVERAGE JPHI')
		line = f1.readline()
		for i in range(npts-2):
			line = f1.readline()
			jphi[i] = line.split()[1]

		count = 0
		target_psi = 0.8
		
		for i in range(npts-2):
			if ((PSIN[i]-target_psi)*(PSIN[i+1]-target_psi) < 0):
				count = i;	
			
		AMAX = max(alpha[105:npts-2])
		AMAX2 = max(alpha2[count:npts-2])
		AMAX3 = min(alpha3[count:npts-2])
		JMAX = max(jphi[count:npts-3])
		JMIN = jphi[-1]
		qedge = qprof[npts-2]
		f1.close()
		return (AMAX, AMAX2, AMAX3, JMAX, qedge, JMIN)

	def read_hel2(self,filename,npts):
		f1 = open(filename,'r')
		finda = -1
		while (finda < 0):
			line = f1.readline()
			if not line: break
			finda = line.find('ALPHA')

		line = f1.readline()
		line = f1.readline()

		alpha = np.zeros(shape=(npts-1,1))
		alpha2 = np.zeros(shape=(npts-1,1))
		ball = np.zeros(shape=(npts-1,1))
		PSIN = np.zeros(shape=(npts-1,1))
		jphi = np.zeros(shape=(npts-2,1))
		qprof = np.zeros(shape=(npts-1,1))
		for i in range(npts-1):
			line = f1.readline()
			alpha[i] = line.split()[6]
			alpha2[i] = line.split()[7]
			try:
				ball[i] = float(line.split()[10])
			except:
				ball[i] = 100.0
			if (ball[i] == 0.0):
				ball[i] = 100
			qprof[i] = line.split()[3]
			PSIN[i] = float(line.split()[1])

		finda = -1
		while (finda < 0):
			line = f1.readline()
			if not line: break
			finda = line.find('AVERAGE JPHI')
		line = f1.readline()
		for i in range(npts-2):
			line = f1.readline()
			jphi[i] = line.split()[1]

		count = 0
		target_psi = 0.8
		
		for i in range(npts-2):
			if ((PSIN[i]-target_psi)*(PSIN[i+1]-target_psi) < 0):
				count = i;	
			
		AMAX = max(alpha[105:npts-2])
		AMAX2 = max(alpha2[105:npts-2])
		ballm = min(ball[count:npts-2])
		JMAX = max(jphi[count:npts-3])
		JMIN = jphi[-1]
		qedge = qprof[npts-2]
		f1.close()
		return (AMAX, AMAX2, JMAX, qedge, ballm, JMIN)
		
	def generate_output_helena(self):
	
		nw = int(self.grid_n1)
		nh = int(self.grid_n2)

		#---
		hel_temp = self.run_dir +'/HEL_TEMP'
		mis_temp = self.run_dir +'/MIS_TEMP'
		eli_temp = self.run_dir +'/ELI_TEMP'

		if (self.run_stab == 'mishka'):
			resname = self.run_dir + '/ja_result_' + self.RUN_ID+'_mis'
		elif (self.run_stab == 'elite'):
			resname = self.run_dir + '/ja_result_' + self.RUN_ID+'_eli'

		resname1 = self.run_dir + '/alp_omega_' + self.RUN_ID

		equ_dir = self.run_dir + '/EQU_HELENA/fort.10'
		equ_out_dir = self.run_dir + '/EQU_HELENA/fort.20'

		if (self.run_stab == 'mishka'):
			nn = np.size(self.modens.split(','))
		elif (self.run_stab == 'elite'):
			nn = np.size(self.modenes.split(','))

		f1 = open(resname,'w')
		f_line = self.RUN_ID + ' '+ str(nw) + ' ' + str(nh) + ' ' + str(nn)+ ' \n'
		f1.write(f_line)

		#---
		hel_shape,hel_prof,hel_phys,hel_num,xiab,bvac,rvac,eps = self.get_HELENA(equ_dir)
		NPT,NPTP = self.get_num_opt(hel_num)
		a1,a2,a3,a4,a5,a6 = self.read_hel1(equ_out_dir,NPT)

		f1.write('%12.6f %12.6f %12.6f %12.6f\n'%(a1,a2,a4,a6))

		for i in range(nw):
			for j in range(nh):
				savename = '/gr_' + self.RUN_ID
				if (self.run_stab == 'mishka'):
					stab_dir = mis_temp +'/'+str(i+1)+'_'+str(j+1)
				elif (self.run_stab == 'elite'):
					stab_dir = eli_temp +'/'+str(i+1)+'_'+str(j+1)
				os.chdir(stab_dir)
				if (self.run_stab == 'mishka'):
					if(self.del_data):
						self.run_wo_out('rm CASPLOT')
						self.run_wo_out('rm fort.12')
						self.run_wo_out('rm eig_str')
						self.run_wo_out('rm fort.20')
						self.run_wo_out('rm fort.22')
				elif (self.run_stab == 'elite'):
					if(self.del_data):
		#			self.run_wo_out('rm *.eq*')
						self.run_wo_out('rm *.out')
						self.run_wo_out('rm *.fun2d')
						self.run_wo_out('rm *.surf')
						self.run_wo_out('rm *.plas')
						self.run_wo_out('rm fort.76')

				os.chdir(self.run_dir)
				savename= stab_dir + savename
				f2 = open(savename,'r')
				linecount = 0
				while True:
					line = f2.readline()
					linecount = linecount + 1
					if not line: break
					if (linecount > nn): break
					if len(line.split())>15: f1.write(line)
					else:
						equ_dir = hel_temp +'/'+str(i+1)+'_'+str(j+1) + '/fort.20'
						equ_out_dir = hel_temp +'/'+str(i+1)+'_'+str(j+1) + '/fort.20'
						hel_shape,hel_prof,hel_phys,hel_num,xiab,bvac,rvac,eps = self.get_HELENA(equ_dir)
						NPT,NPTP = self.get_num_opt(hel_num)
						a1,a2,a3,a4,a5,a6 = self.read_hel1(equ_out_dir,NPT)
						f1.write(line.split('\n')[0]+' %f \n'%(a6/1.e6))
				f2.close()
		f1.close()

		print('>>> mapping r #',NPT,'mapping t #',NPTP)
		if(NPT < 1 or NPTP < 1):
			print('>>> failed to read grid #')
			exit()

		f1 = open(resname1,'w')
		for i in range(nw):
			for j in range(nh):
				if not self.batch_run:
					update_progress(float((i*nh+j+1)/nw/nh))
				else:
					print(i,j)
				work_dir = hel_temp + '/'+str(i+1)+'_'+str(j+1)
				os.chdir(work_dir)
				elitename = work_dir + '/eliteinp2'
	#			print(elitename)
				npsi,psi,dpdpsi,fdf,te,ne,ti,zeff,zimp,alpw = self.read_elite(elitename)
				if(self.del_data):
					self.run_wo_out('rm eliteinp')
					self.run_wo_out('rm fort.12')
					self.run_wo_out('rm RZ_psi')
					self.run_wo_out('rm fort.21')
					self.run_wo_out('rm PCUBEQ')

				try:
					a1,a2,a3,a4,a5,a6 = self.read_hel1('fort.20',NPT)
				except:
					self.run_wo_out('cp tempin fort.10')
					os.system(self.helena_dir)
					a1,a2,a3,a4,a5,a6 = self.read_hel1('fort.20',NPT)
					if(self.del_data):
						self.run_wo_out('rm eliteinp')
						self.run_wo_out('rm fort.12')
						self.run_wo_out('rm RZ_psi')
						self.run_wo_out('rm fort.21')
						self.run_wo_out('rm PCUBEQ')

				f1.write('%i %i %12.6f %12.6f \n'%(i+1,j+1,alpw,a3))
				os.chdir(self.run_dir)
		f1.close()	
	
	def new_helinp(self,filename,xiab,btor,initB):
		f1 = open(filename,'r')
		f2 = open('fort.10','w')
		while True:
			line = f1.readline()
			if not line: break
			line2 = line.split('=')[0].split()[0]
			if ( line2 == 'XIAB'):
				xiab2 = float(line.split('=')[1])
				line =  'XIAB = ' + str(float(xiab)) + '\n'
			elif ( line2 == 'BVAC'):
				btor2 = float(line.split('=')[1])
				line = 'BVAC = ' + str(float(btor)) + '\n'
			elif ( line2 == 'B' or line2 == 'b'):
				initB2 = float(line.split('=')[1])
				line = 'B = ' + str(float(initB)) + '\n'
			f2.write(line)
		f1.close()
		f2.close()
		return (xiab2,btor2,initB2)		
		
	def write_mishka_inp(self,filename,ntor,gridn,psis,xr1,sig1,xr2,sig2,ias):

		f4 = open(filename,'w')

		f4.write('&NEWRUN \n')
		f4.write("EQNAME = 'jet', \n")
		f4.write('MODE=4, NLTORE = .T., \n')
		f4.write('NG=%i, \n'%(gridn))
		f4.write('RFOUR(1)=0, \n')
		f4.write('VFOUR(1)=0, \n')
		f4.write('IFAST =1 \n')
		f4.write('NTOR= -%i \n'%(ntor))
		f4.write('VSHIFT(1) = (2.e-1, 0), \n')
		f4.write('VSHIFT(2) = (2e-2, 0), \n')
		f4.write('VSHIFT(3) = (2e-3, 0), \n')
		f4.write('VSHIFT(4) = (4e-4, 0), \n')
		f4.write('NSHIFT=4, \n')
		f4.write('Q0ZYL = -1, \n')
		f4.write('GAMMA = 0., \n')
		f4.write('ASPECT = 1., \n')
		f4.write('IEQ = 2, \n')
		f4.write('sbegin=%4.3f \n'%(psis))
		f4.write('SIG1 = %4.3f, SIG2 = %4.3f, \n'%(sig1,sig2))
		f4.write('XR1 = %4.3f, XR2 = %4.3f, \n'%(xr1,xr2))
		f4.write('RWALL = 10.0, IVAC=1, \n')
		f4.write('NVPSI = 51, NGV = 51, \n')
		f4.write('SIGV = 0.1, \n')
		if (ias):
			f4.write('IAS=1, \n')
		else:
			f4.write('IAS=0, \n')
			
		f4.write('RMIN=0.890026, ZCNTR=-0.29906, \n')
		f4.write('&END \n')
		f4.write('\n')
		f4.write('&NEWRUN  MODE=0 &END \n')
		f4.write('&NEWLAN    &END \n')

		f4.close()
	
		return

	def write_elite_inp(self,filename,ntor,qdelfix,ngrid,npmap,psise,ndist,ias):
 
		f4 = open(filename,'w')

		f4.write('&equil \n')
		f4.write("shape='eqbm'\n")
		f4.write('delmin=.f.\n')
		f4.write('nnmin=96 \n')
		f4.write('nnmax=110 \n')
		f4.write('setdel=.t.\n')
		f4.write('qafix=.f.\n')
		f4.write('del_fix=%3.2f\n'%(qdelfix))
		f4.write('percenflux=1.0\n')
		f4.write('alpsi=-.98\n')
		f4.write('npts=%i\n'%(npmap))
		f4.write('&end\n')
		f4.write('&vac\n')
		f4.write('npts=%i \n'%(ngrid))
		f4.write('ng=8\n')
		f4.write('&end\n')
		f4.write('&qref_modes\n')
		f4.write('nn=%i\n'%(ntor))
		f4.write('nmlow=3\n')
		f4.write('nmvac=-5\n')
		f4.write('psimin=%4.3f\n'%(psise))
		f4.write('nmwinhalf=-1\n')
		f4.write('nxinterp=155\n')
		f4.write('dens=.true.\n')
		if (self.use_rot):
			f4.write('rotation=.t. \n')
		f4.write('&end\n')
		f4.write('&plas\n')
		f4.write('dx=0.01\n')
		f4.write('nd=20\n')
		f4.write('dw=0.15\n')
		f4.write('ndist=%i\n'%(ndist))
		f4.write('lamdist=0.3\n')
		f4.write('n1term=2\n')
		f4.write('gamsq=0.01\n')
		f4.write('igam=20\n')
		f4.write('dmercier=0.0\n')
		f4.write('ns=2500\n')
		f4.write('meshtype=3\n')
		f4.write('autorun=.f.\n')
		f4.write('funcal=.t.\n')
		if (ias):
			f4.write('updownsym=.f.\n')
		else:
			f4.write('updownsym=.t.\n')
			
		f4.write('vacuum=.t.\n')
		if (self.use_comp):
			f4.write('compression=.t. \n')
		if (self.use_rot):
			f4.write('rot_model=-1 \n')
			f4.write('rot_sep = %f \n'%self.rot_sep)
			f4.write('rot_const = %f \n'%self.rot_const)
			f4.write('rot_pedmid = %f \n'%self.rot_pedmid)
			f4.write('rot_pedwid = %f \n'%self.rot_pedwid)
			f4.write('rot_axis = %f \n'%self.rot_axis)
			f4.write('rot_expin = %f \n'%self.rot_expin)
			f4.write('rot_expout = %f \n'%self.rot_expout)

		f4.write('nowindow=.f.\n')
		f4.write('bug(1)=0.0\n')
		f4.write('&end\n')

		f4.close()
		return

	def run_wo_out(self,command):

		try:
			stat,out = subprocess.getstatusoutput(command)
		except:
			pass
		return
		
	def read_chease_bnd(self,filename):
	
		f4 = open(filename,'r')
		
		bnd_data = []
		
		for i in range(3):
			line = f4.readline()
			bnd_data.append(line)
		line = f4.readline()	
		bnd_num = int(line.split()[0])
		bnd_data.append(line)
		for i in range(bnd_num):
			line = f4.readline()
			bnd_data.append(line)
			
		f4.close()	
		self.bnd_data = bnd_data
		return
		
	def read_chease_out(self,filename):
	
		f4 = open(filename,'r')
		dat_num = int(f4.readline().split()[0])
		self.psin = np.zeros(dat_num)
		self.pprime = np.zeros(dat_num)
		self.q = np.zeros(dat_num)
		self.jav1 = np.zeros(dat_num)
		self.jav2 = np.zeros(dat_num)
		self.alp1 = np.zeros(dat_num)
		self.alp2 = np.zeros(dat_num)
		self.ball = np.zeros(dat_num)
		
		for i in range(2):
			line = f4.readline()
		self.IPEXP = float(f4.readline().split('=')[1].split()[0])
		self.R0EXP = float(f4.readline().split('=')[1].split()[0])
		self.B0EXP = float(f4.readline().split('=')[1].split()[0])
		self.ASPCT = float(f4.readline().split('=')[1].split()[0])
		line = f4.readline()
		for i in range(dat_num):
			line = f4.readline().split()
			self.psin[i] = float(line[1])
			self.pprime[i] = float(line[2])
			self.q[i] = float(line[3])
			self.jav1[i] = float(line[4])
			self.jav2[i] = float(line[5])
			self.alp1[i] = float(line[6])
			self.alp2[i] = float(line[7])
			self.ball[i] = float(line[8])
			
		return
			
	def write_chease_expeq(self,filename,psin,pprime,jav,B0,IP):
	
		if (filename == None):
			f4 = open('EXPEQ','w')
		else:
			f4 = open(filename,'w')
		
		for item in self.bnd_data:
			f4.write(item)

		f4.write('%i \n'%len(psin))
		f4.write('2 0 \n')
		for i in range(len(psin)):
			f4.write('%f \n'%(np.sqrt(psin[i])))
		for i in range(len(psin)):
			f4.write('%f \n'%(-1.0*pprime[i]*abs(pprime[100])/pprime[100]*4*np.pi*1e-7*self.R0EXP*self.R0EXP/B0))
		for i in range(len(psin)):
			f4.write('%f \n'%(abs(jav[i])*4*np.pi*1e-7*self.R0EXP / B0))

	#	f4.write('%f, R0 [M] \n'%self.R0EXP)
#		f4.write('%f, B_T [T] \n'%B0)
#		f4.write('%f, TOTAL CURRENT \n'%(IP*4*np.pi*1e-7/self.R0EXP / B0))
		
		f4.close()
		
		return

	def modify_chease_namelist(self,filename,ns,nt,npsi,nchi):
	
		copyfile(filename,'tempin')
		f4 = open(filename,'r')
		chease_name = []
		count = 0
		while True:
			line = f4.readline()
			if not line: break
			if (line.find('/') == -1):
				chease_name.append(line)
		f4.close()
		
		f4 = open(filename,'w')
		for item in chease_name:
			f4.write(item)
		f4.write('NS=%i,	NT=%i, \n'%(ns,nt));
		f4.write('NPSI=%i, NCHI=%i, \n'%(npsi,nchi));
		f4.write('NIDEAL = 8\n')
		f4.write('/\n');
		
		f4.close()
					
		return

	def modify_chease_namelist2(self,filename):
	
		copyfile(filename,'tempin')
		f4 = open(filename,'r')
		chease_name = []
		count = 0
		while True:
			line = f4.readline()
			if not line: break
			if (line.find('/') == -1):
				chease_name.append(line)
		f4.close()
		
		f4 = open(filename,'w')
		for item in chease_name:
			f4.write(item)
		f4.write('NIDEAL = 8\n')
		f4.write('/\n');
		
		f4.close()
					
		return		
	
	def make_chease_namelist(self):
		# default values
		nideal = 8

		chease_namelist1 = []
		chease_namelist1.append('*************************\n');
		chease_namelist1.append('***    CHEASE NAMELIST FILE\n');
		chease_namelist1.append('***    namelist created by write_namelist_chease.m\n');
		chease_namelist1.append('*************************\n');
		chease_namelist1.append('&EQDATA\n');
		chease_namelist1.append('! --- CHEASE input file option\n');
		chease_namelist1.append('NOPT=0,\n');
		chease_namelist1.append('NSURF=6,\n');
		chease_namelist1.append('NEQDSK=0 \n');
		chease_namelist1.append('TENSBND=    -0.1,\n');
		chease_namelist1.append('NSYM=0, \n');
		chease_namelist1.append('NTCASE=0, \n');
		chease_namelist1.append('RELAX = 0.,\n');
		chease_namelist1.append('NPLOT=0,\n');
		chease_namelist1.append('NSMOOTH=1,\n');
		chease_namelist1.append('NDIAGOP=1,\n');
		chease_namelist1.append('NPROPT=2,\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- Solving opt\n');
		chease_namelist1.append('SIGNB0XP=       1,\n');
		chease_namelist1.append('SIGNIPXP=       1,\n');
		chease_namelist1.append('TREEITM="euitm", "euitm",\n');
		chease_namelist1.append('XI=       0,\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- Profile opt\n');
		chease_namelist1.append('NFUNRHO=0,\n');
		chease_namelist1.append('TENSPROF=    -0.1,\n');
		chease_namelist1.append('NPPFUN = 4, NPP=2,\n');
		chease_namelist1.append('NFUNC=4, NIPR=2, \n');
		chease_namelist1.append('CFNRESS=1.00,\n');
		chease_namelist1.append('NRSCAL=0, NTMF0=0,\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- MESH GRID distributing option\n');
		#chease_namelist1.append('NMESHA=1, NPOIDA=3, SOLPDA=.20,APLACE=0.9, .98, 1.0,\n');
		#chease_namelist1.append('                               AWIDTH=0.1, .03,.001,\n');
		chease_namelist1.append('NMESHPOL=0,SOLPDPOL=0.25,\n');
		#chease_namelist1.append('NMESHC=0, NPOIDC=2, SOLPDC=.70,CPLACE=.95,.98,1.0,\n');
		#chease_namelist1.append('                               CWIDTH=.10,.02,.05,\n');
		width = 2.*(1.001 - self.pedc)	
		chease_namelist1.append('NMESHA=1, NPOIDA=2, SOLPDA=.20,APLACE= %4.3f, 1.000,\n'%(1.-0.5*width)); #####mesh
		chease_namelist1.append('                               AWIDTH= %4.3f, 0.003,\n'%(0.9*width));
		chease_namelist1.append('NMESHC=1, NPOIDC=3, SOLPDC=.20,CPLACE=0.000, %4.3f, 1.00,\n'%(1.-0.5*width));
		chease_namelist1.append('                               CWIDTH=0.005, %4.3f, 0.003,\n'%(0.9*width));
		chease_namelist1.append('NMESHD=0, NPOIDD=2, SOLPDD=.60,DPLACE=-1.80,-1.80,\n');
		chease_namelist1.append('                               DWIDTH=.18,.08,.05,\n');
		chease_namelist1.append('NMESHE=0, NPOIDE=4, SOLPDE=.50,EPLACE=-1.70,-1.70,\n');
		chease_namelist1.append('                               EWIDTH=.18,.08,.18,\n');
		chease_namelist1.append('PSISCL= 1.0, \n');
		chease_namelist1.append('NDIFPS=0,\n');
		chease_namelist1.append('NDIFT=0, \n');
		chease_namelist1.append('NINMAP=50,  NINSCA=50,\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! ---  Define Jacobian\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- Bootstrap current input\n');
		chease_namelist1.append('NBSEXPQ=0,\n');
		chease_namelist1.append('ETAEI= 0.1, RPEOP= 0.5,  RZION=1.5,\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- Ballooning optimization option including\n');
		chease_namelist1.append('NBAL=0, NBLOPT=0, \n');
		chease_namelist1.append('NPPR=30,\n');
		chease_namelist1.append('NTURN=20, NBLC0=16, \n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- Bootstrap current constraint optimization\n');
		chease_namelist1.append('NBSOPT=0, \n');
		chease_namelist1.append('NBSTRP=1, BSFRAC=0, NBSFUN=1, \n');
		chease_namelist1.append('EPSLON=%e, GAMMA=1.6666666667,\n'%self.epslon);
		chease_namelist1.append('CFBAL=10.00,\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- Eqdsk Write option \n');
		chease_namelist1.append('COCOS_IN = 2,\n');
		chease_namelist1.append('COCOS_OUT = 2,\n');
		chease_namelist1.append('NRBOX=129,\n');
		chease_namelist1.append('NZBOX=129,\n');
		chease_namelist1.append('NEQDXTPO=2,\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- Stability input (ELITE & MISHKA, NIDEAL =8 ) \n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- Stability input (MARS, NIDEAL=0,3) \n');
		chease_namelist1.append('NVEXP=1, REXT=10.0, R0W=1., RZ0W=0.,\n');
		chease_namelist1.append('MSMAX=1,  \n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- Stability NOVA & XTOR CODE input \n');
		chease_namelist1.append('NTNOVA=12,\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('! --- ITM opt (not used) \n');
		chease_namelist1.append('NITMRUN=1, 99, NITMSHOT=10, 10,\n');
		chease_namelist1.append('NITMOPT=0,\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('\n');
		chease_namelist1.append('NCSCAL=2,\n');
		chease_namelist1.append('CSSPEC=0,\n');
		chease_namelist1.append('! --- CHEASE plasma opt \n');
		chease_namelist1.append('NEGP=%i, NER=%i, \n'%(0,2));
		chease_namelist1.append('NIDEAL=%i, \n'%(nideal));

		return (chease_namelist1)
		
	def write_chease_namelist(self,filename,ns,nt,npsi,nchi,B0,IP):
	
		namelist = self.make_chease_namelist()
		if (filename == None):
			f4 = open('chease_namelist','w')
		else:
			f4 = open(filename,'w')
		
		for item in namelist:
			f4.write(item)
		f4.write('NS=%i,	NT=%i, \n'%(ns,nt));
		f4.write('NPSI=%i, NCHI=%i, \n'%(npsi,nchi));
		f4.write('ASPCT= %f,\n'%self.ASPCT);
		f4.write('R0EXP= %f, B0EXP=%f, \n'%(self.R0EXP,abs(B0)));
		f4.write('CURRT=%f, QSPEC=%f, \n'%(abs(IP*4*np.pi*1e-7 / self.R0EXP / B0),self.q[0]));
		f4.write('/\n');
		
		f4.close()
		
		return
	
	def run_chease(self):
	
		try:
			stat,out = subprocess.getstatusoutput(self.chease_dir+' > log.chease')
		except:
			print('Failed to run Chease')
			pass
		print(out)
		
	def find_max_ja(self,psint,psin,tarray,max_min):
	
		for i in range(len(psin)-1):
			if (psin[i] > psint): break

		if (max_min == 0):
			maxv = min(tarray[i:len(tarray)+1])
		else:
			maxv = max(tarray[i:len(tarray)+1])
		
		return maxv
		
	def make_scan_equilibrium_chease(self):
		
		# --- run_dir
		work_dir = self.run_dir + '/EQU_CHEASE'
		equ_save_dir = self.run_dir + '/equ_inp'

		self.make_directory(work_dir,'work')
		self.make_directory(equ_save_dir,'equ save')

		print('>>> Calculate given equilibrium')

		copyfile(self.init_HELinp,work_dir+'/EXPEQ')
		copyfile(self.init_CHEinp,work_dir+'/chease_namelist')
		if (self.use_kin_prof):
			if not (self.equ_type ==1):
				copyfile(self.kin_prof_name,work_dir+'/chease_kinprof')
				copyfile(self.kin_prof_name,equ_save_dir+'/chease_kinprof')
			elif (self.adjust_prof):
				copyfile(self.kin_prof_name+'_new',work_dir+'/chease_kinprof')
				copyfile(self.kin_prof_name+'_new',equ_save_dir+'/chease_kinprof')
			else:
				copyfile(self.kin_prof_name,work_dir+'/chease_kinprof')
				copyfile(self.kin_prof_name,equ_save_dir+'/chease_kinprof')
				
		if (self.use_rot):
			copyfile(self.rot_prof_name,work_dir+'/chease_vtor')
			copyfile(self.rot_prof_name,equ_save_dir+'/chease_vtor')

		os.chdir(work_dir)
		self.modify_chease_namelist('chease_namelist',self.CNS,self.CNT,self.CNPSI,self.CNCHI)

		CHEout = work_dir + '/NJA'
		ELIout = work_dir + '/NELITE'
		
		self.read_chease_bnd('EXPEQ')
		
		try:
			self.read_chease_out(CHEout)
			yes_prev = True
		except:
			yes_prev = False

		if not(self.use_prev and yes_prev):
			self.run_chease()
			self.read_chease_out(CHEout)
		else:
			print('>>> use previous CHEASE namelist')

		npsi,psi,dpdpsi,fdf,te,ne,ti,zeff,zimp,alf_over_maxomega = self.read_elite(ELIout)

		if (self.use_ped_find):
			dpdpf = interp1d(self.psin,abs(self.pprime),'cubic')
			psi_temp = np.linspace(0.75,1.0,1001)
			dpdp_temp = dpdpf(psi_temp)
			ind = np.argmax(dpdp_temp)
			print('pedestal_center %5.4f -> %5.4f'%(self.pedc,psi_temp[ind]))
			self.pedc = psi_temp[ind]		

		pedw = 1.01 - self.pedc
		print('>>> pedestal weight option...')
		print('half width',round(pedw,5),'pedestal center',self.pedc)

		print('>>> Generate equilibrium files for scan')

		dpdp1 = np.zeros(len(self.pprime))
		zjz1 = np.zeros(len(self.pprime))
		
		for i in range(self.grid_n1):

			pweight  =((self.phigh-self.plow)*float(i)/float(self.grid_n1-1) + self.plow)
			
			for k in range(len(self.pprime)):
				weig = self.mult_prof2(self.psin[k], self.pedc, pedw, pweight)
				dpdp1[k] = weig * self.pprime[k]
			
			for j in range(self.grid_n2):
				jweight  =((self.jhigh-self.jlow)*float(j)/float(self.grid_n2-1) + self.jlow)
				
				for l in range(len(self.pprime)):
					weig = self.mult_prof2(self.psin[l], self.pedc, pedw, jweight)	
					zjz1[l] = weig * self.jav2[l]
				
				fname = equ_save_dir + '/equ_scan_'+str(i+1)+'_'+str(j+1)		
				fname_expeq = fname + '_expeq'
				fname_name = fname + '_name'

				self.write_chease_namelist(fname_name,self.CNS,self.CNT,self.CNPSI,self.CNCHI,self.B0EXP,self.IPEXP)
				self.write_chease_expeq(fname_expeq,self.psin,dpdp1,zjz1,self.B0EXP,self.IPEXP)

		print('>>> Equilibrium run is finished')
		
		return()
		
	def read_kin_prof(self,filename):
	
		f4 = open(filename,'r')
		dat_num = int(f4.readline().split()[0])
		self.PSIX = np.zeros(dat_num)
		self.NE = np.zeros(dat_num)
		self.TE = np.zeros(dat_num)
		self.NI = np.zeros(dat_num)
		self.TI = np.zeros(dat_num)

		line = f4.readline()
		
		self.zeff  = float(line.split()[0])
		self.zimp  = float(line.split()[1])
		self.amain = float(line.split()[2])
		self.aimp  = float(line.split()[3])
		
		for i in range(dat_num):
			line = f4.readline()
			self.PSIX[i] = float(line.split()[0])
			self.TE[i] = float(line.split()[1])
			self.NE[i] = float(line.split()[2])
			self.TI[i] = float(line.split()[3])
			self.NI[i] = float(line.split()[4])
	
		return
	
	def read_rot_prof(self,filename):
		
		f4 = open(filename,'r')
		dat_num = int(f4.readline().split()[0])
		self.RR = np.zeros(dat_num)
		self.VT = np.zeros(dat_num)
		self.FROT = np.zeros(dat_num)
		self.PSI_ROT = np.zeros(dat_num)

		line = f4.readline()

		for i in range(dat_num):
			line = f4.readline()
			self.PSI_ROT[i] = float(line.split()[0])
			self.RR[i] = float(line.split()[1])
			self.VT[i] = float(line.split()[2])
			self.FROT[i] = self.VT[i] / self.RR[i] / 2.0 / np.pi
		f4.close()

		self.VT = self.VT - self.VT[-1]
		self.FROT = self.FROT - self.FROT[-1]

		return
		
	def generate_output_chease(self):
	
		nw = int(self.grid_n1)
		nh = int(self.grid_n2)

		#---
		chease_temp = self.run_dir +'/CHE_TEMP'
		mis_temp = self.run_dir +'/MIS_TEMP'
		eli_temp = self.run_dir +'/ELI_TEMP'
		
		if (self.run_stab == 'mishka'):
			resname = self.run_dir + '/ja_result_' + self.RUN_ID+'_mis'
		elif (self.run_stab == 'elite'):
			resname = self.run_dir + '/ja_result_' + self.RUN_ID+'_eli'

		resname1 = self.run_dir + '/alp_omega_' + self.RUN_ID

		equ_dir = self.run_dir + '/EQU_CHEASE/'
		equ_out_dir = self.run_dir + '/EQU_CHEASE/NJA'

		if (self.run_stab == 'mishka'):
			nn = np.size(self.modens.split(','))
		elif (self.run_stab == 'elite'):
			nn = np.size(self.modenes.split(','))

		f1 = open(resname,'w')
		f_line = self.RUN_ID + ' '+ str(nw) + ' ' + str(nh) + ' ' + str(nn)+ ' \n'
		f1.write(f_line)

		#---
		self.read_chease_out(equ_out_dir)
		
		a1 = self.find_max_ja(0.8,self.psin,self.alp1,1)
		a2 = self.find_max_ja(0.8,self.psin,self.alp2,1)
		a4 = self.find_max_ja(0.8,self.psin,abs(self.jav1),1)
		a4e= abs(self.jav1[-1])

		f1.write('%12.6f %12.6f %12.6f %12.6f \n'%(a1,a2,a4,a4e))

		for i in range(nw):
			for j in range(nh):
				savename = '/gr_' + self.RUN_ID
				if (self.run_stab == 'mishka'):
					stab_dir = mis_temp +'/'+str(i+1)+'_'+str(j+1)
				elif (self.run_stab == 'elite'):
					stab_dir = eli_temp +'/'+str(i+1)+'_'+str(j+1)
				
				os.chdir(stab_dir)
				if (self.run_stab == 'mishka'):
					if(self.del_data):
						self.run_wo_out('rm CASPLOT')
						self.run_wo_out('rm fort.12')
						self.run_wo_out('rm eig_str')
						self.run_wo_out('rm fort.20')
						self.run_wo_out('rm fort.22')
				elif(self.run_stab == 'elite'):
					if(self.del_data):
						self.run_wo_out('rm *.out')
						self.run_wo_out('rm *.fun2d')
						self.run_wo_out('rm *.surf')
						self.run_wo_out('rm *.plas')
						self.run_wo_out('rm fort.76')

				os.chdir(self.run_dir)
				savename= stab_dir + savename
				f2 = open(savename,'r')
				while True:
					line = f2.readline()
					if not line: break
					if len(line.split())>15: f1.write(line)
					else:
						equ_out_dir = chease_temp + '/' + str(i+1)+'_'+str(j+1) + '/NJA'
						self.read_chease_out(equ_out_dir)
						f1.write(line.split('\n')[0]+' %f \n'%(self.jav1[-1]/1.e6))
				f2.close()
		f1.close()

		f1 = open(resname1,'w')
		for i in range(nw):
			for j in range(nh):
				work_dir = chease_temp + '/'+str(i+1)+'_'+str(j+1)
				os.chdir(work_dir)
				elitename = work_dir + '/NELITE2'
				if not self.batch_run:
					update_progress(float((i*nh+j+1)/nw/nh))
				else:
					print(i,j)
				try:
					npsi,psi,dpdpsi,fdf,te,ne,ti,zeff,zimp,alpw = self.read_elite(elitename)
				except:
					elitename = work_dir + '/NELITE'
					copyfile('NELITE','NELITE2')
					npsi,psi,dpdpsi,fdf,te,ne,ti,zeff,zimp,alpw = self.read_elite(elitename)
					elitename = work_dir + '/NELITE2'
				self.run_wo_out('rm NELITE')
				self.run_wo_out('rm NMISHKA NOUT NDES EQDSK_*')

				try:
					self.read_chease_out(equ_out_dir)
					a3 = self.find_max_ja(0.8,self.psin,self.ball,0)
				except:
					copyfile('EXPEQ_temp','EXPEQ')
					copyfile('chease_namelist_temp','chease_namelist')
					sim.run_chease()
					self.read_chease_out(equ_out_dir)
					a3 = self.find_max_ja(0.8,self.psin,self.ball,0)
					self.run_wo_out('rm NELITE')


				f1.write('%i %i %12.6f %12.6f \n'%(i+1,j+1,alpw,a3))
				os.chdir(self.run_dir)
		f1.close()	

	def qval_change_chease_namelist(self,filename,qloc,qval):

		copyfile(filename,'tempin')
		f4 = open(filename,'r')
		chease_name = []
		count = 0
		while True:
			line = f4.readline()
			if not line: break
			if (line.find('/') == -1):
				chease_name.append(line)
		f4.close()

		f4 = open(filename,'w')
		for item in chease_name:
			f4.write(item)
		f4.write('NCSCAL=1, CSSPEC=%f, \n'%(qloc));
		f4.write('QSPEC=%f , \n'%(qval))
		f4.write('/\n');

		f4.close()
		return

	def read_elite_eq(self,filename_eq,filename_edge):

		f4 = open(filename_eq,'r')
		monotonicq = True
		while True:
			line = f4.readline()
			if not line: break
			if (line.find('monotonic') > -1): 
				monotonicq = False
				break
		f4.close()
			
#		f4 = open(filename_edge,'r')
#		for i in range(7):
#			line = f4.readline()
#		line = f4.readline().split()
#		alpha = float(line[0])
#		jphi = float(line[3])	

		return monotonicq

	def eped_ped(self, x, a1, a2, a3, a4, a5, a6, a7):

		y = a1
		y = y + a2*(np.tanh((1.0 - abs(a3))/(0.01+abs(a4))) - np.tanh((x - abs(a3))/(0.01+abs(a4))))
		yt = (x/(abs(a3)-abs(a4)-0.01)) ** (a6)
		y = y + a5 * (abs((1-yt)) **(abs(a7))) * 0.5 * (1.0 + np.sign(abs(a3)-abs(a4) - 0.01 -x))
	
		return (y)

	def elite_rot_prof(self,psin,f_rot):

		for i in range(len(psin)):
			if psin[i] > min(0.6,(self.psise-0.05)):
				break

		popt, pcov = curve_fit(self.eped_ped,psin[i:len(psin)+1],f_rot[i:len(psin)+1],p0=[f_rot[-1],0.5*f_rot[0],0.95,0.05,0.5*f_rot[0],1.1,1.1],maxfev=30000)
	
		self.rot_sep = popt[0]
		self.rot_const = popt[1] * (np.tanh(1.0) + np.tanh((1-abs(popt[2]))/(0.01+abs(popt[3])))) + self.rot_sep
		self.rot_pedmid = abs(popt[2])
		self.rot_pedwid = abs(popt[3])+0.01
		self.rot_axis = self.eped_ped(0.0,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],abs(popt[6]))
		self.rot_expin = popt[5]
		self.rot_expout = abs(popt[6])

		self.rot_fit = self.eped_ped(psin,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],abs(popt[6]))
	
		return

	def clear_equ_dat(self,clearall=False):

		os.chdir(self.run_dir +'/CHE_TEMP')

		nw = self.grid_n1
		nh = self.grid_n2
		
		for i in range(nw):
			for j in range(nh):
				os.chdir('%i_%i'%(i+1,j+1))
				try:
					os.remove('NDES')
				except:
					pass
				try:
					os.remove('NOUT')
				except:
					pass
				try:
					os.remove('NMISHKA')
				except:
					pass
				try:
					os.remove('NELITE')
				except:
					pass
				try:
					os.remove('EQDSK_COCOS_02.OUT')
				except:
					pass
				try:
					os.remove('EQDSK_COCOS_02_POS.OUT')
				except:
					pass
				if clearall:
					try:
						os.remove('NELITE2')
					except:
						pass					
						
				os.chdir('../')
				update_progress((nh*i+j+1)/nw/nh)
		return
				

	def initialise_variables(self):

		self.init_HELinp = 'EXPEQ'
		self.init_CHEinp = 'chease_namelist'
		self.equ_type = 0
		self.run_equ = 'chease'
		self.RUN_ID = '1'
		self.grid_n1 = 8
		self.grid_n2 = 8
		self.run_dir = '/ja_' + self.RUN_ID
		self.pedc = 0.98
		self.use_ped_find = 'True'
		self.NPT = 301
		self.plow = -0.5
		self.phigh = 0.5
		self.jlow = -0.5
		self.jhigh = 0.5
		self.b = 0.5
		self.coreneo = 0.02
		self.use_prev = 'False'
		self.use_hag = 'False'
		self.use_neo = 'False'
		self.use_chang = 'True'
		self.qsub_dir = 'qsub'
		self.qstat_dir = 'qstat'
		self.nodelist = node_default
		self.python_dir = python3_exec
		self.batch_run = 'False'
		self.helena_dir = helena_exec
		self.chease_dir = chease_exec
		self.mis_dir = mis_exec
		self.elite_dir = elite_dir
		self.qdelfix = 0.3
		self.qdel_crit = 0.003
		self.highq = 'False'
		self.run_stab = 'mishka'
		self.ias = 'True'
		self.psis = 0.75
		self.xr2 = 0.9
		self.sig2 = 0.1
		self.xr1 = 1.0
		self.sig1 = 0.05
		self.ngride = 2000
		self.psise = 0.5
		self.ndist = 50
		self.modens = '5, 10, 15, 20, 25, 30'
		self.modenes = '5, 10, 15, 20, 25, 30'
		self.gridn = 301
		self.CNS = 70
		self.CNT = 70
		self.CNPSI = 200
		self.CNCHI = 150
		self.CNSE = 150
		self.CNTE = 200
		self.CNPSIE = 300
		self.CNCHIE = 512
		self.use_kin_prof = 'True'
		self.kin_prof_name = 'chease_kinprof'
		self.rot_prof_name = 'chease_rot'
		self.use_prescribed_helena = 'True'
		self.use_rot = 'False'
		self.use_comp = 'False'
		self.del_data = 'True'
		self.beta_criterion = 1
		self.adjust_prof = 'True'
		self.line_den = 0.
		self.target_bnd = 0.995
		self.epsilon = 1.e-8
		self.beta_crit = 0.01
		self.li_crit = 0.02
		self.ip_crit = 1.e-5
		self.bs_crit = 1.e-5
		self.iternc = 25
		self.iternl = 10
		self.relax = 0.8
		self.li_target = 1.0
		self.beta_target = 1.0
		self.use_li = 'True'
		self.use_li2 = 'False'
		self.use_core_mod = 'True'
		self.use_core_mod_psin = 0.3
		self.bsmulti = 1.0
		self.apf = 1.0
		self.bpf = 1.0
		self.cpf = 1.1
		self.dpf = 1.5
		self.ajf = 0.2
		self.bjf = 1.0
		self.cjf = 2.0
		self.djf = 2.0
		self.nqa = 27.7
		self.ncut = 1
		self.grcrit1 = 0.03
		self.grcrit2 = 0.25
		self.use_adja = 'True'
		self.use_bilin = 'True'
		self.rmag = 1.84
		self.epslon = 1.e-7

		return

	def make_chease_opt(self,filename='chease_opt'):

		f = open(filename,'w')

		f.write('!-- Plasma geometry.\n')
		f.write('EQDSK = geqdsk \n')
		f.write('BND_PSIN = %f\n'%self.target_bnd)
		f.write('\n')
		f.write('!-- Plasma 0-d parameters \n')
		if (self.beta_criterion == 1):
			f.write('Beta_val = %f \n'%self.beta_target)
		else:
			f.write('Beta_val = %f \n'%self.rmag)
		f.write('LITARGET = %f \n'%self.li_target)
		f.write('		 \n')
		f.write('!-- Fast ion profile option \n')
		f.write('APF = %f \n'%self.apf)
		f.write('BPF = %f \n'%self.bpf)
		f.write('CPF = %f \n'%self.cpf)
		f.write('DPF = %f \n'%self.dpf)
		f.write('		 \n')
		f.write('!-- Current profile option \n')
		if (self.use_neo):	
			f.write('USE_NEO = True \n')
		else:
			f.write('USE_NEO = False \n')
		if (self.use_hag):
			f.write('USE_HAGER = True \n')
		else:
			f.write('USE_HAGER = False \n')
		if (self.use_chang):
			f.write('USE_CHANG = True \n')
		else:
			f.write('USE_CHANG = False \n')
		if (self.use_core_mod):
			f.write('HAG_CORE_MOD=True \n')
		else:
			f.write('HAG_CORE_MOD = False \n')
		f.write('HAG_CORE_MOD_PSIN = %f\n'%self.use_core_mod_psin)
		f.write('Core_neo = %f \n'%self.coreneo)
		f.write('BSMULTI = %f \n'%self.bsmulti)
		f.write('AJF = %f \n'%self.ajf)
		f.write('BJF = %f \n'%self.bjf)
		f.write('CJF = %f \n'%self.cjf)
		f.write('DJF = %f \n'%self.djf)
		f.write('Current_ITERN = %i \n'%self.iternc)
		f.write('RELAX = %f \n'%self.relax)
		f.write('\n')
		f.write('!-- convergence option\n')
		f.write('Beta_crit = %f \n'%self.beta_crit)
		f.write('Li_crit = %f \n'%self.li_crit)
		f.write('Ip_crit = %f \n'%self.ip_crit)
		f.write('Bs_crit = %f \n'%self.bs_crit)
		f.write('Li_ITERN = %i \n'%self.iternl)
		f.write('\n')
		f.write('!-- default option (Do not adjust it!) \n')
		if (self.use_li and self.use_li2):
			f.write('RUN_MODE = EPED3 \n')			
		elif (self.use_li):
			f.write('RUN_MODE = EPED2 \n')
		else:
			f.write('RUN_MODE = NORMAL \n')
		if (self.beta_criterion == 1):
			f.write('Beta_criterion_type = 1 \n')
		else:
			f.write('Beta_criterion_type = 2 \n')
		if (self.adjust_prof):
			f.write('ADJUST_PROF = True \n')
		else:
			f.write('ADJUST_PROF = False \n')
		f.write('kinetic_profile_type = 1\n')
		f.write('chease_kinetic_file = chease_kinprof\n')
		f.write('NIDEAL = 8\n')

		f.close()

		return

	def __init__(self,filename):
	
		print('>>> Start JA-diagram plot tool...')
		print('>>> Read namelist...')
		self.batch_run = 'False'
		self.highq = 'False'
		self.ias = 'True'
		self.use_ped_find = 'True'
		self.use_ped_plot = 'False'
		self.use_rot = 'False'
		self.use_comp = 'False'
		self.del_data = 'True'
		self.beta_criterion = 1
		self.init_CHEinp = 'chease_namelist'
		self.line_den = 0.0
		self.target_bnd = 0.995
		self.adjust_prof = 'False'
		self.initialise_variables()
		self.read_namelist(filename)
		print('>>> Initialization complete...')		
