#!/usr/local/anaconda3/bin/python3

import numpy as np
import os, sys
import subprocess
from shutil import copyfile
import ch_tool as ch
from scipy.interpolate import interp1d
from exec_dirs import chease_exec,elite_dir,mis_exec

savename='growth_list'
currdir = os.getcwd()
equ_save_dir = '/CHEASE'
equ_save_dir = currdir + equ_save_dir
try:
	nw = int(sys.argv[1])
	nh = int(sys.argv[2])
except:
	pass

width = 0.05
try:
	width = float(sys.argv[3])
except:
	pass

class pedscan_stab:

	def str_to_logic(self,s):
		s = s.lower()
		if (s == 'true'):
			return True
		if (s == 'false'):
			return False
			
	def read_namelist(self,sim,filename='stab_input'):
		f4= open(filename,'r')

		while True:
			line = f4.readline()
			if not line: break

			self.CNS = sim.read_namelist_str(line,'NS',self.CNS,1)
			self.CNT = sim.read_namelist_str(line,'NT',self.CNT,1)
			self.CNPSI = sim.read_namelist_str(line,'MAP_NS',self.CNPSI,1)
			self.CNCHI = sim.read_namelist_str(line,'MAP_NT',self.CNCHI,1)
			
			self.CNSE = sim.read_namelist_str(line,'NSE',self.CNSE,1)
			self.CNTE = sim.read_namelist_str(line,'NTE',self.CNTE,1)
			self.CNPSIE = sim.read_namelist_str(line,'MAP_NSE',self.CNPSIE,1)
			self.CNCHIE = sim.read_namelist_str(line,'MAP_NTE',self.CNCHIE,1)
			
			self.gridn = sim.read_namelist_str(line,'MIS_GRIDN',self.gridn,1)
			self.psis = sim.read_namelist_str(line,'MIS_PSISTART',self.psis,2)
			self.xr1 = sim.read_namelist_str(line,'XR1',self.xr1,2)
			self.sig1 = sim.read_namelist_str(line,'SIG1',self.sig1,2)
			self.xr2 = sim.read_namelist_str(line,'XR2',self.xr2,2)
			self.sig2 = sim.read_namelist_str(line,'SIG2',self.sig2,2)
			
			self.ngride = sim.read_namelist_str(line,'ELI_GRIDN',self.ngride,1)
			self.psise = sim.read_namelist_str(line,'ELI_PSISTART',self.psise,2)
			self.ndist = sim.read_namelist_str(line,'ELI_NDIST',self.ndist,1)

			self.use_comp = sim.read_namelist_str(line,'USE_COMP',self.use_comp,4)
			self.use_rot = sim.read_namelist_str(line,'USE_ROT',self.use_rot,4)
			self.run_stab = sim.read_namelist_str(line,'RUNSTAB',self.run_stab,3)

			self.qdelfix = sim.read_namelist_str(line,'QDEL',self.qdelfix,2)
			self.epslon = sim.read_namelist_str(line,'EPSLON',self.epslon,2)
			self.highq = sim.read_namelist_str(line,'HIGHQ',self.highq,4)

			if (line.lower().find('moden') > -1):
				self.mode_line = line.split('=')[1]					

		f4.close()

		moden = self.mode_line.split(',')

		len2 = len(moden)
		self.moden = np.zeros(len2)

		for i in range(len2):
			self.moden[i] = moden[i]
			
		return	

	def make_directory(self):
		try:
			os.mkdir(equ_save_dir)
		except:
			pass
			
		return
		
	def initialise_var(self):
	
		self.CNS = 80
		self.CNT = 80
		self.CNPSI = 200
		self.CNCHI = 200
		
		self.CNSE = 150
		self.CNTE = 200
		self.CNPSIE = 300
		self.CNCHIE = 512
		
		self.gridn = 301
		self.psis = 0.75
		self.xr1 = 1.0
		self.sig1 = 0.07
		self.xr2 = 0.95
		self.sig2 = 0.05
		self.ias = True
		
		self.ngride = 2000
		self.psise = 0.6
		self.ndist = 30
		
		self.chease_dir = chease_exec
		self.rot_prof_name = 'chease_vtor'
		self.use_comp = False
		self.use_kin_prof = True
		self.use_rot = False
		self.elite_dir = elite_dir
		self.mis_dir = mis_exec
		self.run_stab = 'mishka'
		
		self.qdel_crit = 0.003
		self.qdelfix = 0.3
		self.highq = True
		self.epslon = 1.e-8
		
		self.mode_line = '5,7,10,15,20'			#str13

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
		
	def modify_chease_namelist(self,filename,ns,nt,npsi,nchi,nideal):
	
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
		f4.write('NIDEAL=%i, \n'%nideal)
		f4.write('NS=%i,	NT=%i, \n'%(ns,nt));
		f4.write('NPSI=%i, NCHI=%i, \n'%(npsi,nchi));
		f4.write('/\n');
		
		f4.close()
					
		return
		
	def run_chease(self):
		try:
			stat,out = subprocess.getstatusoutput(self.chease_dir+' > log.chease')
		except:
			pass
		print(out)
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
		sign = self.jav1[30] / abs(self.jav1[30])	
		self.jav1 = np.copy(self.jav1 * sign)		
		self.jav2 = np.copy(self.jav2 * sign)	
		
		return
		
	def find_max_ja(self,psint,psin,tarray,max_min):
	
		for i in range(len(psin)-1):
			if (psin[i] > psint): break

		if (max_min == 0):
			maxv = min(tarray[i:len(tarray)+1])
		else:
			maxv = max(tarray[i:len(tarray)+1])
		
		return maxv
		
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

		self.VT = self.VT - self.VT[-1]
		self.FROT = self.FROT - self.FROT[-1]

		f4.close()	
		return	

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
		
		f4.close()
		
		return

	def make_chease_namelist(self):
		# default values

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
		chease_namelist1.append('NMESHA=1, NPOIDA=2, SOLPDA=.20,APLACE= %4.3f, 1.000,\n'%(1.-0.5*self.width)); #####mesh
		chease_namelist1.append('                               AWIDTH= %4.3f, 0.003,\n'%(0.9*self.width));
		chease_namelist1.append('NMESHC=1, NPOIDC=3, SOLPDC=.20,CPLACE=0.000, %4.3f, 1.00,\n'%(1.-0.5*self.width));
		chease_namelist1.append('                               CWIDTH=0.005, %4.3f, 0.003,\n'%(0.9*self.width));
		chease_namelist1.append('NMESHPOL=0,SOLPDPOL=0.25,\n');
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

		return (chease_namelist1)
		
	def write_chease_namelist(self,filename,ns,nt,npsi,nchi,B0,IP,nideal):
	
		namelist = self.make_chease_namelist()
		if (filename == None):
			f4 = open('chease_namelist','w')
		else:
			f4 = open(filename,'w')
		
		for item in namelist:
			f4.write(item)
		f4.write('NIDEAL = %i,\n'%(nideal))
		f4.write('NS=%i,	NT=%i, \n'%(ns,nt));
		f4.write('NPSI=%i, NCHI=%i, \n'%(npsi,nchi));
		f4.write('ASPCT= %f,\n'%self.ASPCT);
		f4.write('R0EXP= %f, B0EXP=%f, \n'%(self.R0EXP,abs(B0)));
		f4.write('CURRT=%f, QSPEC=%f, \n'%(abs(IP*4*np.pi*1e-7 / self.R0EXP / B0),self.q[0]));
		f4.write('/\n');
		
		f4.close()
		
		return
		
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
		rho=ne[-1]*md   ############# ne[-1] -> ne[0]
#		rho=nmainion[0]*md
		omegan=2.9979e10/(ne*1.e-6)/4.8032e-10*(dne*1.e-14)*1.6022e-12*ti
		omegaspi=omegan+2.9979e10/(ne*1.e-6)/4.8032e-10*(dti*1.e-8)*1.6022e-12*(ne*1.e-6)
		omega_alf=np.sqrt(B*B/rho/rav/rav/mu0)#/q[0]/q[0])
		alf_over_maxomega=-omega_alf/np.min(omegaspi)#[0:200])
#		for i in range(len(psi)):
#			print(psi[i],omega_alf,omegaspi[i],rav)
		file.close()

		print('CHECK the definition of rho ne[0] or ne[1]?')

		return(npsi,psi,dpdpsi,ffp,te,ne,ti,zeff,zimp,alf_over_maxomega)

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
			
		f4 = open(filename_edge,'r')
		for i in range(7):
			line = f4.readline()
		line = f4.readline().split()
		alpha = float(line[0])
		jphi = float(line[3])	

		return monotonicq		
	
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

	def run_stability(self):
		self.run_stab = self.run_stab.lower()	
		if (self.run_stab == 'mishka'):
			try:
				stab_dir = currdir + '/MISHKA'
				os.mkdir(stab_dir)
			except:
				pass

		elif (self.run_stab == 'elite'):
			try:
				stabe_dir = currdir + '/ELITE'
				os.mkdir(stabe_dir)
			except:
				pass
	
		f3 = open(savename,'w')
	
	
		os.chdir(equ_save_dir)
		namelist_expeq = equ_save_dir + '/EXPEQ'
		namelist_name = equ_save_dir + '/chease_namelist'
		namelist_prof = currdir + '/chease_kinprof'
		namelist_rot = currdir + '/chease_vtor'

		copyfile(namelist_expeq,'EXPEQ_init')
		copyfile(namelist_name,'chease_namelist_init')

		if(self.use_kin_prof):
			copyfile(namelist_prof,'chease_kinprof')
		if(self.use_rot):
			copyfile(namelist_rot,'chease_vtor')

		self.read_chease_bnd('EXPEQ')
		if (self.run_stab == 'mishka'):
			self.modify_chease_namelist('chease_namelist',self.CNS,self.CNT,self.CNPSI,self.CNCHI,8)
		elif (self.run_stab == 'elite'):
			self.modify_chease_namelist('chease_namelist',self.CNSE,self.CNTE,self.CNPSIE,self.CNCHIE,8)
	
		
		self.run_chease()
		self.read_chease_bnd('EXPEQ.OUT')

		CHEout = equ_save_dir + '/NJA'

		copyfile('EXPEQ','EXPEQ_temp_init')
		copyfile('chease_namelist','chease_namelist_temp')
		copyfile('EXPEQ.OUT','EXPEQ')
		copyfile('EXPEQ.OUT','EXPEQ_temp')
		self.run_chease()
		self.read_chease_out(CHEout)

		qa0 = self.q[-1]
		qf = interp1d(self.psin,self.q,'cubic')
		q95 = qf(0.95)
		psin = np.copy(self.psin)
		jav = np.copy(self.jav2)
		pprime = np.copy(self.pprime)
		IPEXP = abs(self.IPEXP)
		B0EXP = abs(self.B0EXP)

		print('init qa0',qa0)

		jm0 = self.find_max_ja(0.8,self.psin,self.jav1,1)
		am0 = self.find_max_ja(0.8,self.psin,self.alp1,1)
		am20 = self.find_max_ja(0.8,self.psin,self.alp2,1)
		ball0 = self.find_max_ja(0.8,self.psin,self.ball,0)

	
		copyfile('NELITE','NELITE2')
		elitename = equ_save_dir  + '/NELITE2'

		if (self.run_stab.lower() == 'mishka'):
			nn = np.size(self.moden)
			print('SELECTED STABILITY CODE = MISHKA')
		elif (self.run_stab.lower() == 'elite'):
			nn = np.size(self.moden)
			print('SELECTED STABILITY CODE = RUN ELITE')
		if(self.use_rot):
			self.read_rot_prof(self.rot_prof_name)
			self.elite_rot_prof(self.PSI_ROT,self.FROT)

		for i in range(nn):
			os.chdir(equ_save_dir)
			if (self.run_stab == 'mishka'):

				n = self.moden[i]	
				qdel = np.ceil(qa0*n) - qa0*n
				mm = np.ceil(qa0*n)
			
				if (qdel > 0.8):
					mm = mm - 1
					qdel1 = mm - qa0*n
					print('qdel0, new qdel, m0, m1',qdel,qdel1,mm+1,mm)
				else:
					print('qdel, m',qdel,mm)

				qa_target = float(mm - self.qdelfix) / float(n)
				print('qa-target',qa_target)
			
				self.write_chease_namelist(None,self.CNSE,self.CNTE,self.CNPSIE,self.CNCHIE,B0EXP,IPEXP,8)	
				self.write_chease_expeq(None,psin,pprime,jav,B0EXP,IPEXP)
				self.qval_change_chease_namelist('chease_namelist',1.0,qa_target)
				self.run_chease()
				
				self.read_chease_out(CHEout)
				
				jm1 = self.find_max_ja(0.8,self.psin,self.jav1,1)
				am1 = self.find_max_ja(0.8,self.psin,self.alp1,1)
				am21 = self.find_max_ja(0.8,self.psin,self.alp2,1)
				ball1 = self.find_max_ja(0.8,self.psin,self.ball,0)
				qa1 = self.q[-1]
				qf = interp1d(self.psin,self.q,'cubic')
				q95 = qf(0.95)
				qdel1 = np.ceil(qa1*n) - qa1*n

				print('qdelfix,qdel1,qa0,qa1,n,mm')
				print(self.qdelfix,qdel1,qa0,qa1,n,mm)

				if(abs(self.qdelfix-qdel1)>self.qdel_crit):
					print('q converge fail')
					mm = mm + 1
					qa_target = float(mm - self.qdelfix) / float(n)
					print('qa-target',qa_target)
					self.write_chease_namelist(None,self.CNSE,self.CNTE,self.CNPSIE,self.CNCHIE,B0EXP,IPEXP,8)
					self.write_chease_expeq(None,psin,pprime,jav,B0EXP,IPEXP)
					self.qval_change_chease_namelist('chease_namelist',1.0,qa_target)
					self.run_chease()
			
					self.read_chease_out(CHEout)

					jm1 = self.find_max_ja(0.8,self.psin,self.jav1,1)
					am1 = self.find_max_ja(0.8,self.psin,self.alp1,1)
					am21 = self.find_max_ja(0.8,self.psin,self.alp2,1)
					ball1 = self.find_max_ja(0.8,self.psin,self.ball,0)
					qa1 = self.q[-1]
					qf = interp1d(self.psin,self.q,'cubic')
					q95 = qf(0.95)
					qdel1 = np.ceil(qa1*n) - qa1*n

					print('qdelfix,qdel1,qa0,qa1,n,mm')
					print(self.qdelfix,qdel1,qa0,qa1,n,mm)
				
				misinp2 = equ_save_dir + '/NMISHKA'
				os.chdir(stab_dir)
				copyfile(misinp2,'fort.12')
				os.system('rm '+misinp2)

				if not(self.highq):
					if (qa1 > 5.0):
						self.highq = True

				self.write_mishka_inp('fort.10',n,self.gridn,self.psis,self.xr1,self.sig1,self.xr2,self.sig2,self.ias)

				if (self.highq):
					mn = 41
					if (n > 2):
						mn = 41
					if (n > 5):
						mn = 51
					if (n > 10):
						mn = 71
					if (n > 16):
						mn = 91
				else:
					mn = 31
					if (n > 2):
						mn = 31
					if (n > 5):
						mn = 41
					if (n > 10):
						mn = 51
					if (n > 16):
						mn = 71	

				mis_dir1 = self.mis_dir + str(mn)+'_501'
				print(n,mn,mis_dir1)
				os.system(mis_dir1)
				try:
					status, output = subprocess.getstatusoutput('cat fort.20 | grep INST')
					gr_temp = output.split(':')[1].split()
					gr_nn = -1*int(gr_temp[0])
					gr_chk = int(gr_temp[1])
					if (gr_chk < 20):
						gr_n1 = float(gr_temp[2])
					else:
						gr_n1 = 0.0;
					npsi,psi,dpdpsi,ffp,te,ne,ti,zeff,zimp,alpw = self.read_elite(elitename)
					gr_n1 = np.sqrt(gr_n1)
					gr_n2 = gr_n1 * alpw / float(gr_nn)
				except:
					gr_nn = self.moden[i]
					gr_n1 = 0.0
					gr_n2 = 0.0
				print('growth rate',gr_n1)
				f3.write('%i %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n'%(nw,nh,am0,am20,jm0/1.e6,am1,am21,jm1/1.e6,gr_nn,gr_n1,gr_n2,qdel1,q95,ball0,ball1))
				
			elif(self.run_stab == 'elite'):
			
				os.chdir(stabe_dir)
				copyfile(elitename,'eqin')
				for ii in range(1):
					n = int(self.moden[i])
					
					if (n < 5):
						print('mode n is too low for ELITE, try to increase it')

					if (n > 0):
						inp_name = 'work_' + str(n)	
						inp_name2 = 'work_' + str(n) + '.in'	
						elite_eq = self.elite_dir + '/eliteeq5.7 ' + inp_name
						elite_vac = self.elite_dir + '/elitevac5.7 ' + inp_name
						elite_run = self.elite_dir + '/elite5.7 ' + inp_name
						elite_out = stabe_dir + '/' + inp_name + '.gamma'
						elite_eqout = stabe_dir + '/' + inp_name + '.eqout'
						elite_jedge = stabe_dir + '/' + inp_name + '.jedge'

						self.write_elite_inp(inp_name2,n,self.qdelfix,self.ngride,self.CNCHIE+1,self.psise,self.ndist,self.ias)
						
						status, output = subprocess.getstatusoutput(elite_eq)
						status, output = subprocess.getstatusoutput(elite_vac)
						status, output = subprocess.getstatusoutput(elite_run)

						isqmon = self.read_elite_eq(elite_eqout,elite_jedge)
						if (isqmon):
							qmon = 0.0
						else:
							qmon = -1.0
						
						file = open(elite_out)
						line = file.readline()
						line = file.readline()
						line = file.readline()
						file.close()
		#				print(line)	
		#				npsi,psi,dpdpsi,ffp,te,ne,ti,zeff,zimp,alpw = read_elite(elitename)
						gr_nn = int(line.split()[0])
						gr_n1 = float(line.split()[1])
						gr_n2 = float(line.split()[2])/4.0  #### change later
						qdelf0 = float(line.split()[4])
						print(n,'growth rate',gr_n1)
						f3.write('%i %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n'%(nw,nh,am0,am20,jm0/1.e6,am0,am20,jm0/1.e6,gr_nn,gr_n1,gr_n2,qdelf0,q95,ball0,qmon))
		f3.close()	
		
	def __init__(self):

		self.width = width	
		self.initialise_var()
		#self.read_namelist()	
		
		return
		
	if __name__ == "__main__":

		import pedstab_chease as pedstab

		
		peds = pedstab.pedscan_stab()
		
		sim = ch.chease()
		peds.read_namelist(sim)
		sim.nomap = True
		epslon = sim.epslon
		sim.epslon = 5

			
		if (sim.use_eped2):	sim.beta_iteration_eped2()
		elif (sim.use_eped3):	sim.beta_iteration_eped3()
		else:	sim.beta_iteration()
		sim.nomap = False
		sim.epslon = epslon

		peds.run_stability()
		
