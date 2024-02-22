import os,stat
import sys
import subprocess
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import eqdsk
import warnings
import time
from shutil import copyfile, move, rmtree
import matplotlib.pyplot as plt
from exec_dirs import chease_exec, python3_exec, nubeam_dir

currdir = os.getcwd()
e0 = 1.601 * 1.e-19;
mu0 = 4.0 * np.pi * 1.e-7;

class chease:

	def psin_grid(self,xr1,sig1,xr2,sig2,num):
		# from HELENA
		self.num = num;
		num2 = 1001
		S1 = np.zeros(num2)
		F1 = np.zeros(num2)
		F2 = np.zeros(num2)
		F3 = np.zeros(num2)
		F4 = np.zeros(num2)
		FSUM = np.zeros(num2)
		BGF  = 0.3

		DSI = 1. / float(num2-1)
		S1[0] = 0.
		FSUM[0]  = 0.
		SUM = 0.
		self.fgauss(S1[0],BGF,xr1,xr2,sig1,sig2)
		FINT2 = self.fgaus
		
		for i in range(num2-1):
			S1[i+1] = float(i+1) * DSI
			FINT1 = FINT2
			self.fgauss(S1[i+1],BGF,xr1,xr2,sig1,sig2)
			FINT2 = self.fgaus
			SUM = SUM + (FINT1+FINT2)/2. * DSI
			FSUM[i+1] = SUM

		for i in range(num2-1):
			FSUM[i] = FSUM[i]/FSUM[num2-1]
			
		FSUM[-1] = 1.
		ALFA = 0.
		BETA = 0.
		TYP = 2
		
		psinf = interp1d(FSUM,S1,'cubic')
		self.psin = np.zeros(self.num)
		for i in range(self.num-1):
			FI = float(i)/float(self.num-1)
			self.psin[i] = psinf(FI)*psinf(FI)
	
		self.psin[-1] = 1.0
		self.psin[0] = 0.
			
		return

	def fgauss(self,zs,bgf,xr1,xr2,sig1,sig2):
	
		ZNORM1 = 0.39894 / sig1
		ZNORM2 = 0.39894 / sig2
		
		ZEX1   = -0.5 * (zs - xr1)**2 / sig1**2
		ZEX2   = -0.5 * (zs - xr2)**2 / sig2**2
		DEX1   = -(zs-xr1)/sig1**2
		DEX2   = -(zs-xr2)/sig2**2

		F1     = ZNORM1 * np.exp(ZEX1)
		F2     = ZNORM2 * np.exp(ZEX2)
		DF1    = ZNORM1 * DEX1 * np.exp(ZEX1)
		DF2    = ZNORM2 * DEX2 * np.exp(ZEX2)

		self.fgaus = bgf + (1.0 - bgf) * (F1 + F2) 
	
		return
		
	def read_kinprof_chease(self,filename):
	
		f4 = open(filename,'r')
		self.dat_numk = int(f4.readline().split()[0])
		line = f4.readline().split()
		self.zeff = float(line[0])
		self.zimp = float(line[1])
		self.amain = float(line[2])
		self.aimp = float(line[3])
		
		self.psink = np.zeros(self.dat_numk)
		self.tek = np.zeros(self.dat_numk)
		self.nek = np.zeros(self.dat_numk)
		self.tik = np.zeros(self.dat_numk)
		self.nik = np.zeros(self.dat_numk)
		self.vtk = np.zeros(self.dat_numk)
		self.nI1k= np.zeros(self.dat_numk)
		self.isvt = False

		for i in range(self.dat_numk):
		
			line = f4.readline().split()
			self.psink[i] = float(line[0])	
			self.tek[i] = float(line[1])
			self.nek[i] = float(line[2])
			self.tik[i] = float(line[3])
			self.nik[i] = float(line[4])
			if len(line)>5: self.vtk[i] = float(line[5]); self.isvt = True
			if len(line)>6: self.nI1k[i]= float(line[6]);
		f4.close()
		
		return
		
	def read_kinprof_kprofile_fun(self,filename):
	
		f4 = open(filename,'r')
		
		psind = np.linspace(0,1.0,401)

		
		line_num = 0
		
		while True:
			line = f4.readline()
			if not line: break
			line_num = line_num + 1
			
		f4.close()
		
		rhok = np.zeros(line_num-2)
		psink = np.zeros(line_num-2)
		datk = np.zeros(line_num-2)
		
		f4 = open(filename,'r')
		line = f4.readline()
		line = f4.readline()
		
		for i in range(line_num-2):
			
			line = f4.readline().split()
			rhok[i] = float(line[0])
			psink[i] = float(line[1])
			datk[i] = float(line[2])
		
		f4.close()
	
		if not self.use_rho:	
			datf = interp1d(psink,datk,'cubic')
		else:
			datf = interp1d(rhok,datk,'cubic')
		
		dat = datf(psind)
		
		return (dat)
		
	def read_kinprof_kprofile(self):
		self.isvt = False	
		self.psink = np.linspace(0,1.0,401)
		self.dat_numk = len(self.psink)
		self.tek = self.read_kinprof_kprofile_fun(self.te_file)/1.e3
		self.nek = self.read_kinprof_kprofile_fun(self.ne_file)/1.e1
		self.tik = self.read_kinprof_kprofile_fun(self.ti_file)/1.e3
	
		if (self.ni_file == None):
			self.nik = np.copy(self.nek * (1.0 - (self.zeff - 1.0)/self.zimp))
		else:
			self.nik = self.read_kinprof_kprofile_fun(self.ni_file)/1.e1
		
		self.vtk = np.zeros(self.dat_numk)
		if not (self.vt_file == None):
			self.vtk = self.read_kinprof_kprofile_fun(self.vt_file)
			self.isvt = True
		self.nI1k= np.zeros(self.dat_numk)
		return

	def write_kinprof_kprofile(self,prof,name,filename):

		proff = interp1d(self.rho,prof,'cubic')
		psif = interp1d(self.rho,self.psin,'cubic')
		Rf = interp1d(self.psi_map,self.R_map,'cubic')
		
		f4 = open(filename,'w')
		f4.write('R(m)       Z(m)       %s\n'%name)
		f4.write('Fitted Values: X_Norm_Rho, Psi_Norm, Y\n')

		for i in range(101):

			rho = float(i) / 100.0
			val = proff(rho)
			psi = psif(rho)
			f4.write('%9.6f\t%9.6f\t%9.6f\n'%(rho,psi,val))

		f4.close()
		return
		
	def read_kinprof_eped_fun(self,filename):
	
		self.eped_prof = np.zeros(shape=(3,8))

		f4 = open(filename,'r')
		
		for i in range(3):
			line = f4.readline().split()
			
			for j in range(8):
				self.eped_prof[i,j] = float(line[j])
				
		f4.close()

		return
	
	def read_kinprof_eped(self):
	
		self.psink = np.copy(self.psin)
		self.dat_numk = len(self.psink)
		self.tek = np.zeros(self.dat_numk)
		self.nek = np.zeros(self.dat_numk)
		self.tik = np.zeros(self.dat_numk)
		self.nik = np.zeros(self.dat_numk)
		self.vtk = np.zeros(self.dat_numk)
		self.nI1k= np.zeros(self.dat_numk)
		
		if(self.load_eped_file):
			self.read_kinprof_eped_fun(self.eped_file)
		else:
			pass
		
		epedf =np.copy(self.eped_prof)
		
		for i in range(self.num):
			psin = self.psink[i]
			self.nek[i] = epedf[0,1] + epedf[0,0]*(np.tanh(2.0*(1.0-epedf[0,3])/epedf[0,4]) - np.tanh(2.0*(psin-epedf[0,3])/epedf[0,4]))
			if (psin < epedf[0,5]):
				self.nek[i] = self.nek[i] + epedf[0,2] * ((1.0 - (psin/epedf[0,5]) ** epedf[0,6]) ** epedf[0,7])
				
			self.tek[i] = epedf[1,1] + epedf[1,0]*(np.tanh(2.0*(1.0-epedf[1,3])/epedf[1,4]) - np.tanh(2.0*(psin-epedf[1,3])/epedf[1,4]))
			if (psin < epedf[1,5]):
				self.tek[i] = self.tek[i] + epedf[1,2] * ((1.0 - (psin/epedf[1,5]) ** epedf[1,6]) ** epedf[1,7])
	
			self.tik[i] = epedf[2,1] + epedf[2,0]*(np.tanh(2.0*(1.0-epedf[2,3])/epedf[2,4]) - np.tanh(2.0*(psin-epedf[2,3])/epedf[2,4]))
			if (psin < epedf[2,5]):
				self.tik[i] = self.tik[i] + epedf[2,2] * ((1.0 - (psin/epedf[2,5]) ** epedf[2,6]) ** epedf[2,7])
			
			self.nik[i] = self.nek[i] * (1. - (self.zeff-1.0)/self.zimp)
		
		##Regrid with updated pedestal width	
		if not (self.ped_width == epedf[1,4]):
			print('>>> Re-grid for new pedestal width, %3.2f'%epedf[1,4])
			self.ped_width = epedf[1,4]
			self.grid_xr2 = 1.0 - 0.5*self.ped_width
			self.grid_sig2 = self.ped_width
			self.psin_grid(self.grid_xr1,self.grid_sig1,self.grid_xr2,self.grid_sig2,self.grid_n)
		return

	def construct_zeff(self):

		self.zeffk = np.zeros(self.dat_numk)
		if self.nI1k[0] == 0: self.zeffk = np.ones(self.dat_numk) * self.zeff;
		else:
			self.zeffk = 1. + (self.zimp-1.)*(self.zimp)*self.nI1k/self.nek

		#reconstruct nik (main+imp)
		self.nik = (1.-(self.zeffk-1.)/self.zimp) * self.nek;
		return

	def scale_density(self,not_only_print=True,print_density=True):

		len2 = len(self.rc)
		
		RR = np.zeros(2*len2-1)
		NN = np.copy(RR)

		for i in range(len2-1):
			RR[i] = self.rc[len2-1-i] - self.r[len2-1-i]
			RR[2*len2-2-i] = self.rc[len2-1-i] + self.r[len2-1-i]
			NN[i] = self.ne[len2-1-i]
			NN[2*len2-2-i] = self.ne[len2-1-i]

		RR[len2-1] = self.rc[0]
		NN[len2-1] = self.ne[0]

		line_sum = np.trapz(NN,x=RR)
		line_avg = line_sum /(RR[-1]-RR[0]);

		if not (not_only_print):
			self.scaled_density = 1.0
			if print_density:
				print('Line average density(horizontal) = %f 10(19)/m3 '%(line_avg))
			return
	
		scaled_density = self.target_density / line_avg

		if not (self.scaled_density == scaled_density):
			if print_density: #################################################
				print('Line average density(horizontal) = %f >>> %f 10(19)/m3 '%(line_avg,self.target_density))
			self.scaled_density = scaled_density

		return
				
	def write_kinprof(self):
	
		if (self.kinprof_type == 1):
			f = open('chease_kinprof_new','w')

			f.write('%i \n'%self.dat_numk)
			f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.zeff,self.zimp,self.amain,self.aimp))
			for i in range(self.dat_numk):
				line = '%9.6f\t'%(self.psink[i]);
				line+= '%9.6f\t'%(self.tek[i]);
				line+= '%9.6f\t'%(self.nek[i]*self.scaled_density);
				line+= '%9.6f\t'%(self.tik[i]);
				line+= '%9.6f\t'%(self.nik[i]*self.scaled_density);
				line+= '%9.6f\t'%(self.vtk[i]);
				line+= '%9.6f\n'%(self.nI1k[i]);
				f.write(line)

			f.close()

		elif (self.kinprof_type == 2):
			self.write_kinprof_kprofile(self.ne * 1.e1,'Ne [10(18)/m3]','NE.dat_new')
			self.write_kinprof_kprofile(self.te * 1.e3,'Te [eV]',       'TE.dat_new')
			self.write_kinprof_kprofile(self.ti * 1.e3,'Ti [eV]',       'TI.dat_new')
			self.write_kinprof_kprofile(self.vt * 1.e0,'VT [km/s]',     'VT.dat_new')
			self.write_kinprof_kprofile(self.nI1* 1.e1,'NI1[10(18)/m3]','NI1.dat_new')
				

		return
		
	def read_ext_prof(self,filename,isp=False):

		f4 = open(filename,'r')
		line = f4.readline()
		dat_num = int(line.split()[0])
		
		dat1 = np.zeros(dat_num)
		dat2 = np.zeros(dat_num)
		dat3 = np.zeros(dat_num)

		if not(isp):
			for i in range(dat_num):
				line = f4.readline().split()

				dat1[i] = float(line[0])
				dat2[i] = float(line[1])
		else:
			for i in range(dat_num):
				line = f4.readline().split()

				dat1[i] = float(line[0])
				dat2[i] = float(line[1])
				dat3[i] = float(line[2])

		rhop = interp1d(self.rho_map,self.psi_map,'cubic')
		psit = rhop(dat1)

		if (isp):
			return (psit,dat2,dat3)
		else:
			return (psit,dat2)
		
	def read_hagermap(self,filename):
	
		f4 = open(filename,'r')
		line = f4.readline()	
		self.dat_numh = int(f4.readline().split()[0])
		
		self.psinh = np.zeros(self.dat_numh)
		self.ravh = np.zeros(self.dat_numh)
		self.b02avh = np.zeros(self.dat_numh)
		self.bgpbh = np.zeros(self.dat_numh)
		self.fth = np.zeros(self.dat_numh)
		self.rch = np.zeros(self.dat_numh)
		self.rh = np.zeros(self.dat_numh)
		self.drch = np.zeros(self.dat_numh)
		self.qh = np.zeros(self.dat_numh)
		self.fh = np.zeros(self.dat_numh)
		self.bmaxh = np.zeros(self.dat_numh)
		self.bminh = np.zeros(self.dat_numh)
		self.gm1h = np.zeros(self.dat_numh)
		self.gm9h = np.zeros(self.dat_numh)
		self.ah = np.zeros(self.dat_numh)
		self.vh = np.zeros(self.dat_numh)
		
		for i in range(self.dat_numh):
		
			line = f4.readline().split()
			
			self.psinh[i] = float(line[0])
			self.ravh[i] = float(line[1])
			self.b02avh[i] = float(line[2])
			self.bgpbh[i] = float(line[3])
			self.fth[i] = float(line[4])
			self.rch[i] = float(line[5])
			self.rh[i] = float(line[6])
			self.drch[i] = float(line[7])
			self.qh[i] = float(line[8])
			self.fh[i] = float(line[9])
			self.bmaxh[i] = float(line[10])
			self.bminh[i] = float(line[11])
			self.gm1h[i] = float(line[12])
			self.ah[i] = float(line[13])
			self.vh[i] = float(line[14])
			self.gm9h[i] = float(line[15])
			
		
		line = f4.readline().split()
		
		self.psia = abs(float(line[0]))
		
		f4.close()
		
		return
		
	def interpol(self,print_density=True):

		nef = interp1d(self.psink,self.nek,  'cubic')
		tef = interp1d(self.psink,self.tek,  'cubic')
		nif = interp1d(self.psink,self.nik,  'cubic')
		tif = interp1d(self.psink,self.tik,  'cubic')
		vtf = interp1d(self.psink,self.vtk,  'cubic')
		nI1f= interp1d(self.psink,self.nI1k, 'cubic')
		zef = interp1d(self.psink,self.zeffk,'cubic')

		nef2 = interp1d(self.psink[0:3],self.nek[0:3], 'quadratic')
		tef2 = interp1d(self.psink[0:3],self.tek[0:3], 'quadratic')
		nif2 = interp1d(self.psink[0:3],self.nik[0:3], 'quadratic')
		tif2 = interp1d(self.psink[0:3],self.tik[0:3], 'quadratic')
		vtf2 = interp1d(self.psink[0:3],self.vtk[0:3], 'quadratic')
		nI12 = interp1d(self.psink[0:3],self.nI1k[0:3],'quadratic')
		
		if (self.use_ext_pressure):
			pres_exf = interp1d(self.psi_ex1,self.pres_ext,'cubic')
			pres_expf = interp1d(self.psi_ex1,self.pres_extp,'cubic')
		if (self.use_ext_current):
			curr_exf = interp1d(self.psi_ex2,self.curr_ext,'cubic')
		
		if not (self.eped_first_run):

			ravf = interp1d(self.psinh,self.ravh,'cubic')
			b02avf = interp1d(self.psinh,self.b02avh,'cubic')
			bgpbf = interp1d(self.psinh,self.bgpbh,'cubic')
			ftf = interp1d(self.psinh,self.fth,'cubic')
			rcf = interp1d(self.psinh,self.rch,'cubic')
			af = interp1d(self.psinh,self.ah,'cubic')
			vf = interp1d(self.psinh,self.vh,'cubic')
			qf = interp1d(self.psinh,self.qh,'cubic')
			rf = interp1d(self.psinh,self.rh,'cubic')
			ff = interp1d(self.psinh,self.fh,'cubic')
			bmaxf = interp1d(self.psinh,self.bmaxh,'cubic')
			bminf = interp1d(self.psinh,self.bminh,'cubic')
			gm1f = interp1d(self.psinh,self.gm1h,'cubic')
			drcf = interp1d(self.psinh,self.drch,'cubic')
			gm9f = interp1d(self.psinh,self.gm9h,'cubic')

		else:
			self.psinh = np.linspace(0,1,20)
			self.psinh2 = np.linspace(1,2,20)
			self.psia = 1.0
			ravf = interp1d(self.psinh,self.psinh2,'cubic')
			b02avf =  interp1d(self.psinh,self.psinh2,'cubic')
			bgpbf = interp1d(self.psinh,self.psinh2,'cubic')
			ftf =  interp1d(self.psinh,self.psinh2,'cubic')
			rcf = interp1d(self.psinh,self.psinh2,'cubic')	
			af = interp1d(self.psinh,self.psinh2,'cubic')
			vf = interp1d(self.psinh,self.psinh2,'cubic')
			qf = interp1d(self.psinh,self.psinh2,'cubic')
			rf = interp1d(self.psinh,self.psinh2,'cubic')
			ff = interp1d(self.psinh,self.psinh2,'cubic')
			bmaxf = interp1d(self.psinh,self.psinh2,'cubic')
			bminf = interp1d(self.psinh,self.psinh2,'cubic')
			drcf = interp1d(self.psinh,self.psinh2,'cubic')
			gm1f = interp1d(self.psinh,self.psinh2,'cubic')
			gm9f = interp1d(self.psinh,self.psinh2,'cubic')
	
		self.ne = np.zeros(self.num)
		self.te = np.zeros(self.num)
		self.ni = np.zeros(self.num)
		self.ti = np.zeros(self.num)
		self.vt = np.zeros(self.num)
		self.nI1= np.zeros(self.num)

		self.zeffs = np.zeros(self.num)

		self.rav = np.zeros(self.num)
		self.b02av = np.zeros(self.num)
		self.bgpb = np.zeros(self.num)
		self.ft = np.zeros(self.num)
		self.rc = np.zeros(self.num)
		self.q = np.zeros(self.num)
		self.f = np.zeros(self.num)
		self.r = np.zeros(self.num)
		self.pe = np.zeros(self.num)
		self.pt = np.zeros(self.num)
		self.bmax = np.zeros(self.num)
		self.bmin = np.zeros(self.num)
		self.gm1 = np.zeros(self.num)
		self.gm9 = np.zeros(self.num)
		self.area = np.zeros(self.num)
		self.vol = np.zeros(self.num)
		
		self.pres_ex = np.zeros(self.num)
		self.pres_exp = np.zeros(self.num)
		self.curr_ex = np.zeros(self.num)
		
		for i in range(self.num):
			psinn = self.psin[i]
		
			if (psinn > self.psink[1]):	
				self.ne[i] = nef(psinn)
				self.te[i] = tef(psinn)
				self.ni[i] = nif(psinn)
				self.ti[i] = tif(psinn)
			else:
				self.ne[i] = nef2(psinn)
				self.te[i] = tef2(psinn)
				self.ni[i] = nif2(psinn)
				self.ti[i] = tif2(psinn)

			self.vt[i] = vtf(psinn)
			self.nI1[i]= nI1f(psinn)
			self.zeffs[i]= zef(psinn)
			self.rav[i] = ravf(psinn)
			self.b02av[i] = b02avf(psinn)
			self.bgpb[i] = bgpbf(psinn)
			self.ft[i] = ftf(psinn)
			self.rc[i] = rcf(psinn)
			self.q[i] = qf(psinn)
			self.r[i] = rf(psinn)
			self.f[i] = abs(ff(psinn)) ## USE abs |F|
			self.bmax[i] = bmaxf(psinn)
			self.bmin[i] = bminf(psinn)
			self.gm1[i] = gm1f(psinn)
			self.gm9[i] = gm9f(psinn)
			self.area[i] = af(psinn)
			self.vol[i] = vf(psinn)
			
			self.pe[i] = self.ne[i] * self.te[i] * 1.e19 * 1.e3 * e0
			self.pt[i] = (self.ne[i]*self.te[i] + self.ni[i]*self.ti[i]) * 1.e19 * 1.e3 * e0
	
		self.BMAG = self.f[0] / self.rc[0]

		if (self.use_ext_pressure):
			for i in range(self.num):
				psinn = self.psin[i]
				self.pres_ex[i] = pres_exf(psinn)
				self.pres_exp[i] = pres_expf(psinn)

		if (self.use_ext_current):
			for i in range(self.num):
				psinn = self.psin[i]
				self.curr_ex[i] = curr_exf(psinn) #* abs(self.bcent) / self.BMAG

		self.scale_density(self.use_scaled_density,print_density)
		self.ne = self.ne * self.scaled_density
		self.ni = self.ni * self.scaled_density
		self.pe = self.pe * self.scaled_density
		self.pt = self.pt * self.scaled_density
				
		self.dne = np.zeros(self.num)
		self.dte = np.zeros(self.num)
		self.dni = np.zeros(self.num)
		self.dti = np.zeros(self.num)
		self.drc = np.zeros(self.num)
		self.dq = np.zeros(self.num)
		self.dpe = np.zeros(self.num)
		self.dpt = np.zeros(self.num)
		self.dp_ex = np.zeros(self.num)
		
		deps = 1.e-5
		for i in range(self.num-2):
		
			psinn = self.psin[i+1]
			if (psinn > self.psink[1]):		
				self.dne[i+1] = (nef(psinn+deps)-nef(psinn-deps))/2.0/deps/self.psia
				self.dte[i+1] = (tef(psinn+deps)-tef(psinn-deps))/2.0/deps/self.psia
				self.dni[i+1] = (nif(psinn+deps)-nif(psinn-deps))/2.0/deps/self.psia
				self.dti[i+1] = (tif(psinn+deps)-tif(psinn-deps))/2.0/deps/self.psia
			else:
				self.dne[i+1] = (nef2(psinn+deps)-nef2(psinn-deps))/2.0/deps/self.psia
				self.dte[i+1] = (tef2(psinn+deps)-tef2(psinn-deps))/2.0/deps/self.psia
				self.dni[i+1] = (nif2(psinn+deps)-nif2(psinn-deps))/2.0/deps/self.psia
				self.dti[i+1] = (tif2(psinn+deps)-tif2(psinn-deps))/2.0/deps/self.psia
			
			self.dq[i+1] = (qf(psinn+deps)-qf(psinn-deps))/2.0/deps/self.psia

		for i in range(self.num-2):

			if (self.dne[i+1]*self.psia > 0.):	self.dne[i+1] = 0.
			if (self.dte[i+1]*self.psia > 0.):	self.dte[i+1] = 0.
			if (self.dni[i+1]*self.psia > 0.):	self.dni[i+1] = 0.
			if (self.dti[i+1]*self.psia > 0.):	self.dti[i+1] = 0.

		self.dne[0]	 = 2.0*self.dne[1] - self.dne[2]
		self.dte[0]	 = 2.0*self.dte[1] - self.dte[2]
		self.dni[0]	 = 2.0*self.dni[1] - self.dni[2]
		self.dti[0]	 = 2.0*self.dti[1] - self.dti[2]
		self.drc[0]	 = 2.0*self.drc[1] - self.drc[2]
		self.dq[0]	 = 2.0*self.dq[1] - self.dq[2]
		#self.dpe[0]	 = 2.0*self.dpe[1] - self.dpe[2]
		#self.dpt[0]	 = 2.0*self.dpt[1] - self.dpt[2]
		
		
		self.dne[self.num-1] = 2.0*self.dne[self.num-2] - self.dne[self.num-3]
		self.dte[self.num-1] = 2.0*self.dte[self.num-2] - self.dte[self.num-3]
		self.dni[self.num-1] = 2.0*self.dni[self.num-2] - self.dni[self.num-3]
		self.dti[self.num-1] = 2.0*self.dti[self.num-2] - self.dti[self.num-3]
		self.drc[self.num-1] = 2.0*self.drc[self.num-2] - self.drc[self.num-3]
		self.dq[self.num-1] = 2.0*self.dq[self.num-2] - self.dq[self.num-3]
		#self.dpe[self.num-1] = 2.0*self.dpe[self.num-2] - self.dpe[self.num-3]
		#self.dpt[self.num-1] = 2.0*self.dpt[self.num-2] - self.dpt[self.num-3]

		if self.dne[0] * self.psia > 0.: self.dne[0] = 0.
		if self.dte[0] * self.psia > 0.: self.dte[0] = 0.
		if self.dni[0] * self.psia > 0.: self.dni[0] = 0.
		if self.dti[0] * self.psia > 0.: self.dti[0] = 0.

		if self.dne[-1] * self.psia > 0.: self.dne[-1] = 0.
		if self.dte[-1] * self.psia > 0.: self.dte[-1] = 0.
		if self.dni[-1] * self.psia > 0.: self.dni[-1] = 0.
		if self.dti[-1] * self.psia > 0.: self.dti[-1] = 0.
		
		self.drc = drcf(self.psin)
		self.drc[0] = 0.0

		for i in range(self.num):
			
			self.dpe[i] = (self.dne[i]*self.te[i]+self.ne[i]*self.dte[i]) * 1.e19 * 1.e3 * e0
			self.dpt[i] = (self.dni[i]*self.ti[i]+self.ni[i]*self.dti[i]) * 1.e19 * 1.e3 * e0 + self.dpe[i]
			
		if (self.use_ext_pressure):	
			for i in range(self.num-2):
				psinn = self.psin[i+1]
				self.dp_ex[i+1] = (pres_exf(psinn+deps) - pres_exf(psinn-deps))/2.0/deps/self.psia

		#		if self.dp_ex[i+1]*self.psia > 0.:	self.dp_ex[i+1] = 0.

		self.dp_ex[0] 	 = 2.0*self.dp_ex[1] - self.dp_ex[2]
		self.dp_ex[self.num-1] = 2.0*self.dp_ex[self.num-2] - self.dp_ex[self.num-3]

		if self.dp_ex[0] * self.psia > 0.: self.dp_ex[0] = 0.

		return

	def load_kin_profile(self,print_density=True):

		if (self.kinprof_type == 1):
			self.read_kinprof_chease(self.chkin_file)
		elif (self.kinprof_type == 2):
			self.read_kinprof_kprofile()
		elif (self.kinprof_type == 3):
			self.read_kinprof_eped()

		if (self.pedscan_ti):
			psink = np.linspace(0,1.0,401)
			tik = self.read_kinprof_kprofile_fun(self.ti_file)/1.e3
			tikf = interp1d(psink,tik,'cubic')
			self.tik = tikf(self.psink)

		if not (self.vt_file == None):
			psink = np.linspace(0,1.0,401)
			vtk = self.read_kinprof_kprofile_fun(self.vt_file)
			vtkf = interp1d(psink,vtk,'cubic')
			self.vtk = vtkf(self.psink)
		else: self.vtk = np.zeros(self.dat_numk)

		if (self.use_rho and not self.eped_first_run):
			rhok = np.copy(self.psink)
			psif = interp1d(self.rho_map,self.psi_map,'cubic')
			psink = psif(rhok)
			self.psink = np.copy(psink)
			
		if (self.adjust_prof and not self.eped_first_run):
		
			self.nek = self.adjust_kinprof(self.R_mapi,self.psi_mapi,self.R_map,self.psi_map,self.psink,self.nek)
			self.tek = self.adjust_kinprof(self.R_mapi,self.psi_mapi,self.R_map,self.psi_map,self.psink,self.tek)
			self.nik = self.adjust_kinprof(self.R_mapi,self.psi_mapi,self.R_map,self.psi_map,self.psink,self.nik)
			self.tik = self.adjust_kinprof(self.R_mapi,self.psi_mapi,self.R_map,self.psi_map,self.psink,self.tik)
			self.vtk = self.adjust_kinprof(self.R_mapi,self.psi_mapi,self.R_map,self.psi_map,self.psink,self.vtk)
			self.nI1k= self.adjust_kinprof(self.R_mapi,self.psi_mapi,self.R_map,self.psi_map,self.psink,self.nI1k)
				
		if not(self.eped_first_run): self.read_hagermap(self.hager_mapf)	
			
		if (self.use_ext_pressure):
			self.psi_ex1, self.pres_extl, self.pres_extp = self.read_ext_prof(self.pres_file,True)
			self.pres_ext = np.copy(2.0 / 3.0 * (self.pres_extp + self.pres_extl*0.5))
			
		if (self.use_ext_current):
			self.psi_ex2, self.curr_ext = self.read_ext_prof(self.curr_file)

		if not self.vtk[0]==0: self.isvt = True
		else: self.isvt = False
	
		self.construct_zeff()
		self.interpol(print_density)
		
		return
				
	def hager_bs(self):

		f = open('neo_coefs','w')
		f.write('%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%13s\n'%('psi_norm','gamma','nu_i','nu_e','gamma_cs','b02av','sigma'))
		f.write('%i\n'%self.num)	
		f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%13.6e\n'%(0.,0.,0.,0.,0.,self.b02av[0],0.))

		self.hager_const()
		self.alps = np.zeros(self.num)

		if (self.use_neo):# Radl POP 2021
			self.use_hager = False
			self.use_chang = False
		
		if (self.use_hager):
			self.use_chang = False
			self.use_neo   = False

		self.jb_hag = np.zeros(self.num)
		self.jb_hag1 = np.zeros(self.num)
		self.jb_hag2 = np.zeros(self.num)
		
		self.betacol = np.zeros(self.num)
		self.betagradti = np.zeros(self.num)
		self.cb_hag = np.zeros(self.num)
		self.delpsi = np.zeros(self.num)
		self.nuis = np.zeros(self.num)
		self.nues = np.zeros(self.num)
		self.ZLE = np.zeros(self.num)
	
		for i in range(self.num-1):
			i = i + 1
			eps_hag = self.r[i]/self.rc[i]
			if (eps_hag <= 0.0):
				eps_hag = 1.e-8
			
			Z1_hag = self.ne[i]/self.ni[i]
			Z2_hag = self.zeffs[i]
			Z_hag = np.power((self.zmain**2.0)*Z1_hag*Z2_hag,0.25)
			self.Z = Z2_hag
			
			ne_hag = self.ne[i] * 1.e19
			ni_hag = self.ni[i] * 1.e19
			te_hag = self.te[i] * 1.e3
			ti_hag = self.ti[i] * 1.e3
			rpe_hag = ne_hag*te_hag / (ni_hag*ti_hag+ne_hag*te_hag)
			
			vti_hag = np.sqrt(2.0*ti_hag*e0/self.mi_hag)
			DPSI_hag = self.mi_hag * vti_hag * self.f[i] / self.zmain / e0 / self.bmin[i] * np.sqrt(1.0 - (self.bmin[i]/self.bmax[i]))

			psid_hag = abs(self.psia) * (1-self.psin[i])
			
			self.delpsi[i] = DPSI_hag
			
			zle_hag = 31.3 - np.log(np.sqrt(ne_hag)/te_hag)

			if (self.use_hager):		
				zli_hag = 30.0 - np.log((Z1_hag**3.0)*np.sqrt(ni_hag)/(ti_hag**1.5))
			else:
				zli_hag = 30.0 - np.log((Z_hag **3.0)*np.sqrt(ni_hag)/(ti_hag**1.5))
			
			self.ZLE[i] = zle_hag
			
			nues_hag = 6.921*1.e-18 * self.q[i] * self.rav[i] * ne_hag * Z2_hag * zle_hag / (te_hag**2)/(eps_hag**1.5)
			nuis_hag = 4.900*1.e-18 * self.q[i] * self.rav[i] * ni_hag * (Z_hag**4.0) * zli_hag / (ti_hag**2) / (eps_hag**1.5)
			self.nues[i] = nues_hag
			self.nuis[i] = nuis_hag
	
			if (nues_hag >= 1.e8):
				nues_hag = 1.e8
			if (nuis_hag >= 1.e8):
				nuis_hag = 1.e8
				
			if (self.use_chang):
			
				eps1_hag=(eps_hag-0.44) + 0j
		
				beta1=(eps1_hag)**(0.7)
				beta=beta1.real
				if(Z2_hag > 5.0):
					alpha1 = 0
				else:
					alpha1 =(-(Z2_hag*Z2_hag)+5.998*Z2_hag-4.981)/(4.294*Z2_hag*Z2_hag-14.07* Z2_hag+12.61)
				 
				a1 = np.tanh(3.2*beta*((eps_hag**1.5) * nues_hag)**1.4/ (Z2_hag**alpha1))
				a2 = np.tanh(2.2*beta*(eps_hag**2.8) * (nues_hag**0.1)/(Z2_hag**alpha1))
				delta = 1.0+0.55*Z2_hag**0.2*(a1+(1-np.exp(-nues_hag/0.1))*a2)
				deltapsi = self.rav[i] * np.sqrt(eps_hag)*np.sqrt(2*te_hag*self.me_hag/e0)

				H= (0.6/Z2_hag**4)*np.exp(-abs(psid_hag/(3.3*np.log(eps_hag**1.5* nues_hag + 2. )*deltapsi)))

				self.ft[i] = self.ft[i] * (1-H)
								
				
			c1_hag = 0.5 * ((self.a1_hag-self.b1_hag)*np.tanh(2.0*(eps_hag-self.epsc1_hag)/self.w1_hag) + (self.a1_hag+self.b1_hag)) - 0.5 * ((self.a2_hag-self.b2_hag)*np.tanh(2.0*(eps_hag-self.epsc2_hag)/self.w2_hag) + (self.a2_hag-self.b2_hag))
        
			self.betacol[i] = (1.0+c1_hag*self.c2_hag*(nues_hag**2)) / (1.0+c1_hag*(nues_hag**2))
        
			self.betagradti[i] = -1.0 * self.lam1_hag * nuis_hag / (1.+self.lam2_hag*(nuis_hag**2) +self.lam3_hag*(nuis_hag**4)) * ((1.-eps_hag)**self.lam4_hag) * (abs(self.drc[i])**self.lam5_hag) * DPSI_hag / self.ti[i] * abs(self.dti[i])

			self.cb_hag[i] = 3.0 * self.bgpb[i] / self.b02av[i]
			
			ce1_hag = (10./12.)*(289.914+408.*(Z2_hag**2)) / (178.+425.678 * Z2_hag +192.*(Z2_hag**2))
			ce2_hag = (10./12.)*(468.105+912.*Z2_hag)      / (178.+425.678 * Z2_hag +192.*(Z2_hag**2))
			ce3_hag = (10./12.)*(1477.85+2592.*Z2_hag)     / (178.+425.678 * Z2_hag +192 *(Z2_hag**2))
        
			l31h_hag = self.cb_hag[i] * ((Z2_hag*self.rc[-1]*self.q[i])**2) * ((4.*np.sqrt(2.) + 13.*Z2_hag)*ce1_hag + 6.*Z2_hag*ce2_hag) / 4. / (np.sqrt(2.) + Z2_hag) / (eps_hag**3)/(nues_hag**2)
			l32h_hag = self.cb_hag[i] * ((Z2_hag*self.rc[-1]*self.q[i])**2) * ((4.*np.sqrt(2.) + 13.*Z2_hag)*ce2_hag + 6.*Z2_hag*ce3_hag) / 4. / (np.sqrt(2.) + Z2_hag) / (eps_hag**3)/(nues_hag**2)
			l34h_hag = l31h_hag

			l31_hag = 2.*self.ft[i]*Z2_hag*(2.4+Z2_hag)/(1.+Z2_hag)/(1.-0.99*self.ft[i])/nues_hag	
			l32_hag = self.ft[i]*((0.05+0.62*Z2_hag)/(0.18 - 0.0666*self.ft[i])* np.sqrt(Z2_hag) - (0.56+1.93*Z2_hag)/(0.85-0.3145*self.ft[i])/ (1.+Z2_hag)) / (1.+0.44*Z2_hag)/Z2_hag/nues_hag
			l34_hag = self.ft[i] * Z2_hag * (2.4 + Z2_hag)/(1.+Z2_hag)/ (0.5 - 0.25*self.ft[i])/nues_hag

			a_hag = 32.*(ce1_hag**2) + 8.*np.sqrt(2.)*(13.*(ce1_hag**2) + 9.*ce1_hag*ce2_hag + 2.*(ce2_hag**2))*Z2_hag + (((13.*ce1_hag+9.*ce2_hag)**2) + 7.*(ce2_hag**2) + 6.*ce3_hag* (6.*ce1_hag + 4.*ce2_hag))*(Z2_hag**2)
			b_hag = 6.*ce2_hag*Z2_hag*(np.sqrt(2.)+Z2_hag) + ce1_hag* (8.+17.*np.sqrt(2.)*Z2_hag+13.*(Z2_hag**2))
        
			nuesc_hag = 5.*np.sqrt(self.cb_hag[i])*self.q[i]*self.rc[-1]*Z2_hag / (eps_hag**1.5) * np.sqrt(a_hag/b_hag)
			
			c5_hag = l31h_hag/l31_hag/nuesc_hag
			
			gam1_hag = (1.+c5_hag*nues_hag) /(1. + (c5_hag*nues_hag*l31_hag/l31h_hag))
			gam2_hag = (1.+c5_hag*nues_hag) /(1. + (c5_hag*nues_hag*l32_hag/l32h_hag))
			gam4_hag = (1.+c5_hag*nues_hag) /(1. + (c5_hag*nues_hag*l34_hag/l34h_hag))
			gama_hag = (1.+c5_hag*(nuis_hag**2)) /(1. + c5_hag*(nuis_hag**2)*self.alp_hag /self.alph_hag)
		
			f31_hag = self.ft[i] / (1. + ((1.-0.1*self.ft[i])*np.sqrt(nues_hag)) + (0.5*(1.-1.00*self.ft[i])*nues_hag/Z2_hag))
		
			if (self.use_hager):
				f31_hag = self.ft[i] / (1. + ((1.-0.1*self.ft[i])*np.sqrt(nues_hag)) + (0.5*(1.-0.99*self.ft[i])*nues_hag/Z2_hag))
			elif (self.use_chang):
				f31_hag = f31_hag * delta
				
			
			f32_hag = self.ft[i] / (1. + 0.26*(1.-self.ft[i])*np.sqrt(nues_hag) + 0.18*(1.-0.37*self.ft[i])*nues_hag/np.sqrt(Z2_hag))
			f33_hag = self.ft[i] / (1. + (1.+0.6*self.ft[i])*np.sqrt(nues_hag) + 0.85*(1.-0.37*self.ft[i])*nues_hag*(1.+Z1_hag))
			f34_hag = self.ft[i] / (1. + (1.-0.1*self.ft[i])*np.sqrt(nues_hag) + 0.5*(1.-0.5*self.ft[i])*nues_hag/Z2_hag)
			alp0_hag = -1.17*(1.-self.ft[i])/(1.-0.22*self.ft[i]-0.19*(self.ft[i]**2))
			alps_hag = ((alp0_hag + 0.25*(1-(self.ft[i]**2))*np.sqrt(nuis_hag))/ (1.+0.5*np.sqrt(nuis_hag)) + 0.315*(nuis_hag**2)*(self.ft[i]**6)) / (1.+0.15*(nuis_hag**2)*(self.ft[i]**6))
			
			if (self.use_chang):
				f32_hag = f32_hag * delta
				f33_hag = f33_hag * delta
				f34_hag = f34_hag * delta

			if (self.use_neo):
				f31_hag = 1.+(0.67*(1-0.70*self.ft[i])*np.sqrt(nues_hag)/(0.56+0.44*Z2_hag)) + (0.52+0.086*np.sqrt(nues_hag))*(1+0.87*self.ft[i])*nues_hag/(1.+1.13*np.sqrt(Z2_hag-1.))
				f31_hag = self.ft[i]/f31_hag
				f32_hag = 1.+(0.23*(1-0.96*self.ft[i])*np.sqrt(nues_hag)/np.sqrt(Z2_hag))
				f32_hag = f32_hag + 0.13*(1-0.38*self.ft[i])*nues_hag/(Z2_hag**2) * (np.sqrt(1.+2.*np.sqrt(Z2_hag-1))+(self.ft[i]**2)*np.sqrt((0.075+0.25*(Z2_hag-1.)**2)*nues_hag))
				f32_hag = self.ft[i]/f32_hag
				
				f33_hag = 1.+0.87*(1+0.39*self.ft[i])*np.sqrt(nues_hag)/(1.+2.95*(Z2_hag-1.)**2)+1.53*(1-0.37*self.ft[i])*nues_hag*(2.+0.375*(Z2_hag-1.))
				f33_hag = self.ft[i]/f33_hag

				alp0_hag= -(0.62+0.055*(Z2_hag-1.))/(0.53+0.17*(Z2_hag-1.))*(1-self.ft[i])/(1.-(0.31-0.065*(Z2_hag-1.))*self.ft[i]-0.25*self.ft[i]**2.)
				alps_hag = ((alp0_hag+0.7*Z2_hag*np.sqrt(self.ft[i])*np.sqrt(nuis_hag))/(1+0.18*np.sqrt(nuis_hag))-0.002*nuis_hag**2.*self.ft[i]**6.)/(1.+0.004**nuis_hag**2.*self.ft[i]**6.)

			l31s_hag = (1.+1.4/(Z2_hag+1.))*f31_hag - 1.9/(Z2_hag+1.)* (f31_hag**2) + 0.3/(Z2_hag+1.)*(f31_hag**3) + 0.2/(Z2_hag+1.)*(f31_hag**4)
			
			if (self.use_chang):
				l32s_hag = (0.05+0.61*Z2_hag)/Z2_hag/(1.+0.44*Z2_hag)*(f32_hag - (f32_hag**4)) + 1./(1.+0.22*Z2_hag)*((f32_hag**2)-(f32_hag**4) -1.2*(f32_hag**3)+1.2*(f32_hag**4)) + 1.2/(1.+0.5*Z2_hag)* (f32_hag**4)
			else:
				l32s_hag = (0.05+0.62*Z2_hag)/Z2_hag/(1.+0.44*Z2_hag)*(f32_hag - (f32_hag**4)) + 1./(1.+0.22*Z2_hag)*((f32_hag**2)-(f32_hag**4) -1.2*(f32_hag**3)+1.2*(f32_hag**4)) + 1.2/(1.+0.5*Z2_hag)* (f32_hag**4)
				
				
			l32s_hag = l32s_hag - (0.56+1.93*Z2_hag)/Z2_hag/(1.+0.44*Z2_hag)*(f33_hag - (f33_hag**4)) + 4.95/(1.+2.48*Z2_hag)*((f33_hag**2)-(f33_hag**4) -0.55*(f33_hag**3)+0.55*(f33_hag**4)) - 1.2/(1.+0.5*Z2_hag)* (f33_hag**4)

			l34s_hag = (1.+1.4/(Z1_hag+1.))*f34_hag - 1.9/(Z1_hag+1.)* (f34_hag**2) + 0.3/(Z1_hag+1.)*(f34_hag**3) + 0.2/(Z1_hag+1.)*(f34_hag**4)

			if self.use_neo:
				l31s_hag = (1.+0.15/(Z2_hag**1.2-0.71))*f31_hag - 0.22/(Z2_hag**1.2-0.71)*f31_hag**2 + 0.01/(Z2_hag**1.2-0.71)*f31_hag**3 + 0.06/(Z2_hag**1.2-0.71)*f31_hag**4

				l32s_hag = (0.1+0.6*Z2_hag)/Z2_hag/(0.77+0.63*(1+(Z2_hag-1)**1.1))*(f32_hag-f32_hag**4)
				l32s_hag = l32s_hag + 0.7/(1+0.2*Z2_hag)*(f32_hag**2-f32_hag**4-1.2*(f32_hag**3-f32_hag**4))+1.3/(1+0.5*Z2_hag)*f32_hag**4

				l32s_hag = l32s_hag - (0.4+1.93*Z2_hag)/Z2_hag/(0.8+0.6*Z2_hag)*(f33_hag-f33_hag**4) + 5.5/(1.5+2.*Z2_hag)*(f33_hag**2-f33_hag**4-0.8*(f33_hag**3-f33_hag**4))
				l32s_hag = l32s_hag - 1.3/(1+0.5*Z2_hag)*f33_hag**4

				l34s_hag = l31s_hag

  
			LTI_hag = abs(self.dti[i] / self.ti[i])
			LN_hag = abs(self.dne[i] / self.ne[i])
			Lq_hag = abs(self.dq[i] /self.q[i])
			self.alps[i] = alps_hag
			if self.use_hager: gama_hag2 = gama_hag
			else: gama_hag2 = 1.

			F33TEF = self.ft[i] / (1. + (0.55-0.1*self.ft[i])*np.sqrt(self.nues[i]) + 0.45*(1.-self.ft[i])*self.nues[i]/self.Z**1.5)
			if self.use_neo:
				F33TEF = 1+0.25*(1-0.7*self.ft[i])*np.sqrt(self.nues[i])*(1+0.45*(Z2_hag-1.)**2)+0.61*(1-0.41*self.ft[i])*self.nues[i]/np.sqrt(Z2_hag)
				F33TEF = self.ft[i]/F33TEF
	
			ZNZ = 0.58 + 0.74/(0.76+Z2_hag)
			SIGSPITZ = 1.9012*1.e4 * (self.te[i]*1.e3)**1.5 / (self.Z * ZNZ * self.ZLE[i])
			
			X=F33TEF
			SIGNEO = 1. -(1.+0.36/Z2_hag)*X + 0.59/Z2_hag * X*X - 0.23/Z2_hag * X**3
			if self.use_neo:
				SIGNEO = 1. - (1.+0.21/Z2_hag)*X +0.54/Z2_hag*X**2-0.33/Z2_hag*X**3
			SIGNEO = SIGNEO * SIGSPITZ

			f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%13.6e\n'%(self.psin[i],-alps_hag*gama_hag2,self.nuis[i],self.nues[i],-alps_hag,self.b02av[i],SIGNEO))
			
			self.jb_hag1[i] = (gam1_hag*l31s_hag)/self.pe[i]*self.dpt[i] + (gam2_hag*l32s_hag) / self.te[i] * self.dte[i] + (gam4_hag*l34s_hag) * (gama_hag*alps_hag) / Z2_hag / self.te[i] * (1.-2.*DPSI_hag*(-1.5*LTI_hag + LN_hag + Lq_hag))*self.dti[i]

			if (self.use_chang):
				self.jb_hag2[i] = l31s_hag/self.pe[i]*self.dpt[i] + l32s_hag / self.te[i] * self.dte[i] + l34s_hag * alps_hag / Z1_hag / self.te[i] * self.dti[i]
			elif(self.use_hager):
				self.jb_hag2[i] = l31s_hag/self.pe[i]*self.dpt[i] + l32s_hag / self.te[i] * self.dte[i] + l34s_hag * alps_hag / Z2_hag / self.te[i] * self.dti[i]
			else:
				self.jb_hag2[i] = l31s_hag/self.ne[i]*self.dne[i] + rpe_hag*(l31s_hag+l32s_hag)/self.te[i]*self.dte[i] + (1.-rpe_hag)*(l31s_hag+l34s_hag * alps_hag) * self.dti[i]/self.ti[i]
			
			if (self.use_chang):
				self.ft[i] = self.ft[i] / (1-H)
		
		self.jb_hag1[0] = 0.0
		self.jb_hag2[0] = 0.0
		self.betacol[0] = 0.0
		self.betagradti[0] = 0.0
		self.cb_hag[0] = 0.0
		self.ZLE[0] = self.ZLE[1]
		self.nues[0] = self.nues[1]
		self.nuis[0] = self.nuis[1]

		betaw = np.zeros(self.num)	
		betaw0 = np.copy(betaw)
		first_run = True
		psi_damp =10.
		damp_limit = - self.betaw_low_limit
		for i in range(self.num):
			if self.betagradti[i] < damp_limit:
				psi_damp = self.psin[i]
				break

		if (psi_damp < 1.0): psi_damp = psi_damp - 0.03
		self.betagradti2 = np.copy(self.betagradti)
		for i in range(self.num):
			self.betagradti2[i] = (self.betagradti[i]-damp_limit) * 0.5 * (1.- np.tanh((self.psin[i]-psi_damp)/0.03)) + damp_limit
	
		for i in range(self.num):
	
			betaw[i]  = self.betacol[i] + self.betagradti2[i]
			betaw0[i] = self.betacol[i] + self.betagradti[i]

			self.jb_hag1[i] = abs( self.jb_hag1[i] * self.f[i] * self.pe[i] / self.BMAG) * betaw[i]
		
			if (self.use_hager or self.use_chang):
				self.jb_hag2[i] = abs( self.jb_hag2[i] * self.f[i] * self.pe[i] / self.BMAG)
			else:
				self.jb_hag2[i] = abs( self.jb_hag2[i] * self.f[i] * self.pt[i] / self.BMAG)
			
			if (self.use_hager):
				self.jb_hag = np.copy(self.jb_hag1)
			else:
				self.jb_hag = np.copy(self.jb_hag2)
			
		if (self.hag_mod and self.use_hager):

			hag_over_saut = np.zeros(self.num)			
			hag_mod_index = 0
			
			for i in range(self.num):
				if (self.psin[i] > 0.0):
					hag_over_saut[i] = self.jb_hag1[i] / self.jb_hag2[i]
				else:
					hag_over_saut[i] = 1.0

				if (self.psin[i] <= self.hag_mod_psin):
					hag_mod_index = i

			hag_mod_index = hag_mod_index + 1

			temp2_hag = (hag_over_saut[hag_mod_index+1] - hag_over_saut[hag_mod_index])/(self.psin[hag_mod_index+1] - self.psin[hag_mod_index])
			temp1_hag = (1.-hag_over_saut[hag_mod_index]+temp2_hag*self.hag_mod_psin)
			
			for i in range(self.num):

				if (self.psin[i] <= self.hag_mod_psin):
				
					hag_weight =  temp1_hag*(((self.psin[i]/self.hag_mod_psin) - 1.)**2) + temp2_hag*(self.psin[i] - self.hag_mod_psin) + hag_over_saut[hag_mod_index]
				else:
					hag_weight = hag_over_saut[i]

				hag_weight2 = 0.5 * (1.0 + np.tanh((self.psin[i] - self.hag_mod_psin)/self.hag_mod_psinw))

				if (self.psin[i] < 0.03):
					hag_weight2 = hag_weight2 + (0.-hag_weight2)*(1.0 - np.tanh((self.psin[i] - 0.02)/0.01))*0.5
				
				self.jb_hag[i] = self.jb_hag2[i] * (1.0 - hag_weight2) + self.jb_hag1[i] * hag_weight2
		f.close()

		#multiplication
		self.jb_hag = self.jb_hag * self.bsmulti

		self.neoclass_flow3()
		return

	def neoclass_flow(self):

		me = 9.1 * 1.e-31
		mi = 1.672 * 1.e-27
		mI = self.aimp * mi
		e0 = 1.602 * 1.e-19

		R = np.copy(self.r + self.rc)
		nef = interp1d(R,self.ne,'cubic')
		tif = interp1d(R,self.ti,'cubic')
		
		deps = 1.e-6
		Lti = 1.e2	
		Lne = 1.e2
		
		for i in range(self.num):

			ni = self.ne[i] * (self.zimp - self.zeffs[i]) / (self.zimp - 1.0)
			nI = self.ne[i] * (self.zeffs[i] - 1.0) / (self.zimp - 1.0) / self.zimp

			vti = np.sqrt(2.0 * self.ti[i] * 1.e3 * e0 / mi)
			vtI = np.sqrt(2.0 * self.ti[i] * 1.e3 * e0 / mI)

			r = self.r[i] + self.rc[i]

			Bt = self.f[i] / (self.r[i] + self.rc[i])
			rhoi = vti * mi / e0 / Bt
		
			alpha = nI / ni * (self.zimp ** 2)
			beta = ((27./4.)**2) * ((mi / mI)**2) / (15./2. + np.sqrt(2.0*alpha)*vti/vtI)

			g = self.ft[i] / (1.0 - self.ft[i])

			mui00 = g * (alpha + np.sqrt(2.0) - np.log(1.0+np.sqrt(2.0)))
			mui01 = g * (1.5*alpha + 4./np.sqrt(2.0) - 5.0/2.0 * np.log(1.0+np.sqrt(2.0)))
			mui10 = mui01
			mui11 = g * (13./4.*alpha + 39./4./np.sqrt(2.) - 25./4.* np.log(1.0+np.sqrt(2.0)))

			D = mui00 * (mui11 + np.sqrt(2.) + alpha - alpha*beta) - (mui01**2)

			if (D == 0.0):
				D = 1.0

			K1 = mui01*(np.sqrt(2.)+alpha - alpha*beta)/D
			K2 = (mui00*mui11 - (mui01**2))/D

			if (i > 0 and i < (self.num-1)):
				Lti = 2.0 * deps / (tif(r+deps) - tif(r-deps)) * tif(r)
				Lne = 2.0 * deps / (nef(r+deps) - nef(r-deps)) * nef(r)

			Lpi = Lti * Lne / (Lti + Lne)
			LpI = Lpi

			Vthetai = 0.5 * vti * rhoi * K1 / Lti * Bt * Bt / (self.b02av[i]**2)
			VthetaI = 0.5 * vti * rhoi * ((K1+1.5*K2)/Lti - 1./Lpi + 1./self.zimp/LpI) * Bt * Bt / (self.b02av[i]**2)

		return
	
	def neoclass_flow2(self):

		eps0 = 8.85 * 1.e-12

		for i in range(self.num):

			ne = self.ne[i] * 1.e19
			te = self.te[i] * 1.e3
			ni = ne * (self.zimp - self.zeffs[i]) / (self.zimp - 1.0)
			nI = ne * (self.zeffs[i] - 1.0) / (self.zimp - 1.0) / self.zimp
			ti = self.ti[i] * 1.e3

			LneeNRL = 24.-0.5*np.log(ne*1.e-6) + np.log(te)
			LneiNRL = 30.+np.log(self.amain)-2.0*np.log(1.0) - 0.5 * np.log(ni*1.e-6) + 1.5 * np.log(ti)
			LnejNRL = 30.+np.log(self.aimp)-2.0*np.log(self.zimp) - 0.5 * np.log(ne*1.e-6) + 1.5 * np.log(ti)
			LniiNRL = 23.+np.log(2.0*0.5) - 2.0*np.log(1.0) - 0.5 * np.log(ne*1.e-6) + 1.5 * np.log(ti)
			LnijNRL = 23.+np.log((self.aimp+self.amain)/self.amain*0.5) - np.log(self.zimp) - 0.5 * np.log(2.0*ne*1.e-6) + 1.5*np.log(ti)
			LnjjNRL = 23.+np.log(2.0*0.5) - 2.0*np.log(self.zimp) - 0.5 * np.log(ne*1.e-6) + 1.5 * np.log(ti)

			Tauee = 3 * eps0 * np.sqrt(1)

		return

	def neoclass_flow3(self,filename='Vneo.dat'):
		f = open(filename,'w')
		for i in range(self.num): f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.psin[i],self.alps[i],self.dti[i],self.nuis[i]))
		f.close()
		return

	def core_cond(self):
	
		self.sig = np.zeros(self.num)
		for i in range(self.num):

			zeff = self.zeffs[i]
	
			F33TEF = self.ft[i] / (1. + (0.55-0.1*self.ft[i])*np.sqrt(self.nues[i]) + 0.45*(1.-self.ft[i])*self.nues[i]/self.Z**1.5)
			if self.use_neo:
				F33TEF = 1+0.25*(1-0.7*self.ft[i])*np.sqrt(self.nues[i])*(1+0.45*(zeff-1.)**2)+0.61*(1-0.41*self.ft[i])*self.nues[i]/np.sqrt(zeff)
				F33TEF = self.ft[i]/F33TEF
	
			ZNZ = 0.58 + 0.74/(0.76+zeff)
			SIGSPITZ = 1.9012*1.e4 * (self.te[i]*1.e3)**1.5 / (self.Z * ZNZ * self.ZLE[i])
			
			X=F33TEF
			SIGNEO = 1. -(1.+0.36/zeff)*X + 0.59/zeff * X*X - 0.23/zeff * X**3
			if self.use_neo:
				SIGNEO = 1. - (1.+0.21/zeff)*X +0.54/zeff*X**2-0.33/zeff*X**3
			SIGNEO = SIGNEO * SIGSPITZ
			
			self.sig[i] = SIGNEO
	
		return
		
	def core_current(self):
	
		self.jcore = np.zeros(self.num)
		
		for i in range(self.num):
			self.jcore[i] = self.f[i] / 2.0 / np.pi  * self.sig[i] / self.BMAG * self.gm1[i]

		for i in range(self.num):
			if (self.psin[i] > self.core_neo):
				break
		if (i<(self.num-1)):
			self.j_ind = i
			j = i
		else:
			self.j_ind = 0
			j = 0

		for i in range(j):
			self.jcore[i] =  (self.jcore[j+1] - self.jcore[j]) / (self.psin[j+1] - self.psin[j]) * (self.psin[i] - self.psin[j]) + self.jcore[j]

		jcore = self.jcore[0]
		
		if (self.use_beam_current):
			for i in range(self.num):
				if (self.psin[i] < self.bjf):
					if not (self.absolute_ajf):
						self.jcore[i] = self.jcore[i] + jcore * self.ajf * ( 1 - (self.psin[i]/self.bjf)**self.cjf) ** self.djf
						self.ajf0 = self.ajf  * jcore
					else:
						self.jcore[i] = self.jcore[i] + self.ajf * ( 1 - (self.psin[i]/self.bjf)**self.cjf) ** self.djf

				if self.use_eped3:
					self.jcore[i] = jcore * ( 1 - (self.psin[i]/self.bjf)**self.cjf) ** self.djf

		if (self.use_eped2 or self.use_eped3): 
			if (self.core_neo > 0.):
				for i in range(self.j_ind):
					self.jcore[i] =  (self.jcore[j+1] - self.jcore[j]) / (self.psin[j+1] - self.psin[j]) * (self.psin[i] - self.psin[j]) + self.jcore[j]


		return
		
	def init_current(self):
	
		self.bscurt = np.trapz(self.jb_hag,x=self.area)
		
		self.corecurt = np.trapz(self.jcore,x=self.area)
		
		self.excurt = np.trapz(self.curr_ex,x=self.area)

		self.VLOOP = (abs(self.ip) - self.bscurt - self.excurt)/self.corecurt

		if (self.given_j):
			self.corecurt = np.trapz(self.zjzg,x=self.area)
			self.VLOOP = (abs(self.ip) - 0.*self.bscurt - self.excurt)/self.corecurt

		if (self.VLOOP < 0.0):
			self.VLOOP = 5.e-2
		
		return
	
	def make_current(self,type=1):
		if not(len(self.zjz) == 2):
			self.zjz_old = np.copy(self.zjz)
		self.zjz = np.zeros(self.num)
		
		for i in range(self.num):
			
			self.zjz[i] = self.VLOOP * (1.0 - (1.0-self.vloop_core_mod)*((1.0-self.psin[i]**2.0)))* self.jcore[i] + self.jb_hag[i]

			if (self.use_ext_current):
				self.zjz[i] = self.zjz[i] +  self.curr_ex[i] * self.gm9[i] / self.gm1[i] * self.b02av[i] / self.f[i] / self.BMAG	#Change NUBEAM Current to jdotB/Bmag

		if (self.core_neo > 0.):
			j = self.j_ind
			for i in range(self.j_ind):
				self.zjz[i] =  (self.zjz[j+1] - self.zjz[j]) / (self.psin[j+1] - self.psin[j]) * (self.psin[i] - self.psin[j]) + self.zjz[j]

		if type == 1:

			f = open(self.chease_rundir+'/curr_prof','w')
			f.write('%9s\t%9s\t%9s\t%9s\t%9s\t%9s\n'%('PSIN','INDUCTION','BS_SAUTER','BS_TARGET','EXTERNAL','TOTAL'))

			for i in range(self.num):
				ww = self.BMAG / self.rc[-1] / self.f[i] / (self.gm1[i]) / (self.BMAG * self.rc[i]/self.f[i])
				j1 = self.jcore[i]*self.VLOOP
				j2 = self.jb_hag[i]
				j3 = self.jb_hag2[i] * self.bsmulti
				j4 = self.curr_ex[i] * self.gm9[i] / self.gm1[i] * self.b02av[i] / self.f[i] / self.BMAG
				if (len(self.zjz_old) >3):
					f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.psin[i],j1,j2,j3,j4,self.zjz[i],self.zjz_old[i]))
			f.close()

		self.zjzc = np.zeros(self.num)
		self.zjzbsc = np.zeros(self.num)

		if (len(self.zjz_old) == 2):
			self.zjz_old = np.copy(self.zjz)

		mod = self.zjz[0] * self.relax + self.zjz_old[0] * (1.0-self.relax)

		if mod > 1.e5:	dv = 1.e6
		elif mod > 1.e4:	dv = 1.e5
		else:	dv = 1.e4

		print('new_j0','%4.2f'%(self.zjz[0]/dv),'old_j0','%4.2f'%(self.zjz_old[0]/dv),'mod_j0','%4.2f'%(mod/dv))

		self.zjz = self.zjz * self.relax + self.zjz_old * (1.0-self.relax)

		if (self.given_j):
			self.zjz = self.VLOOP * np.copy(self.zjzg)
	
		for i in range(self.num):
			self.zjzc[i] = self.zjz[i] * self.BMAG / self.rc[-1] / self.f[i] / (self.gm1[i])
			self.zjzbsc[i] = self.jb_hag[i] * self.BMAG / self.rc[-1] / self.f[i] / (self.gm1[i])
	
		return

	def make_current2(self):
		self.relax_old = self.relax
		self.relax = 1.0
		self.make_init_current(self.hager_mapf,None,2)
		vloop0 = self.VLOOP
		self.relax = 1.0;

		ip_crit = abs(self.ip)
		errcrit = 1.e-3
		s = np.sign(self.eq.sbdy - self.eq.smag)
		R0 = self.rc[0]
		jtor = np.copy(self.zjz)
		ffp = np.copy(self.zjz)
		for i in range(self.num):	
			jpar = self.zjz[i] * self.BMAG / self.f[i] * R0 / s
			ffp[i] = -(jpar + R0*self.dpp[i]) * mu0  / R0 / self.b02av[i] * (self.f[i]**2)
			jtor[i] = -s * (R0*self.dpp[i] + ffp[i] / mu0 * R0 * self.gm1[i]) / (R0 * self.gm9[i])

		ip0 = np.trapz(jtor,x=self.area)

		err = (ip0 -ip_crit)/ip_crit

		if abs(err) > errcrit:
			if err < 0.:	vloop1 = vloop0 * 1.2
			else:	vloop1 = vloop0 * 0.8
		else: vloop1 = vloop0
		print('>>> Ip = %f [MA] for Vloop = %f [V], Target = %f [MA] err= %f'%(ip0/1.e6,vloop0,ip_crit/1.e6,err))
		while abs(err) > errcrit:

			self.VLOOP = vloop1
			self.make_current(2)

			for i in range(self.num):	
				jpar = self.zjz[i] * self.BMAG / self.f[i] * R0 / s
				ffp[i] = -(jpar + R0*self.dpp[i]) * mu0  / R0 / self.b02av[i] * (self.f[i]**2)
				jtor[i] = -s * (R0*self.dpp[i] + ffp[i] / mu0 * R0 * self.gm1[i]) / (R0 * self.gm9[i])

			ip1 = np.trapz(jtor,x=self.area)

			err = (ip1 -ip_crit)/ip_crit	

			vloop = (vloop1 - vloop0) / (ip1-ip0) * (ip_crit-ip0) + vloop0 
			vloop0 = vloop1
			vloop1 = vloop
			ip0 = ip1
			print('>>> Ip = %f [MA] for Vloop = %f [V], Target = %f [MA] err= %f'%(ip0/1.e6,vloop0,ip_crit/1.e6,err))
		self.make_current(1)
		self.relax = self.relax_old
		return (jtor,ffp,vloop1)

	def read_given_j(self,filename):
		
		f = open(filename,'r')

		line = f.readline()
		num = int(line)

		dat = np.zeros(shape=(num,2))

		for i in range(num):
			line = f.readline().split()
			dat[i,0] = float(line[0])
			dat[i,1] = float(line[1])

			jfit = interp1d(dat[:,0],dat[:,1])
		if not self.use_rho:
			dat2 = jfit(self.psin)
		else:	
			dat2 = jfit(self.rho)

		f.close()

		return dat2

	def make_given_j(self):
		isfile = os.path.isfile('extj')
		if not isfile:
			print('No external given j profile !!!')
			exit()
		if (self.zjzg[0] == 0):
			print('USING EXTERNAL GIVEN J!!')
		dat = self.read_given_j('extj')
		self.zjzg = 1.0 * np.copy(dat)

	def make_init_current(self,filename,vloop=None,type=1):
	
		self.read_hagermap(filename)
		
		self.interpol()
		self.hager_bs()
		self.core_cond()
		self.core_current()
		if (self.given_j):
			self.make_given_j()
		
		if not (self.efit_const):
			self.init_current()
		else:
			if not (self.vloop_ext == 0.0):
				self.VLOOP = self.vloop_ext
				self.bscurt = 0.0
				self.excurt = 0.0
				self.corecurt = 1.0
			else:
				self.init_current()

		if not vloop == None:
			self.VLOOP = vloop
			
		self.make_current(type)
		
		return	
				
	def make_chease_input(self,filename1,filename2,nideal=11):
	
		zjzcf = interp1d(self.psin,self.zjzc,'cubic')
		zjzc = zjzcf(self.eq.psin)
	
		self.eq.jpar = np.copy(zjzc)
		
#		self.eq.target_psin = self.target_psin;
		self.eq.ncscal = 2
		if(self.use_chease_jprof):
			self.eq.jconst = 3
		else:
			self.eq.jconst = 1
		
		self.eq.nideal = nideal
		self.eq.epslon = self.epslon	
	
		self.eq.make_chease_input(filename1,filename2)
	
		return
		
	def run_chease(self,nideal=11,input_only=False):
	
		try:
			os.mkdir(self.chease_rundir)
		except:
			pass
			
		self.make_chease_input(self.chease_rundir+'/EXPEQ',self.chease_rundir+'/chease_namelist',nideal)
		if (input_only):
			return
		os.chdir(self.chease_rundir)
		chease_run = chease_exec + ' > log.chease'
		stat, out = subprocess.getstatusoutput(chease_run)

		f4 = open('log.chease','r')
		linec = 0
		betc = 0
		while True:
			line = f4.readline()
			if not line: break
			
			if (line.find('MKSA') > -1):
				linec = 1
			
			if ((line.find('TOTAL CURRENT -->') > -1) and (linec == 1)):
				ipc = float(line.split()[5])
				
			if ((line.find('POLOIDAL BETA') > -1) and (linec == 1 )):
				bpc = float(line.split()[0])
			if ((line.find('BETA_EXP') > -1) and (linec == 1 )):
				betc = float(line.split()[0])
			if ((line.find('WMHD')>-1) and (linec == 1)):
				self.wmhd = float(line.split()[0])
			if ((line.find('LI')>-1) and (linec == 1)):
				self.li = float(line.split()[0])		
				
		self.bp = bpc

		self.beta = betc * 100.
		a0 = 0.5*(max(self.eq.rzbdy[:,0])-min(self.eq.rzbdy[:,0]))
		R0 = 0.5*(max(self.eq.rzbdy[:,0])+min(self.eq.rzbdy[:,0]))
		B0 = self.f[-1]/R0
		self.betan = self.beta * a0 * abs(B0/self.ip*1.e6)
		
		self.read_rho_psi_R()
		
		os.chdir(self.currdir)
		
		return(ipc,bpc)
		
	def read_rho_psi_R(self,filename='RHO_PSI_R'):
	
		currdir2 = os.getcwd()
		f4 = open(filename,'r')
		line = f4.readline()
		line = f4.readline()
		self.num_map = int(line)
		
		self.rho_map = np.zeros(self.num_map)
		self.psi_map = np.zeros(self.num_map)
		self.R_map = np.zeros(self.num_map)
		
		for i in range(self.num_map):
		
			line = f4.readline().split()
			self.rho_map[i] = float(line[0])
			self.psi_map[i] = float(line[1])
			self.R_map[i] = float(line[2])	
			
		f4.close()

		rhop = interp1d(self.psi_map,self.rho_map,'cubic')
		self.rho = rhop(self.psin)

		return
		
	def adjust_psin(self,R1,psin1,R2,psin2,psint):
	
		psinf1 = interp1d(psin1,R1,'cubic')
		R3 = np.zeros(len(R2))
		
		for i in range(len(R3)):
			R3[i] = (R1[-1] - R1[0])/(R2[-1] - R2[0])*(R2[i] -R2[0]) + R1[0]
			
		Rf2 = interp1d(R3,psin2,'cubic')

		R = psinf1(psint)
		psi = Rf2(R)
		psi[0] = 0.0
		psi[-1] = 1.0
		
		return psi
		
	def adjust_kinprof(self,R1,psin1,R2,psin2,psink,prof):
	
		Rf1 = interp1d(psin1,R1,'quadratic')
		Rf2 = interp1d(psin2,R2,'quadratic')
		Rk1 = Rf1(psink)
		Rk2 = Rf2(psink)
		prof2 = np.copy(psink)
		Rk2 = Rk2 - Rk2[-1] + Rk1[-1]
		
		proff = interp1d(Rk1,prof,'cubic')
		ind = np.where(Rk2>=Rk1[0])
		sindex = len(Rk1) - len(ind[0])
		dp0 = (prof[1]-prof[0])/(psink[1]-psink[0])
		
		
		for i in range(sindex,len(Rk1)):
			prof2[i] = proff(Rk2[i])

		dp1 = (prof2[sindex]-prof2[sindex+1])/(psink[sindex]-psink[sindex+1])
		dp02 = (prof[sindex]-prof[sindex+1])/(psink[sindex]-psink[sindex+1])
		dp0 = dp1 - (dp02-dp0)/(dp02-0.1)*(dp1-0.1)
		
		#for i in range(6):	print(psink[i],prof2[i],dp0,dp1)

		if not sindex == 0:
			for i in range(sindex):
				prof2[i] = (dp1-dp0)/psink[sindex]*0.5*(psink[i]**2 - psink[sindex]**2) + dp0*(psink[i] - psink[sindex]) + prof2[sindex]

		dp0 = (prof2[1]-prof2[0])/(psink[1]-psink[0])
		ind = np.where(psink>0.05)
		prof3 =  prof2[ind]
		psink2 = psink[ind]
		dp1 = (prof3[1]-prof3[0])/(psink2[1]-psink2[0])
		x1 = (psink2[1]+psink2[0])/2.
		dp2 = (prof3[2]-prof3[1])/(psink2[2]-psink2[1])
		x2 = (psink2[2]+psink2[1])/2.
		dp3 = (dp2-dp1)/(x2-x1) * (0.-x1) + dp1

#		if dp0*dp3 < 0.:	return

		if (abs(dp1) < abs(dp0) and abs(dp3) < abs(dp0)):
			for i in range(ind[0][0]):
				x = psink[i]
				x3 = psink[ind[0][0]]
				prof2[i] = prof2[ind[0][0]] + (dp2-dp1)/(x2-x1)*0.5*((x-x1)**2 - (x3-x1)**2) + dp1*(x-x3)
				
		return prof2
		
	def make_param_bnd(self):

		self.chi = np.linspace(0,2*np.pi,512)
		self.rzbdyp = np.zeros(shape=(512,2))

		for i in range(512):

			theta = self.chi[i]
			
			self.rzbdyp[i,0] = self.rcent + self.amin * np.cos(theta + self.triang * np.sin(theta) + self.square*np.sin(2.0*theta))
			self.rzbdyp[i,1] = self.elong * self.amin * np.sin(theta + 0.0*self.square*np.cos(theta))
				
		return

	def read_bnd_shape(self,filename):
	
		f4 = open(filename,'r')
		line = f4.readline()
		bnd_num = int(line)

		self.rzbdyp = np.zeros(shape=(bnd_num,2))

		for i in range(bnd_num):
			line = f4.readline().split()
			self.rzbdyp[i,0] = float(line[0])
			self.rzbdyp[i,1] = float(line[1])
		f4.close()
		
		self.rcent = (max(self.rzbdyp[:,0])+min(self.rzbdyp[:,0]))*0.5

		return
	
	def fast_ion_pressure(self,psin):
	
		len2 = len(psin)
		
		self.fpt = np.zeros(len2)
		self.fptp = np.zeros(len2)
		
		for i in range(len2):
		
			if (psin[i] <= self.bpf):
				self.fpt[i] = self.apf * ((1.-((psin[i]/self.bpf)**self.cpf)) ** self.dpf) * self.pt[0]
				self.fptp[i] = - self.cpf * self.apf * self.dpf  / (self.bpf**self.cpf) * (psin[i] ** (self.cpf-1.0)) * ((1.-((psin[i]/self.bpf)**self.cpf)) ** (self.dpf-1.0)) * self.pt[0]  / self.psia
				
			else:
				self.fpt[i] = 0.0
				self.fptp[i] = 0.0
				
		return
	
	def make_pressure_prof(self):
	
		self.fast_ion_pressure(self.psin)
		self.dpp = np.zeros(len(self.psin))
		self.dpres = np.zeros(len(self.psin))

		ind = np.where(self.psin>0.03)
		ind0 = ind[0][0];	ind1 = ind[0][1]
		for i in range(ind0):
			self.dp_ex[ind0-i-1] = (self.dp_ex[ind1]-self.dp_ex[ind0])/(self.psin[ind1]-self.psin[ind0])*(self.psin[ind0-i-1]-self.psin[ind0])+self.dp_ex[ind0]
		
		for i in range(len(self.psin)):
		
			self.dpp[i] = self.dpt[i] + self.fptp[i]
			self.dpres[i] = self.pt[i] + self.fpt[i]

			
			if (self.use_ext_pressure):	
				self.dpp[i] = self.dpp[i] + self.dp_ex[i]
				self.dpres[i] = self.dpres[i] + self.pres_ex[i]
			if (self.dpp[i] * self.psia > 0.): self.dpp[i] = 0.
	
		self.p_thermal = np.trapz(self.pt,x=self.area)
		if not(self.use_ext_pressure):
			self.p_fast = np.trapz(self.fpt,x=self.area)
		else:
			self.p_fast = np.trapz(self.fpt+self.pres_ex,x=self.area)

		pedloc = 1-self.ped_width
		for i in range(len(self.psin)-1):
			if (self.psin[i]-pedloc)*(self.psin[i+1]-pedloc)<=0: break
		w = (pedloc-self.psin[i])/(self.psin[i+1]-self.psin[i]);
		self.pped = w*self.pt[i+1]+(1.-w)*self.pt[i];

		try:	f = open(self.chease_rundir+'/pres_prof','w')
		except:	f = open(os.getcwd()+'/pres_prof','w')
		f.write('%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\n'%('PSIN','PREST','PRESF','PRESEX','PTOT','DPREST','DPRESF','DPRESEX','DPTOT'))
		for i in range(len(self.psin)):
			f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.psin[i],self.pt[i],self.fpt[i],self.pres_ex[i],self.dpres[i],self.dpt[i],
				self.fptp[i],self.dp_ex[i],self.dpp[i]))
		f.close()

		return

	def current_iteration(self,vloop=None):

		start_time = time.time()	
		ip0 = abs(self.ip)
		self.load_kin_profile(False)
		self.make_init_current(self.hager_mapf,vloop)
		
		psia = self.psia
		
		vloop0 = self.VLOOP
		ip1,self.bp1c = self.run_chease()
		print('iteration #0, IP=%6.3f [MA]'%(ip1/1.e6))
		niter = 0
		
		bs_old = np.copy(self.jb_hag)
		
		err = (ip1-ip0)/ip0
		err1 = 0.1
		while ((abs(err) > self.ip_crit or abs(err1) > self.bs_crit) and niter < self.niterc):
		
			vloop1 = vloop0 - vloop0 * (ip1-ip0)/ip0 * (self.bscurt + self.excurt + self.corecurt) /self.corecurt
			if (vloop1 < 0.0):
				vloop1 = 1.e-2
			
			self.load_kin_profile(False)			########## it can be numerically consuming and vulnerable
			
			self.make_init_current(self.hager_mapf,vloop1)
			bs_new = np.copy(self.jb_hag)
			bs_diff = np.copy(bs_new - bs_old) / max(bs_old)
			
			self.make_pressure_prof()		########## it can be numerically consuming and vulnerable
			self.eq.psin = np.copy(self.psin)
			self.eq.pp = np.copy(self.dpp)
			self.eq.pres = np.copy(self.dpres)

			if not self.ped_width == 0.05:
				self.eq.width = self.ped_width
			
			ip1,self.bp1c = self.run_chease()
			self.bs_tor_calc()
			
			err = (ip1-ip0)/ip0
			err1 = max(abs(bs_diff))
			vloop0 = vloop1
			niter = niter + 1
			print('iteration #%2i,bp = %6.3f, bs= %6.3f [MA], errIP = %9.6f, errBS = %9.6f, li = %6.3f, vloop = %9.6f [V]'%(niter,self.bp1c,self.bs_tot/1.e6,err,err1,self.li,self.VLOOP))
			bs_old = np.copy(bs_new)

		self.VLOOP2 = vloop1
		self.load_kin_profile(False)
		self.make_init_current(self.hager_mapf,vloop1)
		
		fastpv =  np.trapz(self.pres_ex,x=self.vol)
		fastppv = np.trapz(self.pres_exp,x=self.vol)
		
		self.wdia = self.wmhd - 1.5*(fastpv - fastppv)*1.e-3
		self.wfast = 1.5*fastpv*1.e-3
				
		self.ip1 = ip1
		self.zjz_old = np.zeros(2)
		print('Elapsed time = %f [s]'%(time.time()-start_time))
			
		return
		
	def init_current_iteration_eped(self):
			
		print('Initial guess calculation for EPED Equil')
		start_time = time.time()	
		if (self.eqdsk_name.lower() == 'none'):
			ip0 = abs(self.ip)

		self.eped_first_run = True
		self.load_kin_profile(False)
		self.eped_first_run = False
	
		self.make_pressure_prof()
		self.eq = eqdsk.eqdsk(self.eqdsk_name)
		if (self.target_psin == 1.0):
			self.eq.target_psin = 0.0
		else:
			self.eq.target_psin = self.target_psin

		if not (self.eqdsk_name.lower() == 'none'):
			self.eq.read_eqdsk_file()

		if (self.eqdsk_name.lower() == 'none'):
			if(self.use_param_shape):
				self.make_param_bnd()
				self.eq.use_bnd_smooth = True
			else:
				self.read_bnd_shape(self.bnd_file)
				self.eq.use_bnd_smooth = True
		
			self.eq.device_name = 'absolute'
			self.eq.rzbdy = np.copy(self.rzbdyp)
			self.eq.target_psin = 0.0
#			self.eq.use_bnd_smooth = False

		self.zjzc = np.linspace(1,0,self.num)
	
		self.eq.psin = np.copy(self.psin)
		self.eq.pp = np.copy(self.dpp)
		self.eq.pres = np.copy(self.dpres)
		self.eq.ffp = np.copy(self.dpres)

		for i in range(self.num):
			self.eq.pres[i] = self.dpres[0]*0.5*(1.0 - 2.0*self.psin[i] +  self.psin[i]**2)
			self.eq.ffp[i] = -(1.0 - 0.5*self.psin[i] - 0.5*(self.psin[i]**2))

		if (self.eqdsk_name.lower() == 'none'):
			self.eq.q = np.copy(self.q)
			self.eq.ip = self.ip
			self.eq.bcentr = self.bcent
			self.eq.rcentr = self.rcent

		self.eq.use_pres_chease = True
		self.use_chease_jprof = False

		ip1,self.bp1c = self.run_chease()
		self.eq.use_pres_chease = False
		self.use_chease_jprof = True

		self.rho_mapi = np.copy(self.rho_map)
		self.psi_mapi = np.copy(self.psi_map)
		self.R_mapi = np.copy(self.R_map)

		self.read_hagermap(self.hager_mapf)
		self.interpol()

		self.eq.q = np.copy(self.q)

		print('iteration #0, Bp = %5.3f, IP=%5.3f [MA]'%(self.bp1c,ip1/1.e6))
				
		print('Elapsed time = %f [s]'%(time.time()-start_time))
		print('Initial EPED Equil finished')

		return		
		
	def beta_iteration(self,criterion=None,map=True):
		# Beta iteration by adjusting non-thermal component
		try: os.mkdir('OUTPUT')
		except: pass	
		
		start_time = time.time()	
		if (self.apf == 0. and self.beta_criterion > 0.):
			self.apf = 0.3
	
		ap0 = self.apf

		if not(self.eqdsk_name.lower() == 'none'):
			self.load_eqdsk()
			f5 = open('OUTPUT/nubeam_iter_result','w')
			f5.write('iter#       bp              li           Wmhd[kJ]        Wdia[kJ]       Wfast[kJ]        Raxis[m]\n')
			f5.write('%2i\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(0,self.bp,self.li,self.wmhd,self.wdia,0.0,self.eq.rmag))
		else:
			self.init_current_iteration_eped()
			
		if (criterion == None):
			criterion = self.beta_critval	
		
		self.load_kin_profile(False)
		self.make_pressure_prof()
		self.eq.psin = np.copy(self.psin)
		self.eq.pp = np.copy(self.dpp)
		self.eq.pres = np.copy(self.dpres)
		self.current_iteration()
		
		if (self.beta_criterion == 1):
			err = abs(self.bp1c - criterion) /  criterion
			cval0 = self.bp1c
			print('Beta iteration #0 crit = %5.3f, val = %5.3f, err = %9.6f'%(criterion,cval0,err))
			
		elif(self.beta_criterion == 2):
			err = abs(self.rav[0] - criterion) /  criterion
			cval0 = self.rav[0]
			print('Beta iteration #0 crit = %5.3f, val = %5.3f, err = %9.6f'%(criterion,cval0,err))

		elif(self.beta_criterion == 3):
			err = abs(self.betan - criterion) /  criterion
			cval0 = self.betan
			print('Beta iteration #0 crit = %5.3f, val = %5.3f, err = %9.6f'%(criterion,cval0,err))
			
		elif(self.beta_criterion == 0):
			err = 1.0
			cval0 = self.bp1c
			cval = self.ip1
			print('Beta = %5.3f, Ip = %5.3f'%(cval0,cval))
			
		niter = 0
		
		if (self.beta_criterion == 1 or self.beta_criterion == 3):
			ap1 = ap0 + ap0 * (criterion - cval0)/criterion * (self.p_fast + self.p_thermal) / self.p_fast
		elif(self.beta_criterion == 2):
			ap1 = ap0 + ap0 * (criterion - cval0)/criterion * (self.p_fast + self.p_thermal) / self.p_fast		
		if (self.beta_criterion > 0):	
			self.apf = ap1
			self.make_pressure_prof()
			self.eq.psin = np.copy(self.psin)
			self.eq.pp = np.copy(self.dpp)
			self.eq.pres = np.copy(self.dpres)
			self.current_iteration(vloop=self.VLOOP)

		niter = niter + 1	

		if (self.beta_criterion == 1):
			err = abs(self.bp1c - criterion) /  criterion
			cval = self.bp1c
			print('Beta iteration #%2i crit = %f, val = %f, apf = %f, li = %f, err = %f'%(niter,criterion,cval,ap1,self.li,err))
			
		elif(self.beta_criterion == 2):
			err = abs(self.rav[0] - criterion) /  criterion
			cval = self.rav[0]	
			print('Beta iteration #%2i crit = %f, val = %f, apf = %f, li = %f, err = %f'%(niter,criterion,cval,ap1,self.li,err))

		elif (self.beta_criterion == 3):
			err = abs(self.betan - criterion) /  criterion
			cval = self.betan
			print('Beta iteration #%2i crit = %f, val = %f, apf = %f, li = %f, err = %f'%(niter,criterion,cval,ap1,self.li,err))

		while (abs(err) > self.beta_crit and self.beta_criterion > 0 and niter < self.niterb):
			
			ap2 = ap1 + (criterion - cval) / (cval - cval0) * (ap1 - ap0)
			
			cval0 = cval
			ap0 = ap1
			
			self.apf = ap2
			
			self.make_pressure_prof()
			self.eq.psin = np.copy(self.psin)
			self.eq.pp = np.copy(self.dpp)
			self.current_iteration(vloop=self.VLOOP)
			
			if (self.beta_criterion == 1):
				err = abs(self.bp1c - criterion) /  criterion
				cval = self.bp1c
			elif(self.beta_criterion == 2):
				err = abs(self.rav[0] - criterion) /  criterion
				cval = self.rav[0]
			elif(self.beta_criterion == 3):
				err= abs(self.betan - criterion) /  criterion
				cval = self.betan	
			
			ap1 = ap2
			niter = niter + 1
			print('Beta iteration #%2i crit = %f, val = %f, apf = %f, li = %f, err = %f'%(niter,criterion,cval,ap1,self.li,err))
			
		if not (self.nideal == 0):
			if map:
				self.mapping_run()			

		print('Write final kinetic profile')
		kinprof_type_old = self.kinprof_type
		self.kinprof_type = 1.
		self.write_kinprof()
		self.kinprof_type = kinprof_type_old

		copyfile(self.chease_rundir+'/EQDSK_COCOS_02.OUT','OUTPUT/geqdsk_0')
		copyfile(self.chease_rundir+'/EQDSK_COCOS_02.OUT','OUTPUT/geqdsk_1')
		copyfile(self.chease_rundir+'/chease_namelist','OUTPUT/chease_namelist_0')
		copyfile(self.chease_rundir+'/chease_namelist','OUTPUT/chease_namelist_1')
		copyfile(self.chease_rundir+'/EXPEQ','OUTPUT/EXPEQ_0')
		copyfile(self.chease_rundir+'/EXPEQ','OUTPUT/EXPEQ_1')
		copyfile(self.chease_rundir+'/EFIT_JCONST','OUTPUT/EFIT_JCONST_0')
		copyfile(self.chease_rundir+'/EFIT_JCONST','OUTPUT/EFIT_JCONST_1')
		copyfile('Vneo.dat','OUTPUT/Vneo.dat_0')
		copyfile('Vneo.dat','OUTPUT/Vneo.dat_1')

		if not(self.eqdsk_name.lower() == 'none'):
			f5.write('%2i\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(1,self.bp,self.li,self.wmhd,self.wdia,0.0,self.rav[0]))
			f5.close()

		f5 = open('OUTPUT/run_success','w')
		f5.write('1\n')
		f5.close()		
		
		print('Elapsed time = %f [s]'%(time.time()-start_time))
		return

	def beta_iteration_eped(self,criterion=None):
		# beta iteration by adjusting thermal pressure
		start_time = time.time()
		if not (self.eqdsk_name.lower() == 'none'):
			self.load_eqdsk()
		self.init_current_iteration_eped()

		if (criterion == None):
			criterion = self.beta_critval
		
		ap0 = self.eped_prof[1,2]
		self.load_kin_profile(False)
		self.make_pressure_prof()
		self.eq.psin = np.copy(self.psin)
		self.eq.pp = np.copy(self.dpp)
		self.eq.pres = np.copy(self.dpres)
		self.current_iteration()
		
		self.load_eped_file = False
		
		if (self.beta_criterion == 1):
			err = abs(self.bp1c - criterion) /  criterion
			cval0 = self.bp1c
			print('Beta iteration #0 crit = %5.3f, val = %5.3f, atf = %f, li = %f, err = %9.6f'%(criterion,cval0,ap0,self.li,err))

		elif(self.beta_criterion == 2):
			err = abs(self.rav[0] - criterion) /  criterion
			cval0 = self.rav[0]
			print('Beta iteration #0 crit = %5.3f, val = %5.3f, atf = %f, li = %f, err = %9.6f'%(criterion,cval0,ap0,self.li,err))

		elif(self.beta_criterion == 3):
			err = abs(self.betan - criterion) /  criterion
			cval0 = self.betan
			print('Beta iteration #0 crit = %5.3f, val = %5.3f, atf = %f, li = %f, err = %9.6f'%(criterion,cval0,ap0,self.li,err))

		elif(self.beta_criterion == 0):
			err = 1.0
			cval0 = self.bp1c
			cval = self.ip1
			print('Beta = %5.3f, Ip = %5.3f'%(cval0,cval))

		niter = 0

		if (self.beta_criterion == 1 or self.beta_criterion == 3):
			ap1 = ap0 + ap0 * (criterion - cval0)/criterion * (self.p_fast + self.p_thermal) / self.p_thermal
		elif(self.beta_criterion == 2):
			ap1 = ap0 + ap0 * (criterion - cval0)/criterion * (self.p_fast + self.p_thermal) / self.p_thermal

		if (self.beta_criterion > 0):	
			self.eped_prof[1,2] = ap1
			self.eped_prof[2,2] = ap1
			self.load_kin_profile(False)
			self.make_pressure_prof()
			self.eq.psin = np.copy(self.psin)
			self.eq.pp = np.copy(self.dpp)
			self.eq.pres = np.copy(self.dpres)
			self.current_iteration()
		
		if (self.beta_criterion == 1):
			err = abs(self.bp1c - criterion) /  criterion
			cval = self.bp1c
			print('Beta iteration #1 crit = %5.3f, val = %5.3f, atf = %f, li = %f, err = %9.6f'%(criterion,cval,ap1,self.li,err))
		elif(self.beta_criterion == 2):
			err = abs(self.rav[0] - criterion) /  criterion
			cval = self.rav[0]
			print('Beta iteration #1 crit = %5.3f, val = %5.3f, atf = %f, li = %f, err = %9.6f'%(criterion,cval,ap1,self.li,err))
		elif(self.beta_criterion == 3):
			err = abs(self.betan - criterion) /  criterion
			cval = self.betan
			print('Beta iteration #1 crit = %5.3f, val = %5.3f, atf = %f, li = %f, err = %9.6f'%(criterion,cval,ap1,self.li,err))

		niter = niter + 1

		while (abs(err) > self.beta_crit and self.beta_criterion > 0 and niter < self.niterb):
			
			ap2 = ap1 + (criterion - cval) / (cval - cval0) * (ap1 - ap0)
			
			cval0 = cval
			ap0 = ap1
			
			self.eped_prof[1,2] = ap2
			self.eped_prof[2,2] = ap2
		
			self.load_kin_profile(False)	
			self.make_pressure_prof()
			self.eq.psin = np.copy(self.psin)
			self.eq.pp = np.copy(self.dpp)
			self.eq.pres = np.copy(self.dpres)
			self.current_iteration()
			
			if (self.beta_criterion == 1):
				err = abs(self.bp1c - criterion) /  criterion
				cval = self.bp1c
			elif(self.beta_criterion == 2):
				err = abs(self.rav[0] - criterion) /  criterion
				cval = self.rav[0]		
			elif(self.beta_criterion == 3):
				err = abs(self.betan - criterion) /  criterion
				cval = self.betan	
			
			ap1 = ap2
			niter = niter + 1
			print('Beta iteration #%2i crit = %f, val = %f, atf = %f, li = %f, err = %f'%(niter,criterion,cval,ap1,self.li, err))
			
		if not (self.nideal == 0):
			self.mapping_run()

		print('Write final kinetic profile')
		kinprof_type_old = self.kinprof_type
		self.kinprof_type = 1.
		self.write_kinprof()
		self.kinprof_type = kinprof_type_old
			
		print('Elapsed time = %f [s]'%(time.time()-start_time))

		return	

	def beta_iteration_eped2(self):
		#beta iteration by adjusting non-thermal + ext current amp iteration for matching li 
		start_time = time.time()

		niter = 0
		self.ajf = 0.	
		self.bjf = 0.9

		self.beta_iteration(None,False)
		li0 = self.li
		ierr = (li0 - self.li_target) / self.li_target		

		print('Li iteration #-1 crit = %f, val = %f, err = %f'%(self.li_target,li0,ierr))
		self.absolute_ajf = True

		if (self.cjf == 1.0 and self.djf == 1.0): #if cjf/djf preset
			if (ierr > 0.):
				self.cjf = 4.
				self.djf = 2.5
			else:
				self.cjf = 2.
				self.djf = 6.
			print('CJF, DJF changed to %f, %f'%(self.cjf,self.djf))

		self.ajf = max(0.1,abs(ierr)) * self.jcore[0]
		at0 = self.ajf

		self.load_eped_file = False
		self.beta_iteration(None,False)
		li0 = self.li
		err = (li0 - self.li_target) / self.li_target
		print('Li iteration #0 crit = %f, val = %f, ajf = %f, err = %f'%(self.li_target,li0,at0,err))

		self.ajf = 1.2 * self.ajf
		at1 = self.ajf
		
		if (abs(err) > self.li_crit and self.beta_criterion > 0 and niter < self.niterl):

			self.beta_iteration(None,False)
			li1 = self.li
			err = (li1 - self.li_target) / self.li_target
			print('Li iteration #1 crit = %f, val = %f, ajf = %f, err = %f'%(self.li_target,li0,at1,err))

		niter = niter + 1

		while (abs(err) > self.li_crit and self.beta_criterion > 0 and niter < self.niterl):
		
			if (ierr > 0.):
				at2 = at1 + (1./self.li_target - 1./li1) / (1./li1 - 1./li0) * (at1 - at0)
			else:
				at2 = at1 + (self.li_target - li1) / (li1 - li0) * (at1 - at0)
			
			li0 = li1
			at0 = at1
			self.ajf = at2
		
			self.beta_iteration(None,False)
			li1 = self.li
			err = (li1 - self.li_target) / self.li_target

			at1 = at2
			niter = niter + 1
			print('Li iteration #%2i crit = %f, val = %f, ajf = %f, err = %f'%(niter,self.li_target,li1,at1,err))
			
		if not (self.nideal == 0):
			self.mapping_run()
			
		print('Elapsed time = %f [s]'%(time.time()-start_time))

		return

	def beta_iteration_eped3(self):
		#beta iteration by adjusting non-thermal + ext current shape iteration for matching li
		start_time = time.time()

		niter = 0
		if(self.kinprof_type == 3):
			self.read_kinprof_eped()	
			self.bjf = 1.0
			if (self.cjf == 1.0 and self.djf == 1.0): 
				self.cjf = self.eped_prof[1,6]
				self.djf = self.eped_prof[1,7] * 1.5
				print('CJF, DJF changed to %f, %f'%(self.cjf,self.djf))
		else:
			self.bjf = 1.0
			self.cjf = 1.2
			self.djf = 3.5
			print('CJF, DJF changed to %f, %f'%(self.cjf,self.djf))

		self.beta_iteration(None,False)
		li0 = self.li
		ierr = (li0 - self.li_target) / self.li_target		

		if ierr > 0.:
				at0 = self.cjf
				print('Li iteration #0 crit = %f, val = %f, cjf = %f, djf = %f, err = %f'%(self.li_target,li0,self.cjf,self.djf,ierr))
		else:
				at0 = self.djf
				print('Li iteration #0 crit = %f, val = %f, cjf = %f, djf = %f, err = %f'%(self.li_target,li0,self.cjf,self.djf,ierr))
	
		if (ierr > 0.):
			at1 = self.cjf * 1.1
			self.cjf = at1
		else:
			at1 = self.djf * 1.1
			self.djf = at1

		self.load_eped_file = False
		self.beta_iteration(None,False)
		li1 = self.li
		err = (li1 - self.li_target) / self.li_target

		if ierr > 0.:
			print('Li iteration #1 crit = %f, val = %f, cjf = %f, djf = %f, err = %f'%(self.li_target,li1,at1,self.djf,err))
		else:
			print('Li iteration #1 crit = %f, val = %f, cjf = %f, djf = %f, err = %f'%(self.li_target,li1,self.cjf,at1,err))

		niter = niter + 1

		while (abs(err) > self.li_crit and self.beta_criterion > 0 and niter < self.niterl):
		
			if (ierr > 0.):
				at2 = at1 + (1./self.li_target - 1./li1) / (1./li1 - 1./li0) * (at1 - at0)
			else:
				at2 = at1 + (self.li_target - li1) / (li1 - li0) * (at1 - at0)
			
			li0 = li1
			at0 = at1

			if ierr > 0.:
				self.cjf = at2
			else:
				self.djf = at2

			self.beta_iteration(None,False)
			li1 = self.li
			err = (li1 - self.li_target) / self.li_target

			at1 = at2
			niter = niter + 1
			if ierr > 0.:
				print('Li iteration #%2i crit = %f, val = %f, cjf = %f, djf = %f, err = %f'%(niter,self.li_target,li1,at1,self.djf,err))
			else:
				print('Li iteration #%2i crit = %f, val = %f, cjf = %f, djf = %f, err = %f'%(niter,self.li_target,li1,self.cjf,at1,err))
		
		if not (self.nideal == 0):
			self.mapping_run()
			
		print('Elapsed time = %f [s]'%(time.time()-start_time))

		return		
		
	def beta_iteration_nubeam(self):
	
		start_time = time.time()

		try:
			os.mkdir('NUBEAM')
		except:
			pass
			
		try:
			os.mkdir('PROFILES')
		except:
			pass
		try:
			os.mkdir('OUTPUT')
		except:
			pass
					
		self.apf = 0.e0
		self.use_ext_current = False
		self.use_ext_pressure = False
#		self.adjust_prof = True
		self.beta_criterion = 0
		
		if not(self.eqdsk_name.lower() == 'none'):
			self.load_eqdsk()
		else:
			self.init_current_iteration_eped()
		
		if (os.path.isfile('PROFILES/VT_fit.dat')):
			self.vt_file = 'PROFILES/VT_fit.dat'
	
		if (self.kinprof_type == 1):
			kinprof_type = 1
			self.load_kin_profile(False)

		if (kinprof_type == 1):
			self.kinprof_type = 2
			self.write_kinprof()
			self.kinprof_type = kinprof_type
			copyfile('TI.dat_new','NUBEAM/TI.dat')
			copyfile('NE.dat_new','NUBEAM/NE.dat')
			copyfile('TE.dat_new','NUBEAM/TE.dat')
			copyfile('VT.dat_new','NUBEAM/VT.dat')

			move('TI.dat_new','PROFILES/TI.dat_0')
			move('NE.dat_new','PROFILES/NE.dat_0')
			move('TE.dat_new','PROFILES/TE.dat_0')
			move('VT.dat_new','PROFILES/VT.dat_0')

		else:

			copyfile(self.ti_file,'NUBEAM/TI.dat')
			copyfile(self.ne_file,'NUBEAM/NE.dat')
			copyfile(self.te_file,'NUBEAM/TE.dat')
			copyfile(self.vt_file,'NUBEAM/VT.dat')
		
			copyfile(self.ti_file,'PROFILES/TI.dat_0')
			copyfile(self.ne_file,'PROFILES/NE.dat_0')
			copyfile(self.te_file,'PROFILES/TE.dat_0')
			copyfile(self.vt_file,'PROFILES/VT.dat_0')

		copyfile(self.chease_rundir+'/EQDSK_COCOS_02.OUT','OUTPUT/geqdsk_0')
		copyfile(self.chease_rundir+'/chease_namelist','OUTPUT/chease_namelist_0')
		copyfile(self.chease_rundir+'/EXPEQ','OUTPUT/EXPEQ_0')
		copyfile(self.chease_rundir+'/EFIT_JCONST','OUTPUT/EFIT_JCONST_0')
		copyfile(self.chease_rundir+'/EFIT_JCONST','OUTPUT/Vneo.dat_0')
		copyfile(self.chease_rundir+'/EFIT_JCONST','OUTPUT/curr_prof_0')
		
		copyfile('nubeam_opt','NUBEAM/nubeam_opt')

		self.use_ext_current = True
		self.use_ext_pressure = True
		
		var_converge = np.zeros(shape=(self.nubeam_itern+1,6))
		var_converge[0,0] = self.bp
		var_converge[0,1] = self.li
		var_converge[0,2] = self.wmhd
		var_converge[0,3] = self.wdia
		var_converge[0,5] = self.eq.rmag

		self.pres_file = 'NUBEAM/chease_pres'
		self.curr_file = 'NUBEAM/chease_curr'

		for i in range(self.nubeam_itern):
		
			copyfile(self.chease_rundir+'/EQDSK_COCOS_02.OUT','NUBEAM/geqdsk')
			
			os.chdir('NUBEAM')
			os.system(python3_exec + ' ' + nubeam_dir)
			#print(python2_exec + ' ' + nubeam_exec)

			copyfile('chease_pres','../PROFILES/chease_pres_'+str(i+1))
			copyfile('chease_curr','../PROFILES/chease_curr_'+str(i+1))
			copyfile('nubeam_out0d','../PROFILES/nubeam_out0d_'+str(i+1))
			copyfile('nubeam_out1d','../PROFILES/nubeam_out1d_'+str(i+1))
			os.chdir('../')
			
			self.load_kin_profile(False)
			self.make_pressure_prof()
			self.eq.psin = np.copy(self.psin)
			self.eq.pp = np.copy(self.dpp)
			self.eq.pres = np.copy(self.dpres)
			self.current_iteration()
			
			print('Iteration #%2i, bp = %f, li = %f, wmhd = %f, wdia = %f, wfast = %f, rmag = %f \n'%(i,self.bp,self.li,self.wmhd,self.wdia,self.wfast,self.rav[0]))
			var_converge[i+1,0] = self.bp
			var_converge[i+1,1] = self.li
			var_converge[i+1,2] = self.wmhd
			var_converge[i+1,3] = self.wdia
			var_converge[i+1,4] = self.wfast
			var_converge[i+1,5] = self.rav[0]

			if (self.kinprof_type == 1):
				kinprof_type = 1
			elif (self.kinprof_type == 2):
				kinprof_type = 2
	
			self.kinprof_type = 2
			self.write_kinprof()
			self.kinprof_type = 1
			self.write_kinprof()
			self.kinprof_type = kinprof_type

			copyfile('NE.dat_new','PROFILES/NE.dat_'+str(i+1))
			copyfile('TE.dat_new','PROFILES/TE.dat_'+str(i+1))
			copyfile('TI.dat_new','PROFILES/TI.dat_'+str(i+1))
			copyfile('VT.dat_new','PROFILES/VT.dat_'+str(i+1))

			copyfile('chease_kinprof_new','PROFILES/chease_kinprof_'+str(i+1))
			copyfile('neo_coefs','PROFILES/neo_coefs_'+str(i+1))
			move('chease_kinprof_new','PROFILES/chease_kinprof_new')
				
			move('NE.dat_new','NUBEAM/NE.dat')
			move('TE.dat_new','NUBEAM/TE.dat')
			move('TI.dat_new','NUBEAM/TI.dat')
			move('VT.dat_new','NUBEAM/VT.dat')

			copyfile(self.chease_rundir+'/EQDSK_COCOS_02.OUT','OUTPUT/geqdsk_'+str(i+1))
			copyfile(self.chease_rundir+'/chease_namelist','OUTPUT/chease_namelist_'+str(i+1))
			copyfile(self.chease_rundir+'/EXPEQ','OUTPUT/EXPEQ_'+str(i+1))
			copyfile(self.chease_rundir+'/EFIT_JCONST','OUTPUT/EFIT_JCONST_'+str(i+1))
			copyfile(self.chease_rundir+'/curr_prof','OUTPUT/curr_prof_'+str(i+1))
			copyfile('Vneo.dat','OUTPUT/Vneo.dat_'+str(i+1))
			
		if not (self.nideal == 0):
			self.mapping_run()
			
		print('Elapsed time = %f [s]'%(time.time()-start_time))

		var_converge[0,4] = var_converge[1,4]

		f5 = open('OUTPUT/nubeam_iter_result','w')
		f5.write('iter#       bp              li           Wmhd[kJ]        Wdia[kJ]       Wfast[kJ]        Raxis[m]\n')
		for i in range(self.nubeam_itern+1):
			f5.write('%2i\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(i,var_converge[i,0],var_converge[i,1],var_converge[i,2],var_converge[i,3],var_converge[i,4],var_converge[i,5]))
		f5.close()

		os.chdir(self.currdir)

		f5 = open('OUTPUT/run_success','w')
		f5.write('1\n')
		f5.close()

		return
		
	def mapping_run(self):

		if (self.nideal == 8 or self.nideal == 11):
			self.eq.ns = self.ns #150
			self.eq.nt = self.nt #200
			self.eq.npsi = self.map_ns #300
			self.eq.nchi = self.map_nt #512
			self.eq.map_man = True
			self.eq.zeroedgep = True
			if not (self.ped_width == 0.05):
				self.eq.width = self.ped_width

		self.load_kin_profile(False)
		self.make_init_current(self.hager_mapf,self.VLOOP2)
		self.make_pressure_prof()
		self.eq.psin = np.copy(self.psin)
		self.eq.pp = np.copy(self.dpp)
		self.eq.pres = np.copy(self.dpres)
	
		kinprof_type_old  = self.kinprof_type
		self.kinprof_type = 1
		self.write_kinprof()
		self.kinprof_type = kinprof_type_old
		copyfile('chease_kinprof_new','CHEASE/chease_kinprof')
		self.run_chease(self.nideal,self.nomap)
		
		if not (self.nomap):
			print('--- Mapping procedure start with nideal = %i ---'%self.nideal)
			if (self.nideal == 10):
				move(self.chease_rundir+'/hamada.dat','hamada.dat')
			
			elif (self.nideal == 8):
				move(self.chease_rundir+'/NMISHKA','fort.12')
				move(self.chease_rundir+'/NELITE','eliteinp')
		else:
			print('--- Ignore mapping procedure ...')
			return		
		return	
				
	def nubeam_iter_plot(self,var_converge):
	
		fig, ([ax1,ax2,ax3],[ax4,ax5,ax6]) = plt.subplots(2,3)
		
		iter = np.linspace(0,self.nubeam_itern,self.nubeam_itern+1)

		ax1.set_title('Poloidal beta')
		ax1.set_xlabel('Iteration #')
		ax1.set_ylabel('betap [a.u]')
		ax1.scatter(iter,var_converge[:,0],c='r')
	
		ax2.set_title('Internal inductance')
		ax2.set_xlabel('Iteration #')
		ax2.set_ylabel('li [a.u]')
		ax2.scatter(iter,var_converge[:,1],c='r')

		ax3.set_title('WMHD')
		ax3.set_xlabel('Iteration #')
		ax3.set_ylabel('WMHD [kJ]')
		ax3.scatter(iter,var_converge[:,2],c='r')
		
		ax4.set_title('WDIA')
		ax4.set_xlabel('Iteration #')
		ax4.set_ylabel('WDIA [kJ]')
		ax4.scatter(iter,var_converge[:,3],c='g')
		
		ax5.set_title('WFAST')
		ax5.set_xlabel('Iteration #')
		ax5.set_ylabel('WFAST [kJ]')
		ax5.scatter(iter,var_converge[:,4],c='g')

		ax6.set_title('Magnetic axis')
		ax6.set_xlabel('Iteration #')
		ax6.set_ylabel('R [m]')
		ax6.scatter(iter,var_converge[:,5],c='g')
		
		plt.show(block=False)
		input()
		
	def bs_tor_calc(self):
			
		f4 = open(self.chease_rundir+'/JPAR','r')
		
		line= f4.readline()
		
		datnum= int(f4.readline())
		
		zjz = np.zeros(shape=(datnum,3))
		
		for i in range(datnum):
			
			line = f4.readline().split()
			
			zjz[i,0] = float(line[0])
			zjz[i,1] = abs(float(line[1]))
			zjz[i,2] = abs(float(line[2]))
			
		f4.close()
		
		self.zjz_parf = interp1d(zjz[:,0],zjz[:,1],'cubic')
		self.zjz_torf = interp1d(zjz[:,0],zjz[:,2],'cubic')
		self.par_over_tor = np.zeros(self.num)
		self.bs_tor = np.zeros(self.num)

		for i in range(self.num-1):
		
			self.par_over_tor[i] = abs(self.zjz_parf(self.psin[i]) / self.zjz_torf(self.psin[i]))

		self.par_over_tor[-1] = 2*self.par_over_tor[self.num-2] - self.par_over_tor[self.num-3]
	
		for i in range(self.num):
			
			if (self.psin[i] < self.par_over_tor_cut_psin):

				j= i
			else:
				self.par_over_tor[i] = self.par_over_tor[j]

		for i in range(self.num):

			self.bs_tor[i] = self.zjzbsc[i] / self.par_over_tor[i]
			
		self.bs_tot = np.trapz(self.bs_tor,x=self.area)
		if (self.use_ext_current):
			self.ex_tot = np.trapz(self.curr_ex,x=self.area)

		return
		
	def make_efit_constraint(self):
	
		self.load_kin_profile(False)
		self.make_pressure_prof()
		if (self.simple_run):
			#self.load_eqdsk()
			#self.read_hagermap(self.hager_mapf)
			self.make_init_current(self.hager_mapf)
		
		try:
			os.mkdir('OUTPUT')
		except:
			pass

		jtor, ffp, self.vloop2 = self.make_current2()
		self.bs_tor_calc()
		copyfile('Vneo.dat','OUTPUT/Vneo.dat')
		f4 = open('OUTPUT/EFIT_JCONST','w')
		f5 = open('OUTPUT/BS_Profile','w')
		f4.write('%9s\t%9s\t%9s\t%9s\n'%('PSIN[a.u]','JTOR[A/m^2]','JTOR_NORM','FFp[T]  '))
		f5.write('%9s\t%9s\t%9s\n'%('PSIN[a.u]','JBS_PAR[A/m^2]','JBS_TOR[A/m^2]'))
		f4.write('%i\n'%len(self.psin))
		for i in range(len(self.psin)):
			
			val = jtor[i]
			val_norm = val / (abs(self.ip) / self.area[-1])

			f4.write('%9.6f\t%9.6e\t%9.6f\t%9.6f\n'%(self.psin[i],val,val_norm,ffp[i]))
			f5.write('%9.6f\t%9.6e\t%9.6e\n'%(self.psin[i],self.zjzbsc[i],self.bs_tor[i]))

		f5.write('TOT BSCURRENT \t %f [A] \n'%self.bs_tot)
		f5.write('TOT CURRENT \t %f [A] \n'%abs(self.ip))
		f4.close()
		f5.close()
		print('IP = %f, BS = %f [kA]'%(abs(self.ip)/1.e3,self.bs_tot/1.e3))
		if self.use_ext_current:
			print('IPEXT = %f, IPIND = %f [kA]'%(self.ex_tot/1.e3,(abs(self.ip)-self.bs_tot-self.ex_tot)/1e3))
		print('EFIT current profile constraint done')
		return		

	def make_efit_constraint2(self):

		models = ['sauter','csauter','hager','neo']

		self.load_kin_profile(False)
		self.make_pressure_prof()
	
		self.use_hager = False
		self.use_chang = False
		self.use_neo   = False

		self.bs2kstar = dict()
		len2 = len(self.psin)
		self.bs2kstar['opt'] = dict()
		self.bs2kstar['opt']['sauter']  = [False, False, False]
		self.bs2kstar['opt']['csauter'] = [True,  False, False]
		self.bs2kstar['opt']['hager']   = [True,  True,  False]
		self.bs2kstar['opt']['neo']     = [False,  False,  True]

		self.bs2kstar['bstot'] = dict()
		self.bs2kstar['bspar'] = dict()
		self.bs2kstar['bstor'] = dict()
		self.bs2kstar['ffp']   = dict()
		self.bs2kstar['jtor']  = dict()

		self.bs2kstar['psin'] = np.copy(self.psin)
		slegend = [];	plegend = [];
		for i in ['sauter','csauter','hager','neo']:

			print('>>> For the BS model %s...'%i)
			self.use_chang = self.bs2kstar['opt'][i][0]
			self.use_hager = self.bs2kstar['opt'][i][1]
			self.use_neo   = self.bs2kstar['opt'][i][2]

			jtor, ffp, self.vloop2 = self.make_current2()
			self.bs_tor_calc()

			self.bs2kstar['bstot'][i] = self.bs_tot
			self.bs2kstar['bspar'][i] = np.copy(self.zjzbsc)
			self.bs2kstar['bstor'][i] = np.copy(self.bs_tor)
			self.bs2kstar['ffp'][i]   = np.copy(ffp)
			self.bs2kstar['jtor'][i]  = np.copy(jtor)
			
			slegend.append(i)
		fig,[ax1, ax2] = plt.subplots(1,2,figsize=(12,5))
		for i in ['sauter','csauter','hager','neo']:
			ax1.plot(self.psin,self.bs2kstar['bstor'][i]/1.e6)
			ax2.plot(self.psin,self.bs2kstar['jtor'][i]/1.e6)

		print('>>> ======================================== <<<')
		
		print('%10s\t%9s\t%9s'%('Model','BSTOT[MA]','Fraction[%]'))
		for i in ['sauter','csauter','hager','neo']: print('%10s\t%9.7f\t%9.7f'%(i.upper(),self.bs2kstar['bstot'][i]/1.e6,self.bs2kstar['bstot'][i]/self.ip*-1))
		print('>>> ======================================== <<<')
		ax1.set_xlabel('$\psi_N$')
		ax1.set_ylabel('B.S current density [MA/$m^2$]')
		ax1.legend(slegend)
		ax2.set_xlabel('$\psi_N$')
		ax2.set_ylabel('Current density [MA/$m^2$]')
		ax2.legend(slegend)
		plt.tight_layout()
		plt.show(block=False)
		count = 1
		while count == 1:
			model_type = input('>>> Choose desired mode (1:sauter, 2:csauter, 3:hager, 4:neo)\n>>> ')
			count = 0
			try:	model_type = models[int(model_type)-1]
			except:	count = 1

		cur_norm = (abs(self.ip) / self.area[-1])
		f = open('EFIT_JCONST.dat','w')
		f.write('Model type  = %s\n'%model_type)
		f.write('%13s\t%13s\t%13s\t%13s\n'%('psi_norm','rho_norm','JTOR[MA/m2]','JTOR_norm'))
		for i in range(len(self.psin)):	
			f.write('%13.8f\t%13.8f\t%13.8f\t%13.8f\n'%(self.psin[i],self.rho[i],self.bs2kstar['jtor'][model_type][i]/1.e6,self.bs2kstar['jtor'][model_type][i]/cur_norm))
		f.close()
		f = open('BS_PROF.dat','w')
		f.write('Model type  = %s\n'%model_type)
		f.write('%13s\t%13s\t%13s\t%13s\n'%('psi_norm','rho_norm','BSPAR[MA/m2]','BSTOR[MA/m2]'))		
		for i in range(len(self.psin)):
			f.write('%13.8f\t%13.8f\t%13.8f\t%13.8f\n'%(self.psin[i],self.rho[i],self.bs2kstar['bspar'][model_type][i]/1.e6,self.bs2kstar['bstor'][model_type][i]/1.e6))
				
		f.close()
		
		print('>>> BS_PROF.DAT (bootstrap current profile), EFIT_JCONST.DAT (EFIT input) are saved!')

		return		

	def hager_const(self):
	
		self.mi_hag = self.amain * 1.672*1.e-27
		self.me_hag = 9.109 * 1.e-31
		self.a1_hag = 0.0488;
		self.a2_hag = 0.0488;
		self.b1_hag = 0.0
		self.b2_hag = 0.0086
		self.c2_hag = 0.2947
		self.epsc1_hag = 0.15
		self.epsc2_hag = 0.5099
		self.w1_hag = 0.1;
		self.w2_hag = 0.15;
		self.lam1_hag = 19.1702;
		self.lam2_hag = 1.9056;
		self.w2_hag = 0.15;
		self.lam1_hag = 19.1702;
		self.lam2_hag = 1.9056;
		self.lam3_hag = 0.0106;
		self.lam4_hag = 2.994;
		self.lam5_hag = 0.9958;
		self.ci1_hag = 1.357;
		self.ci2_hag = 2.192;
		self.alp_hag = 2.1;
		self.alph_hag = self.ci2_hag/self.ci1_hag;
		
		return
		
	def initialise_var(self):
	
		self.beta_criterion = 0
		self.beta_crit = 1.e-2
		self.li_crit = 2.e-2
		self.beta_critval = np.nan
		self.ip_crit = 1.e-5
		self.bs_crit = 1.e-5
		self.kinprof_type = 0
		self.te_file = 'chease_te'
		self.ti_file = 'chease_ti'
		self.ne_file = 'chease_ne'
		self.ni_file = None
		self.vt_file = None
		self.chkin_file = 'chease_kinprof'
		self.eped_file = 'chease_eped'
		self.bnd_file = 'BND_RZ'
		self.load_eped_file = True
		self.hager_mapf = 'hager_map'
		self.chease_rundir = 'CHEASE'
		self.eqdsk_name = 'geqdsk'
		self.use_rho = False

		self.pres_file = 'chease_pres'
		self.curr_file = 'chease_curr'
		self.isvt = False		
	
		self.target_psin = 0.995
		self.bsmulti = 1.0
		self.ped_width = 0.05
		self.pped    = 0.
			
		self.eped_first_run = False
		self.use_param_shape = False
		self.elong = np.nan
		self.triang = np.nan
		self.square = np.nan
		self.amin = np.nan
		self.rcent = np.nan
		self.ip = np.nan
		self.bcent = np.nan
		
		self.grid_xr1 = 1.0
		self.grid_sig1 = 0.01
		self.grid_xr2 = 1.0 - 0.5*self.ped_width
		self.grid_sig2 = self.ped_width
		self.grid_n = 301
		
		self.use_hager = True
		self.use_chang = True
		self.use_neo   = True
		self.use_chease_jprof = True
		
		self.use_beam_pressure = True
		self.use_ext_pressure = False
		self.apf = 0.5
		self.bpf = 0.8
		self.cpf = 2.0
		self.dpf = 2.0
		
		self.use_beam_current = True
		self.use_ext_current = False
		self.ajf = 0.0
		self.bjf = 0.7
		self.cjf = 1.0
		self.djf = 1.0
		
		self.zmain = 1.0
		self.zeff = 2.0
		self.zimp = 6.0
		self.amain = 2.0
		self.aimp = 12.0
		
		self.niterc = 10
		self.niterb = 10
		self.niterl = 10
		
		self.vloop_ext = 0.0;
		self.vloop_core_mod = 1.0;
		self.core_neo = 0.001;

		self.hag_mod = True
		self.hag_mod_psin = 0.3
		self.hag_mod_psinw = 0.15
		self.par_over_tor_cut_psin = 1.0
		self.betaw_low_limit = 0.4

		self.pedscan_ti = False

		self.remove_temp = False
		
		self.run_mode = 'normal'
		self.map_mode = 'none'
		self.use_eped = False
		self.use_eped2 = False
		self.use_eped3 = False
		self.simple_run = False
		self.nubeam_run = False

		self.nubeam_itern = 6
		
		self.astra_const = False
		self.efit_const = False
		self.efit_const2 = False
		self.adjust_prof = False

		self.use_scaled_density = False
		self.scaled_density = 1.0
		self.target_density = 0.0
		
		self.nideal = 0
		self.nomap = False

		self.zjz = np.zeros(2)
		self.zjz_old = np.zeros(2)
		self.relax = 0.75

		self.ns = 150
		self.nt = 200
		self.map_ns = 300
		self.map_nt = 512
		
		if (self.kinprof_type == 3):
			self.adjust_prof = False

#-- ext
		self.wdia = 0.0
		self.device_name = 'KSTAR'
		self.epslon = 8

		self.given_j = False
		self.zjzg = np.zeros(1)

		self.li_target = 0.
		self.absolute_ajf = False
		
		self.beta = 0.
		self.betan = 0.

		return
		
	def load_eqdsk(self,nideal=11,chease_skip=False):

		start_time = time.time()	
		self.eq = eqdsk.eqdsk(self.eqdsk_name)
		self.eq.nideal = nideal
		self.eq.device_name = self.device_name
		if (self.target_psin == 1.0):
			self.eq.target_psin = 0.0
		else:
			self.eq.target_psin = self.target_psin
		self.eq.read_eqdsk_file()
		
		try:
			os.mkdir(self.chease_rundir)
		except:
			pass
		print('Initial CHEASE run with EQDSK')
		
		self.eq.make_chease_input(self.chease_rundir+'/EXPEQ',self.chease_rundir+'/chease_namelist')
		os.chdir(self.chease_rundir)
		if not chease_skip:	
			stat,out = subprocess.getstatusoutput(chease_exec +' > log.chease')
		print('Run finished, Elapsed time=%f [s]'%(time.time()-start_time))
		
		f4 = open('log.chease','r')

		linec = 0
	
		while True:

			line = f4.readline()
			if not line: break

			if (line.find('MKSA') > -1):
				linec = 1

			if ((line.find('TOTAL CURRENT -->') > -1) and (linec == 1)):
				ipc = float(line.split()[5])
			if ((line.find('POLOIDAL BETA') > -1) and (linec == 1)):
				self.bp =float(line.split()[0])
			if ((line.find('WMHD') > -1) and (linec == 1)):
				self.wmhd = float(line.split()[0])
			if ((line.find('LI') > -1) and (linec == 1)):
				self.li = float(line.split()[0])

		self.read_rho_psi_R()
		
		self.R_mapi = np.copy(self.R_map)
		self.psi_mapi = np.copy(self.psi_map)
		self.rho_mapi = np.copy(self.rho_map)
		self.bcent = self.eq.bcentr
		
		if (self.beta_criterion == 2):
			if (self.beta_critval == 0):
					self.beta_critval = self.eq.rmag
		if (self.beta_criterion == 1):
			if (self.beta_critval == 0):
					self.beta_critval = self.bp
		if (self.beta_criterion == 3):
			if (self.beta_critval == 0):
					self.beta_critval = self.betan

		if (self.li_target == 0):
			self.li_target = self.li
		
		os.chdir(self.currdir)
			
		
		self.ip = self.eq.ip

		return
		
	def read_namelist_str(self,line,var_name,vars,vartype):

		line = line.split('!')[0]
		line = line.split('#')[0]
		line = line.split('\n')[0]	
		try:
			name = line.split('=')[0].split()[0].lower()
		except:
			name = None
		
		if (name == var_name.lower()):
			
			if (vartype < 4):
				if (vartype == 1):
					var = int(float(line.split('=')[1].replace(' ','')))
				elif (vartype == 2):
					var = float(line.split('=')[1].replace(' ',''))
				elif (vartype == 3):
					var = line.split('=')[1].replace(' ','')
					
				exist = 1
			else:
				var2 = line.split('=')[1].replace(' ','').lower()
				
				if (var2 == 'true'):
					var = True
					
				elif (var2 == 'false'):
					var = False
				else:
					var = vars
					
		else:
		
			var = vars;
		line2 = line.split('!')[0].split()
		if(line2 == []):
			var = vars;

		if (var == ''):
			var = None
		
		return (var)
		
	def read_namelist(self,filename):
	
		f4 = open(filename,'r')
		
		while True:
		
			line = f4.readline()
			if not line: break
	
			self.beta_criterion = self.read_namelist_str(line,'Beta_criterion_type',self.beta_criterion,1)
			self.beta_crit = self.read_namelist_str(line,'Beta_crit',self.beta_crit,2)
			self.beta_critval = self.read_namelist_str(line,'Beta_val',self.beta_critval,2)
			self.ip_crit = self.read_namelist_str(line,'Ip_crit',self.ip_crit,2)
			self.bs_crit = self.read_namelist_str(line,'Bs_crit',self.bs_crit,2)
			self.li_crit = self.read_namelist_str(line,'Li_crit',self.li_crit,2)
			self.kinprof_type = self.read_namelist_str(line,'kinetic_profile_type',self.kinprof_type,1)
			self.ne_file = self.read_namelist_str(line,'ne_file',self.ne_file,3)
			self.ni_file = self.read_namelist_str(line,'ni_file',self.ni_file,3)
			self.te_file = self.read_namelist_str(line,'te_file',self.te_file,3)
			self.ti_file = self.read_namelist_str(line,'ti_file',self.ti_file,3)
			self.vt_file = self.read_namelist_str(line,'vt_file',self.vt_file,3)
			
			self.chkin_file = self.read_namelist_str(line,'chease_kinetic_file',self.chkin_file,3)
			self.eped_file = self.read_namelist_str(line,'eped_kinetic_file',self.eped_file,3)
			self.hager_mapf = self.read_namelist_str(line,'hager_map_file',self.hager_mapf,3)
			self.chease_rundir = self.read_namelist_str(line,'chease_run_dir',self.chease_rundir,3)
			self.target_psin = self.read_namelist_str(line,'BND_PSIN',self.target_psin,2)
			self.eqdsk_name = self.read_namelist_str(line,'EQDSK',self.eqdsk_name,3)
			self.pres_file = self.read_namelist_str(line,'PRES_FILE',self.pres_file,3)
			self.curr_file = self.read_namelist_str(line,'CURR_FILE',self.curr_file,3)

			self.bnd_file = self.read_namelist_str(line,'BND_FILE',self.bnd_file,3)
			self.use_param_shape = self.read_namelist_str(line,'USE_PARAM_SHAPE',self.use_param_shape,4)
			
			self.bsmulti = self.read_namelist_str(line,'BSMULTI',self.bsmulti,2)
			self.ped_width = self.read_namelist_str(line,'PED_WIDTH',self.ped_width,2)
			self.elong = self.read_namelist_str(line,'ELONG',self.elong,2)
			self.triang = self.read_namelist_str(line,'TRIA',self.triang,2)
			self.square = self.read_namelist_str(line,'SQUARE',self.square,2)
			self.amin = self.read_namelist_str(line,'AMINOR',self.amin,2)
			self.rcent = self.read_namelist_str(line,'RMAJOR',self.rcent,2)
			self.ip =  self.read_namelist_str(line,'IP',self.ip,2)
			self.bcent = self.read_namelist_str(line,'BVAC',self.bcent,2)
			
			self.grid_xr1 = self.read_namelist_str(line,'XR1',self.grid_xr1,2)
			self.grid_sig1 = self.read_namelist_str(line,'SIG1',self.grid_sig1,2)
			self.grid_n = self.read_namelist_str(line,'GRIDN',self.grid_n,1)
			self.use_hager = self.read_namelist_str(line,'USE_HAGER',self.use_hager,4)
			self.use_neo   = self.read_namelist_str(line,'USE_NEO',self.use_neo,4)
			self.use_chang = self.read_namelist_str(line,'USE_CHANG',self.use_chang,4)
			self.use_beam_pressure = self.read_namelist_str(line,'USE_FAST_P',self.use_beam_pressure,4)
			self.use_ext_pressure = self.read_namelist_str(line,'USE_EXT_P',self.use_ext_pressure,4)
			self.use_beam_current = self.read_namelist_str(line,'USE_FAST_J',self.use_beam_current,4)
			self.use_ext_current = self.read_namelist_str(line,'USE_EXT_J',self.use_ext_current,4)
			
			self.apf = self.read_namelist_str(line,'APF',self.apf,2)
			self.bpf = self.read_namelist_str(line,'BPF',self.bpf,2)
			self.cpf = self.read_namelist_str(line,'CPF',self.cpf,2)
			self.dpf = self.read_namelist_str(line,'DPF',self.dpf,2)
			
			self.ajf = self.read_namelist_str(line,'AJF',self.ajf,2)
			self.bjf = self.read_namelist_str(line,'BJF',self.bjf,2)
			self.cjf = self.read_namelist_str(line,'CJF',self.cjf,2)
			self.djf = self.read_namelist_str(line,'DJF',self.djf,2)			
			
			self.zmain = self.read_namelist_str(line,'ZMAIN',self.zmain,1)
			self.zeff = self.read_namelist_str(line,'ZEFF',self.zeff,2)
			self.zimp = self.read_namelist_str(line,'ZIMP',self.zimp,1)
			self.amain = self.read_namelist_str(line,'AMAIN',self.amain,1)
			self.aimp = self.read_namelist_str(line,'AIMP',self.aimp,1)
			
			self.pedscan_ti = self.read_namelist_str(line,'PEDSCAN_TI',self.pedscan_ti,4)	

			self.niterc = self.read_namelist_str(line,'Current_ITERN',self.niterc,1)
			self.niterl = self.read_namelist_str(line,'Li_ITERN',self.niterl,1)
			
			self.vloop_ext = self.read_namelist_str(line,'EXT_VLOOP',self.vloop_ext,2)
			self.vloop_core_mod = self.read_namelist_str(line,'VLOOP_MOD',self.vloop_core_mod,2)
			
			self.core_neo = self.read_namelist_str(line,'CORE_NEO',self.core_neo,2)

			self.hag_mod = self.read_namelist_str(line,'HAG_CORE_MOD',self.hag_mod,4)
			self.hag_mod_psin = self.read_namelist_str(line,'HAG_CORE_MOD_PSIN',self.hag_mod_psin,2)

			self.remove_temp = self.read_namelist_str(line,'REMOVE_TEMP',self.remove_temp,4)
			self.adjust_prof = self.read_namelist_str(line,'Adjust_prof',self.adjust_prof,4)

			self.use_eped = self.read_namelist_str(line,'USE_EPED',self.use_eped,4)
			self.use_eped2 = self.read_namelist_str(line,'USE_EPED2',self.use_eped2,4)
			self.use_eped3 = self.read_namelist_str(line,'USE_EPED3',self.use_eped3,4)
			
			self.use_rho = self.read_namelist_str(line,'USE_RHO',self.use_rho,4)
			self.astra_const = self.read_namelist_str(line,'ASTRA_CONST',self.astra_const,4)
			self.efit_const = self.read_namelist_str(line,'EFIT_CONST',self.efit_const,4)
			self.efit_const2 = self.read_namelist_str(line,'EFIT_CONST2',self.efit_const2,4)
			self.simple_run = self.read_namelist_str(line,'SIMPLE_RUN',self.simple_run,4)
			self.run_mode = self.read_namelist_str(line,'RUN_MODE',self.run_mode,3)

			self.map_mode = self.read_namelist_str(line,'MAP_MODE',self.map_mode,3)

			self.target_density = self.read_namelist_str(line,'DENSITY_SCALE',self.target_density,2)

			self.nubeam_itern = self.read_namelist_str(line,'NUBEAM_ITERN',self.nubeam_itern,1)

			self.betaw_low_limit = self.read_namelist_str(line,'HAGER_LIMIT',self.betaw_low_limit,2)
			
			self.nideal = self.read_namelist_str(line,'NIDEAL',self.nideal,1)
			self.epslon = self.read_namelist_str(line,'EPSLON',self.epslon,1)

			self.ns = self.read_namelist_str(line,'NS',self.ns,1)
			self.nt = self.read_namelist_str(line,'NT',self.nt,1)
			self.map_ns = self.read_namelist_str(line,'MAP_NS',self.map_ns,1)
			self.map_nt = self.read_namelist_str(line,'MAP_NT',self.map_nt,1)

			self.wdia = self.read_namelist_str(line,'WDIA',self.wdia,2)
			self.relax = self.read_namelist_str(line,'RELAX',self.relax,2)

			self.device_name = self.read_namelist_str(line,'DEVICE',self.device_name,3)
		
			self.given_j = self.read_namelist_str(line,'GIVENJ',self.given_j,4)
			self.li_target = self.read_namelist_str(line,'LITARGET',self.li_target,2)

			self.adjust_variables()
		
		return

	def adjust_variables(self):

		if (self.target_density > 0.0):
			self.use_scaled_density = True	
	
		if (self.run_mode.lower() == 'eped'):
			self.use_eped = True
			self.simple_run = False
		elif (self.run_mode.lower() == 'eped2'):
			self.use_eped = False
			self.use_eped2 = True
			self.use_eped3 = False
			self.simple_run = False
			self.ajf = 1.0
		elif (self.run_mode.lower() == 'eped3'):
			self.use_eped = False
			self.use_eped2 = False
			self.use_eped3 = True
			self.simple_run = False
			self.ajf = 1.0			
		elif (self.run_mode.lower() == 'simple'):
			self.use_eped = False
			self.simple_run = True
		elif (self.run_mode.lower() == 'nubeam'):
			self.use_eped = False
			self.simple_run = False
			self.nubeam_run = True
		else:
			self.use_eped = False
			self.simple_run = False

		if (self.relax > 1.0 or self.relax <=0.0):
			print('Relax should be 0 < relax < 1')
			self.relax = 0.7

		if (self.nomap):
			self.epslon = 5

		if not (self.nomap):
			if self.nideal == 8:
				if self.epslon < 8:
					self.epslon = 8

		return		

		
	def __init__(self,filename='chease_opt'):
		
		self.initialise_var()
		try:
			self.read_namelist(filename)
		except:
			print('No CHEASE extra namelist')
		#	exit()

		self.psin_grid(self.grid_xr1,self.grid_sig1,self.grid_xr1,self.grid_sig1,self.grid_n)

		self.currdir = os.getcwd()

		self.chease_rundir = self.currdir + '/' + self.chease_rundir	
		self.hager_mapf = self.chease_rundir + '/' + self.hager_mapf
		
		return
