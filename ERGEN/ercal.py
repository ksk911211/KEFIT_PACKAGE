#!/usr/bin/env python
## Concept is developed by Y.H.Lee
import os, sys
import numpy as np
import time
from scipy.interpolate import interp1d, interp2d
import subprocess
import matplotlib.pyplot as plt
import eqdsk
from exec_dirs import chease_exec

currdir = os.getcwd()
e0 = 1.601 * 1.e-19;
mu0 = 4.0 * np.pi * 1.e-7;

class Er_profile:

	def read_kinprof_chease(self):
	
		f = open(self.f_profile,'r')
		self.dat_numk = int(f.readline().split()[0])
		line = f.readline().split()
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
		self.isvt = True	
		for i in range(self.dat_numk):
		
			line = f.readline().split()
				
			self.psink[i] = float(line[0])	
			self.tek[i] = float(line[1])
			self.nek[i] = float(line[2])
			self.tik[i] = float(line[3])
			self.nik[i] = float(line[4])
			try: self.vtk[i] = float(line[5])
			except: self.isvt = False

		f.close()
		
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

	def interpol(self):

		nef = interp1d(self.psink,self.nek,'cubic')
		tef = interp1d(self.psink,self.tek,'cubic')
		nif = interp1d(self.psink,self.nik,'cubic')
		tif = interp1d(self.psink,self.tik,'cubic')

		nef2 = interp1d(self.psink[0:3],self.nek[0:3],'quadratic')
		tef2 = interp1d(self.psink[0:3],self.tek[0:3],'quadratic')
		nif2 = interp1d(self.psink[0:3],self.nik[0:3],'quadratic')
		tif2 = interp1d(self.psink[0:3],self.tik[0:3],'quadratic')

		if self.isvt:
			vtf = interp1d(self.psink,self.vtk,'cubic')
		else:
			vtf = interp1d(self.psink,self.psink)

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
	
		self.ne = np.zeros(self.num)
		self.te = np.zeros(self.num)
		self.ni = np.zeros(self.num)
		self.ti = np.zeros(self.num)
		self.vt = np.zeros(self.num)
		self.rav = np.zeros(self.num)
		self.b02av = np.zeros(self.num)
		self.bgpb = np.zeros(self.num)
		self.ft = np.zeros(self.num)
		self.rc = np.zeros(self.num)
		self.q = np.zeros(self.num)
		self.f = np.zeros(self.num)
		self.r = np.zeros(self.num)
		self.pe = np.zeros(self.num)
		self.pi = np.zeros(self.num)
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
			self.pi[i] = self.ni[i] * self.ti[i] * 1.e19 * 1.e3 * e0
			self.pt[i] = (self.ne[i]*self.te[i] + self.ni[i]*self.ti[i]) * 1.e19 * 1.e3 * e0
	
		self.BMAG = self.f[0] / self.rc[0]
				
		self.dne = np.zeros(self.num)
		self.dte = np.zeros(self.num)
		self.dni = np.zeros(self.num)
		self.dti = np.zeros(self.num)
		self.drc = np.zeros(self.num)
		self.dq = np.zeros(self.num)
		self.dpe = np.zeros(self.num)
		self.dpi = np.zeros(self.num)
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
		
		self.dne[self.num-1] = 2.0*self.dne[self.num-2] - self.dne[self.num-3]
		self.dte[self.num-1] = 2.0*self.dte[self.num-2] - self.dte[self.num-3]
		self.dni[self.num-1] = 2.0*self.dni[self.num-2] - self.dni[self.num-3]
		self.dti[self.num-1] = 2.0*self.dti[self.num-2] - self.dti[self.num-3]
		self.drc[self.num-1] = 2.0*self.drc[self.num-2] - self.drc[self.num-3]
		self.dq[self.num-1] = 2.0*self.dq[self.num-2] - self.dq[self.num-3]

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
			self.dpi[i] = (self.dni[i]*self.ti[i]+self.ni[i]*self.dti[i]) * 1.e19 * 1.e3 * e0
			self.dpt[i] = (self.dni[i]*self.ti[i]+self.ni[i]*self.dti[i]) * 1.e19 * 1.e3 * e0 + self.dpe[i]
			
		self.dp_ex[0] 	 = 2.0*self.dp_ex[1] - self.dp_ex[2]
		self.dp_ex[self.num-1] = 2.0*self.dp_ex[self.num-2] - self.dp_ex[self.num-3]

		if self.dp_ex[0] * self.psia > 0.: self.dp_ex[0] = 0.

		return

	def load_eqdsk(self):

		start_time = time.time()	
		self.eq = eqdsk.eqdsk(self.f_eqdsk)
		self.eq.nideal = 11
		self.eq.device_name = self.device_name
		self.eq.target_psin = 0.995
		self.eq.read_eqdsk_file()		
		try: os.mkdir(self.chease_rundir)
		except: pass
		print('Initial CHEASE run with EQDSK')
		
		self.eq.make_chease_input(self.chease_rundir+'/EXPEQ',self.chease_rundir+'/chease_namelist')
		os.chdir(self.chease_rundir)
		stat,out = subprocess.getstatusoutput(chease_exec +' > log.chease')
		print('Run finished, Elapsed time=%f [s]'%(time.time()-start_time))
		
		f = open('log.chease','r')

		linec = 0	
		while True:
			line = f.readline()
			if not line: break

			if (line.find('MKSA') > -1):
				linec = 1

			if ((line.find('TOTAL CURRENT -->') > -1) and (linec == 1)):
				ipc = float(line.split()[5])
			if ((line.find('POLOIDAL BETA') > -1) and (linec == 1)):
				self.eq.bp =float(line.split()[0])
			if ((line.find('WMHD') > -1) and (linec == 1)):
				self.wmhd = float(line.split()[0])
			if ((line.find('LI') > -1) and (linec == 1)):
				self.eq.li = float(line.split()[0])

		os.chdir(self.currdir)
		self.read_rho_psi_R()
		
		self.R_mapi = np.copy(self.R_map)
		self.psi_mapi = np.copy(self.psi_map)
		self.rho_mapi = np.copy(self.rho_map)
		self.bcent = self.eq.bcentr
		self.ip = self.eq.ip
		
		return

	def read_rho_psi_R(self):
	
		currdir2 = os.getcwd()
		f = open(self.f_rhopsir,'r')
		line = f.readline()
		line = f.readline()
		self.num_map = int(line)
		
		self.rho_map = np.zeros(self.num_map)
		self.psi_map = np.zeros(self.num_map)
		self.R_map = np.zeros(self.num_map)
		
		for i in range(self.num_map):
		
			line = f.readline().split()
			self.rho_map[i] = float(line[0])
			self.psi_map[i] = float(line[1])
			self.R_map[i] = float(line[2])	
			
		f.close()

		rhop = interp1d(self.psi_map,self.rho_map,'cubic')
		Rp = interp1d(self.psi_map,self.R_map,'cubic')
		self.rho = rhop(self.psin)
		self.RR  = Rp(self.psin)

		return		

	def read_hagermap(self):
	
		f = open(self.f_hager,'r')
		line = f.readline()	
		self.dat_numh = int(f.readline().split()[0])
		
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
		
			line = f.readline().split()
			
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
			
		
		line = f.readline().split()
		
		self.psia = abs(float(line[0]))
		
		f.close()
		
		return

	def hager_bs(self):

		f = open('hager_coef','w')
	
		self.hager_const()
		self.alps = np.zeros(self.num)

		if (self.use_neo): self.use_hager = True	
		if (self.use_hager): self.use_chang = False
		
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
			Z2_hag = self.zeff
			Z_hag = np.power((self.zmain**2.0)*Z1_hag*Z2_hag,0.25)
			self.Z = Z2_hag
			
			ne_hag = self.ne[i] * 1.e19
			ni_hag = self.ni[i] * 1.e19
			te_hag = self.te[i] * 1.e3
			ti_hag = self.ti[i] * 1.e3
			rpe_hag = ne_hag*te_hag / (ni_hag*ti_hag+ne_hag*te_hag)
			
			vti_hag = np.sqrt(2.0*ti_hag*e0/self.mi_hag)

			if (self.use_neo):
				DPSI_hag = 0.
			else:
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
			
			if (self.use_chang):
				f32_hag = f32_hag * delta
				f33_hag = f33_hag * delta
				f34_hag = f34_hag * delta

			alp0_hag = -1.17*(1.-self.ft[i])/(1.-0.22*self.ft[i]-0.19*(self.ft[i]**2))
        
			alps_hag = ((alp0_hag + 0.25*(1-(self.ft[i]**2))*np.sqrt(nuis_hag))/ (1.+0.5*np.sqrt(nuis_hag)) + 0.315*(nuis_hag**2)*(self.ft[i]**6)) / (1.+0.15*(nuis_hag**2)*(self.ft[i]**6))

			l31s_hag = (1.+1.4/(Z2_hag+1.))*f31_hag - 1.9/(Z2_hag+1.)* (f31_hag**2) + 0.3/(Z2_hag+1.)*(f31_hag**3) + 0.2/(Z2_hag+1.)*(f31_hag**4)
			
			if (self.use_chang):
				l32s_hag = (0.05+0.61*Z2_hag)/Z2_hag/(1.+0.44*Z2_hag)*(f32_hag - (f32_hag**4)) + 1./(1.+0.22*Z2_hag)*((f32_hag**2)-(f32_hag**4) -1.2*(f32_hag**3)+1.2*(f32_hag**4)) + 1.2/(1.+0.5*Z2_hag)* (f32_hag**4)
			else:
				l32s_hag = (0.05+0.62*Z2_hag)/Z2_hag/(1.+0.44*Z2_hag)*(f32_hag - (f32_hag**4)) + 1./(1.+0.22*Z2_hag)*((f32_hag**2)-(f32_hag**4) -1.2*(f32_hag**3)+1.2*(f32_hag**4)) + 1.2/(1.+0.5*Z2_hag)* (f32_hag**4)
				
				
			l32s_hag = l32s_hag - (0.56+1.93*Z2_hag)/Z2_hag/(1.+0.44*Z2_hag)*(f33_hag - (f33_hag**4)) + 4.95/(1.+2.48*Z2_hag)*((f33_hag**2)-(f33_hag**4) -0.55*(f33_hag**3)+0.55*(f33_hag**4)) - 1.2/(1.+0.5*Z2_hag)* (f33_hag**4)
			l34s_hag = (1.+1.4/(Z1_hag+1.))*f34_hag - 1.9/(Z1_hag+1.)* (f34_hag**2) + 0.3/(Z1_hag+1.)*(f34_hag**3) + 0.2/(Z1_hag+1.)*(f34_hag**4)
  
			LTI_hag = abs(self.dti[i] / self.ti[i])
			LN_hag = abs(self.dne[i] / self.ne[i])
			Lq_hag = abs(self.dq[i] /self.q[i])
			if self.use_hager: self.alps[i] = alps_hag*gama_hag

			f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.psin[i],self.betacol[i],self.betagradti[i],-alps_hag,gama_hag,(1.-2.*DPSI_hag*(-1.5*LTI_hag + LN_hag + Lq_hag)),self.nues[i]))
			
			self.jb_hag1[i] = (gam1_hag*l31s_hag)/self.pe[i]*self.dpt[i] + (gam2_hag*l32s_hag) / self.te[i] * self.dte[i] + (gam4_hag*l34s_hag) * (gama_hag*alps_hag) / Z2_hag / self.te[i] * (1.-2.*DPSI_hag*(-1.5*LTI_hag + LN_hag + Lq_hag))*self.dti[i]

			if (self.use_chang):
			    self.jb_hag2[i] = l31s_hag/self.pe[i]*self.dpt[i] + l32s_hag / self.te[i] * self.dte[i] + l34s_hag * alps_hag / Z1_hag / self.te[i] * self.dti[i]
				
			elif(self.use_hager):
				self.jb_hag2[i] = l31s_hag/self.pe[i]*self.dpt[i] + l32s_hag / self.te[i] * self.dte[i] + l34s_hag * alps_hag / Z2_hag / self.te[i] * self.dti[i]
				
			else:
				self.jb_hag2[i] = l31s_hag/self.ne[i]*self.dne[i] + rpe_hag*(l31s_hag+l32s_hag)/self.te[i]*self.dte[i] + (1.-rpe_hag)*(1.+l34s_hag/l31s_hag * alps_hag) * l31s_hag * self.dti[i]/self.ti[i]
			
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
		for i in range(self.num):
			f.write('%f %f %f \n'%(self.psin[i],betaw[i],betaw0[i]))
			
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

		return

	def make_vtheta(self):

		self.vpol = np.copy(self.alps); self.vpoli = np.copy(self.alps); self.vpolp = np.copy(self.alps); self.vpole = np.copy(self.alps); self.vpolE = np.copy(self.alps);
		self.Er = np.copy(self.alps); self.wexb = np.copy(self.alps); self.we = np.copy(self.alps); self.wtor = np.copy(self.alps); self.wi = np.copy(self.alps);
		self.wneo=np.copy(self.alps);
		self.btheta  = np.copy(self.vpol); self.alps[0] = 0.
		for i in range(1,len(self.vpol)-1):
			self.btheta[i] = (self.psin[i+1]-self.psin[i-1])/ (self.RR[i+1]-self.RR[i-1])/ self.RR[i] * self.psia
			if self.psin[i] > 0.993: self.btheta[i] = self.btheta[i-1]

		self.btheta[0] = 2.*self.btheta[1] - self.btheta[2]
		self.btheta[-1] = 2.*self.btheta[-2] - self.btheta[-3]

		f = open('vtheta.dat','w')
		f.write('%9s\t%13s\t%13s\t%13s\t%13s\t%13s\n'%('psin[a.u]','Vi*[m/s]','Ve*[m/s]','V||[m/s]','VE[m/s]','Vneo[m/s]'))
		for i in range(len(self.vpol)): 
			self.vpol[i] =  self.alps[i] * self.eq.fpol[0] * self.btheta[i] / (self.b02av[i]) * self.dti[i] * 1.e3
			self.vpolp[i] = -self.vt[i]  * 1.e3 * self.btheta[i] / self.eq.fpol[0] * self.RR[i]
			self.vpoli[i] = -self.dpi[i] * self.RR[i]*self.btheta[i]/e0/self.ni[i]/1.e19 / self.eq.fpol[0] * self.RR[i]
			self.vpole[i] = +self.dpe[i] * self.RR[i]*self.btheta[i]/e0/self.ne[i]/1.e19 / self.eq.fpol[0] * self.RR[i]
#			self.vpolp[i] = +elf.vpolp[i]* 0.75 #neoclassical toroidal offset
#			self.vpoli[i] = +self.vpoli[i]* 0.6  #Density width correction.
			self.vpolE[i] = +self.vpol[i] - self.vpolp[i] - self.vpoli[i]
			self.Er[i]    = -self.vpolE[i]* self.eq.fpol[0] / self.RR[i]
			self.wexb[i]  = +self.Er[i] / self.RR[i] / self.btheta[i]
			self.we[i]    = -self.dpe[i] * self.RR[i]*self.btheta[i]/e0/self.ne[i]/1.e19 / self.RR[i]/self.btheta[i]
			self.wi[i]    = +self.dpi[i] * self.RR[i]*self.btheta[i]/e0/self.ni[i]/1.e19 / self.RR[i]/self.btheta[i]
			self.wtor[i]  = +self.vt[i] * 1.e3 /self.RR[i]
			
			f.write('%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%f\n'%(self.psin[i],self.vpoli[i],self.vpole[i],self.vpolp[i],self.vpolE[i],self.vpol[i],self.alps[i]))

		plt.plot(self.psin,self.wi,self.psin,self.wexb,self.psin,self.we,self.psin,self.wtor)
		plt.show()
		Erf = interp1d(self.psin,self.wexb); eps = 1.e-4
		raf = interp1d(self.psinh,self.rh);
		qaf = interp1d(self.psinh,self.qh);
		dEr = np.copy(self.wexb*0.)
		for i in range(1,len(self.vpol)-2):
			psint = self.psin[i]
			dr    = raf(psint+eps) - raf(psint-eps)
			rr    = raf(psint)
			dErf  = (Erf(psint+eps)-Erf(psint-eps))/dr
			dErf  = dErf * rr / qaf(psint)
			e1 = 1.602*1.e-19; mi2 = 1.672*2.*1.e-27
			cs    = np.sqrt( 900 * e1 / mi2)
			dEr[i] = dErf/cs*0.35
#			print(self.psin[i],rr,dEr[i],dErf,Erf(psint))
		ind = np.where(self.psin>0.7);
#		print(max(dEr[ind]),dEr[90],self.psin[90])
		f.close()
		return

	def psin_grid(self,xr1,sig1,xr2,sig2,num):

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
		
		self.use_hager = False
		self.use_chang = True
		self.use_neo   = False
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
		self.cjf = 2.0
		self.djf = 1.5
		
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

	def __init__(self):
		
		self.initialise_var()
		self.psin_grid(self.grid_xr1,self.grid_sig1,self.grid_xr1,self.grid_sig1,self.grid_n)

		self.currdir = os.getcwd()
		self.f_hager = 'hager_map'
		self.chease_rundir = 'CTEMP'
		self.f_rhopsir = 'RHO_PSI_R'

		self.chease_rundir = self.currdir + '/' + self.chease_rundir	
		self.f_hager = self.chease_rundir + '/' + self.f_hager
		self.f_rhopsir = self.chease_rundir + '/' + self.f_rhopsir
		
		return

if __name__ == "__main__":

	import ercal
	sim = ercal.Er_profile()

	sim.f_eqdsk = sys.argv[1]
	sim.f_profile = sys.argv[2]
	try: bsm = int(sys.argv[3])
	except: 
		bsm = 1;
	print('>>> 1:satuer 2: hager')
	print('>>> BS model = %i'%bsm)
	if bsm == 2:
		sim.use_hager = True
		
	sim.read_kinprof_chease()
	sim.load_eqdsk()
	sim.read_hagermap()
	sim.interpol()
	sim.hager_bs()
	sim.make_vtheta()
	
