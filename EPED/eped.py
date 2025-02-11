#!/usr/local/anaconda3/bin/python3
import os,sys
import numpy as np
import time
import matplotlib.pyplot as plt
from shutil import copyfile, move
from scipy.interpolate import interp1d
import subprocess
import ch_tool as ch
from progress import update_progress
from nodelist import make_nodelist
import random
from batch_run import *
from exec_dirs import chease_dir,chease_exec,elite_exec,elite_dir,python3_exec,stab_dir,mis_exec,jastab_dir,japlot_dir,version,author,qdel_exec,node_default,node_machine

mu0 = 4. * np.pi * 1.e-7;
e0 = 1.602 * 1.e-19;
class eped:

	def make_directories(self):

		try:
			os.mkdir('equil')
		except:
			pass
			
		try:
			os.mkdir('stab')
		except:
			pass

		try:
			os.mkdir('script')
		except:
			pass

		os.chdir('equil')
		for i in range(self.w_num):
			try:
				os.mkdir('scan_%i'%i)
			except:
				pass
		os.chdir('../')

		os.chdir('stab')
		for i in range(self.w_num):
			try:
				os.mkdir('scan_%i'%i)
			except:
				pass
		os.chdir('../')		
			
		return		

	def get_perimeter(self):

		if not(self.ch.eqdsk_name.lower() == 'none'):
			if (self.ch.target_psin == 1.0):
				self.ch.rzbdyp = np.copy(self.ch.eq.rzbdy)
			else:
				self.ch.rzbdyp = np.copy(self.ch.eq.rzbdyt)

		elif(self.ch.use_param_shape):
			self.ch.make_param_bnd()
		else:
			self.ch.read_bnd_shape(self.ch.bnd_file)

		bndrz = np.copy(self.ch.rzbdyp)

		sum = 0.

		for i in range(len(bndrz)):

			j = i + 1
			if (j==(len(bndrz))):
				j = 0

			dx = abs(bndrz[i,0]-bndrz[j,0])
			dz = abs(bndrz[i,1]-bndrz[j,1])
			dl = np.sqrt(dx**2. + dz**2.)

			sum = sum + dl

		self.perim = sum

		return 

	def make_kbm(self,width):

		bped = (width / self.kbm_coef1) ** (1./self.kbm_coef2)
		bpav = mu0 * self.ch.ip / self.perim

		pped = bped * (bpav**2. / 2. / mu0)

		tped = pped / (self.neped * 1.e19) / e0 / 1.e3 / (2.0 - (self.ch.zeff-1.)/self.ch.zimp) #[ev->keV]

		#print(width,bped,pped,self.neped,tped,self.ch.zeff,self.ch.zimp,self.perim,self.ch.ip)

		return (tped)

	def make_eped_prof(self,width,filename='chease_eped_ref'):

		#fi = open(filename,'r')
		fo = open('chease_eped','w')

		eped_prof = np.zeros(shape=(3,8))

	#	for i in range(3):
#			line = fi.readline()
#			if not line: break
#			line = line.split()
#			eped_prof[i,:] = np.array(line,dtype='float')
#
#		fi.close()

		teped = self.make_kbm(width)

		eped_prof[0,0] = (self.neped - self.neped*self.nesep) / 2. / np.tanh(1.)
		eped_prof[0,1] = self.neped * self.nesep
		eped_prof[0,2] = self.an1
		eped_prof[0,3] = 1. - width*0.5
		eped_prof[0,4] = width
		eped_prof[0,5] = 1. - width
		eped_prof[0,6] = self.alpn1
		eped_prof[0,7] = self.alpn2
		eped_prof[1,0] = (teped - self.tesep) / 2. / np.tanh(1.)
		eped_prof[1,1] = self.tesep
		eped_prof[1,2] = self.at1
		eped_prof[1,3] = 1. - width*0.5
		eped_prof[1,4] = width
		eped_prof[1,5] = 1. - width
		eped_prof[1,6] = self.alpt1
		eped_prof[1,7] = self.alpt2

		eped_prof[2,:] = np.copy(eped_prof[1,:])

		for i in range(3):
			line = ''
			for j in range(8):
				line = line + '%8.6f\t'%eped_prof[i,j]
			line = line + '\n'

			fo.write(line)

		fo.close()

		self.write_eped_prof(eped_prof)

		return
		
	def eped_fun(self,x,m):
	
		y= m[1] + m[0]*(np.tanh(2.*(m[3]-x)/m[4])-np.tanh(2.*(m[3]-1.0)/m[4]))
		
		y = y + m[2] * (abs(1. - (x/m[5])**m[6])**m[7] )* 0.5 * (1.0 + np.sign(m[5]-x))
	
		return y
		
	def write_eped_prof(self,eped_prof,filename='chease_kinprof'):
	
		f = open(filename,'w')
		num = 401
		psin = np.linspace(0,1.,num)
		f.write('%i\n'%num)
		f.write('%f\t%f\t%f\t%f\n'%(self.ch.zeff,self.ch.zimp,self.ch.amain,self.ch.aimp))
		for i in range(num):
		
			ne = self.eped_fun(psin[i],eped_prof[0,:])
			te = self.eped_fun(psin[i],eped_prof[1,:])
			ti = self.eped_fun(psin[i],eped_prof[2,:])
			
			ni = ne * (1.0 - (self.ch.zeff -1.)/self.ch.zimp)
		
			f.write('%f\t%f\t%f\t%f\t%f\n'%(psin[i],te,ne,ti,ni))
			
		f.close()
	
		return

	def make_elite_inp(self,moden,filename2,filename='elite.in'):

		fi = open(filename,'r')
		fo = open(filename2,'w')

		while True:

			line = fi.readline()
			if not line: break

			if (line.find('nn') > -1):
				line = 'nn = %i\n'%moden
			fo.write(line)

		fi.close()
		fo.close()

		return

	def make_elite_run(self):

			f = open('elite_run.sh','w')
			f.write('#!/bin/bash\n')
			f.write('rm result.dat\n')
			line = ''
			for j in range(len(self.mode_n)):
				nn = int(self.mode_n[j])
				f.write(elite_exec + 'eq5.7 work_%i \n'%nn)
				f.write(elite_exec + 'vac5.7 work_%i \n'%nn)
				f.write(elite_exec + '5.7 work_%i \n'%nn)
				f.write('cat work_$i.gamma >> result.dat \n')

				f.write('\n')

			f.write('rm *.fun2d *.plas *.surf \n')

			f.close()

			return

	def make_equ_script(self):

		for i in range(self.w_num):

			filename = self.currdir + '/script/scan_equ_%i'%i
			nodelist = self.node_list
			log_e = filename +'.e'
			log_o = filename +'.o'

			command = 'cd '+self.currdir+'/equil/scan_%i \n'%i
			command = command + chease_dir + ' ' + self.inputfile + ' nomap \n'
			command = command +'mv CHEASE/EXPEQ ' + self.currdir+'/stab/scan_%i/EXPEQ \n'%i 
			command = command +'mv CHEASE/chease_namelist ' + self.currdir+'/stab/scan_%i/chease_namelist \n'%i 
			command = command +'mv chease_kinprof ' + self.currdir+'/stab/scan_%i/chease_kinprof \n'%i 
			command = command +'cp ' + self.inputfile + ' ' + self.currdir+'/stab/scan_%i/'%i + self.inputfile + ' \n'
			command = command +'cp chease_eped ' + self.currdir+'/stab/scan_%i/chease_eped \n'%i 
			if not (self.ch.use_param_shape):
				command = command +'cp ' + self.ch.bnd_file + ' ' + self.currdir+'/stab/scan_%i/'%i + self.ch.bnd_file + ' \n' 
			if not (self.ch.eqdsk_name.lower() == 'none'):
				command = command +'cp ' + self.ch.eqdsk_name + ' ' + self.currdir+'/stab/scan_%i/.'%i + ' \n' 

			command = command +"echo '"+str(i)+"' >> " + self.currdir + '/log.batch_equ \n'
			make_batch_script(filename,nodelist,log_e,log_o,command,'EPED_EQU')	

		return

	def make_stab_script(self):

		for i in range(self.w_num):

			filename = self.currdir + '/script/scan_stab_%i'%i
			nodelist = self.node_list
			log_e = filename +'.e'
			log_o = filename +'.o'

			command = 'cd '+self.currdir+'/stab/scan_%i \n'%i
			command = command +python3_exec + ' ' + stab_dir + ' ' + self.inputfile + ' %i\n'%i
			command = command +"echo '"+str(i)+"' >> " + self.currdir + '/log.batch_stab \n'
			make_batch_script(filename,nodelist,log_e,log_o,command,'EPED_STAB')	

		return		

	def submit_script(self,isequ=True):

		batch_num = []
		batch_index = np.zeros(shape=(self.w_num,2))

		if isequ:
			batch_log_dir = self.currdir + '/log.batch_equ'
		else:
			batch_log_dir = self.currdir + '/log.batch_stab'

		if os.path.isfile(batch_log_dir):
			os.remove(batch_log_dir)

		f = open(batch_log_dir,'w')
		f.close()
	
		script_dir = self.currdir + '/script'
	
		for i in range(self.w_num):
				
			if (isequ):
				filename = script_dir + '/scan_equ_%i'%i
			else:
				filename = script_dir + '/scan_stab_%i'%i
			
			run_id = submit_batch_script(filename)	
			batch_num.append(str(run_id))
			batch_index[i,0] = run_id
			batch_index[i,1] = 0
		
		print('>>> Scripts are submitted')
		if (isequ):
			print('>>> Collect EQU data')
		else:
			print('>>> Collect STAB data')
		
		batch_chk = 1

		bat_count0 = 1
		while ( bat_count0 > 0):
			output = batch_script_status()
			bat_count = 0
			for i in range(self.w_num):
					index1 = i
					batch_chk = output.find(batch_num[index1])
					if (batch_chk > -1):
						bat_count = bat_count + 1;
						batch_index[i,1] = 1

			if not (bat_count == bat_count0):
				bat_count0 = bat_count
				update_progress(1.0 - bat_count0/self.w_num)

			self.check_submit(batch_index,batch_log_dir)
			if (bat_count0 > 0):
				time.sleep(30)
		if (isequ):
			print('>>> EQU scripts are finished')
		else:
			print('>>> STAB scripts are finished')
		print('>>> Finalizing process --> Collect data')
		
		return

	def check_submit(self,index,filename):

		### check the wrong ended & zombie batch script in the schduler!
		f = open(filename,'r')
		while True:
			line = f.readline()
			if not line: break
			ind = int(line)
			if (index[ind,1] == 0):
				status, output = subprocess.getstatusoutput(qdel_exec + ' '+ str(index[ind,0]))
				print('>>> Delete wrong terminated script from Job scheduler -> #%i'%index[ind,0])
				
		f.close()
			
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

	def run_chease(self):
	
		try:
			stat,out = subprocess.getstatusoutput(chease_exec+' > log.chease')
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
			
		return	

	def find_max_ja(self,psint,psin,tarray,max_min):
	
		for i in range(len(psin)-1):
			if (psin[i] > psint): break

		if (max_min == 0):
			maxv = min(tarray[i:len(tarray)+1])
		else:
			maxv = max(tarray[i:len(tarray)+1])
		
		return maxv		
		
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
		rho=ne[-1]*md   ##############################ne[-1] - > ne[0]
#		rho=nmainion[0]*md
		omegan=2.9979e10/ne/4.8032e-10*dne*1.6022e-12*ti
		omegaspi=(omegan+2.9979e10/ne/4.8032e-10*dti*1.6022e-12*ne)/1e8
		omega_alf=np.sqrt(B*B/rho/rav/rav/mu0)#/q[0]/q[0])
		alf_over_maxomega=-omega_alf/np.min(omegaspi)#[0:200])

		file.close()
		
		return(npsi,psi,dpdpsi,ffp,te,ne,ti,zeff,zimp,alf_over_maxomega)		
		
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

		f4.write('nowindow=.f.\n')
		f4.write('bug(1)=0.0\n')
		f4.write('&end\n')

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
			
		f4 = open(filename_edge,'r')
		for i in range(7):
			line = f4.readline()
		line = f4.readline().split()
		alpha = float(line[0])
		jphi = float(line[3])	

		return monotonicq

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
		f4.write('EPSLON = %e\n'%self.epslon)
		f4.write('/\n');
		
		f4.close()
					
		return			
		
	def collect_data(self,filename='output'):

		n_num = len(self.mode_n)
		dat = np.zeros(shape=(self.w_num,n_num,14))
		gra = np.zeros(shape=(self.w_num,2))
		grt = np.zeros(shape=(self.w_num,n_num))
		grd = np.copy(gra)
		temp = np.zeros(shape=(n_num,14))
		self.result = np.zeros(shape=(self.w_num,14,2))
		nofile = False
		for i in range(self.w_num):
		
			outname = self.currdir + '/stab/scan_%i/growth.dat'%i
			try:	
				f = open(outname,'r')
				nofile = False
			except:	nofile = True
			if nofile:
				print('>>> Problems in #%i'%i)
				for j in range(n_num):
					dat[i,j,0] = i
					dat[i,j,:] = 0.
					dat[i,j,7] = 3.
					dat[i,j,8] = 0.
					dat[i,j,9] = 0.
					dat[i,j,11] = 1.
			else:
				for j in range(n_num):
					line = f.readline().split()
					try:
						dat[i,j,:] = np.array(line,dtype=float)
					except:
						dat[i,j,11] = 1.
						pass
					if (np.isnan(dat[i,j,8]) or dat[i,j,7] < self.ncrit):
						dat[i,j,8] = 0.
						dat[i,j,9] = 0.
						dat[i,j,11] = 1.
					
				f.close()	
			
		for i in range(self.w_num):
				
			ind = np.argmax(dat[i,:,8])
		
			gra[i,0] = dat[i,ind,8]
			gra[i,1] = dat[i,ind,7]
			self.result[i,:,0] = np.copy(dat[i,ind,:])
			for j in range(n_num):
				ncrit = self.nqcrit / dat[i,j,11]
				if not (ncrit == 0.0):
					if (ncrit < dat[i,j,7]):
						grt[i,j] = dat[i,j,9] * dat[i,j,7] / ncrit
					else:
						grt[i,j] = dat[i,j,9]
		
		for i in range(self.w_num):
		
			ind = np.argmax(grt[i,:])
			
			grd[i,0] = grt[i,ind]
			grd[i,1] = dat[i,ind,7]
			self.result[i,:,1] = np.copy(dat[i,ind,:])
			self.result[i,9,1] = grt[i,ind]

		w1 = np.array([])
		gra1 = np.array([])
		grd1 = np.array([])

		index = np.linspace(0,len(self.w)-1,len(self.w))
		for i in range(len(self.exclude)):
			index = [x for x in index if x !=self.exclude[i]]

		for i in index:
			i = int(i)
			w1 = np.append(w1,self.w[i])
			gra1 = np.append(gra1,gra[i,0])
			grd1 = np.append(grd1,grd[i,0])
		fit1 = interp1d(w1,gra1,'slinear');
		fit2 = interp1d(w1,grd1,'slinear');
		
		width1 = self.find_target_point(fit1,self.grcrit1)
		width2 = self.find_target_point(fit2,self.grcrit2)

		self.width1 = width1
		self.width2 = width2

		self.tped1 = self.make_kbm(width1)
		self.tped2 = self.make_kbm(width2)		

		self.gra = np.copy(gra)
		self.grd = np.copy(grd)

		self.fit1 = fit1
		self.fit2 = fit2

		if (self.make_equ):
			print('>>> Make result equ... ')
			self.make_result_equ(width1,width2)
		self.print_result(width1, width2)
		self.write_result(width1, width2)

		if (self.plot_only):
			self.plot_result(gra,grd,width1,width2,fit1,fit2)


		print('>>> Run finished')	
						
		return
		
	def find_target_point(self,fit,crit):
	
		w = np.linspace(self.w_min,self.w_max,1.e3)
		gr = fit(w)
		
		width = -1
		
		for i in range(999):
			if (gr[i] > crit):
			
				width = w[i]
				break
		if (width > 0):			
			return width
		else:
			return 0.
		
	def make_result_equ(self,width1,width2):
	
		os.chdir(self.currdir)
	
		try:
			os.mkdir('result')
		except:
			pass
			
		os.chdir('result')	
		try:
			os.mkdir('temp')
		except:
			pass
		copyfile(self.currdir+'/'+self.inputfile,self.inputfile)
		if not (self.ch.use_param_shape):
			copyfile(self.currdir + '/' + self.ch.bnd_file,self.ch.bnd_file)
		if not (self.ch.eqdsk_name.lower() == 'none'):
			copyfile(self.currdir + '/' + self.ch.eqdsk_name,self.ch.eqdsk_name)

		if (width1 > 0.):
			self.make_eped_prof(width1,self.currdir+'/chease_eped')
			os.system(chease_dir + ' ' + self.inputfile + ' nomap')
			move('chease_eped','temp/chease_eped1')
			move('chease_kinprof_new','temp/chease_kinprof1')
			move('OUTPUT/EQDSK_OUT','geqdsk1')
			move('CHEASE/NJA','temp/NJA1')
		
		if (width2 > 0.):
			self.make_eped_prof(width2,self.currdir+'/chease_eped')
			os.system(chease_dir + ' ' + self.inputfile + ' nomap')
			move('chease_eped','temp/chease_eped2')
			move('chease_kinprof_new','temp/chease_kinprof2')
			move('OUTPUT/EQDSK_OUT','geqdsk2')
			move('CHEASE/NJA','temp/NJA2')

		os.system('rm -r OUTPUT')
		os.system('rm -r CHEASE')

		os.chdir(self.currdir)

		return	

	def plot_result(self,gra,grd,width1,width2,fit1,fit2,figv=None):

		if figv == None:
			fig, [ax1, ax2, ax3] = plt.subplots(1,3,figsize=(18,5))
		else:
			fig = figv
			[ax1,ax2, ax3] = fig.axes
			ax1.cla()
			ax2.cla()
			ax3.cla()

		ax1.scatter(self.w,gra[:,0])
		for i in range(self.w_num):
			if (gra[i,0] > 0.):
				ax1.text(self.w[i],gra[i,0]+0.002,str(int(gra[i,1])))
			ax1.text(self.w[i],0.000,str(int(i)),color='lime')
		ax1.plot(self.w,fit1(self.w),color='g',linestyle='--');
		for i in range(len(self.exclude)):
			i = int(self.exclude[i])
			ax1.scatter(self.w[i],gra[i,0],color='r')

		ax1.axhline(y=self.grcrit1,color='magenta',linestyle='--')
		ax1.axvline(x=width1,color='magenta',linestyle='--')
		ax1.text(self.w_min,self.grcrit1+0.002,'$\gamma/\omega_A$ = %3.2f'%self.grcrit1)
		if width1 > 0.:
			ax1.scatter(width1,fit1(width1),color='magenta',s=100,marker='x')
		if width2 > 0:
			ax1.scatter(width2,fit1(width2),color='orange', s=100,marker='x')
		ax1.set_ylabel('growth rate ($\gamma/\omega_A$) [a.u]')
		ax1.set_xlabel('Pedestal Width ($\Delta_{ped}$) [psin]')
		ax1.set_title('EPED results with criterion 1')
		ax1.axis([self.w_min*0.9, self.w_max*1.1, 0., 0.1])
		if width1 > 0.:
			ax1.text(self.w_min,0.01,'$\Delta_{ped}$=%f'%width1)

		ax2.scatter(self.w,grd[:,0])
		for i in range(self.w_num):
			if (grd[i,0]>0.):
				ax2.text(self.w[i],grd[i,0]+0.05,str(int(grd[i,1])))
			ax2.text(self.w[i],0.0,str(int(i)),color='lime')
		ax2.plot(self.w,fit2(self.w),color='g',linestyle='--');
		for i in range(len(self.exclude)):
			i = int(self.exclude[i])
			ax2.scatter(self.w[i],grd[i,0],color='r')

		ax2.axhline(y=self.grcrit2,color='orange',linestyle='--')
		ax2.text(self.w_min,self.grcrit2+0.02,'$\gamma/\omega_i*$ = %3.2f'%self.grcrit2)
		ax2.axvline(x=width2,color='orange',linestyle='--')
		if width1 > 0.:
			ax2.scatter(width1,fit2(width1),color='magenta',s=100,marker='x')
		if width2 > 0.:
			ax2.scatter(width2,fit2(width2),color='orange', s=100,marker='x')

		ax2.set_ylabel('growth rate ($\gamma/\omega_i*$) [a.u]')
		ax2.set_xlabel('Pedestal Width ($\Delta_{ped}$) [psin]')
		ax2.set_title('EPED results with criterion 2')		
		ax2.axis([self.w_min*0.9, self.w_max*1.1, 0., 0.8])
		if width2 > 0.:
			ax2.text(self.w_min,0.1,'$\Delta_{ped}$=%f'%width2)

		tped = self.make_kbm(self.w)
		tped1 = self.make_kbm(width1)
		tped2 = self.make_kbm(width2)

		ax3.scatter(self.w,tped)
		fit3 = interp1d(self.w,tped,'slinear')
		ax3.plot(self.w,fit3(self.w),color='g',linestyle='--')

		legend = []
		legend2 = []

		if (self.teped > 0.):
			sc1 = ax3.scatter(self.tewidth,self.teped,s=100,marker='*')
			legend.append(sc1)
			legend2.append('Exp Ref')


		if width1 > 0.:
			leg1 = ax3.scatter(width1,tped1,color='magenta',s=100,marker='x')
			ax3.axvline(x=width1,color='magenta',linestyle='--')
			ax3.axhline(y=tped1,color='magenta',linestyle='--')			
			legend.append(leg1)
			legend2.append('Criterion1 $\gamma/\omega_A$ = %3.2f'%self.grcrit1)
		if width2 > 0.:
			leg2 = ax3.scatter(width2,tped2,color='orange', s=100,marker='x')
			ax3.axvline(x=width2,color='orange',linestyle='--')
			ax3.axhline(y=tped2,color='orange',linestyle='--')
			legend.append(leg2)
			legend2.append('Criterion2 $\gamma/\omega_{i*}$ = %3.2f'%self.grcrit2)

		ax3.set_ylabel('Pedestal height ($T_{ped}$) [keV]')
		ax3.set_xlabel('Pedestal Width ($\Delta_{ped}$) [psin]')
		ax3.set_title('KBM constraint')	
		ax3.axis([self.w_min*0.9, self.w_max*1.1, tped[0]*0.9, tped[-1]*1.1])
		if (width2 > 0 and width1 > 0):
			line = '$\delta\Delta_{ped}$ [%%]  = %4.2f \n'%((width2/width1-1)*100) + '$\delta  T_{ped}$ [%%]  = %4.2f \n'%((tped2/tped1-1)*100)
			ax3.text(self.w_min,tped[-1]*0.9,line)
		ax3.legend(legend,legend2,loc=4)

		fig.tight_layout()
		#plt.savefig('EPED_result')
		
		if (figv == None):
			plt.show(block=False)
			input('Press Enter...')
		else:
			plt.draw()

		return

	def print_result(self,width1,width2):

		if width1 > 0.:
			print('>>> EPED prediction for criterion1 gamma/w_A = %3.2f'%self.grcrit1)
			tped = self.make_kbm(width1)
			print('>>> Width = %f, T_{ped} = %f [keV], n_{ped} = %f [10(19)/m3] '%(width1,tped,self.neped))
		else:
			print('>>> No available results for criterion1\n')

		if width2 > 0.:
			print('>>> EPED prediction for criterion2 gamma/w_i* = %3.2f'%self.grcrit2)
			tped = self.make_kbm(width2)
			print('>>> Width = %f, T_{ped}$ = %f [keV], n_{ped} = %f [10(19)/m3] '%(width2,tped,self.neped))
		else:
			print('>>> No available results for criterion2\n')
		
		return

	def write_result(self,width1,width2):

		f = open('result.dat','w')

		f.write('id\twidth\t\ttped\t\tALPHA\t\tJPHI[MA/m2]\tq95\t\tnn\tgrowth\n')

		tped = self.make_kbm(self.w)
		for i in range(self.w_num):
			dat = np.copy(self.result[i,:,0])
			f.write('%i\t%f\t%f\t%f\t%f\t%f\t%i\t%f\n'%(dat[0],self.w[i],tped[i],dat[4],dat[6],dat[11],dat[7],dat[8]))
		f.write('\n')
		for i in range(self.w_num):
			dat = np.copy(self.result[i,:,1])
			f.write('%i\t%f\t%f\t%f\t%f\t%f\t%i\t%f\n'%(dat[0],self.w[i],tped[i],dat[4],dat[6],dat[11],dat[7],dat[9]))
		f.write('\n')
		if width1 > 0.:
			f.write('EPED prediction for criterion1 gamma/w_A = %3.2f\n'%self.grcrit1)
			tped = self.make_kbm(width1)
			f.write('Width = %f, Tped = %f [keV], nped = %f [10(19)/m3] \n'%(width1,tped,self.neped))
		else:
			f.write('No available results for criterion1\n')

		if width2 > 0.:
			f.write('EPED prediction for criterion2 gamma/w_i* = %3.2f\n'%self.grcrit2)
			tped = self.make_kbm(width2)
			f.write('Width = %f, Tped = %f [keV], nped = %f [10(19)/m3] \n'%(width2,tped,self.neped))
		else:
			f.write('No available results for criterion2\n')

		f.close()

		return

	def make_eqdsk_opt(self,filename='eqdsk_opt',width=0.05):

		f = open(filename,'w')
		f.write('WIDTH = %f\n'%width)
		f.close()

		return

	def run_eped(self):

		self.make_directories()

		for i in range(self.w_num):
			os.chdir('equil')
			os.chdir('scan_%i'%i)
			copyfile(self.currdir+'/' + self.inputfile,self.inputfile)
			if not (self.ch.use_param_shape):
				copyfile(self.currdir + '/' + self.ch.bnd_file,self.ch.bnd_file)
			if not (self.ch.eqdsk_name.lower() == 'none'):
				copyfile(self.currdir + '/' + self.ch.eqdsk_name,self.ch.eqdsk_name)

			self.make_eped_prof(self.w[i],self.currdir+'/chease_eped')
			self.make_eqdsk_opt('eqdsk_opt',self.w[i])
			os.chdir(self.currdir)
			os.chdir('stab')
			os.chdir('scan_%i'%i)
			os.chdir(self.currdir)

		self.make_equ_script()
		self.make_stab_script()
		self.submit_script(True)
		self.submit_script(False)

		self.clear_files()

		return

	def make_ja_namelist(self,filename = 'ja_namelist_eped',fileid = 'runid'):

		f = open(filename,'w')

		n_min = min(self.mode_n)
		n_max = max(self.mode_n)

		if (n_min < 5):
			mode_n = '%i, 5'%n_min
		else:
			mode_n = '5'

		nn = int(round((n_max-5.)/5.))

		if (nn> 5): nn = 5

		for i in range(nn):
			mode_n = mode_n + ', %i'%(5*(i+2))

		runid = random.randint(1,500)	

		f2 = open(fileid,'w')
		f2.write('%i'%runid)
		f2.close()

		f.write('***** JA diagram namelist ***** (by S.Kim ver.1)\n')
		f.write('--- Run option ---\n')
		f.write('run_id = %i 					# set run id\n'%runid)
		f.write('equ_name = geqdsk                              # input g-file\n')
		f.write('batch_run = True				# use batch run\n')
		f.write('\n')
		f.write('--- Kinetic profile option\n')
		f.write('adjust_prof = False			# Adjust profile in equilibrium construction\n')
		f.write('\n')
		f.write('--- Scan options ---\n')
		f.write('grid_n = 8			# # of radial point\n')
		f.write('pmin = -0.5		# minimum p\n')
		f.write('pmax = +0.5		# maximum p\n')
		f.write('jmin = -0.5		# minimum j\n')
		f.write('jmax = +0.5		# maximum j\n')
		f.write('\n')
		f.write('--- CHEASE options ---\n')
		f.write('target_bnd = 0.999\n')
		f.write('CNS = 80\n')
		f.write('CNT = 80\n')
		f.write('CNPSI = 250\n')
		f.write('CNCHI = 150\n')
		f.write('\n')
		f.write('CNSE = 150\n')
		f.write('CNTE = 200\n')
		f.write('CNPSIE = 300\n')
		f.write('CNCHIE = 512\n')
		f.write('\n')
		f.write('--- Stability options ---\n')
		f.write('run_stab = %s                       # or elite\n'%self.run_stab)
		f.write('qdelfix = 0.3                           # fixed qa_del\n')
		f.write('use_compression = %s\n'%self.use_comp)
		f.write('\n')
		f.write('--- MISHKA options ---\n')
		f.write('moden = %s			# mode n to calculate\n'%mode_n)
		f.write('gridn = 301 							# raidal grid for stab cal\n')
		f.write('psis = 0.75							# radial start pt\n')
		f.write('xr1 = 1.0                             	# grid options xr,sig\n')
		f.write('sig1 = 0.05\n')
		f.write('xr2 = 0.9\n')
		f.write('sig2 = 0.1\n')
		f.write('\n')
		f.write('--- ELITE options ---\n')
		f.write('modene = %s\n'%mode_n)
		f.write('psise = 0.6\n')
		f.write('ndist = 50\n')
		f.write('\n')
		f.write('--- Machine options ---\n')
		f.write('nodelist = %s\n'%self.node_list2)

		f.close()

		return

	def prepare_ja(self):

		if (self.width1 > 0. and self.width2 > 0.):
			equ_type = self.get_input('Which equil type ? (1/2)',[1,2])
		elif (self.width1 > 0.):
			equ_type = 1
			print('>>> Equ_type %i is only available...'%equ_type)
		elif (self.width2 > 0.):
			equ_type = 2
			print('>>> Equ_type %i is only available...'%equ_type)
		else:
			print('>>> No pedestal width is found from EPED')
			exit()

		print('>>> Pedestal width type %i is selected...'%equ_type)

		try:
			os.mkdir('jastab_%i'%equ_type)
		except:
			pass
		try:
			os.mkdir('jastab_%i/input'%equ_type)
		except:
			pass

		copyfile('result/geqdsk%i'%equ_type,'jastab_%i/input/geqdsk'%equ_type)
		copyfile('result/chease_kinprof%i'%equ_type,'jastab_%i/input/chease_kinprof'%equ_type)
		self.make_ja_namelist('jastab_%i/ja_namelist_eped'%equ_type,'jastab_%i/runid'%equ_type)

		print('>>> JA STAB started...')

		os.chdir('jastab_%i'%equ_type)
		os.system(python3_exec + ' ' + jastab_dir + ' ja_namelist_eped')

		return

	def draw_ja(self):

		is1 = os.path.isdir('jastab_1')
		is2 = os.path.isdir('jastab_2')

		if (is1 and is2):
			equ_type = self.get_input('Which equil type ? (1/2)',[1,2])
		elif (is1):
			equ_type = 1
		elif (is2):
			equ_type = 2
		else:
			print('>>> No ja stability results...')
			exit()

		print('>>> Draw ja stability diagram for equil = %i'%equ_type)

		f = open('jastab_%i/runid'%equ_type)
		runid = int(float(f.readline()))
		f.close()

		if (self.run_stab == 'mishka'):
			result = 'jastab_%i/ja_%i/ja_result_%i_mis'%(equ_type,runid,runid)
		elif (self.run_stab == 'elite'):
			result = 'jastab_%i/ja_%i/ja_result_%i_eli'%(equ_type,runid,runid)

		os.system(python3_exec + ' ' + japlot_dir + ' ' + result)

		return

	def get_input(self,message,var_ind,istr=False):

		var_ind = np.array(var_ind,dtype='float')

		print('>>> '+ message)

		chk = True

		while chk:

			var = input('>>> ')
			if istr:
				var = var.lower()
			else:
				var = float(var)	

			for i in range(len(var_ind)):
				if var_ind[i] == var:
					chk = False
					break

		return var

	def clear_files(self):

		print('>>> Remove temporary files')

		for i in range(self.w_num):

			update_progress((i+1)/self.w_num)

			self.remove_file('equil/scan_%i/fort.12'%i)
			self.remove_file('equil/scan_%i/eliteinp'%i)
			self.remove_file('equil/scan_%i/CHEASE/EQDSK_COCOS_02.OUT'%i)
			self.remove_file('equil/scan_%i/CHEASE/EQDSK_COCOS_02_POS.OUT'%i)
			self.remove_file('equil/scan_%i/CHEASE/NDES'%i)
			self.remove_file('equil/scan_%i/CHEASE/NOUT'%i)

			for j in range(len(self.mode_n)):

				nn = int(float(self.mode_n[j]))
				self.remove_file('stab/scan_%i/work_%i.plas'%(i,nn))
				self.remove_file('stab/scan_%i/work_%i.fun2d'%(i,nn))
				self.remove_file('stab/scan_%i/work_%i.eq'%(i,nn))
				self.remove_file('stab/scan_%i/work_%i.eqdat'%(i,nn))
				self.remove_file('stab/scan_%i/work_%i.eqout'%(i,nn))
				self.remove_file('stab/scan_%i/work_%i.jmg'%(i,nn))
				self.remove_file('stab/scan_%i/work_%i.surf'%(i,nn))
				self.remove_file('stab/scan_%i/work_%i.surfdat'%(i,nn))
				self.remove_file('stab/scan_%i/work_%i.vac'%(i,nn))
				self.remove_file('stab/scan_%i/work_%i.xtopsi'%(i,nn))
				self.remove_file('stab/scan_%i/CASPLOT'%i)
				self.remove_file('stab/scan_%i/eig_str'%i)
				self.remove_file('stab/scan_%i/eqin'%i)
				self.remove_file('stab/scan_%i/equ/NELITE'%i)
				self.remove_file('stab/scan_%i/equ/NELITE2'%i)
				self.remove_file('stab/scan_%i/equ/NDES'%i)
				self.remove_file('stab/scan_%i/equ/NMISHKA'%i)
				self.remove_file('stab/scan_%i/equ/NOUT'%i)
				self.remove_file('stab/scan_%i/fort.12'%i)
				self.remove_file('stab/scan_%i/fort.22'%i)
				self.remove_file('stab/scan_%i/fort.76'%i)
		return

	def remove_file(self,filename):

		try:
			os.remove(filename)
		except:
			pass

		return

	def read_namelist(self,filename):
	
		f4 = open(filename,'r')
		
		while True:
		
			line = f4.readline()
			if not line: break
	
			self.kbm_coef1 = self.ch.read_namelist_str(line,'KBM_COEF',self.kbm_coef1,2)
			self.kbm_coef2 = self.ch.read_namelist_str(line,'KBM_POWER',self.kbm_coef2,2)
			self.node_list = self.ch.read_namelist_str(line,'NODE_LIST',self.node_list,3)
			self.moden 	= self.ch.read_namelist_str(line,'MODEN',self.moden,3)
			self.run_stab = self.ch.read_namelist_str(line,'RUN_STAB',self.run_stab,3)

			self.neped = self.ch.read_namelist_str(line,'NEPED',self.neped,2)
			self.nesep = self.ch.read_namelist_str(line,'NESEP',self.nesep,2)
			self.tesep = self.ch.read_namelist_str(line,'TESEP',self.tesep,2)
			self.w_min = self.ch.read_namelist_str(line,'W_min',self.w_min,2)
			self.w_max = self.ch.read_namelist_str(line,'W_max',self.w_max,2)
			self.w_num = self.ch.read_namelist_str(line,'W_STEP',self.w_num,1)

			self.at1 = self.ch.read_namelist_str(line,'AT1',self.at1,2)
			self.an1 = self.ch.read_namelist_str(line,'AN1',self.an1,2)
			self.alpt1 = self.ch.read_namelist_str(line,'ALPT1',self.alpt1,2)
			self.alpt2 = self.ch.read_namelist_str(line,'ALPT2',self.alpt2,2)
			self.alpn1 = self.ch.read_namelist_str(line,'ALPN1',self.alpn1,2)
			self.alpn2 = self.ch.read_namelist_str(line,'ALPN2',self.alpn2,2)
			
			self.nqcrit = self.ch.read_namelist_str(line,'NQA',self.nqcrit,2)

			self.gridn = self.ch.read_namelist_str(line,'MIS_GRIDN',self.gridn,1)
			self.psis = self.ch.read_namelist_str(line,'MIS_PSISTART',self.psis,2)
			self.xr1 = self.ch.read_namelist_str(line,'MIS_XR1',self.xr1,2)
			self.sig1 = self.ch.read_namelist_str(line,'MIS_SIG1',self.sig1,2)
			self.xr2 = self.ch.read_namelist_str(line,'MIS_XR2',self.xr2,2)
			self.sig2 = self.ch.read_namelist_str(line,'MIS_SIG2',self.sig2,2)			

			self.use_comp = self.ch.read_namelist_str(line,'ELI_COMP',self.use_comp,4)
			self.ngride = self.ch.read_namelist_str(line,'ELI_GRIDN',self.ngride,2)
			self.psise = self.ch.read_namelist_str(line,'ELI_PSISTART',self.psise,2)
			self.ndist = self.ch.read_namelist_str(line,'ELI_NDIST',self.ndist,2)

			self.ncrit = self.ch.read_namelist_str(line,'NCUTOFF',self.ncrit,1)
			self.teped = self.ch.read_namelist_str(line,'TEPED',self.teped,2)
			self.tewidth = self.ch.read_namelist_str(line,'TEWIDTH',self.tewidth,2)
			self.grcrit1 = self.ch.read_namelist_str(line,'GRCRIT1',self.grcrit1,2)
			self.grcrit2 = self.ch.read_namelist_str(line,'GRCRIT2',self.grcrit2,2)

			self.epslon = self.ch.read_namelist_str(line,'EPSILON',self.epslon,2)
			self.exclude = self.ch.read_namelist_str(line,'EXCLUDE',self.exclude,3)

			self.qdelfix = self.ch.read_namelist_str(line,'QDELFIX',self.epslon,2)
			self.use_raw_prof = self.ch.read_namelist_str(line,'USERAW',self.use_raw_prof,4)
			self.ti_file = self.ch.read_namelist_str(line,'TIFILE',self.ti_file,3)


		f4.close()
		
		self.adjust_variable()

		return	
		
	
	def initialise_var(self):

		self.kbm_coef1 = 0.076
		self.kbm_coef2 = 0.5
		self.perim = 0.
		self.node_list = node_default
		self.moden = '5,7,10,15'
#		self.mode_n = np.array(self.mode_n,dtype=int)
#		self.mode_n = [3,5,7,10]
		self.run_stab = 'mishka'
		self.vloop_ext = 0.
		
		self.highq = False
		self.mis_dir = mis_exec
		self.qdelfix = 0.3
		self.gridn = 301
		self.psis = 0.75
		self.xr1 = 1.0 
		self.sig1 = 0.05
		self.xr2 = 0.9
		self.sig2 = 0.1
		self.ias = True

		self.elite_dir = elite_dir
		self.use_comp = False
		self.ngride = 2000
		self.psise = 0.55
		self.ndist = 50
		self.nqcrit = 27.7

		self.neped = 1.5
		self.nesep = 0.3
		self.tesep = 0.2
		self.w_min = 0.03
		self.w_max = 0.06
		self.w_num = 15

		self.at1 = 2.0
		self.an1 = 1.0
		self.alpt1 = 1.2
		self.alpt2 = 1.4
		self.alpn1 = 1.1
		self.alpn2 = 1.1

		self.ncrit = 1

		self.teped = 0.
		self.tewidth = 0.

		self.grcrit1 = 0.03
		self.grcrit2 = 0.25

		self.epslon = 1.e-8

		self.exclude = ''

		self.use_raw_prof = False

		self.ti_file = None

		return

	def adjust_variable(self):

		self.run_stab = self.run_stab.lower()
		self.w = np.linspace(self.w_min,self.w_max,self.w_num)
		moden = self.moden.split(',')
		self.mode_n = np.array(moden,dtype='int')
		self.node_list2 = self.node_list
		if node_machine.lower() == 'fusma':	self.node_list = make_nodelist(self.node_list)
		if self.exclude == None:
			self.exclude = ''			
		exclude = str(self.exclude)
		exclude = exclude.replace(" ","")
		self.exclude = np.array([])
		if not exclude=='':
			exclude = exclude.split(',')
			for i in range(len(exclude)):
				self.exclude = np.append(self.exclude,int(exclude[i]))
		return
		
	def __init__(self,inputfile='eped_opt',clearf=False):

		self.currdir = os.getcwd()
		self.ch = ch.chease(inputfile)
		self.inputfile = inputfile
		self.initialise_var()
		self.read_namelist(inputfile)
		if clearf:
			return

		if not(self.ch.eqdsk_name.lower() == 'none'):
			self.ch.use_param_shape = True
			self.ch.load_eqdsk()

		self.get_perimeter()

		return

	if __name__ == "__main__":

		import eped

		plot_only = False
		make_equ  = False
		make_ja   = False
		draw_ja   = False
		clear_f   = False
		inputfile = None

		try:
			if (sys.argv[2] == 'plot'):
				plot_only = True
			if (sys.argv[2] == 'make'):
				make_equ = True
			if (sys.argv[2] == 'jastab'):				
				make_ja = True
			if (sys.argv[2] == 'japlot'):
				draw_ja = True				
			if (sys.argv[2] == 'clear'):
				clear_f = True

		except:
			pass

		try:
			inputfile = sys.argv[1]
		except:
			pass

		if (inputfile==None):
			print('Command should include [input_file] [options]')
			print('Available options: eped(default) plot make jastab japlot clear')
			exit()

		copyfile(inputfile,'input_temp')
		f = open('input_temp','r')
		f2 = open(inputfile,'w')
		while True:
			line = f.readline()
			if not line: break
			line2 = line.lower()
			if (line2.find('nideal') > -1):
				line = 'NIDEAL = 8\n'
			if (line2.find('chease_kinetic_file') > -1):
				line = 'chease_kinetic_file = chease_eped \n' 
			f2.write(line)

		f.close()
		f2.close()
		os.remove('input_temp')
	
		sim = eped.eped(inputfile,clear_f)

		sim.plot_only = plot_only
		sim.make_equ = make_equ
		sim.make_ja  = make_ja

		if not (plot_only or make_equ or make_ja or draw_ja or clear_f):
			sim.run_eped()

		if not (clear_f):
			print('NCUTOFF = %i'%sim.ncrit)
			sim.collect_data()

		if (make_ja):
			sim.prepare_ja()

		if (draw_ja):
			sim.draw_ja()

		if (clear_f):
			sim.clear_files()
