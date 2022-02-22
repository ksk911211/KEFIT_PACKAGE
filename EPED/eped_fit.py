#!/usr/local/anaconda3/bin/python3
import os, sys
import numpy as np
import time
from scipy.interpolate import interp1d, interp2d
import matplotlib.pyplot as plt
import ch_tool
import fit_checkfile
import lmfit
import eqdsk
import copy

class eped_ftool:
						
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
		
		datf = interp1d(psink,datk,'cubic')

		dat = datf(psind)
		
		return (dat)
		
	def read_kinprof_kprofile(self):
	
		self.psink = np.linspace(0,1.0,401)
		self.dat_numk = len(self.psink)

		if not (self.te_file == None):
			self.tek = self.read_kinprof_kprofile_fun(self.te_file)
			if (self.tek[2] > 1.e2):
				self.tek = self.tek/1.e3
		if not (self.ne_file == None):
			self.nek = self.read_kinprof_kprofile_fun(self.ne_file)
			if (self.nek[2] > 1.5e1):
				self.nek = self.nek / 1.e1
		if not (self.ti_file == None):
			self.tik = self.read_kinprof_kprofile_fun(self.ti_file)
			if (self.tik[2] > 1.e2):
				self.tik = self.tik / 1.e3

		return

	def read_raw_kinprof_fun(self,filename):
	
		f4 = open(filename,'r')
		line = f4.readline()
		linec = 0
		while True:
			line = f4.readline()
			if not line: break
			if (float(line.split()[0]) > 0.):
				linec = linec + 1
		f4.close
		
		datR = np.zeros(linec)
		datZ = np.zeros(linec)
		datX = np.zeros(linec)
		datP = np.zeros(linec)
		datS = np.zeros(linec)
	
		f4 = open(filename,'r')
		line = f4.readline()
		i = 0
		while (i < linec):
			line = f4.readline().split()
			if (float(line[0]) > 0.):
				datR[i] = float(line[0])
				datZ[i] = float(line[1])	
				datX[i] = float(line[2])
				datS[i] = float(line[3])

				datP[i] = self.psif(datR[i],datZ[i])
				i = i + 1

		for i in range(linec):
			if not (datX[0] == 0.0):	
				if (datX[i] > 3.0*datX[0] and datX[0] < 2.e3):
					datX[i] = 0.0
		f4.close()
				
		self.datz = datZ[0]		

		return (datP, datX, datS)
		
	def read_raw_kinprof(self):

		if not (self.te_dat_file == None):
			self.te_datx,self.te_daty, self.te_dats = self.read_raw_kinprof_fun(self.te_dat_file)
			self.te_datz = self.datz
			if (self.te_daty[2] > 1.e2):
				self.te_daty = self.te_daty / 1.e3
				self.te_dats = self.te_dats / 1.e3
			self.te_datx2, self.te_daty2, self.te_dats2 = self.make_avg_val(self.te_datx,self.te_daty,self.te_dats,self.te_sep,self.fixed_sep_te,False,'te')
		else:
			self.te_datx2 = np.copy(self.psink)
			self.te_daty2 = np.zeros(shape=(len(self.tek),2))
			self.te_daty2[:,0] = np.copy(self.tek)
			self.te_dats2 = np.ones(len(self.tek))
		if not (self.ti_dat_file == None):
			self.ti_datx,self.ti_daty, self.ti_dats = self.read_raw_kinprof_fun(self.ti_dat_file)
			self.vi_datz = self.datz
			if (self.ti_daty[2] > 1.e2):
				self.ti_daty = self.ti_daty / 1.e3
				self.ti_dats = self.ti_dats / 1.e3
			self.ti_datx2, self.ti_daty2, self.ti_dats2 = self.make_avg_val(self.ti_datx,self.ti_daty,self.ti_dats,self.ti_sep,self.fixed_sep_ti,False,'ti')
		else:
			self.ti_datx2 = np.copy(self.psink)
			self.ti_daty2 = np.zeros(shape=(len(self.tik),2))
			self.ti_daty2[:,0] = np.copy(self.tik)
			self.ti_dats2 = np.ones(len(self.tik))
			
		if not (self.ne_dat_file == None):
			self.ne_datx,self.ne_daty, self.ne_dats = self.read_raw_kinprof_fun(self.ne_dat_file)
			self.ne_datz = self.datz
			if (self.ne_daty[2] > 1.5e1):
				self.ne_daty = self.ne_daty / 1.e1
				self.ne_dats = self.ne_dats / 1.e1
			self.ne_datx2, self.ne_daty2, self.ne_dats2 = self.make_avg_val(self.ne_datx,self.ne_daty,self.ne_dats,self.ne_sep,self.fixed_sep_ne,False,'ne')
		else:
			self.ne_datx2 = np.copy(self.psink)
			self.ne_daty2 = np.zeros(shape=(len(self.nek),2))
			self.ne_daty2[:,0] = np.copy(self.nek)
			self.ne_dats2 = np.ones(len(self.nek))
		return

	def modify_profile(self):

		f = open('PROFILES/chease_kinprof_extended','r')
		dat_num = int(float(f.readline()))
		line = f.readline().split()

		te_datz = float(line[0])
		ne_datz = float(line[1])
		ti_datz = float(line[2])

		psi = np.zeros(dat_num)
		R1 = np.zeros(dat_num)
		R2 = np.zeros(dat_num)
		R3 = np.zeros(dat_num)
		R4 = np.zeros(dat_num)

		te = np.zeros(dat_num)
		ne = np.zeros(dat_num)
		ti = np.zeros(dat_num)

		for i in range(dat_num):
			line = f.readline().split()
			psi[i] = float(line[0])
			R1[i] = float(line[5])
			R2[i] = float(line[6])
			R3[i] = float(line[7])
			R4[i] = float(line[8])
			te[i] = float(line[1])
			ne[i] = float(line[2])
			ti[i] = float(line[3])
		f.close()

		te_psi = self.psif(R1,te_datz)
		ne_psi = self.psif(R2,ne_datz)
		ti_psi = self.psif(R3,ti_datz)

		tefit = self.core_extend(te_psi,te)
		nefit = self.core_extend(ne_psi,ne)
		tifit = self.core_extend(ti_psi,ti)

		try: os.mkdir ('PROFILES')
		except: pass

		f = open('PROFILES/chease_kinprof_mod','w')
		f.write('%i\n'%len(tefit))
		f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.zeff,self.zimp,self.amain,self.aimp))

		for i in range(len(tefit)):
			psit = float(i / (len(tefit)-1))
			f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(psit,tefit[i],nefit[i],tifit[i],nefit[i]*(1.0 - (self.zeff-1.0)/self.zimp)))
		f.close()


		return

	def core_extend(self,psi,dat):

		psi_out = np.linspace(0,1.,401)

		if (psi[0] > min(psi)):
			for i in range(1,len(psi)-3):
				if ((psi[i]-psi[i-1])*(psi[i+1]-psi[i])<=0.0):break

			a =np.array([-psi[i-1],psi[i],psi[i+1]])
			b =np.array([dat[i-1],dat[i],dat[i+1]])

			fitt = interp1d(a,b,'cubic')
			psit = np.zeros(len(psi)-i+2)
			datt = np.copy(psit)

			for j in range(len(psit)-1):
				psit[j+1] = psi[i-1+j]
				datt[j+1] = dat[i-1+j]
			datt[0] = fitt(0)


		else:  
			psitt = np.zeros(len(psi)+1)
			datt = np.zeros(len(psi)+1)

			for j in range(len(psi)):
				psitt[j+1] = psi[j]
				datt[j+1] = dat[j]
			fittp = np.polyfit(psi[0:3],dat[0:3],2)
			fitt = np.poly1d(fittp)

			datt[0] = fitt(0)
			if (datt[0] < datt[1]):
					datt[0] = datt[1]

		fit = interp1d(psitt,datt,'cubic')			

		return fit(psi_out)


	def scale_density(self,type=1,print_only=False):

		if (type ==1):
			nef = interp1d(self.psink,self.nek,'cubic')
		elif (type ==2):
			nef = interp1d(self.psin2,self.ne_fit2,'cubic')

		NN = nef(self.psi_den)	

		
		line_sum = np.trapz(NN,x=self.R_den)
		line_avg = line_sum /1.9 * 2.0;

		if (print_only):
			print('>>> Line average density(horizontal) = %f 10(19)/m3 '%(line_avg))
			return
	
		scaled_density = self.target_density / line_avg

		if not (self.scaled_density == scaled_density):
			if (self.use_density_scale):
				print('>>> Line average density(horizontal) = %f >>> %f 10(19)/m3 '%(line_avg,self.target_density))
				self.scaled_density = scaled_density

		return
		
	def residual(self,p,x,y,s,fun,oli_cut):

		yf = fun(x,**p)
		yt = (yf - y)**2 / (s+0.1) ** 2  

		if (oli_cut == 0.0):
				oli_cut2 = 90
		else:
				oli_cut2 = oli_cut * 100

		yt[np.percentile(yt,oli_cut2)<yt] = 0

		return (yt)

	def lmfit_init(self):

		self.param = dict()

		for i in ['te','ne','ti','vt']:
			self.param[i] = dict()
			self.param[i]['vary'] = [True]*6
			self.param[i]['val'] = np.zeros(6)
			self.param[i]['min'] = np.zeros(6)
			self.param[i]['max'] = np.zeros(6)
			for j in range(6):
				self.lmfit_set_param(self.param[i],j,np.nan,np.nan,np.nan)
		return		

	def lmfit_set_param(self,param,ind,val,minv,maxv,vary=True):

		ind = ind - 1

		param['vary'][ind] = vary

		if not (val==-1):
			param['val'][ind] = val		
		if not (minv==-1):
			param['min'][ind] = minv
		if not (maxv==-1):
			param['max'][ind] = maxv

		nan_val = np.isnan(param['val'][ind])
		nan_min = np.isnan(param['min'][ind])
		nan_max = np.isnan(param['max'][ind])

		if not nan_min and not nan_max:

			if param['min'][ind] == param['max'][ind]:
				if minv == -1:
					minv = param['max'][ind] * 0.8
					param['min'][ind] = minv
				else:
					maxv = param['min'][ind] * 1.2
					param['max'][ind] = maxv

		if nan_val:

			if not nan_min and nan_max:
				param['val'][ind] =param['min'][ind] * 1.1
			elif not nan_min and not nan_max:
				param['val'][ind] = 0.5* (param['max'][ind]+param['min'][ind])
			elif nan_min and not nan_max:
				param['val'][ind] = 0.9 * param['max'][ind]

		else:

			if not nan_min and nan_max:
				if param['val'][ind] <=  param['min'][ind]:
					if val == -1:
						param['val'][ind] = param['min'][ind] * 1.1
					else:
						param['min'][ind] = param['val'][ind] * 0.9			

			elif not nan_min and not nan_max:
				if param['val'][ind] <=  param['min'][ind]:
					if val == -1:
						param['val'][ind] = (param['min'][ind] + param['max'][ind]) * 0.5
					else:
						param['min'][ind] = param['val'][ind] * 0.9

				if param['val'][ind] >=  param['max'][ind]:
					if val == -1:
						param['val'][ind] = (param['min'][ind] + param['max'][ind]) * 0.5
					else:
						param['max'][ind] = param['val'][ind] * 1.1

			elif nan_min and not nan_max:
				if param['val'][ind] >=  param['max'][ind]:
					if val == -1:
						param['val'][ind] = param['max'][ind] * 0.9
					else:
						param['max'][ind] = param['val'][ind] * 1.1

		return

	def lmfit_load_param(self,param,fit_types,filename='fit_opt_param'):

		try:
			f = open(filename,'r')
			print('>>> Load lmfit params...')
		except:
			print('>>> No saved lmfit params...')
			return
		k = -1
		for i in ['te','ne','ti','vt']:
			line = f.readline()
			line = line.split('=')
			fit_type = int(line[1])
			k = k + 1
			if (fit_type == fit_types[k]):
				for j in range(6):
					paramt = param[i]
					line = f.readline().split()
					var1 = line[0]
					var2 = line[1]
					var3 = line[2]
					var4 = line[3]

					if var1.lower() == 'inf': var1 = np.nan
					if var2.lower() == '-inf': var2 = np.nan
					if var3.lower() == 'inf': var3 = np.nan
					if var4.lower() == 'true': var4 = True
					elif var4.lower() == 'false': var4 = False
					else: print('>>> no logic var')
	
					paramt['val'][j] = float(var1)
					paramt['min'][j] = float(var2)
					paramt['max'][j] = float(var3)
					paramt['vary'][j] = var4
			else:
				for j in range(6):
					line = f.readline()

		f.close()

		return		

	def lmfit_write_param(self,param,filename='fit_opt_param'):

		fit_type = dict()
		fit_type['te'] = self.fit_type_te
		fit_type['ne'] = self.fit_type_ne
		fit_type['ti'] = self.fit_type_ti
		fit_type['vt'] = self.fit_type_vt

		f = open(filename,'w')
		for i in ['te','ne','ti','vt']:
			line = i + '_fit_variables(val,min,max,vary), fit_type = %i\n'%fit_type[i]
			f.write(line)
			for j in range(6):
				paramt = param[i]
				var1 = paramt['val'][j]
				var2 = paramt['min'][j]
				var3 = paramt['max'][j]
				var4 = paramt['vary'][j]

				if var1 == np.nan: var1 = np.inf
				if var2 == np.nan: var2 = -np.inf
				if var3 == np.nan: var3 = np.inf

				line = '%9.6f\t%9.6f\t%9.6f\t%s\n'%(var1,var2,var3,var4)
				f.write(line)

		f.close()

		return

	def lmfit_put_param(self,param,p,varn):

		for i in range(varn):

			var_name = 'a%i'%(i+1)
			args = dict()
			if not np.isnan(param['val'][i]):
				args['value'] = param['val'][i]

			if not np.isnan(param['min'][i]):
				args['min'] = param['min'][i]

			if not np.isnan(param['max'][i]):
				args['max'] = param['max'][i]


			args['vary'] = param['vary'][i]

			p.add(var_name,**args)

		return	

	def lmfit_init_param(self,fit_type,flag):

		if (flag == 'te'):
			core = self.te_core
			height = self.te_height
		elif (flag == 'ne'):
			core = self.ne_core
			height = self.ne_height
		elif (flag == 'ti'):
			core = self.ti_core
			height = self.ti_height
		elif (flag == 'vt'):
			core = self.vt_core
			height = self.vt_height

		if (flag == 'te' or flag == 'ti'):
			self.lmfit_set_param(self.param[flag],1,0.1,0.1,0.6,True)
			self.lmfit_set_param(self.param[flag],2,0.5,0.15,5.0,True)
		elif (flag == 'ne'):
			self.lmfit_set_param(self.param[flag],1,0.5,0.3,1.5,True)
			self.lmfit_set_param(self.param[flag],2,1.0,0.15,5.0,True)
		elif (flag == 'vt'):
			self.lmfit_set_param(self.param[flag],1,50,0.0,-1,True)
			self.lmfit_set_param(self.param[flag],2,100,0.0,-1,True)

		self.lmfit_set_param(self.param[flag],3,0.05,0.02,0.2,True)

		if (flag == 'vt'):
			self.lmfit_set_param(self.param[flag],4,200,0.0,-1,True)
		else:
			self.lmfit_set_param(self.param[flag],4,3.0,0.2,10.0,True)

		self.lmfit_set_param(self.param[flag],5,1.1,1.01,7.0,True)
		self.lmfit_set_param(self.param[flag],6,2.0,1.01,7.0,True)

		return

	def lmfit_init_params(self):

		self.lmfit_init()
		self.lmfit_init_param(self.fit_type_te,'te')
		self.lmfit_init_param(self.fit_type_ne,'ne')
		self.lmfit_init_param(self.fit_type_ti,'ti')
		self.lmfit_init_param(self.fit_type_vt,'vt')

		return

	def lmfit_renew_param(self,fit_type,fixed_width,width,pedmid,fixed_sep,sep,flag): #CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		if (flag == 'te'):
			core = self.te_core
			height = self.te_height
		elif (flag == 'ne'):
			core = self.ne_core
			height = self.ne_height
		elif (flag == 'ti'):
			core = self.ti_core
			height = self.ti_height
		elif (flag == 'vt'):
			core = self.vt_core
			height = self.vt_height
	
		if (fixed_sep):
			self.lmfit_set_param(self.param[flag],1,abs(sep),abs(sep),-1,False)
		else:
			if (sep < 0.0):
				self.lmfit_set_param(self.param[flag],1,-1,abs(sep),-1,True)
			elif (sep > 0.0):
				self.lmfit_set_param(self.param[flag],1,-1,-1,abs(sep),True)

		if (height < 0.0):
			self.lmfit_set_param(self.param[flag],2,-1,abs(height),-1,True)
		elif (height > 0.0):
			self.lmfit_set_param(self.param[flag],2,-1,-1,abs(height),True)

		if (fixed_width):
			self.lmfit_set_param(self.param[flag],3,abs(width),abs(width),-1,False)
		else:
			if (width < 0.0):
				self.lmfit_set_param(self.param[flag],3,-1,abs(width),-1,True)
			elif (width > 0.0):
				self.lmfit_set_param(self.param[flag],3,-1,-1,abs(width),True)

		if (core < 0.0):
			self.lmfit_set_param(self.param[flag],4,-1,abs(core),-1,True)
		elif (core > 0.0):
			self.lmfit_set_param(self.param[flag],4,-1,-1,abs(core),True)

		return	
		
	def lmfit_params(self,fit_type,fixed_width,width,pedmid,fixed_sep,sep,flag):

		varn = 6
		
		p = lmfit.Parameters()

		self.lmfit_renew_param(fit_type,fixed_width,width,pedmid,fixed_sep,sep,flag)	
		self.lmfit_put_param(self.param[flag],p,varn)

		return (p)

	def lmfit_params2array(self,p,fit_type):
	
		varn = 6
		
		dat = np.zeros(varn)
		
		for i in range(varn):
		
			s = 'a%i'%(i+1)
			dat[i] = p[s].value
		return dat
		
	def lmfit_fit(self,datx,daty,dats,excn,fit_type,fixed_width,width,pedmid,fixed_sep,sep,psi_end,use_std,s0,use_avg,use_oli,oli_cut,oli_cutn,use_raw,use_raws,flag):

		itime = time.time()

		p = self.lmfit_params(fit_type,fixed_width,width,pedmid,fixed_sep,sep,flag)

		func = self.eped_prof
	
		if (flag.lower() == 'ne'):
			if (self.target_density > 0.):
				var_crit = self.target_density * 1.5
			else:
				var_crit = 3.
		elif (flag.lower() == 'te' or flag.lower() == 'ti'):
			var_crit = 0.8
		elif (flag.lower() == 'vt'):
			var_crit = 1.5e2
		
		if (use_raw):
			datx0 = np.copy(datx)
			daty0 = np.copy(daty)
			dats0 = np.copy(dats)

			for i in range(len(datx0)):
				try:
					if (datx0[0] == datx0[i+10]): break
				except:
					i = 0.
					break

			if not (i==0.0):
				len1 = i+10
				len2 = int(len(daty0)/len1)
			else:
				len1 = len(datx0)
				len2 = 1

			if not (len2 == 1):
				
				for k in range(len1+1):
					if (datx0[k] == datx0[k+1]):break

				if not( k == len1):
					for i in range(len2):
						ratio = daty0[len1*i+k]/daty0[len1*i+k+1]
						if (daty0[len1*i+k] == 0.0 or daty0[len1*i+k+1] == 0.0):
							print('>>> Core/Edge calibration error in channel %s, index# = %i !'%(flag,len1*i+k))
							ratio = 1.

						for j in range(k+1,len1):
							daty0[len1*i+j] = ratio*daty0[len1*i+j]
							dats0[len1*i+j] = ratio*dats0[len1*i+j]
			else:
				k = len1 +1
										
			dat0 = np.zeros(shape=(len1,len2))
			dat1 = np.zeros(shape=(len1,len2))
			for i in range(len2):
				for j in range(len1):
					dat0[j,i] = daty0[len1*i+j]
					dat1[j,i] = dats0[len1*i+j]

			dats = np.copy(daty)

			for i in range(len1):
				sum = []
				for j in range(len2):
					if not (datx[i] > 1.0 and dat0[i,j] > var_crit):
						sum.append(dat0[i,j])
				sum = np.array(sum,dtype=float)
				if not (len(sum) == 0):
					stds = np.std(sum)
					if (stds >= 0):
						dats[i] = stds
					else:
						dats[i] = 0.1
				else:
					dats[i] = 0.1

				for j in range(len2):
					dats[len1*j+i] = dats[i]

			excn = np.array(excn,dtype=int)
			datx2 = []
			daty2 = []
			dats2 = []

			if not (use_avg):
				for i in range(len1):
					chk = excn == i
					if not(chk.any()):
						for j in range(len2):
							if not (datx[i] > 1.0 and dat0[i,j] > var_crit):
								if (np.sum(abs(dat0[i,:]))>0.0 and not i == k+1):
									daty2.append(dat0[i,j])
									datx2.append(datx[i])
									if use_raws:
										dats2.append(dats0[len1*j+i])
									else:
										dats2.append(dats[i])

				if (abs(sep) > 0.0):
						for i in range(len2):
							datx2.append(1.0)
							daty2.append(abs(sep))
							if fixed_sep:
								dats2.append(0.1)
							else:
								dats2.append(abs(sep)/2.)

			else:
				for i in range(len1):
					chk = excn == i
					if not (chk.any()):
							sum = []
							sums = []
							for j in range(len2):
								if not (datx[i] > 1.0 and dat0[i,j] > var_crit):
									if (np.sum(abs(dat0[i,:]))>0.0 and not i == k+1):
										sum.append(dat0[i,j])
										sums.append(dat1[i,j])

							sum = np.array(sum,dtype=float)
							sums = np.array(sums,dtype=float)
							if not (len(sum) == 0):
								daty2.append(np.mean(sum))
								datx2.append(datx[i])
								if use_raws:
									stds = np.sqrt(np.sum(sums))/len(sums)
								else:
									stds = np.std(sum)
								dats2.append(stds)

				if (abs(sep) > 0.0):
					datx2.append(1.0)
					daty2.append(abs(sep))
					if fixed_sep:
						dats2.append(0.1)
					else:
						dats2.append(abs(sep)/2.)

			daty2 = np.array(daty2,dtype=float)
			datx2 = np.array(datx2,dtype=float)
			dats2 = np.array(dats2,dtype=float)

			if not (psi_end == 0.0):
				ind = np.where(datx2 <= psi_end)
				daty2 = daty2[ind]
				datx2 = datx2[ind]
				dats2 = dats2[ind]
								
		else:
				datx2 = np.copy(datx)
				daty2 = np.copy(daty)
				dats2 = np.ones(len(datx))
			
		if not(use_std):
			dats2 = np.ones(len(datx2))

		datx3 = np.copy(datx2)	
		daty3 = np.copy(daty2)
		dats3 = np.copy(dats2)

		if (use_oli and not use_avg and use_raw):
			for i in range(oli_cutn):
				result = lmfit.minimize(self.residual,p,args=(datx3,daty3,dats3,func,oli_cut))
				
				p = result.params
		
				yt = (daty3 - func(datx3,**result.params))**2/(dats3+0.1)**2.0
				if (oli_cut == 0.0):
					oli_cut2 = 90
				else:
					oli_cut2 = min(99,oli_cut*100. + 0.05)

				yn = np.percentile(yt,oli_cut2)>yt

				datx3 = datx3[yn]
				daty3 = daty3[yn]
				dats3 = dats3[yn]
			print('=----------------------------------------------------------------------------=')	
			print('>>> Fitting result for %s, FIT_FUNC = %i'%(flag.upper(),fit_type))			
			print('>>> XI2 = %6.3f, RED_XI2 = %6.4e'%(result.chisqr,result.redchi))
			print('>>> Elapsed time %6.3f(s) '%(time.time()-itime))
			print('=----------------------------------------------------------------------------=')
			result.params.pretty_print()
			print('=----------------------------------------------------------------------------=')

			self.oli_in = np.zeros(shape=(len(daty3),3))
			self.oli_in[:,0] = np.copy(datx3)
			self.oli_in[:,1] = np.copy(daty3)
			self.oli_in[:,2] = np.copy(dats3)
		else:
			result = lmfit.minimize(self.residual,p,args=(datx2,daty2,dats2,func,oli_cut))
			print('=----------------------------------------------------------------------------=')
			print('>>> Fitting result for %s, FIT_FUNC = %i'%(flag.upper(),fit_type))
			print('>>> XI2 = %6.3f, RED_XI2 = %6.4e'%(result.chisqr,result.redchi))
			print('>>> Elapsed time %6.3f(s) '%(time.time()-itime))
			print('=----------------------------------------------------------------------------=')
			result.params.pretty_print()
			print('=----------------------------------------------------------------------------=')

		self.popt = self.lmfit_params2array(result.params,fit_type)
		self.prof1 = func(self.psin1,**result.params)
		self.prof2 = func(self.psin2,**result.params)
		
		self.datx2 = np.copy(datx2)
		self.daty2 = np.copy(daty2)
		self.dats2 = np.copy(dats2)

		try:
			self.chi = result.chisqr
		except:
			self.chi = 1.0
		
		return self.prof2
		
	def eped_prof(self, x, a1, a2, a3, a4, a5, a6):

		y = a1
		y = y + a2*(np.tanh(1) - np.tanh((x - 1.0 + 0.5*(abs(abs(a3)+1.e-7)))/(abs(abs(a3)+1.e-7))*2.0))
		yt = (x/(1.0-abs(abs(a3)+1.e-7))) ** (abs(a5))
		y = y + a4 * (abs((1-yt)) **(abs(a6))) * 0.5 * (1.0 + np.sign(1.0-abs(abs(a3)+1.e-7)-x))
		
		return y

	def make_avg_val(self,datx,daty,dats,sep_val,fixed_sep,outliar=False,flag='ne'):

		len1 = len(datx)
		for i in range(len1-1):
			if (datx[0] == datx[i+1]):
				break
		len2 = i+1
		if (len2 == (len1-1)):
			len2 = len1

		len3 = int(len1/len2)

		datyt = np.zeros(shape=(len2,len3))
		datst = np.zeros(shape=(len2,len3))

		datx2 = np.zeros(len2)
			
		for j in range(len2):
			datx2[j] = datx[j]
		j = 0
		k = 0
			
		for i in range(len1):
			datyt[j,k] = daty[i]
			datst[j,k] = dats[i]

			j = j + 1
			if ( j == len2):
				j = 0
				k = k + 1
		var_crit = 0

		if (flag.lower() == 'ne'):
			if (self.target_density > 0.):
				var_crit = self.target_density*1.5
			else:
				var_crit = 4.
		elif (flag.lower() == 'te' or flag.lower() == 'ti'):
			var_crit = 1.5
		elif (flag.lower() == 'vt'):
			var_crit = 1.5e2
	
		dats3 = []
		daty3 = []
		datx3 = []

		for i in range(len2):
			sum = []
			sums = []
			for j in range(len3):
				if (datx2[i] > 1.0 and datyt[i,j] > var_crit):
					pass
				else:
					sum = np.append(sum,datyt[i,j])
					sums = np.append(sums,datst[i,j])
			
			if (len(sum) > 0.):
				mean = np.mean(sum)
				if (np.isnan(mean)):
					mean = 0.e0

				std = np.std(sum)
				if (np.isnan(std)):
					std = 0.1
				if std < 0.0: std = 0.1
		
				sums = sums ** 2
				snew = np.sum(sums)
				snew = np.sqrt(snew)/len(sum)

				daty3.append([mean,std,len(sum)])
				datx3 = np.append(datx3,datx2[i])
				dats3 = np.append(dats3,snew)

		if not (sep_val == 0.0):
			if (fixed_sep):
				daty3.append([abs(sep_val),0.1,1])
				dats3 = np.append(dats3,0.1)
			else:
				daty3.append([abs(sep_val),abs(sep_val)/2.,1])
				dats3 = np.append(dats3,abs(sep_val)/2.)

			datx3 = np.append(datx3,1.0)

		daty3 = np.array(daty3)

		len4 = len(datx3)
		index = len4-2
		ind = []
		for i in range(len4-1):
			if (datx3[i] == datx3[i+1]):	index = i
			else:
				ind.append(i)
		ind.append(i+1)
		
		if (index < len4-2): ratio = daty3[index,0]/daty3[index+1,0]
		for i in range(len4-2-i):
			daty3[i+index+1] = daty3[i+index+1] * ratio
			dats3[i+index+1] = dats3[i+index+1] * ratio

		datx3 = datx3[ind]
		daty3 = daty3[ind]
		dats3 = dats3[ind]

		return(datx3,daty3,dats3)
				
	def draw_plot_part(self,fig,flag,raw_fit,use_avg,use_oli,use_std,use_raws,fit_type,raw_file,datx,daty,dats,datx2,daty2,dats2,excn,fit_file,fity,oli_in,fity2,popt,xi):

		if (flag == 'te'):
			title = 'TE [keV]'
		elif (flag == 'ne'):
			title = 'NE [10(19)/m3]'
		elif (flag == 'ti'):
			title = 'TI [keV]'
		elif (flag == 'vt'):
			title = 'VT [km/s]'
		
		excn = np.array(excn,dtype=int)
		fig.set_title(title)

		fig.set_xlabel('$\psi_n$ [a.u]')

		fig.set_ylabel(title)

		plegend = []
		llegend = []

		if (fit_type == 4 or fit_type ==6):
			fig.text(max(self.psink),max(fity2),'$W_{ped}$=%4.3f'%popt[2])
		if (fit_type == 2):
			fig.text(max(self.psink),max(fity2),'$W_{ped}$=%4.3f'%(popt[2]*4.))


		if not (fit_file == None):
			line1, = fig.plot(self.psink,fity,'blue',linestyle= '--')
			plegend.append('Pre Fitted prof')
			llegend.append(line1)


		line2, = fig.plot(self.psin2,fity2,'r')
		plegend.append('New Fitted prof')
		llegend.append(line2)
		
		if (flag.lower() == 'ne'):
			if (self.target_density > 0.0):
				var_crit = self.target_density * 3.0
			else:
				var_crit = 10.0
		elif (flag.lower() == 'te' or flag.lower() == 'ti'):
			var_crit = 2.0
		elif (flag.lower() == 'vt'):
			var_crit = 2.e2


		if not (raw_file == None):
			plegend.append('Raw data')
			
			if not (use_avg):
				for i in range(len(datx)):
					try:
						if (datx[0] == datx[i+10]): break
					except:
						i = len(datx) - 1

				if not (i==len(datx)-1):
					len1 = i+10
					len2 = int(len(daty)/len1)
				else:
					len1 = len(datx)
					len2 = 1
					
				if not (len2 == 1):
					for k in range(len1+1):
						if (datx[k] == datx[k+1]):break

					if not( k == len1):
						for i in range(len2):
							ratio = daty[len1*i+k]/daty[len1*i+k+1]
							if (daty[len1*i+k] == 0.0 or daty[len1*i+k+1]):
								ratio = 1.
								print('Core/Edge calibration error in channel %s, index# = %i !'%(flag,len1*i+k))
							for j in range(k+1,len1):
								daty[len1*i+j] = ratio*daty[len1*i+j]
								dats[len1*i+j] = ratio*dats[len1*i+j]

				else:
					k = len1 + 1
											
				dat0 = np.zeros(shape=(len1,len2))
				dat1 = np.zeros(shape=(len1,len2))
				for i in range(len2):
					for j in range(len1):
						dat0[j,i] = daty[len1*i+j]
						dat1[j,i] = dats[len1*i+j]
	
			if not (use_avg):
				
				ind = np.where(datx > 1.0)
				ind2 = np.where(datx <= 1.0)
				datxt = datx[ind]
				datyt = daty[ind]
				datst = dats[ind]

				ind = np.where(datyt < var_crit)
				if (use_std):
					line3 = fig.errorbar(datxt[ind],datyt[ind], yerr = datst[ind],fmt='x',markersize='5',c='gray',ecolor='gray',capthick=2)
					line3 = fig.errorbar(datx[ind2],daty[ind2], yerr = dats[ind2],fmt='x',markersize='5',c='gray',ecolor='gray',capthick=2)
				else:
					line3 = fig.scatter(datxt[ind],datyt[ind],s=20,marker='x',c='gray')
					line3 = fig.scatter(datx[ind2],daty[ind2],s=20,marker='x',c='gray')

				for i in range(len1):
					fig.text(datx[i],min(np.median(dat0[i,:]),0.7*var_crit),str(i))

				if (raw_fit):
					if (use_std):
						fig.errorbar(datx2,daty2,yerr=dats2,fmt='x',markersize='5',c='g',ecolor='g',capthick=2)
					else:
						fig.scatter(datx2,daty2,marker='x',s=20,c='g')
					if (use_oli):
						ind = np.where(oli_in[:,0] > 1.0)
						ind2 = np.where(oli_in[:,0] <= 1.0)
						oli_intx = oli_in[ind,0]
						oli_inty = oli_in[ind,1]
						oli_ints = oli_in[ind,2]
						ind = np.where(oli_inty < var_crit)
						if (use_std):
							fig.errorbar(oli_intx[ind],oli_inty[ind], yerr = oli_ints[ind],fmt='x',markersize='5',c='orange',ecolor='orange',capthick=2)
							fig.errorbar(oli_in[ind2,0],oli_in[ind2,1], yerr = oli_in[ind2,2],fmt='x',markersize='5',c='orange',ecolor='orange',capthick=2)
						else:
							fig.scatter(oli_intx[ind],oli_inty[ind],s=20,marker='x',c='orange')
							fig.scatter(oli_in[ind2,0],oli_in[ind2,1],s=20,marker='x',c='orange')
			else:
				for i in range(len(datx)):
					chk = excn == i
					if not (chk.any() or not raw_fit ):
						if not (use_std):
							line3 = fig.scatter(datx[i],daty[i],s=20,marker='x',c='g')
						else:
							line3 = fig.errorbar(datx[i],daty[i], yerr = dats[i],fmt='*',markersize='5',c='g',ecolor='g',capthick=2)
					else:
						if not (use_std):
							line3 = fig.scatter(datx[i],daty[i],s=20,marker='x',c='gray')
						else:
							line3 = fig.errorbar(datx[i],daty[i], yerr = dats[i],fmt='*',markersize='5',c='gray',ecolor='gray',capthick=2)

					fig.text(datx[i],daty[i],str(i))
					
			llegend.append(line3)		
					
		line4 = fig.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
		plegend.append('xi2 = %4.2e'%xi)
		llegend.append(line4)
		fig.legend(llegend,plegend)
		
		return

	def draw_plot_part2(self,fig,flag,fit_file,fity,fity2):

		if (flag == 'te'):
			title = 'dTE/dpsin'
		elif (flag == 'ne'):
			title = 'dNE/dpsin'
		elif (flag == 'ti'):
			title = 'dTI/dpsin'
		elif (flag == 'vt'):
			title = 'dVT/dpsin'
		
		fig.set_title(title)

		fig.set_xlabel('$\psi_n$ [a.u]')

		#fig.set_ylabel(title)

		plegend = []
		llegend = []
		deps = 1.e-5

		if not (fit_file == None):

			fitf = interp1d(self.psink,fity,'quadratic')
			fitd = np.copy(fity)
			for i in range(len(fity)-2):
				xx = self.psink[i+1]
				fitd[i+1] = (fitf(xx+deps)-fitf(xx-deps))/2./deps
			fitd[0] = (fitf(self.psink[0]+deps)-fitf(self.psink[0]))/deps
			fitd[-1] = (fitf(self.psink[-1])-fitf(self.psink[-1]-deps))/deps

			line1, = fig.plot(self.psink,fitd,'blue',linestyle= '--')
			plegend.append('Pre Fitted prof')
			llegend.append(line1)

		fitf = interp1d(self.psin2,fity2,'quadratic')
		fitd = np.copy(fity2)
		for i in range(len(fity2)-2):
			xx = self.psin2[i+1]
			fitd[i+1] = (fitf(xx+deps)-fitf(xx-deps))/2./deps
		fitd[0] = (fitf(self.psin2[0]+deps)-fitf(self.psin2[0]))/deps
		fitd[-1] = (fitf(self.psin2[-1])-fitf(self.psin2[-1]-deps))/deps		

		line2, = fig.plot(self.psin2,fitd,'r')
		plegend.append('New Fitted prof')
		llegend.append(line2)
							
		line4 = fig.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
		plegend.append('Separatrix')
		llegend.append(line4)
		fig.legend(llegend,plegend)
		
		return

	def draw_plot(self,fig_ex=None,second=False):

		legend4 = []

		len2 = len(plt.get_fignums())
		if not fig_ex == None:
			fig = fig_ex
			[ax1,ax2, ax3, ax4] = fig.axes
			ax1.cla()
			ax2.cla()
			ax3.cla()
			ax4.cla()

			fig.canvas.set_window_title('Profile FIT')

		else:
			if (len(plt.get_fignums()) == 0):
				fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2,2,figsize=(15,9))
				self.fig = fig
			else:
				fig = self.fig
				[ax1,ax2, ax3, ax4] = fig.axes
				ax1.cla()
				ax2.cla()
				ax3.cla()
				ax4.cla()
		if not (self.no_te):
			if (self.te_raw_fit):
				xx = self.te_datx
				yy = self.te_daty
				xx2 = self.datx2_te
				yy2 = self.daty2_te
				ss = self.te_dats
				if not (self.use_raws_te):
					len1 = len(ss)
					len2 = len(self.te_daty2)
					for i in range(len2):
						for j in range(int(len1/len2)):
							ss[len2*j+i] = self.te_daty2[i,1]
				ss2 = self.dats2_te
				if (self.te_use_avg_dat):
					xx = self.te_datx2
					yy = self.te_daty2[:,0]
					if (self.use_raws_te):
						ss = self.te_dats2
					else:
						ss = self.te_daty2[:,1]
				if not (self.te_file == None):
					fity = self.tek
				else:
					fity = []
			else:
				if not (self.te_dat_file == None):
					xx = self.te_datx
					yy = self.te_daty
					ss = self.te_dats
				else:
					xx = self.te_datx2
					yy = self.te_daty2
					ss = self.te_dats2
				xx2 = self.datx2_te
				yy2 = self.daty2_te
				ss2 = self.dats2_te
				if (self.te_use_avg_dat):
					xx = self.te_datx2
					yy = self.te_daty2[:,0]
					if (self.use_raws_te):
						ss = self.te_dats2
					else:
						ss = self.te_daty2[:,1]
				fity = self.tek
			if not (second):
				self.draw_plot_part(ax1,'te',self.te_raw_fit,self.te_use_avg_dat,self.te_use_oli,self.std_te,self.use_raws_te,self.fit_type_te,self.te_dat_file,xx,yy,ss,xx2,yy2,ss2,self.te_excn,self.te_file,fity,self.oli_in_te,self.te_fit2,self.popt_te,self.te_xi)
			else:
				self.draw_plot_part2(ax1,'te',self.te_file,fity,self.te_fit2)

			ax4.plot(self.psin2,self.te_fit2)
			legend4.append('Te [keV]')
		else:
			ax1.text(0.44,0.5,'No te data',color='r')

		if not (self.no_ne):
			if (self.ne_raw_fit):
				xx = self.ne_datx
				yy = self.ne_daty
				xx2 = self.datx2_ne
				yy2 = self.daty2_ne
				ss = self.ne_dats
				if not (self.use_raws_ne):
					len1 = len(ss)
					len2 = len(self.ne_daty2)
					for i in range(len2):
						for j in range(int(len1/len2)):
							ss[len2*j+i] = self.ne_daty2[i,1]
				ss2 = self.dats2_ne
				if (self.ne_use_avg_dat):
					xx = self.ne_datx2
					yy = self.ne_daty2[:,0]
					if (self.use_raws_ne):
						ss = self.ne_dats2
					else:
						ss = self.ne_daty2[:,1]
				if not (self.ne_file == None):
					fity = self.nek/self.scaled_density
				else:
					fity = []
			else:
				if not (self.ne_dat_file == None):
					xx = self.ne_datx
					yy = self.ne_daty
					ss = self.ne_dats
				else:
					xx = self.ne_datx2
					yy = self.ne_daty2
					ss = self.ne_dats2
				xx2 = self.datx2_ne
				yy2 = self.daty2_ne
				ss2 = self.dats2_ne
				if (self.ne_use_avg_dat):
					xx = self.ne_datx2
					yy = self.ne_daty2[:,0]
					if (self.use_raws_ne):
						ss = self.ne_dats2
					else:
						ss = self.ne_daty2[:,1]
				fity = self.nek/self.scaled_density
			
			if not (second):
				self.draw_plot_part(ax2,'ne',self.ne_raw_fit,self.ne_use_avg_dat,self.ne_use_oli,self.std_ne,self.use_raws_ne,self.fit_type_ne,self.ne_dat_file,xx,yy,ss,xx2,yy2,ss2,self.ne_excn,self.ne_file,fity,self.oli_in_ne,self.ne_fit2/self.scaled_density,self.popt_ne,self.ne_xi)
			else:
				self.draw_plot_part2(ax2,'ne',self.ne_file,fity,self.ne_fit2/self.scaled_density)	

			ax4.plot(self.psin2,self.ne_fit2/self.ne_fit2[0]*self.te_fit2[0]*1.1)
			legend4.append('NE [a.u]')
		else:
			ax2.text(0.44,0.5,'No ne data',color='r')	

		if not (self.no_ti):

			if (self.ti_raw_fit):
				xx = self.ti_datx
				yy = self.ti_daty
				xx2 = self.datx2_ti
				yy2 = self.daty2_ti
				ss = self.ti_dats
				if not (self.use_raws_ti):
					len1 = len(ss)
					len2 = len(self.ti_daty2)
					for i in range(len2):
						for j in range(int(len1/len2)):
							ss[len2*j+i] = self.ti_daty2[i,1]
				ss2 = self.dats2_ti
				if (self.ti_use_avg_dat):
					xx = self.ti_datx2
					yy = self.ti_daty2[:,0]
					if (self.use_raws_ti):
						ss = self.ti_dats2
					else:
						ss = self.ti_daty2[:,1]
				if not (self.ti_file == None):
					fity = self.tik
				else:
					fity = []
			else:
				if not (self.ti_dat_file == None):
					xx = self.ti_datx
					yy = self.ti_daty
					ss = self.ti_dats
				else:
					xx = self.ti_datx2
					yy = self.ti_daty2
					ss = self.ti_dats2
				xx2 = self.datx2_ti
				yy2 = self.daty2_ti
				ss2 = self.dats2_ti
				if (self.ti_use_avg_dat):
					xx = self.ti_datx2
					yy = self.ti_daty2[:,0]
					if (self.use_raws_ti):
						ss = self.ti_dats2
					else:
						ss = self.ti_daty2[:,1]
				fity = self.tik
			
			if not (second):
				self.draw_plot_part(ax3,'ti',self.ti_raw_fit,self.ti_use_avg_dat,self.ti_use_oli,self.std_ti,self.use_raws_ti,self.fit_type_ti,self.ti_dat_file,xx,yy,ss,xx2,yy2,ss2,self.ti_excn,self.ti_file,fity,self.oli_in_ti,self.ti_fit2,self.popt_ti,self.ti_xi)
			else:
				self.draw_plot_part2(ax3,'ti',self.ti_file,fity,self.ti_fit2)

			ax4.plot(self.psin2,self.ti_fit2)
			legend4.append('TI [keV]')				
		else:
			ax3.text(0.44,0.5,'No ti data',color='r')


		ax4.set_xlabel('$\psi_n$ [a.u]')
		ax4.set_ylabel('[a.u]')
		ax4.legend(legend4)

		if (len2 == 0):
			plt.show(block=False)
		else:
			plt.draw()
		return

	def write_kinprof(self):

		try:
			os.mkdir('PROFILES')
		except:
			pass
				
		if not (self.no_ne):	f = open('PROFILES/NE_fit.dat','w')
		if not (self.no_ne):	f.write('--- Fitted profile by FUSMA fittool --- \n Fitted Values: X_Norm_Rho, Psi_Norm, NE[1E18/m3]\n')		
		if not (self.no_te):	f2 = open('PROFILES/TE_fit.dat','w')
		if not (self.no_te):	f2.write('--- Fitted profile by FUSMA fittool --- \n Fitted Values: X_Norm_Rho, Psi_Norm, TE[eV]\n')
		if not (self.no_ti):	f3 = open('PROFILES/TI_fit.dat','w')
		if not (self.no_ti):	f3.write('--- Fitted profile by FUSMA fittool --- \n Fitted Values: X_Norm_Rho, Psi_Norm, TI[eV]\n')

		if not (self.no_ne and self.no_te and self.no_ti):	f5 = open('PROFILES/chease_kinprof_fit','w')	
		if not (self.no_ne and self.no_te and self.no_ti):	f5.write('401\n')
		if not (self.no_ne and self.no_te and self.no_ti):	f5.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.zeff,self.zimp,self.amain,self.aimp))
		
		ppp = np.linspace(0,1,101)
		ppp2 = np.linspace(0,1,401)
		pp1 = interp1d(self.psin2,self.ne_fit2,'cubic')
		pp2 = interp1d(self.psin2,self.te_fit2,'cubic')
		pp3 = interp1d(self.psin2,self.ti_fit2,'cubic')

		for i in range(101):
			psi_temp = self.rho_to_psi(ppp[i])
			psi_temp2 = psi_temp
			if not (self.no_ne):	f.write('%9.6f\t%9.6f\t%9.6f\n'%(ppp[i],psi_temp,pp1(psi_temp2)*1.e1))
			if not (self.no_te):	f2.write('%9.6f\t%9.6f\t%9.6f\n'%(ppp[i],psi_temp,pp2(psi_temp2)*1.e3))
			if not (self.no_ti):	f3.write('%9.6f\t%9.6f\t%9.6f\n'%(ppp[i],psi_temp,pp3(psi_temp2)*1.e3))

		if not (self.no_ne and self.no_te and self.no_ti):
			for i in range(401):
				psint = ppp2[i] 
				f5.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(psint,pp2(ppp2[i]),pp1(ppp2[i]),pp3(ppp2[i]),pp1(ppp2[i])* (1.0 - (self.zeff-1.0)/self.zimp)))

		try:	R1 = self.make_psiRZ_extended(self.psin2,self.te_datz)
		except:	R1.ones(len(self.psin2))
		try:	R2 = self.make_psiRZ_extended(self.psin2,self.ne_datz)
		except:	R1.ones(len(self.psin2))
		try:	R3 = self.make_psiRZ_extended(self.psin2,self.ti_datz)
		except:	R1.ones(len(self.psin2))

		if not (self.no_ne):	f.close()
		if not (self.no_te):	f2.close()
		if not (self.no_ti):	f3.close()
		if not (self.no_ne and self.no_te and self.no_ti):	f5.close()

		dpdr = (1.0 - self.psi_to_rho(0.9999))/ (1.0 - 0.9999)


		if not (self.no_ne):
			f4 = open('PROFILES/scaled_factor','w')
			f4.write('%f\n'%self.scaled_density)
			f4.close()
	
		if (not self.no_ne):		
			if (self.use_density_scale):
				scale = self.scaled_density
			else:
				scale = 1.

			f4 = open('PROFILES/NE_non_scaled.dat','w')
			f4.write('--- Fitted profile by FUSMA fittool --- \n Fitted Values: X_Norm_Rho, Psi_Norm, NE[1E18/m3]\n')
			for i in range(101):
				psi_temp = self.rho_to_psi(ppp[i])
				psi_temp2 = psi_temp
				if self.use_rho:
					psi_temp2 = ppp[i]
				f4.write('%9.6f\t%9.6f\t%9.6f\n'%(ppp[i],psi_temp,pp1(psi_temp2)*1.e1/self.scaled_density))
			f4.close()
			
		f4 = open('PROFILES/chease_eped_mod','w')
		
		pedw = abs(self.popt_ne[2])
		self.popt_ne[0] = self.popt_ne[0]*self.scaled_density
		self.popt_ne[1] = self.popt_ne[1]*self.scaled_density
		self.popt_ne[3] = self.popt_ne[3]*self.scaled_density

		if (self.use_rho):
			popt_temp = np.copy(self.popt_ne)
			self.popt_ne = self.modify_rho_to_psi_eped(popt_temp)#,'ne',True)
			pedw = abs(self.popt_ne[2])

		f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.popt_ne[1],self.popt_ne[0],self.popt_ne[3],1-0.5*pedw,pedw,1-pedw,abs(self.popt_ne[4]),abs(self.popt_ne[5])))
		if (self.use_rho):
			self.popt_ne = np.copy(popt_temp)

		repeat = 1
		if not (self.use_ti_eped):
			repeat = 2
		
		for i in range(repeat):
			pedw = abs(self.popt_te[2])
			if (self.use_rho):
				popt_temp = np.copy(self.popt_te)
				self.popt_te = self.modify_rho_to_psi_eped(popt_temp)#,'Te',True)
				pedw = abs(self.popt_te[2])

			f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.popt_te[1],self.popt_te[0],self.popt_te[3],1-0.5*pedw,pedw,1-pedw,abs(self.popt_te[4]),abs(self.popt_te[5])))
			if (self.use_rho):
				self.popt_te = np.copy(popt_temp)
				
		if (self.use_ti_eped):
			pedw = abs(self.popt_ti[2])
			if (self.use_rho):
				popt_temp = np.copy(self.popt_ti)
				self.popt_ti = self.modify_rho_to_psi_eped(popt_temp)#,'Ti',True)
				pedw = abs(self.popt_ti[2])

			f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.popt_ti[1],self.popt_ti[0],self.popt_ti[3],1-0.5*pedw,pedw,1-pedw,abs(self.popt_ti[4]),abs(self.popt_ti[5])))
			if (self.use_rho):
				self.popt_ti = np.copy(popt_temp)

		f4.close()
		
		return
			
	def read_namelist_str(self,line,var_name,vars,vartype):
	
		line = line.split('!')[0]
	
		try:
			name = line.split('=')[0].split()[0].lower()
		except:
			name = None
		
		if (name == var_name.lower()):
			
			if (vartype < 4):
				try:
					if (vartype == 1):
						var = int(float(line.split('=')[1].split()[0]))
					elif (vartype == 2):
						var = float(line.split('=')[1].split()[0])
					elif (vartype == 3):
						#var = line.split('=')[1].split()[0]
						var = line.split('=')[1].strip()
				except:
					var = vars;
				exist = 1
			else:
				var2 = line.split('=')[1].split()[0].lower()
				
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

		if (vartype ==3):	
			try:
				if (var.lower() == 'none'):
					var = None
			except:
				pass
		
		return (var)

	def read_namelist(self,filename='fit_opt'):
	
		f4 = open(filename,'r')
		
		while True:
		
			line = f4.readline()
			if not line: break
	
			self.psi_end_te = self.read_namelist_str(line,'PSI_END_TE',self.psi_end_te,2)
			self.psi_end_ti = self.read_namelist_str(line,'PSI_END_TI',self.psi_end_ti,2)		
			self.psi_end_ne = self.read_namelist_str(line,'PSI_END_NE',self.psi_end_ne,2)
			self.psi_end_vt = self.read_namelist_str(line,'PSI_END_VT',self.psi_end_vt,2)

			self.te_file = self.read_namelist_str(line,'TE_FILE',self.te_file,3)
			self.ti_file = self.read_namelist_str(line,'TI_FILE',self.ti_file,3)
			self.ne_file = self.read_namelist_str(line,'NE_FILE',self.ne_file,3)
			self.vt_file = self.read_namelist_str(line,'VT_FILE',self.vt_file,3)
			
			self.te_dat_file = self.read_namelist_str(line,'TE_DAT_FILE',self.te_dat_file,3)
			self.ti_dat_file = self.read_namelist_str(line,'TI_DAT_FILE',self.ti_dat_file,3)
			self.ne_dat_file = self.read_namelist_str(line,'NE_DAT_FILE',self.ne_dat_file,3)
			self.vt_dat_file = self.read_namelist_str(line,'VT_DAT_FILE',self.vt_dat_file,3)
			
			self.target_density = self.read_namelist_str(line,'DENSITY_SCALE',self.target_density,2)
			
			self.use_ti_eped = self.read_namelist_str(line,'USE_TI_EPED',self.use_ti_eped,4)
			self.use_ti_width = self.read_namelist_str(line,'USE_TI_WIDTH',self.use_ti_width,4)
		
			self.fixed_sep_ne = self.read_namelist_str(line,'FIX_NE_sep',self.fixed_sep_ne,4)
			self.fixed_sep_te = self.read_namelist_str(line,'FIX_TE_sep',self.fixed_sep_te,4)
			self.fixed_sep_ti = self.read_namelist_str(line,'FIX_TI_sep',self.fixed_sep_ti,4)
			self.fixed_sep_vt = self.read_namelist_str(line,'FIX_VT_sep',self.fixed_sep_vt,4)
	
			self.ne_sep = self.read_namelist_str(line,'NE_sep',self.ne_sep,2)
			self.te_sep = self.read_namelist_str(line,'TE_sep',self.te_sep,2)
			self.ti_sep = self.read_namelist_str(line,'TI_sep',self.ti_sep,2)
			self.vt_sep = self.read_namelist_str(line,'VT_sep',self.vt_sep,2)
			
			self.ped_scan_fit = self.read_namelist_str(line,'Ped_scan_fit',self.ped_scan_fit,4)

			self.fixed_width_ne = self.read_namelist_str(line,'FIX_NE_width',self.fixed_width_ne,4)
			self.fixed_width_te = self.read_namelist_str(line,'FIX_TE_width',self.fixed_width_te,4)
			self.fixed_width_ti = self.read_namelist_str(line,'FIX_TI_width',self.fixed_width_ti,4)
			self.fixed_width_vt = self.read_namelist_str(line,'FIX_VT_width',self.fixed_width_vt,4)

			self.te_width = self.read_namelist_str(line,'TE_width',self.te_width,2)
			self.ti_width = self.read_namelist_str(line,'TI_width',self.ti_width,2)
			self.ne_width = self.read_namelist_str(line,'NE_width',self.ne_width,2)
			self.vt_width = self.read_namelist_str(line,'VT_width',self.vt_width,2)

			self.te_height = self.read_namelist_str(line,'TE_height',self.te_height,2)
			self.ne_height = self.read_namelist_str(line,'NE_height',self.ne_height,2)
			self.ti_height = self.read_namelist_str(line,'TI_height',self.ti_height,2)
			self.vt_height = self.read_namelist_str(line,'VT_height',self.vt_height,2)
		
			self.te_core = self.read_namelist_str(line,'TE_core',self.te_core,2)
			self.ne_core = self.read_namelist_str(line,'NE_core',self.ne_core,2)
			self.ti_core = self.read_namelist_str(line,'TI_core',self.ti_core,2)
			self.vt_core = self.read_namelist_str(line,'VT_core',self.vt_core,2)

			self.te_use_avg_dat = self.read_namelist_str(line,'TE_AVG_DAT',self.te_use_avg_dat,4)
			self.ne_use_avg_dat = self.read_namelist_str(line,'NE_AVG_DAT',self.ne_use_avg_dat,4)
			self.ti_use_avg_dat = self.read_namelist_str(line,'TI_AVG_DAT',self.ti_use_avg_dat,4)
			self.vt_use_avg_dat = self.read_namelist_str(line,'VT_AVG_DAT',self.vt_use_avg_dat,4)

			self.te_exclude = self.read_namelist_str(line,'TE_EXC',self.te_exclude,3)
			self.ne_exclude = self.read_namelist_str(line,'NE_EXC',self.ne_exclude,3)
			self.ti_exclude = self.read_namelist_str(line,'TI_EXC',self.ti_exclude,3)
			self.vt_exclude = self.read_namelist_str(line,'VT_EXC',self.vt_exclude,3)
		
			self.te_raw_fit = self.read_namelist_str(line,'TE_RAW_FIT',self.te_raw_fit,4)
			self.ne_raw_fit = self.read_namelist_str(line,'NE_RAW_FIT',self.ne_raw_fit,4)
			self.ti_raw_fit = self.read_namelist_str(line,'TI_RAW_FIT',self.ti_raw_fit,4)
			self.vt_raw_fit = self.read_namelist_str(line,'VT_RAW_FIT',self.vt_raw_fit,4)
			
			self.te_use_oli = self.read_namelist_str(line,'TE_USE_OUT',self.te_use_oli,4)
			self.ne_use_oli = self.read_namelist_str(line,'NE_USE_OUT',self.ne_use_oli,4)
			self.ti_use_oli = self.read_namelist_str(line,'TI_USE_OUT',self.ti_use_oli,4)
			self.vt_use_oli = self.read_namelist_str(line,'VT_USE_OUT',self.vt_use_oli,4)

			self.te_oli_cut = self.read_namelist_str(line,'TE_OLI_CUT',self.te_oli_cut,2)
			self.ne_oli_cut = self.read_namelist_str(line,'NE_OLI_CUT',self.ne_oli_cut,2)
			self.ti_oli_cut = self.read_namelist_str(line,'TI_OLI_CUT',self.ti_oli_cut,2)
			self.vt_oli_cut = self.read_namelist_str(line,'VT_OLI_CUT',self.vt_oli_cut,2)
			
			self.te_oli_cutn = self.read_namelist_str(line,'TE_OLI_CUTN',self.te_oli_cutn,1)
			self.ne_oli_cutn = self.read_namelist_str(line,'NE_OLI_CUTN',self.ne_oli_cutn,1)
			self.ti_oli_cutn = self.read_namelist_str(line,'TI_OLI_CUTN',self.ti_oli_cutn,1)
			self.vt_oli_cutn = self.read_namelist_str(line,'VT_OLI_CUTN',self.vt_oli_cutn,1)

			self.fit_type_te = self.read_namelist_str(line,'FIT_TYPE_TE',self.fit_type_te,1)
			self.fit_type_ne = self.read_namelist_str(line,'FIT_TYPE_NE',self.fit_type_ne,1)
			self.fit_type_ti = self.read_namelist_str(line,'FIT_TYPE_TI',self.fit_type_ti,1)
			self.fit_type_vt = self.read_namelist_str(line,'FIT_TYPE_VT',self.fit_type_vt,1)

			self.s_ne = self.read_namelist_str(line,'S_NE',self.s_ne,2)
			self.s_te = self.read_namelist_str(line,'S_TE',self.s_te,2)
			self.s_ti = self.read_namelist_str(line,'S_TI',self.s_ti,2)
			self.s_vt = self.read_namelist_str(line,'S_VT',self.s_vt,2)
	
			self.std_te = self.read_namelist_str(line,'STD_TE',self.std_te,4)
			self.std_ne = self.read_namelist_str(line,'STD_NE',self.std_ne,4)
			self.std_ti = self.read_namelist_str(line,'STD_TI',self.std_ti,4)
			self.std_vt = self.read_namelist_str(line,'STD_VT',self.std_vt,4)
			
			self.use_raws_te = self.read_namelist_str(line,'RAW_STD_TE',self.use_raws_te,4)
			self.use_raws_ne = self.read_namelist_str(line,'RAW_STD_NE',self.use_raws_ne,4)
			self.use_raws_ti = self.read_namelist_str(line,'RAW_STD_TI',self.use_raws_ti,4)
			self.use_raws_vt = self.read_namelist_str(line,'RAW_STD_VT',self.use_raws_vt,4)
		
			self.amain = self.read_namelist_str(line,'AMAIN',self.amain,2)
			self.zeff = self.read_namelist_str(line,'ZEFF',self.zeff,2)
			self.aimp = self.read_namelist_str(line,'AIMP',self.aimp,2)
			self.zimp = self.read_namelist_str(line,'ZIMP',self.zimp,2)

			self.eqdsk_name = self.read_namelist_str(line,'EQDSK',self.eqdsk_name,3)
			self.use_chease = self.read_namelist_str(line,'USE_CHEASE',self.use_chease,4)
			self.use_rho = self.read_namelist_str(line,'USE_RHO',self.use_rho,4)

		f4.close()

		self.input_file_check()
		self.adjust_variables()
		
		return

	def input_file_check(self):

		if not (self.te_file == None):
			if not (os.path.isfile(self.te_file)): self.te_file = None
		if not (self.ne_file == None):
			if not (os.path.isfile(self.ne_file)): self.ne_file = None

		if not (self.ti_file == None):	
			if not (os.path.isfile(self.ti_file)): self.ti_file == None
		if not (self.te_dat_file == None):
			if not (os.path.isfile(self.te_dat_file)): self.te_dat_file = None
		if not (self.ne_dat_file == None):
			if not (os.path.isfile(self.ne_dat_file)): self.ne_dat_file = None
		if not (self.ti_dat_file == None):
			if not (os.path.isfile(self.ti_dat_file)): self.ti_dat_file = None

		return
		
	def adjust_variables(self):

		if (self.target_density > 0.0):
				self.use_density_scale = True

		self.fit_type_te = 4
		self.fit_type_ne = 4
		self.fit_type_ti = 4

		try:
			self.ne_excn = self.ne_exclude.split(',')
			self.ne_exclude = ''
		except:
			self.ne_excn = []
		try:
			self.te_excn = self.te_exclude.split(',')
			self.te_exclude = ''
		except:
			self.te_excn = []
		try:
			self.ti_excn = self.ti_exclude.split(',')
			self.ti_exclude = ''
		except:
			self.ti_excn = []
		try:
			self.vt_excn = self.ti_exclude.split(',')
			self.vt_exclude = ''
		except:
			self.vt_excn = []			

		self.no_te = False
		self.no_ne = False
		self.no_ti = False
			
		if (self.te_dat_file == None and self.te_file == None): self.no_te = True #exit('No te data')		
		if (self.ne_dat_file == None and self.ne_file == None): self.no_ne = True #exit('No ne data')
		if (self.ti_dat_file == None and self.ti_file == None): self.no_ti = True #exit('No ti data')

		if (self.te_file == None): self.te_raw_fit = True
		if (self.ne_file == None): self.ne_raw_fit = True
		if (self.ti_file == None): self.ti_raw_fit = True
		
		if(self.te_dat_file == None): self.te_raw_fit = False
		if(self.ne_dat_file == None): self.ne_raw_fit = False
		if(self.ti_dat_file == None): self.ti_raw_fit = False

		if (self.use_ti_width):
			self.fixed_width_te = True
			self.fixed_width_ne = True
			
		return	
			
	def initialise_variables(self):
	
		self.use_ti_width = False
	
		self.psi_end_te = 0.00
		self.psi_end_ti = 0.00
		self.psi_end_ne = 0.00
		self.psi_end_vt = 0.00

		self.te_file = None
		self.ti_file = None
		self.ne_file = None
		self.vt_file = None
		
		self.te_dat_file = None
		self.ti_dat_file = None
		self.ne_dat_file = None
		self.vt_dat_file = None
		
		self.target_density = 0.0
		self.use_density_scale = False
		self.scaled_density = 1.0
		self.use_ti_eped = True

		self.fixed_sep_te = False
		self.fixed_sep_ne = False
		self.fixed_sep_ti = False
		self.fixed_sep_vt = False
		
		self.ne_sep = 0.0
		self.te_sep = 0.0
		self.ti_sep = 0.0
		self.vt_sep = 0.0

		self.te_height = 0.0
		self.ne_height = 0.0
		self.ti_height = 0.0
		self.vt_height = 0.0

		self.te_core = 0.0
		self.ne_core = 0.0
		self.ti_core = 0.0
		self.vt_core = 0.0		
		
		self.ped_scan_fit  = True

		self.fixed_width_te = False
		self.fixed_width_ne = False
		self.fixed_width_ti = False
		self.fixed_width_vt = False

		self.te_width = 0.0
		self.ti_width = 0.0
		self.ne_width = 0.0
		self.vt_width = 0.0

		self.te_pedmid = 0.97
		self.ti_pedmid = 0.97
		self.ne_pedmid = 0.97

		self.te_use_avg_dat = True
		self.ti_use_avg_dat = True
		self.ne_use_avg_dat = True
		self.vt_use_avg_dat = True

		self.te_exclude = ''
		self.ne_exclude = ''
		self.ti_exclude = ''
		self.vt_exclude = ''

		self.te_raw_fit = False
		self.ne_raw_fit = False
		self.ti_raw_fit = False
		self.vt_raw_fit = False
		
		self.te_use_oli = False
		self.ne_use_oli = False
		self.ti_use_oli = False
		self.vt_use_oli = False

		self.te_oli_cut = 0.0
		self.ne_oli_cut = 0.0
		self.ti_oli_cut = 0.0
		self.vt_oli_cut = 0.0

		self.te_oli_cutn = 4
		self.ne_oli_cutn = 4
		self.ti_oli_cutn = 4
		self.vt_oli_cutn = 4

		self.fit_type_te = 4
		self.fit_type_ne = 4
		self.fit_type_ti = 4
		self.fit_type_vt = 4

		self.s_te = -1.
		self.s_ne = -1.
		self.s_ti = -1.
		self.s_vt = -1.

		self.std_te = False
		self.std_ne = False
		self.std_ti = False
		self.std_vt = False

		self.amain = 2.0
		self.zeff = 2.0
		self.zimp = 6.0
		self.aimp = 12.0

		self.te_datz = 0.
		self.ne_datz = 0.
		self.ti_datz = 0.
		self.vt_datz = 0.

		self.psin1 = np.linspace(0,1.0,401)
		self.psin2 = np.linspace(0,1.2,401)
		self.oli_in = np.ones(shape=(3,2))

		self.te_fit2 = np.ones(len(self.psin2))
		self.ne_fit2 = np.ones(len(self.psin2))
		self.ti_fit2 = np.ones(len(self.psin2))
		self.vt_fit2 = np.ones(len(self.psin2))

		self.use_chease = False
		self.eqdsk_name = 'geqdsk'
		self.use_rho = False

		self.popt_te = np.ones(9)
		self.popt_ne = np.ones(9)
		self.popt_ti = np.ones(9)
		self.popt_vt = np.ones(9)

		self.fit_options = [[0]*15 for i in range(5)]

		self.same_tefit = False
		self.same_nefit = False
		self.same_tifit = False
		self.same_vtfit = False	

		self.use_raws_te = False
		self.use_raws_ne = False
		self.use_raws_ti = False
		self.use_raws_vt = False	

		return
			
	def den_scale(self,type=1,konly=False):	
		self.read_kinprof_kprofile()
		self.read_raw_kinprof()

		if (self.target_density <= 0.0):
			self.use_density_scale = False
			self.scale_density(type,True)
			self.scaled_density = 1.0

		if (self.use_density_scale):

			self.scale_density(type)
			try:	self.nek = self.nek * self.scaled_density
			except:	pass
			if not (konly):
				try:	self.ne_fit2 = self.ne_fit2 * self.scaled_density
				except:	pass
			
		return
		
	def write_avg_dat(self,datx,daty,flag):
	
		f = open('PROFILES/'+flag.upper()+'_avg.dat','w')

		x_type = 'Psi_Norm'
		x_type2 = 'Rho_Norm'
		datx2 = self.psi_to_rho(datx)
	
		if (flag == 'ne'):
			f.write(x_type+'  '+x_type2+', NE[1E18/m3], STD \n')
		elif (flag == 'te'):
			f.write(x_type+'  '+x_type2+', TE[keV], 	STD \n')
		elif (flag == 'ti'):
			f.write(x_type+'  '+x_type2+', TI[kev],	    STD \n')
		
		for i in range(len(datx)):
			f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(datx[i],datx2[i],daty[i,0],daty[i,1]))
		f.close()

		return
				
	def main_run(self):

		try:	os.mkdir('PROFILES')
		except:	pass

		if not self.eqdsk_name_temp == self.eqdsk_name:
			self.eqdsk_name_temp = self.eqdsk_name
			self.mapping_variables()
			print('>>> Re-equilibrium mapping!')

		self.read_kinprof_kprofile()
		self.read_raw_kinprof()
		if not (self.fit_type_te_temp == self.fit_type_te): 
			self.fit_type_te_temp = self.fit_type_te
			self.lmfit_init_param(self.fit_type_te,'te')
		if not (self.fit_type_ne_temp == self.fit_type_ne): 
			self.fit_type_ne_temp = self.fit_type_ne
			self.lmfit_init_param(self.fit_type_ne,'ne')

		if not (self.fit_type_ti_temp == self.fit_type_ti): 
			self.fit_type_ti_temp = self.fit_type_ti
			self.lmfit_init_param(self.fit_type_ti,'ti')
		
		if (self.te_use_avg_dat and not self.no_te):	self.write_avg_dat(self.te_datx2,self.te_daty2,'te')
		if (self.ne_use_avg_dat and not self.no_ne):	self.write_avg_dat(self.ne_datx2,self.ne_daty2,'ne')
		if (self.ti_use_avg_dat and not self.no_ti):	self.write_avg_dat(self.ti_datx2,self.ti_daty2,'ti')

		if (self.use_raws_te and np.sum(self.te_dats) == 0.):	
			self.use_raws_te = False
			print('>>> No Te Raw STD')
		if (self.use_raws_ne and np.sum(self.ne_dats) == 0.):        
			self.use_raws_ne = False
			print('>>> No Ne Raw STD')
		if (self.use_raws_ti and np.sum(self.ti_dats) == 0.):        
			self.use_raws_ti = False
			print('>>> No Ti Raw STD')
	
		print('=----------------------------------------------------------------------------=')
		print('=                             FITTING STARTED                                =')
		print('=----------------------------------------------------------------------------=')
		
		self.check_options_change()

		if not (self.no_ti or self.same_tifit):

			if (self.ti_raw_fit):
				xx = self.ti_datx
				yy = self.ti_daty
				ss = self.ti_dats
			else:
				xx = self.psink
				yy = np.copy(self.tik)
				ss = np.ones(len(self.psink))

			self.ti_fit2 = self.lmfit_fit(xx,yy,ss,self.ti_excn,self.fit_type_ti,self.fixed_width_ti,self.ti_width,1.e4,self.fixed_sep_ti,self.ti_sep,self.psi_end_ti,self.std_ti,self.s_ti,self.ti_use_avg_dat,self.ti_use_oli,self.ti_oli_cut,self.ti_oli_cutn,self.ti_raw_fit,self.use_raws_ti,'ti')

			self.popt_ti = np.copy(self.popt)
			self.oli_in_ti = np.copy(self.oli_in)
			self.datx2_ti = np.copy(self.datx2)
			self.daty2_ti = np.copy(self.daty2)
			self.dats2_ti = np.copy(self.dats2)
			self.ti_xi = self.chi

			if (self.use_ti_width):
				self.te_width = self.popt[2]
				self.ne_width = self.popt[2]
				if (self.fit_type_ti == 6):
					self.te_pedmid = self.popt[6]
					self.ne_pedmid = self.popt[6]

		if not (self.no_ne or self.same_nefit):

			if (self.ne_raw_fit):
				xx = self.ne_datx
				yy = self.ne_daty
				ss = self.ne_dats
			else:
				xx = self.psink
				yy = np.copy(self.nek)
				ss = np.ones(len(self.nek))

			self.ne_fit2 = self.lmfit_fit(xx,yy,ss,self.ne_excn,self.fit_type_ne,self.fixed_width_ne,self.ne_width,self.ne_pedmid,self.fixed_sep_ne,self.ne_sep,self.psi_end_ne,self.std_ne,self.s_ne,self.ne_use_avg_dat,self.ne_use_oli,self.ne_oli_cut,self.ne_oli_cutn,self.ne_raw_fit,self.use_raws_ne,'ne')

			self.popt_ne = np.copy(self.popt)
			self.oli_in_ne = np.copy(self.oli_in)
			self.datx2_ne = np.copy(self.datx2)
			self.daty2_ne = np.copy(self.daty2)
			self.dats2_ne = np.copy(self.dats2)
			self.ne_xi = self.chi

		if not (self.no_te or self.same_tefit):

			if (self.te_raw_fit):
				xx = self.te_datx
				yy = self.te_daty
				ss = self.te_dats
			else:
				xx = self.psink
				yy = np.copy(self.tek)
				ss = np.ones(len(self.tek))

			self.te_fit2 = self.lmfit_fit(xx,yy,ss,self.te_excn,self.fit_type_te,self.fixed_width_te,self.te_width,self.te_pedmid,self.fixed_sep_te,self.te_sep,self.psi_end_te,self.std_te,self.s_te,self.te_use_avg_dat,self.te_use_oli,self.te_oli_cut,self.te_oli_cutn,self.te_raw_fit,self.use_raws_te,'te')

			self.popt_te = np.copy(self.popt)
			self.oli_in_te = np.copy(self.oli_in)
			self.datx2_te = np.copy(self.datx2)
			self.daty2_te = np.copy(self.daty2)
			self.dats2_te = np.copy(self.dats2)
			self.te_xi = self.chi

		print('=----------------------------------------------------------------------------=')
		print('=                             FITTING IS DONE                                =')
		print('=----------------------------------------------------------------------------=')


		return

	def make_psiRZ_extended(self,psi,Z):

		ind = np.where(self.R > self.rmag)
		Rt = self.R[ind]
		R = np.linspace(self.rmag,max(Rt),301)
		psi_map = self.psif(R,Z)
		psi_map[0] = 0.0

		psif = interp1d(psi_map,R,'cubic')

		return(psif(psi))

	def mapping_variables(self):

		if not (self.use_chease):

			self.eq = eqdsk.eqdsk(self.eqdsk_name,False)
			self.eq.read_eqdsk_file()

			self.rho_map = np.copy(self.eq.prhoR[:,1])
			self.psi_map = np.copy(self.eq.prhoR[:,0])
			self.R = np.copy(self.eq.R)
			self.Z = np.copy(self.eq.Z)
			self.rmag = self.eq.rmag

			self.psif = interp2d(self.R,self.Z,(self.eq.psirz-self.eq.smag)/(self.eq.sbdy-self.eq.smag))

			psin = self.psif(self.R,0.0)
			rmin = min(self.eq.rzbdy[:,0])
			rmax = max(self.eq.rzbdy[:,0])

			ind = np.where(self.R >= (rmin-0.1))
			R = self.R[ind]
			psin = psin[ind]

			ind = np.where(R <= (rmax+0.1))
			R = R[ind]
			psin = psin[ind]

			ind = np.where(R<=(self.rmag-0.02))
			rpf = interp1d(psin[ind],R[ind],'cubic')
			self.Rin = rpf(1.0)
			ind = np.where(R>=(self.rmag+0.02))
			rpf = interp1d(psin[ind],R[ind],'cubic')
			self.Rout = rpf(1.0)		
			self.dat_numh = 201
			self.R_den = np.linspace(self.Rin,self.Rout,self.dat_numh)

			prf = interp1d(R,psin,'cubic')
			self.psi_den = prf(self.R_den)
			self.psi_den[0] = 1.0
			self.psi_den[-1] = 1.0

		else:
			ch = ch_tool.chease(self.chease_name)
			self.ch = ch
			self.ch.eqdsk_name = self.eqdsk_name
			ch.load_eqdsk()
			ch.read_hagermap(ch.hager_mapf)
			ch.read_rho_psi_R(ch.chease_rundir+'/RHO_PSI_R')

			self.R = np.copy(self.ch.eq.R)
			self.Z = np.copy(self.ch.eq.Z)
			self.rmag = self.ch.eq.rmag

			self.dat_numh = self.ch.dat_numh
			self.rho_map = np.copy(self.ch.rho_map)
			self.psi_map = np.copy(self.ch.psi_map)
		
			self.psif = interp2d(ch.eq.R,ch.eq.Z,(ch.eq.psirz-ch.eq.smag)/(ch.eq.sbdy-ch.eq.smag))


			self.R_den = np.zeros(2*self.ch.dat_numh-1)
			self.psi_den = np.copy(self.R_den)		
			
			for i in range(self.ch.dat_numh-1):
				self.R_den[i] = self.ch.rch[self.ch.dat_numh-1-i] - self.ch.rh[self.ch.dat_numh-1-i]
				self.R_den[2*self.ch.dat_numh-2-i] = self.ch.rch[self.ch.dat_numh-1-i] + self.ch.rh[self.ch.dat_numh-1-i]
				self.psi_den[i] = self.ch.psinh[self.ch.dat_numh-1-i]
				self.psi_den[2*self.ch.dat_numh-2-i] =  self.ch.psinh[self.ch.dat_numh-1-i]

			self.R_den[self.ch.dat_numh-1] = self.ch.rch[0]
		
		psi_map_extend = np.zeros(len(self.psi_map)+100)
		rho_map_extend = np.copy(psi_map_extend)
	
		self.psi_to_rho = interp1d(self.psi_map,self.rho_map,'cubic')

		for i in range(len(self.psi_map)):
			psi_map_extend[i] = self.psi_map[i]
			rho_map_extend[i] = self.rho_map[i]
		drdp = (1.0-self.psi_to_rho(0.9999))/(1.0-0.9999)

		for i in range(100):
			psi_map_extend[len(self.psi_map)+i] = 1.0 + 0.2/30.*float(i+1)
			rho_map_extend[len(self.psi_map)+i] = 1.0 + drdp * 0.2/30.*float(i+1)
		
		self.rho_to_psi = interp1d(rho_map_extend,psi_map_extend,'cubic')
		self.psi_to_rho = interp1d(psi_map_extend,rho_map_extend,'cubic')

		return

	def put_options(self):

		option = [[0]*15 for i in range(5)]

		option[1][0] = self.ti_raw_fit
		option[1][1] = self.ti_excn
		option[1][2] = self.fit_type_ti
		option[1][3] = self.fixed_width_ti
		option[1][4] = self.ti_width
		option[1][5] = self.fixed_sep_ti
		option[1][6] = self.ti_sep
		option[1][7] = self.psi_end_ti
		option[1][8] = self.std_ti
		option[1][9] = self.s_ti
		option[1][10] = self.ti_use_avg_dat
		option[1][11] = self.ti_use_oli
		option[1][12] = self.ti_oli_cut
		option[1][13] = self.ti_oli_cutn
		option[1][14] = self.use_raws_ti

		option[2][0] = self.vt_raw_fit
		option[2][1] = self.vt_excn
		option[2][2] = self.fit_type_vt
		option[2][3] = self.fixed_width_vt
		option[2][4] = self.vt_width,
		option[2][5] = self.fixed_sep_vt
		option[2][6] = self.vt_sep
		option[2][7] = self.psi_end_vt
		option[2][8] = self.std_vt
		option[2][9] = self.s_vt
		option[2][10] = self.vt_use_avg_dat
		option[2][11] = self.vt_use_oli
		option[2][12] = self.vt_oli_cut
		option[2][13] = self.vt_oli_cutn
		option[2][14] = self.use_raws_vt

		option[3][0] = self.te_raw_fit
		option[3][1] = self.te_excn
		option[3][2] = self.fit_type_te
		option[3][3] = self.fixed_width_te
		option[3][4] = self.te_width,
		option[3][5] = self.fixed_sep_te
		option[3][6] = self.te_sep
		option[3][7] = self.psi_end_te
		option[3][8] = self.std_te
		option[3][9] = self.s_te
		option[3][10] = self.te_use_avg_dat
		option[3][11] = self.te_use_oli
		option[3][12] = self.te_oli_cut
		option[3][13] = self.te_oli_cutn
		option[3][14] = self.use_raws_te

		option[4][0] = self.ne_raw_fit
		option[4][1] = self.ne_excn
		option[4][2] = self.fit_type_ne
		option[4][3] = self.fixed_width_ne
		option[4][4] = self.ne_width,
		option[4][5] = self.fixed_sep_ne
		option[4][6] = self.ne_sep
		option[4][7] = self.psi_end_ne
		option[4][8] = self.std_ne
		option[4][9] = self.s_ne
		option[4][10] = self.ne_use_avg_dat
		option[4][11] = self.ne_use_oli
		option[4][12] = self.ne_oli_cut
		option[4][13] = self.ne_oli_cutn
		option[4][14] = self.use_raws_ne

		option[0][0] = self.use_rho
		option[0][1] = self.target_density
		option[0][2] = self.zeff
		option[0][3] = self.zimp
		option[0][4] = self.amain
		option[0][5] = self.aimp

		for i in range(15):
			for j in range(5):
				self.fit_options[j][i] = option[j][i]

		return

	def check_options_change(self):

		option = [[0]*15 for i in range(5)]

		for i in range(15):
			for j in range(5):
				option[j][i] = self.fit_options[j][i]
	
		self.put_options()

		self.same_tefit = True
		self.same_nefit = True
		self.same_tifit = True
		self.same_vtfit = True

		if not (self.fit_options[0][0] == option[0][0]):
			self.same_tefit = False
			self.same_nefit = False
			self.same_tifit = False
			self.same_vtfit = False
		
		for i in range(10):
			if not(self.fit_options[0][i+1] == option[0][i+1]):
				self.same_nefit = False

		for i in range(15):
			if not (option[1][i] == self.fit_options[1][i]):
				self.same_tifit = False
				break
		for i in range(15):
			if not (option[2][i] == self.fit_options[2][i]):
				self.same_vtfit = False
				break
		for i in range(15):
			if not (option[3][i] == self.fit_options[3][i]):
				self.same_tefit = False
				break
		for i in range(15):
			if not (option[4][i] == self.fit_options[4][i]):
				self.same_nefit = False
				break
		
		if self.param_change:	
			if not self.compare_dict(self.param2['ne'],self.param2_old['ne']):	self.same_nefit = False
			if not self.compare_dict(self.param2['ti'],self.param2_old['ti']):	self.same_tifit = False
			if not self.compare_dict(self.param2['te'],self.param2_old['te']):	self.same_tefit = False
			if not self.compare_dict(self.param2['vt'],self.param2_old['vt']):	self.same_vtfit = False
			self.param_change = False
		return

	def compare_dict(self,a,b):

		ans = True

		for i in ['vary','val','min','max']:
			for j in range(len(a[i])):
				if not np.isnan(a[i][j]):
					if not a[i][j] == b[i][j]:
						ans = False

		return ans

	def __init__(self,filename):

		self.initialise_variables()
		try:
			self.read_namelist('fit_opt')
		except:
			print('>>> No fit_opt...')

		chk = fit_checkfile.chk_file()
		print('>>> Check Raw data files...')

		ischk1 = False
		ischk2 = False
		ischk3 = False

		if not (self.te_dat_file == None):
			ischk1 = chk.writefiles(self.te_dat_file,'Te')
		if not (self.ne_dat_file == None):
			ischk2 = chk.writefiles(self.ne_dat_file,'Ne')
		if not (self.ti_dat_file == None):
			ischk3 = chk.writefiles(self.ti_dat_file,'Ti')
		if (ischk1 or ischk2 or ischk3):
			print('>>> Source files are saved to xxx_source...')
	
		self.chease_name = filename

		self.eqdsk_name_temp = self.eqdsk_name

		self.mapping_variables()

		self.lmfit_init_params()
	
		fit_types = [self.fit_type_te,self.fit_type_ne,self.fit_type_ti,self.fit_type_vt]
		self.lmfit_load_param(self.param,fit_types)
		self.param2 = copy.deepcopy(self.param)
		self.param2_old = copy.deepcopy(self.param)
		self.param_change = False

		self.fit_type_te_temp = self.fit_type_te
		self.fit_type_ne_temp = self.fit_type_ne
		self.fit_type_ti_temp = self.fit_type_ti
		self.fit_type_vt_temp = self.fit_type_vt


		return		

	if __name__ == "__main__":

		import fittool
		fit = fittool.fit_tool(sys.argv[1])
