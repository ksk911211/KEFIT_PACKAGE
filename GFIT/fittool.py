#!/usr/local/anaconda3/bin/python3
import os, sys
import numpy as np
import time
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import UnivariateSpline as smoothspline
from scipy.interpolate import LSQUnivariateSpline as lspline
from scipy.interpolate import make_interp_spline as mspline
from progress import *
from exec_dirs import dummy_dir
import matplotlib.pyplot as plt
import fit_checkfile
from csaps import csaps
import pickle
import lmfit
import eqdsk
import copy

class fit_tool:

	def print(self,lines):
		if not self.nolog: os.system("echo '%s' >> log.gfit"%lines)
		os.system("echo '%s'"%lines)
		return

	def pretty_print(self,lm,err,colwidth=8, precision=4, fmt='g',columns=['value', 'min', 'max', 'stderr', 'vary', 'expr','brute_step']):

	    name_len = max(len(s) for s in lm)
	    allcols = ['name'] + columns
	    title = '{:{name_len}} ' + len(columns) * ' {:>{n}}'
	    self.print(title.format(*allcols, name_len=name_len, n=colwidth).title())
	    numstyle = '{%s:>{n}.{p}{f}}'  # format for numeric columns
	    otherstyles = dict(name='{name:<{name_len}} ', stderr='{stderr!s:>{n}}',
	                       vary='{vary!s:>{n}}', expr='{expr!s:>{n}}',
	                       brute_step='{brute_step!s:>{n}}')
	    line = ' '.join([otherstyles.get(k, numstyle % k) for k in allcols])

	    count = 0
	    for name, values in sorted(lm.items()):
	        pvalues = {k: getattr(values, k) for k in columns}
	        pvalues['name'] = name
	        # stderr is a special case: it is either numeric or None (i.e. str)
	        if 'stderr' in columns and pvalues['stderr'] is not None:
	            pvalues['stderr'] = (numstyle % '').format(
	                pvalues['stderr'], n=colwidth, p=precision, f=fmt)
	        elif 'brute_step' in columns and pvalues['brute_step'] is not None:
	            pvalues['brute_step'] = (numstyle % '').format(
	                pvalues['brute_step'], n=colwidth, p=precision, f=fmt)
	        self.print(line.format(name_len=name_len, n=colwidth, p=precision,
	                          f=fmt, **pvalues))
	        err[count] = pvalues['stderr']
	        count = count + 1
	       
	    return

	def read_kinprof_kprofile_fun(self,filename,use_rho):
	
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
		
		if not(use_rho):
			datf = interp1d(psink,datk,'cubic')
		else:
			datf = interp1d(rhok,datk,'cubic')
		
		dat = datf(psind)
		
		return (dat)
		
	def read_kinprof_kprofile(self):

		for i in self.prof_list:
			filename = self.fit_opt['file'][i]['kfile']
			use_rho  = self.fit_opt['use_rho'][i]
			if not (filename == None):
				self.__dict__['%s_prof'%i]['kfile']['xxr']  = np.copy(self.fit_eq['psin1'])
				self.__dict__['%s_prof'%i]['kfile']['raw']  = self.read_kinprof_kprofile_fun(filename,use_rho)
				self.__dict__['%s_prof'%i]['kfile']['raw']  = self.__dict__['%s_prof'%i]['kfile']['raw'] * self.factor[i]
				self.__dict__['%s_prof'%i]['kfile']['raws'] = np.copy(self.__dict__['%s_prof'%i]['kfile']['raw'])*1.e-1
				self.__dict__['%s_prof'%i]['kfile']['xxa']  = np.copy(self.fit_eq['psin1'])
				self.__dict__['%s_prof'%i]['kfile']['avg']  = np.copy(self.__dict__['%s_prof'%i]['kfile']['raw'])*1.e+0
				self.__dict__['%s_prof'%i]['kfile']['avgs'] = np.copy(self.__dict__['%s_prof'%i]['kfile']['raw'])*1.e-1

		return

	def read_raw_kinprof_fun(self,filename,shift,use_rho):
	
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
	
		if (use_rho):
			drdp = (1.0 - self.fit_eq['psi_to_rho'](0.9999))/(1.0 - 0.9999)
	
		f4 = open(filename,'r')
		line = f4.readline()
		i = 0
		while (i < linec):
			line = f4.readline().split()
			if (float(line[0]) > 0.):
				datR[i] = float(line[0])+shift
				datZ[i] = float(line[1])	
				datX[i] = float(line[2])
				datS[i] = float(line[3])
				if datS[i] <=0.: datS[i] = datX[i]*0.1
				
				if (use_rho):
					if (self.fit_eq['psif'](datR[i],datZ[i]) < 1.0):
						datP[i] = self.fit_eq['psi_to_rho'](self.fit_eq['psif'](datR[i],datZ[i]))
					else:
						datP[i] = 1.0 + drdp*(self.fit_eq['psif'](datR[i],datZ[i])-1.)
				else:
					datP[i] = self.fit_eq['psif'](datR[i],datZ[i])
				i = i + 1

		#CHECK BETTER FILTER METHOD
		for i in range(linec):
			if not (datX[1] <= 10.0): #Default for ECE false value 10.	
				if (datX[i] > 3.0*datX[1] and datX[1] < 2.e3):
					datX[i] = 0.0
		f4.close()

		return (datP, datX, datS, datR, datZ)

	def get_ch_number(self,datx):

		len1 = len(datx)
		i=0
		for i in range(len1-1):
			if (datx[0] == datx[i+1]):
				break
		len2 = i+1

		if (len2 == (len1-1)): len2 = len1

		return len1,len2

	def make_avg_val(self,datx,daty,dats,flag):

		len1, len2 = self.get_ch_number(datx)
		len3 = int(len1/len2)

		datyt = np.reshape(daty,(len3,len2))
		datst = np.reshape(dats,(len3,len2))
		datxt = datx[0:len2]

		var_crit = self.device['var_crit'][flag.lower()]
		if ((flag.lower() == 'ne') and (self.fit_opt['target_density'] > 0.)):     var_crit = self.fit_opt['target_density'] * 2.

		for i in range(len2):
			maxval = max(datyt[:,i])
			for j in range(len3):
				if datyt[j,i] == 0.:
					if maxval>0.: datyt[j,i] = maxval
					else:         datyt[j,i] = var_crit

		daty2 = np.zeros(len2)
		dats2 = np.zeros(len2)
		dats3 = np.zeros(len2)
		dats4 = np.zeros(len1)

		for i in range(len2):
			sum2 = []
			sums = []
			for j in range(len3):
				if not (datxt[i] > 1.0 and datyt[j,i] > var_crit):
					sum2 = np.append(sum2,datyt[j,i])
					sums = np.append(sums,datst[j,i])
			if (len(sum2) > 0.):
				mean = np.mean(sum2)
				if (np.isnan(mean)):
					mean = 0.e0

				std = np.std(sum2)
				if (np.isnan(std)):
					std = 0.1*mean
				if std <= 0.0: std = 0.1*mean
				if std == 0.0: std = 0.1;
		
				sums = sums ** 2
				snew = np.sum(sums)
				snew = np.sqrt(snew)/len(sum2) #std_err = std_eff / sqrt(n); 
				if snew <= 0.: snew = mean*0.1 #std_eff = sum(std**2) / sqrt(n);
				if snew == 0.0: snew = 0.1     #averaged std of independent groups

				daty2[i] = mean
				dats2[i] = snew
				dats3[i] = std
				for j in range(len3): dats4[len2*j+i]= std
			else:
				daty2[i] = var_crit
				dats2[i] = var_crit * 0.1
				dats3[i] = var_crit * 0.1
				for j in range(len3): dats4[len2*j+i] = var_crit*0.1

		return(datxt,daty2,dats2,dats3,dats4,len2,len3)		
		
	def read_raw_kinprof(self):

		for i in self.prof_list:
			use_rho  = self.fit_opt['use_rho'][i]
			for j in self.__dict__['%s_list'%i]:
				filename = self.fit_opt['file'][i][j]
				shift    = self.fit_opt['shift'][i][j]
				if j=='tse': shift = self.fit_opt['shift'][i]['ts'];

				if not filename == None:
					datx,daty,dats,datR,datZ = self.read_raw_kinprof_fun(filename,shift,use_rho)
					self.__dict__['%s_prof'%i][j]['xxr']   = np.copy(datx)
					self.__dict__['%s_prof'%i][j]['raw']   = np.copy(daty * self.factor[i])
					self.__dict__['%s_prof'%i][j]['raws']  = np.copy(dats * self.factor[i])

					self.__dict__['%s_prof'%i][j]['datR']  = datR
					self.__dict__['%s_prof'%i][j]['datZ']  = datZ

		self.te_ts_raw_ori = np.copy(self.te_prof['ts']['raw'])
		self.te_tse_raw_ori = np.copy(self.te_prof['tse']['raw'])
		self.ne_ts_raw_ori = np.copy(self.ne_prof['ts']['raw'])
		self.ne_tse_raw_ori = np.copy(self.ne_prof['tse']['raw'])

		if self.fit_opt['file']['te']['tse'] == None: self.fit_opt['scale']['te']['ts']['ch_cal'] = False
		if self.fit_opt['file']['ne']['tse'] == None: self.fit_opt['scale']['ne']['ts']['ch_cal'] = False

		for i in self.prof_list:
			use_rho  = self.fit_opt['use_rho'][i]
			
			for j in self.__dict__['%s_list'%i]:
				filename = self.fit_opt['file'][i][j]
				shift    = self.fit_opt['shift'][i][j]
				if j=='tse': shift = self.fit_opt['shift'][i]['ts'];

				if not filename == None:
					datx,daty,dats,dats1,dats2,len2,len3 = self.make_avg_val(self.__dict__['%s_prof'%i][j]['xxr'],self.__dict__['%s_prof'%i][j]['raw'],self.__dict__['%s_prof'%i][j]['raws'],i)
					lenx = len(datx);
					lent = int(len(self.__dict__['%s_prof'%i][j]['xxr'])/lenx)
					for k in range(lenx):
						for m in range(lent):
							factor = self.make_scale_factor(i,j,m)
							for n in ['raw','raws']:
								self.__dict__['%s_prof'%i][j][n][lenx*m+k] =  self.__dict__['%s_prof'%i][j][n][lenx*m+k] * factor

					datx,daty,dats,dats1,dats2,len2,len3   = self.make_avg_val(self.__dict__['%s_prof'%i][j]['xxr'],self.__dict__['%s_prof'%i][j]['raw'],self.__dict__['%s_prof'%i][j]['raws'],i)
					self.__dict__['%s_prof'%i][j]['xxa']   = np.copy(datx)
					self.__dict__['%s_prof'%i][j]['avg']   = np.copy(daty)
					self.__dict__['%s_prof'%i][j]['avgs']  = np.copy(dats)
					self.__dict__['%s_prof'%i][j]['avgs2'] = np.copy(dats1)
					self.__dict__['%s_prof'%i][j]['raws2'] = np.copy(dats2)
					self.post['prof_dim'][i][j] = [len2,len3]
					
		return

	def make_scale_factor(self,prof,flag,m):

		if flag in ['ces','refl','ece']: return self.fit_opt['scale'][prof][flag]['core']
		if not prof in ['te','ne']: return 1.

		if self.fit_opt['scale'][prof]['ts']['ch_cal']:

			lents1, lents2 = self.get_ch_number(self.__dict__['%s_prof'%prof]['ts']['xxr'])
			lentse1,lentse2= self.get_ch_number(self.__dict__['%s_prof'%prof]['tse']['xxr'])			
			cind = lents2*m  + min(lents2,int(self.fit_opt['scale'][prof]['ts']['core_ch']))-1
			eind = lentse2*m + min(lentse2,int(self.fit_opt['scale'][prof]['ts']['edge_ch']))-1	
			vratio = self.__dict__['%s_tse_raw_ori'%prof][eind]/self.__dict__['%s_ts_raw_ori'%prof][cind]
			if flag == 'ts':
				if self.fit_opt['scale'][prof]['ts']['fix_core']: factor = self.fit_opt['scale'][prof]['ts']['core']
				else:
					factor = vratio * self.fit_opt['scale'][prof]['ts']['edge'] / self.fit_opt['scale'][prof]['ts']['ratio']
			else:
				if self.fit_opt['scale'][prof]['ts']['fix_core']: 
					factor =  self.fit_opt['scale'][prof]['ts']['ratio'] / vratio * self.fit_opt['scale'][prof]['ts']['core']
				else: factor = self.fit_opt['scale'][prof]['ts']['edge']
			#if prof in ['te','ne']: print(prof,flag,m,factor)

		else:
			if flag == 'ts': factor = self.fit_opt['scale'][prof]['ts']['core'];
			else: factor = self.fit_opt['scale'][prof]['ts']['edge'];


		return factor

	def eped_fun(self,x,a1,a2,a3,a4,a5,a6):

		y = a1
		y = y + a2*(np.tanh(1) - np.tanh((x - 1.0 + 0.5*(abs(abs(a3)+1.e-7)))/(abs(abs(a3)+1.e-7))*2.0))
		yt = (x/(1.0-abs(abs(a3)+1.e-7))) ** (abs(a5)+1.01)
		y = y + a4 * (abs((1-yt)) **(abs(a6)+1.01)) * 0.5 * (1.0 + np.sign(1.0-abs(abs(a3)+1.e-7)-x))
		
		return y

	def modify_rho_to_psi_eped(self,coeff,flag='ne',flagout=False):

		x = np.linspace(0,1.0,201)
		y = self.eped_prof(x,coeff[0],coeff[1],abs(coeff[2]),coeff[3],abs(coeff[4]),abs(coeff[5]))

		xx = self.fit_eq['rho_to_psi'](x)

		pp = [coeff[0],coeff[1],abs(coeff[2]),coeff[3],abs(coeff[4])-0.8,abs(coeff[5])-0.8]

		popt, pcov = curve_fit(self.eped_fun,xx,y,p0=pp,maxfev=300000)

		new_coef = np.copy(popt)
		new_coef[4] = abs(popt[4]) + 1.01
		new_coef[5] = abs(popt[5]) + 1.01
		if flagout:
			line = flag + ' eped coeffs changed from --- to --- for rho -> psin mapping ...'
			self.print(line)
			line = str(coeff)
			self.print(line)
			line = str(new_coef)
			self.print(line)

		return new_coef

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

	def scale_density(self,use_rho=False,print_only=False,noprint=False):

		nef = interp1d(self.fit_eq['psin2'],self.ne_prof['fit2'],'cubic')
		for flag in ['intp1','intp2']:
			for i in range(len(self.device[flag])):
				if self.device[flag][i]>self.fit_eq['psin2'][-1]: self.device[flag][i] = self.fit_eq['psin2'][-1]

		if not (use_rho):
			NN  = nef(self.device['intp1'])
			NN2 = nef(self.device['intp2'])
		else:
			xx  = self.fit_eq['psi_to_rho'](self.device['intp1'])
			NN  = nef(xx)
			xx  = self.fit_eq['psi_to_rho'](self.device['intp2'])
			NN2 = nef(xx)

		line_sum = np.trapz(NN, x=self.device['intr1'])
		self.post['int01'] = 2.0* line_sum / self.device['inter1L']
		line_sum = np.trapz(NN2,x=self.device['intz2'])
		self.post['int02'] = 2.0* line_sum / self.device['inter2L']

		for i in ['tci01','tci02','tci03','tci04','tci05']:
			self.post[i] = self.get_tci_line(self.device['tci_Rmid'][i],self.device['tci_Rend'][i],nef,self.device['tci_L'][i]/4.,use_rho)

		for i in self.inter_list:
			if self.fit_opt[i]['val']> 0.:
				self.post['den_diff'][i] = self.fit_opt[i]['val'] - self.post[i]
			else: self.post['den_diff'][i] = 0.;

		if (print_only or self.fit_opt['line_type'] == 0):
			if not noprint: self.print_scale_density(0.,False)
			return

		if self.fit_opt['line_type'] > 0:	line_avg = self.post[self.inter_list[self.fit_opt['line_type']-1]]
		scaled_density = self.fit_opt['target_density'] / line_avg

		if not (self.post['scaled_density'] == scaled_density):
			if (self.fit_opt['use_density_scale']):
				self.print_scale_density(self.fit_opt['target_density'],True)
				self.post['scaled_density'] = scaled_density

		return

	def den_scale(self,konly=False):	
		self.read_kinprof_kprofile()
		self.read_raw_kinprof()

		if (self.fit_opt['target_density'] <= 0.0 or self.fit_opt['line_type'] == 0):
			self.fit_opt['use_density_scale'] = False
			self.scale_density(self.fit_opt['use_rho']['ne'],True,True)
			self.post['scaled_density'] = 1.0

		if (self.fit_opt['use_density_scale']):

			self.scale_density(self.fit_opt['use_rho']['ne'])
			try:	self.ne_prof['fit_old'] = self.ne_prof['fit_old']#* self.post['scaled_density']
			except:	pass
			try:    self.ne_prof['kfile']['raw']   = self.ne_prof['kfile']['raw'] * self.post['scaled_density']
			except: pass
			if not (konly):
				try:	
						self.ne_prof['fit2'] = self.ne_prof['fit2'] * self.post['scaled_density']
						self.ne_prof['fit1'] = self.ne_prof['fit1'] * self.post['scaled_density']
						self.ne_prof['fit2p']= self.ne_prof['fit2p']* self.post['scaled_density'] 			
				except:	pass
			
		return		
	
	def print_scale_density(self,line_avgt,use_scale):
		if self.noprint: return
		for i in range(len(self.inter_list)):
			ind = self.inter_list[i]
			line_avg  = self.post[ind]
			line_avge = self.fit_opt[ind]['val']
			if (i==(self.fit_opt['line_type']-1) and use_scale):	
				self.print('>>> Line average density(%s)      = %5.3f  >>> %5.3f [10(19)/m3]'%(ind,line_avg,line_avgt))
			else:
				if (line_avge > 0):
					self.print('>>> Line average density(%s)      = %5.3f, EXP %5.3f [10(19)/m3] (Diff %+5.3f)'%(ind,line_avg,line_avge,line_avg-line_avge))
				else:	
					self.print('>>> Line average density(%s)      = %5.3f, EXP -N/A- [10(19)/m3]'%(ind,line_avg))
		if use_scale:	self.print('=----------------------------------------------------------------------------=')

		return

	def get_tci_line(self,Rmid,Rend,nef,LL,use_rho):
		
		L_traj = np.sqrt(Rend**2 - Rmid**2)
		R1 = np.linspace(0.,L_traj,200)
		R2 = np.copy(R1); NN = np.copy(R2)
		for i in range(len(R2)): R2[i] = np.sqrt(Rmid**2+ R1[i]**2)
		psi_den = self.fit_eq['psif'](R2,0.)

		if use_rho: xx = self.fit_eq['psi_to_rho'](psi_den)
		else: xx = np.copy(psi_den)
		
		for i in range(len(xx)):
			if xx[i] >=self.fit_eq['psin2'][-1]: xx[i] = self.fit_eq['psin2'][-1]
			NN[i] = nef(xx[i])
			NN[i] = max(NN[i],0.01)

		line_sum = np.trapz(NN,x=R1)
		line_avg = line_sum /L_traj
		return line_avg
		
	def residual(self,p,x,y,s,ch_list,fun,flag):
		yt = np.array([])
		nch = len(ch_list)
		for i in range(nch):
			yf = fun(x[i],**p)
			yf = np.round(yf,self.fit_trunc)
			ytt = (yf - y[i]) / (abs(s[i])+1.e-8)

			if not ch_list[0] == 'kfile':
				if (self.fit_opt['oli'][flag][ch_list[i]]['use'] and not self.fit_opt['avg'][flag][ch_list[i]]):
					oli_cut = self.fit_opt['oli'][flag][ch_list[i]]['per']
				else:
					oli_cut = 1.0

				if (oli_cut == 0.0):
					oli_cut2 = 95
				else:
					oli_cut2 = oli_cut * 100

				yn = np.percentile(abs(ytt),oli_cut2)>=abs(ytt)
				yt = np.hstack([yt,ytt[yn]*self.fit_opt['weight'][flag][ch_list[i]]/np.sqrt(len(ytt[yn]))])
			else: yt = np.hstack([yt,ytt/np.sqrt(len(yyt))])

		if flag == 'ne':
			use_rho = self.fit_opt['use_rho'][flag]

			tci_n = [self.fit_opt['tci01']['val'], self.fit_opt['tci02']['val'], self.fit_opt['tci03']['val'], self.fit_opt['tci04']['val'], self.fit_opt['tci05']['val']]
			int_n = [self.fit_opt['int01']['val'], self.fit_opt['int02']['val']]
			tci_s = [self.fit_opt['tci01']['sig'], self.fit_opt['tci02']['sig'], self.fit_opt['tci03']['sig'], self.fit_opt['tci04']['sig'], self.fit_opt['tci05']['sig']]
			int_s = [self.fit_opt['int01']['sig'], self.fit_opt['int02']['sig']]

			tci_n = np.round(tci_n,self.line_trunc)
			tci_s = np.round(tci_s,self.line_trunc)
			int_n = np.round(int_n,self.line_trunc)
			int_s = np.round(int_s,self.line_trunc)

			if use_rho:
				xx = self.fit_eq['psi_to_rho'](self.device['tcip1'])
				tci1_n = fun(xx,**p)
				xx = self.fit_eq['psi_to_rho'](self.device['tcip2'])
				tci2_n = fun(xx,**p)
				xx = self.fit_eq['psi_to_rho'](self.device['tcip3'])
				tci3_n = fun(xx,**p)
				xx = self.fit_eq['psi_to_rho'](self.device['tcip4'])
				tci4_n = fun(xx,**p)
				xx = self.fit_eq['psi_to_rho'](self.device['tcip5'])
				tci5_n = fun(xx,**p)								
				xx = self.fit_eq['psi_to_rho'](self.device['intp1'])
				int1_n = fun(xx,**p)
				xx = self.fit_eq['psi_to_rho'](self.device['intp2'])
				int2_n = fun(xx,**p)
			else:
				tci1_n = fun(self.device['tcip1'],**p)
				tci2_n = fun(self.device['tcip2'],**p)
				tci3_n = fun(self.device['tcip3'],**p)
				tci4_n = fun(self.device['tcip4'],**p)
				tci5_n = fun(self.device['tcip5'],**p)				
				int1_n = fun(self.device['intp1'],**p)
				int2_n = fun(self.device['intp2'],**p)

			tci1_a = np.trapz(tci1_n,x=self.device['tcir1'])
			#tci1_a = tci1_a / self.device['tci_L']['tci01'] * 4.
			tci2_a = np.trapz(tci2_n,x=self.device['tcir2'])
			#tci2_a = tci2_a / self.device['tci_L']['tci02'] * 4.
			tci3_a = np.trapz(tci3_n,x=self.device['tcir3'])
			#tci3_a = tci3_a / self.device['tci_L']['tci03'] * 4.
			tci4_a = np.trapz(tci4_n,x=self.device['tcir4'])
			#tci4_a = tci4_a / self.device['tci_L']['tci04'] * 4.
			tci5_a = np.trapz(tci5_n,x=self.device['tcir5'])
			#tci5_a = tci5_a / self.device['tci_L']['tci05'] * 4.

			tci1_a = tci1_a / self.device['tcir1'][-1]
			tci2_a = tci2_a / self.device['tcir2'][-1]
			tci3_a = tci3_a / self.device['tcir3'][-1]
			tci4_a = tci4_a / self.device['tcir4'][-1]
			tci5_a = tci5_a / self.device['tcir5'][-1]

			int1_a = np.trapz(int1_n,x=self.device['intr1'])
			int1_a = int1_a/self.device['inter1L']*2.0
			int2_a = np.trapz(int2_n,x=self.device['intz2'])
			int2_a = int2_a/self.device['inter2L']*2.0

			tci1_a = round(tci1_a,self.line_trunc)
			tci2_a = round(tci2_a,self.line_trunc)
			tci3_a = round(tci3_a,self.line_trunc)
			tci4_a = round(tci4_a,self.line_trunc)
			tci5_a = round(tci5_a,self.line_trunc)			
			int1_a = round(int1_a,self.line_trunc)
			int2_a = round(int2_a,self.line_trunc)
			
			tci1_d = abs(tci1_a-tci_n[0])/tci_s[0]
			tci2_d = abs(tci2_a-tci_n[1])/tci_s[1]
			tci3_d = abs(tci3_a-tci_n[2])/tci_s[2]
			tci4_d = abs(tci4_a-tci_n[3])/tci_s[3]
			tci5_d = abs(tci5_a-tci_n[4])/tci_s[4]

			if self.exp_lsq:
				tci1_d = tci1_d**3
				tci2_d = tci2_d**3
				tci3_d = tci3_d**3
				tci4_d = tci4_d**3
				tci5_d = tci5_d**3

			if (tci_s[0] > 0.and  tci_n[0]>0.): yt = np.append(yt,round(tci1_d,self.line_trunc))
			if (tci_s[1] > 0.and  tci_n[1]>0.): yt = np.append(yt,round(tci2_d,self.line_trunc))
			if (tci_s[2] > 0.and  tci_n[2]>0.): yt = np.append(yt,round(tci3_d,self.line_trunc))
			if (tci_s[3] > 0.and  tci_n[3]>0.): yt = np.append(yt,round(tci4_d,self.line_trunc))
			if (tci_s[4] > 0.and  tci_n[4]>0.): yt = np.append(yt,round(tci5_d,self.line_trunc))			

			if (int_s[0] > 0. and int_n[0]): yt = np.append(yt,round((int1_a-int_n[0])/int_s[0],self.line_trunc))
			if (int_s[1] > 0. and int_n[1]): yt = np.append(yt,round((int2_a-int_n[1])/int_s[1],self.line_trunc))
			
		return (yt)

	def lmfit_init(self,flag='none'):

		if not flag == 'none':
			i = flag.lower()
			self.param[i] = dict()
			self.param[i]['vary'] = [True]*10
			self.param[i]['val'] = np.zeros(10)
			self.param[i]['min'] = np.zeros(10)
			self.param[i]['max'] = np.zeros(10)
			for j in range(10):
				self.lmfit_set_param(self.param[i],j,np.nan,-np.inf,np.inf)

			return

		self.param = dict()
		for i in ['te','ne','ti','vt']:
			self.param[i] = dict()
			self.param[i]['vary'] = [True]*10
			self.param[i]['val'] = np.zeros(10)
			self.param[i]['min'] = np.zeros(10)
			self.param[i]['max'] = np.zeros(10)
			for j in range(10):
				self.lmfit_set_param(self.param[i],j,np.nan,-np.inf,np.inf)
		
		return		

	def lmfit_set_param(self,param,ind,val,minv,maxv,vary=True):

		ind = ind - 1

		param['vary'][ind] = vary

		if not (minv == None):  param['min'][ind] = minv
		if not (maxv == None):	param['max'][ind] = maxv
		if not (val == None):	param['val'][ind] = val		

		nan_val = np.isnan(param['val'][ind])
		inf_min = np.isinf(param['min'][ind])
		inf_max = np.isinf(param['max'][ind])

		if param['min'][ind] == param['max'][ind]:
			if inf_min:
				param['max'][ind] = np.inf
				param['min'][ind] = -np.inf
			else:
				if not ((param['min'][ind]==0.) and  (param['max'][ind] ==0.)):
					param['min'][ind] = 0.9 * param['min'][ind]
					param['max'][ind] = 1.1 * param['max'][ind]
				else:
					param['min'][ind] = +0.
					param['max'][ind] = +0.1

		if nan_val:
			if (inf_min and inf_max):	param['val'][ind] = 0.
			elif (inf_min):			param['val'][ind] = 0.9 * param['max'][ind]
			elif (inf_max):			param['val'][ind] = 1.1 * param['min'][ind]
			else:				param['val'][ind] = 0.5 * (param['min'][ind] + param['max'][ind])

		else:
			if (vary):
				if (param['val'][ind] < param['min'][ind]):	param['val'][ind] = param['min'][ind]
				if (param['val'][ind] > param['max'][ind]):	param['val'][ind] = param['max'][ind]
			else:
				if (param['val'][ind] < param['min'][ind]):   param['min'][ind] = param['val'][ind]
				if (param['val'][ind] > param['max'][ind]):   param['max'][ind] = max(param['val'][ind]*1.1,param['val'][ind]+0.01)

		return

	def lmfit_load_param(self,param,fit_types,filename='fit_opt_param'):

		try:
			f = open(filename,'r')
			self.print('>>> Load lmfit params...')
		except:
			self.print('>>> No saved lmfit params...')
			return
		k = -1
		for i in ['te','ne','ti','vt']:
			line = f.readline()
			line = line.split('=')
			fit_type = int(line[1])
			k = k + 1
			if (fit_type == fit_types[k]):
				for j in range(10):
					paramt = param[i]
					line = f.readline().split()
					var1 = line[0]
					var2 = line[1]
					var3 = line[2]
					var4 = line[3]

					if var1.lower() == 'inf': var1 = np.nan
					if var2.lower() == '-inf': var2 = -np.inf
					if var3.lower() == 'inf': var3 = np.inf
					if var4.lower() == 'true': var4 = True
					elif var4.lower() == 'false': var4 = False
					else: self.print('>>> no logic var')
	
					paramt['val'][j] = float(var1)
					paramt['min'][j] = float(var2)
					paramt['max'][j] = float(var3)
					paramt['vary'][j] = var4
			else:
				self.print('>>> %s fit_opt [type=%i] is not matched with saved option [type=%i]'%(i,fit_types[k],fit_type))
				for j in range(10):
					line = f.readline()

		f.close()

		return		

	def lmfit_write_param(self,param,filename='fit_opt_param.out'):

		fit_type = dict()
		fit_type['te'] = self.fit_opt['func_type']['te']
		fit_type['ne'] = self.fit_opt['func_type']['ne']
		fit_type['ti'] = self.fit_opt['func_type']['ti']
		fit_type['vt'] = self.fit_opt['func_type']['vt']

		f = open(filename,'w')
		for i in ['te','ne','ti','vt']:
			line = i + '_fit_variables(val,min,max,vary), fit_type = %i\n'%fit_type[i]
			f.write(line)
			for j in range(10):
				paramt = param[i]
				var1 = paramt['val'][j]
				var2 = paramt['min'][j]
				var3 = paramt['max'][j]
				var4 = paramt['vary'][j]

				if var1 == np.nan: var1 = np.inf
				if var2 == -np.inf: var2 = -np.inf
				if var3 == np.inf: var3 = np.inf

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

			try:
				if args['min'] == args['max']:
					if np.isinf(args['max']):
						args['min'] = -np.inf
						args['max'] = +np.inf
					else:
						args['min'] = 0.99*args['min']
						args['max'] = 1.01*args['max']
			except: pass

			p.add(var_name,**args)

		return	

	def lmfit_init_param(self,fit_type,flag):

		if (fit_type == 1):
			if (flag == 'te' or flag == 'ti'):
				self.lmfit_set_param(self.param[flag],1,0.15,0.05,0.3,True)
			elif (flag == 'ne'):
				self.lmfit_set_param(self.param[flag],1,0.5,0.3,1.0,True)
			elif (flag == 'vt'):
				self.lmfit_set_param(self.param[flag],1,50,-np.inf,np.inf,True)

			self.lmfit_set_param(self.param[flag],2,2.0,0.0,np.inf,True)	

			self.lmfit_set_param(self.param[flag],3,1.3,0.5,6.0,True)
			self.lmfit_set_param(self.param[flag],4,2.0,0.5,6.0,True)

		if (fit_type == 2):
			if (flag == 'te' or flag == 'ti'):
				self.lmfit_set_param(self.param[flag],1,0.15,0.05,0.3,True)
			elif (flag == 'ne'):
				self.lmfit_set_param(self.param[flag],1,0.5,0.3,1.0,True)
			elif (flag == 'vt'):
				self.lmfit_set_param(self.param[flag],1,50,0.0,np.inf,True)
				

			self.lmfit_set_param(self.param[flag],2,1.0,0.0,np.inf,True)
			if (flag == 'ti' or flag == 'vt'): self.lmfit_set_param(self.param[flag],3,0.05,0.04,0.15,True)
			else: self.lmfit_set_param(self.param[flag],3,0.05,0.035,0.1,True)
			self.lmfit_set_param(self.param[flag],4,0.96,0.85,1.02,True)

			self.lmfit_set_param(self.param[flag],5,1.0,0.0,2.0,True)

			if (flag == 'vt'):
				self.lmfit_set_param(self.param[flag],6,100.0,0.0,np.inf,True)
			else:
				self.lmfit_set_param(self.param[flag],6,3.0,0.0,10.0,True)

			self.lmfit_set_param(self.param[flag],7,0.1,0.0,0.5,True)
			self.lmfit_set_param(self.param[flag],8,1.5,1.01,10.0,True)

		if (fit_type == 3):
			if (flag == 'te' or flag == 'ti'):
				self.lmfit_set_param(self.param[flag],2,1.,0.3,20,True)
				self.lmfit_set_param(self.param[flag],1,0.15,0.05,0.2,True)
			elif (flag == 'ne'):
				self.lmfit_set_param(self.param[flag],2,1.,0.5,10,True)
				self.lmfit_set_param(self.param[flag],1,0.5,0.2,1.0,True)
			elif (flag == 'vt'):
				self.lmfit_set_param(self.param[flag],2,100,0,500,True)
				self.lmfit_set_param(self.param[flag],1,50,0,100,True)
			
			self.lmfit_set_param(self.param[flag],3,-0.1,-np.inf,0.,True)
			self.lmfit_set_param(self.param[flag],4,1,-np.inf,np.inf,True)				
			self.lmfit_set_param(self.param[flag],5,1,-np.inf,np.inf,True)
			if (flag == 'ti' or flag == 'vt'): self.lmfit_set_param(self.param[flag],6,0.05,0.04,0.15,True)
			else: self.lmfit_set_param(self.param[flag],6,0.05,0.03,0.1,True)	
			self.lmfit_set_param(self.param[flag],7,0.98,0.88,1.04,True)

		if (fit_type == 4 or fit_type == 9):
			if (flag == 'te' or flag == 'ti'):
				self.lmfit_set_param(self.param[flag],1,0.15,0.05,0.3,True)
				self.lmfit_set_param(self.param[flag],2,0.5,0.15,5.0,True)
			elif (flag == 'ne'):
				self.lmfit_set_param(self.param[flag],1,0.5,0.3,1.0,True)
				self.lmfit_set_param(self.param[flag],2,1.0,0.15,5.0,True)
			elif (flag == 'vt'):
				self.lmfit_set_param(self.param[flag],1,50,0.0,100,True)
				self.lmfit_set_param(self.param[flag],2,100,0.0,np.inf,True)

			if (flag == 'ti' or flag == 'vt'): self.lmfit_set_param(self.param[flag],3,0.05,0.04,0.15,True)
			else: self.lmfit_set_param(self.param[flag],3,0.05,0.03,0.1,True)

			if (flag == 'vt'):
				self.lmfit_set_param(self.param[flag],4,200,0.0,np.inf,True)
			else:
				self.lmfit_set_param(self.param[flag],4,3.0,0.2,10.0,True)
			coef5 = 4.
			if flag=='ne': coef5 = 2.5
			self.lmfit_set_param(self.param[flag],5,1.3,1.01,coef5,True)
			self.lmfit_set_param(self.param[flag],6,2.0,1.01,coef5,True)

		if (fit_type == 6):
			if (flag == 'te' or flag == 'ti'):
				self.lmfit_set_param(self.param[flag],1,0.15,0.05,0.3,True)
				self.lmfit_set_param(self.param[flag],2,0.5,0.15,5.0,True)
			elif (flag == 'ne'):
				self.lmfit_set_param(self.param[flag],1,0.5,0.3,1.0,True)
				self.lmfit_set_param(self.param[flag],2,1.0,0.15,5.0,True)
			elif (flag == 'vt'):
				self.lmfit_set_param(self.param[flag],1,50,0.0,100.,True)
				self.lmfit_set_param(self.param[flag],2,100,0.0,np.inf,True)

			if (flag == 'ti' or flag == 'vt'): self.lmfit_set_param(self.param[flag],3,0.05,0.04,0.15,True)
			else: self.lmfit_set_param(self.param[flag],3,0.05,0.03,0.1,True)

			if (flag == 'vt'):
				self.lmfit_set_param(self.param[flag],4,200,0.0,np.inf,True)
			else:
				self.lmfit_set_param(self.param[flag],4,3.0,0.2,10.0,True)

			self.lmfit_set_param(self.param[flag],5,1.3,1.01,4.0,True)
			self.lmfit_set_param(self.param[flag],6,2.0,1.01,4.0,True)
			self.lmfit_set_param(self.param[flag],7,0.98,0.9,0.99,True)		

		return

	def lmfit_init_params(self,flag='none'):

		if flag == 'none':
			self.lmfit_init()
			self.lmfit_init_param(self.fit_opt['func_type']['te'],'te')
			self.lmfit_init_param(self.fit_opt['func_type']['ne'],'ne')
			self.lmfit_init_param(self.fit_opt['func_type']['ti'],'ti')
			self.lmfit_init_param(self.fit_opt['func_type']['vt'],'vt')
			return
		else:
			self.lmfit_init(flag.lower())

		if not flag=='none':
			self.lmfit_init_param(self.fit_opt['func_type'][flag.lower()],flag.lower())
		return

	def lmfit_renew_param(self,fit_type,fixed_width,width,pedmid,fixed_sep,sep,flag):

		if (fit_type == 1):
			if (fixed_sep):
				self.lmfit_set_param(self.param[flag],1,sep, None,None,False)
			else:
				self.lmfit_set_param(self.param[flag],1,None,None,None,True)

		elif (fit_type == 2):
			if fixed_sep:
				self.lmfit_set_param(self.param[flag],1,sep, None,None,False)
			else:
				self.lmfit_set_param(self.param[flag],1,None,None,None,True)
					
			if (fixed_width):
				self.lmfit_set_param(self.param[flag],3,width,None,None,False)
#			else:
#				self.lmfit_set_param(self.param[flag],3,None, None,None,True)

		elif (fit_type == 3):
			if (fixed_sep):
				self.lmfit_set_param(self.param[flag],1,sep, None,None,False)
			else:
				self.lmfit_set_param(self.param[flag],1,None,None,None,True)

			if (fixed_width):
				self.lmfit_set_param(self.param[flag],6,width,None,None,False)
#			else:
#				self.lmfit_set_param(self.param[flag],6,None, None,None,True)

		elif (fit_type == 4 or fit_type == 9):
			if (fixed_sep):
				self.lmfit_set_param(self.param[flag],1,sep, None,None,False)
			else:
				self.lmfit_set_param(self.param[flag],1,None,None,None,True)

			if (fixed_width):
				self.lmfit_set_param(self.param[flag],3,width,None,None,False)
#			else:
#				self.lmfit_set_param(self.param[flag],3,None, None,None,True)

		elif (fit_type == 6):
			if (fixed_sep):
				self.lmfit_set_param(self.param[flag],1,sep, None,None,False)
			else:
				self.lmfit_set_param(self.param[flag],1,None,None,None,True)

			if (fixed_width):
				self.lmfit_set_param(self.param[flag],3,width,None,None,False)
				if (self.fit_opt['func_type']['ti'] == 6 and self.fit_opt['use_ti_width'] and pedmid < 10.):
					self.lmfit_set_param(self.param[flag],7,pedmid,None,None,False)
#			else:
#				self.lmfit_set_param(self.param[flag],3,None, None,None,True)

		return	

	def lmfit_params_trunc(self,flag):

		for i in range(10):
			if not np.isinf(self.param[flag]['val'][i]):
				self.param[flag]['val'][i] = round(self.param[flag]['val'][i],4)
			if not np.isinf(self.param[flag]['min'][i]):
				self.param[flag]['min'][i] = round(self.param[flag]['min'][i],4)
			if not np.isinf(self.param[flag]['max'][i]):
				self.param[flag]['max'][i] = round(self.param[flag]['max'][i],4)				
		return
		
	def lmfit_params(self,fit_type,fixed_width,width,pedmid,fixed_sep,sep,flag):

		if not (fit_type == 5):
			varn = self.func_varn[fit_type-1]
		else:
			return	
		
		p = lmfit.Parameters()
		self.lmfit_renew_param(fit_type,fixed_width,width,pedmid,fixed_sep,sep,flag)	
		#self.lmfit_params_trunc(flag)
		self.lmfit_put_param(self.param[flag],p,varn)

		return (p)

	def lmfit_params2array(self,p,fit_type):
	
		varn = self.func_varn[fit_type-1]	
	
		dat = np.zeros(varn)
		
		for i in range(varn):
		
			s = 'a%i'%(i+1)
			dat[i] = p[s].value
		return dat

	def pwidth2rwidth(self,width,pedmid):

		pedtop = pedmid - 0.5*width
		pedend = pedmid + 0.5*width
		pedmid2 = self.fit_eq['rho_to_psi'](pedmid)
		pedtop2 = self.fit_eq['rho_to_psi'](pedtop)
		pedend2 = self.fit_eq['rho_to_psi'](pedend)

		width2 = (pedend2 - pedtop2)
	
		return width2, pedmid2
		
	def lmfit_fit(self,flag):

		itime = time.time()
		fit_type    = self.fit_opt['func_type'][flag]
		fixed_width = self.fit_opt['width_fix'][flag]
		width       = self.fit_opt['width_val'][flag]
		pedmid      = self.post['pedmid'][flag]
		fixed_sep   = self.fit_opt['sep_fix'][flag]
		sep         = self.fit_opt['sep_val'][flag]
		use_raw     = self.fit_opt['raw_fit'][flag]
		sspline     = self.fit_opt['sspline_order'][flag]
		oli_cutn    = self.fit_opt['oli'][flag]['n']
		use_rho     = self.fit_opt['use_rho'][flag]

		if use_rho:
			if (width >0.):
				width,pedmid = self.pwidth2rwidth(self.fit_opt['width_val'][flag],self.post['pedmid'][flag])
		
		p = self.lmfit_params(fit_type,fixed_width,width,pedmid,fixed_sep,sep,flag)

		if   (fit_type==1):func = self.core_prof
		elif (fit_type==2):func = self.mtanh_prof
		elif (fit_type==3):func = self.tanh_prof
		elif (fit_type==4):func = self.eped_prof
		elif (fit_type==6):func = self.eped_prof2
		elif (fit_type==9):func = self.eped_prof3

		if not use_raw:
			if not self.post['iskfile']: use_raw = True
		else:
			if not self.post['isrfile']: use_raw = False		

		datx = dict()
		daty = dict()
		dats = dict()
		ch_list = []
		use_avg = False
		not_use_avg_all = False
		if (use_raw):
			nch  = -1
			for i in self.__dict__['%s_list'%flag]:
				if (self.fit_opt['file'][flag][i] == None): continue
				nch = nch + 1
				use_avg = self.fit_opt['avg'][flag][i]
				excn    = self.fit_opt['exclude'][flag][i]
				if not excn == '':
					if excn[-1] == ',':	excn = excn[:-1].split(',')
					else:	excn = excn.split(',')
				else: excn = []
				use_std = self.fit_opt['std'][flag][i]['use']
				raw_std = self.fit_opt['std'][flag][i]['raw']
				psi_end = self.fit_opt['psi_end'][flag][i]
				excn    = np.array(excn,dtype=int)
	
				if use_avg:
					datx2 = np.copy(self.__dict__['%s_prof'%flag][i]['xxa'])
					daty2 = np.copy(self.__dict__['%s_prof'%flag][i]['avg'])
					len2  = len(datx2)
					if (use_std and raw_std): dats2 = np.copy(self.__dict__['%s_prof'%flag][i]['avgs']) 
					elif use_std:             dats2 = np.copy(self.__dict__['%s_prof'%flag][i]['avgs2']) 
					else:					  dats2 = np.ones(len2)*0.1

				else:
					datx2 = np.copy(self.__dict__['%s_prof'%flag][i]['xxr'])
					daty2 = np.copy(self.__dict__['%s_prof'%flag][i]['raw'])
					len2  = len(datx2)					
					if (use_std and raw_std): dats2 = np.copy(self.__dict__['%s_prof'%flag][i]['raws'])
					elif use_std:             dats2 = np.copy(self.__dict__['%s_prof'%flag][i]['raws2']) 
					else:					  dats2 = np.ones(len2)*0.1
					excn2 = np.copy(excn)
					for j in range(self.post['prof_dim'][flag][i][1]-1):
						excn2 = excn2+self.post['prof_dim'][flag][i][0]
						excn = np.hstack([excn,excn2])

				excn    = np.array(excn,dtype=int)
				dati2 = np.array(np.linspace(0,len2-1,len2),dtype='int')
				dati2 = np.delete(dati2,excn)

				if not (psi_end == 0.0):
					ind = np.where(datx2[dati2] <= psi_end)
					dati2 = dati2[ind]

				if len(dati2) ==0: 
					self.post['fit_ind'][flag][i] = []
					self.post['oli_ind'][flag][i] = [[],[],[]]
					continue
				datx[nch] = np.copy(datx2[dati2])
				daty[nch] = np.copy(daty2[dati2])
				dats[nch] = np.copy(dats2[dati2])
				len2      = len(dats[nch])
				self.post['fit_ind'][flag][i] = dati2
				self.post['oli_ind'][flag][i] = [datx[nch],daty[nch],dats[nch]]
				ch_list.append(i)

				if not use_avg: not_use_avg_all = True

		else:
				nch     = 0
				ch_list =['kfile']
				datx[0] = np.copy(self.__dict__['%s_prof'%flag]['kfile']['xxr'])
				daty[0] = np.copy(self.__dict__['%s_prof'%flag]['kfile']['raw'])
				len2    = len(datx)
				dats[0] = np.ones(len2)
				self.post['fit_ind'][flag]['kfile'] =np.linspace(0,len2-1,len2)

		nch = nch + 1

		if (fit_type == 5 or fit_type == 7 or fit_type == 8):
			if (sspline == -1): 
				if fit_type==5: sspline = len2
				else: sspline = 5

			xx = np.copy(datx[0])
			yy = np.copy(daty[0])
			ss = np.copy(dats[0])/(1.e-4+self.fit_opt['weight'][flag][ch_list[0]])
			for i in range(nch-1):
				xx = np.hstack([xx,datx[i+1]])
				yy = np.hstack([yy,daty[i+1]])
				ss = np.hstack([ss,dats[i+1]/(1.e-4+self.fit_opt['weight'][flag][ch_list[i+1]])])

			if self.fit_opt['sep_fix'][flag]:
				xx = np.append(xx,1.)
				yy = np.append(yy,self.fit_opt['sep_val'][flag])
				ss = np.append(ss,1.e-3)


			for i in range(len(xx)):
				if (ss[i] == 0.0):    tempv = 1.e-4 * max(ss)
				else:		          tempv = ss[i]
				if tempv==0.0:        tempv = 0.1
					
				ss[i] = 1./tempv

			ind = np.argsort(xx)
			xx = xx[ind]
			yy = yy[ind]
			ss = ss[ind]
			if fit_type == 7: sspline = min(5,int(sspline))
			kth=3
			if fit_type == 5: sfit = smoothspline(xx,yy,w=ss,s=sspline)		# Do spline
			if fit_type == 7: sfit = self.spline_fit(xx,yy,ss,max(sspline,1),k=kth)
			if fit_type == 8: sfit = self.csaps_fit(xx,yy,ss,max(sspline,1))			

			self.post['popt'][flag] = np.ones(6)
			self.__dict__['%s_prof'%flag]['fit1'] = sfit(self.fit_eq['psin1'])
			self.__dict__['%s_prof'%flag]['fit2'] = sfit(self.fit_eq['psin2'])

			prof_ind = np.where(~np.isnan(self.__dict__['%s_prof'%flag]['fit2'] ))
			len3     = len(self.__dict__['%s_prof'%flag]['fit2'])
			if not len(prof_ind) == 0:
				psin = self.fit_eq['psin2'][prof_ind]
				prof = self.__dict__['%s_prof'%flag]['fit2'][prof_ind]
				sfit2 = interp1d(psin,prof,'cubic')
				for i in range(len3):
					try:    self.__dict__['%s_prof'%flag]['fit2'] [i] = sfit2(self.fit_eq['psin2'][i])
					except: self.__dict__['%s_prof'%flag]['fit2']  = 0.e0
		
			ind = np.where(self.fit_eq['psin2'] <= 1.0)
			ind = np.where(self.__dict__['%s_prof'%flag]['fit2'][ind] < 1.e4)

			for i in range(len3):
				try:
					if (self.__dict__['%s_prof'%flag]['fit2'][i] > 1.1 * max(self.__dict__['%s_prof'%flag]['fit2'][ind]) or self.__dict__['%s_prof'%flag]['fit2'][i] < 0.0):
						self.__dict__['%s_prof'%flag]['fit2'][i] = 0.e0
				except: pass

			if fit_type < 8: self.post['chi'][flag] = sfit.get_residual()/len2
			else: self.post['chi'][flag] = 1.
			self.post['popte'][flag] = np.zeros(10)
			if fit_type == 7:
				line = '%7.3f'%0.
				for k in range(len(self.spline_knots)): line = line+'%7.3f'%self.spline_knots[k]
				line = line + '%7.3f'%1.0

			self.print('=----------------------------------------------------------------------------=')	
			if fit_type == 5: self.print('>>> Fitting result for %s, SMOOTHED_SPLINE'%(flag.upper()))
			else:  self.print('>>> Fitting result for %s, SKNOTS_SPLINE'%(flag.upper()))
			if fit_type == 7: 
				self.print('>>> Fitting with %i-th order, %i knots'%(kth,max(sspline,1)+2))
				self.print('>>> Internal knots - %s'%line)
			self.print('>>> XI2 = %6.3f'%(self.post['chi'][flag]))
			self.print('>>> Elapsed time %6.3f(s) '%(time.time()-itime))
			self.print('=----------------------------------------------------------------------------=')
			self.print('=----------------------------------------------------------------------------=')		

			return

		FIT_FUNC = 'CORE'
		if   fit_type==2: FIT_FUNC = 'MTANH'
		elif fit_type==3: FIT_FUNC = 'PTANH'
		elif fit_type==4: FIT_FUNC = 'EPED'
		elif fit_type==6: FIT_FUNC = 'EPED2'
		elif fit_type==9: FIT_FUNC = 'EPED3'

		self.post['popte'][flag] = np.zeros(10)

		if (use_raw and not_use_avg_all):
			for i in range(oli_cutn):
				result = lmfit.minimize(self.residual,p,args=(datx,daty,dats,ch_list,func,flag))
					
				p = result.params
				
				for i in datx.keys():
					yf = func(datx[i],**result.params)
					ytt = (yf - daty[i])**1 / (dats[i]+1.e-8) ** 1

					if (self.fit_opt['oli'][flag][ch_list[i]]['use'] and not self.fit_opt['avg'][flag][ch_list[i]]):
						oli_cut = self.fit_opt['oli'][flag][ch_list[i]]['per']
					else:
						oli_cut = 1.0

					if (oli_cut == 0.0):
						oli_cut2 = 95
					else:
						oli_cut2 = oli_cut * 100

					yn = np.percentile(ytt,oli_cut2)>=ytt

					datx[i] = datx[i][yn]
					daty[i] = daty[i][yn]
					dats[i] = dats[i][yn]
					self.post['oli_ind'][flag][ch_list[i]] = [datx[i],daty[i],dats[i]]
			if not self.noprint:
				self.print('=----------------------------------------------------------------------------=')	
				self.print('>>> Fitting result for %s, FIT_FUNC = %s'%(flag.upper(),FIT_FUNC))			
				self.print('>>> XI2 = %6.3f, RED_XI2 = %6.4e'%(result.chisqr,result.redchi))
				self.print('>>> Elapsed time %6.3f(s) '%(time.time()-itime))
				self.print('=----------------------------------------------------------------------------=')
				self.pretty_print(result.params,self.post['popte'][flag])
				self.print('=----------------------------------------------------------------------------=')

		else:
			result = lmfit.minimize(self.residual,p,args=(datx,daty,dats,ch_list,func,flag))
			if not self.noprint:
				self.print('=----------------------------------------------------------------------------=')
				self.print('>>> Fitting result for %s, FIT_FUNC = %s'%(flag.upper(),FIT_FUNC))
				self.print('>>> XI2 = %6.3f, RED_XI2 = %6.4e'%(result.chisqr,result.redchi))
				self.print('>>> Elapsed time %6.3f(s) '%(time.time()-itime))
				self.print('=----------------------------------------------------------------------------=')
				self.pretty_print(result.params,self.post['popte'][flag])
				self.print('=----------------------------------------------------------------------------=')

		self.post['popt'][flag] = self.lmfit_params2array(result.params,fit_type)
		self.__dict__['%s_prof'%flag]['fit1'] = func(self.fit_eq['psin1'],**result.params)
		self.__dict__['%s_prof'%flag]['fit2'] = func(self.fit_eq['psin2'],**result.params)		
		
		try:
			self.post['chi'][flag] = result.chisqr
		except:
			self.post['chi'][flag] = 1.0
		
		return
		
	def eped_prof(self, x, a1, a2, a3, a4, a5, a6):

		y = a1
		y = y + a2*(np.tanh(1) - np.tanh((x - 1.0 + 0.5*(a3))/(a3)*2.0))
		yt = (x/(1.0-a3)) ** (a5)
		y = y + a4 * ((abs(1-yt)) **(a6)) * 0.5 * (1.0 + np.sign(1.0-a3-x)) #(abs)
		
		return y

	def eped_prof3(self, x, a1, a2, a3, a4, a5, a6):
			
		y = a2*(np.tanh(1) - np.tanh((x - 1.0 + 0.5*(a3))/(a3)*2.0))
		alp1 = 1.1 + abs(a5) #3.*0.5*(1.+np.tanh((1-a5)/10.))
		alp2 = 1.1 + abs(a6) #3.*0.5*(1.+np.tanh((1-a6)/10.))
		yt = abs(x/(1.0-a3)) ** (alp1)
		y = y + a1 * a4 * ((abs(1-yt)) **(alp2)) * 0.5 * (1.0 + np.sign(1.0-a3-x))
		y = y + a1 * a2

		return y 

	def eped_prof2(self, x, a1, a2, a3, a4, a5, a6, a7):

		y = a1
		y = y + a2*(np.tanh((1. - a7)/(a3)*2.0) - np.tanh((x - a7)/(a3)*2.0))
		yt = (x/(a7-0.5*a3)) ** (a5)
		y = y + a4 * ((abs(1-yt)) **(a6)) * 0.5 * (1.0 + np.sign(a7-0.5*a3-x))
		
		return y

	def core_prof(self, x, a1, a2, a3, a4): # Core gaussian

		y = a1
		yt = x ** (a3)
		y = y + a2 * abs(1 - yt)**(a4)

		return y

	def mtanh(self,x,b): # mtanh
		
		y = ((1.0 + b*x)*np.exp(x) - np.exp(-x)) / (np.exp(x) + np.exp(-x))
		return y

	def mtanh_prof(self, x, a1, a2, a3, a4, a5, a6, a7, a8): #mtanh_prof

		y = (a2 - a1)/ 2.0 * (self.mtanh((a4-x)/2/(a3/4+1.e-7),a5) + 1.0) + a1

		y = y + (a6 - y) * np.exp(-1.*(x/(a7+0.2))**(a8))

		return y

	def tanh_prof(self,x,a1,a2,a3,a4,a5,a6,a7): #tanh_prof

		y = (a2 -  a1) * (1.0 + a3*x + a4 * x**2 + a5 * x**3) * 0.5 * (1.0 - np.tanh((x-a7)/a6*2.)) + a1

		return y

	def spline_knot_params(self,knots_num):

		p = lmfit.Parameters()
		for i in range(knots_num+1):
				val = 1./float(knots_num+1)
				p.add('d%i'%i,value=val,vary=True)
		return p

	def spline_residual(self,p,x,y,w,knots_num,delx,k=2):

		t = []
		t2 = []
		sum1 = 0.
		for i in range(knots_num+1): sum1 = sum1 + abs(p['d%i'%i].value)
		cmul = 0.8

		sum2 = 0.
		dels = sum1 * delx / (1. - knots_num*delx - cmul*delx)
		sum1 = sum1 + dels*(knots_num+cmul)
		for i in range(knots_num):
			if i == 0: sum2 = sum2 + abs(p['d%i'%i])+dels*(cmul)
			else: sum2 = sum2 + abs(p['d%i'%i])+dels
			t.append(sum2/sum1)
			t2.append(p['d%i'%i])

		ind = np.where(x<1.05)
		yff = lspline(x[ind],y[ind],t,w=w[ind],bbox=[0.,1.05],ext=3,k=k)
		dy0 = yff.derivatives(0.)[1]
		yf = yff(x)
		yf = np.round(yf,self.fit_trunc)
		yt = ((yf-y)**2) * (w**2) 
		yt = yt * (1.+ 5.e2*(1+max(-1,abs(dy0+0.2)-1)))

		self.spline_delx = dels
		self.spline_knots= t
		return (yt)

	def spline_fit(self,datx,daty,datw,knots_num,k=2):

		p = self.spline_knot_params(knots_num)

		result = lmfit.minimize(self.spline_residual,p,args=(datx,daty,datw,knots_num,0.1,k)) # 0.1 min delx b/w knots
		p = result.params
		ind = np.where(datx<1.05)	
		yff = lspline(datx[ind],daty[ind],self.spline_knots,w=datw[ind],bbox=[0.,1.05],ext=3,k=k)

		return yff

	def csaps_fit(self,datx,daty,datw,sfactor=5): #new smooth spline
		smooth = 1. - 10**(-sfactor)

		xx2 = np.linspace(0.,1.4,1401)
		yy2 = np.copy(xx2)
		xx  = np.linspace(min(datx),max(datx),401)
		yy  = csaps(datx,daty,xx,weights=datw, smooth=smooth)
		yyf = interp1d(xx,yy,'cubic')
		d1  = (yy[1]-yy[0])/(xx[1]-xx[0])
		d2  = (yy[-1]-yy[-2])/(xx[-1]-xx[-2])
		for i in range(1401):
			if xx2[i] <= xx[0]:  yy2[i] = d1*(xx2[i]-xx[0])+yy[0]
			elif xx2[i]>=xx[-1]: yy2[i] = d2*(xx2[i]-xx[-1])+yy[-1]
			else: yy2[i] = yyf(xx2[i])
			if yy2[i] <0.01: yy2[i] = 0.01
		yyf2 = interp1d(xx2,yy2,'cubic')
		return yyf2
		
	def write_chease_rot(self):
		
		f4 = open('PROFILES/chease_vtor','w')

		ppp = np.linspace(0,1.0,301)	
		vtf = interp1d(self.fit_eq['psin2'],self.vt_prof['fit2'],'cubic')

		f4.write('%i\n'%(len(ppp)))
		f4.write('PSIN RR[m] VTOR[km/s]\n')
		for i in range(len(ppp)):
			psit = ppp[i]
			if (self.fit_opt['use_rho']['vt']):
				psit = self.fit_eq['rho_to_psi'](psit)
			f4.write('%9.6f\t%9.6f\t%9.6f\n'%(psit,self.fit_eq['Rf'](ppp[i]),vtf(ppp[i])))

		f4.close()
		
		return
				
	def draw_plot_part(self,fig,flag):

		colormap1 = ['orange','purple','magenta']
		colormap2 = ['green','yellowgreen','navy']

		if (flag == 'te'):
			title = '$T_e$ [keV]'
		elif (flag == 'ne'):
			title = '$n_e$ [$10^{19}$/m3]'
		elif (flag == 'ti'):
			title = '$T_i$ [keV]'
		elif (flag == 'vt'):
			title = '$V_\phi$ [km/s]'

		WPED = 0.
		if (self.fit_opt['func_type'][flag] == 4 or self.fit_opt['func_type'][flag] == 2):
			WPED = self.post['popt'][flag][2]
			PMID = 1.0 - 0.5*WPED
			if self.fit_opt['use_rho'][flag]: WPED,PMID = self.pwidth2rwidth(WPED,PMID)

		if (self.fit_opt['func_type'][flag] ==6):
			WPED = self.post['popt'][flag][2]
			PMID = self.post['popt'][flag][6]
			if self.fit_opt['use_rho'][flag]: WPED,PMID = self.pwidth2rwidth(WPED,PMID)

		if (self.fit_opt['func_type'][flag] == 3):
			WPED = self.post['popt'][flag][5]
			PMID = 1.0-0.5*WPED		
			if self.fit_opt['use_rho'][flag]: WPED,PMID = self.pwidth2rwidth(WPED,PMID)	

		fig.set_ylabel(title)
		if (WPED > 0.): 
			if self.post['time'] > 0: title = title + ' %i[ms] - $W_{\psi,ped}$ = %4.3f'%(self.post['time'],WPED)
			else: title = title + ' - $W_{\psi,ped}$ = %4.3f'%WPED
		
		fig.set_title(title)
		if not (self.fit_opt['use_rho'][flag]):
			fig.set_xlabel('$\psi_n$ [a.u]')
		else:
			fig.set_xlabel('$\\rho_t$ [a.u]')

		plegend = []
		llegend = []

		if not (self.fit_opt['file'][flag]['kfile'] == None):
			line1, = fig.plot(self.__dict__['%s_prof'%flag]['kfile']['xxr'],self.__dict__['%s_prof'%flag]['kfile']['raw'],'gray',linestyle= '--')
			plegend.append('Saved Fitted prof')
			llegend.append(line1)

		if np.sum(self.__dict__['%s_prof'%flag]['fit_old']) > 0:
			if self.fit_opt['use_rho'][flag]:
				line2, = fig.plot(self.fit_eq['psi_to_rho'](self.fit_eq['psin2']),self.__dict__['%s_prof'%flag]['fit_old'],'blue',linestyle= '--')
			else:
				line2, = fig.plot(self.fit_eq['psin2'],self.__dict__['%s_prof'%flag]['fit_old'],'blue',linestyle= '--')
			plegend.append('Pre Fitted prof')
			llegend.append(line2)		

		line2, = fig.plot(self.fit_eq['psin2'],self.__dict__['%s_prof'%flag]['fit2'],'red')
		plegend.append('New Fitted prof')
		llegend.append(line2)

		raw_fit  = self.fit_opt['raw_fit'][flag]
		var_crit = self.device['var_crit'][flag]
	
		nch = -1
		for k in self.__dict__['%s_list'%flag]:
			if self.fit_opt['file'][flag][k] == None: continue
			nch = nch + 1
			plegend.append('%s-data'%k.upper())
			use_avg = self.fit_opt['avg'][flag][k]
			use_std = self.fit_opt['std'][flag][k]['use']
			raw_std = self.fit_opt['std'][flag][k]['raw']
			use_oli = self.fit_opt['oli'][flag][k]['use']

			if not use_avg:

				datx = self.__dict__['%s_prof'%flag][k]['xxr']
				daty = self.__dict__['%s_prof'%flag][k]['raw']
				if raw_std: dats = self.__dict__['%s_prof'%flag][k]['raws']
				else:		dats = self.__dict__['%s_prof'%flag][k]['raws2']

				len1 = len(datx)
				ind1 = np.array(np.linspace(0,len1-1,len1),dtype='int')
				if self.fit_opt['raw_fit'][flag]:
					ind2 = np.delete(ind1,self.post['fit_ind'][flag][k])
					ind3 = self.post['fit_ind'][flag][k]
				else:
					ind2 = ind1
					ind3 = []

				if (use_std): line3 = fig.errorbar(datx[ind2],daty[ind2], yerr = dats[ind2],fmt=self.plot_shape1[nch],markersize='5',c='gray',ecolor='gray',capthick=2)
				else:		  line3 = fig.scatter(datx[ind2],daty[ind2],s=20,marker=self.plot_shape1[nch],c='gray')

				for i in range(self.post['prof_dim'][flag][k][0]):
					dat = np.reshape(daty,(self.post['prof_dim'][flag][k][1],self.post['prof_dim'][flag][k][0]))
					fig.text(self.__dict__['%s_prof'%flag][k]['xxa'][i],self.__dict__['%s_prof'%flag][k]['avg'][i],str(i+1),color=colormap1[nch],clip_on=True)

				if (raw_fit):
					if use_oli: cmap = copy.deepcopy(colormap2)
					else: cmap = copy.deepcopy(colormap1)

					if (use_std): line3 = fig.errorbar(datx[ind3],daty[ind3],yerr=dats[ind3],fmt=self.plot_shape1[nch],markersize='5',c=cmap[nch],ecolor=cmap[nch],capthick=2)
					else:		  line3 = fig.scatter(datx[ind3],daty[ind3],marker=self.plot_shape1[nch],s=20,c=cmap[nch])

					if (use_oli):
						oli_intx = self.post['oli_ind'][flag][k][0]
						oli_inty = self.post['oli_ind'][flag][k][1]
						oli_ints = self.post['oli_ind'][flag][k][2]

						if (use_std): line3 = fig.errorbar(oli_intx,oli_inty,yerr = oli_ints,fmt=self.plot_shape1[nch],markersize='5',c=colormap1[nch],ecolor=colormap1[nch],capthick=2)
						else:		  line3 = fig.scatter(oli_intx, oli_inty,s=20,marker=self.plot_shape1[nch],c=colormap1[nch])
			else:

				datx = self.__dict__['%s_prof'%flag][k]['xxa']
				daty = self.__dict__['%s_prof'%flag][k]['avg']

				if raw_std: dats = self.__dict__['%s_prof'%flag][k]['avgs']
				else:		dats = self.__dict__['%s_prof'%flag][k]['avgs2']			

				len1 = len(datx)
				ind1 = np.array(np.linspace(0,len1-1,len1),dtype='int')
				if self.fit_opt['raw_fit'][flag]:
					ind2 = np.delete(ind1,self.post['fit_ind'][flag][k])
					ind3 = self.post['fit_ind'][flag][k]						
				else:
					ind2 = ind1
					ind3 = []

				if not (use_std):
					line3 = fig.scatter(datx[ind2],daty[ind2],s=15,marker=self.plot_shape2[nch],c='gray')
					line3 = fig.scatter(datx[ind3],daty[ind3],s=15,marker=self.plot_shape2[nch],c='g')
					if (raw_fit):	
						line3 = fig.scatter(datx[ind2],daty[ind2],s=15,marker=self.plot_shape2[nch],c='gray')
						line3 = fig.scatter(datx[ind3],daty[ind3],s=15,marker=self.plot_shape2[nch],c='g')
	
				else:
					line3 = fig.errorbar(datx[ind2],daty[ind2], yerr = dats[ind2],fmt=self.plot_shape2[nch],markersize='3',c='gray',ecolor='gray',capthick=2)
					line3 = fig.errorbar(datx[ind3],daty[ind3], yerr = dats[ind3],fmt=self.plot_shape2[nch],markersize='3',c=colormap1[nch],ecolor=colormap1[nch],capthick=2)
					if (raw_fit):	
						line3 = fig.errorbar(datx[ind2],daty[ind2], yerr = dats[ind2],fmt=self.plot_shape2[nch],markersize='3',c='gray',ecolor='gray',capthick=2)
						line3 = fig.errorbar(datx[ind3],daty[ind3], yerr = dats[ind3],fmt=self.plot_shape2[nch],markersize='3',c=colormap1[nch],ecolor=colormap1[nch],capthick=2)

				for i in range(len(datx)):
					fig.text(datx[i],daty[i],str(i+1),color=colormap1[nch],clip_on=True)
				
			llegend.append(line3)		
					
		line4 = fig.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
		if not flag.lower() == 'vt':	line5 = fig.axhline(y=0.0,color='goldenrod',linestyle='--',linewidth=1.0)
		plegend.append('$\chi^2$ = %4.2e'%self.post['chi'][flag])

		fig.set_xlim(self.fit_opt['plot'][flag]['xmin'],self.fit_opt['plot'][flag]['xmax'])
		ind = np.where(datx<0.3)
		if min(datx) > 0.3: yu = max(daty)*1.1
		else:	yu = max(daty[ind])*1.1
		if flag.lower()=='vt':			
			yl = max(min(daty),-50.)
		else:
			if min(daty) > -0.5: yl = -0.5
			else:	yl = max(min(daty),-0.5)
		yu = max(yu,2.0)
		if self.fit_opt['plot'][flag]['ymin'] == -1:
			self.fit_opt['plot'][flag]['ymin'] = yl
		if self.fit_opt['plot'][flag]['ymax'] == -1:
			self.fit_opt['plot'][flag]['ymax'] = yu

		fig.set_ylim(self.fit_opt['plot'][flag]['ymin'],self.fit_opt['plot'][flag]['ymax'])

		llegend.append(line4)
		fig.legend(llegend,plegend)
		
		return

	def draw_plot_part2(self,fig,flag):

		if (flag == 'te'):
			title = 'd$T_e$/d$\psi_n$'
		elif (flag == 'ne'):
			title = 'd$n_e$/d$\psi_n$'
		elif (flag == 'ti'):
			title = 'd$T_i$/d$\psi_n$'
		elif (flag == 'vt'):
			title = 'd$V_{\phi}$/d$\psi_n$'
		

		if not (self.fit_opt['use_rho'][flag]):
			fig.set_xlabel('$\psi_n$ [a.u]')
		else:
			fig.set_xlabel('$\\rho_t$ [a.u]')	
		if self.post['time'] > 0: title = title + ' %i[ms] '%(self.post['time'])
		fig.set_title(title)

		plegend = []
		llegend = []
		deps = 1.e-5

		if (self.post['iskfile'][flag]):

			fitf = interp1d(self.__dict__['%s_prof'%flag]['kfile']['xxr'],self.__dict__['%s_prof'%flag]['kfile']['raw'],'quadratic')
			fitd = np.copy(self.__dict__['%s_prof'%flag]['kfile']['raw'])
			for i in range(len(fitd)-2):
				xx = self.__dict__['%s_prof'%flag]['kfile']['xxr'][i+1]
				fitd[i+1] = (fitf(xx+deps)-fitf(xx-deps))/2./deps
			fitd[0] = (fitf(self.__dict__['%s_prof'%flag]['kfile']['xxr'][0]+deps)-fitf(self.__dict__['%s_prof'%flag]['kfile']['xxr'][0]))/deps
			fitd[-1] = (fitf(self.__dict__['%s_prof'%flag]['kfile']['xxr'][-1])-fitf(self.__dict__['%s_prof'%flag]['kfile']['xxr'][-1]-deps))/deps

			line1, = fig.plot(self.__dict__['%s_prof'%flag]['kfile']['xxr'],fitd,'gray',linestyle= '--')
			plegend.append('Saved Fitted prof')
			llegend.append(line1)

		if np.sum(self.__dict__['%s_prof'%flag]['fit_old']) > 0:
			if not self.fit_opt['use_rho'][flag]:
				fitf = interp1d(self.fit_eq['psin2'],self.__dict__['%s_prof'%flag]['fit_old'],'quadratic')
			else:
				fitf = interp1d(self.fit_eq['psi_to_rho'](self.fit_eq['psin2']),self.__dict__['%s_prof'%flag]['fit_old'],'quadratic')
			fitd = np.copy(self.fit_eq['psin2'])
			for i in range(len(self.fit_eq['psin2'])-2):
				xx = self.fit_eq['psin2'][i+1]
				fitd[i+1] = (fitf(xx+deps)-fitf(xx-deps))/2./deps
			fitd[0] = (fitf(self.fit_eq['psin2'][0]+deps)-fitf(self.fit_eq['psin2'][0]))/deps
			fitd[-1] = (fitf(self.fit_eq['psin2'][-1])-fitf(self.fit_eq['psin2'][-1]-deps))/deps		

			line2, = fig.plot(self.fit_eq['psin2'],fitd,'blue',linestyle= '--')
			plegend.append('Pre Fitted prof')
			llegend.append(line2)

		if (self.post['didfit'][flag]):
			fitf = interp1d(self.fit_eq['psin2'],self.__dict__['%s_prof'%flag]['fit2p'],'quadratic')
			fitd = np.copy(self.fit_eq['psin2'])
			for i in range(len(self.fit_eq['psin2'])-2):
				xx = self.fit_eq['psin2'][i+1]
				fitd[i+1] = (fitf(xx+deps)-fitf(xx-deps))/2./deps
			fitd[0] = (fitf(self.fit_eq['psin2'][0]+deps)-fitf(self.fit_eq['psin2'][0]))/deps
			fitd[-1] = (fitf(self.fit_eq['psin2'][-1])-fitf(self.fit_eq['psin2'][-1]-deps))/deps		

			line2, = fig.plot(self.fit_eq['psin2'],fitd,'r')
			plegend.append('New Fitted prof')
			llegend.append(line2)
							
		line4 = fig.axvline(x=1.0,color='goldenrod',linestyle='--',linewidth=1.0)
		line5 = fig.axhline(y=0.0,color='goldenrod',linestyle='--',linewidth=1.0)
		plegend.append('Separatrix')
		llegend.append(line4)
		fig.legend(llegend,plegend)
		
		return

	def draw_plot(self,fig_ex=None,second=False):

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
				ax1.cla(); ax2.cla(); ax3.cla(); ax4.cla();

		for k in range(4):
			flag = self.prof_list[k]
			ax   = fig.axes[k]

			if self.post['didfit'][flag]:
				if not second: self.draw_plot_part(ax,flag)
				else:          self.draw_plot_part2(ax,flag)
			else: ax.text(0.44,0.5,'No %s data'%flag.upper(),color='r')

		fig.tight_layout()

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

		zeff  = self.fit_opt['zeff']
		zimp  = self.fit_opt['zimp']
		aimp  = self.fit_opt['aimp']
		amain = self.fit_opt['amain']
		ni_ne = 1.0 - (zeff-1.0)/zimp

		ppp  = np.linspace(0,1,101)
		ppp2 = np.linspace(0,1,401)
		ppp3 = np.linspace(0,1.2,101)	

		pp1 = interp1d(self.fit_eq['psin2'],self.ne_prof['fit2p'],'cubic')
		pp2 = interp1d(self.fit_eq['psin2'],self.te_prof['fit2p'],'cubic')
		pp3 = interp1d(self.fit_eq['psin2'],self.ti_prof['fit2p'],'cubic')
		pp4 = interp1d(self.fit_eq['psin2'],self.vt_prof['fit2p'],'cubic')

		f = open('PROFILES/chease_kinprof_fit','w')	
		f.write('401\n')
		f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(zeff,zimp,amain,aimp))

		for i in range(401):
			#if no_vt: f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(ppp2[i],pp2(ppp2[i]),pp1(ppp2[i]),pp3(ppp2[i]),pp1(ppp2[i])*ni_ne))
			f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(ppp2[i],pp2(ppp2[i]),pp1(ppp2[i]),pp3(ppp2[i]),pp1(ppp2[i])*ni_ne,pp4(ppp2[i])))
		f.close()


		f = open('PROFILES/scaled_factor','w')
		f.write('%f\n'%self.post['scaled_density'])
		f.close()

		f = open('PROFILES/NE_non_scaled.dat','w')
		f.write('--- Fitted profile by FUSMA fittool --- \n Fitted Values: X_Norm_Rho, Psi_Norm, NE[1E18/m3]\n')
		for i in range(101):
			psi_temp = self.fit_eq['rho_to_psi'](ppp[i])
			f.write('%9.6f\t%9.6f\t%9.6f\n'%(ppp[i],psi_temp,pp1(psi_temp)/self.factor['ne']/self.post['scaled_density']))
		f.close()			

		self.write_kprofile(ppp)
		self.write_extended_kinprof()
		if (self.fit_opt['ped_scan_fit']): self.write_epedprof()

		self.write_raw_kinprof()
		self.write_avg_dat()

		#write psi-rho mapping
		f = open('PROFILES/psi_rho.dat','w')
		f.write('%i\n'%len(ppp3))
		f.write('psin[a.u]\trho[a.u]\n')
		for i in range(len(ppp3)):
			f.write('%9.6f\t%9.6f\n'%(ppp3[i],self.fit_eq['psi_to_rho'](ppp3[i])))
		f.close()
		
		#write fit_coefs
		self.write_fit_coefs()
		return

	def write_kprofile(self,ppp):

		for i in self.prof_list:
			no_fit = not self.post['didfit'][i]
			use_rho= self.fit_opt['use_rho'][i]
			filename = 'PROFILES/%s_fit.dat'%i.upper()

			if no_fit: continue

			f = open(filename,'w')
			f.write('--- Fitted profile by FUSMA fittool --- \n')
			f.write('Fitted Values: X_Norm_Rho, Psi_Norm, %s\n'%self.flag[i])

			pp = interp1d(self.fit_eq['psin2'],self.__dict__['%s_prof'%i]['fit2p'],'cubic')

			for j in range(101):
				psi_temp = self.fit_eq['rho_to_psi'](ppp[j])
				f.write ('%9.6f\t%9.6f\t%9.6f\n'%(ppp[j],psi_temp,pp(psi_temp)/self.factor[i]))

			f.close()
		return

	def write_extended_kinprof(self):

		try:	R1 = self.make_psiRZ_extended(self.fit_eq['psin2'],0.,True)
		except:	R1 = np.ones(len(self.fit_eq['psin2']))
		try:	R2 = self.make_psiRZ_extended(self.fit_eq['psin2'],0.,False)
		except:	R2 = np.ones(len(self.fit_eq['psin2']))

		dpdr = (1.0 - self.fit_eq['psi_to_rho'](0.9999))/ (1.0 - 0.9999)

		f = open('PROFILES/chease_kinprof_extended','w')
		f.write('%i\n'%len(self.fit_eq['psin2']))
		f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(0.,0.,0.,0.))
		for i in range(len(self.fit_eq['psin2'])):
			psint = self.fit_eq['psin2'][i]
			f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(psint,self.te_prof['fit2p'][i],self.ne_prof['fit2p'][i],self.ti_prof['fit2p'][i],\
				self.vt_prof['fit2p'][i],R1[i],R2[i]))
		f.close()		

		return

	def write_epedprof(self):

		f4 = open('PROFILES/chease_eped_mod','w')
		
		pedw = abs(self.post['popt']['ne'][2])
		self.post['popt']['ne'][0] = self.post['popt']['ne'][0]*self.post['scaled_density']
		self.post['popt']['ne'][1] = self.post['popt']['ne'][1]*self.post['scaled_density']
		self.post['popt']['ne'][3] = self.post['popt']['ne'][3]*self.post['scaled_density']

		popt_ne = self.post['popt']['ne']
		popt_ti = self.post['popt']['ti']
		popt_te = self.post['popt']['te']

		if (self.fit_opt['use_rho']['ne']):
			popt_temp = np.copy(popt_ne)
			popt_ne = self.modify_rho_to_psi_eped(popt_temp)
			pedw = abs(popt_ne[2])
		if self.fit_opt['func_type']['ne'] ==6:
			f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(popt_ne[1],popt_ne[0],popt_ne[3],popt_ne[6],pedw,popt_ne[6]-0.5*pedw,abs(popt_ne[4]),abs(popt_ne[5])))
		else:
			f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(popt_ne[1],popt_ne[0],popt_ne[3],1-0.5*pedw,pedw,1-pedw,abs(popt_ne[4]),abs(popt_ne[5])))
		if (self.fit_opt['use_rho']['ne']): popt_ne= np.copy(popt_temp)

		repeat = 1
		if not (self.fit_opt['use_ti_eped']): repeat = 2
		
		for i in range(repeat):
			pedw = abs(popt_te[2])
			if (self.fit_opt['use_rho']['te']):
				popt_temp = np.copy(popt_te)
				popt_te = self.modify_rho_to_psi_eped(popt_temp)
				pedw = abs(popt_te[2])
			if self.fit_opt['func_type']['te'] ==6:
				f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(popt_te[1],popt_te[0],popt_te[3],popt_te[6],pedw,popt_te[6]-0.5*pedw,abs(popt_te[4]),abs(popt_te[5])))
			else:
				f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(popt_te[1],popt_te[0],popt_te[3],1-0.5*pedw,pedw,1-pedw,abs(popt_te[4]),abs(popt_te[5])))
			if (self.fit_opt['use_rho']['te']): popt_te = np.copy(popt_temp)
				
		if (self.fit_opt['use_ti_eped']):
			pedw = abs(popt_ti[2])
			if (self.fit_opt['use_rho']['ti']):
				popt_temp = np.copy(popt_ti)
				popt_ti = self.modify_rho_to_psi_eped(popt_temp)
				pedw = abs(popt_ti[2])
			if self.fit_opt['func_type']['ti'] == 6:
				f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(popt_ti[1],popt_ti[0],popt_ti[3],1-0.5*pedw,popt_ti[6],popt_ti[6]-0.5*pedw,abs(popt_ti[4]),abs(popt_ti[5])))
			else:
				f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(popt_ti[1],popt_ti[0],popt_ti[3],1-0.5*pedw,pedw,1-pedw,abs(popt_ti[4]),abs(popt_ti[5])))
			if (self.fit_opt['use_rho']['ti']): popt_ti = np.copy(popt_temp)

		f4.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(self.fit_opt['zeff'],self.fit_opt['zimp'],self.fit_opt['amain'],self.fit_opt['aimp']))
		f4.close()	
		return

	def write_raw_kinprof(self):

		for i in self.prof_list:
			datx = []; daty = []; dats = []; datR = [];

			raw_fit = self.fit_opt['raw_fit'][i]
			no_fit  =not self.post['didfit'][i]

			if not (raw_fit and not no_fit): continue

			use_rho = self.fit_opt['use_rho'][i]
			for j in self.__dict__['%s_list'%i]:
				if not self.fit_opt['file'][i][j] == None:
					try:
						datx = np.hstack([datx,self.__dict__['%s_prof'%i][j]['xxr']])
						daty = np.hstack([daty,self.__dict__['%s_prof'%i][j]['raw']])
						dats = np.hstack([dats,self.__dict__['%s_prof'%i][j]['raws']])
						datR = np.hstack([datR,self.__dict__['%s_prof'%i][j]['datR']])
					except: pass
					f = open('PROFILES/%s_raw_%s.dat'%(i.upper(),j.upper()),'w')			
					try: lend = len(self.__dict__['%s_prof'%i][j]['xxr'])
					except: f.close(); continue;
					f.write('%i\n'%lend)
					datx2 = self.__dict__['%s_prof'%i][j]['xxr']
					daty2 = self.__dict__['%s_prof'%i][j]['raw']
					dats2 = self.__dict__['%s_prof'%i][j]['raws']
					datR2 = self.__dict__['%s_prof'%i][j]['datR']
					if use_rho: 
						xx1 = self.fit_eq['rho_to_psi'](datx2)
						xx2 = self.__dict__['%s_prof'%i][j]['xxr']
					else: 
						xx1 = self.__dict__['%s_prof'%i][j]['xxr']
						xx2 = self.fit_eq['psi_to_rho'](datx2)
					
					for k in range(lend): f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(xx1[k],daty2[k],dats2[k],datR2[k],xx2[k]))
					f.close()

					if self.fit_opt['oli'][i][j]['use']:
						f = open('PROFILES/%s_raw_%s_oli.dat'%(i.upper(),j.upper()),'w')
						oli_intx = self.post['oli_ind'][i][j][0]
						oli_inty = self.post['oli_ind'][i][j][1]
						oli_ints = self.post['oli_ind'][i][j][2]
						for k in range(len(oli_intx)): f.write('%9.6f\t%9.6f\t%9.6f\n'%(oli_intx[k],oli_inty[k],oli_ints[k]))
						f.close()
				
			lend = len(datx)
			filename = 'PROFILES/%s_raw.dat'%i.upper()
			f = open(filename,'w')
			f.write('%i\n'%lend)
			if use_rho: 
				xx1 = self.fit_eq['rho_to_psi'](datx)
				xx2 = datx
			else:
				xx1 = datx
				xx2 = self.fit_eq['psi_to_rho'](datx)

			for j in range(lend):	f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(xx1[j],daty[j],dats[j],datR[j],xx2[j]))

			f.close()
		return

	def write_avg_dat(self):
	
		for flag in self.prof_list:
			for k in self.__dict__['%s_list'%flag]:

				filename = 'PROFILES/%s_avg_%s.dat'%(flag.upper(),k.upper())
				nofile = self.fit_opt['file'][flag][k] == None
				use_rho = self.fit_opt['use_rho'][flag]
				if nofile: continue 

				f = open(filename,'w')

				x_type  = 'Psi_Norm'
				x_type2 = 'Rho_Norm'
				datx2 = self.fit_eq['psi_to_rho'](self.__dict__['%s_prof'%flag][k]['xxa'])
				if use_rho:
					x_type  = 'Rho_Norm'		
					x_type2 = 'Psi_Norm'
					datx2 = self.fit_eq['rho_to_psi'](self.__dict__['%s_prof'%flag][k]['xxa'])
		
				if (flag == 'ne'):
					f.write(x_type+'  '+x_type2+', NE[1E18/m3], STD \n')
				elif (flag == 'te'):
					f.write(x_type+'  '+x_type2+', TE[keV], 	STD \n')
				elif (flag == 'ti'):
					f.write(x_type+'  '+x_type2+', TI[kev],	    STD \n')
				elif (flag == 'vt'):
					f.write(x_type+'  '+x_type2+', VT[km/s], 	STD \n')

				datx = self.__dict__['%s_prof'%flag][k]['xxa']
				daty = self.__dict__['%s_prof'%flag][k]['avg']
				dats = self.__dict__['%s_prof'%flag][k]['avgs']

				for i in range(len(datx)):
					f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(datx[i],datx2[i],daty[i],dats[i]))
				f.close()

		return

	def write_fit_coefs(self):

		with open('PROFILES/fit_coefs.dat','w') as f:
			for flag in self.prof_list:
				fflag = self.func_list[self.fit_opt['func_type'][flag]-1];
				f.write('%s FUNC_TYPE: %s\n'%(flag.upper(),fflag.upper()))
				f.write('COEFS[VAL,ERR]\n');
				for i in range(len(self.post['popt'][flag])):
					f.write('a%02i: %13.7e %13.7e\n'%(i+1,self.post['popt'][flag][i],self.post['popte'][flag][i]))
		return

	def make_psiRZ_extended(self,psi,Z,islfs=True):

		if islfs:
			ind = np.where(self.fit_eq['R'] > self.fit_eq['rmag'])
			Rt = self.fit_eq['R'][ind]
			R = np.linspace(self.fit_eq['rmag'],max(Rt),301)
		else:
			ind = np.where(self.fit_eq['R'] < self.fit_eq['rmag'])
			Rt = self.fit_eq['R'][ind]
			R = np.linspace(self.fit_eq['rmag'],min(Rt),301)

		psi_map = self.fit_eq['psif'](R,Z)
		psi_map[0] = 0.0

		psif = interp1d(psi_map,R,'cubic')

		return psif(psi)

	def make_fitp(self,flag):

		self.__dict__['%s_prof'%flag]['fit2p'] = np.copy(self.__dict__['%s_prof'%flag]['fit2'])
		if self.fit_opt['use_rho'][flag]:
			self.__dict__['%s_prof'%flag]['fit2p'] = np.copy(self.__dict__['%s_prof'%flag]['fit2'])
			pf = interp1d(self.fit_eq['psin2'],self.__dict__['%s_prof'%flag]['fit2'])
			xx = self.fit_eq['psi_to_rho'](self.fit_eq['psin2'])
			for i in range(len(xx)):
				if xx[i] <= self.fit_eq['psin2'][-1]: val = pf(xx[i])
				self.__dict__['%s_prof'%flag]['fit2p'][i] = val

		return

	def adjust_sepval(self):

		min_pn = self.fit_opt['ashift']['min']; 
		max_pn = self.fit_opt['ashift']['max'];
		if not (self.post['isdat']['ti']): return
		if not self.fit_opt['ashift']['is']: return

		xp = self.ti_prof['ces']['xxa']
		yy = self.ti_prof['ces']['avg']
		
		max_pn2 = min(max_pn,max(xp))

		xr   = self.make_psiRZ_extended(xp,0.)
		sepr = self.make_psiRZ_extended(1.,0.)
		rp   = np.linspace(min_pn,max_pn2,401)

		yf   = interp1d(xp,yy)
		yt   = yf(rp) - 0.003
		ind  = np.argmin(yt)
		sepp = rp[ind]
		sepr2= self.make_psiRZ_extended(sepp,0.)

		delx = round(sepr - sepr2,4) + self.fit_opt['shift']['ti']['ces']
		self.print('>>> Auto-shift CES channel for %5.2f cm'%float(100.*delx))
		self.fit_opt['shift']['ti']['ces'] = delx
		self.fit_opt['shift']['vt']['ces'] = delx
		self.read_kinprof_kprofile()
		self.read_raw_kinprof()		
		return

	def adjust_ts_scale(self):

		if not self.fit_opt['ascale']: return
		self.print('>>> Start TS-NE auto scale')

		noprint_old = self.noprint
		cmulti_old = self.fit_opt['scale']['ne']['ts']['core']
		emulit_old = self.fit_opt['scale']['ne']['ts']['edge']
		self.noprint = True
		cn = self.fit_opt['scale']['ne']['ts']['cn']
		en = self.fit_opt['scale']['ne']['ts']['en']
		maxc = self.fit_opt['scale']['ne']['ts']['maxc']
		minc = self.fit_opt['scale']['ne']['ts']['minc']
		maxe = self.fit_opt['scale']['ne']['ts']['maxe']
		mine = self.fit_opt['scale']['ne']['ts']['mine']
		if self.fit_opt['ascale1d']: en = 1;

		if not self.fit_opt['ascale1d']: self.print('>>> #CN %i [%4.2f-%4.2f] #EN %i [%4.2f-%4.2f]'%(cn,minc,maxc,en,mine,maxe))
		else: self.print('>>> #CN %i [%4.2f-%4.2f] CORE/EDGE = %5.3f'%(cn,minc,maxc,ratio))

		xis = np.zeros((cn,en))
		core_multi = np.linspace(minc,maxc,cn)
		edge_multi = np.linspace(mine,maxe,en)

		cmulti = np.linspace(minc,maxc,41)
		emulti = np.linspace(mine,maxe,51)

		for i in range(cn):
			self.fit_opt['scale']['ne']['ts']['core'] = core_multi[i]
			for j in range(en):
				if self.fit_opt['ascale1d'] : self.fit_opt['scale']['ne']['ts']['edge'] = core_multi[i] / ratio
				else: self.fit_opt['scale']['ne']['ts']['edge'] = edge_multi[j]

				self.read_kinprof_kprofile()
				self.read_raw_kinprof()
				self.lmfit_fit('ne');
				xis[i,j] = self.post['chi']['ne']
				update_progress(float((i*en+j+1)/cn/en))	

		if not self.fit_opt['ascale1d']:		
			xf = interp2d(edge_multi,core_multi,xis)
			xif= xf(emulti,cmulti)
			ind = np.argmin(xif)
			ind1 = int((ind+1)/51)
			ind2 = int(ind - 51*ind1)
			if ind1 == 41: ind1 = 40; ind2 = 50;

			self.fit_opt['scale']['ne']['ts']['core'] = round(cmulti[ind1],3)
			self.fit_opt['scale']['ne']['ts']['edge'] = round(emulti[ind2],3)

		else:
			xf = interp1d(core_multi,xis[:,0])
			xif= xf(cmulti)
			ind = np.argmin(xif)

			self.fit_opt['scale']['ne']['ts']['core'] = round(cmulti[ind],3)
			self.fit_opt['scale']['ne']['ts']['edge'] = round(cmulti[ind]/ratio,3)			

		self.print('-> Auto Scaling factor, core-> %5.3f, edge-> %.3f '%(self.fit_opt['scale']['ne']['ts']['core'],self.fit_opt['scale']['ne']['ts']['edge']))
		self.read_kinprof_kprofile()
		self.read_raw_kinprof()		

		self.noprint = False
		return

	def mapping_variables(self):

		self.fit_eq['eq'] = eqdsk.eqdsk(self.fit_opt['file']['gfile'],False)
		self.fit_eq['eq'].read_eqdsk_file()
		self.fit_eq['eq'].contour = None

		self.fit_eq['rho_map'] = np.copy(self.fit_eq['eq'].prhoR[:,1])
		self.fit_eq['psi_map'] = np.copy(self.fit_eq['eq'].prhoR[:,0])
		self.fit_eq['R'] = np.copy(self.fit_eq['eq'].R)
		self.fit_eq['Z'] = np.copy(self.fit_eq['eq'].Z)
		self.fit_eq['rmag'] = self.fit_eq['eq'].rmag
		self.fit_eq['zmag'] = self.fit_eq['eq'].zmag

		self.fit_eq['psif'] = interp2d(self.fit_eq['R'],self.fit_eq['Z'],(self.fit_eq['eq'].psirz-self.fit_eq['eq'].smag)/(self.fit_eq['eq'].sbdy-self.fit_eq['eq'].smag))


		psin  = self.fit_eq['psif'](self.fit_eq['R'],self.device['inter1Z'])						
		psin2 = self.fit_eq['psif'](self.device['inter2R'],self.fit_eq['Z'])			
		psin2 = psin2.reshape(len(psin2))
		rmin = min(self.fit_eq['eq'].rzbdy[:,0])
		rmax = max(self.fit_eq['eq'].rzbdy[:,0])
		zmin = min(self.fit_eq['eq'].rzbdy[:,1])
		zmax = max(self.fit_eq['eq'].rzbdy[:,1])

		if (self.device['minwR'] > rmin):	self.device['minwR'] = rmin-0.1
		if (self.device['maxwR'] < rmax):	self.device['maxwR'] = rmax+0.1
		if (self.device['minwZ'] > zmin):	self.device['minwZ'] = zmin-0.2
		if (self.device['maxwZ'] < zmax):	self.device['maxwZ'] = zmax+0.2

		ind = np.where(self.fit_eq['R'] > (self.fit_eq['rmag']+0.01))
		R = self.fit_eq['R'][ind]
		psin0 = psin[ind]
		prf = interp1d(psin0,R,'cubic')
		self.Rsep = prf(1.0)
		try:	Rout0 = prf(1.08)
		except:	Rout0 = max(self.fit_eq['R'])

		ind = np.where(self.fit_eq['R'] < (self.fit_eq['rmag']-0.01))
		R = self.fit_eq['R'][ind]
		psin0 = psin[ind]
		prf = interp1d(psin0,R,'cubic')
	
		if not self.noprint: self.print('>>> R at Sep. midplane %3.2f [m]'%self.Rsep)	
	
		try:	Rin0 = prf(1.08)
		except:	Rin0 = min(self.fit_eq['R'])

		ind = np.where(self.fit_eq['Z'] > (self.fit_eq['zmag']+0.01))
		Z = self.fit_eq['Z'][ind]
		psin0 = psin2[ind]
		prf = interp1d(psin0,Z,'cubic')
		try:	Zout0 = prf(1.08)
		except:	Zout0 = max(self.fit_eq['Z'])

		ind = np.where(self.fit_eq['Z'] < (self.fit_eq['zmag']-0.01))
		Z = self.fit_eq['Z'][ind]
		psin0 = psin2[ind]
		prf = interp1d(psin0,Z,'cubic')		
		try:	Zin0 = prf(1.08)
		except:	Zin0 = min(self.fit_eq['Z'])			

		self.device['Rin']  = max(self.device['minwR'],Rin0)
		self.device['Rout'] = min(self.device['maxwR'],Rout0)
		self.device['Zin']  = max(self.device['minwZ'],Zin0)
		self.device['Zout'] = min(self.device['maxwZ'],Zout0)

		self.device['intr1']    = np.linspace(self.device['Rin'],self.device['Rout'],51)
		self.device['intz2']    = np.linspace(self.device['Zin'],self.device['Zout'],51)
		self.device['intp1']    = self.fit_eq['psif'](self.device['intr1'],self.device['inter1Z'])
		self.device['intp2']    = self.fit_eq['psif'](self.device['inter2R'],self.device['intz2'])
		self.device['intp2']    = self.device['intp2'].reshape(len(self.device['intp2']))
	
		for j in range(5):
			Rm = self.device['tci_Rmid'][self.inter_list[j+2]]
			Re = self.device['tci_Rend'][self.inter_list[j+2]]
#			Re = self.Rsep + 0.001
			L_traj = np.sqrt(Re**2 - Rm**2)
			R1 = np.linspace(0,L_traj,50)
			R2 = np.zeros(50)
			for i in range(50): R2[i] = np.sqrt(Rm**2+ R1[i]**2)
			self.device['tcip%i'%(j+1)] = self.fit_eq['psif'](R2,0.)
			self.device['tcir%i'%(j+1)] = np.copy(R1)
			
		psin = self.fit_eq['psif'](self.fit_eq['R'],self.fit_eq['zmag'])
		ind = np.where(self.fit_eq['R'] > (self.fit_eq['rmag']+0.02))
		R = self.fit_eq['R'][ind]
		psin = psin[ind]
		ind = np.where(R < (rmax+0.1))
		R = R[ind]
		psin = psin[ind]
		prf = interp1d(psin,R,'cubic')
		RO = np.linspace(self.fit_eq['rmag'],prf(1.0),51)
		psin = self.fit_eq['psif'](RO,self.fit_eq['zmag'])
		psin[0] = 0.
		psin[-1] = 1.
		self.fit_eq['Rf'] = interp1d(psin,RO,'cubic')
		self.fit_eq['eq'].construct_volume()
		self.fit_eq['Vf'] = interp1d(self.fit_eq['eq'].avolp[:,0],self.fit_eq['eq'].avolp[:,2],'cubic')

		psi_map_extend = np.zeros(len(self.fit_eq['psi_map'])+100)
		rho_map_extend = np.copy(psi_map_extend)

		self.fit_eq['psi_to_rho'] = interp1d(self.fit_eq['psi_map'],self.fit_eq['rho_map'],'cubic')

		for i in range(len(self.fit_eq['psi_map'])):
			psi_map_extend[i] = self.fit_eq['psi_map'][i]
			rho_map_extend[i] = self.fit_eq['rho_map'][i]
		drdp = (1.0-self.fit_eq['psi_to_rho'](0.9999))/(1.0-0.9999)

		for i in range(100):
			psi_map_extend[len(self.fit_eq['psi_map'])+i] = 1.0 + 0.6/60.*float(i+1)
			rho_map_extend[len(self.fit_eq['psi_map'])+i] = 1.0 + drdp * 0.6/60.*float(i+1)
		
		self.fit_eq['rho_to_psi'] = interp1d(rho_map_extend,psi_map_extend,'cubic')
		self.fit_eq['psi_to_rho'] = interp1d(psi_map_extend,rho_map_extend,'cubic')
		

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
	
			for flag in self.prof_list:
				for k in self.__dict__['%s_list'%flag]:
					self.fit_opt['psi_end'][flag][k] = self.read_namelist_str(line,'PSI_END_%s'%flag.upper(),self.fit_opt['psi_end'][flag][k],2)

			self.fit_opt['file']['te']['ts']  = self.read_namelist_str(line,'TE_DAT_FILE',self.fit_opt['file']['te']['ts'],3)
			self.fit_opt['file']['ne']['ts']  = self.read_namelist_str(line,'NE_DAT_FILE',self.fit_opt['file']['ne']['ts'],3)
			self.fit_opt['file']['ti']['ces'] = self.read_namelist_str(line,'TI_DAT_FILE',self.fit_opt['file']['ti']['ces'],3)
			self.fit_opt['file']['vt']['ces'] = self.read_namelist_str(line,'VT_DAT_FILE',self.fit_opt['file']['vt']['ces'],3)

			self.fit_opt['file']['te']['kfile'] = self.read_namelist_str(line,'TE_FILE',self.fit_opt['file']['te']['kfile'],3)
			self.fit_opt['file']['ne']['kfile'] = self.read_namelist_str(line,'NE_FILE',self.fit_opt['file']['ne']['kfile'],3)
			self.fit_opt['file']['ti']['kfile'] = self.read_namelist_str(line,'TI_FILE',self.fit_opt['file']['ti']['kfile'],3)
			self.fit_opt['file']['vt']['kfile'] = self.read_namelist_str(line,'VT_FILE',self.fit_opt['file']['vt']['kfile'],3)
			
			self.fit_opt['target_density'] = self.read_namelist_str(line,'DENSITY_SCALE',self.fit_opt['target_density'],2)
			
			self.fit_opt['use_ti_eped']  = self.read_namelist_str(line,'USE_TI_EPED',self.fit_opt['use_ti_eped'],4)
			self.fit_opt['use_ti_width'] = self.read_namelist_str(line,'USE_TI_WIDTH',self.fit_opt['use_ti_width'],4)
		
			self.fit_opt['sep_fix']['te'] = self.read_namelist_str(line,'FIX_TE_sep',self.fit_opt['sep_fix']['te'],4)
			self.fit_opt['sep_fix']['ne'] = self.read_namelist_str(line,'FIX_NE_sep',self.fit_opt['sep_fix']['ne'],4)
			self.fit_opt['sep_fix']['ti'] = self.read_namelist_str(line,'FIX_TI_sep',self.fit_opt['sep_fix']['ti'],4)
			self.fit_opt['sep_fix']['vt'] = self.read_namelist_str(line,'FIX_VT_sep',self.fit_opt['sep_fix']['vt'],4)

			self.fit_opt['sep_val']['te'] = self.read_namelist_str(line,'TE_sep',self.fit_opt['sep_val']['te'],2)
			self.fit_opt['sep_val']['ne'] = self.read_namelist_str(line,'NE_sep',self.fit_opt['sep_val']['ne'],2)
			self.fit_opt['sep_val']['ti'] = self.read_namelist_str(line,'TI_sep',self.fit_opt['sep_val']['ti'],2)
			self.fit_opt['sep_val']['vt'] = self.read_namelist_str(line,'VT_sep',self.fit_opt['sep_val']['vt'],2)
			
			self.fit_opt['ped_scan_fit'] = self.read_namelist_str(line,'Ped_scan_fit',self.fit_opt['ped_scan_fit'],4)

			self.fit_opt['width_fix']['te'] = self.read_namelist_str(line,'FIX_TE_width',self.fit_opt['width_fix']['te'],4)
			self.fit_opt['width_fix']['ne'] = self.read_namelist_str(line,'FIX_NE_width',self.fit_opt['width_fix']['ne'],4)
			self.fit_opt['width_fix']['ti'] = self.read_namelist_str(line,'FIX_TI_width',self.fit_opt['width_fix']['ti'],4)
			self.fit_opt['width_fix']['vt'] = self.read_namelist_str(line,'FIX_VT_width',self.fit_opt['width_fix']['vt'],4)

			self.fit_opt['width_val']['te'] = self.read_namelist_str(line,'TE_width',self.fit_opt['width_val']['te'],2)
			self.fit_opt['width_val']['ne'] = self.read_namelist_str(line,'NE_width',self.fit_opt['width_val']['ne'],2)
			self.fit_opt['width_val']['ti'] = self.read_namelist_str(line,'TI_width',self.fit_opt['width_val']['ti'],2)
			self.fit_opt['width_val']['vt'] = self.read_namelist_str(line,'VT_width',self.fit_opt['width_val']['vt'],2)

			self.fit_opt['avg']['te']['ts']  = self.read_namelist_str(line,'TE_AVG_DAT',self.fit_opt['avg']['te']['ts'],4)
			self.fit_opt['avg']['ne']['ts']  = self.read_namelist_str(line,'NE_AVG_DAT',self.fit_opt['avg']['ne']['ts'],4)
			self.fit_opt['avg']['ti']['ces'] = self.read_namelist_str(line,'TI_AVG_DAT',self.fit_opt['avg']['ti']['ces'],4)
			self.fit_opt['avg']['vt']['ces'] = self.read_namelist_str(line,'VT_AVG_DAT',self.fit_opt['avg']['vt']['ces'],4)

			self.fit_opt['exclude']['te']['ts']  = self.read_namelist_str(line,'TE_EXC',self.fit_opt['exclude']['te']['ts'],3)
			self.fit_opt['exclude']['ne']['ts']  = self.read_namelist_str(line,'NE_EXC',self.fit_opt['exclude']['ne']['ts'],3)
			self.fit_opt['exclude']['ti']['ces'] = self.read_namelist_str(line,'TI_EXC',self.fit_opt['exclude']['ti']['ces'],3)
			self.fit_opt['exclude']['vt']['ces'] = self.read_namelist_str(line,'VT_EXC',self.fit_opt['exclude']['vt']['ces'],3)
		
			self.fit_opt['raw_fit']['te'] = self.read_namelist_str(line,'TE_RAW_FIT',self.fit_opt['raw_fit']['te'],4)
			self.fit_opt['raw_fit']['ne'] = self.read_namelist_str(line,'NE_RAW_FIT',self.fit_opt['raw_fit']['ne'],4)
			self.fit_opt['raw_fit']['ti'] = self.read_namelist_str(line,'TI_RAW_FIT',self.fit_opt['raw_fit']['ti'],4)
			self.fit_opt['raw_fit']['vt'] = self.read_namelist_str(line,'VT_RAW_FIT',self.fit_opt['raw_fit']['vt'],4)
			
			self.fit_opt['oli']['te']['ts']['use']  = self.read_namelist_str(line,'TE_USE_OUT',self.fit_opt['oli']['te']['ts']['use'],4)
			self.fit_opt['oli']['ne']['ts']['use']  = self.read_namelist_str(line,'NE_USE_OUT',self.fit_opt['oli']['ne']['ts']['use'],4)
			self.fit_opt['oli']['ti']['ces']['use'] = self.read_namelist_str(line,'TI_USE_OUT',self.fit_opt['oli']['ti']['ces']['use'],4)
			self.fit_opt['oli']['vt']['ces']['use'] = self.read_namelist_str(line,'VT_USE_OUT',self.fit_opt['oli']['vt']['ces']['use'],4)

			self.fit_opt['oli']['te']['ts']['per']  = self.read_namelist_str(line,'TE_OLI_CUT',self.fit_opt['oli']['te']['ts']['per'],2)
			self.fit_opt['oli']['ne']['ts']['per']  = self.read_namelist_str(line,'NE_OLI_CUT',self.fit_opt['oli']['ne']['ts']['per'],2)
			self.fit_opt['oli']['ti']['ces']['per'] = self.read_namelist_str(line,'TI_OLI_CUT',self.fit_opt['oli']['ti']['ces']['per'],2)
			self.fit_opt['oli']['vt']['ces']['per'] = self.read_namelist_str(line,'VT_OLI_CUT',self.fit_opt['oli']['vt']['ces']['per'],2)
		
			self.fit_opt['oli']['te']['n'] = self.read_namelist_str(line,'TE_OLI_CUTN',self.fit_opt['oli']['te']['n'],1)
			self.fit_opt['oli']['ne']['n'] = self.read_namelist_str(line,'NE_OLI_CUTN',self.fit_opt['oli']['ne']['n'],1)
			self.fit_opt['oli']['ti']['n'] = self.read_namelist_str(line,'TI_OLI_CUTN',self.fit_opt['oli']['ti']['n'],1)
			self.fit_opt['oli']['vt']['n'] = self.read_namelist_str(line,'VT_OLI_CUTN',self.fit_opt['oli']['vt']['n'],1)

			self.fit_opt['func_type']['te'] = self.read_namelist_str(line,'FIT_TYPE_TE',self.fit_opt['func_type']['te'],1)
			self.fit_opt['func_type']['ne'] = self.read_namelist_str(line,'FIT_TYPE_NE',self.fit_opt['func_type']['ne'],1)
			self.fit_opt['func_type']['ti'] = self.read_namelist_str(line,'FIT_TYPE_TI',self.fit_opt['func_type']['ti'],1)
			self.fit_opt['func_type']['vt'] = self.read_namelist_str(line,'FIT_TYPE_VT',self.fit_opt['func_type']['vt'],1)

			self.fit_opt['sspline_order']['te'] = self.read_namelist_str(line,'S_TE',self.fit_opt['sspline_order']['te'],2)
			self.fit_opt['sspline_order']['ne'] = self.read_namelist_str(line,'S_NE',self.fit_opt['sspline_order']['ne'],2)
			self.fit_opt['sspline_order']['ti'] = self.read_namelist_str(line,'S_TI',self.fit_opt['sspline_order']['ti'],2)
			self.fit_opt['sspline_order']['vt'] = self.read_namelist_str(line,'S_VT',self.fit_opt['sspline_order']['vt'],2)
	
			self.fit_opt['std']['te']['ts']['use']  = self.read_namelist_str(line,'STD_TE',self.fit_opt['std']['te']['ts']['use'],4)
			self.fit_opt['std']['ne']['ts']['use']  = self.read_namelist_str(line,'STD_NE',self.fit_opt['std']['ne']['ts']['use'],4)
			self.fit_opt['std']['ti']['ces']['use'] = self.read_namelist_str(line,'STD_TI',self.fit_opt['std']['ti']['ces']['use'],4)
			self.fit_opt['std']['vt']['ces']['use'] = self.read_namelist_str(line,'STD_VT',self.fit_opt['std']['vt']['ces']['use'],4)
			
			self.fit_opt['std']['te']['ts']['raw']  = self.read_namelist_str(line,'RAW_STD_TE',self.fit_opt['std']['te']['ts']['raw'],4)
			self.fit_opt['std']['ne']['ts']['raw']  = self.read_namelist_str(line,'RAW_STD_NE',self.fit_opt['std']['ne']['ts']['raw'],4)
			self.fit_opt['std']['ti']['ces']['raw'] = self.read_namelist_str(line,'RAW_STD_TI',self.fit_opt['std']['ti']['ces']['raw'],4)
			self.fit_opt['std']['vt']['ces']['raw'] = self.read_namelist_str(line,'RAW_STD_VT',self.fit_opt['std']['vt']['ces']['raw'],4)
		
			self.fit_opt['amain'] = self.read_namelist_str(line,'AMAIN',self.fit_opt['amain'],2)
			self.fit_opt['zeff']  = self.read_namelist_str(line,'ZEFF', self.fit_opt['zeff'],2)
			self.fit_opt['aimp']  = self.read_namelist_str(line,'AIMP', self.fit_opt['aimp'],2)
			self.fit_opt['zimp']  = self.read_namelist_str(line,'ZIMP', self.fit_opt['zimp'],2)

			self.fit_opt['file']['gfile'] = self.read_namelist_str(line,'EQDSK',self.fit_opt['file']['gfile'],3)
			self.fit_opt['use_rho']['te'] = self.read_namelist_str(line,'USE_RHO',self.fit_opt['use_rho']['te'],4)
			self.fit_opt['use_rho']['ne'] = self.read_namelist_str(line,'USE_RHO',self.fit_opt['use_rho']['ne'],4)
			self.fit_opt['use_rho']['ti'] = self.read_namelist_str(line,'USE_RHO',self.fit_opt['use_rho']['ti'],4)
			self.fit_opt['use_rho']['vt'] = self.read_namelist_str(line,'USE_RHO',self.fit_opt['use_rho']['vt'],4)

			self.fit_opt['shift']['te']['ts']  = self.read_namelist_str(line,'SHIFT_TE',self.fit_opt['shift']['te']['ts'],2)
			self.fit_opt['shift']['ne']['ts']  = self.read_namelist_str(line,'SHIFT_NE',self.fit_opt['shift']['ne']['ts'],2)
			self.fit_opt['shift']['ti']['ces'] = self.read_namelist_str(line,'SHIFT_TI',self.fit_opt['shift']['ti']['ces'],2)
			self.fit_opt['shift']['vt']['ces'] = self.read_namelist_str(line,'SHIFT_VT',self.fit_opt['shift']['vt']['ces'],2)
	
			self.fit_opt['int01']['val'] = self.read_namelist_str(line,'INT_N1',self.fit_opt['int01']['val'],2)
			self.fit_opt['int02']['val'] = self.read_namelist_str(line,'INT_N2',self.fit_opt['int02']['val'],2)
			self.fit_opt['int01']['sig'] = self.read_namelist_str(line,'INT_S1',self.fit_opt['int01']['sig'],2)
			self.fit_opt['int02']['sig'] = self.read_namelist_str(line,'INT_S2',self.fit_opt['int02']['sig'],2)			

			self.fit_opt['tci01']['val'] = self.read_namelist_str(line,'TCI_N1',self.fit_opt['tci01']['val'],2)
			self.fit_opt['tci02']['val'] = self.read_namelist_str(line,'TCI_N2',self.fit_opt['tci02']['val'],2)
			self.fit_opt['tci03']['val'] = self.read_namelist_str(line,'TCI_N3',self.fit_opt['tci03']['val'],2)
			self.fit_opt['tci04']['val'] = self.read_namelist_str(line,'TCI_N4',self.fit_opt['tci04']['val'],2)
			self.fit_opt['tci05']['val'] = self.read_namelist_str(line,'TCI_N5',self.fit_opt['tci05']['val'],2)

			self.fit_opt['tci01']['sig'] = self.read_namelist_str(line,'TCI_S1',self.fit_opt['tci01']['sig'],2)
			self.fit_opt['tci02']['sig'] = self.read_namelist_str(line,'TCI_S2',self.fit_opt['tci02']['sig'],2)
			self.fit_opt['tci03']['sig'] = self.read_namelist_str(line,'TCI_S3',self.fit_opt['tci03']['sig'],2)
			self.fit_opt['tci04']['sig'] = self.read_namelist_str(line,'TCI_S4',self.fit_opt['tci04']['sig'],2)
			self.fit_opt['tci05']['sig'] = self.read_namelist_str(line,'TCI_S5',self.fit_opt['tci05']['sig'],2)			

			self.fit_opt['line_type'] = self.read_namelist_str(line,'LINE_TYPE',self.fit_opt['line_type'],1)

		f4.close()
		self.input_file_check()
		self.adjust_variables()
		return

	def input_file_check(self):

		for flag in self.prof_list:
			for k in self.__dict__['%s_list'%flag]:
				filename = self.fit_opt['file'][flag][k]
				if not filename == None:
					if not os.path.isfile(filename): self.fit_opt['file'][flag][k] = None				

			k = 'kfile'
			filename = self.fit_opt['file'][flag][k]
			if not filename == None:
				if not os.path.isfile(filename): self.fit_opt['file'][flag][k] = None
		return
		
	def adjust_variables(self):

		if (self.fit_opt['target_density'] > 0.0):
				self.fit_opt['use_density_scale'] = True
		else: self.fit_opt['use_density_scale'] = False

		if (self.fit_opt['ped_scan_fit']):
			self.fit_opt['func_type']['te'] = 4
			self.fit_opt['func_type']['ne'] = 4

		self.post['isdat'] = dict()
		self.post['iskfile'] = dict()
		self.post['isrfile'] = dict()

		for flag in self.prof_list:	
			if (self.fit_opt['func_type'][flag] == 5 or self.fit_opt['func_type'][flag] == 7 ):
				for k in self.__dict__['%s_list'%flag]:
					self.fit_opt['avg'][flag][k] = True

			for k in self.__dict__['%s_list'%flag]:
				if self.fit_opt['exclude'][flag][k] == None: self.fit_opt['exclude'][flag][k] =''
				len1 = len(self.fit_opt['exclude'][flag][k].split())
				if len1 == 0: self.fit_opt['exclude'][flag][k] = ''
				else: self.fit_opt['exclude'][flag][k] = self.fit_opt['exclude'][flag][k].split()[0]

			self.post['iskfile'][flag] = False
			if not self.fit_opt['file'][flag]['kfile'] == None:	self.post['iskfile'][flag] = True
			
			self.post['isrfile'][flag] = False
			for k in self.__dict__['%s_list'%flag]:
				if not self.fit_opt['file'][flag][k] == None: self.post['isrfile'][flag] = True

			if flag == 'ne':
				self.ne_tcifit = False
				tci_count = 0
				for i in self.inter_list:
					if (self.fit_opt[i]['val']>.0 and self.fit_opt[i]['sig'] > 0.):
						tci_count = tci_count + 1		
	
				if (not self.post['isrfile'][flag] and (tci_count >= self.ne_tcifit_ch_nmin)):
					self.fit_opt['file'][flag]['ts'] = dummy_dir
					self.post['isrfile'][flag] = True
					self.print('No TS and REFL profiles and use #%g INTER CH.'%tci_count)
					self.ne_tcifit = True

			self.post['isdat'][flag]= True
			if self.fit_opt['raw_fit'][flag]:	
				if not self.post['isrfile'][flag]:self.post['isdat'][flag]= False
			else:
				if not self.post['iskfile'][flag]:self.post['isdat'][flag]= False

		if self.fit_opt['use_ti_eped']:	self.fit_opt['func_type']['ti'] = 4

		if self.fit_opt['use_ti_width']:
			self.fit_opt['width_fix']['te'] = True
			self.fit_opt['width_fix']['ne'] = True

		if self.fit_opt['file']['gfile'] == None:
			self.print('>>> No eq file!')
			exit()
			
		if not self.fit_opt['file']['gfile'] == self.post['eq_old']:	self.post['eq_change'] = True

		return	

	def declare_variable(self):
		# -- variable list
		self.prof_list  = ['te','ne','ti','vt']
		self.inter_list = ['int01','int02','tci01','tci02','tci03','tci04','tci05']
		self.func_list  = ['core','mtanh','ptanh','eped','spline','eped2','nspline','spline2','eped3']
		self.func_varn  = [     4,      8,      7,     6,       0,      7,        0,        0,      6]
		self.func       = dict()
		self.ne_list    = ['ts','tse','refl']
		self.te_list    = ['ts','tse','ece']
		self.ti_list    = ['ces']
		self.vt_list    = ['ces']
		self.factor     = dict()
		self.flag 	= dict()

		self.fit_opt = dict()
		self.fit_eq  = dict()
		self.post    = dict()
		self.device  = dict()

		self.plot_shape1 = ['x','>','D']
		self.plot_shape2 = ['o','s','p']

		self.first_run = True

		return
			
	def initialise_variables(self):
		
		self.declare_variable()
		self.initialise_fitopt()
		self.force_fit = False
		self.noprint = False
		self.forcefit = False
		self.forcefit2 = False
		self.hmodefit = False
		self.ne_tcifit= False
		self.exp_lsq  = False

		self.ne_tcifit_ch_nmin = 3

		self.fit_trunc = 12
		self.line_trunc = 12

		self.factor['te'] = 1.e-3
		self.factor['ne'] = 1.e-1
		self.factor['ti'] = 1.e-3
		self.factor['vt'] = 1.e+0		

		self.func['list'] = self.func_list
		self.func['varn'] = dict()
		self.func['num'] = dict()
		for i in range(len(self.func_list)):
			flag = self.func_list[i]
			varn = self.func_varn[i]
			self.func['varn'][flag.lower()] = varn
			self.func['num'][flag.lower()] = i+1

		self.flag['te'] = 'TE[eV]'
		self.flag['ne'] = 'NE[1E18/m3]'
		self.flag['ti'] = 'TI[eV]'
		self.flag['vt'] = 'VT[km/s]'		

		for i in self.prof_list: self.__dict__['%s_prof'%i] = dict()	

		self.fit_eq['psin1'] = np.linspace(0,1.0,401)
		self.fit_eq['psin2'] = np.linspace(0,1.3,401)

		for i in self.prof_list:
			self.__dict__['%s_prof'%i]['fit1']     = np.zeros(401)	
			self.__dict__['%s_prof'%i]['fit2']     = np.zeros(401)
			self.__dict__['%s_prof'%i]['fit2p']    = np.zeros(401)
			self.__dict__['%s_prof'%i]['fit_old']  = np.zeros(401)
			for j in self.__dict__['%s_list'%i]:
				self.__dict__['%s_prof'%i][j] = dict()
				self.__dict__['%s_prof'%i][j]['file'] = None
				self.__dict__['%s_prof'%i][j]['opt']  = dict()
				for k in ['xxr','raw','raws','raws2','datR','datz','xxa','avg','avgs','avgs2']: self.__dict__['%s_prof'%i][j][k] = None
			j = 'kfile'
			self.__dict__['%s_prof'%i][j] = dict()
			self.__dict__['%s_prof'%i][j]['file'] = None
			self.__dict__['%s_prof'%i][j]['opt']  = dict()
			for k in ['xxr','raw','raws','raws2','datR','datz','xxa','avg','avgs','avgs2']: self.__dict__['%s_prof'%i][j][k] = None

		for i in self.inter_list:	self.post[i] = 0.0
		
		self.post['time']   = 0.
		self.post['scaled_density'] = 1.0
		self.post['eq_change']   = True
		self.post['eq_old']      = None
		self.post['pedmid'] 	 = dict()
		self.post['popt']   	 = dict()
		self.post['popte']  	 = dict()
		self.post['same_opt']    = dict()
		self.post['fit_ind']     = dict()
		self.post['oli_ind']     = dict()
		self.post['prof_dim']	 = dict()
		self.post['opt_change']  = dict()
		self.post['chi']	     = dict()
		self.post['err']	     = dict()		
		self.post['didfit']	     = dict()
		self.post['den_diff']    = dict()
		self.post['width']       = dict()
		self.post['func_history']= dict()
		self.post['spline_knots']= dict()

		for i in self.inter_list: self.post['den_diff'][i] = 0.0
		
		for i in self.prof_list:
			self.post['pedmid'][i]    = 0.97
			self.post['width'][i]     = 0.97
			self.post['popt'][i]      = np.ones(10)
			self.post['popte'][i]     = np.ones(10)
			self.post['same_opt'][i]  = False
			self.post['didfit'][i]    = False
			self.post['fit_ind'][i]   = dict()
			self.post['oli_ind'][i]   = dict()
			self.post['prof_dim'][i]  = dict()
			self.post['func_history'][i]=dict()

			self.post['opt_change'][i]= True
			for j in self.__dict__['%s_list'%i]:
				self.post['fit_ind'][i][j]  = []
				self.post['oli_ind'][i][j]  = []
				self.post['prof_dim'][i][j] = [1,1]

			for j in range(len(self.func_list)):
				self.post['func_history'][i][j+1] = False

		self.post['bpped'] = 0.
		self.post['bppede'] = 0.
		self.post['wmhd'] = 0.
		self.post['wkin'] = 0.
		self.post['width'] = np.zeros(4)
		self.post['widthe'] = np.zeros(4)
		self.post['shift']  = dict()
		self.post['shift']['ts']  = 0.
		self.post['shift']['ces'] = 0.
		self.post['shift']['ref'] = 0.
		self.post['shift']['ece'] = 0.

		#Device variables
		self.device['inter1Z']  = 0.0				
		self.device['inter2R']  = 1.8
		self.device['inter1L']  = 1.9
		self.device['inter2L']  = 2.75
		self.device['tci_Rmid'] = dict()
		self.device['tci_Rend'] = dict()
		self.device['tci_Rmid']['tci01'] = 1.343 
		self.device['tci_Rmid']['tci02'] = 1.784
		self.device['tci_Rmid']['tci03'] = 1.912
		self.device['tci_Rmid']['tci04'] = 2.039
		self.device['tci_Rmid']['tci05'] = 2.163
		self.device['tci_Rend']['tci01'] = 2.253
		self.device['tci_Rend']['tci02'] = 2.253
		self.device['tci_Rend']['tci03'] = 2.253
		self.device['tci_Rend']['tci04'] = 2.253
		self.device['tci_Rend']['tci05'] = 2.253
		self.device['ne_sol_decayL'] = 0.08
		self.device['tci_L'] = dict()
		self.device['tci_L']['tci01'] = 7.23
		self.device['tci_L']['tci02'] = 5.51
		self.device['tci_L']['tci03'] = 4.76
		self.device['tci_L']['tci04'] = 3.80
		self.device['tci_L']['tci05'] = 2.52

		self.device['maxwR'] =  2.31
		self.device['minwR'] =  1.27
		self.device['maxwZ'] =  1.0
		self.device['minwZ'] = -1.0

		self.device['var_crit'] = dict()
		self.device['var_crit']['te'] = 1.
		self.device['var_crit']['ti'] = 1.
		self.device['var_crit']['vt'] = 200
		self.device['var_crit']['ne'] = 4.

		return

	def initialise_fitopt(self):

		for i in self.inter_list:
			self.fit_opt[i] = dict()
			self.fit_opt[i]['val'] = 0.0
			self.fit_opt[i]['sig'] = 0.0

		self.fit_opt['use_ti_width']      = False
		self.fit_opt['use_ti_eped']       = False
		self.fit_opt['ped_scan_fit']      = False
		self.fit_opt['line_type']         = 0
		self.fit_opt['target_density']    = 0.
		self.fit_opt['use_density_scale'] = False
		self.fit_opt['ashift']            = dict()
		self.fit_opt['ashift']['is']	  = False
		self.fit_opt['ashift']['min']	  = 0.9
		self.fit_opt['ashift']['max']	  = 1.05
		self.fit_opt['ascale']            = False
		self.fit_opt['ascale1d']          = False

		self.fit_opt['dacrit']            = 1.5
		self.fit_opt['duty']              = 30.

		self.fit_opt['use_rho']  	 	  = dict()
		self.fit_opt['psi_end'] 		  = dict()
		self.fit_opt['sep_fix']  		  = dict()
		self.fit_opt['sep_val'] 		  = dict()
		self.fit_opt['width_fix']		  = dict()
		self.fit_opt['width_val']		  = dict()
		self.fit_opt['raw_fit']  		  = dict()
		self.fit_opt['func_type']		  = dict()
		self.fit_opt['sspline_order']	  = dict()

		self.fit_opt['file']			  = dict()
		self.fit_opt['avg']				  = dict()
		self.fit_opt['shift']			  = dict()
		self.fit_opt['exclude']			  = dict()
		self.fit_opt['oli']				  = dict()
		self.fit_opt['std']			      = dict()
		self.fit_opt['weight']			  = dict()
		
		self.fit_opt['file']['gfile']     = None

		for i in self.prof_list:	
			self.fit_opt['use_rho'][i]       = False
			self.fit_opt['psi_end'][i]       = dict()
			self.fit_opt['sep_fix'][i]       = False
			self.fit_opt['sep_val'][i]       = 0.
			self.fit_opt['width_fix'][i]     = False
			self.fit_opt['width_val'][i]     = 0.
			self.fit_opt['raw_fit'][i]       = True
			self.fit_opt['sspline_order'][i] = -1

			self.fit_opt['file'][i]		 = dict()
			self.fit_opt['avg'][i]		 = dict()
			self.fit_opt['shift'][i]	 = dict()
			self.fit_opt['exclude'][i]	 = dict()
			self.fit_opt['oli'][i]		 = dict()
			self.fit_opt['std'][i]		 = dict()
			self.fit_opt['weight'][i]	 = dict()

			self.fit_opt['oli'][i]['n']  = 2
			self.fit_opt['file'][i]['kfile'] = None

			for j in self.__dict__['%s_list'%i]:

				self.fit_opt['file'][i][j]		 = None
				self.fit_opt['avg'][i][j]		 = True
				self.fit_opt['shift'][i][j]	     = 0.
				self.fit_opt['exclude'][i][j]	 = ''
				self.fit_opt['oli'][i][j]		 = dict()
				self.fit_opt['std'][i][j]		 = dict()
				self.fit_opt['oli'][i][j]['use'] = False
				self.fit_opt['oli'][i][j]['per'] = 0.95
				self.fit_opt['std'][i][j]['use'] = True
				self.fit_opt['std'][i][j]['raw'] = True
				self.fit_opt['weight'][i][j]     = 1.
				self.fit_opt['psi_end'][i][j]    = 0.

		self.fit_opt['exclude']['te']['ece'] = '8,13,14,15,16,17'
		self.fit_opt['weight']['te']['ece'] = 0.3
		self.fit_opt['weight']['ne']['ts']  = 0.5
		self.fit_opt['weight']['ne']['tse'] = 0.5

		self.fit_opt['psi_end']['te']['ece'] = 0.55
		self.fit_opt['psi_end']['te']['ts']  = 1.0
		self.fit_opt['psi_end']['te']['tse'] = 1.0
		self.fit_opt['psi_end']['ne']['ts']  = 1.0
		self.fit_opt['psi_end']['ne']['tse'] = 1.0
		
		self.fit_opt['func_type']['te']   = 4
		self.fit_opt['func_type']['ne']   = 4
		self.fit_opt['func_type']['ti']   = 3
		self.fit_opt['func_type']['vt']   = 3

		self.fit_opt['amain'] = 2.0
		self.fit_opt['zeff']  = 2.0
		self.fit_opt['zimp']  = 6.0
		self.fit_opt['aimp']  = 12.0

		self.fit_opt['plot'] = dict()
		for flag in self.prof_list:
			self.fit_opt['plot'][flag] = dict()
			self.fit_opt['plot'][flag]['xmin'] = -0.03
			self.fit_opt['plot'][flag]['xmax'] = 1.25
			self.fit_opt['plot'][flag]['ymin'] = -1
			self.fit_opt['plot'][flag]['ymax'] = -1

		self.fit_opt['mds']	 = dict()
		self.fit_opt['mds']['shot'] = 0
		self.fit_opt['mds']['time'] = 0
		self.fit_opt['mds']['times'] = ''
		self.fit_opt['mds']['cmse'] = False
		self.fit_opt['mds']['dt'] = dict()
		self.fit_opt['mds']['dt']['ts']  = 0
		self.fit_opt['mds']['dt']['ces'] = 0
		self.fit_opt['mds']['dt']['tci'] = 0
		self.fit_opt['mds']['dt']['ece'] = 0
		self.fit_opt['mds']['dt']['ref'] = 0
		self.fit_opt['mds']['dt']['ref2'] = 0
		self.fit_opt['mds']['efit'] = 1

		self.fit_opt['scale'] = dict()
		for flag in self.prof_list:
			self.fit_opt['scale'][flag] = dict()
			for kk in self.__dict__['%s_list'%flag]:
				self.fit_opt['scale'][flag][kk] = dict()
				self.fit_opt['scale'][flag][kk]['core'] = 1.
				self.fit_opt['scale'][flag][kk]['edge'] = 1.
				self.fit_opt['scale'][flag][kk]['mine'] = 0.5
				self.fit_opt['scale'][flag][kk]['maxe'] = 2.
				self.fit_opt['scale'][flag][kk]['minc'] = 0.5
				self.fit_opt['scale'][flag][kk]['maxc'] = 2.
				self.fit_opt['scale'][flag][kk]['cn'] = 7
				self.fit_opt['scale'][flag][kk]['en'] = 7
				
				self.fit_opt['scale'][flag][kk]['ch_cal']  = False
				self.fit_opt['scale'][flag][kk]['core_ch'] = 14
				self.fit_opt['scale'][flag][kk]['edge_ch'] = 1
				self.fit_opt['scale'][flag][kk]['fix_core']= True
				self.fit_opt['scale'][flag][kk]['ratio']= 1.
		return

	def wash_exclude(self):

		for flag in ['te','ne','ti','vt']:
			for kk in self.__dict__['%s_list'%flag]:
				if not self.fit_opt['exclude'][flag][kk] == '':
					last_ind = int(self.fit_opt['exclude'][flag][kk].split(',')[-1])+1
					if last_ind > self.post['prof_dim'][flag][kk][0]:
						self.print('>>> Wash Ch. exlude for %s-%s'%(flag,kk))
						self.print('>>> %s'%self.fit_opt['exclude'][flag][kk])
						self.fit_opt['exclude'][flag][kk] = ''
		return

	def main_run(self):

		self.input_file_check()
		self.adjust_variables()

		try:	os.mkdir('PROFILES')
		except:	pass
		if not (self.post['eq_change']):
			if not self.post['eq_old'] == self.fit_opt['file']['gfile']:
				self.post['eq_change'] = True

		if (self.post['eq_change'] or self.first_run):
			self.post['eq_old'] = self.fit_opt['file']['gfile']
			self.mapping_variables()
			self.post['eq_change'] = False
			if not self.noprint:
				self.print('>>> New equilibrium is given -> Re-equilibrium mapping!')
				self.print('>>> '+self.fit_opt['file']['gfile'])

		self.read_kinprof_kprofile()
		self.read_raw_kinprof()
		self.wash_exclude()
		if not self.noprint:
			self.print('=----------------------------------------------------------------------------=')
			self.print('=                             FITTING STARTED                                =')
			self.print('=----------------------------------------------------------------------------=')

		self.adjust_sepval()

		if not self.first_run: self.check_options_change()
		else: 
			self.param_old = copy.deepcopy(self.param)
			for k in self.prof_list:        self.post['opt_change'][k] = True

		if (self.forcefit or self.forcefit2): self.forced_fit_params()

		if (self.post['isdat']['ti'] and self.post['opt_change']['ti']):
			if self.post['didfit']['ti']: self.ti_prof['fit_old'] = np.copy(self.ti_prof['fit2p'])
			self.lmfit_fit('ti'); self.post['didfit']['ti'] = True;

			if (self.fit_opt['ped_scan_fit'] and self.fit_opt['use_ti_eped'] and self.fit_opt['use_rho']['ti']):
				temp_coef = self.modify_rho_to_psi_eped(self.post['popt']['ti'],'ti',True)
				temp_psi  = self.fit_eq['rho_to_psi'](self.fit_eq['psin2'])
				self.ti_prof['fit2'] = self.eped_prof(temp_psi,temp_coef[0],temp_coef[1],temp_coef[2],temp_coef[3],temp_coef[4],temp_coef[5])
			self.make_fitp('ti')
		else:	
			if self.post['opt_change']['ti']:
				self.post['didfit']['ti'] = False

		if ((self.forcefit or self.forcefit2) and self.post['isdat']['ti'] and self.hmodefit): 
			self.param['te']['val'][1] = self.post['popt']['ti'][1]
			self.param['te']['min'][1] = self.post['popt']['ti'][1]*0.8
			self.param['te']['max'][1] = self.post['popt']['ti'][1]*1.2
			self.param['te']['vary'][1]= True

		if (self.post['isdat']['ti']):
			if (self.fit_opt['use_ti_width']):
				if (self.fit_opt['func_type']['ti'] == 2 or self.fit_opt['func_type']['ti'] == 4):
					width = self.post['popt']['ti'][2]; pedmid = 1 - 0.5*self.post['popt']['ti'][2];
					if self.fit_opt['use_rho']['ti']: width, pedmid = self.pwidth2rwidth(width,pedmid)
					self.fit_opt['width_val']['te'] = round(width,3)
					self.fit_opt['width_val']['ne'] = round(width,3)
					self.fit_opt['width_val']['ti'] = round(width,3)
					self.post['pedmid']['te'] = round(pedmid,3)
					self.post['pedmid']['te'] = round(pedmid,3)
				if (self.fit_opt['func_type']['ti'] == 6):
					width = self.post['popt']['ti'][2]; pedmid = self.post['popt']['ti'][6];
					if self.fit_opt['use_rho']['ti']: width, pedmid = self.pwidth2rwidth(width,pedmid)
					self.fit_opt['width_val']['te'] = round(width,3)
					self.fit_opt['width_val']['ne'] = round(width,3)
					self.fit_opt['width_val']['ti'] = round(width,3)				
					self.post['pedmid']['te'] = round(pedmid,3)
					self.post['pedmid']['ne'] = round(pedmid,3)
					self.post['pedmid']['ti'] = round(pedmid,3)
				if (self.fit_opt['func_type']['ti'] == 3):
					width = self.post['popt']['ti'][5]; pedmid = 1-0.5*self.post['popt']['ti'][5];
					if self.fit_opt['use_rho']['ti']: width, pedmid = self.pwidth2rwidth(width,pedmid)
					self.fit_opt['width_val']['te'] = round(width,3)
					self.fit_opt['width_val']['ne'] = round(width,3)

		if self.ne_tcifit:
			for k in self.ne_list: self.fit_opt['weight']['ne'][k] = 0.

		if (self.forcefit or self.forcefit2):
			for k in self.ne_list: self.fit_opt['weight']['ne'][k] = 0. 
			self.fit_opt['use_ti_width']      = False; self.fit_opt['width_fix']['ne'] = False
			self.param['ne']['max'][2] = 0.1; self.param['ne']['min'][2] = 0.04; self.param['ne']['val'][2] = 0.05
				

		if (self.post['isdat']['ne'] and self.post['opt_change']['ne']):
			if self.post['didfit']['ne']: self.ne_prof['fit_old'] = np.copy(self.ne_prof['fit2p'])
			self.adjust_ts_scale()
			self.lmfit_fit('ne'); self.post['didfit']['ne'] = True;

			if (self.fit_opt['ped_scan_fit'] and self.fit_opt['use_rho']['ne']):
				temp_coef = self.modify_rho_to_psi_eped(self.post['popt']['ne'],'ne',True)
				temp_psi  = self.fit_eq['rho_to_psi'](self.fit_eq['psin2'])
				self.ne_prof['fit2'] = self.eped_prof(temp_psi,temp_coef[0],temp_coef[1],temp_coef[2],temp_coef[3],temp_coef[4],temp_coef[5])
			self.make_fitp('ne')	
		else: 
			if self.post['opt_change']['ne']:
				self.post['didfit']['ne'] = False

		if (self.forcefit or self.forcefit2):
			self.fit_opt['use_ti_width']      = True; self.fit_opt['width_fix']['te'] = True


		if (self.post['isdat']['te'] and self.post['opt_change']['te']):
			if self.post['didfit']['te']: self.te_prof['fit_old'] = np.copy(self.te_prof['fit2p'])			
			self.lmfit_fit('te'); self.post['didfit']['te'] = True;
			if (self.fit_opt['ped_scan_fit'] and self.fit_opt['use_rho']['te']):
				temp_coef = self.modify_rho_to_psi_eped(self.post['popt']['te'],'te',True)
				temp_psi  = self.fit_eq['rho_to_psi'](self.fit_eq['psin2'])
				self.te_prof['fit2'] = self.eped_prof(temp_psi,temp_coef[0],temp_coef[1],temp_coef[2],temp_coef[3],temp_coef[4],temp_coef[5])
			self.make_fitp('te')

		else:   
			if self.post['opt_change']['te']:
				self.post['didfit']['te'] = False	
		

		if (self.post['isdat']['vt'] and self.post['opt_change']['vt']):	
			if self.post['didfit']['vt']: self.vt_prof['fit_old'] = np.copy(self.vt_prof['fit2p'])
			self.lmfit_fit('vt'); self.post['didfit']['vt'] = True;
			self.make_fitp('vt')

		else:   
			if self.post['opt_change']['vt']:
				self.post['didfit']['vt'] = False
		
		if (self.post['isdat']['ne'] and self.post['opt_change']['ne']):
			if np.sum(self.ne_prof['fit_old'])>0: self.ne_prof['fit_old'] / self.post['scaled_density']
			if self.post['iskfile']['ne']: self.ne_prof['kfile']['raw'] = self.ne_prof['kfile']['raw']/self.post['scaled_density']
			self.den_scale()

		elif (self.post['isdat']['ne']):
			self.ne_prof['fit1']    = self.ne_prof['fit1']    / self.post['scaled_density']
			self.ne_prof['fit2']    = self.ne_prof['fit2']    / self.post['scaled_density']
			self.ne_prof['fit2p']   = self.ne_prof['fit2p']   / self.post['scaled_density']
			if np.sum(self.ne_prof['fit_old'])>0: self.ne_prof['fit_old'] / self.post['scaled_density']
			if self.post['iskfile']['ne']: self.ne_prof['kfile']['raw'] = self.ne_prof['kfile']['raw']/self.post['scaled_density']
			self.den_scale()
		if (self.post['isdat']['vt']): self.write_chease_rot()

		self.post['fit_opt'] = copy.deepcopy(self.fit_opt)

		self.first_run = False
		for flag in self.prof_list:
			self.post['func_history'][self.fit_opt['func_type'][flag]] = True

		self.post_process()
		if not self.noprint:
			self.print('=----------------------------------------------------------------------------=')
			self.print('=                              FITTING  DONE                                 =')
			self.print('=----------------------------------------------------------------------------=')
		return

	def check_options_change(self):

		for flag in self.prof_list: self.post['opt_change'][flag] = False

		for i in self.inter_list:
			if not self.fit_opt[i]['sig'] == 0.:
				if not self.fit_opt[i]['val'] == self.post['fit_opt'][i]['val']: self.post['opt_change']['ne'] = True
			if not self.fit_opt[i]['sig'] == self.post['fit_opt'][i]['sig']: self.post['opt_change']['ne'] = True

		if not self.fit_opt['use_ti_width']  == self.post['fit_opt']['use_ti_width']: 
			self.post['opt_change']['te'] = True; self.post['opt_change']['ne'] = True;

		if not self.fit_opt['file']['gfile'] == self.post['fit_opt']['file']['gfile']: 
			for flag in self.prof_list: self.post['opt_change'][flag] = True

		for i in self.prof_list:
			if not self.fit_opt['use_rho'][i] == self.post['fit_opt']['use_rho'][i]: self.post['opt_change'][i] = True
			if not self.fit_opt['sep_fix'][i] == self.post['fit_opt']['sep_fix'][i]: self.post['opt_change'][i] = True
			if not self.fit_opt['sep_val'][i] == self.post['fit_opt']['sep_val'][i]: self.post['opt_change'][i] = True
			if not self.fit_opt['width_fix'][i] == self.post['fit_opt']['width_fix'][i]: self.post['opt_change'][i] = True
			if not self.fit_opt['width_val'][i] == self.post['fit_opt']['width_val'][i]: self.post['opt_change'][i] = True
			if not self.fit_opt['raw_fit'][i] == self.post['fit_opt']['raw_fit'][i]: self.post['opt_change'][i] = True
			if not self.fit_opt['sspline_order'][i] == self.post['fit_opt']['sspline_order'][i]: self.post['opt_change'][i] = True

			if not self.fit_opt['oli'][i]['n']  == self.post['fit_opt']['oli'][i]['n']: self.post['opt_change'][i] = True
			if not self.fit_opt['file'][i]['kfile'] == self.post['fit_opt']['file'][i]['kfile']: self.post['opt_change'][i] = True

			for j in self.__dict__['%s_list'%i]:
				if not self.fit_opt['file'][i][j] == self.post['fit_opt']['file'][i][j]: self.post['opt_change'][i] = True
				if not self.fit_opt['avg'][i][j] == self.post['fit_opt']['avg'][i][j]: self.post['opt_change'][i] = True
				if not self.fit_opt['shift'][i][j] == self.post['fit_opt']['shift'][i][j]: self.post['opt_change'][i] = True
				if not self.fit_opt['exclude'][i][j] == self.post['fit_opt']['exclude'][i][j]: self.post['opt_change'][i] = True
				if not self.fit_opt['oli'][i][j]['use'] == self.post['fit_opt']['oli'][i][j]['use']: self.post['opt_change'][i] = True
				if not self.fit_opt['oli'][i][j]['per'] == self.post['fit_opt']['oli'][i][j]['per']: self.post['opt_change'][i] = True
				if not self.fit_opt['std'][i][j]['use'] == self.post['fit_opt']['std'][i][j]['use']: self.post['opt_change'][i] = True
				if not self.fit_opt['std'][i][j]['raw'] == self.post['fit_opt']['std'][i][j]['raw']: self.post['opt_change'][i] = True
				if not self.fit_opt['weight'][i][j] == self.post['fit_opt']['weight'][i][j]: self.post['opt_change'][i] = True
				if not self.fit_opt['psi_end'][i][j] == self.post['fit_opt']['psi_end'][i][j]: self.post['opt_change'][i] = True

				for k in ['core','edge','ch_cal','core_ch','edge_ch','ratio','fix_core']:
					if not self.fit_opt['scale'][i][j][k] == self.post['fit_opt']['scale'][i][j][k]:
						self.post['opt_change'][i] = True

			if not self.fit_opt['func_type'][i] == self.post['fit_opt']['func_type'][i]: 
				self.post['opt_change'][i] = True
				try: 
					if not self.post['func_history'][i][self.fit_opt['func_type'][i]]: self.lmfit_init_param(self.fit_opt['func_type'][i],i)
				except: pass

			if not self.compare_dict(self.param[i],self.param_old[i]):self.post['opt_change'][i] = True

		for flag in self.prof_list:
			if (self.fit_opt['func_type'][flag] == 4 and self.post['fit_opt']['func_type'][flag] == 6):
				for i in ['vary','val','min','max']:
					for j in range(6): self.param[flag][i][j] = self.param_old[flag][i][j]

			if (self.fit_opt['func_type'][flag] == 6 and self.post['fit_opt']['func_type'][flag] == 4):
				for i in ['vary','val','min','max']:
					for j in range(6): self.param[flag][i][j] = self.param_old[flag][i][j]

		self.param_old = copy.deepcopy(self.param)

		for flag in self.prof_list:
			if self.post['opt_change'][flag]: self.print('>>> %s fit option is changed'%flag.upper())

		if not self.fit_opt['ascale'] == self.post['fit_opt']['ascale']: self.post['opt_change']['ne'] = True
		try: 
			if not self.fit_opt['ascale1d'] == self.post['fit_opt']['ascale1d']: self.post['opt_change']['ne'] = True
		except: pass
		if self.fit_opt['ascale']:
			for flag in ['minc','maxc','mine','maxe','cn','en']:
				if not self.fit_opt['scale']['ne']['ts'][flag] == self.post['fit_opt']['scale']['ne']['ts'][flag]: 
					self.post['opt_change']['ne'] = True

		return

	def compare_dict(self,a,b):

		ans = True

		for i in ['vary','val','min','max']:
			for j in range(len(a[i])):
				if not np.isnan(a[i][j]):
					if not a[i][j] == b[i][j]:
						ans = False
				else:
					if not np.isnan(b[i][j]):
						ans = False

		return ans

	def post_process(self):

		if not self.post['isdat']['ne']: return
		ne = np.copy(self.ne_prof['fit2'])
		ne2p= np.copy(self.ne_prof['fit2p'])
		xx = np.copy(self.fit_eq['psin2'])
		self.scale_density(use_rho=self.fit_opt['use_rho']['ne'],print_only=True)
		if not (self.post['isdat']['te'] and self.post['isdat']['ti']): return
		elif not (self.post['isdat']['te']):
			te = np.copy(self.ti_prof['fit2'])
			ti = np.copy(self.ti_prof['fit2'])
			te2p= np.copy(self.ti_prof['fit2p'])
			ti2p= np.copy(self.ti_prof['fit2p'])
		elif not (self.post['isdat']['ti']):
			te = np.copy(self.te_prof['fit2'])
			ti = np.copy(self.te_prof['fit2'])
			te2p= np.copy(self.te_prof['fit2p'])
			ti2p= np.copy(self.te_prof['fit2p'])
		else:
			te = np.copy(self.te_prof['fit2'])
			ti = np.copy(self.ti_prof['fit2'])
			te2p= np.copy(self.te_prof['fit2p'])
			ti2p= np.copy(self.ti_prof['fit2p'])

		if not (self.post['isdat']['vt']):
			vt = np.copy(self.fit_eq['psin2'])
		else:
			vt = np.copy(self.vt_prof['fit2'])

		factore = 1.602*1.e-19*1.e19*1.e3
		factori = (1.-(self.fit_opt['zeff']-1.)/self.fit_opt['zimp'])*1.602*1.e-19*1.e19*1.e3
		nef = interp1d(xx,ne,'cubic')
		tef = interp1d(xx,te,'cubic')
		tif = interp1d(xx,ti,'cubic')
		vtf = interp1d(xx,vt,'cubic')
		#pef = interp1d(xx,pe,'cubic')
		#pif = interp1d(xx,pi,'cubic')

		tiw = 0.05; tew = 0.05; new = 0.05; vtw = 0.05;	tiwe = 0.;	vtwe = 0;	newe = 0.;	tewe = 0.
		tip = 0.95; tep = 0.95; nep = 0.95; vtp = 0.95;	tihe = 0.;	vthe = 0.;	nehe = 0.;	tehe = 0.

		if (self.fit_opt['func_type']['ti'] == 3):
			tiw = self.post['popt']['ti'][5];	tip = self.post['popt']['ti'][6] - 0.5*tiw;   tiwe = self.post['popte']['ti'][2];
			s1  = self.post['popte']['ti'][2] * np.sign(self.post['popt']['ti'][2]);	
			s2  = self.post['popte']['ti'][3] * np.sign(self.post['popt']['ti'][3]); 
			s3 = self.post['popte']['ti'][4]* np.sign(self.post['popt']['ti'][4]);
			tihe = (self.post['popt']['ti'][1]-self.post['popt']['ti'][0])*(1+np.tanh(1))*0.5*(s1*tip+s2*tip**2+s3*tip**3)
		elif (self.fit_opt['func_type']['ti'] == 6):
			tiw = self.post['popt']['ti'][2];	tip = self.post['popt']['ti'][6] - 0.5*tiw;	tiwe = self.post['popte']['ti'][2];
			tihe = 2.*np.tanh(1)*self.post['popte']['ti'][1]	
		elif (self.fit_opt['func_type']['ti'] == 4):
			tiw = self.post['popt']['ti'][2];	tip = 1.0 - tiw;					
			tiwe = self.post['popte']['ti'][2];	tihe = 2.*np.tanh(1)*self.post['popte']['ti'][1]
		elif (self.fit_opt['func_type']['ti'] == 2):
			tiw = self.post['popt']['ti'][2];	tip = self.post['popt']['ti'][3] - 0.5*tiw;	tiwe = self.post['popte']['ti'][2];	tihe = self.post['popte']['ti'][1]
		tihe = abs(tihe)

		if (self.fit_opt['func_type']['vt'] == 3):
			vtw = self.post['popt']['vt'][5];	vtp = self.post['popt']['vt'][6] - 0.5*vtw;	vtwe = self.post['popte']['vt'][5];	
			s1 = self.post['popte']['vt'][2]* np.sign(self.post['popt']['vt'][2]);	
			s2 = self.post['popte']['vt'][3]* np.sign(self.post['popt']['vt'][3]); 
			s3 = self.post['popte']['vt'][4]* np.sign(self.post['popt']['vt'][4]);
			vthe = (self.post['popt']['vt'][1]-self.post['popt']['vt'][0])*(1+np.tanh(1))*0.5*(s1*vtp+s2*vtp**2+s3*vtp**3)
		elif (self.fit_opt['func_type']['vt'] == 6):
			vtw = self.post['popt']['vt'][2];	vtp = self.post['popt']['vt'][6] - 0.5*vtw;	vtwe = self.post['popte']['vt'][2];	vthe = 2.*np.tanh(1)*self.post['popte']['vt'][1]			
		elif (self.fit_opt['func_type']['vt'] == 4):
			vtw = self.post['popt']['vt'][2];	vtp = 1.0 - vtw;			vtwe = self.post['popte']['vt'][2];	vthe = 2.*np.tanh(1)*self.post['popte']['vt'][1]
		elif (self.fit_opt['func_type']['vt'] == 2):
			vtw = self.post['popt']['vt'][2];	vtp = self.post['popt']['vt'][3] - 0.5*vtw;	vtwe = self.post['popte']['vt'][2];	vthe = self.post['popte']['vt'][1]
		vthe = abs(vthe)

		if (self.fit_opt['func_type']['te'] == 3):
			tew = self.post['popt']['te'][5];	tep = self.post['popt']['te'][6] - 0.5*tew;	tewe = self.post['popte']['te'][5];	
			s1 = self.post['popte']['te'][2]* np.sign(self.post['popt']['te'][2]);	
			s2 = self.post['popte']['te'][3]* np.sign(self.post['popt']['te'][3]); 
			s3 = self.post['popte']['te'][4]* np.sign(self.post['popt']['te'][4]);
			tehe = (self.post['popt']['te'][1]-self.post['popt']['te'][0])*(1+np.tanh(1))*0.5*(s1*tep+s2*tep**2+s3*tep**3)
		elif (self.fit_opt['func_type']['te'] == 6):
			tew = self.post['popt']['te'][2];	tep = self.post['popt']['te'][6] - 0.5*tew;	tewe = self.post['popte']['te'][2];	tehe = 2.*np.tanh(1)*self.post['popte']['te'][1]		
		elif (self.fit_opt['func_type']['te'] == 4):
			tew = self.post['popt']['te'][2];	tep = 1.0 - tew;					tewe = self.post['popte']['te'][2];	tehe = 2.*np.tanh(1)*self.post['popte']['te'][1]
		elif (self.fit_opt['func_type']['te'] == 2):
			tew = self.post['popt']['te'][2];	tep = self.post['popt']['te'][3] - 0.5*tew;	tewe = self.post['popte']['te'][2];	tehe = self.post['popte']['te'][1]
		tehe = abs(tehe)

		if (self.fit_opt['func_type']['ne'] == 3):
			new = self.post['popt']['ne'][5];	nep = self.post['popt']['ne'][6] - 0.5*new;	newe = self.post['popte']['ne'][5];	
			s1 = self.post['popte']['ne'][2]* np.sign(self.post['popt']['ne'][2]);	
			s2 = self.post['popte']['ne'][3]*np.sign(self.post['popt']['ne'][3]); 
			s3 = self.post['popte']['ne'][4]*np.sign(self.post['popt']['ne'][4]);
			nehe = (self.post['popt']['ne'][1]-self.post['popt']['ne'][0])*(1+np.tanh(1))*0.5*(s1*nep+s2*nep**2+s3*nep**3)
		elif (self.fit_opt['func_type']['ne'] == 6):
			new = self.post['popt']['ne'][2];	nep = self.post['popt']['ne'][6] - 0.5*new;	newe = self.post['popte']['ne'][2];	nehe = 2.*np.tanh(1)*self.post['popte']['ne'][1]			
		elif (self.fit_opt['func_type']['ne'] == 4):
			new = self.post['popt']['ne'][2];	nep = 1.0 - new;					newe = self.post['popte']['ne'][2];	nehe = 2.*np.tanh(1)*self.post['popte']['ne'][1]
		elif (self.fit_opt['func_type']['ne'] == 2):
			new = self.post['popt']['ne'][2];	nep = self.post['popt']['ne'][3] - 0.5*new;	newe = self.post['popte']['ne'][2];	nehe = self.post['popte']['ne'][1]
		nehe = abs(nehe)

		neped = nef(nep); teped = tef(tep);	tiped = tif(tip); vtped = vtf(vtp);
		necor = nef(0.0); tecor = tef(0.0); ticor = tif(0.0); vtcor = vtf(0.0);
		nesep = nef(1.0); tesep = tef(1.0); tisep = tif(1.0); vtsep = vtf(1.0);
		nesol = nef(1.2); tesol = tef(1.2); tisol = tif(1.2); vtsol = vtf(1.2);
		peped = nef(0.5*(nep+tep)) * tef(0.5*(nep+tep))*factore;	
		piped = nef(0.5*(nep+tep)) * tef(0.5*(nep+tep))*factori; 
		ptped = peped + piped
		#pef(0.25*(2.*nep+tip+tep)) + pif(0.25*(2.*nep+tip+tep));

		if self.fit_opt['use_rho']['te']: tew = self.fit_eq['rho_to_psi'](tep+tew) - self.fit_eq['rho_to_psi'](tep)
		if self.fit_opt['use_rho']['ne']: new = self.fit_eq['rho_to_psi'](nep+new) - self.fit_eq['rho_to_psi'](nep)
		if self.fit_opt['use_rho']['ti']: tiw = self.fit_eq['rho_to_psi'](tip+tiw) - self.fit_eq['rho_to_psi'](tip)
		if self.fit_opt['use_rho']['vt']: vtw = self.fit_eq['rho_to_psi'](vtp+vtw) - self.fit_eq['rho_to_psi'](vtp)

		pepede = peped*(tehe/teped + nehe/neped);	pipede = piped*(tihe/tiped + nehe/neped); ptpede = pepede+pipede
		we1 = 0.5*(newe+tewe);	we2 = 0.5*(newe+tiwe);	we3 = 0.25*(2*newe+tewe+tiwe);

		psin2 = np.linspace(0,1.,201)
		
		vol = self.fit_eq['Vf'](psin2)
		nepf = interp1d(self.fit_eq['psin2'],ne2p)
		tepf = interp1d(self.fit_eq['psin2'],te2p)
		tipf = interp1d(self.fit_eq['psin2'],ti2p)
		pt =  nepf(psin2)*(tepf(psin2)*factore + tipf(psin2)*factori)

		wkin = np.trapz(pt,x=vol) * 1.5
		wratio = (1.0-self.fit_eq['psi_to_rho'](1.0-tew))/tew
		
		bpped = 4*np.pi*1.e-7 * abs(self.fit_eq['eq'].ip)/self.fit_eq['eq'].perim
		bpped = ptped * 2. * 4*np.pi*1.e-7 / (bpped**2)
		bppede = ptpede/ptped * bpped
		coef = 0.5*(new+tew) / (bpped**0.5)
		coefe = coef * (we1*2/(new+tew) + 0.5*bppede/bpped)

		self.print('=----------------------------------------------------------------------------=')
		self.print('                       TE [keV]    NE [10(19)/m3]   TI [keV]    VT [km/s]    ')
		self.print('>>> PED  HEIGHT         %-5.3f         %-5.3f          %-5.3f        %-6.3f'%(teped,neped,tiped,vtped))
		self.print('>>> PED  HEIGHT(ERR)    %-5.3f         %-5.3f          %-5.3f        %-6.3f'%(tehe,nehe,tihe,vthe))

		self.print('>>> PED  WIDTH (PSI)    %-5.3f         %-5.3f          %-5.3f        %-6.3f'%(tew,new,tiw,vtw))
		self.print('>>> PED  WIDTH (ERR)    %-5.3f         %-5.3f          %-5.3f        %-6.3f'%(tewe,newe,tiwe,vtwe))			
		self.print('>>> PED    LOC (PSI)    %-5.3f         %-5.3f          %-5.3f        %-6.3f'%(tep,nep,tip,vtp))

		self.print('>>> SEP  VALUE          %-5.3f         %-5.3f          %-5.3f        %-6.3f'%(tesep,nesep,tisep,vtsep))
		self.print('>>> SOL  VALUE          %-5.3f         %-5.3f          %-5.3f        %-6.3f'%(tesol,nesol,tisol,vtsol))
		self.print('>>> CORE VALUE          %-5.3f         %-5.3f          %-5.3f        %-6.3f'%(tecor,necor,ticor,vtcor))
		self.print('=----------------------------------------------------------------------------=')

		self.print('>>> PE   PED HEIGHT [kPa]   %-6.3f +/- [%-5.3f]   at PSI = %-5.3f  (width = %-5.3f +/- [%-5.3f] )'%(peped/1.e3,pepede/1.e3,0.5*(nep+tep),0.5*(new+tew),we1))
		self.print('>>> PI   PED HEIGHT [kPa]   %-6.3f +/- [%-5.3f]   at PSI = %-5.3f  (width = %-5.3f +/- [%-5.3f] )'%(piped/1.e3,pipede/1.e3,0.5*(nep+tip),0.5*(new+tiw),we2))
		self.print('>>> PTOT PED HEIGHT [kPa]   %-6.3f +/- [%-5.3f]   at PSI = %-5.3f  (width = %-5.3f +/- [%-5.3f] )'%(ptped/1.e3,ptpede/1.e3,0.25*(2.*nep+tip+tep),0.25*(2.*new+tiw+tew),we3))
		self.print('>>> RHO_WIDTH/PSI_WIDTH     %-6.3f  '%wratio)
		self.print('>>> W_THERMAL/MHD [kJ]      %-6.3f / %-6.3f (%-4.2f[%%])  W_EXT [kJ] -> %-6.3f '%(wkin/1.e3,self.fit_eq['eq'].wmhd/1.e3,wkin/self.fit_eq['eq'].wmhd*1.e2,(self.fit_eq['eq'].wmhd-wkin)/1.e3))
		self.print('>>> BP_PED / WIDTH          %-6.3f +/- [%-5.3f]  /  %-6.3f +/- [%-5.3f] '%(bpped,bppede,0.5*(new+tew),we1))
		self.print('>>> KBM_COEFF               %-6.3f +/- [%-5.3f]  (Standard coeff-GA -> (0.076 with c2 = 0.5]) '%(coef,coefe))
		self.post['bpped'] = bpped
		self.post['bppede'] = bppede
		self.post['wmhd'] = self.fit_eq['eq'].wmhd/1.e3
		self.post['wkin'] = wkin/1.e3
		self.post['width'] = [tew,new,tiw,vtw]
		self.post['widthe'] = [tewe,newe,tiwe,vtwe]

		return

	def forced_fit(self,ishmode = True):

		self.hmodefit = ishmode
		self.param['te']['min'][4] = 0.5
		self.param['te']['min'][5] = 0.5
		self.param['ti']['min'][4] = 0.5
		self.param['ti']['min'][5] = 0.5
		if ishmode:
			self.fit_opt['ashift']['min']	  = 0.9
			self.fit_opt['ashift']['max']	  = 1.05			
			self.fit_opt['use_ti_width']      = True
			self.fit_opt['use_ti_eped']       = True
			self.fit_opt['ped_scan_fit']      = True			
			self.fit_opt['func_type']['vt']   = 3
			self.param['ne']['min'][4] = 1.05
			self.param['ne']['max'][4] = 2.6
			self.param['ne']['min'][5] = 1.05
			self.param['ne']['max'][5] = 2.6
		else:
			self.fit_opt['func_type']['te']   = 7
			self.fit_opt['sspline_order']['te'] = 1
			self.fit_opt['func_type']['ne']   = 4
			self.fit_opt['func_type']['ti']   = 7
			self.fit_opt['sspline_order']['ti'] = 3
			self.fit_opt['psi_end']['ti']['ces'] = 1.05
			self.fit_opt['func_type']['vt']   = 3


		for i in ['te','ti']:
			self.fit_opt['sep_fix'][i]       = True	

		for i in ['ti','vt','te']:
			for j in self.__dict__['%s_list'%i]:
				self.fit_opt['avg'][i][j]		 = True

		for i in ['ne']:
			self.fit_opt['oli'][i]['n']  = 2			
			for j in self.__dict__['%s_list'%i]:
				self.fit_opt['avg'][i][j]		 = False
				self.fit_opt['oli'][i][j]['use'] = True
				self.fit_opt['oli'][i][j]['per'] = 0.95
				self.fit_opt['std'][i][j]['use'] = True
				self.fit_opt['std'][i][j]['raw'] = True				

		for i in self.prof_list:	
			self.fit_opt['raw_fit'][i]       = True
			for j in self.__dict__['%s_list'%i]:
				self.fit_opt['psi_end'][i][j]    = 1.05

		return

	def forced_fit_params(self):

		self.lmfit_init_params()
		if self.hmodefit:
			self.param['ti']['min'][2] = self.force_width
			if ((self.param['ne']['min'][0] == 0.30) and (self.param['ne']['max'][0] == 1.0) and (self.param['ne']['val'][0] == 0.50)):
				self.fit_opt['sep_fix']['ne'] = True
				self.fit_opt['sep_val']['ne'] = 0.4
				self.param['ne']['min'][0] = 0.30
				self.param['ne']['max'][0] = 0.60
				self.param['ne']['val'][0] = 0.50
				self.param['ne']['vary'][0] = False
			self.param['ti']['min'][4] = 0.5
			self.param['ti']['min'][5] = 0.5

		else:
			self.param['ne']['min'][2] = 0.15
			self.param['ne']['max'][2] = 0.3
			self.param['ne']['min'][0] = 0.1
			self.param['ne']['max'][0] = 0.3
			self.param['ne']['val'][0] = 0.2

		return

	def __init__(self,nolog=False):

		self.initialise_variables()
		self.nolog = nolog

		try:
			if os.path.isfile('fit_opt'): self.read_namelist('fit_opt')
		except: self.print('>>> No fit_opt...')
		else: self.print('>>> No fit_opt...')

		chk = fit_checkfile.chk_file()
		self.print('>>> Check Raw data files...')

		for flag in self.prof_list:
			for k in self.__dict__['%s_list'%flag]:
				filename = self.fit_opt['file'][flag][k]
				ischk = False
				if not (filename == None): ischk = chk.writefiles(filename,flag)
				if ischk:	self.print('>>> Source file for %s [%s] are saved to xxx_source...'%(flag,k.upper()))

		self.lmfit_init_params()
	
		fit_types = [self.fit_opt['func_type']['te'],self.fit_opt['func_type']['ne'],self.fit_opt['func_type']['ti'],self.fit_opt['func_type']['vt']]
		self.param2 = copy.deepcopy(self.param)
		self.param2_old = copy.deepcopy(self.param)
		self.param_change = False

		return		

	if __name__ == "__main__":

		import fittool
		fit = fittool.fit_tool(sys.argv[1])
