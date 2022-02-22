import os, sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.interpolate import LSQUnivariateSpline as sspline
import lmfit

class knots:

	def read_file(self,filename):


		f = open(filename,'r')
		dat_num = int(float(f.readline()))

		datx = np.zeros(dat_num)
		daty = np.zeros(dat_num)

		for i in range(dat_num):
			line = f.readline().split()
			datx[i] = float(line[0])
			daty[i] = float(line[1])

		f.close()

		return (datx,daty)


	def lmfit_params(self,datx,daty,type=2):

		p = lmfit.Parameters()
		self.ppmin_psi = 1.0
		if (type == 2):
			s = np.sign(daty[-10])
			datx = datx
			daty = daty * s

			ind = np.where(datx > 0.8)
			datx = datx[ind]
			daty = daty[ind]
			
			ind = np.argmax(daty)
			xmax = datx[ind]
			datx = datx[0:ind]
			daty = daty[0:ind]

			ind = np.argmin(daty)
			xmin = datx[ind]

			core_end = xmin - 0.2
			if self.minloc == 0.:
				ped_start = min(xmin+0.05,0.92)	
			else:
				ped_start = self.minloc
			self.ppmin_psi = ped_start
			self.xmin = xmin
			self.xmax = xmax

			self.varn = self.coren + self.edgen

			if self.spline_end == 1.:	a = np.linspace(ped_start+0.00,self.spline_end,self.edgen+1)
			else:	a = np.linspace(ped_start+0.00,self.spline_end,self.edgen)
			a1 = min(0.3,(0.2+core_end)*0.5)
			a1 = max(self.spline_start,a1)
			a2 = np.linspace(a1,ped_start,self.coren+1)

			self.init = np.zeros(self.varn+2)
			self.init[-1] = 1.0
			for i in range(self.coren):
				self.init[i+1] = a2[i]
			for i in range(self.edgen):
				self.init[i+self.coren+1] = a[i]

			p.add('d1',value=self.init[1]-self.spline_start,vary=True)

			for i in range(1,self.varn+1):
				val = self.init[i+1] - self.init[i]
				if i == self.coren + 1:	val = (1.0-ped_start)/(self.edgen-1)*2
				if i == self.coren + 2:	val = (1.0-ped_start)/(self.edgen-1)*1
				p.add('d%i'%i,value=val,vary=True)

		else:
			self.init = np.zeros(self.coren+2)
			self.init[1:self.coren+2] = np.linspace(self.spline_start,1.0,self.coren+1)

			self.varn = self.coren
			for i in range(1,self.varn+1):
				val = self.init[i+1] - self.init[i]
				p.add('d%i'%i,value=val,vary=True)

		return p

	def residual(self,p,x,y,init_ind=0,delx=0.01):

		t = []
		t2 = []
		sum1 = 0.
		for i in range(init_ind,self.varn2+init_ind):

			var_name = 'd%i'%(i+1)
			sum1 = sum1 + abs(p[var_name].value)

		len2 = self.knot_end - self.knot_start
		sum2 = 0.
		dels = sum1 * delx / (len2 - self.varn2*delx)
		sum1 = sum1 + dels*self.varn2
		for i in range(init_ind,self.varn2+init_ind):
			var_name = 'd%i'%(i+1)
			sum2 = sum2 + abs(p[var_name])+dels
			t.append(sum2/sum1*len2+self.knot_start)
			t2.append(p[var_name])

		ind = np.where(x>=self.spline_start)
		x0 = x[ind]
		y0 = y[ind]
		ind = np.where(x0<=self.spline_end)
		x0 = x0[ind]
		y0 = y0[ind]

		try:
			#yff = sspline(x,y,t,bbox=[self.spline_start2+1.e-3,self.spline_end],k=2)
			yff = sspline(x0,y0,t,k=2)
		except:
			print(x[0],t,x[-1])

		ind = np.where(x>self.lmfit_start)
		x2 = x[ind]
		y2 = y[ind]
		ind = np.where(x2<self.lmfit_end)
		x2 = x2[ind]
		y2 = y2[ind]

		yf = yff(x2)
		yt = (yf - y2)**2

		self.dels = dels

		return (yt)

	def lmfit_params2array(self,p,init_ind,var_n):
	
		t = []
		for i in range(init_ind,init_ind+var_n):
			var_name = 'd%i'%(i+1)
			t.append(p[var_name].value)
		sum1  = np.zeros(var_n)
		sum1[0] = abs(float(t[0]))+self.dels
		for i in range(1,var_n):
			sum1[i] = sum1[i-1] + abs(float(t[i])) + self.dels

		sum1 = self.knot_start + (sum1)/sum1[-1]*(self.knot_end-self.knot_start)
		return sum1

	def lmfit_fit(self,datx,daty,type=2,isp=True):

		self.print_input()
		datx2 = np.copy(datx)
		daty2 = np.copy(daty)

		datx = np.linspace(self.spline_start,1.0,self.input_datan)
		datf = interp1d(datx2,daty2,'cubic')		
		daty = datf(datx)
		p = self.lmfit_params(datx,daty,type)

		knots = [0.]

		if type == 2:

			delx = self.delmin_core
			self.knot_start = self.start_knot			#first do the core
			self.knot_end = self.ppmin_psi-self.knots_shift
			self.spline_end = min(1.0,self.ppmin_psi+0.1)
			self.lmfit_start = self.spline_start
			self.lmfit_end = self.ppmin_psi
			self.varn2 = self.coren
			self.spline_start2 = self.spline_start
		
			delx2 = min(delx,(self.knot_end-self.knot_start)/float(self.varn2)-0.001)
	
			result = lmfit.minimize(self.residual,p,args=(datx,daty,0,delx2))
			p = result.params
			t2 = self.lmfit_params2array(p,0,self.coren)
			for i in range(self.coren):	knots.append(t2[i])

			delx = self.delmin_edge
			self.knot_start = self.ppmin_psi			#second do the edge
			self.knot_end = self.end_knot
			self.spline_end = 1.0
			self.lmfit_start = self.ppmin_psi
			self.lmfit_end = 1.0
			self.spline_start2 =self.ppmin_psi-self.knots_shift
			self.varn2 = self.edgen

			datx = np.linspace(self.spline_start2,1.0,self.input_datan)
			datf = interp1d(datx2,daty2,'cubic')		
			daty = datf(datx)			
			p = self.lmfit_params(datx,daty)

			delx2 = min(delx,(self.knot_end-self.knot_start)/float(self.varn2)-0.001)
			if not delx2 == delx:	print('DELX %s -> %s'%(delx,delx2))
			
			result = lmfit.minimize(self.residual,p,args=(datx,daty,self.coren,delx2))
			p = result.params
			t2 = self.lmfit_params2array(p,self.coren,self.edgen)
			for i in range(self.edgen):	knots.append(t2[i])

		else:
		
			delx = self.delmin_core
			self.knot_start = self.start_knot			#first do the core
			self.knot_end = self.end_knot
			self.spline_end = 1.0
			self.lmfit_start = self.spline_start
			self.lmfit_end = 1.0
			self.varn2 = self.coren
			self.spline_start2 = self.spline_start
			
			delx = min(delx,(self.knot_end-self.knot_start)/float(self.varn2)-0.001)

			result = lmfit.minimize(self.residual,p,args=(datx,daty,0,delx))
			p = result.params
			t2 = self.lmfit_params2array(p,0,self.coren)
			for i in range(self.coren):	knots.append(t2[i])

		knots.append(1.)
		self.print_result(result,knots,type,isp)

		if (type == 2):
			if not (isp):
				if self.ppmin_psi >= 0.905:	
					knots[self.coren] = 0.85
				else:
					knots[self.coren] = 0.85 + (self.ppmin_psi - 0.905)

		return knots

	def print_result(self,result,knots,type,isp):

		if (type == 1):
			flag = 'L-mode'
			ll = self.varn+2
		elif (type == 2):
			flag = 'H-mode'
			ll = self.varn+2

		if isp:
			flag2 = "P'"
		else:
			flag2 = "FF'"

		s0 = ''
		s1 = ''
		for i in range(ll):
			s0 = s0 + '%4.3f,'%self.init[i]
			s1 = s1 + '%4.3f,'%knots[i]

		print('=----------------------------------------------------------------------------=')
		print('=                             EFIT KNOTS FINDER                              =')	
		print('=----------------------------------------------------------------------------=')	
		print('>>> Knots finding for %s in %s, ped_start = %f'%(flag2,flag,self.ppmin_psi))
		print('>>> XI2 = %6.3f, RED_XI2 = %6.4e'%(result.chisqr,result.redchi))
		print('=----------------------------------------------------------------------------=')
		result.params.pretty_print()
		print('=----------------------------------------------------------------------------=')
		print('>>> Initial guess = %s'%s0)		
		print('>>> Final knots   = %s'%s1)
		print('=----------------------------------------------------------------------------=')
		return

	def set_default(self):

		self.spline_start = 0.2			# spline start
		self.end_knot = 0.985			# knots end
		self.start_knot = 0.3			# knots start
		self.knots_shift = 0.05			# knots shift
		self.coren = 2					# core knots #
		self.edgen = 5					# edge knots #
		self.spline_end = 1.
		self.input_datan =501			# input_datan

		self.delmin_core = 0.1
		self.delmin_edge = 0.01

		self.minloc = 0.

		return

	def print_input(self):
	
		label = ['SPLINE_START','END_KNOT','START_KNOT','KNOT_SHIFT','CORE_KNOT','EDGE_KNOT','SPLINE_PT','DXMIN_CORE','DXMIN_EDGE']
		
		line = '--------- INPUT -------->\n'
		line = line + '%13s\t%s\n'%(label[0],self.spline_start)
		line = line + '%13s\t%s\n'%(label[1],self.end_knot)
		line = line + '%13s\t%s\n'%(label[2],self.start_knot)
		line = line + '%13s\t%s\n'%(label[3],self.knots_shift)
		line = line + '%13s\t%s\n'%(label[4],self.coren)
		line = line + '%13s\t%s\n'%(label[5],self.edgen)
		line = line + '%13s\t%s\n'%(label[6],self.input_datan)
		line = line + '%13s\t%s\n'%(label[7],self.delmin_core)
		line = line + '%13s\t%s'%(label[8],self.delmin_edge)
		print(line)
		return

	def __init__(self):

		self.set_default()

		return

	if __name__ == "__main__":

		import knots_tool3
		ch = knots_tool3.knots()

		#d = np.load('sample.npz')
		#datx = np.copy(d['datx'])
		#temp = np.copy(d['daty'])
		#daty = np.copy(datx)	

		f = open('pre_prof','r')
		line = f.readline()
		num = int(float(line))
		datx = np.zeros(num);	temp = np.copy(datx);	daty = np.copy(datx)

		for i in range(num):
			line = f.readline().split()
			datx[i] = float(line[0])
			temp[i] = float(line[1])
		f.close()

		yy = interp1d(datx,temp,'cubic')
		eps = 1.e-5

		for i in range(len(datx)-2):
			xx = datx[i+1]
			daty[i+1] = (yy(xx+eps) - yy(xx-eps))/2./eps
		daty[0] = 2.*daty[1] - daty[2]
		daty[-1] = 2.*daty[-2] - daty[-3]

		daty = daty * -1.

		xx = np.linspace(0,1.,501)
		yf = interp1d(datx,daty,'cubic')
		yy = yf(xx)

		ch.lmfit_fit(datx,daty,2,True)





