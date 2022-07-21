#!/usr/local/anaconda3/bin/python3
import os, sys
import numpy as np
import tkinter as tk
import tkinter.font
import tkinter.scrolledtext as tkst
from shutil import move, copyfile, copytree, rmtree
from tkinter.filedialog import askopenfilename,asksaveasfilename, askdirectory
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import matplotlib.pyplot as plt
import eqdsk
import subprocess
from gefit_tool import *
from exec_dirs import gfitp_dir,chease_dir,nubeam_dir2,python2_exec,python3_exec,efit_dir,version,author,comment,efit_rmp
from exec_dirs import mse_good_ch, years, shotk
from gefit_mse import *
from get_mse_data import _get_mse, _make_mse
import pfile_convert

class gefitk:

	def declare_vars(self):
		#Vars
		for i in range(1,40):
			self.__dict__['CheckVar%d'%i] = tk.IntVar()
			self.__dict__['CheckVar%d'%i].set(0)

		for i in range(1,250):
			self.__dict__['CoilVar%d'%i] = tk.IntVar()
			self.__dict__['CoilVar%d'%i].set(0)			

		for i in range(1,200):
			self.__dict__['StrVar%d'%i] = tk.StringVar()
			self.__dict__['StrVar%d'%i].set('')

		for i in range(1,20):
			self.__dict__['MenuVar%d'%i] = tk.StringVar()
			self.__dict__['MenuVar%d'%i].set('')			

		for i in range(1,10):
			self.__dict__['DoubleVar%d'%i] = tk.DoubleVar()
			self.__dict__['DoubleVar%d'%i].set(0)

		return

	def save_opt(self,filename):

		for i in range(1,14):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%i].get())

		for i in range(45,46):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%i].get())			

		f = open(filename,'w')
		for i in range(1,40):
			f.write('%i\t%s%i\n'%(self.__dict__['CheckVar%d'%i].get(),'!CHECK',i))

		for i in range(1,250):
			f.write('%i\t%s%i\n'%(self.__dict__['CoilVar%d'%i].get(),'!COIL',i))	

		for i in range(1,200):
			f.write('%s\t%s%i\n'%(self.__dict__['StrVar%d'%i].get(),'!STR',i))

		for i in range(3,20):
			f.write('%s\t%s%i\n'%(self.__dict__['MenuVar%d'%i].get(),'!MENU',i))	

		for i in range(1,10):
			f.write('%f\t%s%i\n'%(self.__dict__['DoubleVar%d'%i].get(),'!DOUBLE',i))
		f.close()

	def load_opt(self,filename):

		if not os.path.isfile(filename):
			print('>>> No saved params...')
			return True

		f = open(filename,'r')
		for i in range(1,40):
			line = f.readline().split('!')[0]
			self.__dict__['CheckVar%d'%i].set(int(float(line)))

		for i in range(1,250):
			line = f.readline().split('!')[0]
			self.__dict__['CoilVar%d'%i].set(int(float(line)))			

		for i in range(1,200):
			try:	line = f.readline().split('\n')[0].split('!')[0].split()[0]
			except:	line = ''
			self.__dict__['StrVar%d'%i].set(line)

		for i in range(3,20):
			try:	line = f.readline().split('\n')[0].split('!')[0].split()[0]
			except:	line = ''
			self.__dict__['MenuVar%d'%i].set(line)			

		for i in range(1,10):
			line = f.readline().split('!')[0]
			self.__dict__['DoubleVar%d'%i].set(float(line))

		f.close()

		for i in range(1,14):

			self.__dict__['e%d'%i].delete(0,'end')
			self.__dict__['e%d'%i].insert(10,self.__dict__['StrVar%d'%i].get())

		for i in range(45,46):
			self.__dict__['e%d'%i].delete(0,'end')
			self.__dict__['e%d'%i].insert(10,self.__dict__['StrVar%d'%i].get())
		
		self.bdiff_index = find_bdiff_index(self,1)
		self.bdiff_index = find_bdiff_index(self,2)
		self.update_status()
		if self.StrVar99.get() == '': self.StrVar99.set('2.0');
		print('>>> Ran with Ver.%s'%self.StrVar99.get())
		if not int(float(self.StrVar99.get())) == int(float(version['gefit'])):
			print('>>> Version is not matached, try with version %s'%self.StrVar99.get())
			return False
		return True

	def reset_opt(self):

		self.initialise_vars()
		self.bdiff_index = find_bdiff_index(self,1)
		if float(self.MenuVar1.get()) > 0.:
			self.firstrun = False

		for i in range(1,14):

			self.__dict__['e%d'%i].delete(0,'end')
			self.__dict__['e%d'%i].insert(10,self.__dict__['StrVar%d'%i].get())

		for i in range(45,46):
			self.__dict__['e%d'%i].delete(0,'end')
			self.__dict__['e%d'%i].insert(10,self.__dict__['StrVar%d'%i].get())		

		self.textpad.delete('1.0',self.END)

		return

	def detect_close(self,ind):

		if ind == 4:
			try:	self.t5.destroy()
			except:	pass
			self.t5_close = True
	
		if ind == 6:
			try:	self.t14.destroy()
			except:	pass
			self.t14_close = True

		if ind == 8:
			self.efit_first_plot = True
			self.window_bnd = False
			self.window_coil = False
			self.window_pres = False
			self.window_j = False
			self.window_riter = False
			for i in [9,10,12,13,18]:	
				self.__dict__['t%i_close'%i] = True
				try:	self.__dict__['t%i'%i].destroy()
				except:	pass

		if ind == 9:	self.efit_first_plot = True
		if ind == 9:	self.window_bnd = False
		if ind == 10:	self.window_coil = False
		if ind == 12:	self.window_pres = False
		if ind == 13:	self.window_j = False
		if ind == 18:   self.window_riter = False

		
		self.__dict__['t%i_close'%ind] = True
		self.__dict__['t%i'%ind].destroy()

		return

	def scroll_func1(self):

		self.target_bnd	= self.DoubleVar1.get()

		if (self.iseqdsk):
			self.eq.target_psin = self.target_bnd
			self.eq.get_target_bnd()
		try:
			self.rzbdy = np.copy(self.eq.rzbdyt)
		except:
			print('>>> Have to load equilibrium first')

		update_psirz(self.iseqdsk,self.fig,self.eq,self.rzbdy,self)

		ind = np.argmax(self.rzbdy[:,0])
#		self.StrVar51.set(round(float(self.rzbdy[ind,0]),3))
#		self.StrVar52.set(round(float(self.rzbdy[ind,1]),3))
		temp = get_midRZ(self.eq.rzbdy,1.9)
		self.StrVar51.set('%6.4f'%temp)
		self.StrVar52.set('0.000')

		return

	def check_func1(self):		#Check option change of fitopt

		if self.change_fitopt:	self.change_fitopt = False
		else:	self.change_fitopt = True
		return

	def check_func2(self,ind):

		if ind == 1:	
			if self.CheckVar31.get() == 1:
				self.CheckVar32.set(0)
				self.CheckVar33.set(0)
				self.CheckVar34.set(0)
		if ind == 2:	
			if self.CheckVar32.get() == 1:
				self.CheckVar31.set(0)
				self.CheckVar33.set(0)
				self.CheckVar34.set(0)
		if ind == 3:	
			if self.CheckVar33.get() == 1:
				self.CheckVar31.set(0)
				self.CheckVar32.set(0)
				self.CheckVar34.set(0)				
		if ind == 4:	
			if self.CheckVar34.get() == 1:
				self.CheckVar31.set(0)
				self.CheckVar32.set(0)
				self.CheckVar33.set(0)	

		return

	def check_func3(self):

		if not self.ismse:	self.CheckVar3.set(1)
		return	

	def button_func1a(self):	#Make new case

		if self.run_index[-1] == '--':	self.ind = 0;
		else: self.ind = int(float(self.run_index[-1])) + 1

		try:	os.mkdir(str(self.ind))
		except:	pass
		try:	os.mkdir(str(self.ind)+'/INPUT')			#inputs
		except:	pass
		try:	os.mkdir(str(self.ind)+'/PROFILE')			#gfit or profiles
		except:	pass
		try:	os.mkdir(str(self.ind)+'/CHEASE')			#Chease solver 
		except:	pass
		try:	os.mkdir(str(self.ind)+'/EFIT')				#EFIT
		except:	pass		
		try:	os.mkdir(str(self.ind)+'/EFIT/RESULT')				#EFIT
		except:	pass								
		try:	os.mkdir('MDS')						#MDS save
		except:	pass	
		try:	os.mkdir(str(self.ind)+'/CSOLVE')			#current result
		except:	pass	
		try:	os.mkdir(str(self.ind)+'/NUBEAM')			#nubeam result
		except:	pass										

		if not (self.run_index[-1] == self.run_index2[-1]):
			self.run_index2.append(self.run_index[-1])
		self.run_index.append(str(self.ind))
		menu = self.m1["menu"]
		menu.delete(0, "end")
		self.MenuVar1.set(str(self.ind))
		for string in self.run_index:
			menu.add_command(label=string, command=lambda value=string: option_menu_update(self.MenuVar1,value,self.button_func1b,2))		

		menu = self.m2["menu"]
		menu.delete(0, "end")
		for string in self.run_index2:
			menu.add_command(label=string, command=lambda value=string: self.MenuVar2.set(value))	

		menu = self.m4["menu"]
		menu.delete(0, "end")
		for string in self.run_index2:
			menu.add_command(label=string, command=lambda value=string: self.MenuVar4.set(value))	

		if int(self.MenuVar1.get()) == 1: 
			self.MenuVar1.set(self.run_index[-1])
			self.MenuVar4.set('0')

		self.reset_opt()
		f= open(self.currdir+'/0/last_index','w')
		f.write('%s'%self.MenuVar1.get())
		f.close()	

		return

	def button_func1b(self,type=1):	#Load prev input x

		if self.firstrun:	return
		if self.MenuVar1.get() == '--':	return
		if (type == 1):	
			load_dir = self.currdir+'/%s/save.dat'%self.MenuVar2.get()
			load_dir2 = self.currdir+'/0/comment.dat'
		elif (type == 2 or type ==3):	
			load_dir = self.currdir+'/%s/save.dat'%self.MenuVar1.get()
			load_dir2 = self.currdir+'/0/comment.dat'
		verok = True
		try:
			verok = self.load_opt(load_dir)
			if verok:
				read_comment(self,load_dir2)
				print('>>> Parameters are loaded from %s'%load_dir)
		except:
			print('>>> No available save opt')
		if not verok: exit()
		if type==3:	
			f= open(self.currdir+'/0/last_index','w')
			f.write('%s'%self.MenuVar1.get())
			f.close()

		if os.path.isfile(self.e3.get()):
			self.button_func3a(True,self.e3.get())

		if os.path.isfile(self.e4.get()):
			self.button_func4a(True,True,self.e4.get())

		if os.path.isfile(self.e45.get()):
			self.button_func4c(True,True,self.e45.get())			

		if os.path.isfile(self.e5.get()):
			self.button_func5(True,self.e5.get(),False)		

		if ( not self.MenuVar1.get() == '0' and self.MenuVar4.get() == '--'): self.MenuVar4.set('0')
		
		self.bdiff_index = find_bdiff_index(self,1)
		self.update_status()
		return

	def button_func2(self):		#Run MDS

		mds_dir = self.currdir + '/MDS/'
		if not os.path.isdir(mds_dir): os.mkdir(mds_dir)
		mdsr_dir = self.currdir + '/MDS/result.dat'
		gtime,ktime,isefit = call_mse_data(mds_dir,int(self.e1.get()),int(self.e2.get()),self)

		if not isefit:	print('>>> There is no efit data!')

		if not float(self.e2.get()) == ktime:	
			print('>>> Target time is changed from %s to %i [ms]'%(self.e2.get(),ktime))
			self.e2.delete(0,'end')
			self.e2.insert(10,str(ktime))

		if not float(self.e2.get()) == gtime:
			print('>>> Equil time is changed from %s to %i [ms]'%(self.e2.get(),gtime))

		os.chdir(self.currdir)
		dat = read_mse_data(mdsr_dir)

		#dat[0] = mds_dir + 'g%06i.%06i'%(int(self.e1.get()),gtime)
		#dat[1] = mds_dir + 'k%06i.%06i'%(int(self.e1.get()),ktime)

		dat[0] = mds_dir + dat[0];
		dat[1] = mds_dir + dat[1];

		self.e3.delete(0,'end')
		try:	copyfile(dat[0],'%s/INPUT/%s'%(self.MenuVar1.get(),dat[0].split('/')[-1]))
		except:	print('>>> There is no proper gfile...');pass
		self.StrVar3.set('%s/INPUT/%s'%(self.MenuVar1.get(),dat[0].split('/')[-1]))
		self.e3.insert(10,self.StrVar3.get())

		self.e5.delete(0,'end')
		try:	copyfile(dat[1],'%s/INPUT/%s'%(self.MenuVar1.get(),dat[1].split('/')[-1]))
		except:	print('>>> There is no proper kfile...');pass
		self.StrVar5.set('%s/INPUT/%s'%(self.MenuVar1.get(),dat[1].split('/')[-1]))
		self.e5.insert(10,self.StrVar5.get())		

		self.wdia = float(dat[8])/1.e3
		wdia = str(round(self.wdia,1))
		self.e6.delete(0,'end')
		self.StrVar6.set(wdia)
		self.e6.insert(10,wdia)

		for i in range(9,21):
			if np.mod(i,2) == 1:	
				val = float(dat[i])
				if val>5.e5: val = val/1.e6
				elif val>5.3:val = val/1.e4
				else:        val = 0.
				if val < 0.1:	val = '0'
				else:	val = str(round(val,1))
				self.__dict__['StrVar%i'%(i+142)].set(val)

			else:	
				val = float(dat[i])
				if val < 10.:	val = '80'
				else:	val = str(round(val,1))
				self.__dict__['StrVar%i'%(i+142)].set(val)

		self.button_func3a(True,self.e3.get())
		self.button_func5(True,self.e5.get(),True)

		self.te_edge_file = 'MDS/'+dat[2].split('/')[-1]
		self.ne_edge_file = 'MDS/'+dat[3].split('/')[-1]
		self.te_file = 'MDS/'+dat[4].split('/')[-1]
		self.ne_file = 'MDS/'+dat[5].split('/')[-1]
		self.ti_file = 'MDS/'+dat[6].split('/')[-1]
		self.vt_file = 'MDS/'+dat[7].split('/')[-1]

		try:	copyfile(self.te_edge_file,'%s/INPUT/te_edge_file.dat'%self.MenuVar1.get())
		except: print('>>> te_edge_file no exists');pass
		try:	copyfile(self.ne_edge_file,'%s/INPUT/ne_edge_file.dat'%self.MenuVar1.get())
		except: print('>>> ne_edge_file no exists');pass
		try:    copyfile(self.ne_file,'%s/INPUT/ne_file.dat'%self.MenuVar1.get())
		except: print('>>> ne_file no exists');pass
		try:    copyfile(self.te_file,'%s/INPUT/te_file.dat'%self.MenuVar1.get())
		except: print('>>> te_file no exists');pass
		try:    copyfile(self.ti_file,'%s/INPUT/ti_file.dat'%self.MenuVar1.get())
		except: print('>>> ti_file no exists');pass
		try:    copyfile(self.vt_file,'%s/INPUT/vt_file.dat'%self.MenuVar1.get())
		except: print('>>> vt_file no exists');pass

		return

	def button_func3a(self,skipin=False,filename=None):	#Load Equ

		self.iseqdsk = False
		if not skipin:
			self.input = askopenfilename()
			if (len(self.input) == 0):	return
			try:	copyfile(self.input,'%s/INPUT/%s'%(self.MenuVar1.get(),self.input.split('/')[-1]))	
			except:	print('>>> Cannot copy gfile...');pass
			self.input = '%s/INPUT/%s'%(self.MenuVar1.get(),self.input.split('/')[-1])

		else:
			self.input = filename

		self.eq = eqdsk.eqdsk(self.input)
		self.eq.read_eqdsk_file()
		self.iseqdsk = True	
		draw_psirz(self.iseqdsk,self.fig,self.eq,self)
		self.scroll_func1()
		update_psirz(self.iseqdsk,self.fig,self.eq,self.rzbdy,self)

		self.rmag = self.eq.rmag
		self.q1 = self.eq.q[0]

		self.e3.delete(0, 'end')
		self.e3.insert(10,self.input)

		self.e9.configure(state='normal')
		self.e9.delete(0,'end')
		self.e9.insert(10,str(round(self.eq.wmhd/1.e3,1)))
		self.e9.configure(state='readonly')		

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
		self.update_status()				
		return

	def button_func3b(self):	#Check Equ

		if not self.iseqdsk:
			print('>> EQDSK must be selected first...')
			return
		if not self.t1_close:
			print('>>> EQU Param Window is already opened...')
			return			

		self.t1 = tk.Toplevel(self.root)
		self.t1.wm_title("EFIT Inform")
		self.t1_close = False		
		self.t1.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(1))

		self.t1.resizable(0,0)
		self.fig2, [ax1, ax2] = plt.subplots(1,2,figsize=(12,5))		

		self.canvas2 = FigureCanvasTkAgg(self.fig2,master=self.t1)
		self.plot_widget2 = self.canvas2.get_tk_widget()
		self.plot_widget2.grid(rowspan=40,row=2,column=0,columnspan=40)

		toolbar_frame = tk.Frame(self.t1)
		toolbar_frame.grid(column=0,row=0)
		
		self.toolbar2 = NavigationToolbar2Tk(self.canvas2,toolbar_frame)
		self.fig2.canvas.draw_idle()

		draw_eqdsk(ax1,ax2,self.eq)

		return

	def button_func4a(self,skip=False,skipin=False,filename=None):	#Load kinprof

		self.iskin = False
		
		if not skipin:
			input2 = askopenfilename()
			if (len(input2) == 0):	return
			try:	copyfile(input2,'%s/INPUT/%s'%(self.MenuVar1.get(),input2.split('/')[-1]))
			except:	print('>>> Cannot copy kinetic profile...');pass
			input2 = '%s/INPUT/%s'%(self.MenuVar1.get(),input2.split('/')[-1])
		else:
			input2 = filename

		self.button_func4b(skip,input2)

		self.iskin = True
		self.e4.delete(0,'end')
		self.e4.insert(10,input2)	
		self.update_status()
		return

	def button_func4b(self,skip=False,ifile = None):	#Check kinprof

		if not self.t2_close:
			print('>>> Kin profile Window is already opened...')
			return			

		if ifile == None:	filename = self.e4.get()
		else:	filename = ifile

		
		if filename.strip() == '':	
			print('>>> No kinetic files...')	
			return

		if not skip:
			self.t2 = tk.Toplevel(self.root)
			self.t2.wm_title("PROFILE Inform")
			self.t2_close = False		
			self.t2.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(2))

			self.t2.resizable(0,0)
			self.fig3, [ax1, ax2, ax3] = plt.subplots(1,3,figsize=(13,5))		

			self.canvas3 = FigureCanvasTkAgg(self.fig3,master=self.t2)
			self.plot_widget3 = self.canvas3.get_tk_widget()
			self.plot_widget3.grid(rowspan=40,row=2,column=0,columnspan=40)

			toolbar_frame = tk.Frame(self.t2)
			toolbar_frame.grid(column=0,row=0)
			
			self.toolbar3 = NavigationToolbar2Tk(self.canvas3,toolbar_frame)
			self.fig3.canvas.draw_idle()
		else:
			self.fig3 = None;	ax1 = None; ax2 = None; ax3 = None;

		if self.iseqdsk:
			self.zeff,self.zimp,self.amain,self.aimp,self.lineden2,self.wkin = draw_kin_file(filename,self.fig3,ax1,ax2,ax3,self.eq,skip)
		else:
			self.zeff,self.zimp,self.amain,self.aimp,self.lineden2,self.wkin = draw_kin_file(filename,self.fig3,ax1,ax2,ax3,None,skip)

		self.e7.delete(0,'end')
		self.e7.insert(10,str(self.zeff))

		self.StrVar14.set(self.trans_vars(self.zimp,2))

		self.e10.configure(state='normal')
		self.e10.delete(0,'end')
		self.e10.insert(10,str(round(self.wkin/1.e3,1)))
		self.e10.configure(state='readonly')	

		self.e11.configure(state='normal')
		self.e11.delete(0,'end')
		self.e11.insert(10,str(round(self.lineden2,2)))
		self.e11.configure(state='readonly')				
		return

	def button_func4c(self,skip=False,skipin=False,filename=None):	#Load rot_file

		self.isrfile = False
		
		if not skipin:
			input2 = askopenfilename()
			if (len(input2) == 0):	return
			try:	copyfile(input2,'%s/INPUT/%s'%(self.MenuVar1.get(),input2.split('/')[-1]))
			except:	print('>>> Cannot copy rot file...');pass
			input2 = '%s/INPUT/%s'%(self.MenuVar1.get(),input2.split('/')[-1])
		else:
			input2 = filename

		self.button_func4d(skip,input2)

		self.isrfile = True
		self.e45.delete(0,'end')
		self.e45.insert(10,input2)	
		self.update_status()
		return

	def button_func4d(self,skip=False,ifile=None):	#Check rot file

		if not self.t15_close:
			print('>>> Rot profile Window is already opened...')
			return			
	
		if ifile == None:	filename = self.e45.get()
		else:	filename = ifile

		if filename.strip() == '':
			print('>>> No rotation file...')
			return

		if not skip:
			self.t15 = tk.Toplevel(self.root)
			self.t15.wm_title("PROFILE Inform")
			self.t15_close = False		
			self.t15.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(15))

			self.t15.resizable(0,0)
			self.fig16, ax1 = plt.subplots(1,1,figsize=(5,5))		

			self.canvas16 = FigureCanvasTkAgg(self.fig16,master=self.t15)
			self.plot_widget16 = self.canvas16.get_tk_widget()
			self.plot_widget16.grid(rowspan=40,row=2,column=0,columnspan=40)

			toolbar_frame = tk.Frame(self.t15)
			toolbar_frame.grid(column=0,row=0)
			
			self.toolbar16 = NavigationToolbar2Tk(self.canvas16,toolbar_frame)
			self.fig16.canvas.draw_idle()
		else:
			return

		draw_rot_file(filename,self.fig16,ax1,skip)

		return

	def button_func5(self,skip=False,filename=None,draw_plot=True):	#Read kfile

		self.iskfile_already = False
		if not skip:
			self.iskfile = False
			input2 = askopenfilename()
			if (len(input2) == 0):	return	
			try:	copyfile(input2,'%s/INPUT/%s'%(self.MenuVar1.get(),input2.split('/')[-1]))	
			except:	print('>>> Cannot copy kfile...');pass
			input2 = '%s/INPUT/%s'%(self.MenuVar1.get(),input2.split('/')[-1])
			self.CheckVar23.set(0)
		else:
			input2 = filename
		
		self.e5.delete(0,'end')
		self.e5.insert(10,input2)
		try:
			read_kfile(self,input2)
			load_kfile(self,False)
			self.iskfile = True
		except:
			self.iskfile = False
			self.iskfile_already = False
			self.CheckVar23.set(0)
	
		if not self.ismse:
			mse_dat = _make_mse(int(float(self.e1.get())),float(float(self.e2.get()))/1000.,0.05)
			self.ismse = mse_dat['is']
			if self.ismse:
				self.tgamma = np.copy(mse_dat['tgam'])
				self.sgamma = np.copy(mse_dat['sgam'])
				self.fwtgam = np.copy(mse_dat['fwt'])
				self.rrrgam = np.copy(mse_dat['RRRGAM'])
				self.zzzgam = np.copy(mse_dat['ZZZGAM'])
				self.aa1gam = np.copy(mse_dat['AA1GAM'])
				self.aa2gam = np.copy(mse_dat['AA1GAM'])
				self.aa3gam = np.copy(mse_dat['AA1GAM'])
				self.aa4gam = np.copy(mse_dat['AA1GAM'])
				self.aa5gam = np.copy(mse_dat['AA1GAM'])
				self.aa6gam = np.copy(mse_dat['AA1GAM'])
				for yy in years:					
					if int(float(self.e1.get())) in shotk[yy]['shot']:
						self.fwtgam = np.copy(mse_good_ch[yy]); print('>>> Discharge Year',yy)
				for i in range(len(self.fwtgam)):
					if self.rrrgam[i] < 0: 
						self.fwtgam[i] = 0.
						self.sgamma[i] = 0.;
		self.dtgamma = np.copy(self.sgamma)
		for i in range(61,61+len(self.sgamma)):
			if self.__dict__['StrVar%i'%i].get() == '': self.__dict__['StrVar%i'%i].set('0.0') 
			self.dtgamma[i-61] = float(self.__dict__['StrVar%i'%i].get())

		if not self.ismse:	self.MenuVar7.set('sMSE')
		if not self.ismse:
			print('>>> No MSE data in K-FILE, please use sMSE option!')
			self.use_smse = True
			self.CheckVar3.set(1)

		self.iskfile_already = True
		self.CheckVar23.set(1)
		self.update_status()
		return


	def button_func5a(self):

		if (not self.ismse or not self.iskfile):
			print('>>> No MSE data')
			return

		if not self.t17_close:
			print('>>> K-FILE window is already opened...')
			return

		self.t17 = tk.Toplevel(self.root)
		self.t17.wm_title("RAW MSE Profile")
		self.t17_close = False
		self.t17.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(17))

		self.t17.resizable(0,0)
		self.fig14, self.eax18 = plt.subplots(1,1,figsize=(5.,5.))
		
		self.canvas17 = FigureCanvasTkAgg(self.fig14,master=self.t17)
		self.plot_widget17 = self.canvas17.get_tk_widget()
		self.plot_widget17.grid(row=1,column=1)
		
		toolbar_frame = tk.Frame(self.t17)
		toolbar_frame.grid(column=1,row=0)

		self.toolbar17 = NavigationToolbar2Tk(self.canvas17,toolbar_frame)
		self.fig14.canvas.draw_idle()

		if float(self.StrVar151.get()) > 0.:   mseb = 'NB1A'; ind = 152
		elif float(self.StrVar153.get()) > 0.: mseb = 'NB1B'; ind = 154
		else: mseb = 'NB1C'; ind=156
		sshot = int(float(self.e1.get()))

		if ((sshot >= 20756) and (sshot < 21759)): 
			mseb = 'NB1B'; ind = 154
		if ((sshot >= 23057) and (sshot < 23134)): 
			mseb = 'NB1C'; ind = 156
		if (sshot >= 25283):
			if self.aa1gam[0] < 0.83: mseb = 'NB1A'; ind = 152
			else: mseb = 'NB1B'; ind = 154
		print('>>> Beam used in MSE', mseb)
		draw_efit_mse(self,self.eax18,mse_rad_err[mseb]/1.e3,mseb)

		return

	def button_func6(self):		#Run GFIT

		if not self.iseqdsk:
			print('>>> Eq file is not selected ..')
			return
		
		try:	os.mkdir('%s/PROFILE'%self.MenuVar1.get())
		except:	print('>>> Error in making PROFILE dir'); pass

		if (self.MenuVar4.get() == '--' and int(self.MenuVar1.get()) > 0.):	return

		self.ind = int(float(self.MenuVar1.get()))

		gfit_dir = self.currdir + '/%i/PROFILE/'%self.ind
		gfit_inp = self.currdir + '/%s/PROFILE/gfitp_history.dat'%self.ind
		prof_dir = self.currdir + '/%i/PROFILE/PROFILES/chease_kinprof.out'%self.ind
		vt_dir   = self.currdir + '/%s/PROFILE/PROFILES/VT_fit.dat'%self.ind

		isrefl = 0; isece = 0;
		if self.ind > 0:
			dirs = self.currdir + '/%s/INPUT/refl_ece.dat'%self.MenuVar4.get()
			if os.path.isfile(dirs):
				f = open(dirs,'r')
				isrefl = int(float(f.readline()))
				isece  = int(float(f.readline()))
				f.close()

		os.chdir(gfit_dir)
		if not os.path.isfile(gfit_inp):
			os.system("echo 'zeff %s ' >> gfitp_history.dat"%self.e7.get())			
			os.system("echo 'zimp %s ' >> gfitp_history.dat"%self.e14.get())	
			os.system("echo 'lineden %s ' >> gfitp_history.dat"%self.e8.get())			
			os.system("echo 'eq_file ../../%s ' >> gfitp_history.dat"%self.e3.get())			
			os.system("echo 'bs_model %s ' >> gfitp_history.dat"%self.MenuVar3.get().lower())
			os.system("echo 'bsmulti %s ' >> gfitp_history.dat"%self.e13.get())
			if not self.te_edge_file == '':
				os.system("echo 'te_edge_file ../../%s ' >> gfitp_history.dat"%self.te_edge_file)
			if not self.ne_edge_file == '':
				os.system("echo 'ne_edge_file ../../%s ' >> gfitp_history.dat"%self.ne_edge_file)
			if not self.te_file == '':
				os.system("echo 'te_file ../../%s ' >> gfitp_history.dat"%self.te_file)
			if not self.ne_file == '':
				os.system("echo 'ne_file ../../%s ' >> gfitp_history.dat"%self.ne_file)
			if not self.ti_file == '':
				os.system("echo 'ti_file ../../%s ' >> gfitp_history.dat"%self.ti_file)
			if not self.vt_file == '':
				os.system("echo 'vt_file ../../%s ' >> gfitp_history.dat"%self.vt_file)

			os.system("echo 'shot %i ' >> gfitp_history.dat"%int(float(self.e1.get())))
			os.system("echo 'time %i ' >> gfitp_history.dat"%int(float(self.e2.get())))
			os.system("echo 'mds_dir ../../../MDS/ ' >> gfitp_history.dat")

		if self.ind > 0.:
			try:    copyfile('../../%s/INPUT/ne_file.dat'%(self.MenuVar4.get()),'../../%s/INPUT/ne_file.dat'%(self.MenuVar1.get()))
			except: print('>>> Error in copying ne_file'); pass
			try:    copyfile('../../%s/INPUT/ne_edge_file.dat'%(self.MenuVar4.get()),'../../%s/INPUT/ne_edge_file.dat'%(self.MenuVar1.get()))
			except: print('>>> Error in copying ne_edge_file'); pass
			try:    copyfile('../../%s/INPUT/te_file.dat'%(self.MenuVar4.get()),'../../%s/INPUT/te_file.dat'%(self.MenuVar1.get()))
			except: print('>>> Error in copying te_file'); pass
			try:    copyfile('../../%s/INPUT/te_edge_file.dat'%(self.MenuVar4.get()),'../../%s/INPUT/te_edge_file.dat'%(self.MenuVar1.get()))
			except: print('>>> Error in copying te_edge_file'); pass
			try:    copyfile('../../%s/INPUT/ti_file.dat'%(self.MenuVar4.get()),'../../%s/INPUT/ti_file.dat'%(self.MenuVar1.get()))
			except: print('>>> Error in copying ti_file'); pass
			try:    copyfile('../../%s/INPUT/vt_file.dat'%(self.MenuVar4.get()),'../../%s/INPUT/vt_file.dat'%(self.MenuVar1.get()))
			except: print('>>> Error in copying vt_file'); pass
			try:	os.mkdir(self.currdir + '/%i/PROFILE/FITPROF'%(self.ind))
			except:	print('>>> Error in making FITPROF dir');pass
			if isrefl == 1:
				try:    copyfile('../../%s/INPUT/refl_file.dat'%(self.MenuVar4.get()),'../../%s/INPUT/refl_file.dat'%(self.MenuVar1.get()))
				except: print('>>> Error in copying refl_file');pass
			if isece == 1:
				try:    copyfile('../../%s/INPUT/ece_file.dat'%(self.MenuVar4.get()),'../../%s/INPUT/ece_file.dat'%(self.MenuVar1.get()))
				except: print('>>> Error in copying ece_file');pass

			gfit_opt_old1 = self.currdir + '/%s/PROFILE/FITPROF/fit_opt.save'%(self.MenuVar4.get())
			gfit_opt_old2 = self.currdir + '/%s/PROFILE/FITPROF/fit_opt.save_param'%(self.MenuVar4.get())
			gfit_opt1 = self.currdir + '/%i/PROFILE/FITPROF/fit_opt.save'%(self.ind)
			gfit_opt2 = self.currdir + '/%i/PROFILE/FITPROF/fit_opt.save_param'%(self.ind)

			if not os.path.isfile(gfit_opt1): self.change_fitopt = True

			if self.change_fitopt:
				self.change_fitopt = False	
				try: os.remove(gfit_opt2)
				except:	print('>>> No pre-existing fit_opt');pass
				try:	copyfile(gfit_opt_old1,gfit_opt1)
				except:	print('>>> Cannot copy pre-existing fit_opt');pass
				if self.CheckVar1.get() == 0:
					try:	copyfile(gfit_opt_old2,gfit_opt2)
					except:	print('>>> Cannot copy initial fit_opt');pass
				try: copytree(self.currdir+'/%s/PROFILE/FITPROF/PROFILES'%self.MenuVar4.get(),self.currdir+'/%i/PROFILE/FITPROF/PROFILES'%self.ind)
				except:	print('>>> Cannot copy pre-existing fittings');pass

			ece=''; refl='';
			if self.CheckVar1.get() == 1:
				te = '%s/PROFILE/PROFILES/TE_pre.dat'%(self.MenuVar4.get())
				ne = '%s/PROFILE/PROFILES/NE_pre.dat'%(self.MenuVar4.get())
				ti = '%s/PROFILE/PROFILES/TI_pre.dat'%(self.MenuVar4.get())
				vt = '%s/PROFILE/PROFILES/VT_pre.dat'%(self.MenuVar4.get())
				modify_fit_opt(gfit_opt1,2,self.e3.get(),te,ne,None,None,ti,vt)
				modify_gfitp_opt(self,gfit_inp,te,ne,'','',ti,vt)

			else:
				te = '%s/INPUT/te_file.dat'%self.MenuVar4.get()
				ne = '%s/INPUT/ne_file.dat'%self.MenuVar4.get()
				te_edge = '%s/INPUT/te_edge_file.dat'%self.MenuVar4.get()
				ne_edge = '%s/INPUT/ne_edge_file.dat'%self.MenuVar4.get()
				ti = '%s/INPUT/ti_file.dat'%self.MenuVar4.get()
				vt = '%s/INPUT/vt_file.dat'%self.MenuVar4.get()
				if isrefl==1: refl = '%s/INPUT/refl_file.dat'%self.MenuVar4.get()
				if isece ==1: ece  = '%s/INPUT/ece_file.dat'%self.MenuVar4.get()
				modify_fit_opt(gfit_opt1,1,self.e3.get(),te,ne,te_edge,ne_edge,ti,vt,refl,ece)
				modify_gfitp_opt(self,gfit_inp,te,ne,te_edge,ne_edge,ti,vt,refl,ece)

		else:
			modify_gfitp_opt(self,gfit_inp,self.te_file,self.ne_file,self.te_edge_file,self.ne_edge_file,self.ti_file,self.vt_file)			
		
		if not self.ts: os.system(gfitp_dir)
		else:	os.system(gfitp_dir+'  -plare')
		os.chdir(self.currdir)

		if os.path.isfile(prof_dir):	
			try:	copyfile(prof_dir,'%s/INPUT/%s'%(self.MenuVar1.get(),prof_dir.split('/')[-1]))
			except:	print('>>> Error in copying fitted profiles');pass
			self.e4.delete(0,'end')
			self.e4.insert(10,'%s/INPUT/%s'%(self.MenuVar1.get(),prof_dir.split('/')[-1]))
			self.button_func4a(True,True,self.e4.get())
		else:
			print('>>> No kinetic profiles..')

		if os.path.isfile(vt_dir):	
			try:	copyfile(vt_dir,'%s/INPUT/%s'%(self.MenuVar1.get(),vt_dir.split('/')[-1]))
			except:	print('>>> Error in copying fitted rotation profile');pass
			self.e45.delete(0,'end')
			self.e45.insert(10,'%s/INPUT/%s'%(self.MenuVar1.get(),vt_dir.split('/')[-1]))
			self.button_func4c(True,True,self.e45.get())
		else:
			print('>>> No rotation profiles..')			

		f = open(gfit_inp,'r')
		ece_file = None; refl_file = None;
		while True:
			line = f.readline()
			if not line: break
			line = line.split()
			if line[0] == 'refl_file': refl_file = line[1]
			if line[0] == 'ece_file' : ece_file = line[1]
		f.close()

		if self.ind==0 :dirs1 = self.currdir + '/%i/PROFILE/'%self.ind
		else: dirs1 = self.currdir + '/%i/PROFILE/'%self.ind
		dirs2 = self.currdir + '/%i/INPUT/'%self.ind

		if not refl_file==None: copyfile(dirs1+refl_file,dirs2+'refl_file.dat')
		if not ece_file==None: copyfile(dirs1+ece_file,dirs2+'ece_file.dat')
		f = open(dirs2+'/refl_ece.dat','w')
		if refl_file == None: f.write('0\n')
		else: f.write('1\n')
		if ece_file == None: f.write('0\n')
		else: f.write('1\n')
		f.close()

		self.update_status()

		return

	def button_func6a(self):

		self.nogfitp = False

		self.button_func6b()	
		if not self.nogfitp:	return	

		input2 = askdirectory()
		if len(input2) == 0:	
			return

		dir1 = '%s/PROFILE'%self.MenuVar1.get()

		os.system('ln -s %s %s'%(input2,dir1))

		dir2 = dir1 + '/gfitp_history.dat'
		read_gfit_opt(self,dir2)
		try:	copyfile(dir1+'/'+self.ne_file,'%s/INPUT/ne_file.dat'%self.MenuVar1.get())
		except:	print('>>> Error in copying ne_file');pass
		try:	copyfile(dir1+'/'+self.te_file,'%s/INPUT/te_file.dat'%self.MenuVar1.get())
		except: print('>>> Error in copying te_file');pass
		try:    copyfile(dir1+'/'+self.te_edge_file,'%s/INPUT/te_edge_file.dat'%self.MenuVar1.get())
		except: print('>>> Error in copying te_edge_file');pass
		try:    copyfile(dir1+'/'+self.ne_edge_file,'%s/INPUT/ne_edge_file.dat'%self.MenuVar1.get())
		except: print('>>> Error in copying ne_edge_file');pass
		try:	copyfile(dir1+'/'+self.ti_file,'%s/INPUT/ti_file.dat'%self.MenuVar1.get())
		except: print('>>> Error in copying ti_file');pass
		try:	copyfile(dir1+'/'+self.vt_file,'%s/INPUT/vt_file.dat'%self.MenuVar1.get())
		except: print('>>> Error in copying vt_file');pass
		
		os.system('diff %s %s'%(self.currdir+'/'+self.e3.get(),dir1+'/'+self.eq_file))

		mds_dir = '%s/MDS'%input2
		if os.path.isdir(mds_dir):
			try:	rmtree('MDS')
			except:	pass
			try:	os.remove('MDS')
			except:	pass
			os.system('ln -s %s %s'%(mds_dir,'MDS'))	

		return

	def button_func6b(self):
		
		self.t11 = tk.Toplevel(self.root)
		self.t11.wm_title('DELETE FILE')
		self.t11.resizable(0,0)

		dirs = '%s/PROFILE'%self.MenuVar1.get()
		
		self.l2 = tk.Label(self.t11, text = ' Do you want to remove previous GFIT #%s ? '%self.MenuVar1.get(),justify='center')
		self.l2.grid(row=0,column=0,columnspan=4,pady=15)
		
		self.b1 = tk.Button(self.t11, text="YES", bg = "lightgray",command=lambda: self.button_func6c(True),height = 1,width = 5)
		self.b1.grid(row=1, column=1, columnspan=1)
		self.b1 = tk.Button(self.t11, text="NO", bg = "lightgray",command=lambda: self.button_func6c(False),height = 1,width = 5)
		self.b1.grid(row=1, column=2, columnspan=1)

		self.t11.mainloop()

		return

	def button_func6c(self,delete):

		if not delete:	
			self.t11.quit()
			self.t11.destroy()
			return

		try:   rmtree('%s/PROFILE'%self.MenuVar1.get())
		except: pass
		try:   os.remove('%s/PROFILE'%self.MenuVar1.get())
		except: pass

		self.nogfitp = delete

		self.t11.quit()
		self.t11.destroy()

		return

	def button_func7aa(self,val):

		input2 = askopenfilename()
		if (len(input2) == 0):	return	
		try:	copyfile(input2,'%s/INPUT/%s'%(self.MenuVar1.get(),input2.split('/')[-1]))
		except:	pass
		self.iskfile = True
		val.delete(0,'end')
		val.insert(10,'%s/INPUT/%s'%(self.MenuVar1.get(),input2.split('/')[-1]))
		return

	def button_func7ab(self):

		for i in range(101,122):
			self.__dict__['StrVar%d'%(i)].set(self.__dict__['e%d'%(i)].get())

		for i in range(151,163):
			self.__dict__['StrVar%d'%(i)].set(self.__dict__['e%d'%(i)].get())			

		self.bpower, self.benergy, self.nbeam = make_nubeam_config(self)

		print('>>> BPOWER',self.bpower)
		print('>>> BENERGY',self.benergy)
		print('>>> NBEAM',self.nbeam)

		self.t3_close = True
		self.t3.destroy()

		return

	def button_func7a(self):	#RUN PARAMS

		if not self.t3_close:
			print('>>> Run Param Window is already opened...')
			return			

		self.t3 = tk.Toplevel(self.root)
		self.t3.wm_title("Run Params")
		self.t3_close = False		
		self.t3.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(3))

		self.t3.resizable(0,0)

		self.l1 = tk.Label(self.t3, text="======= CHEASE =======",justify='center')
		self.l1.grid(row=0, column=0,columnspan=4)
		self.l1 = tk.Label(self.t3, text="------- Current -------",justify='center')
		self.l1.grid(row=1, column=0,columnspan=4)

		label = ['Core_neo [0-1]','VLOOP_MOD [0-1]','VLOOP_EXT [0-1]','HAG_MOD_PSIN [0-1]','HAG_CORE_MOD',]
		for i in range(4):

			self.l1 = tk.Label(self.t3, text=label[i],justify='left')
			self.l1.grid(row=(i+2), column=0,columnspan=2)
			self.__dict__['e%d'%(i+101)] = tk.Entry(self.t3,width=6,justify='center')
			self.__dict__['e%d'%(i+101)].insert(10,self.__dict__['StrVar%d'%(i+101)].get())
			self.__dict__['e%d'%(i+101)].grid(row=(i+2), column=2,columnspan=2)	

		self.l1 = tk.Label(self.t3, text=label[4],justify='left')
		self.l1.grid(row=7, column=0,columnspan=2)
		self.c1 = tk.Checkbutton(self.t3,variable=self.CheckVar21)
		self.c1.grid(row=7, column=3)		

		self.l1 = tk.Label(self.t3, text="------ Numeric ------",justify='center')
		self.l1.grid(row=8, column=0,columnspan=4)

		label = ['ITERC [#]','ITERB [#]','RELAX [0-1]','NS [#]','NT [#]','MAP_NS [#]','MAP_NT[#]','Ip_crit','Bs_crit']
		for i in range(9):

			self.l1 = tk.Label(self.t3, text=label[i],justify='left')
			self.l1.grid(row=(i+9), column=0,columnspan=2)
			self.__dict__['e%d'%(i+105)] = tk.Entry(self.t3,width=6,justify='center')
			self.__dict__['e%d'%(i+105)].insert(10,self.__dict__['StrVar%d'%(i+105)].get())
			self.__dict__['e%d'%(i+105)].grid(row=(i+9), column=2,columnspan=2)		

		self.l1 = tk.Label(self.t3, text="======= NUBEAM =======",justify='center')
		self.l1.grid(row=0, column=4,columnspan=4)
		self.l1 = tk.Label(self.t3, text="-------- NBI1 --------",justify='center')
		self.l1.grid(row=1, column=4,columnspan=4)
		self.l1 = tk.Label(self.t3, text="-------- NBI2 --------",justify='center')
		self.l1.grid(row=5, column=4,columnspan=4)			
		self.l1 = tk.Label(self.t3, text="  Power[MW]",justify='center')
		self.l1.grid(row=3, column=4)
		self.l1 = tk.Label(self.t3, text="  Power[MW]",justify='center')
		self.l1.grid(row=6, column=4)		
		self.l1 = tk.Label(self.t3, text="  Energy[keV]",justify='center')
		self.l1.grid(row=4, column=4)
		self.l1 = tk.Label(self.t3, text="  Energy[keV]",justify='center')
		self.l1.grid(row=7, column=4)		
		label = ['A','B','C']
		for i in range(3):
			self.l1 = tk.Label(self.t3, text=label[i],justify='center')
			self.l1.grid(row=2, column=(i+5))
		count = 151
		for i in range(3):
			for j in range(2):				
					self.__dict__['e%d'%(count)] = tk.Entry(self.t3,width=4,justify='center')
					self.__dict__['e%d'%(count)].insert(10,self.__dict__['StrVar%d'%(count)].get())
					self.__dict__['e%d'%(count)].grid(row=(j+3), column=(i+5),columnspan=1)		
					self.__dict__['e%d'%(count+6)] = tk.Entry(self.t3,width=4,justify='center')
					self.__dict__['e%d'%(count+6)].insert(10,self.__dict__['StrVar%d'%(count+6)].get())
					self.__dict__['e%d'%(count+6)].grid(row=(j+6), column=(i+5),columnspan=1)									
					count = count + 1

		self.l1 = tk.Label(self.t3, text="------ Numeric ------",justify='center')
		self.l1.grid(row=8, column=4,columnspan=4)

		label = ['RUN STEP [#]','RUN AVG [#]','RUN DT [s]','AVG DT [s]','NPROC [#]','MFILE','SFILE','IFILE']
		for i in range(5):

			self.l1 = tk.Label(self.t3, text=label[i],justify='left')
			self.l1.grid(row=(i+9), column=4,columnspan=2)
			self.__dict__['e%d'%(i+114)] = tk.Entry(self.t3,width=6,justify='center')
			self.__dict__['e%d'%(i+114)].insert(10,self.__dict__['StrVar%d'%(i+114)].get())
			self.__dict__['e%d'%(i+114)].grid(row=(i+9), column=6,columnspan=2)

		self.l1 = tk.Label(self.t3, text=label[5],justify='left')
		self.l1.grid(row=14, column=4,columnspan=2)
		b1 = tk.Button(self.t3, text="OPEN", bg = "lightgray",command=lambda: self.button_func7aa(self.e119),height = 1,width = 4)
		b1.grid(row=14, column=5,columnspan=2)
		self.e119 = tk.Entry(self.t3,width=25,justify='center')
		self.e119.insert(10,self.StrVar120.get())
		self.e119.grid(row=15,column=4,columnspan=4)

		self.l1 = tk.Label(self.t3, text=label[6],justify='left')
		self.l1.grid(row=16, column=4,columnspan=2)
		b1 = tk.Button(self.t3, text="OPEN", bg = "lightgray",command=lambda: self.button_func7aa(self.e120),height =1,width = 4)
		b1.grid(row=16, column=5,columnspan=2)
		self.e120 = tk.Entry(self.t3,width=25,justify='center')
		self.e120.insert(10,self.StrVar120.get())
		self.e120.grid(row=17,column=4,columnspan=4)		

		self.l1 = tk.Label(self.t3, text=label[7],justify='left')
		self.l1.grid(row=18, column=4,columnspan=2)
		b1 = tk.Button(self.t3, text="OPEN", bg = "lightgray",command=lambda: self.button_func7aa(self.e121),height =1,width = 4)
		b1.grid(row=18, column=5,columnspan=2)
		self.e121 = tk.Entry(self.t3,width=25,justify='center')
		self.e121.insert(10,self.StrVar121.get())
		self.e121.grid(row=19,column=4,columnspan=4)		

		self.l1 = tk.Label(self.t3, text='USE_BEAM',justify='left')
		self.l1.grid(row=20, column=4,columnspan=2)
		self.c1 = tk.Checkbutton(self.t3,variable=self.CheckVar22)
		self.c1.grid(row=20, column=6)			

		b1 = tk.Button(self.t3, text="SAVE", bg = "lightgray",command=lambda: self.button_func7ab(),height =1,width = 4)
		b1.grid(row=21, column=6,columnspan=2)		

		return		

	def button_func7b(self):	#RUN CHEASE

		if not (self.iseqdsk):	return
		if not (self.iskin):	return
		if not (self.isrfile):	return

		self.bpower, self.benergy, self.nbeam = make_nubeam_config(self)

		self.bdiff_index =  find_bdiff_index(self,1)
		bdiff = round(float(self.e12.get()),2)
		len2 = len(self.bdiff_index)
		for i in range(len2):
			if self.bdiff_index[i] == bdiff:	break

		if len2==0.:	i = 0;
		else:
			if self.bdiff_index[i] == bdiff:
				print('>>> There is previous run! Choose different Beam diffusivity')
				print('>>> Diffusivities ->',self.bdiff_index)
				self.quit = False;
				self.button_func7ba(bdiff)
				if self.quit:	return

		self.bdiff_index =  find_bdiff_index(self,1)

		if self.nbeam == 0.:	
			self.StrVar106.set(1)
			self.CheckVar22.set(0)

		self.button_func10a()

		if not self.ismse:
			if self.MenuVar5.get().lower() == 'emse':
				print('>>> No MSE data, use sMSE')
				return

		if self.MenuVar5.get().lower() == 'smse':
			run_smse(self,bdiff,chease_dir)

		if self.MenuVar5.get().lower() == 'emse':
			run_emse(self,bdiff,nubeam_dir2,python2_exec)

		self.update_status()
		return

	def button_func7ba(self,bdiff):

		self.t11 = tk.Toplevel(self.root)
		self.t11.wm_title('DELETE FILE')
		self.t11.resizable(0,0)

		if	self.MenuVar5.get().lower() == 'smse':	dirs = self.currdir + '/%s/CSOLVE/Bdiff_%03i'%(self.MenuVar1.get(),round(bdiff*100.,0))
		else:	dirs = self.currdir + '/%s/NUBEAM/Bdiff_%03i'%(self.MenuVar1.get(),round(bdiff*100.,0))
		
		self.l2 = tk.Label(self.t11, text = ' Do you want to remove previous run for B-Diff = %s ?    '%bdiff,justify='center')
		self.l2.grid(row=0,column=0,columnspan=4,pady=15)

		self.b1 = tk.Button(self.t11, text="YES", bg = "lightgray",command=lambda: self.button_func7bc(True,dirs),height = 1,width = 5)
		self.b1.grid(row=1, column=1, columnspan=1)	
		self.b1 = tk.Button(self.t11, text="NO", bg = "lightgray",command=lambda: self.button_func7bc(False,dirs),height = 1,width = 5)
		self.b1.grid(row=1, column=2, columnspan=1)	

		self.t11.mainloop()

		return

	def button_func7bc(self,flag,dirs):

		if flag:
			try:	rmtree(dirs)
			except:	pass
		else:
			self.quit = True
		self.t11.quit()	
		self.t11.destroy()
		return

	def button_func7c(self):	#OPEN CHEASE

		if not self.t4_close:
			print('>>> CSOLVE Window is already opened...')
			return			

		self.bdiff_index = find_bdiff_index(self,1)
		if len(self.bdiff_index) == 0.:
			return			

		self.t4 = tk.Toplevel(self.root)
		self.t4.wm_title("CSOLVE Inform #%s %s[ms]"%(self.e1.get(),self.e2.get()))
		self.t4_close = False		
		self.t4.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(4))

		self.t4.resizable(0,0)
		self.fig4, ([ax1, ax2], [ax3, ax4]) = plt.subplots(2,2,figsize=(13,7))		

		self.canvas4 = FigureCanvasTkAgg(self.fig4,master=self.t4)
		self.plot_widget4 = self.canvas4.get_tk_widget()
		self.plot_widget4.grid(rowspan=30,row=2,column=10,columnspan=40)

		toolbar_frame = tk.Frame(self.t4)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar4 = NavigationToolbar2Tk(self.canvas4,toolbar_frame)
		self.fig4.canvas.draw_idle()



		list1 = [];	list2 = [];
		for i in range(int(float(self.StrVar106.get()))):
			list1.append(str(i+1))

		for i in range(len(self.bdiff_index)):
			list2.append(str(round(float(self.bdiff_index[i]),2)))			

		ind = int(float(list1[-1]))

		self.l2 = tk.Label(self.t4, text = '========= Diffusivity scan =========',justify='center')
		self.l2.grid(row=1,column=0,columnspan=8)

		self.MenuVar11.set(list1[-1])
		self.MenuVar12.set(list2[0])
		
		if self.MenuVar5.get().lower() == 'smse':
			self.l1 = tk.Label(self.t4, text=" BEAM ITER [#]",justify='center')
			self.l1.grid(row=3, column=0,columnspan=3,sticky='e')
			self.m11 = tk.OptionMenu(self.t4,self.MenuVar11,*list1, command=lambda value: self.button_func7ca())
			self.m11.config(width=3)
			self.m11.grid(row=3, column=4,columnspan=2,sticky='e')

			self.l1 = tk.Label(self.t4, text=" ITER CHECK [#]",justify='center')
			self.l1.grid(row=2, column=0,columnspan=3,sticky='e')
			self.m12 = tk.OptionMenu(self.t4,self.MenuVar12,*list2)
			self.m12.config(width=3)
			self.m12.grid(row=2, column=4,columnspan=2,sticky='e')
			b1 = tk.Button(self.t4, text="PLOT", bg = "lightgray",command=lambda: self.button_func7cb(),height = 1,width = 4)
			b1.grid(row=2, column=6,columnspan=2,sticky='w')		

		count = 4;
		if self.MenuVar5.get().lower() == 'smse':
			self.l1 = tk.Label(self.t4, text=" BPOL Crit ",justify='center')
			self.l1.grid(row=count, column=0,columnspan=3,sticky='e')
			self.e201 = tk.Entry(self.t4,width=5,justify='center')
			self.e201.insert(0,'0')
			self.e201.grid(row=count, column=4,columnspan=2)
			self.e201.configure(state='readonly')
			self.c1 = tk.Checkbutton(self.t4,variable=self.CheckVar31,command=lambda: self.check_func2(1))
			self.c1.grid(row=count, column=6)
			count = count + 1

		self.l1 = tk.Label(self.t4, text=" WMHD Crit ",justify='center')
		self.l1.grid(row=count, column=0,columnspan=3,sticky='e')
		self.e202 = tk.Entry(self.t4,width=5,justify='center')
		self.e202.insert(0,'0')
		self.e202.grid(row=count, column=4,columnspan=2)
		self.e202.configure(state='readonly')
		self.c1 = tk.Checkbutton(self.t4,variable=self.CheckVar32,command=lambda: self.check_func2(2))
		self.c1.grid(row=count, column=6)
		count = count + 1		

		self.l1 = tk.Label(self.t4, text=" WDIA Crit ",justify='center')
		self.l1.grid(row=count, column=0,columnspan=3,sticky='e')
		self.e203 = tk.Entry(self.t4,width=5,justify='center')
		self.e203.insert(0,'0')
		self.e203.grid(row=count, column=4,columnspan=2)
		self.e203.configure(state='readonly')
		self.c1 = tk.Checkbutton(self.t4,variable=self.CheckVar33,command=lambda: self.check_func2(3))
		self.c1.grid(row=count, column=6)		
		count = count + 1		

		if self.MenuVar5.get().lower() == 'smse':
			self.l1 = tk.Label(self.t4, text=" RMAG Crit ",justify='center')
			self.l1.grid(row=count, column=0,columnspan=3,sticky='e')
			self.e204 = tk.Entry(self.t4,width=5,justify='center')
			self.e204.insert(0,'0')
			self.e204.grid(row=count, column=4,columnspan=2)
			self.e204.configure(state='readonly')
			self.c1 = tk.Checkbutton(self.t4,variable=self.CheckVar34,command=lambda: self.check_func2(4))
			self.c1.grid(row=count, column=6)			

		b1 = tk.Button(self.t4, text="DONE", bg = "lightgray",command=lambda: self.button_func7cc(),height = 1,width = 8)
		b1.grid(row=9, column=0,columnspan=8)			

		if self.MenuVar5.get().lower() == 'smse':
			read_smse_result(self,self.bdiff_index,self.fig4,ind,1)
		elif self.MenuVar5.get().lower() == 'emse':
			read_emse_result(self,self.bdiff_index,self.fig4)
		self.update_status()
		return

	def button_func7ca(self):

		ind = int(float(self.MenuVar11.get()))
		read_smse_result(self,self.bdiff_index,self.fig4,ind,1)
		self.update_status()
		return

	def button_func7cb(self):

		self.t5 = tk.Toplevel(self.t4)
		self.t5.wm_title("Iteration Inform #%s %s[ms]"%(self.e1.get(),self.e2.get()))

		self.t5.resizable(0,0)
		self.fig5, ([ax1, ax2, ax3], [ax4, ax5, ax6]) = plt.subplots(2,3,figsize=(13,7))		

		self.canvas5 = FigureCanvasTkAgg(self.fig5,master=self.t5)
		self.plot_widget5 = self.canvas5.get_tk_widget()
		self.plot_widget5.grid(rowspan=30,row=2,column=0,columnspan=40)

		toolbar_frame = tk.Frame(self.t5)
		toolbar_frame.grid(column=0,row=0)
		
		self.toolbar5 = NavigationToolbar2Tk(self.canvas5,toolbar_frame)
		self.fig5.canvas.draw_idle()

		bdiff = float(self.MenuVar12.get());
		bdiff = round(bdiff,2)
		for i in range(len(self.bdiff_index)):
			if self.bdiff_index[i] == bdiff:	break

		self.bdiff_index = find_bdiff_index(self,1)
		read_smse_result(self,self.bdiff_index,self.fig5,i,2)
		self.update_status()
		return

	def button_func7cc(self):

		self.t4_close = True
		line = None
		if self.CheckVar31.get() == 1:	line = self.e201.get()
		if self.CheckVar32.get() == 1:	line = self.e202.get()
		if self.CheckVar33.get() == 1:	line = self.e203.get()
		if self.CheckVar34.get() == 1:	line = self.e204.get()

		if not line == None:
			self.e12.delete(0,'end')
			self.e12.insert(10,line)

		self.ischease = True
		self.t4.destroy()
		self.update_status()
		return

	def button_func8a(self):
		self.CheckVar3.set(1)
		if (self.ismse and self.MenuVar7.get().lower() == 'emse'):	self.CheckVar3.set(0)
		if not self.ismse:	self.CheckVar3.set(1)
		
		if not self.iskfile:
			print('>>> KFILE is not selected...')
			return	
		if not self.t6_close:
			print('>>> EFIT constraint is already opened...')
			return			

		if not self.t8_close:
			print('>>> EFIT run window is already opened...')
			return	

		self.bdiff_index = find_bdiff_index(self,2)
		if len(self.bdiff_index)==0.:	
			print('>>> No chease run')
			self.ischease = False
			return

		self.t6 = tk.Toplevel(self.root)
		self.t6.wm_title("EFIT constraint Inform #%s %s[ms]"%(self.e1.get(),self.e2.get()))
		self.t6_close = False		
		self.t6.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(6))

		self.t6.resizable(0,0)
		self.fig6, ([ax1, ax2], [ax3, ax4]) = plt.subplots(2,2,figsize=(13,7))		

		self.canvas6 = FigureCanvasTkAgg(self.fig6,master=self.t6)
		self.plot_widget6 = self.canvas6.get_tk_widget()
		self.plot_widget6.grid(rowspan=30,row=2,column=10,columnspan=40)

		toolbar_frame = tk.Frame(self.t6)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar6 = NavigationToolbar2Tk(self.canvas6,toolbar_frame)
		self.fig6.canvas.draw_idle()

		self.l2 = tk.Label(self.t6, text = '========= Singal constraints =========',justify='center')
		self.l2.grid(row=1,column=0,columnspan=10)		

		b1 = tk.Button(self.t6, text="Select Signals", bg = "lightgray",command=lambda: self.button_func8b(),height = 1,width = 15)
		b1.grid(row=2, column=0,columnspan=10)

		self.l1 = tk.Label(self.t6, text="MSE",justify='center')
		self.l1.grid(row=3, column=0,columnspan=2)
#		self.e16 = tk.Entry(self.t6,width=7,justify='center')
#		self.e16.insert(10,self.StrVar16.get())
#		self.e16.grid(row=3, column=4,columnspan=2,sticky='w')	

		b1 = tk.Button(self.t6, text="SGAM", bg = "lightgray",command=lambda: self.button_func11(1),height = 1,width = 5)
		b1.grid(row=3,column=2,columnspan=2)

		b1 = tk.Button(self.t6, text="DTGAM", bg = "lightgray",command=lambda: self.button_func12(1),height = 1,width = 5)
		b1.grid(row=3,column=4,columnspan=2)

		self.l1 = tk.Label(self.t6, text="sMSE",justify='center')
		self.l1.grid(row=3, column=6,columnspan=3,sticky='e')
		self.c1 = tk.Checkbutton(self.t6,variable=self.CheckVar3,command=lambda: self.check_func3())
		self.c1.grid(row=3, column=9)		
		if (self.ismse and self.MenuVar7.get().lower() == 'emse'): self.c1.config(state='disabled')

		self.l2 = tk.Label(self.t6, text = '========= Knots constraint =========',justify='center')
		self.l2.grid(row=4,column=0,columnspan=10)	
		self.l1 = tk.Label(self.t6, text="Pressure",justify='center')
		self.l1.grid(row=5, column=0,columnspan=3)
		self.e17 = tk.Entry(self.t6,width=30,justify='center')
		self.e17.insert(10,self.StrVar17.get())
		self.e17.grid(row=6, column=0,columnspan=10)	

		self.l1 = tk.Label(self.t6, text="FFprime",justify='center')
		self.l1.grid(row=7, column=0,columnspan=3)
		self.e18 = tk.Entry(self.t6,width=30,justify='center')
		self.e18.insert(10,self.StrVar18.get())
		self.e18.grid(row=8, column=0,columnspan=10)

		self.l1 = tk.Label(self.t6, text="Current",justify='center')
		self.l1.grid(row=9, column=0,columnspan=3)
		self.e19 = tk.Entry(self.t6,width=30,justify='center')
		self.e19.insert(10,self.StrVar19.get())
		self.e19.grid(row=10, column=0,columnspan=10)

		b1 = tk.Button(self.t6, text="Plot", bg = "lightgray",command=lambda: draw_pp_constraint(self),height = 1,width = 4)
		b1.grid(row=5, column=4,columnspan=2)	
		b1 = tk.Button(self.t6, text="Plot", bg = "lightgray",command=lambda: draw_ff_constraint(self),height = 1,width = 4)
		b1.grid(row=7, column=4,columnspan=2)	
		b1 = tk.Button(self.t6, text="Plot", bg = "lightgray",command=lambda: draw_j_constraint(self),height = 1,width = 4)
		b1.grid(row=9, column=4,columnspan=2)	

		b1 = tk.Button(self.t6, text="FIND KNOTS", bg = "lightgray",command=lambda: self.button_func8c(1),height = 1,width = 8)
		b1.grid(row=5, column=6,columnspan=4)	
		b1 = tk.Button(self.t6, text="FIND KNOTS", bg = "lightgray",command=lambda: self.button_func8c(2),height = 1,width = 8)
		b1.grid(row=7, column=6,columnspan=4)	
		b1 = tk.Button(self.t6, text="MAKE PNTS", bg = "lightgray",command=lambda: self.button_func8d(),height = 1,width = 8)
		b1.grid(row=9, column=6,columnspan=4)			

		self.l1 = tk.Label(self.t6, text="USE_JCONST",justify='center')
		self.l1.grid(row=11, column=0,columnspan=4)
		self.c1 = tk.Checkbutton(self.t6,variable=self.CheckVar4)
		self.c1.grid(row=11, column=5)	

		self.l2 = tk.Label(self.t6, text = '=========   EFIT options   =========',justify='center')
		self.l2.grid(row=12,column=0,columnspan=10)	

		self.l1 = tk.Label(self.t6, text="MXITER [1-999]",justify='center')
		self.l1.grid(row=13, column=0,columnspan=4)
		self.e38 = tk.Entry(self.t6,width=6,justify='center')
		self.e38.insert(10,self.StrVar38.get())
		self.e38.grid(row=13, column=4,columnspan=3)

		self.l1 = tk.Label(self.t6, text="RELAX [0-1]",justify='center')
		self.l1.grid(row=14, column=0,columnspan=4)
		self.e39 = tk.Entry(self.t6,width=6,justify='center')
		self.e39.insert(10,self.StrVar39.get())
		self.e39.grid(row=14, column=4,columnspan=3)

		self.l1 = tk.Label(self.t6, text="CONVERG [0-1]",justify='center')
		self.l1.grid(row=15, column=0,columnspan=4)
		self.e40 = tk.Entry(self.t6,width=6,justify='center')
		self.e40.insert(10,self.StrVar40.get())
		self.e40.grid(row=15, column=4,columnspan=3)

		self.l1 = tk.Label(self.t6, text="FWTCUR [0-10]",justify='center')
		self.l1.grid(row=16, column=0,columnspan=4)
		self.e41 = tk.Entry(self.t6,width=6,justify='center')
		self.e41.insert(10,self.StrVar41.get())
		self.e41.grid(row=16, column=4,columnspan=3)

		self.l1 = tk.Label(self.t6, text="QCONST",justify='center')
		self.l1.grid(row=17, column=0,columnspan=4)
		self.e42 = tk.Entry(self.t6,width=6,justify='center')
		self.e42.insert(10,self.StrVar42.get())
		self.e42.grid(row=17, column=4,columnspan=3)
		self.c1 = tk.Checkbutton(self.t6,variable=self.CheckVar5)
		self.c1.grid(row=17, column=7)

		self.l1 = tk.Label(self.t6, text="SUPP-P0",justify='center')
		self.l1.grid(row=19, column=0,columnspan=4)
		self.e43 = tk.Entry(self.t6,width=6,justify='center')
		self.e43.insert(10,self.StrVar43.get())
		self.e43.grid(row=19, column=4,columnspan=3)		

		self.l1 = tk.Label(self.t6, text="SHIFT_MSE [cm]",justify='center')
		self.l1.grid(row=20, column=0,columnspan=4)
		self.e44 = tk.Entry(self.t6,width=6,justify='center')
		self.e44.insert(10,self.StrVar44.get())
		self.e44.grid(row=20, column=4,columnspan=3)
		b1 = tk.Button(self.t6, text="Plot", bg = "lightgray",command=lambda: draw_mse_constraint(self),height = 1,width = 4)
		b1.grid(row=20, column=7,columnspan=2)

		self.l1 = tk.Label(self.t6, text="",justify='center')
		self.l1.grid(row=21, column=0,columnspan=4)		

		b1 = tk.Button(self.t6, text="GENERATE INPUT", bg = "lightgray",command=lambda: self.button_func8e(),height = 1,width = 15)
		b1.grid(row=22, column=0,columnspan=10)

		if not self.MenuVar8.get().lower() == 'hmode':
			self.use_jconst = True
			self.StrVar17.set('0.0,0.5,1.0')
			self.StrVar18.set('0.0,0.6,1.0')
			self.StrVar19.set('0.80.85,0.9,0.95,1.0')
			self.CheckVar4.set(1)
		else:
			self.use_jconst = True
			self.CheckVar4.set(1)

		print('>>> Generate Field information....')
		get_efit_constraint_eq(self)

		print('>>> Done!')
		draw_efit_constraint(self,self.fig6)

		self.l1 = tk.Label(self.t6, text="Q0[REF/MOD]",justify='center')
		self.l1.grid(row=18, column=0,columnspan=4)
		self.e46 = tk.Entry(self.t6,width=6,justify='center')
		self.e46.insert(10,str(round(self.q1,3)))
		self.e46.grid(row=18, column=4,columnspan=3)
		self.e46.configure(state='readonly')

		self.e47 = tk.Entry(self.t6,width=6,justify='center')
		self.e47.insert(10,str(round(self.q2,3)))
		self.e47.grid(row=18, column=7,columnspan=3)
		self.e47.configure(state='readonly')		

		return

	def button_func8b(self,type=1):

		if not self.t7_close:
			print('>>> Signal Window is already opened...')
			return			

		self.t7 = tk.Toplevel(self.root)
		self.t7.wm_title("Probe signal Inform")
		self.t7_close = False		
		self.t7.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(7))

		self.t7.resizable(0,0)
		
		self.l2 = tk.Label(self.t7, text = '===========  Mangetic coil ===========',justify='center')
		self.l2.grid(row=0,column=0,columnspan=10)
		self.l2 = tk.Label(self.t7, text = '=========== PSI LOOP coil ============',justify='center')
		self.l2.grid(row=0,column=11,columnspan=10)		
		self.l2 = tk.Label(self.t7, text = '=======  MSE =======',justify='center')
		self.l2.grid(row=0,column=22,columnspan=5)				

		i = 1; j = 0; k = 1; 
		while k <= self.clen1:
			self.__dict__['cc%d'%k] = tk.Checkbutton(self.t7,variable=self.__dict__['CoilVar%d'%k],bg='lime')
			self.__dict__['cc%d'%k].grid(row=i, column=j)
			j = j + 1
			k = k + 1
			if j==10:
				i = i +1
				j = 0		

		i = 1; j = 0; k = 1; 
		while k <= self.clen2:
			self.__dict__['cc%d'%(k+self.clen1)] = tk.Checkbutton(self.t7,variable=self.__dict__['CoilVar%d'%(k+150)],bg='magenta')
			self.__dict__['cc%d'%(k+self.clen1)].grid(row=i, column=(j+11))
			j = j + 1
			k = k + 1
			if j==10:
				i = i +1
				j = 0	

		i = 1; j = 0; k = 1; 
		while k <= self.clen3:
			if type == 1:
				self.__dict__['cc%d'%(k+self.clen1+self.clen2)] = tk.Checkbutton(self.t7,variable=self.__dict__['CoilVar%d'%(k+200)],bg='yellow',command=lambda: draw_mse_constraint(self))
				if (self.rrrgam[k-1] < 0.): 
					self.__dict__['CoilVar%d'%(k+200)].set(0)
					self.__dict__['cc%d'%(k+self.clen1+self.clen2)].config(state='disabled')
				else: self.__dict__['cc%d'%(k+self.clen1+self.clen2)].config(state='normal')
			else:
				self.__dict__['cc%d'%(k+self.clen1+self.clen2)] = tk.Checkbutton(self.t7,variable=self.__dict__['CoilVar%d'%(k+200)],bg='yellow')
				if (self.rrrgam[k-1] < 0.): 
					self.__dict__['CoilVar%d'%(k+200)].set(0)
					self.__dict__['cc%d'%(k+self.clen1+self.clen2)].config(state='disabled')
				else: self.__dict__['cc%d'%(k+self.clen1+self.clen2)].config(state='normal')
			self.__dict__['cc%d'%(k+self.clen1+self.clen2)].grid(row=i, column=(j+22))
			j = j + 1
			k = k + 1
			if j==5:
				i = i +1
				j = 0				

		b1 = tk.Button(self.t7, text="CLOSE", bg = "lightgray",command=lambda: self.button_func8ba(),height = 1,width = 5)
		b1.grid(row=11, column=23,columnspan=4, rowspan=2)								

		return

	def button_func8ba(self):
		self.t7_close = True
		self.t7.destroy()
		return

	def button_func8c(self,type):

		if not self.t14_close:
			print('>>> Knots finder is already opened...')
			return			

		if self.MenuVar8.get().lower() == 'lmode':
			if float(self.StrVar21.get()) == 0.985: self.StrVar21.set('0.6')
			if float(self.StrVar46.get()) == 0.1:   self.StrVar46.set('0.2')
			if float(self.StrVar28.get()) == 0.985: self.StrVar28.set('0.6')
			if float(self.StrVar48.get()) == 0.1:   self.StrVar48.set('0.2')

		self.t14 = tk.Toplevel(self.t6)
		self.t14.wm_title("KNOTS FINDER OPT")
		self.t14_close = False		
		self.t14.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(14))

		self.t14.resizable(0,0)

		self.l1 = tk.Label(self.t14, text="=========== KNOTS FINDER ===========",justify='center')
		self.l1.grid(row=0, column=0,columnspan=10)

		if type == 1:	count = 20;	count2 = 46
		elif type == 2:	count = 27;	count2 = 48

		labels = ['SPLINE START [0-1]','END KNOT [0-1]','START KNOT [0-1]','KNOT SHIFT [0-1]','CORE KNOT [#]','EDGE KNOT [#]','SPLINE POINT [#]','DXMIN CORE [0-1]','DXMIN EDGE [0-1]']

		self.l1 = tk.Label(self.t14, text="=========== KNOTS FINDER ===========",justify='center')
		self.l1.grid(row=0, column=0,columnspan=10)

		for i in range(count,count+7):

			self.l1 = tk.Label(self.t14, text=labels[i-count],justify='center')
			self.l1.grid(row=(i-count+1), column=0,columnspan=5)
			self.__dict__["e%i"%i] = tk.Entry(self.t14,width=7,justify='center')
			self.__dict__["e%i"%i].insert(10,self.__dict__["StrVar%i"%i].get())
			self.__dict__["e%i"%i].grid(row=(i-count+1), column=6,columnspan=3)			

		for i in range(count2,count2+2):

			self.l1 = tk.Label(self.t14, text=labels[i-count2+7],justify='center')
			self.l1.grid(row=(i-count2+8), column=0,columnspan=5)
			self.__dict__["e%i"%(i+2)] = tk.Entry(self.t14,width=7,justify='center')
			self.__dict__["e%i"%(i+2)].insert(10,self.__dict__["StrVar%i"%i].get())
			self.__dict__["e%i"%(i+2)].grid(row=(i-count2+8), column=6,columnspan=3)

		b1 = tk.Button(self.t14, text="FIND", bg = "lightgray",command=lambda: self.button_func8ca(type),height = 1,width = 5)
		b1.grid(row=11, column=0,columnspan=4)	

		b1 = tk.Button(self.t14, text="CLOSE", bg = "lightgray",command=lambda: self.button_func8cb(),height = 1,width = 5)
		b1.grid(row=11, column=5,columnspan=4)			
		return

	def button_func8ca(self,type):

		if type == 1:	count = 20;	count2 = 7;	count3 = 46;
		elif type == 2:	count = 27;	count2 = 7;	count3 = 48;
		elif type == 3: count = 34; 	count2 = 4;

		for i in range(count,count+count2):	self.__dict__["StrVar%i"%i].set(self.__dict__["e%i"%i].get())
		if (type == 1 or type == 2):	
			for i in range(count3,count3+2):	self.__dict__["StrVar%i"%i].set(self.__dict__["e%i"%(i+2)].get())
		
		if type == 3:	sknot = make_current_knots(self,int(float(self.e34.get())),int(float(self.e35.get())),
			float(self.e36.get()),float(self.e37.get()))
		else:	sknot = make_pf_knots(self,type)

		if type == 1:
			self.e17.delete(0,'end')
			self.e17.insert(10,sknot)
			draw_pp_constraint(self)

		if type == 2:
			self.e18.delete(0,'end')
			try:	self.e18.insert(10,sknot)
			except:	print('>>> Check error message..')
			draw_ff_constraint(self)

		if type == 3:
			self.e19.delete(0,'end')
			self.e19.insert(10,sknot)
			draw_j_constraint(self)		

		return

	def button_func8cb(self):
		self.t14_close = True
		self.t14.destroy()
		return

	def button_func8d(self):

		if not self.t14_close:
			print('>>> Knots finder is already opened...')
			return			

		self.t14 = tk.Toplevel(self.t6)
		self.t14.wm_title("KNOTS FINDER OPT")
		self.t14_close = False		
		self.t14.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(14))

		self.t14.resizable(0,0)

		self.l1 = tk.Label(self.t14, text="=========== KNOTS FINDER ===========",justify='center')
		self.l1.grid(row=0, column=0,columnspan=10)


		count = 34

		labels = ['CORE KNOT [#]','EDGE KNOT [#]','START KNOT [0-1]','END KNOT [0-1]']

		self.l1 = tk.Label(self.t14, text="=========== KNOTS FINDER ===========",justify='center')
		self.l1.grid(row=0, column=0,columnspan=10)

		for i in range(count,count+4):

			self.l1 = tk.Label(self.t14, text=labels[i-count],justify='center')
			self.l1.grid(row=(i-count+1), column=0,columnspan=5)
			self.__dict__["e%i"%i] = tk.Entry(self.t14,width=7,justify='center')
			self.__dict__["e%i"%i].insert(10,self.__dict__["StrVar%i"%i].get())
			self.__dict__["e%i"%i].grid(row=(i-count+1), column=6,columnspan=3)			

		b1 = tk.Button(self.t14, text="FIND", bg = "lightgray",command=lambda: self.button_func8ca(3),height = 1,width = 5)
		b1.grid(row=9, column=0,columnspan=4)	

		b1 = tk.Button(self.t14, text="CLOSE", bg = "lightgray",command=lambda: self.button_func8cb(),height = 1,width = 5)
		b1.grid(row=9, column=5,columnspan=4)		

		return

	def button_func8e(self):

		if self.is_sgam <= 0.:
			print('>>> SGAM configuration is not ready!')
			return

		efitdir = self.currdir +'/%s/EFIT/'%(self.MenuVar1.get())
		shot = int(float(self.e1.get()));	time = int(float(self.e2.get()));
		kfile_dir = efitdir + 'k%06i.%06i_kin'%(shot,time)

		self.StrVar17.set(self.e17.get())
		self.StrVar18.set(self.e18.get())
		self.StrVar19.set(self.e19.get())

		self.DoubleVar2.set(self.q1)
		self.DoubleVar3.set(self.q2)

		for i in range(38,45):
			self.__dict__['StrVar%i'%i].set(self.__dict__['e%i'%i].get())

		make_kfile(self)
		write_efit_input(self,kfile_dir)

		bdiff = int(round(100*float(self.MenuVar6.get()),0))
		if self.MenuVar7.get().lower() == 'emse': copyfile(self.e4.get(),efitdir+'p%06i.%06i'%(shot,time))
		else:
			pdir = self.currdir + '/%s/CSOLVE/Bdiff_%03i/PROFILES/chease_kinprof_new'%(self.MenuVar1.get(),bdiff)
			copyfile(pdir,efitdir+'p%06i.%06i'%(shot,time))

		if self.CheckVar22.get() == 1:

			if self.MenuVar7.get().lower() == 'emse':
				pf_dir = self.currdir + '/%s/NUBEAM/Bdiff_%03i/chease_pres'%(self.MenuVar1.get(),bdiff)
				jf_dir = self.currdir + '/%s/NUBEAM/Bdiff_%03i/chease_curr'%(self.MenuVar1.get(),bdiff)
			else:
				index = int(self.MenuVar9.get())
				pf_dir = self.currdir + '/%s/CSOLVE/Bdiff_%03i/PROFILES/chease_pres_%i'%(self.MenuVar1.get(),bdiff,index)
				jf_dir = self.currdir + '/%s/CSOLVE/Bdiff_%03i/PROFILES/chease_curr_%i'%(self.MenuVar1.get(),bdiff,index)
	
			copyfile(pf_dir,efitdir+'/chease_pres')
			copyfile(jf_dir,efitdir+'/chease_curr')	
			
		make_chease_opt(self,efitdir+'/chease_opt')
		print('>>> K-FILE is save to %s'%kfile_dir)
		self.button_func10a()
		self.t6_close = True
		self.t6.destroy()

		return

	def button_func9a(self):

		if not self.t8_close:
			print('>>> EFIT run is already opened...')
			return		
		if not self.t6_close:
			print('>>> EFIT constraint is already opened...')
			return
	
		self.window_coil = False
		self.window_bnd  = False
		self.window_pres = False
		self.window_j = False
		self.window_riter = False

		self.MenuVar14.set(self.MenuVar1.get())
		efitdir = self.currdir +'/%s/EFIT/'%(self.MenuVar1.get())
		shot = int(float(self.e1.get()));	time = int(float(self.e2.get()));
		kfile_dir = efitdir + 'k%06i.%06i_kin'%(shot,time)

		init_efit = efitdir+'RESULT/g%06i.%06i_kin_0'%(shot,time)

		init_pfile= efitdir+'p%06i.%06i'%(shot,time)
		init_pf   = efitdir+'chease_pres'
		init_jf   = efitdir+'chease_curr'
		if not os.path.isfile(init_pfile):
			print('>>> No PFILE, copy the reference...')
			if self.MenuVar7.get().lower() == 'emse': 
				copyfile(self.e4.get(),init_pfile)
			else:
				bdiff = int(round(100*float(self.MenuVar6.get()),0))
				pdir = self.currdir + '/%s/CSOLVE/Bdiff_%03i/PROFILES/chease_kinprof_new'%(self.MenuVar1.get(),bdiff)
				copyfile(pdir,init_pfile)
		bdiff = int(round(100*float(self.MenuVar6.get()),0))
		if self.CheckVar22.get() == 1:
			if self.MenuVar7.get().lower() == 'emse':
				pf_dir = self.currdir + '/%s/NUBEAM/Bdiff_%03i/chease_pres'%(self.MenuVar1.get(),bdiff)
				jf_dir = self.currdir + '/%s/NUBEAM/Bdiff_%03i/chease_curr'%(self.MenuVar1.get(),bdiff)
			else:
				index = int(self.MenuVar9.get())
				pf_dir = self.currdir + '/%s/CSOLVE/Bdiff_%03i/PROFILES/chease_pres_%i'%(self.MenuVar1.get(),bdiff,index)
				jf_dir = self.currdir + '/%s/CSOLVE/Bdiff_%03i/PROFILES/chease_curr_%i'%(self.MenuVar1.get(),bdiff,index)

			if not os.path.isfile(init_pf): copyfile(pf_dir,efitdir+'/chease_pres')
			if not os.path.isfile(init_jf): copyfile(jf_dir,efitdir+'/chease_curr')
		
		if not os.path.isfile(kfile_dir):
			print('>>> No KFILE....')
			return

		try:		
			a=np.copy(self.tgamma2)
			get_efit_constraint_eq(self,True)
		except:	get_efit_constraint_eq(self,False)

		currdir = os.getcwd()
		if (int(self.MenuVar1.get()) == 0 and not os.path.isfile(init_efit)):	
			write_efit_input(self,efitdir+'kfile_run',0)
			efit_exec = efitdir +'efit129' #257
			efit_exec2 = efit_dir+'/efit129'
			if not os.path.isfile(efit_exec):       
				try:	copyfile(efit_exec2,efit_exec)
				except:	
					print('>>> Error in copying efit exe')
					return
			os.system('chmod 777 %s'%efit_exec)

			os.chdir(efitdir)	
			subprocess.run([efit_exec])
			move('kfile_run','RESULT/k%06i.%06i_kin_0'%(shot,time))
			try:
				copyfile(currdir+'/'+self.e3.get(),'RESULT/g%06i.%06i_kin_0'%(shot,time))
				copyfile('m%06i.%06i'%(shot,time),'RESULT/m%06i.%06i_kin_0'%(shot,time))
				copyfile(init_pfile,'RESULT/p%06i.%06i_kin_0'%(shot,time))
				copyfile('g%06i.%06i'%(shot,time),'RESULT/g%06i.%06i_kin_00'%(shot,time))
				copyfile('fitout.dat','RESULT/fitout.dat_0')
			except:
				pass
			os.remove('g%06i.%06i'%(shot,time))
			self.button_func9ba(0,True)	
			w_dir2 = efitdir + 'RESULT/wmhd.dat_0'
			mse_dir2 = efitdir + 'RESULT/mse.dat_0'
			pres_dir2 = efitdir + 'RESULT/pres.dat_0'
			j_dir2 = efitdir + 'RESULT/jconst.dat_0'
			map_dir2  = efitdir + 'RESULT/map.dat_0'
			efit_post_process(self,0,currdir+'/'+self.e3.get(),w_dir2,mse_dir2,pres_dir2,j_dir2,map_dir2,True)

		os.chdir(currdir)
		self.efit_index = find_efit_runs(self)
		if self.MenuVar10.get() == '--':	self.MenuVar10.set(self.efit_index[-1])
		else:
			if self.efit_index[-1] == '--':
				self.MenuVar10.set('--')
			else:
				if float(self.MenuVar10.get()) > float(self.efit_index[-1]):
					self.MenuVar10.set(self.efit_index[-1])

		self.efit_index = find_efit_runs(self)

		self.t8 = tk.Toplevel(self.root)
		self.t8.wm_title("EFIT RUN #%s %s[ms]-%s"%(self.e1.get(),self.e2.get(),self.MenuVar1.get()))
		self.t8_close = False		
		self.t8.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(8))

		self.t8.resizable(0,0)
		self.fig8, ([self.eax1, self.eax2], [self.eax3, self.eax4]) = plt.subplots(2,2,figsize=(13,7))		

		self.canvas8 = FigureCanvasTkAgg(self.fig8,master=self.t8)
		self.plot_widget8 = self.canvas8.get_tk_widget()
		self.plot_widget8.grid(rowspan=30,row=2,column=10,columnspan=40)

		toolbar_frame = tk.Frame(self.t8)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar8 = NavigationToolbar2Tk(self.canvas8,toolbar_frame)
		self.fig8.canvas.draw_idle()

		self.l2 = tk.Label(self.t8, text = '========= Singal constraints =========',justify='center')
		self.l2.grid(row=1,column=0,columnspan=10)		

		b1 = tk.Button(self.t8, text="Select Signals", bg = "lightgray",command=lambda: self.button_func8b(2),height = 1,width = 15)
		b1.grid(row=2, column=0,columnspan=10)

		self.l1 = tk.Label(self.t8, text="MSE",justify='center')
		self.l1.grid(row=3, column=0,columnspan=2)

		b1 = tk.Button(self.t8, text="SGAM",  bg = "lightgray",command=lambda: self.button_func11(),height = 1,width = 5)
		b1.grid(row=3,column=2,columnspan=2)

#		b1 = tk.Button(self.t8, text="DTGAM", bg = "lightgray",command=lambda: self.button_func12(),height = 1,width = 5)
#		b1.grid(row=3,column=4,columnspan=2)

#		self.e16 = tk.Entry(self.t8,width=7,justify='center')
#		self.e16.insert(10,self.StrVar16.get())
#		self.e16.grid(row=3, column=4,columnspan=2,sticky='w')	
		
		self.l1 = tk.Label(self.t8, text="sMSE",justify='center')
		self.l1.grid(row=3, column=6,columnspan=2,sticky='e')
		self.c1 = tk.Checkbutton(self.t8,variable=self.CheckVar3,command=lambda: self.check_func3())
		self.c1.grid(row=3, column=8)
		if (self.ismse and self.MenuVar7.get().lower() == 'emse'): self.c1.config(state='disabled')

		self.l2 = tk.Label(self.t8, text = '========= Knots constraint =========',justify='center')
		self.l2.grid(row=4,column=0,columnspan=10)	
		self.l1 = tk.Label(self.t8, text="Pressure",justify='center')
		self.l1.grid(row=5, column=0,columnspan=3)
		self.e17 = tk.Entry(self.t8,width=30,justify='center')
		self.e17.insert(10,self.StrVar17.get())
		self.e17.grid(row=6, column=0,columnspan=10)	

		self.l1 = tk.Label(self.t8, text="FFprime",justify='center')
		self.l1.grid(row=7, column=0,columnspan=3)
		self.e18 = tk.Entry(self.t8,width=30,justify='center')
		self.e18.insert(10,self.StrVar18.get())
		self.e18.grid(row=8, column=0,columnspan=10)

		self.l2 = tk.Label(self.t8, text = '=========   EFIT options   =========',justify='center')
		self.l2.grid(row=12,column=0,columnspan=10)	

		self.l1 = tk.Label(self.t8, text="MXITER [1-999]",justify='center')
		self.l1.grid(row=13, column=0,columnspan=4)
		self.e38 = tk.Entry(self.t8,width=6,justify='center')
		self.e38.insert(10,self.StrVar38.get())
		self.e38.grid(row=13, column=4,columnspan=3)

		self.l1 = tk.Label(self.t8, text="RELAX [0-1]",justify='center')
		self.l1.grid(row=14, column=0,columnspan=4)
		self.e39 = tk.Entry(self.t8,width=6,justify='center')
		self.e39.insert(10,self.StrVar39.get())
		self.e39.grid(row=14, column=4,columnspan=3)

		self.l1 = tk.Label(self.t8, text="CONVERG [0-1]",justify='center')
		self.l1.grid(row=15, column=0,columnspan=4)
		self.e40 = tk.Entry(self.t8,width=6,justify='center')
		self.e40.insert(10,self.StrVar40.get())
		self.e40.grid(row=15, column=4,columnspan=3)

		self.l1 = tk.Label(self.t8, text="FWTCUR [0-10]",justify='center')
		self.l1.grid(row=16, column=0,columnspan=4)
		self.e41 = tk.Entry(self.t8,width=6,justify='center')
		self.e41.insert(10,self.StrVar41.get())
		self.e41.grid(row=16, column=4,columnspan=3)

		self.l1 = tk.Label(self.t8, text="QCONST",justify='center')
		self.l1.grid(row=17, column=0,columnspan=4)
		self.e42 = tk.Entry(self.t8,width=6,justify='center')
		self.e42.insert(10,self.StrVar42.get())
		self.e42.grid(row=17, column=4,columnspan=3)
		self.c1 = tk.Checkbutton(self.t8,variable=self.CheckVar5)
		self.c1.grid(row=17, column=7)

		if self.ts:
			self.l1 = tk.Label(self.t8, text="Er-Correct [#]",justify='center')
			self.l1.grid(row=18, column=0,columnspan=4)
			self.c1 = tk.Checkbutton(self.t8,variable=self.CheckVar35)
			self.c1.grid(row=18, column=7)
			self.m15 = tk.OptionMenu(self.t8,self.MenuVar15,'1','2','3','4')
			self.m15.config(width=2)
			self.m15.grid(row=18, column=4,columnspan=3)
			self.l1 = tk.Label(self.t8, text="RITER[#]",justify='center')
			self.l1.grid(row=23, column=0,columnspan=1)
			self.m16 = tk.OptionMenu(self.t8,self.MenuVar16,'0','1','2','3','4')
			self.m16.config(width=2,bg='lightyellow')
			self.m16.grid(row=23, column=1,columnspan=3)

			self.l1 = tk.Label(self.t8, text="FIX",justify='center')
			self.l1.grid(row=5, column= 6, columnspan=2, sticky='e')
			self.c1 = tk.Checkbutton(self.t8,variable=self.CheckVar36)
			self.c1.grid(row=5, column=8)
	
			self.l1 = tk.Label(self.t8, text="FIX",justify='center')
			self.l1.grid(row=7, column= 6, columnspan=2, sticky='e')
			self.c1 = tk.Checkbutton(self.t8,variable=self.CheckVar37)
			self.c1.grid(row=7, column=8)
				
		self.pres_multi = 1.
		self.l1 = tk.Label(self.t8, text="SUPP-P0[VAL/REF]",justify='center')
		self.l1.grid(row=20, column=0,columnspan=4)
		self.e43 = tk.Entry(self.t8,width=6,justify='center')
		self.e43.insert(10,self.StrVar43.get())
		self.e43.grid(row=20, column=4,columnspan=3)		

		self.l1 = tk.Label(self.t8, text="RBDY [m]",justify='center')
		self.l1.grid(row=21, column=0,columnspan=4)
		self.e51 = tk.Entry(self.t8,width=12,justify='center')
		self.e51.insert(10,self.StrVar51.get())
		self.e51.grid(row=21, column=4,columnspan=6,sticky='w')

		self.l1 = tk.Label(self.t8, text="ZBDY [m]",justify='center')
		self.l1.grid(row=22, column=0,columnspan=4)
		self.e52 = tk.Entry(self.t8,width=12,justify='center')
		self.e52.insert(10,self.StrVar52.get())
		self.e52.grid(row=22, column=4,columnspan=6,sticky='w')

		self.c1 = tk.Checkbutton(self.t8,variable=self.CheckVar6)
		self.c1.grid(row=21, column=8)

		self.l1 = tk.Label(self.t8, text="ITER[#]",justify='center')
		self.l1.grid(row=26, column=0,columnspan=1)
		self.m10 = tk.OptionMenu(self.t8,self.MenuVar10,*self.efit_index,command=lambda value: self.button_func9b(True))
		self.m10.config(width=2)
		self.m10.grid(row=26, column=1,columnspan=3)

		self.l1 = tk.Label(self.t8, text="RUN[#]",justify='center')
		self.l1.grid(row=25, column=0,columnspan=1)
		self.m14 = tk.OptionMenu(self.t8,self.MenuVar14,*self.run_index,command=lambda value: self.button_func9bc())
		self.m14.config(width=2)
		self.m14.grid(row=25, column=1,columnspan=3)		

		self.l1 = tk.Label(self.t8, text="NW[#]",justify='center')
		self.l1.grid(row=24, column=0,columnspan=1)
		self.m13 = tk.OptionMenu(self.t8,self.MenuVar13,'65','129','257')
		self.m13.config(width=2,bg='mistyrose')
		self.m13.grid(row=24, column=1,columnspan=3)				

		b1 = tk.Button(self.t8, text="PLOT RUN", bg = "lightgray",command=lambda: self.button_func9b(),height = 1,width = 15)
		b1.grid(row=24, column=4,columnspan=5)		

		b1 = tk.Button(self.t8, text="RUN EFIT", bg = "mistyrose",command=lambda: self.button_func9c(),height = 1,width = 15)
		b1.grid(row=23, column=4,columnspan=5)

		b1 = tk.Button(self.t8, text="RESET PLOT", bg = "lightgray",command=lambda: delete_efit_runs(self),height = 1,width = 15)
		b1.grid(row=27, column=4,columnspan=5)		

		b1 = tk.Button(self.t8, text="PARAMS CHECK", bg = "lightgray",command=lambda: self.button_func9d(),height = 1,width = 15)
		b1.grid(row=26, column=4,columnspan=5)		

		b1 = tk.Button(self.t8, text="FITTING ERRORS", bg = "lightgray",command=lambda: self.button_func9e(),height = 1,width = 15)
		b1.grid(row=25, column=4,columnspan=5)	

		b1 = tk.Button(self.t8, text="PRESSURE", bg = "lightgray",command=lambda: self.button_func9eb(),height = 1,width = 12)
		b1.grid(row=28, column=0,columnspan=4)		

		b1 = tk.Button(self.t8, text="CURRENT", bg = "lightgray",command=lambda: self.button_func9h(),height = 1,width = 12)
		b1.grid(row=29, column=0,columnspan=4)								

		b1 = tk.Button(self.t8, text="REMOVE FILES", bg = "aliceblue",command=lambda: self.button_func9f(),height = 1,width = 15)
		b1.grid(row=28, column=4,columnspan=5)			

		b1 = tk.Button(self.t8, text="SAVE", bg = "aliceblue",command=lambda: self.button_func9g(),height = 1,width = 15)
		b1.grid(row=29, column=4,columnspan=5)

		if self.ts:
			b1 = tk.Button(self.t8, text="RITER VIEW", bg = "lightyellow",command=lambda: self.button_func9i(),height = 1,width = 12)
			b1.grid(row=27, column=0,columnspan=4)	

		if (float(self.MenuVar1.get()) == 0):
			shot = int(float(self.e1.get()));        time = int(float(self.e2.get()));
			gfile0 = self.currdir +'/%s/EFIT/RESULT/g%06i.%06i_kin_00'%(0,shot,time)
			if not os.path.isfile(gfile0): gfile0 = self.currdir +'/%s/EFIT/RESULT/g%06i.%06i_kin_0'%(0,shot,time)
			eq2 = eqdsk.eqdsk(gfile0)
			eq2.read_eqdsk(gfile0)
			q0 = eq2.q[0]

		self.l1 = tk.Label(self.t8, text="Q0[REF/MOD]",justify='center')
		self.l1.grid(row=19, column=0,columnspan=4)
		self.e46 = tk.Entry(self.t8,width=6,justify='center')
		if (float(self.MenuVar1.get()) > 0): self.e46.insert(10,str(round(self.DoubleVar2.get(),3)))
		else:	self.e46.insert(10,str(round(q0,3)))	
		self.e46.grid(row=19, column=4,columnspan=3)
		self.e46.configure(state='readonly')

		self.e47 = tk.Entry(self.t8,width=6,justify='center')
		self.e47.insert(10,str(round(self.DoubleVar3.get(),3)))
		self.e47.grid(row=19, column=7,columnspan=3)
		self.e47.configure(state='readonly')				


		presf = interp1d(self.psin,self.pres,'cubic')
		ratio = self.pres[0]/presf(0.0375)
		self.e48 = tk.Entry(self.t8,width=6,justify='center')
		self.e48.insert(10,str(round(ratio,2)))
		self.e48.grid(row=20, column=7,columnspan=3)
		self.e48.configure(state='readonly')

		if not self.efit_index[-1] == '--':	self.button_func9bb(int(float(self.MenuVar10.get())))
		delete_efit_runs(self)
		return

	def button_func9b(self,skip=False):

		index = self.MenuVar10.get()
		if index == '--':	return
		index = int(float(index))
		self.button_func9bb(index)
			
		efitdir = self.currdir +'/%s/EFIT/'%(self.MenuVar14.get())	
		save_dir = efitdir + 'RESULT/'
		fit_dir2 = save_dir + 'fitout.dat_%i'%index
		map_dir2 = save_dir + 'map.dat_%i'%index
		mer_dir2 = save_dir + 'mse_er.dat_%i'%index

		if os.path.isfile(mer_dir2):
			f = open(mer_dir2,'r')
			line = f.readline()
			line_n = len(line.split()) - 2;
			print('>>> MSE-Er correction is used with #%i iterations'%line_n)
			print('>>> '+line.split('\n')[0])
			while True:
				line = f.readline().split('\n')[0]
				if not line: break
				print('>>> '+line)
			f.close()

		draw_efit_runs(self,fit_dir2,map_dir2,self.eax1,self.eax2,self.eax3,self.eax4,index,skip)
		self.MenuVar13.set(str(self.efit_num))
		if not self.t9_close:	draw_efit_boundary(self,self.eax5,self.eax6,self.eax7,self.eax8,self.eax9,self.eax16,self.eax17)
		if not self.t10_close:	draw_efit_coils(self,self.eax10,self.eax11,self.eax12,self.eax13)
		if not self.t12_close:	draw_efit_press(self,self.eax14)
		if not self.t13_close:	draw_efit_j(self,self.eax15)
		return

	def button_func9ba(self,index,skip=False):


		filename = self.currdir +'/%s/EFIT/RESULT/save_opt_%i'%(self.MenuVar14.get(),index)
		f = open(filename,'w')

		list1 = [3,5]

		for i in list1:
			f.write('%i\t%s%i\n'%(self.__dict__['CheckVar%d'%i].get(),'!CHECK',int(i)))

		list1 = [17,18,38,39,40,41,42,43]

		i = 16
		f.write('%s\t%s%i\n'%(self.__dict__['StrVar%d'%i].get(),'!STR',int(i)))

		for i in list1:
			if not skip:	self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%i].get())
			f.write('%s\t%s%i\n'%(self.__dict__['StrVar%d'%i].get(),'!STR',int(i)))		

		for i in range(1,250):
			f.write('%i\t%s%i\n'%(self.__dict__['CoilVar%d'%i].get(),'!COIL',int(i)))

		list1 = [51,52]
		for i in list1:
			if not skip:    self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%i].get())
			f.write('%s\t%s%i\n'%(self.__dict__['StrVar%d'%i].get(),'!STR',int(i)))

		list1 = [6]
		for i in list1:
			f.write('%i\t%s%i\n'%(self.__dict__['CheckVar%d'%i].get(),'!CHECK',int(i)))

		for i in range(61,90):
			f.write('%s\t%s%i\n'%(self.__dict__['StrVar%d'%i].get(),'!STR',int(i)))

		for i in range(171,200):
			f.write('%s\t%s%i\n'%(self.__dict__['StrVar%d'%i].get(),'!STR',int(i)))

		list1 = [35]
		for i in list1:
			f.write('%i\t%s%i\n'%(self.__dict__['CheckVar%d'%i].get(),'!CHECK',int(i)))

		list1 = [15,16]
		for i in list1:
			f.write('%s\t%s%i\n'%(self.__dict__['MenuVar%d'%i].get(),'!MENU',int(i)))

		list1 = [36,37]
		for i in list1:
			f.write('%i\t%s%i\n'%(self.__dict__['CheckVar%d'%i].get(),'!CHECK',int(i)))

		f.close()

		return

	def button_func9bb(self,index):

		filename = self.currdir +'/%s/EFIT/RESULT/save_opt_%i'%(self.MenuVar14.get(),index)
		f = open(filename,'r')

		list1 = [3,5]

		for i in list1:
			line = f.readline().split('!')[0]
			self.__dict__['CheckVar%d'%i].set(int(float(line)))

		list1 = [17,18,38,39,40,41,42,43]
		
		try:    line = f.readline().split('\n')[0].split('!')[0].split()[0]
		except: line = ''
		self.StrVar16.set(line)

		for i in list1:
			try:	line = f.readline().split('\n')[0].split('!')[0].split()[0]
			except:	line = ''
			self.__dict__['StrVar%d'%i].set(line)
			self.__dict__['e%d'%i].delete(0,'end')
			self.__dict__['e%d'%i].insert(10,line)

		for i in range(1,250):
			line = f.readline().split('!')[0]
			try:	self.__dict__['CoilVar%d'%i].set(int(float(line)))
			except:	print('>>> '+line)

		list1 = [51,52]
		for i in list1:
			try:    line = f.readline().split('\n')[0].split('!')[0].split()[0]
			except: line = ''
			self.__dict__['StrVar%d'%i].set(line)
			self.__dict__['e%d'%i].delete(0,'end')
			self.__dict__['e%d'%i].insert(10,line)

		list1 = [6]
		for i in list1:
			try:	line = f.readline().split('!')[0]
			except:	line = '0'
			try:	self.__dict__['CheckVar%d'%i].set(int(float(line)))
			except:	self.__dict__['CheckVar%d'%i].set(int(float('0')))

		for i in range(61,90):
			try:    line = f.readline().split('\n')[0].split('!')[0].split()[0]
			except: line = '0.0'
			self.__dict__['StrVar%d'%i].set(line)
	
		for i in range(171,200):
			try:	line = f.readline().split('\n')[0].split('!')[0].split()[0]
			except:	line = self.StrVar16.get()
			self.__dict__['StrVar%d'%i].set(line)

		list1 = [35]
		for i in list1:
			try:	line = f.readline().split('!')[0]
			except: line = '0'
			try:	self.__dict__['CheckVar%d'%i].set(int(float(line)))
			except: self.__dict__['CheckVar%d'%i].set(int(float('0')))
	
		list1 = [15]
		for i in list1:
			try:	line = f.readline().split('\n')[0].split('!')[0].split()[0]
			except: line = '1'
			nn = int(float(line))
			self.__dict__['MenuVar%d'%i].set(str(nn))

		list1 = [16]
		for i in list1:
			try:	line = f.readline().split('\n')[0].split('!')[0].split()[0]
			except: line = '0'
			nn = int(float(line))
			self.__dict__['MenuVar%d'%i].set(str(nn))

		list1 = [36,37]
		for i in list1:
			try:    line = f.readline().split('!')[0]
			except: line = '0'
			try:    self.__dict__['CheckVar%d'%i].set(int(float(line)))
			except: self.__dict__['CheckVar%d'%i].set(int(float('0')))


		f.close()	
		return

	def button_func9bc(self):

		self.efit_index = find_efit_runs(self)
		menu = self.m10["menu"]
		menu.delete(0, "end")
		for string in self.efit_index:
			menu.add_command(label=string, command=lambda value=string: option_menu_update(self.MenuVar10,value,self.button_func9b,True))	
		self.MenuVar10.set(str(string))

		self.button_func9b(True)
		return

	def button_func9c(self):

		riter = int(self.MenuVar16.get())
		if riter > 0:
			self.button_func9c2()
			return

		self.MenuVar14.set(self.MenuVar1.get())
		self.efit_index = find_efit_runs(self)
		
		if self.efit_index[-1] == '--':
			index = 0
		else:
			index = self.efit_index[-1] + 1

		efitdir = self.currdir +'/%s/EFIT/'%(self.MenuVar14.get())
		save_dir = efitdir + 'RESULT/'
		try:	os.mkdir(save_dir)
		except:	pass

		shot = int(float(self.e1.get()));	time = int(float(self.e2.get()));
		kfile_dir = efitdir + 'kfile_run'
		kfile_dir2 = save_dir + 'k%06i.%06i_kin_%i'%(shot,time,index)
		gfile_dir  = efitdir + 'g%06i.%06i'%(shot,time)
		gfile_dir2 = save_dir + 'g%06i.%06i_kin_%i'%(shot,time,index)
		mfile_dir  = efitdir + 'm%06i.%06i'%(shot,time)
		mfile_dir2 = save_dir + 'm%06i.%06i_kin_%i'%(shot,time,index)
		afile_dir  = efitdir + 'a%06i.%06i'%(shot,time)
		afile_dir2 = save_dir + 'a%06i.%06i_kin_%i'%(shot,time,index)
		pfile_dir  = efitdir +'p%06i.%06i'%(shot,time)
		pfile_dir2 = save_dir + 'p%06i.%06i_kin_%i'%(shot,time,index)
		fit_dir  = efitdir  + 'fitout.dat'
		fit_dir2 = save_dir + 'fitout.dat_%i'%index
		w_dir2 = save_dir + 'wmhd.dat_%i'%index
		mse_dir2 = save_dir + 'mse.dat_%i'%index
		pres_dir2 = save_dir + 'pres.dat_%i'%index
		j_dir2 = save_dir + 'jconst.dat_%i'%index
		map_dir2 = save_dir + 'map.dat_%i'%index

		if self.MenuVar13.get() == '257':
			efit_exec = efitdir +'efit257'
			efit_exec2 = efit_dir+'/efit257'
		elif self.MenuVar13.get() == '129':
			efit_exec = efitdir +'efit129'
			efit_exec2 = efit_dir+'/efit129'
		else:
			efit_exec = efitdir +'efit65'
			efit_exec2 = efit_dir+'/efit65'

		if not os.path.isfile(efit_exec):	
			try:	copyfile(efit_exec2,efit_exec)
			except:	
				print('>>> Error in copying efit exe')
				pass
		os.system('chmod 777 %s'%efit_exec)
		if self.CheckVar35.get() == 0:
			run_efit(self)
			make_kfile(self,True)
			write_efit_input(self,kfile_dir)
			os.chdir(efitdir)
			subprocess.run([efit_exec])
			self.ts_run = False
			if not os.path.isfile(gfile_dir):	
				print('>>> RUN FAIL..')
				os.chdir(self.currdir)
				return
			else:
				print('>>> RUN FINISHED..')
	
		elif(self.CheckVar3.get() == 0):
			print('>>> Er correction (by Y.H.LEE) is started')
			if float(self.StrVar151.get()) > 0.:   mseb = 'NB1A'; ind = 152
			elif float(self.StrVar153.get()) > 0.: mseb = 'NB1B'; ind = 154
			else: mseb = 'NB1C'; ind=156
			if (int(float(self.e1.get())) >= 25283):
				if self.aa1gam[0] < 0.83: mseb = 'NB1A'; ind = 152
				else: mseb = 'NB1B'; ind = 154
			print('>>> Beam used in MSE', mseb)
			self.ts_run = False
			run_efit(self)
			a1 = self.StrVar38.get()
			a2 = self.StrVar39.get()
			a3 = self.StrVar40.get()
			self.StrVar38.set('201')
			self.StrVar39.set('0.5')
			self.StrVar40.set('1.e-4')
			make_kfile(self,True)
			write_efit_input(self,kfile_dir)
			os.chdir(efitdir)
			subprocess.run([efit_exec])
			if not os.path.isfile(gfile_dir):
				print('>>> RUN FAIL..')
				return
			self.ts_run = True
			mse_history = dict()
			nn = int(float(self.MenuVar15.get()))
			dtgamma = np.zeros(len(self.dtgamma))
			if self.CheckVar7.get() == 0: dtgamma = np.copy(self.dtgamma)
			mse_history[0] = np.copy(self.tgamma+dtgamma)
			mse_history['q'] = np.zeros(nn+1)
			mse_history['xi']= np.zeros(nn+1)
			read_efit_result(self,'fitout.dat')
			mse_history['q'][0] = self.efit_q[0]
			mse_history['xi'][0]= self.efit_chi
			print('>>> MSE #%i/%i iteration'%(0.,nn))
			for i in range(nn):
				self.tgamma3 = Er_mse_corr(gfile=gfile_dir,pfile='../../'+self.e4.get(),vfile='../../'+self.e45.get(),tgamma=self.tgamma+dtgamma,rgamma=self.rrrgam,ebeam=float(self.__dict__['StrVar%d'%ind].get()))
				if i== nn-1:
					self.StrVar38.set(a1)
					self.StrVar39.set(a2)
					self.StrVar40.set(a3)
				make_kfile(self,True)
				write_efit_input(self,kfile_dir)
				mse_history[i+1] = np.copy(self.tgamma3)
				os.chdir(efitdir)
				subprocess.run([efit_exec])
				if not os.path.isfile(gfile_dir):
					print('>>> RUN FAIL..')
					os.chdir(self.currdir)
					return
				else:
					print('>>> MSE #%i/%i iteration'%(i+1,nn))
	
				read_efit_result(self,'fitout.dat')
				mse_history['q'][i+1] = self.efit_q[0]
				mse_history['xi'][i+1]= self.efit_chi

			f = open('RESULT/mse_er.dat_%i'%index,'w')
			line = '%9s'%'R[m]'
			for i in range(nn+1):line = line +'\t%9s'%('tgam#%1i'%i)
			print('>>> '+line);	f.write(line+'\n');
			for i in range(len(self.rrrgam)):
				line = '%9.6f'%self.rrrgam[i]
				for j in range(nn+1):line = line +'\t%9.6f'%(mse_history[j][i])
				print('>>> '+line);    f.write(line+'\n');

			line = '%9s'%'xisq'
			for i in range(nn+1):line = line +'\t%9.6f'%mse_history['xi'][i]
			print('>>> '+line);    f.write(line+'\n');
			line = '%9s'%'q0'
			for i in range(nn+1):line = line +'\t%9.6f'%mse_history['q'][i]
			print('>>> '+line);    f.write(line+'\n');
			line = '>>>---'
			for i in range(nn+3):   line = line + '---------'
			line = line + '---<<<'
			print('>>> '+line);    f.write(line+'\n');
			f.close()
			print('>>> RUN FINISHED..')	
		else:
			print('>>> Er correction is available with eMSE')
			return

		move(kfile_dir,kfile_dir2)
		move(gfile_dir,gfile_dir2)
		move(mfile_dir,mfile_dir2)
		move(fit_dir,fit_dir2)
		move(afile_dir,afile_dir2)
		copyfile(pfile_dir,pfile_dir2)

		self.button_func9ba(index)

		self.efit_index.append(index)

		menu = self.m10["menu"]
		menu.delete(0, "end")

		efit_post_process(self,index,gfile_dir2,w_dir2,mse_dir2,pres_dir2,j_dir2,map_dir2)
		
		for string in self.efit_index:
			menu.add_command(label=string, command=lambda value=string: option_menu_update(self.MenuVar10,value,self.button_func9b,True))	
					
		draw_efit_runs(self,fit_dir2,map_dir2,self.eax1,self.eax2,self.eax3,self.eax4,index)
		self.efit_list3 = []
		if not self.t9_close:	draw_efit_boundary(self,self.eax5,self.eax6,self.eax7,self.eax8,self.eax9,self.eax16,self.eax17)
		if not self.t10_close:	draw_efit_coils(self,self.eax10,self.eax11,self.eax12,self.eax13)		
		if not self.t12_close:	draw_efit_press(self,self.eax14)
		if not self.t13_close:	draw_efit_j(self,self.eax15)

		os.chdir(self.currdir)

		return

	def button_func9c2(self):

		self.MenuVar14.set(self.MenuVar1.get())
		self.efit_index = find_efit_runs(self)
	
		if self.efit_index[-1] == '--':	 index = 0
		else: index = self.efit_index[-1] + 1

		riter = int(self.MenuVar16.get())

		efitdir  = self.currdir +'/%s/EFIT/'%(self.MenuVar14.get())
		save_dir = efitdir + 'RESULT/'
		temp_dir = save_dir+ '/EFIT_ITER_%i/'%index 
		try:    os.mkdir(save_dir)
		except: pass
		try:	os.mkdir(temp_dir)
		except: pass

		f = open(temp_dir+'/niter','w')
		f.write('%i\n'%(riter))
		f.close()

		shot = int(float(self.e1.get()));       time = int(float(self.e2.get()));
		kfile_dir  = efitdir + 'kfile_run'
		kfile_dir2 = save_dir + 'k%06i.%06i_kin_%i'%(shot,time,index)
		gfile_dir  = efitdir + 'g%06i.%06i'%(shot,time)
		gfile_dir2 = save_dir + 'g%06i.%06i_kin_%i'%(shot,time,index)
		mfile_dir  = efitdir + 'm%06i.%06i'%(shot,time)
		mfile_dir2 = save_dir + 'm%06i.%06i_kin_%i'%(shot,time,index)
		afile_dir  = efitdir + 'a%06i.%06i'%(shot,time)
		afile_dir2 = save_dir + 'a%06i.%06i_kin_%i'%(shot,time,index)
		pfile_dir  = efitdir +'p%06i.%06i'%(shot,time)
		pfile_dir2 = save_dir + 'p%06i.%06i_kin_%i'%(shot,time,index)
		fit_dir  = efitdir  + 'fitout.dat'
		fit_dir2 = save_dir + 'fitout.dat_%i'%index
		w_dir2 = save_dir + 'wmhd.dat_%i'%index
		mse_dir2 = save_dir + 'mse.dat_%i'%index
		pres_dir2 = save_dir + 'pres.dat_%i'%index
		j_dir2 = save_dir + 'jconst.dat_%i'%index
		map_dir2 = save_dir + 'map.dat_%i'%index
		
		if self.MenuVar13.get() == '257':
			efit_exec = efitdir +'efit257'
			efit_exec2 = efit_dir+'/efit257'
		elif self.MenuVar13.get() == '129':
			efit_exec = efitdir +'efit129'
			efit_exec2 = efit_dir+'/efit129'
		else:
			efit_exec = efitdir +'efit65'
			efit_exec2 = efit_dir+'/efit65'
	
		if not os.path.isfile(efit_exec):
			try:    copyfile(efit_exec2,efit_exec)
			except:
				print('>>> Error in copying efit exe')
				pass
		os.system('chmod 777 %s'%efit_exec)

		if (self.CheckVar35.get() > 0 and self.CheckVar3.get() == 0):
			print('>>> Er correction (by Y.H.LEE) is started')

		run_efit(self)
		a1 = self.StrVar38.get()
		a2 = self.StrVar39.get()
		a3 = self.StrVar40.get()
		self.StrVar38.set('201')
		self.StrVar39.set('0.5')
		self.StrVar40.set('1.e-4')

		s1 = self.StrVar17.get()
		s2 = self.StrVar18.get()
		s3 = self.StrVar19.get()

		psin_temp = deepcopy(self.psin)
		pres_temp = deepcopy(self.pres)
		pp_temp   = deepcopy(self.pp)
		ffp_temp  = deepcopy(self.ffp)
		jconst_temp = deepcopy(self.jconst)
		f1 = open(temp_dir+'EFIT_PCONST_0','w')
		f2 = open(temp_dir+'EFIT_JCONST_0','w')
		f1.write('%i\n'%len(psin_temp))
		f2.write('----\n')
		f2.write('%i\n'%len(jconst_temp[:,0]))
		for kk in range(len(psin_temp)): f1.write('%9.6f\t%9.6f\n'%(psin_temp[kk],pres_temp[kk]))
		for kk in range(len(jconst_temp[:,0])): f2.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(jconst_temp[kk,0],1.,jconst_temp[kk,1],1.))
		f1.close(); f2.close();

		for riteri in range(riter+1):
			kfile_dirt = temp_dir + 'k%06i.%06i_kin_%i'%(shot,time,riteri)
			gfile_dirt = temp_dir + 'g%06i.%06i_kin_%i'%(shot,time,riteri)
			mfile_dirt = temp_dir + 'm%06i.%06i_kin_%i'%(shot,time,riteri)
			afile_dirt = temp_dir + 'a%06i.%06i_kin_%i'%(shot,time,riteri)
			pfile_dirt = temp_dir + 'p%06i.%06i_kin_%i'%(shot,time,riteri)
			fit_dirt   = temp_dir + 'fitout.dat_%i'%riteri
			w_dirt     = temp_dir + 'wmhd.dat_%i'%riteri
			mse_dirt   = temp_dir + 'mse.dat_%i'%riteri
			pres_dirt  = temp_dir + 'pres.dat_%i'%riteri
			j_dirt     = temp_dir + 'jconst.dat_%i'%riteri
			map_dirt   = temp_dir + 'map.dat_%i'%riteri

			if riteri > 0:
				os.chdir(temp_dir)
				if not os.path.isfile(efitdir+'/chease_opt'): make_chease_opt(self,'chease_opt')
				else: copyfile(efitdir+'/chease_opt',temp_dir+'/chease_opt')
				gfile0 = self.currdir + '/' + self.e3.get()
				gfile1  = temp_dir + 'g%06i.%06i_kin_%i'%(shot,time,riteri-1)
				renew_efit_constraint_kin(self,gfile0,gfile1)
				pfile_dir = pfile_dirt
				beam_pf   = efitdir + 'chease_pres'
				beam_jf   = efitdir + 'chease_curr'
				copyfile('chease_kinprof_new',pfile_dir)	
				renew_efit_constraint_pj(gfile1,pfile_dir,beam_pf,beam_jf)
				copyfile('EFIT_PCONST','EFIT_PCONST_%i'%riteri)
				copyfile('OUTPUT/EFIT_JCONST','EFIT_JCONST_%i'%riteri)
				get_efit_constraint_cu(self,'OUTPUT/EFIT_JCONST',True)
				f = open('EFIT_PCONST','r')
				datn = int(f.readline())
				self.psin = np.zeros(datn)
				self.pres = np.zeros(datn)
				self.pp   = np.zeros(datn)
				self.ffp  = np.zeros(datn)
				ffpf = interp1d(self.jconst[:,0],self.jconst[:,2])
				for k in range(datn):
					line = f.readline().split()
					self.psin[k] = float(line[0])
					self.pres[k] = float(line[1])
					self.ffp[k]  = ffpf(self.psin[k])
				presf = interp1d(self.psin,self.pres)
				for k in range(1,datn-1):
					ppp = self.psin[k]
					self.pp[k] = (presf(ppp+1.e-5)-presf(ppp-1.e-5))/2.e-5
				self.pp[0] = 2.0 * self.pp[1] - self.pp[2]
				self.pp[-1] = 2.0 * self.pp[-2] - self.pp[-3]

				f.close()
				os.chdir(self.currdir)
				sknot = make_current_knots(self,int(float(self.StrVar34.get())),int(float(self.StrVar35.get())), float(self.StrVar36.get()),float(self.StrVar37.get()))
				self.StrVar19.set(sknot)

				if self.CheckVar36.get() == 0:
					sknot = make_pf_knots(self,1)
					self.StrVar17.set(sknot)
				if self.CheckVar37.get() == 0:
					sknot = make_pf_knots(self,2)
					self.StrVar18.set(sknot)

			else:
				copyfile(pfile_dir,pfile_dirt)
				

			if self.CheckVar35.get() == 0:
				if riteri == (riter):
					self.StrVar38.set(a1)
					self.StrVar39.set(a2)
					self.StrVar40.set(a3)

#				run_efit(self)
				make_kfile(self,True)
				write_efit_input(self,kfile_dir)
				os.chdir(efitdir)
				subprocess.run([efit_exec])
				self.ts_run = False
				if not os.path.isfile(gfile_dir):
					print('>>> RUN FAIL..')
					os.chdir(self.currdir)
					return
				else:
					print('>>> ITER #%i/%i RUN FINISHED..'%(riteri+1,riter+1))
	
			elif(self.CheckVar3.get() == 0):
				if float(self.StrVar151.get()) > 0.:   mseb = 'NB1A'; ind = 152
				elif float(self.StrVar153.get()) > 0.: mseb = 'NB1B'; ind = 154
				else: mseb = 'NB1C'; ind=156
				if (int(float(self.e1.get())) >= 25283):
					if self.aa1gam[0] < 0.83: mseb = 'NB1A'; ind = 152
					else: mseb = 'NB1B'; ind = 154
				print('>>> Beam used in MSE', mseb)
				self.ts_run = False
#				run_efit(self)
				make_kfile(self,True)
				write_efit_input(self,kfile_dir)
				os.chdir(efitdir)
				subprocess.run([efit_exec])
				if not os.path.isfile(gfile_dir):
					print('>>> RUN FAIL..')
					return
				self.ts_run = True
				mse_history = dict()
				nn = int(float(self.MenuVar15.get()))
				dtgamma = np.zeros(len(self.dtgamma))
				if self.CheckVar7.get() == 0: dtgamma = np.copy(self.dtgamma)
				mse_history[0] = np.copy(self.tgamma+dtgamma)
				mse_history['q'] = np.zeros(nn+1)
				mse_history['xi']= np.zeros(nn+1)
				read_efit_result(self,'fitout.dat')
				mse_history['q'][0] = self.efit_q[0]
				mse_history['xi'][0]= self.efit_chi
				print('>>> ITER #%i/%i > MSE #%i/%i iteration'%(riteri+1,riter+1,0.,nn))
				for i in range(nn):
					self.tgamma3 = Er_mse_corr(gfile=gfile_dir,pfile='../../'+self.e4.get(),vfile='../../'+self.e45.get(),tgamma=self.tgamma+dtgamma,rgamma=self.rrrgam,ebeam=float(self.__dict__['StrVar%d'%ind].get()))
					if ((i== nn-1) and (riteri == (riter))):
						self.StrVar38.set(a1)
						self.StrVar39.set(a2)
						self.StrVar40.set(a3)
					make_kfile(self,True)
					write_efit_input(self,kfile_dir)
					mse_history[i+1] = np.copy(self.tgamma3)
					os.chdir(efitdir)
					subprocess.run([efit_exec])
					if not os.path.isfile(gfile_dir):
						print('>>> RUN FAIL..')
						os.chdir(self.currdir)
						return
					else:
						print('>>> ITER #%i/%i > MSE #%i/%i iteration'%(riteri+1,riter+1,i+1,nn))
	
					read_efit_result(self,'fitout.dat')
					mse_history['q'][i+1] = self.efit_q[0]
					mse_history['xi'][i+1]= self.efit_chi
				
				f = open(temp_dir+'/mse_er.dat_%i'%riteri,'w')
				line = '%9s'%'R[m]'
				for i in range(nn+1):line = line +'\t%9s'%('tgam#%1i'%i)
				print('>>> '+line);    f.write(line+'\n');
				for i in range(len(self.rrrgam)):
					line = '%9.6f'%self.rrrgam[i]
					for j in range(nn+1):line = line +'\t%9.6f'%(mse_history[j][i])
					print('>>> '+line);    f.write(line+'\n');
			
				line = '%9s'%'xisq'
				for i in range(nn+1):line = line +'\t%9.6f'%mse_history['xi'][i]
				print('>>> '+line);    f.write(line+'\n');
				line = '%9s'%'q0'
				for i in range(nn+1):line = line +'\t%9.6f'%mse_history['q'][i]
				print('>>> '+line);    f.write(line+'\n');
				line = '>>>---'
				for i in range(nn+3):   line = line + '---------'
				line = line + '---<<<'
				print(line);    f.write(line+'\n');
				f.close()
				print('>>> ITER #%i/%i RUN FINISHED..'%(riteri+1,riter+1))
			else:
				print('>>> Er correction is available with eMSE')
				return

			copyfile(kfile_dir,kfile_dirt)
			copyfile(gfile_dir,gfile_dirt)
			copyfile(mfile_dir,mfile_dirt)
			copyfile(fit_dir,fit_dirt)
			copyfile(afile_dir,afile_dirt)
			print('>>> POST PROCESS for ITER #%i'%(riteri+1))
			eq = eqdsk.eqdsk(gfile_dirt)
			eq.read_eqdsk_file()
			eq.make_grid()
			eq.get_flux_contour()
			eq.make_rho_R_psin()
			eq.surface_zjz()

			f1 = open(temp_dir+'prho_%i.dat'%riteri,'w')
			f2 = open(temp_dir+'qpres_%i.dat'%riteri,'w')
			f3 = open(temp_dir+'javg_%i.dat'%riteri,'w')
			f1.write('%i\n'%len(eq.prhoR[:,0]))
			f2.write('%i\n'%len(eq.psin))
			f3.write('%i\n'%len(eq.psinh))			

			for kk in range(len(eq.psin)): f2.write('%9.6f\t%9.6f\t%9.6f\n'%(eq.psin[kk],eq.pres[kk],eq.q[kk]))
			for kk in range(len(eq.prhoR[:,0])): f1.write('%9.6f\t%9.6f\t%9.6f\n'%(eq.prhoR[kk,0],eq.prhoR[kk,1],eq.prhoR[kk,2]))
			for kk in range(len(eq.psinh)): f3.write('%9.6f\t%9.6f\n'%(eq.psinh[kk],eq.zjzh[kk]/1.e6))
			f1.close()
			f2.close()
			f3.close()
			efit_post_process(self,index,gfile_dirt,w_dirt,mse_dirt,pres_dirt,j_dirt,map_dirt)

			self.psin = deepcopy(psin_temp)
			self.pres = deepcopy(pres_temp)
			self.jconst = deepcopy(jconst_temp)
			self.pp   = deepcopy(pp_temp)
			self.ffp  = deepcopy(ffp_temp)

		move(kfile_dir,kfile_dir2)
		move(gfile_dir,gfile_dir2)
		move(mfile_dir,mfile_dir2)
		move(fit_dir,fit_dir2)
		move(afile_dir,afile_dir2)
		copyfile(pfile_dir,pfile_dir2)
		if (self.CheckVar35.get() > 0 and self.CheckVar3.get() == 0):
			copyfile(temp_dir+'/mse_er.dat_%i'%riteri,'RESULT/mse_er.dat_%i'%index)
		
		self.button_func9ba(index)
		
		self.efit_index.append(index)
		
		menu = self.m10["menu"]
		menu.delete(0, "end")

		copyfile(w_dirt,w_dir2); copyfile(mse_dirt,mse_dir2); copyfile(pres_dirt,pres_dir2); copyfile(j_dirt,j_dir2); copyfile(map_dirt,map_dir2)
#		efit_post_process(self,index,gfile_dir2,w_dir2,mse_dir2,pres_dir2,j_dir2,map_dir2)

		for string in self.efit_index:
			menu.add_command(label=string, command=lambda value=string: option_menu_update(self.MenuVar10,value,self.button_func9b,True))
	
		draw_efit_runs(self,fit_dir2,map_dir2,self.eax1,self.eax2,self.eax3,self.eax4,index)
		self.efit_list3 = []
		if not self.t9_close:   draw_efit_boundary(self,self.eax5,self.eax6,self.eax7,self.eax8,self.eax9,self.eax16,self.eax17)
		if not self.t10_close:  draw_efit_coils(self,self.eax10,self.eax11,self.eax12,self.eax13)
		if not self.t12_close:  draw_efit_press(self,self.eax14)
		if not self.t13_close:  draw_efit_j(self,self.eax15)

		os.chdir(self.currdir)
		return

	def button_func9d(self):

		if not self.t9_close:
			print('>>> Boundary window is already opened...')
			return					
		self.window_bnd = True	
		get_efit_constraint_eq(self,True)

		
		self.efit_index = find_efit_runs(self)
		index = self.MenuVar10.get()
		if index == '--':	return		

		self.t9 = tk.Toplevel(self.t8)
		self.t9.wm_title("EFIT Params")
		self.t9_close = False		
		self.t9.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(9))
		
		self.t9.resizable(0,0)
		self.fig9 = plt.figure(figsize=(17,7))	

		self.eax5 = plt.subplot2grid((2,4),(0,0),rowspan=2,fig=self.fig9)
		self.eax6 = plt.subplot2grid((2,4),(0,1),fig=self.fig9)
		self.eax7 = plt.subplot2grid((2,4),(0,2),fig=self.fig9)
		self.eax8 = plt.subplot2grid((2,4),(1,1),fig=self.fig9)
		self.eax9 = plt.subplot2grid((2,4),(1,2),fig=self.fig9)
		self.eax16 = plt.subplot2grid((2,4),(0,3),fig=self.fig9)
		self.eax17 = plt.subplot2grid((2,4),(1,3),fig=self.fig9)
		
		self.canvas9 = FigureCanvasTkAgg(self.fig9,master=self.t9)
		self.plot_widget9 = self.canvas9.get_tk_widget()
		self.plot_widget9.grid(rowspan=30,row=2,column=10,columnspan=40)

		toolbar_frame = tk.Frame(self.t9)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar9 = NavigationToolbar2Tk(self.canvas9,toolbar_frame)

		draw_efit_boundary(self,self.eax5,self.eax6,self.eax7,self.eax8,self.eax9,self.eax16,self.eax17)

		return

	def button_func9e(self):

		if not self.t10_close:
			print('>>> Fitting error window is already opened...')
			return					
		self.window_coil = True	
		get_efit_constraint_eq(self,True)
		self.efit_index = find_efit_runs(self)
		index = self.MenuVar10.get()
		if index == '--':	return		

		self.t10 = tk.Toplevel(self.t8)
		self.t10.wm_title("Fitting errors")
		self.t10_close = False		
		self.t10.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(10))		

		self.fig10, ([self.eax10, self.eax11], [self.eax12, self.eax13]) = plt.subplots(2,2,figsize=(13,7))		

		self.canvas10 = FigureCanvasTkAgg(self.fig10,master=self.t10)
		self.plot_widget10 = self.canvas10.get_tk_widget()
		self.plot_widget10.grid(rowspan=30,row=2,column=10,columnspan=40)

		toolbar_frame = tk.Frame(self.t10)
		toolbar_frame.grid(column=10,row=0)
		
		self.toolbar10 = NavigationToolbar2Tk(self.canvas10,toolbar_frame)
		self.fig10.canvas.draw_idle()

		draw_efit_coils(self,self.eax10,self.eax11,self.eax12,self.eax13)

		return

	def button_func9eb(self):

		if not self.t12_close:
			print('>>> Pressure window is already opened...')
			return					
		self.window_pres = True	
		get_efit_constraint_eq(self,True)
		self.efit_index = find_efit_runs(self)
		index = self.MenuVar10.get()
		if index == '--':	return		

		self.t12 = tk.Toplevel(self.t8)
		self.t12.wm_title("PRESSURE profile")
		self.t12_close = False		
		self.t12.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(12))		

#		self.t12.resizable(0,0)
		self.fig12, self.eax14 = plt.subplots(1,1,figsize=(5,5))		

		self.canvas12 = FigureCanvasTkAgg(self.fig12,master=self.t12)
		self.plot_widget12 = self.canvas12.get_tk_widget()
		self.plot_widget12.grid(row=1,column=1)

		toolbar_frame = tk.Frame(self.t12)
		toolbar_frame.grid(column=1,row=0)
		
		self.toolbar12 = NavigationToolbar2Tk(self.canvas12,toolbar_frame)
		self.fig12.canvas.draw_idle()

		draw_efit_press(self,self.eax14)

		return		

	def button_func9f(self):

		self.t11 = tk.Toplevel(self.root)
		self.t11.wm_title('DELETE FILE')
		self.t11.resizable(0,0)
		
		self.l2 = tk.Label(self.t11, text = ' Do you want to remove previous EFIT runs ?',justify='center')
		self.l2.grid(row=0,column=0,columnspan=4,pady=15)

		self.b1 = tk.Button(self.t11, text="YES", bg = "lightgray",command=lambda: self.button_func9fa(True),height = 1,width = 5)
		self.b1.grid(row=1, column=1, columnspan=1)	
		self.b1 = tk.Button(self.t11, text="NO", bg = "lightgray",command=lambda: self.button_func9fa(False),height = 1,width = 5)
		self.b1.grid(row=1, column=2, columnspan=1)	

		self.t11.mainloop()

		return

	def button_func9fa(self,input):
		if input:
			efitdir = self.currdir +'/%s/EFIT/RESULT'%(self.MenuVar14.get())
			rmtree(efitdir)
			os.mkdir(efitdir)	

			menu = self.m10["menu"]
			menu.delete(0,"end")
			for string in ['--']:
				menu.add_command(label=string, command=lambda value=string: option_menu_update(self.MenuVar10,value,self.button_func9b,True))
			self.MenuVar10.set('--')

		self.t11.destroy()
		self.t11.quit()
		return

	def button_func9g(self):

		index = self.MenuVar10.get()
		if index == '--':	return
		index = int(float(index))

		efitdir = self.currdir +'/%s/EFIT/'%(self.MenuVar14.get())
		save_dir = efitdir + 'RESULT/'
		shot = int(float(self.e1.get()));	time = int(float(self.e2.get()));
		kfile_dir2 = save_dir + 'k%06i.%06i_kin_%i'%(shot,time,index)
		gfile_dir2 = save_dir + 'g%06i.%06i_kin_%i'%(shot,time,index)
		afile_dir2 = save_dir + 'a%06i.%06i_kin_%i'%(shot,time,index)
		pfile_dir2 = save_dir + 'p%06i.%06i_kin_%i'%(shot,time,index)

		try: bdiff = round(float(self.MenuVar6.get()),2)
		except: print('>>> Select the Bdiff first!'); return
		if self.MenuVar7.get().lower() == 'smse':
			nfile_dir2 = self.currdir + '/%s/CSOLVE/Bdiff_%03i/PROFILES/neo_coefs_%s'%(self.MenuVar1.get(),round(bdiff*100.,0),self.MenuVar9.get())
			bfile_dir2 = self.currdir + '/%s/CSOLVE/Bdiff_%03i/PROFILES/nubeam_out1d_%s'%(self.MenuVar1.get(),round(bdiff*100.,0),self.MenuVar9.get())
			Bfile_dir2 = self.currdir + '/%s/CSOLVE/Bdiff_%03i/PROFILES/nubeam_out0d_%s'%(self.MenuVar1.get(),round(bdiff*100.,0),self.MenuVar9.get())
		else:
			nfile_dir2 = self.currdir + '/%s/NUBEAM/Bdiff_%03i/neo_coefs'%(self.MenuVar1.get(),round(bdiff*100.,0))
			bfile_dir2 = self.currdir + '/%s/NUBEAM/Bdiff_%03i/nubeam_out1d'%(self.MenuVar1.get(),round(bdiff*100.,0))
			Bfile_dir2 = self.currdir + '/%s/NUBEAM/Bdiff_%03i/nubeam_out0d'%(self.MenuVar1.get(),round(bdiff*100.,0))

		outdir = self.currdir+'/OUTPUT'

		try:	os.mkdir(outdir)
		except:	pass

		pfile = pfile_convert.GENpFile()
		pfile.filename = outdir+'/p%06i.%06i_kin_%s'%(shot,time,self.MenuVar1.get())
		pfile._make_pfile(gfile_dir2,pfile_dir2,bfile_dir2,nfile_dir2)

		try:
			copyfile(gfile_dir2,outdir+'/g%06i.%06i_kin_%s'%(shot,time,self.MenuVar1.get()))
			copyfile(kfile_dir2,outdir+'/k%06i.%06i_kin_%s'%(shot,time,self.MenuVar1.get()))
			copyfile(afile_dir2,outdir+'/a%06i.%06i_kin_%s'%(shot,time,self.MenuVar1.get()))
			copyfile(pfile_dir2,outdir+'/c%06i.%06i_kin_%s'%(shot,time,self.MenuVar1.get()))
		except:
			pass
		try: 
			copyfile(bfile_dir2,outdir+'/f%06i.%06i_kin_%s'%(shot,time,self.MenuVar1.get()))
			copyfile(Bfile_dir2,outdir+'/F%06i.%06i_kin_%s'%(shot,time,self.MenuVar1.get()))
			copyfile(nfile_dir2,outdir+'/n%06i.%06i_kin_%s'%(shot,time,self.MenuVar1.get()))
		except: print('>>> Beam file is not ready, re-run the iteration to get full output!')

		try: pscale = float(self.StrVar43.get().split(',')[1])
		except: pscale = 1.
		f = open(outdir+'/efit_opt_%s'%(self.MenuVar1.get()),'w')
		f.write('pscale  = %s\n'%pscale)
		f.write('q0fix   = %s\n'%self.CheckVar5.get())
		f.write('q0val   = %s\n'%self.StrVar42.get())
		f.write('Zeff    = %s\n'%self.StrVar7.get())
		f.write('bsmodel = %s\n'%self.MenuVar3.get())
		f.write('bs_mult = %s\n'%self.StrVar13.get())
		f.write('Beam_dif= %s\n'%self.MenuVar6.get())
		f.write('syn_MSE = %s\n'%self.CheckVar3.get())
		f.close()

		print('>>> Result is saved to OUTPUT dirs with g-file, k-file and p-file (profile)')

		self.t8.destroy()
		self.t8_close = True
		self.button_func10a()
		return

	def button_func9h(self):

		if not self.t13_close:
			print('>>> Current window is already opened...')
			return					
		self.window_j = True	
		get_efit_constraint_eq(self,True)
		self.efit_index = find_efit_runs(self)
		index = self.MenuVar10.get()
		if index == '--':	return		

		self.t13 = tk.Toplevel(self.t8)
		self.t13.wm_title("Current profile")
		self.t13_close = False		
		self.t13.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(13))		

#		self.t13.resizable(0,0)
		self.fig13, self.eax15 = plt.subplots(1,1,figsize=(5.,5.))		

		self.canvas13 = FigureCanvasTkAgg(self.fig13,master=self.t13)
		self.plot_widget13 = self.canvas13.get_tk_widget()
		self.plot_widget13.grid(row=1,column=1)

		toolbar_frame = tk.Frame(self.t13)
		toolbar_frame.grid(column=1,row=0)
		
		self.toolbar13 = NavigationToolbar2Tk(self.canvas13,toolbar_frame)
		self.fig13.canvas.draw_idle()

		draw_efit_j(self,self.eax15)
		return

	def button_func9i(self):

		if not self.t18_close:
			print('>>> RITER window is already opened...')

		self.window_riter = True

		index = self.MenuVar10.get()
		if index == '--':       return

		efitdir  = self.currdir +'/%s/EFIT/'%(self.MenuVar14.get())
		temp_dir = efitdir+ '/RESULT/EFIT_ITER_%i/'%int(self.MenuVar10.get())

		if not os.path.isdir(temp_dir):
			print('>>> No RIteration history...')
			return

		self.t18 = tk.Toplevel(self.t8)

		self.t18.wm_title("RITER overview")
		self.t18_close = False
		self.t18.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(18))

		self.t18.resizable(0,0)

		self.fig18, ([self.eax19, self.eax20], [self.eax21, self.eax22]) = plt.subplots(2,2,figsize=(11.,7.))
		self.canvas18 = FigureCanvasTkAgg(self.fig18,master=self.t18)
		self.plot_widget18 = self.canvas18.get_tk_widget()
		self.plot_widget18.grid(row=1,column=1)

		toolbar_frame = tk.Frame(self.t18)
		toolbar_frame.grid(column=1,row=0)

		self.toolbar18 = NavigationToolbar2Tk(self.canvas18,toolbar_frame)
		self.fig18.canvas.draw_idle()

		draw_efit_riter(self,self.eax19,self.eax20,self.eax21,self.eax22)
	
		return

	def button_func10a(self):

		self.ind = int(float(self.MenuVar1.get()))

		save_dir = self.currdir+'/%i/save.dat'%self.ind
		self.save_opt(save_dir)
		save_dir2 = self.currdir+'/0/comment.dat'
		save_comment(self,save_dir2)

		print('>>> Parameters are saved in %s'%save_dir)
		self.update_status()
		return		

	def button_func10b(self):

		self.root.destroy()
		return

	def button_func11(self,rtype=2):

		if (self.ismse and float(self.StrVar16.get()) == 0.):	self.StrVar16.set('0.03')

		if not self.t16_close:
			print('>>> Current SGAMMA is already opened...')
			return

		self.t16 = tk.Toplevel(self.root)
		self.t16.wm_title("SGAMMA")
		self.t16_close = False
		self.t16.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(16))

		self.t16.resizable(0,0)

		sgam_sum = 0.

		for i in range(171,171+len(self.sgamma)):
			if self.__dict__['StrVar%d'%i].get() == '':     self.__dict__['StrVar%d'%i].set('0')
			sgam_sum = sgam_sum + float(self.__dict__['StrVar%d'%i].get())
		if sgam_sum == 0.:
			for i in range(171,171+len(self.sgamma)):
				if (float(self.StrVar16.get()) == 0. and self.ismse):	self.__dict__['StrVar%d'%i].set('%5.4e'%self.sgamma[i-171])
				else:	self.__dict__['StrVar%d'%i].set(self.StrVar16.get())

		self.l2 = tk.Label(self.t16, text = '============  MSE-SGAMMA ============',justify='center')
		self.l2.grid(row=0,column=0,columnspan=5)

		i = 1;	j=0; k = 0;

		for k in range(len(self.sgamma)):

			self.__dict__['e%d'%(171+k)] = tk.Entry(self.t16,width=10,justify='center')
			self.__dict__['e%d'%(171+k)].insert(10,self.__dict__['StrVar%d'%(171+k)].get())
			self.__dict__['e%d'%(171+k)].grid(row=i,column=j)

			k = k + 1
			j = j + 1
			if j == 5:
				j = 0;
				i = i + 1
			
		self.l1 = tk.Label(self.t16, text="UNIFORM MSE-SGAM ",justify='center')
		self.l1.grid(row=i+1, column=0,columnspan=2,sticky='e')
		self.e16 = tk.Entry(self.t16,width=8,justify='center')
		self.e16.insert(10,self.StrVar16.get())
		self.e16.grid(row=i+1, column=2,columnspan=1)	
	
		b1 = tk.Button(self.t16, text="SET", bg = "lightgray",command=lambda: self.button_func11a(),height = 1,width = 5)
		b1.grid(row=i+1,column=3)
		b1 = tk.Button(self.t16, text="RESET", bg = "lightgray",command=lambda: self.button_func11b(),height = 1,width = 5)
		b1.grid(row=i+1,column=4)
		b1 = tk.Button(self.t16, text="SAVE", bg = "lightgray",command=lambda: self.button_func11c(rtype),height = 1,width = 5)		
		b1.grid(row=i+2,column=0,columnspan=5)

		return
	
	def button_func11a(self):
		self.StrVar16.set(self.e16.get())
		if float(self.StrVar16.get()) <= 0.:	return
		for i in range(171,171+len(self.sgamma)):
			self.__dict__['StrVar%d'%i].set(self.StrVar16.get())
			self.__dict__['e%d'%i].delete(0,'end')	
			self.__dict__['e%d'%i].insert(10,self.__dict__['StrVar%d'%i].get())
		return

	def button_func11b(self):

		for i in range(171,171+len(self.sgamma)):
			self.__dict__['StrVar%d'%i].set('%5.4e'%self.sgamma[i-171])
			self.__dict__['e%d'%i].delete(0,'end')
			self.__dict__['e%d'%i].insert(10,self.__dict__['StrVar%d'%i].get())
		return

	def button_func11c(self,rtype=1):

		for i in range(171,171+len(self.sgamma)):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%i].get())
			self.sgamma2[i-171] = float(self.__dict__['e%d'%i].get())

		self.StrVar16.set(self.e16.get())
		if rtype==1:	draw_mse_constraint(self)
		self.t16_close = True
		self.is_sgam = True
		self.t16.destroy()
		return

	def button_func12(self,rtype=1):

		if not self.t19_close:
			print('>>> Current DTGAMMA is already opened...')
			return

		self.t19 = tk.Toplevel(self.root)
		self.t19.wm_title("DTGAMMA")
		self.t19_close = False
		self.t19.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(19))

		self.t19.resizable(0,0)

		for i in range(61,61+len(self.sgamma)):
			if self.__dict__['StrVar%d'%i].get() == '':     self.__dict__['StrVar%d'%i].set('0.0')

		self.l2 = tk.Label(self.t19, text = '============  MSE-DTGAMMA ============',justify='center')
		self.l2.grid(row=0,column=0,columnspan=5)

		i = 1;  j=0; k = 0;

		for k in range(len(self.sgamma)):
			self.__dict__['e%d'%(61+k)] = tk.Entry(self.t19,width=10,justify='center')
			self.__dict__['e%d'%(61+k)].insert(10,self.__dict__['StrVar%d'%(61+k)].get())
			self.__dict__['e%d'%(61+k)].grid(row=i,column=j)

			k = k + 1
			j = j + 1
			if j == 5: j = 0; i = i + 1

		b1 = tk.Button(self.t19, text="SAVE", bg = "lightgray",command=lambda: self.button_func12a(rtype),height = 1,width = 5)
		b1.grid(row=i+2,column=0,columnspan=5)

		self.l2 = tk.Label(self.t19, text = 'sMSE',justify='center')
		self.l2.grid(row=i+2,column=0,columnspan=1)

		self.c1 = tk.Checkbutton(self.t19,variable=self.CheckVar7)
		self.c1.grid(row=i+2, column=1,sticky='w')


		return

	def button_func12a(self,rtype=1):

		for i in range(61,61+len(self.sgamma)):
			self.__dict__['StrVar%d'%i].set(self.__dict__['e%d'%i].get())
			self.dtgamma[i-61] = float(self.__dict__['e%d'%i].get())

		if rtype==1:    draw_mse_constraint(self)
		self.t19_close = True
		self.t19.destroy()
		return

		
	def gui_efit(self):

		self.fig, ax1 = plt.subplots(1,1,figsize=(5,9))
		ax1.set_title('Equilibrium')
		ax1.set_xlabel('R [m]')
		ax1.set_ylabel('Z [m]')

		self.fig.tight_layout()
		self.root.resizable(0,0)
		self.canvas = FigureCanvasTkAgg(self.fig,master=self.root)
		self.plot_widget = self.canvas.get_tk_widget()
		self.plot_widget.grid(rowspan=50,row=1,column=12,columnspan=10)

		toolbar_frame = tk.Frame(self.root)
		toolbar_frame.grid(column=12,row=0)
		
		self.toolbar = NavigationToolbar2Tk(self.canvas,toolbar_frame)

		self.initialise_vars(True)	

		self.l1 = tk.Label(self.root, text="==================== G-EFIT toolkit  ====================",justify='center')
		self.l1.grid(row=0, column=0,columnspan=12)
		self.l1 = tk.Label(self.root, text="RUN [#]",justify='center')
		self.l1.grid(row=1, column=0,columnspan=2,sticky='e')
		self.m1 = tk.OptionMenu(self.root,self.MenuVar1,*self.run_index, command=lambda value: self.button_func1b(3))
		self.m1.config(width=2)
		self.m1.grid(row=1, column=2,columnspan=2,sticky='e')
		b1 = tk.Button(self.root, text="GEN", bg = "lightgray",command=lambda: self.button_func1a(),height = 1,width = 4)
		b1.grid(row=1, column=4,columnspan=2,sticky='w')		

		self.l1 = tk.Label(self.root, text="LOAD [#]",justify='center')
		self.l1.grid(row=1, column=6,columnspan=2)
		self.m2 = tk.OptionMenu(self.root,self.MenuVar2,*self.run_index2)
		self.m2.config(width=2)
		self.m2.grid(row=1, column=8,columnspan=2)
		b1 = tk.Button(self.root, text="LOAD", bg = "lightgray",command=lambda: self.button_func1b(),height = 1,width = 4)
		b1.grid(row=1, column=10,columnspan=2,sticky='w')


		self.l1 = tk.Label(self.root, text="=====================   INPUTS   ======================",justify='center')
		self.l1.grid(row=2, column=0,columnspan=12)


		self.l1 = tk.Label(self.root, text="Shot [#]",justify='center')
		self.l1.grid(row=3, column=0,columnspan=2)
		self.e1 = tk.Entry(self.root,width=6,justify='center')
		self.e1.insert(10,self.StrVar1.get())
		self.e1.grid(row=3, column=2,columnspan=2)

		self.l1 = tk.Label(self.root, text="Time [ms]",justify='center')
		self.l1.grid(row=3, column=4,columnspan=2)
		self.e2 = tk.Entry(self.root,width=7,justify='center')
		self.e2.insert(10,self.StrVar2.get())
		self.e2.grid(row=3, column=6,columnspan=2)		

		b2 = tk.Button(self.root, text="Open MDS", bg = "lightgray",command=lambda: self.button_func2(),height = 1,width = 8)
		b2.grid(row=3, column=8,columnspan=4)

		self.l1 = tk.Label(self.root, text="EQU",justify='center')
		self.l1.grid(row=4, column=0,columnspan=2)
		self.e3 = tk.Entry(self.root,width=30,justify='center')
		self.e3.insert(10,self.StrVar3.get())
		self.e3.grid(row=4, column=1,columnspan=7,sticky='e')
		b3 = tk.Button(self.root, text="LOAD", bg = "lightgray",command=lambda: self.button_func3a(),height = 1,width = 4)
		b3.grid(row=4, column=8,columnspan=2,sticky='w')
		b3 = tk.Button(self.root, text="CHECK", bg = "lightgray",command=lambda: self.button_func3b(),height = 1,width = 4)
		b3.grid(row=4, column=10,columnspan=2,sticky='w')		

		self.l1 = tk.Label(self.root, text="KIN",justify='center')
		self.l1.grid(row=5, column=0,columnspan=2)
		self.e4 = tk.Entry(self.root,width=30,justify='center')
		self.e4.insert(10,self.StrVar4.get())
		self.e4.grid(row=5, column=1,columnspan=7,sticky='e')
		b4 = tk.Button(self.root, text="LOAD", bg = "lightgray",command=lambda: self.button_func4a(),height = 1,width = 4)
		b4.grid(row=5, column=8,columnspan=2,sticky='w')
		b4 = tk.Button(self.root, text="CHECK", bg = "lightgray",command=lambda: self.button_func4b(),height = 1,width = 4)
		b4.grid(row=5, column=10,columnspan=2,sticky='w')

		self.l1 = tk.Label(self.root, text="RFILE",justify='center')
		self.l1.grid(row=6, column=0,columnspan=2)
		self.e45 = tk.Entry(self.root,width=30,justify='center')
		self.e45.insert(10,self.StrVar45.get())
		self.e45.grid(row=6, column=1,columnspan=7,sticky='e')
		b4 = tk.Button(self.root, text="LOAD", bg = "lightgray", command=lambda: self.button_func4c(),height = 1,width = 4)
		b4.grid(row=6, column=8,columnspan=2,sticky='w')	
		b4 = tk.Button(self.root, text="CHECK", bg = "lightgray",command=lambda: self.button_func4d(),height = 1,width = 4)
		b4.grid(row=6, column=10,columnspan=2,sticky='w')			

		self.l1 = tk.Label(self.root, text="KFILE",justify='center')
		self.l1.grid(row=7, column=0,columnspan=2)
		self.e5 = tk.Entry(self.root,width=30,justify='center')
		self.e5.insert(10,self.StrVar5.get())
		self.e5.grid(row=7, column=1,columnspan=7,sticky='e')
		b4 = tk.Button(self.root, text="LOAD", bg = "lightgray",command=lambda: self.button_func5(),height = 1,width = 4)
		b4.grid(row=7, column=8,columnspan=2,sticky='w')		

		self.l1 = tk.Label(self.root, text=" Wdia[kJ]",justify='center')
		self.l1.grid(row=8, column=0,columnspan=3,sticky='w')
		self.e6 = tk.Entry(self.root,width=6,justify='center')
		self.e6.insert(10,self.StrVar6.get())
		self.e6.grid(row=8, column=2,columnspan=2)

		self.l1 = tk.Label(self.root, text="BSmulti",justify='center')
		self.l1.grid(row=8, column=4,columnspan=2)
		self.e13 = tk.Entry(self.root,width=5,justify='center')
		self.e13.insert(7,self.StrVar13.get())
		self.e13.grid(row=8, column=6,columnspan=2)

		self.l1 = tk.Label(self.root, text="Zeff",justify='center')
		self.l1.grid(row=8, column=8,columnspan=2)
		self.e7 = tk.Entry(self.root,width=6,justify='center')
		self.e7.insert(10,self.StrVar7.get())
		self.e7.grid(row=8, column=10,columnspan=2)

		self.l1 = tk.Label(self.root, text="Zimp",justify='center')
		self.l1.grid(row=9, column=0,columnspan=3)
		self.e14 = tk.Entry(self.root,width=6,justify='center')
		self.e14.insert(10,self.StrVar14.get())
		self.e14.grid(row=9, column=2,columnspan=2)

		self.l1 = tk.Label(self.root, text="NL[(19)/m3]",justify='center')
		self.l1.grid(row=9, column=4,columnspan=2)
		self.e8 = tk.Entry(self.root,width=6,justify='center')
		self.e8.insert(10,self.StrVar8.get())
		self.e8.grid(row=9, column=6,columnspan=2)	

		self.l1 = tk.Label(self.root, text="BStype",justify='center')
		self.l1.grid(row=9, column=8,columnspan=2)
		self.m3 = tk.OptionMenu(self.root,self.MenuVar3,'neo','hager','csauter','sauter')
		self.m3.config(width=5)
		self.m3.grid(row=9, column=8,columnspan=4,sticky='e')					

		self.l1 = tk.Label(self.root, text="====================  PROFILES  =====================",justify='center')
		self.l1.grid(row=11, column=0,columnspan=12)		

		self.l1 = tk.Label(self.root, text="Use_preK",justify='center')
		self.l1.grid(row=12, column=0,columnspan=3)
		self.c1 = tk.Checkbutton(self.root,variable=self.CheckVar1,command=lambda: self.check_func1())
		self.c1.grid(row=12, column=3)

		self.l1 = tk.Label(self.root, text="Load[#] ",justify='center')
		self.l1.grid(row=12, column=4,columnspan=2,sticky='e')		
		self.m4 = tk.OptionMenu(self.root,self.MenuVar4,*self.run_index2)
		self.m4.config(width=2)
		self.m4.grid(row=12, column=6,columnspan=2)	
		b5 = tk.Button(self.root, text="GFIT", bg = "lightgray",command=lambda: self.button_func6(),height = 1,width = 4)
		b5.grid(row=12, column=8,columnspan=2,sticky='w')

		b5 = tk.Button(self.root, text="CHECK", bg = "lightgray",command=lambda: self.button_func5a(),height = 1,width = 4)
		b5.grid(row=7, column=10,columnspan=2,sticky='w')

		self.l1 = tk.Label(self.root, text=" Wmhd[kJ]",justify='center')
		self.l1.grid(row=13, column=0,columnspan=3,sticky='w')
		self.e9 = tk.Entry(self.root,width=5,justify='center')
		self.e9.insert(0,self.StrVar9.get())
		self.e9.grid(row=13, column=2,columnspan=2)
		self.e9.configure(state='readonly')

		self.l1 = tk.Label(self.root, text="Wkin[kJ] ",justify='center')
		self.l1.grid(row=13, column=4,columnspan=2,sticky='e')
		self.e10 = tk.Entry(self.root,width=5,justify='center')
		self.e10.insert(0,self.StrVar10.get())
		self.e10.grid(row=13, column=6,columnspan=2)
		self.e10.configure(state='readonly')	

		self.l1 = tk.Label(self.root, text="NL[19/m3]",justify='center')
		self.l1.grid(row=13, column=8,columnspan=3,sticky='w')
		self.e11 = tk.Entry(self.root,width=5,justify='center')
		self.e11.insert(0,self.StrVar11.get())
		self.e11.grid(row=13, column=10,columnspan=2)
		self.e11.configure(state='readonly')

		self.l1 = tk.Label(self.root, text="================   Pressure & Current   =================",justify='center')
		self.l1.grid(row=14, column=0,columnspan=12)

		self.l1 = tk.Label(self.root, text="RTYPE",justify='center')
		self.l1.grid(row=15, column=0,columnspan=2)
		self.m5 = tk.OptionMenu(self.root,self.MenuVar5,'sMSE','eMSE',command=lambda value: find_bdiff_index(self,1))
		self.m5.config(width=4)
		self.m5.grid(row=15, column=2,columnspan=2)			

		self.l1 = tk.Label(self.root, text="BND_PSI ",anchor='s')	
		self.l1.grid(row=15, column=4,columnspan=2,sticky="e")
		self.s1 = tk.Scale(self.root, variable=self.DoubleVar1, command=lambda x: self.scroll_func1(), orient='horizontal', showvalue=True, from_=0.980,to=0.999,resolution=0.001,length=200)
		self.s1.grid(row=15,column=6,columnspan=6,sticky='n')			

		self.l1 = tk.Label(self.root, text="Adjust_Prof",justify='center')
		self.l1.grid(row=16, column=0,columnspan=3)
		self.c2 = tk.Checkbutton(self.root,variable=self.CheckVar2)
		self.c2.grid(row=16, column=3)	

		self.l1 = tk.Label(self.root, text="B-Diff",justify='center')
		self.l1.grid(row=16, column=4,columnspan=2)
		self.e12 = tk.Entry(self.root,width=5,justify='center')
		self.e12.insert(10,self.StrVar12.get())
		self.e12.grid(row=16, column=6,columnspan=2)

		b6 = tk.Button(self.root, text="RUN OPTION", bg = "lightgray",command=lambda: self.button_func7a(),height = 1,width = 10)
		b6.grid(row=16, column=8,columnspan=4)

		b7 = tk.Button(self.root, text="RUN CHEASE", bg = "lightgray",command=lambda: self.button_func7b(),height = 1,width = 10)
		b7.grid(row=17, column=1,columnspan=4)

		b8 = tk.Button(self.root, text="OPEN RESULT", bg = "lightgray",command=lambda: self.button_func7c(),height = 1,width = 10)
		b8.grid(row=17, column=5,columnspan=4)				

		self.l1 = tk.Label(self.root, text="================   GEN EFIT Constraint  =================",justify='center')
		self.l1.grid(row=18, column=0,columnspan=12)

		self.l1 = tk.Label(self.root, text="B-Diff",justify='center')
		self.l1.grid(row=19, column=0,columnspan=2)
		self.m6 = tk.OptionMenu(self.root,self.MenuVar6,*self.bdiff_index)
		self.m6.config(width=4)
		self.m6.grid(row=19, column=2,columnspan=2)	

		self.l1 = tk.Label(self.root, text="MSE",justify='center')
		self.l1.grid(row=19, column=4,columnspan=2)
		self.m7 = tk.OptionMenu(self.root,self.MenuVar7,'sMSE','eMSE',command=lambda value: find_bdiff_index(self,2,True))
		self.m7.config(width=4)
		self.m7.grid(row=19, column=6,columnspan=2)	

		self.l1 = tk.Label(self.root, text="Mode",justify='center')
		self.l1.grid(row=19, column=8,columnspan=2)
		self.m8 = tk.OptionMenu(self.root,self.MenuVar8,'Hmode','Lmode')
		self.m8.config(width=4)
		self.m8.grid(row=19, column=10,columnspan=2)	

		self.l1 = tk.Label(self.root, text="Iter [#]",justify='center')
		self.l1.grid(row=20, column=0,columnspan=2)
		self.m9 = tk.OptionMenu(self.root,self.MenuVar9,*self.iter_index)
		self.m9.config(width=4)
		self.m9.grid(row=20, column=2,columnspan=2)	

		b8 = tk.Button(self.root, text="GENERATE EFIT INPUT", bg = "lightgray",command=lambda: self.button_func8a(),height = 1,width = 21)
		b8.grid(row=20, column=4,columnspan=6)			

		self.l1 = tk.Label(self.root, text="=================        EFIT RUNS       ==================",justify='center')
		self.l1.grid(row=21, column=0,columnspan=12)	

		b9 = tk.Button(self.root, text="RUN EFIT", bg = "lightgray",command=lambda: self.button_func9a(),height = 2,width = 35)
		b9.grid(row=22, column=0,columnspan=12,pady=5)		

		self.l1 = tk.Label(self.root, text="====================     NOTES      =====================",justify='center')
		self.l1.grid(row=25, column=0,columnspan=12)
		
		self.textpad = tkst.ScrolledText(self.root,width=50,height=8)
		self.textpad.grid(row=26,column=0,columnspan=12)

		self.l1 = tk.Label(self.root, text="======================================================",justify='center')
		self.l1.grid(row=27, column=0,columnspan=12)		

		b9 = tk.Button(self.root, text="SAVE", bg = "lightgray",command=lambda: self.button_func10a(),height = 1,width = 10)
		b9.grid(row=28, column=0,columnspan=5,sticky='e')	

		b10 = tk.Button(self.root, text="EXIT", bg = "lightgray",command=lambda: self.button_func10b(),height = 1,width = 10)
		b10.grid(row=28, column=7,columnspan=5,sticky='w')			

		self.l101 = tk.Label(self.root, text="NONE",fg='red',anchor='center',bg='white',width=7,justify='left')
		self.l101.grid(row=12, column=10,columnspan=2)		
		self.l103 = tk.Label(self.root, text="NONE",fg='red',anchor='center',bg='white',width=7,justify='left')
		self.l103.grid(row=17, column=10,columnspan=2)		
		self.l104 = tk.Label(self.root, text="NONE",fg='red',anchor='center',bg='white',width=7,justify='left')
		self.l104.grid(row=20, column=10,columnspan=2)

		font = tkinter.font.Font(weight='bold',size=10)#,font='Courier')
		self.l105 = tk.Label(self.root, text="-Ver %s "%(version['gefit']),anchor='w',width=20,justify='left',fg='dodgerblue',font=font)
		self.l105.grid(row=46, column=0,columnspan=12,sticky='sw')
	
		self.l106 = tk.Label(self.root, text="-%s"%(author['gefit']),anchor='w',width=56,justify='left',fg='dodgerblue',font=font)
		self.l106.grid(row=47, column=0,columnspan=12,sticky='sw')

#		self.l107 = tk.Label(self.root, text="-%s"%(author['gefit']),anchor='w',width=56,justify='left',fg='dodgerblue',font=font)
#		self.l107.grid(row=48, column=0,columnspan=12,sticky='sw')
		
		self.l108 = tk.Label(self.root, text="-%s"%(comment['gefit']),anchor='w',width=50,justify='left',fg='magenta',font=font)
		self.l108.grid(row=49, column=0,columnspan=12,sticky='sw')


		self.firstrun = False
		self.button_func1b(2)
		self.bdiff_index = find_bdiff_index(self,1)
		self.update_status()
		self.root.mainloop()
		return

	def initialise_vars(self,use_last_index=False):

		self.firstrun = True
		self.iseqdsk = False
		self.iskin = False
		self.iskfile = False
		self.isrfile = False
		self.iskfile_already = False
		self.change_fitopt = False
		self.ischease = False
		self.ismse = False
		self.isrefit = False
		self.rzbdy = np.zeros(2)
		self.bdiff_index = np.array([''])
		self.iter_index =  np.array([''])
		self.efit_first_plot = True
		self.is_sgam = False
		self.END = tk.END
		
		self.bdyrc = ''
		self.bdyzc = ''

		for i in range(25):
			self.__dict__['t%i_close'%i] = True

		self.run_index = ['--']
		self.run_index2 = ['--']
		self.efit_index = ['--']
		self.index = 0  
		self.run_index, self.run_index2, self.last_index = find_run_index(self.run_index, self.run_index2,use_last_index)

		self.shot = ''
		self.time = ''
		self.equ_name = None
		self.target_bnd = 0.997
		self.kin_name = None
		self.kfile_name = None

		self.wdia = 280
		self.zeff = 2.0
		self.zimp = 6.0
		self.lineden = 0.0

		self.bs_model = 'csauter'
		self.use_prev = False

		self.wmhd = 0.0
		self.lineden2 = 0.0
		self.wkin = 0.0
		self.rmag = 1.8

		self.run_type = 'sMSE'
		self.beamdiff = 0.0
		self.adjust_prof = False
		self.bsmulti = 1.0


		self.te_file = ''
		self.ne_file = ''
		self.ti_file = ''
		self.vt_file = ''
		self.te_edge_file = ''
		self.ne_edge_file = ''

		#------RUN OPTION (Str100~, Check20~)
		self.zimp = 6		
		self.core_neo = 0.05
		self.vloop_mod = 1.0						#101
		self.ext_vloop = 0.0						#102

		self.hag_core_mod = True
		self.hag_core_mod_psin = 0.3

		self.iternc = 25
		self.iternb = 5
		self.relax = 0.6
		self.ns = 70
		self.nt = 70
		self.map_ns = 200
		self.map_nt = 200
		self.ip_crit = 1.e-4
		self.bs_crit = 1.e-4		

		self.run_step = 20
		self.run_avg = 0
		self.run_dt = 0.01
		self.avg_dt = 0.01
		self.nproc = 28		
		self.mdescrf = ''
		self.sconfigf = ''
		self.iconfigf = ''
		self.use_beam = True

		self.b1ap = 0.	#BEAM option 151
		self.b1bp = 0.
		self.b1cp = 0.
		self.b1ae = 90.
		self.b1be = 90.
		self.b1ce = 90.
		self.b2ap = 0.
		self.b2bp = 0.
		self.b2cp = 0.
		self.b2ae = 90.
		self.b2be = 90.
		self.b2ce = 90.

		self.nbeam = 3
		self.bpower = '0.0,1.5,1.7,0.,0.,0.'		
		self.benergy = '100.,80.,90.,0.,0.,0.'

		self.bdiff_target = ''
		self.mse_target = 'sMSE'
		self.mode_type = 'Hmode'
		self.iter_type = self.iternb

		self.sgam = 0.0
		self.use_smse = False
		self.pnots = ''
		self.ffnots = ''
		self.jnots = ''

		self.s_startp = 0.2
		self.e_knotp = 0.985
		self.s_knotp = 0.3
		self.knot_shiftp = 0.05
		self.corenp = 2
		self.edgenp = 5
		self.inputnp = 501
		self.dmin_corep = 0.1
		self.dmin_edgep = 0.01	

		self.s_startf = 0.2
		self.e_knotf = 0.985
		self.s_knotf = 0.3
		self.knot_shiftf = 0.05
		self.corenf = 1
		self.edgenf = 5
		self.inputnf = 501
		self.dmin_coref = 0.1
		self.dmin_edgef = 0.01

		self.corenj = 14
		self.edgenj = 20
		self.knotsj = 0.7
		self.knotej = 1.0

		self.use_jconst = True
		self.mxiter = 101
		self.efit_relax = 0.5
		self.efit_conv = 1.e-4
		self.fwtcur = 2.
		self.use_q_const = False
		self.use_bnd_const = True
		self.qval = 1.0
		self.push_corep = 0.
		self.mse_shift = 0.

		self.q1 = 1.
		self.q2 = 1.

		self.efitr_n = 0
		self.efitv = '257'
		self.efit_plg1 = []
		self.efit_slg1 = []
		self.efit_plg2 = []
		self.efit_slg2 = []
		self.efit_slg12 = []	

		self.StrVar1.set(self.trans_vars(self.shot,1))
		self.StrVar2.set(self.trans_vars(self.time,1))
		self.StrVar3.set(self.trans_vars(self.equ_name,3))
		self.StrVar4.set(self.trans_vars(self.kin_name,3))
		self.StrVar5.set(self.trans_vars(self.kfile_name,3))
		self.StrVar6.set(self.trans_vars(self.wdia,2))
		self.StrVar7.set(self.trans_vars(self.zeff,2))
		self.StrVar8.set(self.trans_vars(self.lineden,2))
		self.StrVar9.set(self.trans_vars(self.wmhd,2))
		self.StrVar10.set(self.trans_vars(self.wkin,2))	
		self.StrVar11.set(self.trans_vars(self.lineden2,2))		
		self.StrVar12.set(self.trans_vars(self.beamdiff,2))
		self.StrVar13.set(self.trans_vars(self.bsmulti,2))
		self.StrVar14.set(self.trans_vars(self.zimp,2))
		self.StrVar15.set(self.trans_vars(self.mode_type,3))
		self.StrVar16.set(self.trans_vars(self.sgam,2))
		self.StrVar17.set(self.trans_vars(self.pnots,3))
		self.StrVar18.set(self.trans_vars(self.ffnots,3))
		self.StrVar19.set(self.trans_vars(self.jnots,3))

		self.StrVar20.set(self.trans_vars(self.s_startp,2))
		self.StrVar21.set(self.trans_vars(self.e_knotp,2))
		self.StrVar22.set(self.trans_vars(self.s_knotp,2))
		self.StrVar23.set(self.trans_vars(self.knot_shiftp,2))
		self.StrVar24.set(self.trans_vars(self.corenp,1))
		self.StrVar25.set(self.trans_vars(self.edgenp,1))
		self.StrVar26.set(self.trans_vars(self.inputnp,1))

		self.StrVar27.set(self.trans_vars(self.s_startf,2))
		self.StrVar28.set(self.trans_vars(self.e_knotf,2))
		self.StrVar29.set(self.trans_vars(self.s_knotf,2))
		self.StrVar30.set(self.trans_vars(self.knot_shiftf,2))
		self.StrVar31.set(self.trans_vars(self.corenf,1))
		self.StrVar32.set(self.trans_vars(self.edgenf,1))
		self.StrVar33.set(self.trans_vars(self.inputnf,1))		

		self.StrVar34.set(self.trans_vars(self.corenj,1))
		self.StrVar35.set(self.trans_vars(self.edgenj,1))
		self.StrVar36.set(self.trans_vars(self.knotsj,2))
		self.StrVar37.set(self.trans_vars(self.knotej,2))
		self.StrVar46.set(self.trans_vars(self.dmin_corep,2))
		self.StrVar47.set(self.trans_vars(self.dmin_edgep,2))
		self.StrVar48.set(self.trans_vars(self.dmin_coref,2))
		self.StrVar49.set(self.trans_vars(self.dmin_edgef,2))
		
		self.StrVar50.set('')

		self.StrVar51.set(self.bdyrc)
		self.StrVar52.set(self.bdyzc)
	
		self.StrVar38.set(self.trans_vars(self.mxiter,1))
		self.StrVar39.set(self.trans_vars(self.efit_relax,2))
		self.StrVar40.set('%3.1e'%self.efit_conv)
		self.StrVar41.set(self.trans_vars(self.fwtcur,2))
		self.StrVar42.set(self.trans_vars(self.qval,2))
		self.StrVar43.set(self.trans_vars(self.push_corep,2))
		self.StrVar44.set(self.trans_vars(self.mse_shift,2))		
		self.StrVar45.set(self.trans_vars(self.vt_file,3))
				
		if use_last_index: self.MenuVar1.set(self.last_index)
		else:	self.MenuVar1.set(self.run_index[-1])
		self.MenuVar2.set(self.run_index2[-1])
		self.MenuVar3.set(self.bs_model)
		self.MenuVar4.set(self.run_index2[-1])
		self.MenuVar5.set(self.trans_vars(self.run_type,3))
		self.MenuVar6.set(self.trans_vars(self.bdiff_target,2))
		self.MenuVar7.set(self.mse_target)
		self.MenuVar8.set(self.mode_type)
		self.MenuVar9.set(self.trans_vars(self.iter_type,2))
		self.MenuVar10.set('--')
		self.MenuVar13.set('257')
		self.MenuVar15.set('1')
		self.MenuVar16.set('0')

		self.CheckVar1.set(self.trans_vars(self.use_prev,4))
		self.CheckVar2.set(self.trans_vars(self.adjust_prof,4))
		self.CheckVar3.set(self.trans_vars(self.use_smse,4))
		self.CheckVar4.set(self.trans_vars(self.use_jconst,4))
		self.CheckVar5.set(self.trans_vars(self.use_q_const,4))
		self.CheckVar6.set(self.trans_vars(self.use_bnd_const,4))
		self.CheckVar7.set(0)

		self.DoubleVar1.set(self.target_bnd)
		self.DoubleVar2.set(self.q1)
		self.DoubleVar3.set(self.q2)

		self.StrVar99.set(version['gefit'])

		self.StrVar101.set(self.trans_vars(self.core_neo,2))
		self.StrVar102.set(self.trans_vars(self.vloop_mod,2))					
		self.StrVar103.set(self.trans_vars(self.ext_vloop,2))
		self.StrVar104.set(self.trans_vars(self.hag_core_mod_psin,2))
		self.StrVar105.set(self.trans_vars(self.iternc,1))
		self.StrVar106.set(self.trans_vars(self.iternb,1))
		self.StrVar107.set(self.trans_vars(self.relax,2))
		self.StrVar108.set(self.trans_vars(self.ns,1))
		self.StrVar109.set(self.trans_vars(self.nt,1))
		self.StrVar110.set(self.trans_vars(self.map_ns,1))
		self.StrVar111.set(self.trans_vars(self.map_nt,1))
		self.StrVar112.set('%e'%self.ip_crit)
		self.StrVar113.set('%e'%self.bs_crit)

		self.StrVar114.set(self.trans_vars(self.run_step,1))
		self.StrVar115.set(self.trans_vars(self.run_avg,1))
		self.StrVar116.set(self.trans_vars(self.run_dt,2))
		self.StrVar117.set(self.trans_vars(self.avg_dt,2))		
		self.StrVar118.set(self.trans_vars(self.nproc,1))	
		self.StrVar119.set(self.mdescrf)
		self.StrVar120.set(self.sconfigf)
		self.StrVar121.set(self.iconfigf)

		self.CheckVar21.set(self.trans_vars(self.hag_core_mod,4))
		self.CheckVar22.set(self.trans_vars(self.use_beam,4))	
		self.CheckVar23.set(self.trans_vars(self.iskfile_already,4))
		self.CheckVar24.set(self.trans_vars(self.ismse,4))	

		self.StrVar151.set(self.trans_vars(self.b1ap,2))
		self.StrVar152.set(self.trans_vars(self.b1ae,2))
		self.StrVar153.set(self.trans_vars(self.b1bp,2))
		self.StrVar154.set(self.trans_vars(self.b1be,2))
		self.StrVar155.set(self.trans_vars(self.b1cp,2))
		self.StrVar156.set(self.trans_vars(self.b1ce,2))
		self.StrVar157.set(self.trans_vars(self.b2ap,2))
		self.StrVar158.set(self.trans_vars(self.b2ae,2))
		self.StrVar159.set(self.trans_vars(self.b2bp,2))
		self.StrVar160.set(self.trans_vars(self.b2be,2))
		self.StrVar161.set(self.trans_vars(self.b2cp,2))
		self.StrVar162.set(self.trans_vars(self.b2ce,2))


		self.MenuVar11.set('--')
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

	def update_status(self):

		if len(self.bdiff_index)>0.:	self.ischease = True
		else:	self.ischease = False

		if (self.ischease and self.iskfile and self.isrfile):	self.isrefit = True
		else:	self.isrefit = False

		if (self.iseqdsk and self.iskin and self.iskfile and self.isrfile):	self.l101.config(text='READY',fg='lime')
		else:	self.l101.config(text='NONE',fg='red')

		prof_dir = self.currdir + '/%s/PROFILE/PROFILES/chease_kinprof.out'%self.MenuVar1.get()
		#if os.path.isfile(prof_dir):	self.l102.config(text='READY',fg='lime')
		#else:	self.l102.config(text='NONE',fg='red')

		if self.ischease:	self.l103.config(text='READY',fg='lime')
		else:	self.l103.config(text='NONE',fg='red')		

		if self.isrefit:	self.l104.config(text='READY',fg='lime')
		else:	self.l104.config(text='NONE',fg='red')		

		return

	def __init__(self):

		self.currdir = os.getcwd()

		return

if __name__ == "__main__":

	import gefit
	gefitk = gefit.gefitk()
	gefitk.root = tk.Tk()
	gefitk.root.title('G-EFITK')
	gefitk.ts = True
	gefitk.ts_run = False
#	try:
#		if (sys.argv[1].lower() == '-plare' or sys.argv[1].lower() == '-test') : 
#			gefitk.ts = True
#			gefitk.ts_run = False
#			print('>>> START TEST MODE')
#	except: pass
	gefitk.declare_vars()
	gefitk.gui_efit()
