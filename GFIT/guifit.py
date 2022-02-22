#!/usr/local/anaconda3/bin/python3
import os,sys
import tkinter as tk
from tkinter import ttk
import numpy as np
import fittool
from shutil import move, copyfile, copytree, rmtree
from tkinter.filedialog import askopenfilename,asksaveasfilename
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import matplotlib.pyplot as plt
import copy
import time
from get_efit import *
from gfit_ece_multi import *
from progress import *
from exec_dirs import version,author,comment,mds_dir,mds_tci,pythonc_exec,mds_over,mds_ref,mds_ts
from exec_dirs import mds_ces,mds_lit,mse_corr,mse_dir,python2_exec,gzip_dir,mds_da,dena_dir
import pickle

currdir = os.getcwd()

class guifittool:

	def gfit(self):
		self.fit = fittool.fit_tool()

		if os.path.isfile('fit_opt.save'):
			f = open('fit_opt.save','rb')
			fit_opt =pickle.load(f)
			self.load_fitopt_dict(fit_opt)
			f.close()	

		if not self.fit.fit_opt['file']['gfile'] == None:
			self.fit.post['eq_old'] = self.fit.fit_opt['file']['gfile']
			self.fit.input_file_check()
			self.fit.adjust_variables()				
			self.fit.mapping_variables()
			self.fit.post['eq_change'] = False
			self.fit.read_kinprof_kprofile()
			self.fit.read_raw_kinprof()	
		
		self.root.protocol('WM_DELETE_WINDOW')
		self.declare_variables()
		self.initialise_variables()
		if self.single_only: self.read_gfitp()
		if self.multi_only: self.read_multi_fitopt()
		self.get_fit_opt()
		print('>>> MDS directory is set to "%s"'%self.mds_dir)

		if os.path.isfile('fit_opt.save_param'):
			print('>>> Load recent save param')
			f = open('fit_opt.save_param','rb')
			temp = pickle.load(f)
			f.close()
			for flag in self.fit.prof_list:
				for key in temp[flag].keys(): self.post_opt['param'][flag][key] = copy.deepcopy(temp[flag][key])
		if os.path.isfile('fit_opt.save_param_m'):
			f = open('fit_opt.save_param_m','rb')
			self.post_opt['param_m'] = pickle.load(f)
			f.close()
			

		self.l1 = tk.Label(self.root, text="============ KSTAR Profile fitting tool ============",justify='right')
		self.l1.grid(row=0,column=0,columnspan=4)

		self.make_home_canvas()
		self.make_note_frame()
		self.pedestal_option()
		self.fitting_func_option()
		
		self.input_option()
		endline = dict()
		for k in range(4): 
			flag = self.fit.prof_list[k]
			self.fit_param_option(self.home['page%i'%(k+3)],flag)
			endline[flag] = self.fit_constraint_option(self.home['page%i'%(k+3)],13,flag)
			endline[flag] = self.fit_variable_option(self.home['page%i'%(k+3)],endline[flag]+1,flag)

		self.fit_inter(self.home['page4'],endline['ne']+1,'ne')

		self.fit_fuction_plot(14)
		self.make_mds()
		self.make_plot()
		self.make_etc()

		self.make_plot_canvas()

		for k in range(4): self.root.grid_columnconfigure(k,weight=1)
		b1 = tk.Button(self.root,  text="DO FIT", height = 2, width = 8,command=lambda: self.run_fit())
		b1.grid(row=2, column=0)
		self.b2 = tk.Button(self.root,  text="SAVE-S",  height = 2,width = 8,command=lambda: self.run_save())
		self.b2.grid(row=2, column=2)
		b3 = tk.Button(self.root,  text="EXIT",  height = 2,width = 8,command=lambda: self.run_exit())
		b3.grid(row=2, column=3)
		b1b = tk.Button(self.root, text="SYNC FIT", height = 2,width = 8,command=lambda: self.make_mds_f())
		b1b.grid(row=2,column=1)

		if self.single_only: b1b.configure(state='disabled')

		b4 = tk.Button(self.root, text="Derivative", height = 1, width = 7,command=lambda: self.make_plot_b())	
		b4.grid(row=0, column=4,padx=70,sticky='e')		

		b5 = tk.Button(self.root, text="Restore", height = 1, width = 5,command=lambda: self.make_plot_a())	
		b5.grid(row=0, column=4,padx=5,sticky='e')		

		b5 = tk.Button(self.root, text="Back", height = 1, width = 5,command=lambda: self.back_())	
		b5.grid(row=0, column=4,padx=270,sticky='e')	

		b5 = tk.Button(self.root, text="Next", height = 1, width = 5,command=lambda: self.next_())	
		b5.grid(row=0, column=4,padx=200,sticky='e')					

		b6 = tk.Button(self.root, text="Restore Opt", height = 1, width = 9,command=lambda: self.restore_fitopt())	
		b6.grid(row=0, column=4,sticky='w')		

		l1 = tk.Label(self.root, text='Ver. %s'%version['gfit'],fg='magenta')
		l1.grid(row=1,column=3,sticky='en')

		self.initialise_setup()
		self.check_runmode()

		if os.path.isfile('fit_opt.save_param'):
			for flag in self.fit.prof_list: self.fitting_func_option_a(flag)
	
		self.root.resizable(0,0)
		if self.multi_only:
			self.fit.forcefit2= True
			self.fit.hmodefit = True
			if os.path.isfile('result_multi.save'): b1.invoke() 

		self.root.mainloop()

		return

	def make_home_canvas(self):

		self.figure['name']['home'] = plt.figure(1,figsize=self.figure['size']['home'])
		ax1 = self.figure['name']['home'].add_subplot(2,2,1)
		ax2 = self.figure['name']['home'].add_subplot(2,2,2)
		ax3 = self.figure['name']['home'].add_subplot(2,2,3)
		ax4 = self.figure['name']['home'].add_subplot(2,2,4)

		axes = self.figure['name']['home'].axes
		titles = ['$T_e$ [keV]','$n_e$ [$10^{19}$/m3]','$T_i$ [keV]','$V_\phi$ [km/s]']

		for i in range(4):
			ax = axes[i]
			flag = self.fit.prof_list[i]
			
			ax.set_title(titles[i])
			if self.fit.fit_opt['use_rho'][flag]: ax.set_xlabel('$\psi_n$ [a.u]')
			else: ax.set_xlabel('$\\rho_t$ [a.u]')
			ax.set_ylabel(titles[i])

		self.figure['name']['home'].tight_layout()
		self.figure['canvas']['home'] = FigureCanvasTkAgg(self.figure['name']['home'],master=self.root)

		self.figure['widget']['home'] = self.figure['canvas']['home'].get_tk_widget()
		self.figure['widget']['home'].grid(row=1,column=4,rowspan=2,sticky='w')

		toolbar_frame = tk.Frame(self.root)
		toolbar_frame.grid(row=0,column=4)
		
		self.figure['toolbar']['home'] = NavigationToolbar2Tk(self.figure['canvas']['home'],toolbar_frame)		
		return

	def make_plot_canvas(self):
		self.figure['name']['plot'] = plt.figure(0,figsize=self.figure['size']['plot'])
		ax1 = self.figure['name']['plot'].add_subplot(1,1,1)

		self.figure['canvas']['plot'] = FigureCanvasTkAgg(self.figure['name']['plot'],master=self.home['page2'])
		self.figure['widget']['plot'] = self.figure['canvas']['plot'].get_tk_widget()
		self.figure['widget']['plot'].grid(row=16,column=0,columnspan=9)
		[ax] = self.figure['name']['plot'].axes
		ax.set_facecolor('b')
		ax.axis('off')
		return
		
	def leftclick(self,event):
		page = self.home['nb'].tk.call(self.home['nb']._w,'identify','tab',event.x,event.y)
		if not page == '':
			if not self.curr_page == page:
				self.curr_page = page
				flag = 'ne'
				self.sync_gui_opt()				
				self.put_fit_opt()							
				self.open_page()
				
		return

	def make_note_frame(self):

		titles = ['INPUT','FUNC',' TE ',' NE ',' TI ',' VT ','MDS','PLOT','ETC']
		self.home['nb'] = ttk.Notebook(self.root,width=400,height=710) #665
		self.home['nb'].bind('<Button-1>',self.leftclick)
		for i in range(1,len(titles)+1):
			self.home['page%i'%i] = ttk.Frame(self.home['nb'])
			self.home['nb'].add(self.home['page%i'%i],text=titles[i-1])
			for k in range(9): self.home['page%i'%i].grid_columnconfigure(k,weight=1)

		self.home['nb'].grid(row=1,column=0,columnspan=4,sticky='n')

		return

	def find_filename(self,entry,readkin=False,readeq=False):

		inputf = askopenfilename()
		if (len(inputf) == 0):      return		
		self.note_in['e%i'%entry].delete(0,'end')
		self.note_in['e%i'%entry].insert(10,inputf)

		if (readkin or readeq):
			self.sync_gui_opt()
			self.put_fit_opt()
			if readeq: 
				self.fit.post['eq_old'] = self.fit.fit_opt['file']['gfile']
				self.fit.mapping_variables()
				self.fit.post['eq_change'] = True

			if readkin: self.fit.read_raw_kinprof()

		return

	def input_option(self):

		self.l1 = tk.Label(self.home['page1'], text="=============== EQ INPUT FILES ================",justify='center')
		self.l1.grid(row=0,column=0,columnspan=9,pady=5)

		self.l1 = tk.Label(self.home['page1'], text='g-file')
		self.l1.grid(row=1,column=0)

		self.note_in['e1'] = tk.Entry(self.home['page1'],width=35,justify='center')
		self.note_in['e1'].grid(row=1, column=1,columnspan=7)
		self.note_in['e1'].insert(10,self.note_in['file']['gfile'].get())	

		self.note_in['b1'] = tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(1,False,True))
		self.note_in['b2'] = tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(2,True))
		self.note_in['b3'] = tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(3,True))
		self.note_in['b4'] = tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(4,True))
		self.note_in['b5'] = tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(5,True))
		self.note_in['b6'] = tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(6,True))
		self.note_in['b7'] = tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(7,True))
		self.note_in['b8'] = tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(8,True))
		self.note_in['b9'] = tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(9,True))
		self.note_in['b19']= tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(19))
		self.note_in['b20']= tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(20))
		self.note_in['b21']= tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(21))
		self.note_in['b22']= tk.Button(self.home['page1'],width=3,text='OPEN',command=lambda: self.find_filename(22))
		
		self.note_in['b1'].grid(row=1, column=8)			
		count1 = 2;	count2 = 1;
		for flag in self.fit.prof_list:
			self.l1 = tk.Label(self.home['page1'], text="=============== %s INPUT FILES ================"%flag.upper(),justify='center')
			self.l1.grid(row=count1,column=0,columnspan=9,pady=5)
			count1 = count1 + 1;	
			for k in self.fit.__dict__['%s_list'%flag]:				
				count2 = count2 + 1;

				self.l1 = tk.Label(self.home['page1'], text='%s'%(k.upper()))
				self.l1.grid(row=count1,column=0)

				self.note_in['e%i'%count2] = tk.Entry(self.home['page1'],width=35,justify='center')
				self.note_in['e%i'%count2].grid(row=count1, column=1,columnspan=7)
				self.note_in['e%i'%count2].insert(10,self.note_in['file'][flag][k].get())				
				self.note_in['b%i'%count2].grid(row=count1, column=8)			

				count1 = count1 + 1;

		inter_list = ['int01','int02','tci01','tci02','tci03','tci04','tci05']
		self.l1 = tk.Label(self.home['page1'], text="============== INTERF PARARMS ===============",justify='center')
		self.l1.grid(row=count1,column=0,columnspan=9,pady=5)

		for i in range(0,4):
			flag= inter_list[i]
			count2 = count2 + 1;
			self.l1 = tk.Label(self.home['page1'], text=flag.upper())
			self.l1.grid(row=count1+1,column=2*i)			
			self.note_in['e%i'%count2] = tk.Entry(self.home['page1'],width=5,justify='center')
			self.note_in['e%i'%count2].grid(row=count1+1, column=2*i+1,sticky='w',pady=3)
			self.note_in['e%i'%count2].insert(10,self.note_in[flag].get())	

		for i in range(0,3):
			flag= inter_list[i+4]
			count2 = count2 + 1;
			self.l1 = tk.Label(self.home['page1'], text=flag.upper())
			self.l1.grid(row=count1+2,column=2*i)			
			self.note_in['e%i'%count2] = tk.Entry(self.home['page1'],width=5,justify='center')
			self.note_in['e%i'%count2].grid(row=count1+2, column=2*i+1,sticky='w',pady=3)
			self.note_in['e%i'%count2].insert(10,self.note_in[flag].get())	

		count1 = count1 + 2

		self.l1 = tk.Label(self.home['page1'], text="=============== INPUT PARARMS ================",justify='center')
		self.l1.grid(row=count1+1,column=0,columnspan=9,pady=5)

		self.l1 = tk.Label(self.home['page1'], text='ZEFF')
		self.l1.grid(row=count1+2,column=0); count2 += 1
		self.note_in['e%i'%count2] = tk.Entry(self.home['page1'],width=5,justify='center')
		self.note_in['e%i'%count2].grid(row=count1+2, column=1,sticky='w')
		self.note_in['e%i'%count2].insert(10,self.note_in['zeff'].get())	

		self.l1 = tk.Label(self.home['page1'], text='ZIMP')
		self.l1.grid(row=count1+2,column=2); count2 += 1	
		self.note_in['e%i'%count2] = tk.Entry(self.home['page1'],width=5,justify='center')
		self.note_in['e%i'%count2].grid(row=count1+2, column=3,sticky='w')
		self.note_in['e%i'%count2].insert(10,self.note_in['zimp'].get())

		self.l1 = tk.Label(self.home['page1'], text='AMAIN')
		self.l1.grid(row=count1+2,column=4); count2 += 1	
		self.note_in['e%i'%count2] = tk.Entry(self.home['page1'],width=5,justify='center')
		self.note_in['e%i'%count2].grid(row=count1+2, column=5,sticky='w')
		self.note_in['e%i'%count2].insert(10,self.note_in['amain'].get())

		self.l1 = tk.Label(self.home['page1'], text='AIMP')
		self.l1.grid(row=count1+2,column=6); count2 += 1	
		self.note_in['e%i'%count2] = tk.Entry(self.home['page1'],width=5,justify='center')
		self.note_in['e%i'%count2].grid(row=count1+2, column=7,sticky='w')
		self.note_in['e%i'%count2].insert(10,self.note_in['aimp'].get())				

		count1 = count1 + 3
		self.l1 = tk.Label(self.home['page1'], text="=============== KPROFILE FILES ================",justify='center')
		self.l1.grid(row=count1,column=0,columnspan=9,pady=5)

		count2 += 1
		for i in range(4):
			flag = self.fit.prof_list[i]
			count1 = count1 + 1
			self.l1 = tk.Label(self.home['page1'], text='%s'%(flag.upper()))
			self.l1.grid(row=count1,column=0)

			self.note_in['e%i'%(i+count2)] = tk.Entry(self.home['page1'],width=35,justify='center')
			self.note_in['e%i'%(i+count2)].grid(row=count1, column=1,columnspan=7)
			self.note_in['e%i'%(i+count2)].insert(10,self.note_in['file']['kfile'][flag].get())
	
			self.note_in['b%i'%(i+19)].grid(row=count1, column=8)			


		b1 = tk.Button(self.home['page1'], text="LOAD OPT", height = 1, width = 10,command=lambda: self.load_fitopt())
		b1.grid(row=count1+2, column=0,columnspan=3,pady=10)
		b2 = tk.Button(self.home['page1'], text="LOAD PARAM", height = 1,width = 10,command=lambda: self.load_fitopt_param())
		b2.grid(row=count1+2, column=3,columnspan=3,pady=10)
		b3 = tk.Button(self.home['page1'], text="RESET OPT", height = 1,width = 10,command=lambda: self.reset_fitopt())
		b3.grid(row=count1+2, column=6,columnspan=3,pady=10)

		return

	def pedestal_option(self):

		self.l1 = tk.Label(self.home['page2'], text="=================== Pedestal option ===================",justify='center')
		self.l1.grid(row=0, column=0,columnspan=9,pady=5)
		 
		self.l2 = tk.Label(self.home['page2'], text="  ped_scan_fit")
		self.note_fn['c1'] = tk.Checkbutton(self.home['page2'],variable=self.ped_scan_fit, command=lambda: self.pedestal_option_a())
		self.l2.grid(row=1, column=0, columnspan=2)
		self.note_fn['c1'].grid(row=1, column=2)
		
		self.l2 = tk.Label(self.home['page2'], text="  TI_ped_scan")
		self.note_fn['c2'] = tk.Checkbutton(self.home['page2'],variable=self.ti_ped_scan, command=lambda: self.pedestal_option_a())
		self.l2.grid(row=1, column=3, columnspan = 2)
		self.note_fn['c2'].grid(row=1, column=5)
		self.note_fn['c2'].configure(state='disabled')

		self.l2 = tk.Label(self.home['page2'], text="  TI_width")
		self.note_fn['c3'] = tk.Checkbutton(self.home['page2'],variable=self.ti_ped_width,command=lambda: self.pedestal_option_a())
		self.l2.grid(row=1, column=6, columnspan = 2)
		self.note_fn['c3'].grid(row=1, column=8)		

		return

	def pedestal_option_a(self):

		if self.ped_scan_fit.get() == 1:	self.note_fn['c2'].configure(state='normal')
		else:	
			self.note_fn['c2'].configure(state='disabled')
			self.ti_ped_scan.set(0)

		if self.ped_scan_fit.get() == 1:
			self.post_opt['func_type']['te'] = self.note_fn['func_type']['te'].get()
			self.post_opt['func_type']['ne'] = self.note_fn['func_type']['ne'].get()

			self.note_fn['func_type']['ne'].set('EPED')
			self.note_fn['func_type']['te'].set('EPED')
			self.fitting_func_option_a('te')
			self.fitting_func_option_a('ne')

			self.note_fn['m1'].config(state='disabled')
			self.note_fn['m2'].config(state='disabled')
		else:
			if not self.post_opt['func_type']['te'] == None: 
				self.note_fn['func_type']['te'].set(self.post_opt['func_type']['te'])
				self.post_opt['func_type']['te'] = None
			if not self.post_opt['func_type']['ne'] == None: 
				self.note_fn['func_type']['ne'].set(self.post_opt['func_type']['ne'])
				self.post_opt['func_type']['ne'] = None
			self.fitting_func_option_a('te')
			self.fitting_func_option_a('ne')

			self.note_fn['m1'].config(state='normal')
			self.note_fn['m2'].config(state='normal')			

		if self.ti_ped_scan.get() == 1:
			self.post_opt['func_type']['ti'] = self.note_fn['func_type']['ti'].get()

			self.note_fn['func_type']['ti'].set('EPED')
			self.note_fn['m3'].config(state='disabled')
			self.fitting_func_option_a('ti')
		else:
			if not self.post_opt['func_type']['ti'] == None: 
				self.note_fn['func_type']['ti'].set(self.post_opt['func_type']['ti'])
				self.post_opt['func_type']['ti'] = None
			self.note_fn['m3'].config(state='normal')
			self.fitting_func_option_a('ti')


		if self.ti_ped_width.get() == 1:
			self.post_opt['width_fix']['te'] = self.note_fn['width_fix']['te'].get()
			self.post_opt['width_fix']['ne'] = self.note_fn['width_fix']['ne'].get()
			self.note_fn['width_fix']['te'].set(1)
			self.note_fn['width_fix']['ne'].set(1)
			self.note_fn['c13'].configure(state='disabled')
			self.note_fn['c14'].configure(state='disabled')

			self.post_opt['width_val']['te'] = self.note_fn['e5'].get()
			self.post_opt['width_val']['ne'] = self.note_fn['e6'].get()
			self.note_fn['e5'].delete(0,'end')
			self.note_fn['e6'].delete(0,'end')
			self.note_fn['e5'].insert(10,'0.0')
			self.note_fn['e6'].insert(10,'0.0')
			self.note_fn['e5'].configure(state='readonly')
			self.note_fn['e6'].configure(state='readonly')

		else:	
			if not self.post_opt['width_val']['te']==None: 
				self.note_fn['e5'].delete(0,'end')
				self.note_fn['e5'].insert(10,self.post_opt['width_val']['te'])
				self.post_opt['width_val']['te']=None
			
			if not self.post_opt['width_val']['ne']==None: 
				self.note_fn['e6'].delete(0,'end')	
				self.note_fn['e6'].insert(10,self.post_opt['width_val']['ne'])
				self.post_opt['width_val']['ne']=None
			

			self.note_fn['e5'].configure(state='normal')
			self.note_fn['e6'].configure(state='normal')

			if not self.post_opt['width_fix']['te'] == None: 
				self.note_fn['width_fix']['te'].set(self.post_opt['width_fix']['te'])
				self.post_opt['width_fix']['te'] = None
			if not self.post_opt['width_fix']['ne'] == None: 
				self.note_fn['width_fix']['ne'].set(self.post_opt['width_fix']['ne'])
				self.post_opt['width_fix']['ne'] = None
			self.note_fn['c13'].configure(state='normal')
			self.note_fn['c14'].configure(state='normal')

		return

	def fitting_func_option(self):

		self.l1 = tk.Label(self.home['page2'], text="================== Fitting Func Option ==================",justify='center')
		self.l1.grid(row=2, column=0,columnspan=9,pady=5)

		titles = [" TE [keV]"," NE [(19)]"," TI [keV]"," VT [km/s]"]

		self.l2 = tk.Label(self.home['page2'], text="RAW")
		self.l2.grid(row=3, column=1)
		self.l2 = tk.Label(self.home['page2'], text="SEP")
		self.l2.grid(row=3, column=2)	
		self.l2 = tk.Label(self.home['page2'], text="VAL")
		self.l2.grid(row=3, column=3)			
		self.l2 = tk.Label(self.home['page2'], text="WIDTH")
		self.l2.grid(row=3, column=4)		
		self.l2 = tk.Label(self.home['page2'], text="VAL")
		self.l2.grid(row=3, column=5)				
		self.l3 = tk.Label(self.home['page2'], text = '   FUNC_TYPE   ',justify='center')
		self.l3.grid(row=3, column=6, columnspan=3)			

		self.note_fn['m1'] = tk.OptionMenu(self.home['page2'],self.note_fn['func_type']['te'],*self.func_list,command=lambda value: self.fitting_func_option_a('te'))
		self.note_fn['m2'] = tk.OptionMenu(self.home['page2'],self.note_fn['func_type']['ne'],*self.func_list,command=lambda value: self.fitting_func_option_a('ne'))
		self.note_fn['m3'] = tk.OptionMenu(self.home['page2'],self.note_fn['func_type']['ti'],*self.func_list,command=lambda value: self.fitting_func_option_a('ti'))
		self.note_fn['m4'] = tk.OptionMenu(self.home['page2'],self.note_fn['func_type']['vt'],*self.func_list,command=lambda value: self.fitting_func_option_a('vt'))

		for k in range(4):
			flag = self.fit.prof_list[k]
			self.l2 = tk.Label(self.home['page2'], text=titles[k])
			self.l2.grid(row=4+k, column=0)
			self.note_fn['c%i'%(k+5)] = tk.Checkbutton(self.home['page2'],variable=self.note_fn['raw_fit'][flag])
			self.note_fn['c%i'%(k+5)].grid(row=4+k, column=1)

			self.note_fn['c%i'%(k+9)] = tk.Checkbutton(self.home['page2'],variable=self.note_fn['sep_fix'][flag])
			self.note_fn['c%i'%(k+9)].grid(row=4+k, column=2)			

			self.note_fn['e%i'%(k+1)] = tk.Entry(self.home['page2'],width=6,justify='center')
			self.note_fn['e%i'%(k+1)].grid(row=4+k, column=3)
			self.note_fn['e%i'%(k+1)].insert(10,self.note_fn['sep_val'][flag].get())

			self.note_fn['c%i'%(k+13)] = tk.Checkbutton(self.home['page2'],variable=self.note_fn['width_fix'][flag])
			self.note_fn['c%i'%(k+13)].grid(row=4+k, column=4)			

			self.note_fn['e%i'%(k+5)] = tk.Entry(self.home['page2'],width=7,justify='center')
			self.note_fn['e%i'%(k+5)].grid(row=4+k, column=5)
			self.note_fn['e%i'%(k+5)].insert(10,self.note_fn['width_val'][flag].get())

			#Fitting Fun type			
			self.note_fn['m%i'%(k+1)].grid(row=4+k,column=6,columnspan=3)
			self.note_fn['m%i'%(k+1)].config(width=7)

		self.l1 = tk.Label(self.home['page2'], text="================= Fitting Option ETC. =================",justify='center')
		self.l1.grid(row=8, column=0,columnspan=9,pady=5)

		self.l1 = tk.Label(self.home['page2'], text="USE_RHO",justify='center')
		self.l1.grid(row=10, column=0)
		self.l1 = tk.Label(self.home['page2'], text="V-SPLINE",justify='center')
		self.l1.grid(row=11, column=0)		

		for k in range(4):
			i = self.fit.prof_list[k]
			self.l1 = tk.Label(self.home['page2'], text="%s"%i.upper(),justify='center')
			self.l1.grid(row=9, column=2*k+1,columnspan=2)			
			self.note_fn['c%i'%(k+17)] = tk.Checkbutton(self.home['page2'],variable=self.note_fn['use_rho'][i])
			self.note_fn['c%i'%(k+17)].grid(row=10, column=2*k+1,columnspan=2,pady=3)		

			self.note_fn['e%i'%(k+9)] = tk.Entry(self.home['page2'],width=7,justify='center')
			self.note_fn['e%i'%(k+9)].grid(row=11, column=2*k+1,columnspan=2)
			self.note_fn['e%i'%(k+9)].insert(10,self.note_fn['smooth'][i].get())

		return

	def fitting_func_option_a(self,flag,isupdate=True):

		func = self.note_fn['func_type'][flag].get().lower()
		if isupdate: self.put_param(flag); self.put_param_gui(flag);
		for i in range(self.fit.func['varn'][func],10):
			self.__dict__['note_%s'%flag]['e%i'%(3*i+1)].configure(state='readonly')
			self.__dict__['note_%s'%flag]['e%i'%(3*i+2)].configure(state='readonly')
			self.__dict__['note_%s'%flag]['e%i'%(3*i+3)].configure(state='readonly')
			self.__dict__['note_%s'%flag]['c%i'%(i+1)].configure(state='disabled')
		return

	def put_param(self,flag):

		func = self.note_fn['func_type'][flag].get().lower()
		if not self.post_opt['param'][flag][func]['is']:
			self.fit.lmfit_init_param(self.fit.func['num'][func],flag )
			self.fit.post['func_history'][flag][self.fit.func['num'][func]] = True
			self.post_opt['param'][flag][func]['val'] = copy.deepcopy(self.fit.param[flag])
			self.post_opt['param'][flag][func]['is'] = True

		return

	def put_param_gui(self,flag):

		func = self.note_fn['func_type'][flag].get().lower()
		count1 = 1; count2 = 1;
		for i in range(10):
			self.__dict__['note_%s'%flag]['e%i'%count1].configure(state='normal')
			self.__dict__['note_%s'%flag]['e%i'%(count1+1)].configure(state='normal')
			self.__dict__['note_%s'%flag]['e%i'%(count1+2)].configure(state='normal')
			self.__dict__['note_%s'%flag]['c%i'%count2].configure(state='normal')	

			self.__dict__['note_%s'%flag]['e%i'%count1].delete(0,'end')
			self.__dict__['note_%s'%flag]['e%i'%(count1+1)].delete(0,'end')
			self.__dict__['note_%s'%flag]['e%i'%(count1+2)].delete(0,'end')

			self.__dict__['note_%s'%flag]['e%i'%count1].insert(10,'%5.3f'%self.post_opt['param'][flag][func]['val']['val'][i])
			self.__dict__['note_%s'%flag]['e%i'%(count1+1)].insert(10,'%5.3f'%self.post_opt['param'][flag][func]['val']['min'][i])
			self.__dict__['note_%s'%flag]['e%i'%(count1+2)].insert(10,'%5.3f'%self.post_opt['param'][flag][func]['val']['max'][i])

			count1 = count1 + 3
			count2 = count2 + 1

			if self.post_opt['param'][flag][func]['val']['vary'][i]:
				self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['fix'].set(0)
			else:
				self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['fix'].set(1)


		return

	def fit_param_option(self,frame,flag):

		self.l1 = tk.Label(frame, text="================== Fitting Param Option ==================",justify='center')
		self.l1.grid(row=0, column=0,columnspan=9,pady=2)

		b1 = tk.Button(frame, text="RESET",height = 1,width = 5,command=lambda: self.reset_param_option(flag))
		b1.grid(row=0, column=7,columnspan=2)		

		count1 = 1; count2 = 1;
		titles = ['val','min','max','fix']
		index  = [1,3,5,7]
		for i in range(4):
			self.l1 = tk.Label(frame,text=titles[i].upper(),justify='center')
			self.l1.grid(row=1, column=index[i],columnspan=2,pady=2)			

		for i in range(10):
			self.l1 = tk.Label(frame,text='a%i'%(i+1),justify='center')
			self.l1.grid(row=i+2, column=0,pady=2)

			self.__dict__['note_%s'%flag]['e%i'%count1] = tk.Entry(frame,width=7,justify='center')
			self.__dict__['note_%s'%flag]['e%i'%count1].grid(row=i+2, column=1,columnspan=2)
			self.__dict__['note_%s'%flag]['e%i'%count1].insert(10,self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['val'].get())
			self.__dict__['note_%s'%flag]['e%i'%count1].configure(state='readonly')

			self.__dict__['note_%s'%flag]['e%i'%(count1+1)] = tk.Entry(frame,width=7,justify='center')
			self.__dict__['note_%s'%flag]['e%i'%(count1+1)].grid(row=i+2, column=3,columnspan=2)
			self.__dict__['note_%s'%flag]['e%i'%(count1+1)].insert(10,self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['min'].get())
			self.__dict__['note_%s'%flag]['e%i'%(count1+1)].configure(state='readonly')

			self.__dict__['note_%s'%flag]['e%i'%(count1+2)] = tk.Entry(frame,width=7,justify='center')
			self.__dict__['note_%s'%flag]['e%i'%(count1+2)].grid(row=i+2, column=5,columnspan=2)
			self.__dict__['note_%s'%flag]['e%i'%(count1+2)].insert(10,self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['max'].get())
			self.__dict__['note_%s'%flag]['e%i'%(count1+2)].configure(state='readonly')

			self.__dict__['note_%s'%flag]['c%i'%count2] = tk.Checkbutton(frame,variable=self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['fix'])
			self.__dict__['note_%s'%flag]['c%i'%count2].grid(row=i+2, column=7,columnspan=2)
			self.__dict__['note_%s'%flag]['c%i'%count2].configure(state='disabled')

			count1 = count1 + 3
			count2 = count2 + 1

		return

	def reset_param_option(self,flag):
		func = self.note_fn['func_type'][flag].get().lower()
		self.fit.lmfit_init_param(self.fit.func['num'][func],flag )
		self.fit.post['func_history'][flag][self.fit.func['num'][func]] = True
		self.post_opt['param'][flag][func]['is'] = True
		TEMP = copy.deepcopy(self.fit.param[flag])

		count1 = 1; count2 = 1;
		for i in range(10):
			self.__dict__['note_%s'%flag]['e%i'%count1].configure(state='normal')
			self.__dict__['note_%s'%flag]['e%i'%(count1+1)].configure(state='normal')
			self.__dict__['note_%s'%flag]['e%i'%(count1+2)].configure(state='normal')
			self.__dict__['note_%s'%flag]['c%i'%count2].configure(state='normal')	

			self.__dict__['note_%s'%flag]['e%i'%count1].delete(0,'end')
			self.__dict__['note_%s'%flag]['e%i'%(count1+1)].delete(0,'end')
			self.__dict__['note_%s'%flag]['e%i'%(count1+2)].delete(0,'end')

			self.__dict__['note_%s'%flag]['e%i'%count1].insert(10,'%5.3f'%TEMP['val'][i])
			self.__dict__['note_%s'%flag]['e%i'%(count1+1)].insert(10,'%5.3f'%TEMP['min'][i])
			self.__dict__['note_%s'%flag]['e%i'%(count1+2)].insert(10,'%5.3f'%TEMP['max'][i])

			count1 = count1 + 3
			count2 = count2 + 1

			if TEMP['vary'][i]:
				self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['fix'].set(0)
			else:
				self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['fix'].set(1)

		self.fitting_func_option_a(flag=flag,isupdate=False)
		return

	def fit_constraint_option(self,frame,start_row,flag):

		self.l1 = tk.Label(frame, text="================== Fitting Constraint ==================",justify='center')
		self.l1.grid(row=start_row, column=0,columnspan=9,pady=2)
		count1 = start_row; count2 = 1;

		for k in self.fit.__dict__['%s_list'%flag]:
			count1 = count1 + 1
			self.l1 = tk.Label(frame, text='%s'%k.upper(),justify='center')
			self.l1.grid(row=count1+1, column=0,pady=2)
		count1 = start_row;

		titles = ['AVG','OLI','STD','ESTD']
		self.l1 = tk.Label(frame, text='-CH-',justify='center')
		self.l1.grid(row=count1+1, column=0,pady=2)		
		len2 = len(self.fit.__dict__['%s_list'%flag])
		len3 = len(titles)
		for k in range(len2):
			for j in range(len3):
				kk = self.fit.__dict__['%s_list'%flag][k]
				if k== 0:
					self.l1 = tk.Label(frame, text='%s'%titles[j].upper(),justify='center')
					self.l1.grid(row=count1+1, column=2*j+1,columnspan=2,pady=2)

				self.__dict__['note_%s'%flag]['c%i'%(len3*k+j+11)] = tk.Checkbutton(frame,variable=self.__dict__['note_%s'%flag][titles[j].lower()][kk],command=lambda: self.fit_constraint_option_a())
				self.__dict__['note_%s'%flag]['c%i'%(len3*k+j+11)].grid(row=count1+2+k, column=2*j+1,columnspan=2)			
		return		count1+2+k

	def fit_constraint_option_a(self,init=False):

		titles = ['AVG','OLI','STD','ESTD']
		len3 = len(titles)
		for flag in self.fit.prof_list:
			len2 = len(self.fit.__dict__['%s_list'%flag])
			for k in range(len2):
				kk = self.fit.__dict__['%s_list'%flag][k]
				if self.__dict__['note_%s'%flag]['avg'][kk].get()==1:
					if not init: self.post_opt['oli'][flag][kk] = self.__dict__['note_%s'%flag]['oli'][kk].get()
					self.__dict__['note_%s'%flag]['oli'][kk].set(0)
					self.__dict__['note_%s'%flag]['c%i'%(len3*k+1+11)].configure(state='disabled')
				else:
					if not self.post_opt['oli'][flag][kk] == None: 
						self.__dict__['note_%s'%flag]['oli'][kk].set(self.post_opt['oli'][flag][kk])
						self.post_opt['oli'][flag][kk] = None
					self.__dict__['note_%s'%flag]['c%i'%(len3*k+1+11)].configure(state='normal')

				if self.__dict__['note_%s'%flag]['std'][kk].get()==0:
					self.post_opt['estd'][flag][kk] = self.__dict__['note_%s'%flag]['estd'][kk].get()
					self.__dict__['note_%s'%flag]['estd'][kk].set(0)
					self.__dict__['note_%s'%flag]['c%i'%(len3*k+3+11)].configure(state='disabled')
				else:
					if not self.post_opt['estd'][flag][kk] == None: 
						self.__dict__['note_%s'%flag]['estd'][kk].set(self.post_opt['estd'][flag][kk])
						self.post_opt['estd'][flag][kk] = None
					self.__dict__['note_%s'%flag]['c%i'%(len3*k+3+11)].configure(state='normal')
		return

	def fit_fuction_plot(self,start_row):
		self.l1 = tk.Label(self.home['page2'], text="================== Fitting Functions ==================",justify='center')
		self.l1.grid(row=start_row, column=0,columnspan=9,pady=5)
		self.l1 = tk.Label(self.home['page2'], text='FUNC_TYPE')
		self.l1.grid(row=start_row+1,column=0,pady=5,padx=5)
		self.m5 = tk.OptionMenu(self.home['page2'],self.note_fn['func_type']['etc'],*self.func_list,command=lambda value: self.fit_function_plot_a()) 
		self.m5.config(width=6)
		self.m5.grid(row=start_row+1,column=1,columnspan=3,pady=1,sticky='w')
		return

	def fit_function_plot_a(self):
	
		self.figure['name']['plot'].canvas.draw_idle()
		[ax] = self.figure['name']['plot'].axes
		ax.cla()
		ax.set_facecolor('g')
		ax.axis('off')
		ax.set_xlim(0,1.)
		ax.set_ylim(0,2.)
		if self.note_fn['func_type']['etc'].get().lower() == 'core':
			line = 'y = $a_1 +  a_2(1 - x^{a_3})^{a_4}$'
			line = line + '\n\n\n'
		elif self.note_fn['func_type']['etc'].get().lower() == 'mtanh':
			line = '$mtanh(x,b) = \\frac{(1.0 + bx)e^x - e^{-x}}{e^x + e^{-x}}$\n'
			line = line + '$F(x) = \\frac{a_2 - a_1}{2}{\\times}(mtanh[\\frac{a_4-x}{a_3/2},a_5] + 1) + a_1$\n'
			line = line + '$y = F(x) + (a_6-F(x)){\\times}Exp[-(\\frac{x}{a_7+0.2})^{a_8}]$'
		elif self.note_fn['func_type']['etc'].get().lower() == 'ptanh':
			line = '$y = \\frac{a_2 -  a_1}{2}{\\times}(1 + a_3x + a_4x^2 + a_5x^3)$\n'
			line = line +'$y = y{\\times}(1 - tanh[\\frac{x-a_7}{a_6/2}]) + a_1$'
			line = line + '\n'
		elif self.note_fn['func_type']['etc'].get().lower() == 'eped':
			line = '$y = a_1 + a_2(tanh[1] - tanh[\\frac{x - 1 + a_3/2}{a_3/2}])$\n'
			line = line + '$y = y + a_4(1 - (\\frac{x}{1-a3})^{a_5})^{a_6}$'
			line = line + '\n'
		elif self.note_fn['func_type']['etc'].get().lower() == 'eped2':
			line = '$y = a_1 + a_2(tanh[\\frac{1-a7}{a_3/2}] - tanh[\\frac{x - a7}{a3/2}] $\n'
			line = line + '$y = y + a_4(1 - (\\frac{x}{a7-a_3/2})^{a_5})^{a_6}$'
			line = line + '\n'
		elif self.note_fn['func_type']['etc'].get().lower() == 'spline':
			line = 'Smoothed Cubic SPLINE\n\n\n'
		else:   
			line = 'Cubic SPLINE with \n'
			line = line + 'Optimized knots \n\n'

		ax.text(-0.12,0.,line,color='b')

		return

	def fit_variable_option(self,frame,start_row,flag):

		self.l1 = tk.Label(frame, text="================== Fitting Variables ==================",justify='center')
		self.l1.grid(row=start_row, column=0,columnspan=9,pady=5)
		count1 = start_row; count2 = 1;

		b1 = tk.Button(frame, text="SELEC",height = 1,width = 5,command=lambda: self.fit_exclude(flag))
		b1.grid(row=start_row, column=7,columnspan=2)		

		for k in self.fit.__dict__['%s_list'%flag]:
			count1 = count1 + 1
			self.l1 = tk.Label(frame, text='%s'%k.upper(),justify='center')
			self.l1.grid(row=count1+1, column=0,pady=3)
		count1 = start_row;

		titles   = ['RANGE','OLCN','OLC','WEIGHT']
		var_name = ['RANGE','OLCN','OLC','WEIGHT']
		self.l1 = tk.Label(frame, text='-CH-',justify='center')
		self.l1.grid(row=count1+1, column=0,pady=3)		
		len2 = len(self.fit.__dict__['%s_list'%flag])
		len3 = len(titles)
		for k in range(len2):
			for j in range(len3):
				kk = self.fit.__dict__['%s_list'%flag][k]
				if k== 0:
					self.l1 = tk.Label(frame, text='%s'%titles[j],justify='center')
					self.l1.grid(row=count1+1, column=2*j+1,columnspan=2,pady=3)

				self.__dict__['note_%s'%flag]['e%i'%(len3*k+j+31)] = tk.Entry(frame,width=7,justify='center')
				if (k<=1 or not j==1):
					self.__dict__['note_%s'%flag]['e%i'%(len3*k+j+31)].grid(row=count1+2+k, column=2*j+1,columnspan=2)	
					self.__dict__['note_%s'%flag]['e%i'%(len3*k+j+31)].insert(10,self.__dict__['note_%s'%flag][var_name[j].lower()][kk].get())		
		return		count1+2+k

	def fit_inter(self,frame,start_row,flag):

		self.l1 = tk.Label(frame, text="================== Interf Constraint ==================",justify='center')
		self.l1.grid(row=start_row, column=0,columnspan=9,pady=5)
		self.l1 = tk.Label(frame,text='-CH-')
		self.l1.grid(row=start_row+1,column=0,pady=2)


		titles = ['SIG','EXP','DIFF','SCALE']

		len2 = len(self.fit.inter_list)
		for i in range(4):
			self.l1 = tk.Label(frame,text='%s'%titles[i].upper())
			self.l1.grid(row=start_row+1,column=2*i+1,columnspan=2,pady=2)

		for k in range(len2):
			chan = self.fit.inter_list[k]
			self.l1 = tk.Label(frame,text='%s'%chan.upper())
			self.l1.grid(row=start_row+k+2,column=0,pady=2)

			self.note_ne['e%i'%(k+100)] = tk.Entry(frame,width=7,justify='center')
			self.note_ne['e%i'%(k+100)].grid(row=start_row+k+2, column=1,columnspan=2)
			self.note_ne['e%i'%(k+100)].insert(10,self.note_in[chan.lower()+'s'].get())

			self.note_ne['e%i'%(k+107)] = tk.Entry(frame,width=7,justify='center')
			self.note_ne['e%i'%(k+107)].grid(row=start_row+k+2, column=3,columnspan=2)
			self.note_ne['e%i'%(k+107)].insert(10,self.note_in[chan.lower()].get())
			self.note_ne['e%i'%(k+107)].configure(state='readonly')

			self.note_ne['e%i'%(k+114)] = tk.Entry(frame,width=7,justify='center')
			self.note_ne['e%i'%(k+114)].grid(row=start_row+k+2, column=5,columnspan=2)
			if self.note_in[chan.lower()] == 0.:
				self.note_ne['e%i'%(k+114)].insert(10,'-N/A-')
			else:
				self.note_ne['e%i'%(k+114)].insert(10,self.note_in[chan.lower()+'d'].get())
			self.note_ne['e%i'%(k+114)].configure(state='readonly')

		self.note_ne['c101'] = tk.Checkbutton(frame,variable=self.note_in['int01'.lower()+'f'], command=lambda: self.fit_inter_a('int01'))
		self.note_ne['c102'] = tk.Checkbutton(frame,variable=self.note_in['int02'.lower()+'f'], command=lambda: self.fit_inter_a('int02'))
		self.note_ne['c103'] = tk.Checkbutton(frame,variable=self.note_in['tci01'.lower()+'f'], command=lambda: self.fit_inter_a('tci01'))
		self.note_ne['c104'] = tk.Checkbutton(frame,variable=self.note_in['tci02'.lower()+'f'], command=lambda: self.fit_inter_a('tci02'))
		self.note_ne['c105'] = tk.Checkbutton(frame,variable=self.note_in['tci03'.lower()+'f'], command=lambda: self.fit_inter_a('tci03'))
		self.note_ne['c106'] = tk.Checkbutton(frame,variable=self.note_in['tci04'.lower()+'f'], command=lambda: self.fit_inter_a('tci04'))
		self.note_ne['c107'] = tk.Checkbutton(frame,variable=self.note_in['tci05'.lower()+'f'], command=lambda: self.fit_inter_a('tci05'))
		for i in range(7): self.note_ne['c%i'%(i+101)].grid(row=start_row+2+i, column=7,columnspan=2)

		return

	def fit_inter_a(self,chan):

		sum0 = 0
		for i in self.fit.inter_list:
			sum0 = sum0 + self.note_in[i.lower()+'f'].get()
			if (self.note_in[i.lower()+'f'].get() > 0 and not chan == i): flag = i

		if sum0 > 1: 
			self.note_in[flag.lower()+'f'].set(0)

		return

	def fit_exclude(self,flag):

		if not self.t_close1:
			print('>>> Window is already opened...')
			return		

		if not 	self.fit.post['isrfile'][flag]:
			print('>>> No raw data files...')
			return

		frame = tk.Toplevel(self.root)
		frame.wm_title("Profile Ch. Select")
		self.t_close1 = False
		frame.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(1,frame))
		frame.grab_set()
		l = tk.Label(frame, text="Channel Selection")
		len2 = len(self.fit.__dict__['%s_list'%flag])

		n_channel = 0
		count1 = -1;
		for k in range(len2):
			kk = self.fit.__dict__['%s_list'%flag][k]
			if self.fit.fit_opt['file'][flag][kk]==None: continue
			self.exclude[kk] = dict()

			count1 = count1 + 1

			self.l1 = tk.Label(frame, text="= Ch. %s-%s ="%(flag.upper(),kk.upper()),justify='center')
			self.l1.grid(row=0, column=2*count1+1,columnspan=1,pady=3)
			self.l1 = tk.Label(frame, text="= # =",justify='center')
			self.l1.grid(row=0, column=2*count1,columnspan=1,pady=3)			

			len3 = self.fit.post['prof_dim'][flag][kk][0]
			for i in range(len3):
				self.exclude[kk][i] = tk.IntVar()
				self.exclude[kk][i].set(1)

			excn    = self.fit.fit_opt['exclude'][flag][kk]
			if not excn == '':
				if excn[-1] == ',':	excn = excn[:-1].split(',')
				else:	excn = excn.split(',')
			else: excn = []
			excn = np.array(excn,dtype=int)
			len4 = len(excn)

			for i in range(len4):
				self.exclude[kk][excn[i]].set(0)

			count2 = 0
			for i in range(len3):
				count2 = count2 + 1
				self.l1 = tk.Label(frame, text="#%2i"%(i+1),justify='center')
				self.l1.grid(row=count2, column=2*count1,pady=3)				
				self.__dict__['exc%i'%(i+1+n_channel)] = tk.Checkbutton(frame,variable=self.exclude[kk][i])
				self.__dict__['exc%i'%(i+1+n_channel)].grid(row=count2, column=2*count1+1,pady=3)

				if  count2==10:
					count1 = count1 + 1
					count2 = 0
					self.l1 = tk.Label(frame, text="= Ch. %s-%s ="%(flag.upper(),kk.upper()),justify='center')
					self.l1.grid(row=0, column=2*count1+1,columnspan=1,pady=3)
					self.l1 = tk.Label(frame, text="= # =",justify='center')
					self.l1.grid(row=0, column=2*count1,columnspan=1,pady=3)					

			n_channel = n_channel + len3

		self.exclude_flag = flag
		return

	def make_exclude(self):

		flag = self.exclude_flag
		for kk in self.fit.__dict__['%s_list'%flag]:
			excn = []
			if self.fit.fit_opt['file'][flag][kk]==None: continue
			len2 = self.fit.post['prof_dim'][flag][kk][0]
			for i in range(len2):
				if self.exclude[kk][i].get() == 0: excn = np.append(excn,i)

			excn = np.array(excn,dtype='int')
			len3 = len(excn)
			if len3 == 0: self.fit.fit_opt['exclude'][flag][kk] = ''
			else:
				excn2 = '%i'%excn[0]
				for i in range(len3-1): excn2 = excn2 + ',%i'%excn[i+1]
				self.fit.fit_opt['exclude'][flag][kk] = excn2
			if len3 > 0:
				print('>>> Ch. %s of %s - %s are excluded'%(self.fit.fit_opt['exclude'][flag][kk],flag.upper(),kk.upper()))
			else:	print('>>> All Ch.of %s - %s are included'%(flag.upper(),kk.upper()))

		return

	def make_mds(self):

		try: os.mkdir('GFILES')
		except: pass
		try: os.mkdir('INPUTS')
		except: pass
		try: os.mkdir('MFIT')
		except: pass		
		try: os.mkdir(self.mds_dir)
		except: pass

		self.l1 = tk.Label(self.home['page7'], text="===================== LOAD GFILE =====================",justify='center')
		self.l1.grid(row=0, column=0,columnspan=9,pady=5)
		b1 = tk.Button(self.home['page7'],text="LOAD", height = 1, width = 5,command=lambda: self.make_mds_a())	
		b1.grid(row=1, column=8,columnspan=1,pady=5)

		self.l1 = tk.Label(self.home['page7'], text='SHOT [#]',justify='center')
		self.l1.grid(row=1, column=0,pady=3)		
		self.note_md['e1'] = tk.Entry(self.home['page7'],width=9,justify='center')
		self.note_md['e1'].insert(10,self.note_md['shot'].get())
		self.note_md['e1'].grid(row=1,column=1,columnspan=2)		

		self.l1 = tk.Label(self.home['page7'], text='TIME [ms]',justify='center')
		self.l1.grid(row=1, column=3,pady=3)		
		self.note_md['e2'] = tk.Entry(self.home['page7'],width=9,justify='center')	
		self.note_md['e2'].insert(10,self.note_md['time'].get())			
		self.note_md['e2'].grid(row=1,column=4,columnspan=2)			

#		self.l1 = tk.Label(self.home['page7'], text='EFIT[#]',justify='center')
#		self.l1.grid(row=1, column=6,pady=3)
		self.note_md['m1'] = tk.OptionMenu(self.home['page7'],self.note_md['efit_type'],'1','2','3','4','5')
		self.note_md['m1'].grid(row=1, column=6)
		self.note_md['m1'].config(width=1)

		self.l1 = tk.Label(self.home['page7'], text="==================== LOAD TS-CES ====================",justify='center')
		self.l1.grid(row=2, column=0,columnspan=9,pady=5)
		b2 = tk.Button(self.home['page7'],text="LOAD TS & CES", height = 2, width = 30,command=lambda: self.make_mds_b())	
		b2.grid(row=3, column=3,columnspan=6,pady=5,rowspan=2)
		self.l1 = tk.Label(self.home['page7'], text='TS-dt [ms]',justify='center')
		self.l1.grid(row=3, column=0,pady=3)		
		self.l1 = tk.Label(self.home['page7'], text='CES-dt [ms]',justify='center')
		self.l1.grid(row=4, column=0,pady=3)			
		self.note_md['e3'] = tk.Entry(self.home['page7'],width=5,justify='center')				
		self.note_md['e3'].grid(row=3,column=1,columnspan=2)	
		self.note_md['e3'].insert(10,self.note_md['dt']['ts'].get())
		self.note_md['e4'] = tk.Entry(self.home['page7'],width=5,justify='center')				
		self.note_md['e4'].grid(row=4,column=1,columnspan=2)			
		self.note_md['e4'].insert(10,self.note_md['dt']['ces'].get())			

		self.l1 = tk.Label(self.home['page7'], text="====================== LOAD TCI ======================",justify='center')
		self.l1.grid(row=5, column=0,columnspan=9,pady=5)
		b3 = tk.Button(self.home['page7'],text="LOAD TCI", height = 1, width = 30,command=lambda: self.make_mds_c())	
		b3.grid(row=6, column=3,columnspan=6,pady=5)
		self.l1 = tk.Label(self.home['page7'], text='dt [ms]',justify='center')
		self.l1.grid(row=6, column=0,pady=3)		
		self.note_md['e5'] = tk.Entry(self.home['page7'],width=5,justify='center')				
		self.note_md['e5'].grid(row=6,column=1,columnspan=2)			
		self.note_md['e5'].insert(10,self.note_md['dt']['tci'].get())

		self.l1 = tk.Label(self.home['page7'], text="============== LOAD ECE (by Y.H.Lee, ~1min) ==============",justify='center')
		self.l1.grid(row=7, column=0,columnspan=9,pady=5)
		b4 = tk.Button(self.home['page7'],text="LOAD ECE", height = 1, width = 19,command=lambda: self.make_mds_d())	
		b4.grid(row=8, column=4,columnspan=5,padx=12,sticky='e')
		self.l1 = tk.Label(self.home['page7'], text='dt [ms]',justify='center')
		self.l1.grid(row=8, column=0,pady=3)
		self.l1 = tk.Label(self.home['page7'], text='HFS',justify='center')
		self.l1.grid(row=8, column=3)

		self.note_md['c2'] = tk.Checkbutton(self.home['page7'],variable=self.note_md['ece_hfs'])
		self.note_md['c2'].grid(row=8, column=4,sticky='w')
		self.note_md['c2'].configure(state='disabled')
	
		self.note_md['e6'] = tk.Entry(self.home['page7'],width=5,justify='center')				
		self.note_md['e6'].grid(row=8,column=1,columnspan=2)	
		self.note_md['e6'].insert(10,self.note_md['dt']['ece'].get())		

		self.l1 = tk.Label(self.home['page7'], text="============ LOAD REFLEC (SINGLE FIT ONLY) ============",justify='center')
		self.l1.grid(row=9, column=0,columnspan=9,pady=5)
		b5 = tk.Button(self.home['page7'],text="LOAD REFLECTOMETER", height = 2, width = 30,command=lambda: self.make_mds_e())	
		b5.grid(row=10, column=3,columnspan=6,rowspan=2,pady=5)		
		self.l1 = tk.Label(self.home['page7'], text='dt [ms]',justify='center')
		self.l1.grid(row=11, column=0,pady=3)		
		self.l1 = tk.Label(self.home['page7'], text='time [ms]',justify='center')
		self.l1.grid(row=10, column=0,pady=3)
		self.note_md['e7'] = tk.Entry(self.home['page7'],width=5,justify='center')				
		self.note_md['e7'].grid(row=11,column=1,columnspan=2)			
		self.note_md['e7'].insert(10,self.note_md['dt']['ref'].get())
		self.note_md['e9'] = tk.Entry(self.home['page7'],width=5,justify='center')
		self.note_md['e9'].grid(row=10,column=1,columnspan=2)
		self.note_md['e9'].insert(10,self.note_md['dt']['ref2'].get())
		self.l1 = tk.Label(self.home['page7'], text="================== MDS LOAD STATUS ==================",justify='center')
		self.l1.grid(row=12, column=0,columnspan=9,pady=5)

		self.make_mds_a0()

		self.l1 = tk.Label(self.home['page7'], text="=================== MULTI-TIME OPT. ===================",justify='center')
		self.l1.grid(row=16, column=0,columnspan=9,pady=5)
		self.note_md['e8'] = tk.Entry(self.home['page7'],width=40,justify='center')				
		self.note_md['e8'].grid(row=17,column=1,columnspan=8)	
		self.note_md['e8'].insert(10,self.note_md['times'].get())
		self.l1 = tk.Label(self.home['page7'], text='TIMEs [ms]',justify='center')
		self.l1.grid(row=17, column=0,pady=3)			

		self.note_md['e10'] = tk.Entry(self.home['page7'],width=7,justify='center')
		self.note_md['e10'].grid(row=18,column=1,columnspan=2)
		self.note_md['e10'].insert(10,self.note_md['dacrit'].get())
		self.l1 = tk.Label(self.home['page7'], text='DACRIT[#]',justify='center')
		self.l1.grid(row=18, column=0,pady=3)

		self.note_md['e11'] = tk.Entry(self.home['page7'],width=7,justify='center')
		self.note_md['e11'].grid(row=18,column=4,columnspan=2)
		self.note_md['e11'].insert(10,self.note_md['duty'].get())
		self.l1 = tk.Label(self.home['page7'], text='DUTY[%] ',justify='center')
		self.l1.grid(row=18, column=3,pady=3)
	
		b10 = tk.Button(self.home['page7'],text="ELM CHECK", height = 1, width = 10, command=lambda: self.make_mds_i())
		b10.grid(row=18, column=6,columnspan=4,pady=5,padx=10,sticky='e')

		b6 = tk.Button(self.home['page7'],text="CHANNEL CHECK", height = 1, width = 23,command=lambda: self.make_mds_g())
		b6.grid(row=19, column=3,columnspan=9,pady=5,padx=10,sticky='e')
		b7 = tk.Button(self.home['page7'],text="2D-RT VIEW", height = 1, width = 23,command=lambda: self.make_mds_h())
		b7.grid(row=19, column=0,columnspan=9,pady=5,padx=10,sticky='w')

		b8 = tk.Button(self.home['page7'],text="RESTORE SINGLE FIT", height = 1, width = 23,command=lambda: self.restore_singlefit())
		b8.grid(row=21, column=0,columnspan=9,pady=5,padx=10,sticky='w')
		b9 = tk.Button(self.home['page7'],text="RESTORE MULTI FIT", height = 1, width = 23, command=lambda: self.restore_multifit())
		b9.grid(row=21, column=0,columnspan=9,pady=5,padx=10,sticky='e')

		self.note_md['b10'] = tk.Button(self.home['page7'],text="NE CROSS CHECK", height = 1, width = 23, command=lambda: self.make_mds_j())
		self.note_md['b10'].grid(row=20, column=0,columnspan=9,pady=5,padx=10,sticky='w')

		if (self.single_only or self.multi_only):
			self.note_md['e8'].configure(state='disabled')
			b1.configure(state='disabled')
			b2.configure(state='disabled')
			b6.configure(state='disabled')
			b6.configure(state='disabled')
			b7.configure(state='disabled')
			b8.configure(state='disabled')
			b9.configure(state='disabled')

		return

	def make_mds_a0(self):

		for k in range(4): self.note_md['s%i'%(k+1)] = tk.Entry(self.home['page7'],width=12,justify='center')
		for k in range(5,9): self.note_md['s%i'%(k)] = tk.Entry(self.home['page7'],width=26,justify='center')
		
		self.note_md['s1'].grid(row=13,column=0,columnspan=9,pady=5,padx=10,sticky='w')
		self.note_md['s2'].grid(row=13,column=0,columnspan=9,pady=5,padx=107,sticky='w')
		self.note_md['s3'].grid(row=13,column=0,columnspan=9,pady=5,padx=107,sticky='e')
		self.note_md['s4'].grid(row=13,column=0,columnspan=9,pady=5,padx=10,sticky='e')
		self.note_md['s5'].grid(row=14,column=0,columnspan=9,pady=5,padx=10,sticky='w')
		self.note_md['s6'].grid(row=14,column=0,columnspan=9,pady=5,padx=10,sticky='e')
		self.note_md['s7'].grid(row=15,column=0,columnspan=9,pady=5,padx=10,sticky='w')
		self.note_md['s8'].grid(row=15,column=0,columnspan=9,pady=5,padx=10,sticky='e')

		for k in range(8): self.note_md['s%i'%(k+1)].configure(state='readonly')
		return
	
	def check_mds_in(self):
		isok = True
		try:
			shot = int(self.note_md['e1'].get())
			time = int(self.note_md['e2'].get())
		except: print('>>> No Shot # and time given'); isok = False
		if not shot>0: isok = False
		if not time>0: isok = False
		return isok

	def make_mds_a(self):
		if not self.check_mds_in(): return
		shot = int(self.note_md['e1'].get())
		time = int(self.note_md['e2'].get())
		efit_no = int(self.note_md['efit_type'].get())
		self.note_md['shot'].set(shot)
		self.note_md['time'].set(time)
		efit_list = get_efit_list2(shot)
		if not efit_list['isefit'][efit_no]: print('>>> No EFIT'); return
		tind = np.argmin(abs(efit_list['times'][efit_no]-time))
		ttime= efit_list['times'][efit_no][tind]
		filename = 'GFILES/g%06i.%06i_%i'%(shot,time,efit_no)
		if not os.path.isfile(filename):
				filename1 = efit_list['dirs'][efit_no]+'g%06i.%06i'%(shot,ttime)
				copyfile(filename1,filename)
		copyfile(filename,self.mds_dir+'/'+filename.split('/')[1])
		self.note_in['file']['gfile'].set(filename)		
		self.sync_gui_opt(g2v=False,skip_mds=True)
		self.make_mds_e2()
		return

	def make_mds_b(self):
		if not self.check_mds_in(): return
		shot = int(self.note_md['e1'].get())
		time = int(self.note_md['e2'].get())
		try:
			dt1  = int(float(self.note_md['e3'].get()))
			dt2  = int(float(self.note_md['e4'].get()))
		except: print('>>> No averaging time is given'); return
		
		if (dt1==0 and dt2==0): return
		#dt   = max(dt1,dt2)
		self.note_md['shot'].set(shot)
		self.note_md['time'].set(time)
		self.note_md['dt']['ts'].set(dt1)
		self.note_md['dt']['ces'].set(dt2)
		cdir = os.getcwd()
		os.chdir(self.mds_dir)

		command = '%s %i %i %i %i'%(mds_lit,shot,time,max(dt2,10),max(dt1,10))
		os.system(command)
		os.chdir(cdir)

		if os.path.isfile(self.mds_dir+'/result.dat'):
			f = open(self.mds_dir+'/result.dat','r')
			while True:
				line = f.readline()
				if not line: break
				for i in range(4):
					flag = self.fit.prof_list[i]
					if line.split()[0] == '%s_file'%flag:
						if len(line.split()) > 1:
							line2 = self.mds_dir+'/'+line.split()[1].split('/')[-1]
							if (flag=='ne' or flag=='te'):							
								try:self.note_in['file'][flag]['ts'].set(line2)
								except: print('>>> No %s input file'%flag); pass
							else:							
								try:self.note_in['file'][flag]['ces'].set(line2)
								except: print('>>> No %s input file'%flag); pass
					if line.split()[0] == '%s_edge_file'%flag:
						if len(line.split()) > 1:
							line2 = self.mds_dir+'/'+line.split()[1].split('/')[-1]
							if (flag=='ne' or flag=='te'):
								try:self.note_in['file'][flag]['tse'].set(line2)
								except: print('>>> No %s input file'%flag); pass

			f.close()
		self.sync_gui_opt(g2v=False,skip_mds=True)
		self.make_mds_e2()
		return

	def make_mds_c(self):
		if not self.check_mds_in(): return
		shot = int(self.note_md['e1'].get())
		time = int(self.note_md['e2'].get())
		try: dt  = int(self.note_md['e5'].get())
		except: print('>>> No averaging time is given'); return
		if (dt == 0): return
		self.note_md['shot'].set(shot)
		self.note_md['time'].set(time)
		self.note_md['dt']['tci'].set(dt)

		cdir = os.getcwd()
		os.chdir(self.mds_dir)
		command = '%s %s %i %i %i %s'%(pythonc_exec, mds_tci,shot,time,dt,'tci_result.dat')
		os.system(command)
		os.chdir(cdir)

		if os.path.isfile(self.mds_dir+'/tci_result.dat'):
			f = open(self.mds_dir+'/tci_result.dat','r')
			for i in range(7):
				line = f.readline()	
				flag = self.fit.inter_list[i]
				if float(line.split()[0]) == 0: self.note_in[flag].set(line.split()[1])
				else: 
					self.note_in[flag].set(line.split()[1])
					if len(line.split()) > 2: self.note_in[flag.lower()+'s'].set(line.split()[2])
	
			f.close()
		self.sync_gui_opt(g2v=False,skip_mds=True)
		self.make_mds_e2()
		return

	def make_mds_d(self):
		if not self.check_mds_in(): return
		shot = int(self.note_md['e1'].get())
		time = int(self.note_md['e2'].get())

		try: dt   = int(self.note_md['e6'].get())
		except: print('>>> No averaging time is given'); return

		if (dt == 0): return
		self.note_md['shot'].set(shot)
		self.note_md['time'].set(time)
		self.note_md['dt']['ece'].set(dt)
		cdir = os.getcwd()
		os.chdir(self.mds_dir)
		if not os.path.isfile('DATASAVE/%i/ECE_size'%shot):
			if self.note_md['ece_hfs'].get() == 0: gefit_ece(shot)
			else:  gefit_ece(shot, lfs_option ='n')
		#load_ece(shot,[time/1.e3],dt/1.e3,'',None)
		self.make_mds_d1(shot,time,dt)
		os.chdir(cdir)
		filename = self.mds_dir+'/te_%06i_%ims_ECE_2nd.dat'%(shot,time)
		if os.path.isfile(filename):
			self.note_in['e4'].delete(0,'end')
			self.note_in['e4'].insert(10,filename)
		self.make_mds_e2()
		return

	def make_mds_e(self):
		if not self.check_mds_in(): return
		shot = int(self.note_md['e1'].get())

		try: time = int(self.note_md['e9'].get())
		except: print('>>> No time slice is given'); time = int(self.note_md['e2'].get())

		try: dt  = int(self.note_md['e7'].get())
		except: print('>>> No averaging time is given'); return
		if (dt == 0 or shot == 0): return
		self.note_md['dt']['ref'].set(dt)
		self.note_md['dt']['ref2'].set(time)
		cdir = os.getcwd()
		os.chdir(self.mds_dir)
		command = '%s %s %i %f %f %s'%(pythonc_exec, mds_ref,shot,time/1.e3,dt/1.e3,'ne_%06i_%ims_reflec.dat'%(shot,time))
		os.system(command)
		os.chdir(cdir)
		if os.path.isfile(self.mds_dir+'/ne_%06i_%ims_reflec.dat'%(shot,time)):
			self.note_in['e7'].delete(0,'end')
			self.note_in['e7'].insert(10,'%s/ne_%06i_%ims_reflec.dat'%(self.mds_dir,shot,time))
		self.make_mds_e2()
		return						

	def make_mds_d1(self,shot,time,dt):

		if not self.t_close2: return

		frame = tk.Toplevel(self.root)
		frame.wm_title("ECE Time")
		self.t_close2 = False
		frame.protocol("WM_DELETE_WINDOW", lambda: self.detect_close(2,frame))
#		frame.resizable(0,0)
	
		fig, ax = plt.subplots(1,1,figsize=(8,5))
		canvas = FigureCanvasTkAgg(fig,master=frame)
		plot_widget = canvas.get_tk_widget()
		plot_widget.grid(rowspan=2,row=1,column=0)

		toolbar_frame = tk.Frame(frame)
		toolbar_frame.grid(column=0,row=0)
		
		toolbar =  NavigationToolbar2Tk(canvas,toolbar_frame)
		fig.canvas.draw_idle()

		load_ece(shot,[time/1.e3],dt/1.e3,'',fig)
		frame.grab_set()

		return

	def make_mds_e2(self):
		
		try: shot = int(float(self.note_md['e1'].get()))
		except: shot = 0
		try: time = int(float(self.note_md['e2'].get()))
		except: time = 0
		try: times = np.array(self.note_md['e8'].get().split(','),dtype='int')
		except: times = [0]

		dirs = '%s/DATASAVE/%i/'%(self.mds_dir,shot)
		efit_no = int(self.note_md['efit_type'].get())
		gdirs = 'GFILES/g%06i.%06i_%i'%(shot,time,efit_no)
		list1 = ['tst','tsn','cest','cesv','tci','ece','ref','g']
		if shot == 0: 
			for k in list1: self.mds_load[k] = False
		else:
		#CHECK TS
			if not os.path.isfile(dirs+'TE_size'): self.mds_load['tst'] = False
			elif os.path.isfile(dirs+'TS_err'): self.mds_load['tst'] = None
			else: self.mds_load['tst'] = True
			if not os.path.isfile(dirs+'NE_size'): self.mds_load['tsn'] = False
			elif os.path.isfile(dirs+'TS_err'): self.mds_load['tsn'] = None
			else: self.mds_load['tsn'] = True
			if not os.path.isfile(dirs+'TI_size'): self.mds_load['cest'] = False
			elif os.path.isfile(dirs+'CES_err'): self.mds_load['cest'] = None
			else: self.mds_load['cest'] = True
			if not os.path.isfile(dirs+'VT_size'): self.mds_load['cesv'] = False
			elif os.path.isfile(dirs+'CES_err'): self.mds_load['cesv'] = None
			else: self.mds_load['cesv'] = True
			if not os.path.isfile(dirs+'TCI_size'): self.mds_load['tci'] = False
			else: self.mds_load['tci'] = True
			if not os.path.isfile(dirs+'ECE_size'): self.mds_load['ece'] = False
			else: self.mds_load['ece'] = True
			if not os.path.isfile(dirs+'REF.npz.gz'): self.mds_load['ref'] = False
			else: self.mds_load['ref'] = True
			if self.runmode == 'single':
				if not os.path.isfile(gdirs): self.mds_load['g'] = False
				else: self.mds_load['g'] = True
			else:
				self.mds_load['g'] = True
				for m in range(len(times)):
					gdirs = 'GFILES/g%06i.%06i_%i'%(shot,times[m],efit_no)
					if not os.path.isfile(gdirs): self.mds_load['g'] = False
	
		list2 = ['TS-TE','TS-NE','CES-TI','CES-VT','TCI','ECE','REFLEC','GFILE']
		for k in range(8):
			self.note_md['s%i'%(k+1)].configure(state='normal')
			self.note_md['s%i'%(k+1)].delete(0,'end')
			flag = list1[k]
			if self.mds_load[flag] == True: 
				self.note_md['s%i'%(k+1)].insert(0,'[%s]-OK'%list2[k])
				self.note_md['s%i'%(k+1)].configure(fg='white',readonlybackground='deepskyblue')
			elif self.mds_load[flag] == False:
				self.note_md['s%i'%(k+1)].insert(0,'[%s]-NO'%list2[k])
				self.note_md['s%i'%(k+1)].configure(fg='white',readonlybackground='tomato')
			else:
				self.note_md['s%i'%(k+1)].insert(0,'[%s]-N/A'%list2[k])
				self.note_md['s%i'%(k+1)].configure(fg='white',readonlybackground='silver')
			
			self.note_md['s%i'%(k+1)].configure(state='readonly')

		if (self.single_only or self.multi_only): self.note_md['e8'].configure(state='disabled')	
		if (self.mds_load['tsn'] and self.mds_load['tci']): self.note_md['b10'].configure(state='normal')
		else: self.note_md['b10'].configure(state='disabled')
	
		return

	def make_mds_f(self):
		if not os.path.isfile('result_single.save'):
			print('>>> You must save the single time fit first!')
			return

		try: shot = int(float(self.note_md['e1'].get()))
		except: print('>>> No Shot #'); return

		self.filename = dict()
		for flag in self.fit.prof_list: self.filename[flag] = dict()
		if len(self.note_md['e8'].get()) == 0: 
			print('>>> Times are not given...')
			return
		self.times = np.array(self.note_md['e8'].get().split(','),dtype='int')
		self.shot  = int(self.note_md['e1'].get())
		efit_list  = get_efit_list2(self.shot)
		efit_no    = int(self.note_md['efit_type'].get())
		tlen = len(self.times)
		self.files = dict()
		try: os.mkdir('GFILES')
		except: pass
		try: os.mkdir('INPUTS')
		except: pass
		self.filename['gfile'] = dict()
		self.filename['tci']   = dict()
		print('>>> Load gfiles')
		if not efit_list['isefit'][efit_no]: print('>>> No EFIT'); return
		for i in range(tlen):
			filename = 'GFILES/g%06i.%06i_%i'%(self.shot,self.times[i],efit_no)
			ttime    = self.times[i]
			tind     = np.argmin(abs(efit_list['times'][efit_no]-ttime))
			ttime2   = efit_list['times'][efit_no][tind]
			if not os.path.isfile(filename):
				filename1 = efit_list['dirs'][efit_no]+'g%06i.%06i'%(self.shot,ttime2)
				copyfile(filename1,filename)
			self.filename['gfile'][i] = filename
			self.filename['tci'][i] = None
			update_progress(float((i+1)/tlen))

		chan = ['TS','TSE','CES','TCI','ECE']
		for flag in self.fit.prof_list:
			self.filename[flag]['kfile'] = dict()
			for kk in self.fit.__dict__['%s_list'%flag]:
				self.filename[flag][kk] = dict()
				for ll in range(tlen):
					self.filename[flag]['kfile'][ll] = None
					self.filename[flag][kk][ll] = None
	
		tas_ind = [3,3,4,5,6]
		dirs  = '%s/DATASAVE/%i'%(self.mds_dir,self.shot)
		times = np.array(self.note_md['e8'].get().split(','),dtype='int');
		self.fit.fit_opt['dacrit'] = float(self.note_md['e10'].get())
		self.fit.fit_opt['duty']   = float(self.note_md['e11'].get())
		self.get_elm_peak(dirs,times);
		for k in range(5):
			flag = chan[k].lower()
			if float(self.note_md['e%i'%(tas_ind[k])].get()) == 0: continue
			print('>>> Load %s raw profiles'%flag.upper())
			if (flag == 'ts' or flag =='tse' or flag == 'ces'): 
				ftemp,flist = self.make_ts_ces_raw_files(float(self.note_md['e%i'%tas_ind[k]].get()),flag,dirs,'INPUTS')
				for kk in flist:
					if ftemp[kk]['isfile']:
						print('>>> There is %s - %s raw file'%(flag.upper(),kk.upper()))
						for ll in range(tlen):
							if (kk=='ne' or kk=='te'): 
								if flag=='ts': self.filename[kk]['ts'][ll] = ftemp[kk][ll]
								else:          self.filename[kk]['tse'][ll]= ftemp[kk][ll]
							else:	self.filename[kk]['ces'][ll] = ftemp[kk][ll]

			if (flag == 'tci'):
				dirs = '%s/DATASAVE/%i'%(self.mds_dir,self.shot)
				ftemp,isfile = self.make_tci_raw_files(float(self.note_md['e%i'%(tas_ind[k])].get()),dirs,'INPUTS')
				if isfile: print('>>> There is %s - %s raw file'%(flag.upper(),'NE'))
				for ll in range(tlen):
					self.filename['tci'][ll] = ftemp[ll]

			if (flag == 'ece'):
				dirs = '%s/DATASAVE/%i'%(self.mds_dir,self.shot)
				ftemp = self.make_ece_raw_files(float(self.note_md['e%i'%(tas_ind[k])].get())/1.e3,dirs)
				if ftemp['isfile']: 
					print('>>> There is ECE raw file')
					for ll in range(tlen):
						self.filename['te']['ece'][ll] = ftemp[ll]

		self.make_mds_e2()
		self.multi_time_fit()
		self.b2.configure(text='SAVE-M')
		print('>>> Run mode -> MULTI')
		return	

	def make_mds_g(self):

		try: self.shot  = int(self.note_md['e1'].get())
		except: print('>>> No shot #'); return
		os.system('%s %s %i %s'%(pythonc_exec,mds_over,self.shot,self.note_md['e8'].get()))
	
		return

	def make_mds_h(self):
		
		try: self.shot = int(self.note_md['e1'].get())
		except: print('>>> No shot #'); return
		dts = '%i,%i,%i'%(int(self.note_md['e3'].get()),int(self.note_md['e4'].get()),int(self.note_md['e6'].get()))
		os.system('%s %s %i %s 2d %s'%(pythonc_exec,mds_over,self.shot,self.note_md['e8'].get(),dts))


		return

	def make_mds_i(self):

		try: self.shot = int(self.note_md['e1'].get())
		except: print('>>> No shot #'); return
		dirs1 = '%s/DATASAVE/%i/'%(self.mds_dir,self.shot)
		times = np.array(self.note_md['e8'].get().split(','),dtype='int')
		self.fit.fit_opt['dacrit'] = float(self.note_md['e10'].get())
		self.get_elm_peak(dirs1,times)
		dirs2 = '%s/DATASAVE/%i/elm_peaks.npz'%(self.mds_dir,self.shot)
		tmin = min(times)-400; tmax = max(times)+400;
		os.system('%s %s %s %i %i %s'%(pythonc_exec,mds_da,dirs1+'/DA.npz',tmin,tmax,dirs2))

		return

	def make_mds_j(self):
		os.system('%s %i r &'%(dena_dir,int(self.note_md['e1'].get())))
		time.sleep(5)
		return

	def get_elm_peak(self,dirs,times):

		nfile = dirs+'/DA.npz'
		self.peaks = [];
		if (not os.path.isfile(nfile) and not os.path.isfile(nfile+'.gz')): print('>>> No DA signal'); return
		if not os.path.isfile(nfile):
			comm=gzip_dir+' -d '+nfile+'.gz'
			os.system(comm)
		npzfile  = np.load(nfile, mmap_mode='r')
		self.da = [];
		for i in range(2): self.da.append(npzfile['arr_%i'%i])
		comm=gzip_dir+' '+nfile
		os.system(comm)

		tmin = min(times)/1.e3 - 0.3
		tmax = max(times)/1.e3 + 0.3

		ind1 = np.where(self.da[0]>tmin); 		
		ind2 = np.where(self.da[0][ind1]<tmax);
		ind3 = int(round(0.003/(self.da[0][4]-self.da[0][3])))
		ind4 = int(round(0.00003/(self.da[0][4]-self.da[0][3])))

		tt = self.da[0][ind1][ind2]
		yy = self.da[1][ind1][ind2]
		lent = len(tt)

		mean_sig = np.mean(yy);

		if mean_sig>0: base = mean_sig*0.9; dbase = mean_sig*0.3;
		else: base = mean_sig*1.1; dbase = -mean_sig*0.3;
		dbase2 = np.max(yy)-base;
		dbase  = max(dbase,dbase2*0.2)

		ii=0
		while ii<lent:
			if (yy[ii]-base) > self.fit.fit_opt['dacrit']*dbase:
				self.peaks.append(tt[ii]*1.e3);
				ii += ind3
			else: ii += ind4
	
		np.savez(dirs+'/elm_peaks',self.peaks)
		return

	def get_duty(self,time):

		lenp = len(self.peaks)
		if  lenp< 3: return 1.

		for i in range(lenp-1):
			left = time - self.peaks[i]
			right = time - self.peaks[i+1]
			if left*right <=0.: break

		if i==(lenp-2): return 1.1

		duty = -100.*(right)/(self.peaks[i+1]-self.peaks[i]);

		if duty <5.: return 1.2;
		else: return duty
		
	def make_ts_ces_raw_files(self,dtime1,tflag,dirs,save_dir):

		scale = dict()
		scale['te'] = 1.e0
		scale['ne'] = 1.e-18
		scale['ti'] = 1.e0
		scale['vt'] = 1.e0
		filename = dict()
		if (tflag=='ts' or tflag=='tse'): flag_list=self.fit.prof_list[0:2]
		else: flag_list=self.fit.prof_list[2:4]
		for flag in flag_list:
			filename[flag] = dict()
			filename[flag]['isfile'] = True
			filen = dirs + '/%s_size'%flag.upper()
			if not os.path.isfile(filen):
				filename[flag]['isfile'] = False

				continue

			f = open(filen,'r')	
			stmp = f.readline().split()
			f.close()
			
			size = int(stmp[0])
			size2 = int(stmp[1])
			nbprobe = int(stmp[2])
			try: ncore   = int(stmp[3])
			except: ncore = nbprobe
			time = np.zeros(size); val = np.zeros(size); err = np.zeros(size); 
			tail = -1; tail2 = -1
			PROFV = dict(); PROFE = dict(); PROFR = dict();
			print('>>> %s -> '%flag.upper())
			for i in range(1,nbprobe+1):
				filen = dirs+'/%s%i.npz'%(flag.upper(),i)
				try:os.system('gzip -d %s.gz'%filen)
				except: pass
				npzfile = np.load(filen,mmap_mode='r')
				time = npzfile['arr_0']
				val  = npzfile['arr_1']
				err  = npzfile['arr_2']
				rr   = npzfile['arr_3']
				os.system('gzip %s'%filen)
				update_progress(float((i/nbprobe)))

				PROFV[i] = dict(); PROFE[i] = dict(); PROFR[i] = dict();
				for k in range(len(self.times)):
					time1 = self.times[k]
					ltime = time1 - dtime1 * 0.5
					utime = time1 + dtime1 * 0.5					
					ind1 = np.where((time*1000 >= ltime))
					ind2 = np.where((time[ind1]*1000 <= utime))
					tt   = time[ind1][ind2]
					ind3 = []; ind4 = [];
					for l in range(len(tt)):
						if tflag=='ces':
							duty = self.get_duty(tt[l]*1000.)
							if duty <= self.fit.fit_opt['duty']: ind3.append(l)
							ind4.append(l)
						else: ind3.append(l)
					if len(ind3) == 0: ind3 = ind4;
					PROFR[i][k] = rr * 1.e-3
					PROFV[i][k] = np.nan_to_num(val[ind1][ind2][ind3]) * scale[flag]
					PROFE[i][k] = np.nan_to_num(err[ind1][ind2][ind3]) * scale[flag]
				
			for k in range(len(self.times)):
				time1 = self.times[k]
				if (flag.lower()=='te' or flag.lower()=='ne'): 
					if (tflag=='ts'):    filename[flag][k] = save_dir+'/TS_%s_%ims.dat'%(flag.upper(),time1); istart = 1; iend = ncore+1;
					elif (tflag=='tse'): filename[flag][k] = save_dir+'/TSE_%s_%ims.dat'%(flag.upper(),time1); istart = ncore+1;  iend = nbprobe+1
					else:                filename[flag][k] = None; istart = 1; iend = -1;
				else:
					filename[flag][k] = save_dir+'/CES_%s_%ims.dat'%(flag.upper(),time1); istart = 1; iend = nbprobe+1
				if iend<0: continue
				f = open(filename[flag][k],'w')
				f.write('R[m]\tZ[m]\tVAL\tERR\n')
				tlen = len(PROFV[1])
				for i in range(len(PROFV[1][k])):
					for j in range(istart,iend):
						f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(PROFR[j][k],0.,PROFV[j][k][i],PROFE[j][k][i]))
				f.close()

		return filename,flag_list

	def make_tci_raw_files(self,dtime,dirs,save_dir):

	    NE_conv=[1.9,2.75,1.,1.,1.,1.,1.]
	    isdata =[0,0,0,0,0,0,0] 
	    size   =[0,0,0,0,0,0,0]
	    names2 =['ne_int1','ne_int2','ne_tci1','ne_tci2','ne_tci3','ne_tci4','ne_tci5']
	    tlen = len(self.times)
	    tciavg =np.zeros((tlen,7))

	    filename = dict()
	    for k in range(tlen): filename[k]   = ''

	    if not os.path.isfile(dirs+'/TCI_size'):
	        for t in range(tlen):
	            filename[t] = save_dir+'/TCI_%ims.dat'%self.times[t]
	            f = open(filename[t],'w')
	            for k in range(7):
	            	f.write('0. 0. \n')
	            f.close()	    	
	        return filename, False

	    f=open(dirs+'/TCI_size','r')
	    for k in range(7):
	        stmp=f.readline()
	        isdata[k]=int(stmp.split()[0])
	        size[k]=int(stmp.split()[1])
	    f.close()    
	    print('>>> %s -> '%'TCI')
	    for k in range(7):
	        if isdata[k] == 0:
	        	tciavg[:,k] = 0.
	        	continue
	        filen = dirs+'/'+names2[k]+'.npz'
	        comm='gzip -d '+filen+'.gz'
	        os.system(comm)
	        npzfile  = np.load(filen, mmap_mode='r')
	        time=npzfile['arr_0']
	        ne=npzfile['arr_1']
	        comm='gzip '+filen
	        os.system(comm)
	        lent = len(time); lenn = len(ne); lens = min(lent,lenn)
	        time = time[:lens]; ne = ne[:lens]

	        for t in range(tlen):
	        	utime = (self.times[t]+0.5*dtime)/1000.
	        	ltime = (self.times[t]-0.5*dtime)/1000.
	        	ind1 = np.where(time>ltime)
	        	ind2 = np.where(time[ind1]<utime)
	        	tciavg[t,k] = np.sum(ne[ind1][ind2])/len(ne[ind1][ind2])/NE_conv[k]
	        	ind3 = np.where(abs(ne[ind1][ind2]/NE_conv[k]-tciavg[t,k]) < 0.35)
	        	tciavg[t,k] = np.sum(ne[ind1][ind2][ind3])/len(ne[ind1][ind2][ind3])/NE_conv[k]

	        for t in range(tlen):
	            filename[t] = save_dir+'/TCI_%ims.dat'%self.times[t]
	            f = open(filename[t],'w')
	            for k in range(7):
	            	f.write('%i %4.2f\n'%(isdata[k],max(0.,tciavg[t,k])))
	            f.close()
	    update_progress(float(((k+1)/7)))
	    return filename,True

	def make_ece_raw_files(self,dtime1,dirs):

		try: shot = int(float(self.note_md['e1'].get()))
		except: shot = self.shot
		filen = dirs+'/ECE_size'
		ttime = self.times / 1.e3
		filename = dict()
		filename['isfile'] = False	
		if not os.path.isfile(filen):
			return filename
		filename['isfile'] = True
		tlen = len(self.times)
		times = self.times / 1.e3
		cdir = os.getcwd()
		os.chdir(self.mds_dir)
		load_ece(shot,times,dtime1,'../INPUTS/',None)
		for k in range(tlen):
			filename[k] = 'INPUTS/te_%06i_%ims_ECE_2nd.dat'%(shot,self.times[k])
			update_progress(float((k+1)/tlen))
				
		os.chdir(cdir)
		return	filename
	
	def multi_time_fit(self):

		try: os.mkdir('MPROFILES')
		except: pass

		if (self.didfit and (self.runmode=='single')):
			f = open('result_single.save','wb')
			pickle.dump([self.fit.fit_opt,self.fit.te_prof,self.fit.ne_prof,self.fit.ti_prof,self.fit.vt_prof,self.fit.post,self.fit.fit_eq],f)
			f.close()
		else:
			f = open('result_single.save','rb')
			##################################
			f.close()
		self.didfit = True

		self.run_fit_page()
		self.sync_gui_opt()
		self.put_fit_opt()

		self.mfit_opt = dict()
		self.mfit_opt['times'] = copy.deepcopy(self.times)
		self.mfit_opt['post'] = dict()

		for flag in self.fit.prof_list: self.fit.__dict__['%s_prof'%flag]['fit_old'] = []
		tlen = len(self.times)
		nel = np.zeros(7); bpped = np.zeros(2); width = np.zeros(8); ww = np.zeros(2); tci = np.zeros(7); tcie = np.copy(tci)
		for i in range(tlen):
			self.reset_input_file()
			self.fit.fit_opt['file']['gfile'] = self.filename['gfile'][i]
			if not self.filename['tci'][i] == None:
				f = open(self.filename['tci'][i],'r')
				for tflag in self.fit.inter_list:
					line = f.readline()
					self.fit.fit_opt[tflag]['val'] = round(float(line.split()[1]),2)
			for flag in self.fit.prof_list:
				for kk in self.fit.__dict__['%s_list'%flag]:
					self.fit.fit_opt['file'][flag][kk] = self.filename[flag][kk][i]
			print('>>> Fitting for %05i [ms]'%self.times[i])
			for flag in self.fit.prof_list: self.fit.post['didfit'][flag] = False
			self.fit.first_run = True
			self.fit.fit_opt['mds']['time'] = self.times[i]
			self.fit.main_run()
			self.fit.write_kinprof()
			self.fit.read_kinprof_kprofile()
			self.fit.post['time'] = self.times[i]
			self.post_opt['param_m'][self.times[i]] = copy.deepcopy(self.post_opt['param'])
			nel = [self.fit.post['int01'],self.fit.post['int02'],self.fit.post['tci01'],self.fit.post['tci02'],self.fit.post['tci03'],self.fit.post['tci04'],self.fit.post['tci05']]
			bpped = [self.fit.post['bpped'],min(self.fit.post['bppede'],self.fit.post['bpped'])]
			for k in range(4):
				width[k]   = self.fit.post['width'][k]
				width[k+4] = min(self.fit.post['width'][k],self.fit.post['widthe'][k])
			ww = [self.fit.post['wmhd'],self.fit.post['wkin']]
			for kk in range(7):
				tci[kk] = self.fit.fit_opt[self.fit.inter_list[kk]]['val']
				tcie[kk] = self.fit.fit_opt[self.fit.inter_list[kk]]['sig']
			copyfile('PROFILES/chease_kinprof_fit','MFIT/chease_kinprof_fit_%05i_temp'%self.times[i])
			try: copyfile('PROFILES/VT_fit.dat','MFIT/VT_fit.dat_%05i_temp'%self.times[i])
			except: pass
			self.mfit_opt[self.times[i]] = copy.deepcopy([self.fit.fit_opt,self.fit.te_prof,self.fit.ne_prof,self.fit.ti_prof,self.fit.vt_prof,self.fit.post,self.fit.fit_eq])
			self.mfit_opt['post'][self.times[i]] = copy.deepcopy([nel,bpped,width,ww,tci,tcie])
			try: rmtree('MFIT/%i'%self.times[i])
			except: pass
			copytree('PROFILES','MFIT/%i'%self.times[i])

		if self.plotpage > (tlen-1): self.plotpage = tlen-1
		[self.fit.fit_opt,self.fit.te_prof,self.fit.ne_prof,self.fit.ti_prof,self.fit.vt_prof,self.fit.post,self.fit.fit_eq] = copy.deepcopy(self.mfit_opt[self.times[self.plotpage]])
		self.restore_gui()
		for flag in self.fit.prof_list: self.fit.fit_opt['file'][flag]['kfile'] = None
		self.draw_plot()
		self.didmfit = True
		self.runmode = 'multi'
		self.fit.first_run = True
		return

	def restore_gui(self):
	
		self.get_fit_opt()
		self.sync_gui_opt(g2v=False)
		self.pedestal_option_a()
		for flag in self.fit.prof_list:
			self.fitting_func_option_a(flag,False)
		self.fit_constraint_option_a()
		self.run_post()
		return

	def restore_singlefit(self):
		if self.runmode == 'single': return
		if not os.path.isfile('result_single.save'):	
			print('>>> No prev. single run')
			return
		self.runmode = 'single'
		self.didmfit = False
		self.load_run_save('result_single.save',False)

		temp = self.note_md['e8'].get()
		self.restore_gui()
		self.plot_type = 0
		self.note_md['e8'].delete(0,'end')
		self.note_md['e8'].insert(10,temp)
		self.note_md['times'].set(temp)
		self.fit.fit_opt['mds']['times'] = temp
		self.draw_plot()
		self.b2.configure(text='SAVE-S')
		print('>>> Run mode MULTI -> SINGLE')
		return

	def restore_multifit(self):
		if self.runmode == 'multi': return
		if len(self.note_md['e8'].get()) == 0: return
		self.times = np.array(self.note_md['e8'].get().split(','),dtype='int')
		isfile = True
		if not os.path.isfile('result_multi.save'): isfile = False
		if not isfile: 
			print('>>> No saved mutli run')
			return
		self.plot_type = 0	
		self.didmfit = True
		self.runmode = 'multi'
		self.load_run_save('result_multi.save',True)	
		self.restore_gui()
		self.times = np.array(self.note_md['e8'].get().split(','),dtype='int')
		for flag in self.fit.prof_list: self.fit.fit_opt['file'][flag]['kfile'] = None
		self.draw_plot()
		self.b2.configure(text='SAVE-M')
		print('>>> Run mode SINGLE -> MULTI')
		return

	def restore_fitopt(self):

		if self.runmode=='multi':
			if not os.path.isfile('result_multi.save'): return
			else: 
				self.load_run_save('result_multi.save',True)
				self.restore_gui()
				self.plot_type = 0
				self.draw_plot()
		else:
			if not os.path.isfile('result_single.save'): return
			else:
				self.load_run_save('result_single.save',False)
				temp = self.note_md['e8'].get()
				self.restore_gui()
				self.plot_type = 0
				self.note_md['e8'].delete(0,'end')
				self.note_md['e8'].insert(10,temp)
				self.note_md['times'].set(temp)
				self.fit.fit_opt['mds']['times'] = temp
				self.draw_plot()			

		return

	def check_runmode(self):
		if not os.path.isfile('run_mode'): 
			self.runmode = 'single'
			print('>>> Initial run, run mode -> SINGLE')
			return
	
		f = open('run_mode','r')
		line=f.readline().split()
		f.close()
		if line[0].upper() == 'SINGLE': 
			if not os.path.isfile('result_single.save'):
				self.runmode = 'single'
				print('>>> Initial run, run mode -> SINGLE')
				return				
			self.didfit = True
			self.load_run_save('result_single.save',False)
			self.runmode = 'single'
			self.restore_gui()
			self.b2.configure(text='SAVE-S')

		elif line[0].upper() == 'MULTI':
			if not os.path.isfile('result_multi.save'):
				self.runmode = 'single'
				print('>>> Initial run, run mode -> SINGLE')
				return		

			self.didmfit = True
			self.load_run_save('result_multi.save',True)	
			for flag in self.fit.prof_list: self.fit.fit_opt['file'][flag]['kfile'] = None
			self.runmode = 'multi'
			self.restore_gui()
			self.times = np.array(self.note_md['e8'].get().split(','),dtype='int')

			self.b2.configure(text='SAVE-M')

		print('>>> Current run mode = %s'%self.runmode.upper())
		self.draw_plot()
		return

	def load_run_save(self,file,ismult=False):

		f = open(file,'rb')
		if not ismult: 
			[fit_opt,self.fit.te_prof,self.fit.ne_prof,self.fit.ti_prof,self.fit.vt_prof,post,self.fit.fit_eq] = pickle.load(f)
		else:
			self.mfit_opt = dict()
			temp = pickle.load(f)
			for key in temp.keys(): self.mfit_opt[key] = temp[key]
		f.close()

		if ismult:
			[fit_opt,self.fit.te_prof,self.fit.ne_prof,self.fit.ti_prof,self.fit.vt_prof,post,self.fit.fit_eq] = self.mfit_opt[self.mfit_opt['times'][0]]
		for key in post.keys(): self.fit.post[key] = post[key]
		self.load_fitopt_dict(fit_opt)

		return 

	def back_(self):
		if not self.didmfit: return
		if self.plotpage == 0: 
			if self.plot_type == 0: return
			else: self.plotpage = self.plotpage + 1
		if len(self.times) == 1: return
		if self.plot_type > 2: self.plot_type = 0
		self.plotpage = self.plotpage - 1
		[fit_opt,self.fit.te_prof,self.fit.ne_prof,self.fit.ti_prof,self.fit.vt_prof,post,self.fit.fit_eq] = copy.deepcopy(self.mfit_opt[self.mfit_opt['times'][self.plotpage]])
		for key in post.keys(): self.fit.post[key] = post[key]
		self.load_fitopt_dict(fit_opt)		

#		self.restore_gui()
		if self.mfit_opt['times'][self.plotpage] in self.post_opt['param_m'].keys():
			self.post_opt['param'] = copy.deepcopy(self.post_opt['param_m'][self.mfit_opt['times'][self.plotpage]])
			for flag in self.fit.prof_list: self.put_param_gui(flag)
		self.restore_gui()
		self.draw_plot()
		return

	def next_(self):
		if not self.didmfit: return
		if self.plotpage == (len(self.times)-1): 
			if self.plot_type == 0: return
			else: self.plotpage = self.plotpage - 1
		if len(self.times) == 1: return
		if (self.plot_type>2): self.plot_type = 0
		self.plotpage = self.plotpage + 1
		[fit_opt,self.fit.te_prof,self.fit.ne_prof,self.fit.ti_prof,self.fit.vt_prof,post,self.fit.fit_eq] = copy.deepcopy(self.mfit_opt[self.mfit_opt['times'][self.plotpage]])
		for key in post.keys(): self.fit.post[key] = post[key]
		self.load_fitopt_dict(fit_opt)		
#		self.restore_gui()
		if self.mfit_opt['times'][self.plotpage] in self.post_opt['param_m'].keys():
			self.post_opt['param'] = copy.deepcopy(self.post_opt['param_m'][self.mfit_opt['times'][self.plotpage]])
			for flag in self.fit.prof_list: self.put_param_gui(flag)
		self.restore_gui()
		self.draw_plot()
		return		

	def reset_input_file(self):

		self.fit.fit_opt['file']['gfile']     = None
		for i in self.fit.inter_list:
			self.fit.fit_opt[i]['val'] = 0.0
			#self.fit.fit_opt[i]['sig'] = 0.0		

		for i in self.fit.prof_list:	
			self.fit.fit_opt['file'][i]['kfile'] = None
			for j in self.fit.__dict__['%s_list'%i]:
				self.fit.fit_opt['file'][i][j] = None		

		return

	def make_plot(self):

		self.l1 = tk.Label(self.home['page8'], text="================== Plot Axis Option ==================",justify='center')
		self.l1.grid(row=0, column=0,columnspan=9,pady=5)

		titles = [" TE "," NE "," TI "," VT "]
		count = 1
		for flag in ['xmin','xmax','ymin','ymax']:
			self.l2 = tk.Label(self.home['page8'], text=flag.upper())
			self.l2.grid(row=1,column=count,columnspan=2)
			count = count + 2

		count = 1; count1 = 1;
		for flag in self.fit.prof_list:
			self.l2 = tk.Label(self.home['page8'], text=' '+flag.upper()+' ')
			self.l2.grid(row=count+1,column=0)
			count2 = 1
			for kk in ['xmin','xmax','ymin','ymax']:
				self.note_pl['e%i'%count1] = tk.Entry(self.home['page8'],width=7,justify='center')
				self.note_pl['e%i'%count1].grid(row=count+1,column=count2,columnspan=2)
				self.note_pl['e%i'%count1].insert(10,self.note_pl[flag][kk].get())
				count2 = count2 + 2
				count1 = count1 + 1

			count = count + 1

		b1 = tk.Button(self.home['page8'],text="Update", height = 1, width = 5,command=lambda: self.make_plot_a(False))	
		b1.grid(row=count1+1, column=7,columnspan=2,pady=5)

		self.l1 = tk.Label(self.home['page8'], text="==================== Multi Plots ====================",justify='center')
		self.l1.grid(row=count1+2, column=0,columnspan=9,pady=5)
		b1 = tk.Button(self.home['page8'],text="Ti & Te overlap", height = 1, width = 40,command=lambda: self.make_plot_c())	
		b1.grid(row=count1+3, column=0,columnspan=9,pady=5)
		b2 = tk.Button(self.home['page8'],text="Multi Time Fit overlap", height = 1, width = 40,command=lambda: self.make_plot_d())	
		b2.grid(row=count1+4, column=0,columnspan=9,pady=5)		
		b3 = tk.Button(self.home['page8'],text="Multi Time Fit 0D params", height = 1, width = 40,command=lambda: self.make_plot_e())	
		b3.grid(row=count1+5, column=0,columnspan=9,pady=5)			

		count1 = count1+6

		if self.single_only:
			b2.configure(state='disabled')
			b3.configure(state='disabled')
	
		return

	def make_plot_a(self,restore=True):
	
		if restore: 
			if self.plot_type > 0:
				self.plot_type = 0 
				self.fit.draw_plot(fig_ex=self.figure['name']['home'],second=False)

		axes = self.figure['name']['home'].axes
		self.figure['name']['home'].canvas.draw_idle()
		count1 = 1
		for i in range(4):
			ax = axes[i]
			flag = self.fit.prof_list[i]
			for kk in ['xmin','xmax','ymin','ymax']:
				self.fit.fit_opt['plot'][flag][kk] = float(self.note_pl['e%i'%count1].get())
				count1 = count1 +1

			ax.set_xlim(self.fit.fit_opt['plot'][flag]['xmin'],self.fit.fit_opt['plot'][flag]['xmax'])
			if not ((self.fit.fit_opt['plot'][flag]['ymin'] == -1) and (self.fit.fit_opt['plot'][flag]['ymax'] == -1)):
				ax.set_ylim(self.fit.fit_opt['plot'][flag]['ymin'],self.fit.fit_opt['plot'][flag]['ymax'])

		if self.runmode =='multi':
			for time in self.mfit_opt['times']:
				self.mfit_opt[time][0]['plot'] = copy.deepcopy(self.fit.fit_opt['plot'])


		return	

	def make_plot_b(self):

		self.plot_type = 1
		self.figure['name']['home'].canvas.draw_idle()
		self.fit.draw_plot(fig_ex=self.figure['name']['home'],second=True)		

		return

	def make_plot_c(self):

		if self.plot_type >0: return
		self.plot_type = 2
		self.figure['name']['home'].canvas.draw_idle()
		axes = self.figure['name']['home'].axes

		if (self.fit.post['didfit']['te'] and self.fit.post['didfit']['ti']):
			xx = self.fit.fit_eq['psin2']
			if self.fit.fit_opt['use_rho']['te']: xx = self.fit.fit_eq['psi_to_rho'](xx)
			axes[0].plot(xx,self.fit.ti_prof['fit2p'],color='steelblue',linestyle='--')

		return

	def make_plot_d(self):

		if not self.didmfit: return
		self.plot_type = 3
		self.figure['name']['home'].canvas.draw_idle()
		axes = self.figure['name']['home'].axes
		tlen = len(self.times); plegend = []; slegend = [];

		for i in range(tlen): slegend.append('t = %i [ms]'%self.times[i])
		for k in range(4): axes[k].cla()
		for time in self.mfit_opt['times']:
			[fit_opt,te_prof,ne_prof,ti_prof,vt_prof,post,fit_eq] = self.mfit_opt[time]
			for k in range(4):
				flag = self.fit.prof_list[k]
				if (flag == 'te'):
					title = '$T_e$ [keV]'
				elif (flag == 'ne'):
					title = '$n_e$ [$10^{19}$/m3]'
				elif (flag == 'ti'):
					title = '$T_i$ [keV]'
				elif (flag == 'vt'):
					title = '$V_\phi$ [km/s]'				

				WPED = 0.
				if (fit_opt['func_type'][flag] == 4 or fit_opt['func_type'][flag] == 2):
					WPED = post['popt'][flag][2]
					PMID = 1.0 - 0.5*WPED
					if fit_opt['use_rho'][flag]: WPED,PMID = self.fit.pwidth2rwidth(WPED,PMID)

				if (fit_opt['func_type'][flag] ==6):
					WPED = post['popt'][flag][2]
					PMID = post['popt'][flag][6]
					if fit_opt['use_rho'][flag]: WPED,PMID = self.fit.pwidth2rwidth(WPED,PMID)

				if (fit_opt['func_type'][flag] == 3):
					WPED = post['popt'][flag][5]
					PMID = 1.0-0.5*WPED		
					if fit_opt['use_rho'][flag]: WPED,PMID = self.fit.pwidth2rwidth(WPED,PMID)	

				axes[k].set_ylabel(title)
				if (WPED > 0.): 
					title = title + ' - $W_{\psi,ped}$ = %4.3f'%WPED
				
				axes[k].set_title(title)
				if not (fit_opt['use_rho'][flag]):
					axes[k].set_xlabel('$\psi_n$ [a.u]')
				else:
					axes[k].set_xlabel('$\\rho_t$ [a.u]')

				if post['didfit'][flag]: 
					if (flag == 'te'): axes[k].plot(fit_eq['psin2'],te_prof['fit2'],linestyle='--')
					if (flag == 'ne'): axes[k].plot(fit_eq['psin2'],ne_prof['fit2'],linestyle='--')
					if (flag == 'ti'): axes[k].plot(fit_eq['psin2'],ti_prof['fit2'],linestyle='--')
					if (flag == 'vt'): axes[k].plot(fit_eq['psin2'],vt_prof['fit2'],linestyle='--')

				else: axes[k].text(0.5,0.5,'No data',color='red')

				axes[k].legend(slegend)

		return

	def make_plot_e(self):
		colorm = ['blue','orange','green','red','purple','brown','gray']

		if not self.didmfit: return
		self.plot_type = 4
		self.figure['name']['home'].canvas.draw_idle()
		axes = self.figure['name']['home'].axes
		tlen = len(self.times); plegend = []; slegend = [];

		for k in range(4): 
			axes[k].cla()
			axes[k].set_xlabel('Times [ms]')
			axes[k].set_xlim(min(self.times)-100,max(self.times)+100)
		axes[0].set_title('$<n_e> [10^{19}m^{-3}]$')
		axes[1].set_title('$W_{stored}$ [kJ]')
		axes[2].set_title('$\\beta_{ped}$ [a.u]')
		axes[3].set_title('$W_{ped,\psi}$ [a.u]')

		tlen = len(self.mfit_opt['times'])
		nel = np.zeros((tlen,7)); tci = np.zeros((tlen,7)); tcie = np.zeros((tlen,7))
		bpped = np.zeros((tlen,2)); width = np.zeros((tlen,8)); ww = np.zeros((tlen,2));
		for k in range(tlen):
			[nel[k,:],bpped[k,:],width[k,:],ww[k,:],tci[k,:],tcie[k,:]]=self.mfit_opt['post'][self.times[k]]

		no_tci = True
		plegend = []
		tlegend = []
		for k in range(7):
			if not tci[0,k]*tcie[0,k] == 0.:
				line1, = axes[0].plot(self.times,nel[:,k],'x-',color=colorm[k])
				line2  = axes[0].errorbar(self.times,tci[:,k],tcie[:,k],fmt='o--',c=colorm[k])
				tlegend.append('%s - FIT'%self.fit.inter_list[k].upper())
				tlegend.append('%s - EXP.'%self.fit.inter_list[k].upper())
				plegend.append(line1)
				plegend.append(line2)
				no_tci = False
		if no_tci:
			axes[0].plot(self.times,nel,'o--')
			tlegend = copy.deepcopy(self.fit.inter_list)

			axes[0].legend(tlegend)

		else: axes[0].legend(plegend,tlegend)
		axes[1].plot(self.times,ww,'o--')
		axes[1].legend(['$W_{MHD}$','$W_{th}$'])
		axes[2].errorbar(self.times,bpped[:,0],bpped[:,1],fmt='o--')
		for i in range(4):
			axes[3].errorbar(self.times,width[:,i],width[:,i+4],fmt='o--')

		llegend = []
		for i in range(4):llegend.append('%s'%self.fit.prof_list[i].upper())

		axes[3].legend(llegend)
		return

	def make_etc(self):
		count1 = 0
		self.l1 = tk.Label(self.home['page9'], text="==================== Raw Scale Opt. ====================",justify='center')
		self.l1.grid(row=count1, column=0,columnspan=9,pady=5)
		titles = ["TCORE","TEDGE","NCORE","NEDGE"]
		for i in range(4):
			self.l2 = tk.Label(self.home['page9'],text=titles[i].upper())
			self.l2.grid(row=count1+1,column=2*i+1,columnspan=2,pady=5)

		self.l2 = tk.Label(self.home['page9'], text='TS')
		self.l2.grid(row=count1+2,column=0)

		count = 1; i = 0
		for flag in ['te','ne']:
			for k in ['core','edge']:
				self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
				self.note_ec['e%i'%count].grid(row=count1+2,column=2*i+1,columnspan=2)
				self.note_ec['e%i'%count].insert(10,self.note_ec[flag]['scale']['ts'][k].get())
				count = count + 1			
				i = i + 1

		count1 = 2
		titles = ["TCORE","TEDGE","VCORE","VEDGE"]
		for i in range(4):
			self.l2 = tk.Label(self.home['page9'],text=titles[i].upper())
			self.l2.grid(row=count1+1,column=2*i+1,columnspan=2,pady=5)		
		self.l2 = tk.Label(self.home['page9'], text='CES')
		self.l2.grid(row=count1+2,column=0)

		i = 0;
		for flag in ['ti','vt']:
			for k in ['core','edge']:
				self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
				self.note_ec['e%i'%count].grid(row=count1+2,column=2*i+1,columnspan=2)
				self.note_ec['e%i'%count].insert(10,self.note_ec[flag]['scale']['ces'][k].get())
				count = count + 1			
				i = i + 1

		count1 = 4
		titles = ["ECE","REFL"]
		for i in range(2):
			self.l2 = tk.Label(self.home['page9'],text=titles[i].upper())
			self.l2.grid(row=count1+1,column=2*i+1,columnspan=2,pady=5)		
		self.l2 = tk.Label(self.home['page9'], text='ETC')
		self.l2.grid(row=count1+2,column=0)

		for i in range(2):
				self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
				self.note_ec['e%i'%count].grid(row=count1+2,column=2*i+1,columnspan=2)
				if i==0: self.note_ec['e%i'%count].insert(10,self.note_ec['te']['scale']['ece']['core'].get())
				else: self.note_ec['e%i'%count].insert(10,self.note_ec['ne']['scale']['refl']['core'].get())
				count = count + 1			


		count1 = 7;
		self.l2 = tk.Label(self.home['page9'], text='AUTO-S')
		self.l2.grid(row=count1+2,column=0)		
		titles = ["MINC","MAXC","MINE","MAXE"]
		for i in range(4):
			self.l2 = tk.Label(self.home['page9'],text=titles[i].upper())
			self.l2.grid(row=count1+1,column=2*i+1,columnspan=2,pady=5)

		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*0+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['tscale']['minc'].get())
		count = count + 1
		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*1+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['tscale']['maxc'].get())
		count = count + 1
		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*2+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['tscale']['mine'].get())
		count = count + 1
		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*3+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['tscale']['maxe'].get())		
		count = count + 1

		count1 = 9;
		titles = ["CN#","EN#","USE","1D"]
		for i in range(4):
			self.l2 = tk.Label(self.home['page9'],text=titles[i].upper())
			self.l2.grid(row=count1+1,column=2*i+1,columnspan=2,pady=5)

		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*0+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['tscale']['cn'].get())
		count = count + 1
		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*1+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['tscale']['en'].get())
		count = count + 1

		self.note_ec['c1'] = tk.Checkbutton(self.home['page9'],variable=self.note_ec['tscale']['is'])
		self.note_ec['c1'].grid(row=count1+2,column=2*2+1,columnspan=2)
		self.note_ec['c3'] = tk.Checkbutton(self.home['page9'],variable=self.note_ec['tscale']['1d'])
		self.note_ec['c3'].grid(row=count1+2,column=2*3+1,columnspan=2)		

		count1 = 12; 
		self.l1 = tk.Label(self.home['page9'], text="==================== Radial Shift Opt. ====================",justify='center')
		self.l1.grid(row=count1, column=0,columnspan=9,pady=5)
		self.l2 = tk.Label(self.home['page9'], text='SHIFT[m] →')
		self.l2.grid(row=count1+2,column=0)		
		titles = ["TS","CES","ECE","REFL"]
		for i in range(4):
			self.l2 = tk.Label(self.home['page9'],text=titles[i].upper())
			self.l2.grid(row=count1+1,column=2*i+1,columnspan=2,pady=5)

		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*0+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['shift']['ts'].get())
		count = count + 1
		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*1+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['shift']['ces'].get())
		count = count + 1
		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*2+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['shift']['ece'].get())
		count = count + 1
		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*3+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['shift']['ref'].get())
		count = count + 1				


		count1 = 14;
		self.l2 = tk.Label(self.home['page9'], text='AUTO-S')
		self.l2.grid(row=count1+2,column=0)		
		titles = ["MINP","MAXP","USE"]
		for i in range(3):
			self.l2 = tk.Label(self.home['page9'],text=titles[i].upper())
			self.l2.grid(row=count1+1,column=2*i+1,columnspan=2,pady=5)

		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*0+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['ashift']['min'].get())
		count = count + 1
		self.note_ec['e%i'%count] = tk.Entry(self.home['page9'],width=7,justify='center')
		self.note_ec['e%i'%count].grid(row=count1+2,column=2*1+1,columnspan=2)
		self.note_ec['e%i'%count].insert(10,self.note_ec['ashift']['max'].get())
		count = count + 1
		self.note_ec['c2'] = tk.Checkbutton(self.home['page9'],variable=self.note_ec['ashift']['is'])
		self.note_ec['c2'].grid(row=count1+2,column=2*2+1,columnspan=2)

		return

	def detect_close(self,ind,frame):

		frame.destroy()
		self.__dict__['t_close%i'%ind] = True
		if ind == 1: self.make_exclude()

		return

	def print_tab(self):

		print('aaaaaa')
		return

	def initialise_setup(self):

		self.pedestal_option_a()
		print('>>> Initialising env.')
		for flag in self.fit.prof_list:
			self.fitting_func_option_a(flag)

		self.fit_constraint_option_a(True)
		return

	def declare_variables(self):

		self.figure     = dict()
		self.figure['name']   = dict()
		self.figure['size']   = dict()
		self.figure['canvas'] = dict()
		self.figure['widget'] = dict()
		self.figure['toolbar'] = dict()

		self.exclude    = dict()

		self.home       = dict()
		self.note_in    = dict()
		self.note_fn    = dict()
		self.note_md    = dict()
		self.note_te    = dict()
		self.note_ne    = dict()
		self.note_ti    = dict()
		self.note_vt    = dict()
		self.note_pl    = dict()
		self.note_ec    = dict()

		self.post_opt   = dict()

		self.fit.prof_list = ['te','ne','ti','vt']
		self.func_list = ['CORE','MTANH','PTANH','EPED','SPLINE','EPED2','NSPLINE']
		self.initialise_opt = True
		self.didfit = False
		self.plot_type = 0
		self.post_times = ''
		self.mds_load = dict()

		for i in range(4): self.__dict__['t_close%i'%i] = True

		return

	def initialise_variables(self):

		self.curr_page = 0
		self.prev_page = 0
		self.figure['size']['home'] = (11,7.7)
		self.figure['size']['plot'] = (3.8,0.9)

		for k in ['tst','tsn','cest','cesv','tci','ece','ref']: self.mds_load[k] = False

		self.ped_scan_fit = tk.IntVar()
		self.ti_ped_scan  = tk.IntVar()
		self.ti_ped_width = tk.IntVar()		

		self.note_fn['raw_fit'] = dict()
		self.note_fn['sep_fix'] = dict()
		self.note_fn['sep_val'] = dict()
		self.note_fn['width_fix'] = dict()
		self.note_fn['width_val'] = dict()
		self.note_fn['func_type'] = dict()

		self.note_fn['func_type']['etc'] = tk.StringVar()
		self.note_fn['func_type']['etc'].set('----')
		self.note_in['file'] = dict()
		self.note_in['file']['gfile']       = tk.StringVar()
		self.note_in['file']['kfile']       = dict()
		self.note_in['zeff'] = tk.StringVar()
		self.note_in['zimp'] = tk.StringVar()
		self.note_in['amain'] = tk.StringVar()
		self.note_in['aimp']  = tk.StringVar()

		for k in self.fit.inter_list:
			self.note_in[k.lower()] = tk.StringVar()
			self.note_in[k.lower()+'s'] = tk.StringVar()
			self.note_in[k.lower()+'d'] = tk.StringVar()
			self.note_in[k.lower()+'f'] = tk.IntVar()

		self.note_fn['use_rho'] = dict()
		self.note_fn['smooth'] = dict()
		self.note_fn['olcn'] = dict()
		
		for flag in self.fit.prof_list:
			self.note_fn['raw_fit'][flag]   = tk.IntVar()
			self.note_fn['sep_fix'][flag]   = tk.IntVar()
			self.note_fn['sep_val'][flag]   = tk.StringVar()
			self.note_fn['width_fix'][flag] = tk.IntVar()
			self.note_fn['width_val'][flag] = tk.StringVar()
			self.note_fn['func_type'][flag] = tk.StringVar()

			self.note_in['file'][flag]   = dict()
			self.note_in['file']['kfile'][flag] = tk.StringVar()
			self.note_fn['use_rho'][flag] = tk.IntVar()
			self.note_fn['smooth'][flag] = tk.StringVar()
			self.note_fn['olcn'][flag] = tk.StringVar()

			for k in self.fit.__dict__['%s_list'%flag]:
				self.note_in['file'][flag][k] = tk.StringVar()

		for flag in self.fit.prof_list:
			self.__dict__['note_%s'%flag]['param'] = dict()
			for k in range(10):
				self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)] = dict()
				self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['min'] = tk.StringVar()
				self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['max'] = tk.StringVar()
				self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['val'] = tk.StringVar()
				self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['fix'] = tk.IntVar()

			self.__dict__['note_%s'%flag]['avg'] = dict()
			self.__dict__['note_%s'%flag]['oli'] = dict()
			self.__dict__['note_%s'%flag]['std'] = dict()
			self.__dict__['note_%s'%flag]['estd'] = dict()

			self.__dict__['note_%s'%flag]['olcn'] = dict()
			self.__dict__['note_%s'%flag]['range'] = dict()
			self.__dict__['note_%s'%flag]['olc'] = dict()
			self.__dict__['note_%s'%flag]['weight'] = dict()

			for k in self.fit.__dict__['%s_list'%flag]:
				self.__dict__['note_%s'%flag]['avg'][k] = tk.IntVar()
				self.__dict__['note_%s'%flag]['oli'][k] = tk.IntVar()
				self.__dict__['note_%s'%flag]['std'][k] = tk.IntVar()
				self.__dict__['note_%s'%flag]['estd'][k] = tk.IntVar()				

				self.__dict__['note_%s'%flag]['olcn'][k] = tk.StringVar()
				self.__dict__['note_%s'%flag]['range'][k] = tk.StringVar()
				self.__dict__['note_%s'%flag]['olc'][k]   = tk.StringVar()
				self.__dict__['note_%s'%flag]['weight'][k]= tk.StringVar()

		self.post_opt['func_type'] = dict()
		self.post_opt['width_fix'] = dict()
		self.post_opt['width_val'] = dict()
		self.post_opt['param']     = dict()
		self.post_opt['param_m']   = dict()
		self.post_opt['oli']       = dict()
		self.post_opt['estd']      = dict()
		for flag in self.fit.prof_list:
			self.post_opt['func_type'][flag] = None
			self.post_opt['width_fix'][flag] = None
			self.post_opt['width_val'][flag] = None
			self.post_opt['oli'][flag]       = dict()
			self.post_opt['estd'][flag]      = dict()
			for kk in self.fit.__dict__['%s_list'%flag]:
				self.post_opt['oli'][flag][kk] = None
				self.post_opt['estd'][flag][kk]= None

		for flag in self.fit.prof_list:
			self.post_opt['param'][flag] = dict()
			for ff in self.fit.func_list:
				self.post_opt['param'][flag][ff.lower()]        = dict()
				self.post_opt['param'][flag][ff.lower()]['is']  = False
				self.post_opt['param'][flag][ff.lower()]['val'] = [] 

		for flag in self.fit.prof_list:
			self.note_pl[flag] = dict()
			self.note_pl[flag]['xmin'] = tk.StringVar()
			self.note_pl[flag]['xmax'] = tk.StringVar()
			self.note_pl[flag]['ymin'] = tk.StringVar()
			self.note_pl[flag]['ymax'] = tk.StringVar()

		self.note_ec = dict()
		for flag in self.fit.prof_list:
			self.note_ec[flag] = dict()
			self.note_ec[flag]['scale'] = dict()
			for kk in self.fit.__dict__['%s_list'%flag]:
				self.note_ec[flag]['scale'][kk] = dict()
				self.note_ec[flag]['scale'][kk]['core'] = tk.StringVar()
				self.note_ec[flag]['scale'][kk]['edge'] = tk.StringVar()
		
		self.note_ec['tscale'] = dict()
		self.note_ec['tscale']['minc'] = tk.StringVar()
		self.note_ec['tscale']['maxc'] = tk.StringVar()
		self.note_ec['tscale']['mine'] = tk.StringVar()
		self.note_ec['tscale']['maxe'] = tk.StringVar()
		self.note_ec['tscale']['cn'] = tk.StringVar()
		self.note_ec['tscale']['en'] = tk.StringVar()		
		self.note_ec['tscale']['is'] = tk.IntVar()
		self.note_ec['tscale']['1d'] = tk.IntVar()

		self.note_ec['shift'] = dict()
		self.note_ec['shift']['ces']   = tk.StringVar()
		self.note_ec['shift']['ts']    = tk.StringVar()
		self.note_ec['shift']['ref']   = tk.StringVar()
		self.note_ec['shift']['ece']   = tk.StringVar()
		self.note_ec['ashift'] = dict()
		self.note_ec['ashift']['is'] = tk.IntVar()
		self.note_ec['ashift']['min'] = tk.StringVar()
		self.note_ec['ashift']['max'] = tk.StringVar()

		self.note_md['cmse'] = tk.IntVar()
		self.note_md['efit_type'] = tk.StringVar()
		self.note_md['ece_hfs'] = tk.IntVar()
		self.note_md['shot'] = tk.StringVar()
		self.note_md['time'] = tk.StringVar()
		self.note_md['times']= tk.StringVar()
		self.note_md['dt']   = dict()
		self.note_md['dt']['ts'] = tk.StringVar()
		self.note_md['dt']['ces']= tk.StringVar()
		self.note_md['dt']['tci']= tk.StringVar()
		self.note_md['dt']['ece']= tk.StringVar()
		self.note_md['dt']['ref']= tk.StringVar()
		self.note_md['dt']['ref2']=tk.StringVar()

		self.note_md['dacrit'] = tk.StringVar()
		self.note_md['duty']   = tk.StringVar()

		self.didmfit = False
		self.runmode = 'single'
		self.plotpage = 0

		self.mds_dir  = './MDS'
		return

	def transfer_logic(self,var1,g2f,inv=False):

		if g2f:
			if var1.get() == 1: var2 = True
			else:	var2 = False
			if inv: var2 = not var2
		else:
			if var1: var2 = 1
			else: var2 = 0

		return var2

	def transfer_int(self,var1,g2f):

		if g2f: var2 = int(float(var1.get()))
		else: var2 = '%i'%var1
		return var2

	def transfer_float(self,var1,g2f):
		if g2f: var2 = float(var1.get())
		else: var2 = '%5.2f'%var1
		return var2		

	def transfer_name(self,var1,g2f):
		if g2f: 
			if len(var1.get().split()) == 0: var2 = None
			else: var2 = var1.get()

		else: 
			if var1 == None: var2 = ''
			else: var2 = var1
		return var2		

	def transfer_func(self,var1,g2f):

		len2 = len(self.func_list)
		if g2f:
			for i in range(len2):
				flag = self.func_list[i].lower()
				if var1.get().lower() == flag:
					var2 = i+1
		else:
			var2 = self.func_list[var1-1].upper()
		return var2

	def transfer_param(self,var1,g2f):

		if g2f:
			var3 = var1.get().lower()
			if  (var3 == '-inf'): var2 = np.inf
			elif(var3 == 'inf'):  var2 = -np.inf
			else: var2 = float(var1.get())

		else:
			var2 = str(var1)

		return var2

	def put_fit_opt(self):

		self.fit.fit_opt['use_ti_width']      = self.transfer_logic(self.ti_ped_width,True)
		self.fit.fit_opt['ped_scan_fit']      = self.transfer_logic(self.ped_scan_fit,True)
		self.fit.fit_opt['use_ti_eped']       = self.transfer_logic(self.ti_ped_scan,True)
		self.fit.fit_opt['line_type']         = 0.
		self.fit.fit_opt['target_density']    = 0.
		self.fit.fit_opt['use_density_scale'] = False
		for i in range(7): 
			flag = self.fit.inter_list[i]
			if self.note_in[flag+'f'].get() == 1:	
				self.fit.fit_opt['line_type'] = i + 1
				self.fit.fit_opt['target_density'] = self.transfer_float(self.note_in[flag],True)
				self.fit.fit_opt['use_density_scale'] = True

		self.fit.fit_opt['file']['gfile']        = self.transfer_name(self.note_in['file']['gfile'],True)

		for i in self.fit.prof_list:	
			self.fit.fit_opt['use_rho'][i]       = self.transfer_logic(self.note_fn['use_rho'][i],True)
			
			self.fit.fit_opt['sep_fix'][i]       = self.transfer_logic(self.note_fn['sep_fix'][i],True)
			self.fit.fit_opt['sep_val'][i]       = self.transfer_float(self.note_fn['sep_val'][i],True)
			self.fit.fit_opt['width_fix'][i]     = self.transfer_logic(self.note_fn['width_fix'][i],True)
			self.fit.fit_opt['width_val'][i]     = self.transfer_float(self.note_fn['width_val'][i],True)
			self.fit.fit_opt['raw_fit'][i]       = self.transfer_logic(self.note_fn['raw_fit'][i],True)
			self.fit.fit_opt['sspline_order'][i] = self.transfer_int(self.note_fn['smooth'][i],True)

			self.fit.fit_opt['file'][i]['kfile'] = self.transfer_name(self.note_in['file']['kfile'][i],True)
			count = 0;
			for j in self.fit.__dict__['%s_list'%i]:

				self.fit.fit_opt['file'][i][j]		 = self.transfer_name(self.note_in['file'][i][j],True)
				self.fit.fit_opt['avg'][i][j]		 = self.transfer_logic(self.__dict__['note_%s'%i]['avg'][j],True)
				if count==0: self.fit.fit_opt['oli'][i]['n']      = self.transfer_int(self.__dict__['note_%s'%i]['olcn'][j],True)
				self.fit.fit_opt['oli'][i][j]['use'] = self.transfer_logic(self.__dict__['note_%s'%i]['oli'][j],True)
				self.fit.fit_opt['oli'][i][j]['per'] = self.transfer_float(self.__dict__['note_%s'%i]['olc'][j],True)
				self.fit.fit_opt['std'][i][j]['use'] = self.transfer_logic(self.__dict__['note_%s'%i]['std'][j],True)
				self.fit.fit_opt['std'][i][j]['raw'] = self.transfer_logic(self.__dict__['note_%s'%i]['estd'][j],True)
				self.fit.fit_opt['weight'][i][j]     = self.transfer_float(self.__dict__['note_%s'%i]['weight'][j],True)
				self.fit.fit_opt['psi_end'][i][j]    = self.transfer_float(self.__dict__['note_%s'%i]['range'][j],True)
				count = count + 1
	
		for flag in self.fit.prof_list:
			self.fit.fit_opt['func_type'][flag]   = self.transfer_func(self.note_fn['func_type'][flag],True)

		self.fit.fit_opt['amain'] = self.transfer_float(self.note_in['amain'],True)
		self.fit.fit_opt['zeff']  = self.transfer_float(self.note_in['zeff'],True)
		self.fit.fit_opt['zimp']  = self.transfer_float(self.note_in['zimp'],True)
		self.fit.fit_opt['aimp']  = self.transfer_float(self.note_in['aimp'],True)

		for i in self.fit.inter_list:
			self.fit.fit_opt[i]['val'] = self.transfer_float(self.note_in[i],True)
			self.fit.fit_opt[i]['sig'] = self.transfer_float(self.note_in[i+'s'],True)

		for flag in self.fit.prof_list:
			for k in range(10):
				self.fit.param[flag]['min'][k]    = self.transfer_param(self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['min'],True)
				self.fit.param[flag]['max'][k]    = self.transfer_param(self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['max'],True)
				self.fit.param[flag]['val'][k]    = self.transfer_param(self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['val'],True)
				self.fit.param[flag]['vary'][k]   = self.transfer_logic(self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['fix'],True,True)

		for flag in self.fit.prof_list:
			for k in ['xmin','xmax','ymin','ymax']:
				self.fit.fit_opt['plot'][flag][k] = self.transfer_float(self.__dict__['note_pl'][flag][k],True)

			for kk in self.fit.__dict__['%s_list'%flag]:
				self.fit.fit_opt['scale'][flag][kk]['core'] = self.transfer_float(self.note_ec[flag]['scale'][kk]['core'],True)
				self.fit.fit_opt['scale'][flag][kk]['edge'] = self.transfer_float(self.note_ec[flag]['scale'][kk]['edge'],True)

		self.fit.fit_opt['scale']['ne']['ts']['minc'] = self.transfer_float(self.note_ec['tscale']['minc'],True)
		self.fit.fit_opt['scale']['ne']['ts']['maxc'] = self.transfer_float(self.note_ec['tscale']['maxc'],True)
		self.fit.fit_opt['scale']['ne']['ts']['mine'] = self.transfer_float(self.note_ec['tscale']['mine'],True)
		self.fit.fit_opt['scale']['ne']['ts']['maxe'] = self.transfer_float(self.note_ec['tscale']['maxe'],True)
		self.fit.fit_opt['scale']['ne']['ts']['cn']   = self.transfer_int(self.note_ec['tscale']['cn'],True)
		self.fit.fit_opt['scale']['ne']['ts']['en']   = self.transfer_int(self.note_ec['tscale']['en'],True)		
		self.fit.fit_opt['ascale']   				  = self.transfer_logic(self.note_ec['tscale']['is'],True)
		self.fit.fit_opt['ascale1d']   				  = self.transfer_logic(self.note_ec['tscale']['1d'],True)

		self.fit.fit_opt['mds']['shot']   = self.transfer_int(self.note_md['shot'],True)
		self.fit.fit_opt['mds']['time']   = self.transfer_int(self.note_md['time'],True)
		self.fit.fit_opt['mds']['times']  = self.note_md['times'].get()
		self.fit.fit_opt['mds']['cmse']   = self.transfer_logic(self.note_md['cmse'],True)
		self.fit.fit_opt['mds']['dt']['ts']   = self.transfer_int(self.note_md['dt']['ts'],True)
		self.fit.fit_opt['mds']['dt']['ces']  = self.transfer_int(self.note_md['dt']['ces'],True)
		self.fit.fit_opt['mds']['dt']['tci']  = self.transfer_int(self.note_md['dt']['tci'],True)
		self.fit.fit_opt['mds']['dt']['ece']  = self.transfer_int(self.note_md['dt']['ece'],True)
		self.fit.fit_opt['mds']['dt']['ref']  = self.transfer_int(self.note_md['dt']['ref'],True)			
		self.fit.fit_opt['mds']['dt']['ref2'] = self.transfer_int(self.note_md['dt']['ref2'],True)
		self.fit.fit_opt['mds']['efit']   = self.transfer_int(self.note_md['efit_type'],True)

		self.fit.fit_opt['dacrit']            = self.transfer_float(self.note_md['dacrit'],True)
		self.fit.fit_opt['duty']              = self.transfer_float(self.note_md['duty'],True)

		self.fit.fit_opt['shift']['ti']['ces']	= self.transfer_float(self.note_ec['shift']['ces'],True)
		self.fit.fit_opt['shift']['vt']['ces']	= self.transfer_float(self.note_ec['shift']['ces'],True)
		self.fit.fit_opt['shift']['ne']['ts']	= self.transfer_float(self.note_ec['shift']['ts'],True)
		self.fit.fit_opt['shift']['te']['ts']	= self.transfer_float(self.note_ec['shift']['ts'],True)
		self.fit.fit_opt['shift']['te']['ece']	= self.transfer_float(self.note_ec['shift']['ece'],True)
		self.fit.fit_opt['shift']['ne']['refl']	= self.transfer_float(self.note_ec['shift']['ref'],True)

		self.fit.fit_opt['ashift']['is']        = self.transfer_logic(self.note_ec['ashift']['is'],True)
		self.fit.fit_opt['ashift']['min'] 	= self.transfer_float(self.note_ec['ashift']['min'],True)
		self.fit.fit_opt['ashift']['max'] 	= self.transfer_float(self.note_ec['ashift']['max'],True)

		return

	def get_fit_opt(self):

		self.ped_scan_fit.set(self.transfer_logic(self.fit.fit_opt['ped_scan_fit'],False))
		self.ti_ped_scan.set(self.transfer_logic(self.fit.fit_opt['use_ti_eped'],False))
		self.ti_ped_width.set(self.transfer_logic(self.fit.fit_opt['use_ti_width'],False))

		self.note_in['file']['gfile'].set(self.transfer_name(self.fit.fit_opt['file']['gfile'],False))
		self.note_in['zeff'].set(self.transfer_float(self.fit.fit_opt['zeff'],False))
		self.note_in['zimp'].set(self.transfer_float(self.fit.fit_opt['zimp'],False))
		self.note_in['amain'].set(self.transfer_float(self.fit.fit_opt['amain'],False))
		self.note_in['aimp'].set(self.transfer_float(self.fit.fit_opt['aimp'],False))

		for k in self.fit.inter_list:
			self.note_in[k.lower()].set(self.transfer_float(self.fit.fit_opt[k]['val'],False))
			self.note_in[k.lower()+'s'].set(self.transfer_float(self.fit.fit_opt[k]['sig'],False))
			if self.fit.post['den_diff'][k] > 0.:
				self.note_in[k.lower()+'d'].set(self.transfer_float(self.fit.post['den_diff'][k],False))
			else:
				self.note_in[k.lower()+'d'].set('0.0')
			self.note_in[k.lower()+'f'].set(0)

		if self.fit.fit_opt['line_type'] > 0:
			flag = self.fit.inter_list[self.fit.fit_opt['line_type']-1]
			self.note_in[flag.lower()+'f'].set(1)

		for flag in self.fit.prof_list:
			self.note_fn['raw_fit'][flag].set(self.transfer_logic(self.fit.fit_opt['raw_fit'][flag],False))
			self.note_fn['sep_fix'][flag].set(self.transfer_logic(self.fit.fit_opt['sep_fix'][flag],False))
			self.note_fn['sep_val'][flag].set(self.transfer_float(self.fit.fit_opt['sep_val'][flag],False))
			self.note_fn['width_fix'][flag].set(self.transfer_logic(self.fit.fit_opt['width_fix'][flag],False))
			self.note_fn['width_val'][flag].set(self.transfer_float(self.fit.fit_opt['width_val'][flag],False))
			self.note_fn['func_type'][flag].set(self.transfer_func(self.fit.fit_opt['func_type'][flag],False))

			self.note_in['file']['kfile'][flag].set(self.transfer_name(self.fit.fit_opt['file'][flag]['kfile'],False))
			self.note_fn['use_rho'][flag].set(self.transfer_logic(self.fit.fit_opt['use_rho'][flag],False))
			self.note_fn['smooth'][flag].set(self.transfer_int(self.fit.fit_opt['sspline_order'][flag],False))			

			for k in self.fit.__dict__['%s_list'%flag]:
				self.note_in['file'][flag][k].set(self.transfer_name(self.fit.fit_opt['file'][flag][k],False))

		for flag in self.fit.prof_list:
			for k in range(10):
				self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['min'].set(self.transfer_float(self.fit.param[flag]['min'][k],False))
				self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['max'].set(self.transfer_float(self.fit.param[flag]['max'][k],False))
				self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['val'].set(self.transfer_float(self.fit.param[flag]['val'][k],False))
				self.__dict__['note_%s'%flag]['param']['a%i'%(k+1)]['fix'].set(self.transfer_logic(not self.fit.param[flag]['vary'][k],False))

			for k in self.fit.__dict__['%s_list'%flag]:
				self.__dict__['note_%s'%flag]['avg'][k].set(self.transfer_logic(self.fit.fit_opt['avg'][flag][k],False))
				self.__dict__['note_%s'%flag]['oli'][k].set(self.transfer_logic(self.fit.fit_opt['oli'][flag][k]['use'],False))
				self.__dict__['note_%s'%flag]['std'][k].set(self.transfer_logic(self.fit.fit_opt['std'][flag][k]['use'],False))
				self.__dict__['note_%s'%flag]['estd'][k].set(self.transfer_logic(self.fit.fit_opt['std'][flag][k]['raw'],False))			

				self.__dict__['note_%s'%flag]['range'][k].set(self.transfer_param(self.fit.fit_opt['psi_end'][flag][k],False))
				self.__dict__['note_%s'%flag]['olc'][k].set(self.transfer_param(self.fit.fit_opt['oli'][flag][k]['per'],False))
				self.__dict__['note_%s'%flag]['olcn'][k].set(self.transfer_int(self.fit.fit_opt['oli'][flag]['n'],False))
				self.__dict__['note_%s'%flag]['weight'][k].set(self.transfer_float(self.fit.fit_opt['weight'][flag][k],False))

		for flag in self.fit.prof_list:
			for k in ['xmin','xmax','ymin','ymax']:
				self.note_pl[flag][k].set(self.transfer_float(self.fit.fit_opt['plot'][flag][k],False))

			for kk in self.fit.__dict__['%s_list'%flag]:
				for k in ['core','edge']:
					self.note_ec[flag]['scale'][kk][k].set(self.transfer_float(self.fit.fit_opt['scale'][flag][kk][k],False))
				

		self.note_md['shot'].set(self.transfer_int(self.fit.fit_opt['mds']['shot'],False))
		self.note_md['time'].set(self.transfer_int(self.fit.fit_opt['mds']['time'],False))
		self.note_md['times'].set(self.fit.fit_opt['mds']['times'])
		self.note_md['cmse'].set(self.transfer_logic(self.fit.fit_opt['mds']['cmse'],False))
		self.note_md['dt']['ts'].set(self.transfer_int(self.fit.fit_opt['mds']['dt']['ts'],False))
		self.note_md['dt']['ces'].set(self.transfer_int(self.fit.fit_opt['mds']['dt']['ces'],False))
		self.note_md['dt']['tci'].set(self.transfer_int(self.fit.fit_opt['mds']['dt']['tci'],False))
		self.note_md['dt']['ece'].set(self.transfer_int(self.fit.fit_opt['mds']['dt']['ece'],False))
		self.note_md['dt']['ref'].set(self.transfer_int(self.fit.fit_opt['mds']['dt']['ref'],False))
		self.note_md['efit_type'].set(self.transfer_int(self.fit.fit_opt['mds']['efit'],False))
		try:self.note_md['dt']['ref2'].set(self.transfer_int(self.fit.fit_opt['mds']['dt']['ref2'],False))
		except: self.note_md['dt']['ref2'].set('0')

		self.note_md['dacrit'].set(self.transfer_float(self.fit.fit_opt['dacrit'],False))
		self.note_md['duty'].set(self.transfer_float(self.fit.fit_opt['duty'],False))

		self.note_ec['ashift']['is'].set(self.transfer_logic(self.fit.fit_opt['ashift']['is'],False))
		self.note_ec['ashift']['min'].set(self.transfer_param(self.fit.fit_opt['ashift']['min'],False))
		self.note_ec['ashift']['max'].set(self.transfer_param(self.fit.fit_opt['ashift']['max'],False))

		shift_ces = 0.5* (self.fit.fit_opt['shift']['ti']['ces']+self.fit.fit_opt['shift']['vt']['ces'])
		shift_ts  = 0.5* (self.fit.fit_opt['shift']['te']['ts']+self.fit.fit_opt['shift']['ne']['ts'])
		shift_ref = self.fit.fit_opt['shift']['ne']['refl']
		shift_ece = self.fit.fit_opt['shift']['te']['ece']

		self.note_ec['shift']['ces'].set(self.transfer_param(shift_ces,False))
		self.note_ec['shift']['ts'].set(self.transfer_param(shift_ts,False))
		self.note_ec['shift']['ece'].set(self.transfer_param(shift_ece,False))
		self.note_ec['shift']['ref'].set(self.transfer_param(shift_ref,False))	

		self.note_ec['tscale']['minc'].set(self.transfer_param(self.fit.fit_opt['scale']['ne']['ts']['minc'],False))
		self.note_ec['tscale']['maxc'].set(self.transfer_param(self.fit.fit_opt['scale']['ne']['ts']['maxc'],False))
		self.note_ec['tscale']['mine'].set(self.transfer_param(self.fit.fit_opt['scale']['ne']['ts']['mine'],False))
		self.note_ec['tscale']['maxe'].set(self.transfer_param(self.fit.fit_opt['scale']['ne']['ts']['maxe'],False))
		self.note_ec['tscale']['cn'].set(self.transfer_int(self.fit.fit_opt['scale']['ne']['ts']['cn'],False))
		self.note_ec['tscale']['en'].set(self.transfer_int(self.fit.fit_opt['scale']['ne']['ts']['en'],False))		
		
		self.note_ec['tscale']['is'].set(self.transfer_logic(self.fit.fit_opt['ascale'],False))
		self.note_ec['tscale']['1d'].set(self.transfer_logic(self.fit.fit_opt['ascale1d'],False))

		return

	def sync_gui_opt(self,g2v=True,skip_mds=False):

		# INPUT PAGE
		if g2v: self.note_in['file']['gfile'].set(self.note_in['e1'].get())
		else: self.put_vars_to_gui(['note_in','e1'],self.note_in['file']['gfile'].get())

		count2 = 1;
		for flag in self.fit.prof_list:
			for k in self.fit.__dict__['%s_list'%flag]:				
				count2 = count2 + 1;		
				if g2v: self.note_in['file'][flag][k].set(self.note_in['e%i'%count2].get())
				else: self.put_vars_to_gui(['note_in','e%i'%count2],self.note_in['file'][flag][k].get())

		len2 = len(self.fit.inter_list)
		for i in range(0,len2):
			count2 = count2 + 1;
			flag = self.fit.inter_list[i]
			if g2v: self.note_in[flag].set(self.note_in['e%i'%count2].get())
			else: self.put_vars_to_gui(['note_in','e%i'%count2],self.note_in[flag].get())
			if g2v: self.note_in[flag.lower()+'s'].set(self.note_ne['e%i'%(i+100)].get())
			else: self.put_vars_to_gui(['note_ne','e%i'%(i+100)],self.note_in[flag.lower()+'s'].get())

		if g2v:
			self.note_in['zeff'].set( self.note_in['e%i'%(count2+1)].get())
			self.note_in['zimp'].set( self.note_in['e%i'%(count2+2)].get())
			self.note_in['amain'].set(self.note_in['e%i'%(count2+3)].get())
			self.note_in['aimp'].set( self.note_in['e%i'%(count2+4)].get())
		else:
			self.put_vars_to_gui(['note_in','e%i'%(count2+1)],self.note_in['zeff'].get())
			self.put_vars_to_gui(['note_in','e%i'%(count2+2)],self.note_in['zimp'].get())
			self.put_vars_to_gui(['note_in','e%i'%(count2+3)],self.note_in['amain'].get())
			self.put_vars_to_gui(['note_in','e%i'%(count2+4)],self.note_in['aimp'].get())

		count2 += 5
		for i in range(4):
			flag = self.fit.prof_list[i]
			if g2v: self.note_in['file']['kfile'][flag].set(self.note_in['e%i'%(i+count2)].get())
			else:self.put_vars_to_gui(['note_in','e%i'%(i+count2)],self.note_in['file']['kfile'][flag].get())

		for k in range(4):
			flag = self.fit.prof_list[k]
			if g2v:
				self.note_fn['sep_val'][flag].set(self.note_fn['e%i'%(k+1)].get())
				self.note_fn['width_val'][flag].set(self.note_fn['e%i'%(k+5)].get())
				self.note_fn['smooth'][flag].set(self.note_fn['e%i'%(k+9)].get())
			else:
				self.put_vars_to_gui(['note_fn','e%i'%(k+1)], self.note_fn['sep_val'][flag].get())
				self.put_vars_to_gui(['note_fn','e%i'%(k+5)], self.note_fn['width_val'][flag].get())
				self.put_vars_to_gui(['note_fn','e%i'%(k+9)], self.note_fn['smooth'][flag].get())

		# OTHER PAGES
		for flag in self.fit.prof_list:
			count1 = 1;
			for i in range(10):
				if g2v:
					self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['val'].set(self.__dict__['note_%s'%flag]['e%i'%count1].get())
					self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['min'].set(self.__dict__['note_%s'%flag]['e%i'%(count1+1)].get())
					self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['max'].set(self.__dict__['note_%s'%flag]['e%i'%(count1+2)].get())
				else:
					self.put_vars_to_gui(['note_%s'%flag,'e%i'%(count1+0)],self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['val'].get())
					self.put_vars_to_gui(['note_%s'%flag,'e%i'%(count1+1)],self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['min'].get())
					self.put_vars_to_gui(['note_%s'%flag,'e%i'%(count1+2)],self.__dict__['note_%s'%flag]['param']['a%i'%(i+1)]['max'].get())

				count1 = count1 + 3

			titles   = ['RANGE','OLCN','OLC','WEIGHT']
			var_name = ['RANGE','OLCN','OLC','WEIGHT']
			len2 = len(self.fit.__dict__['%s_list'%flag])
			len3 = len(titles)
			for k in range(len2):
				for j in range(len3):
					kk = self.fit.__dict__['%s_list'%flag][k]
					if g2v: self.__dict__['note_%s'%flag][var_name[j].lower()][kk].set(self.__dict__['note_%s'%flag]['e%i'%(len3*k+j+31)].get())
					else: self.put_vars_to_gui(['note_%s'%flag,'e%i'%(len3*k+j+31)],self.__dict__['note_%s'%flag][var_name[j].lower()][kk].get())

		count1 = 1;
		for flag in self.fit.prof_list:
			for kk in ['xmin','xmax','ymin','ymax']:
				if g2v: self.note_pl[flag][kk].set(self.note_pl['e%i'%count1].get())
				else: self.put_vars_to_gui(['note_pl','e%i'%count1],self.note_pl[flag][kk].get())

				count1 = count1 + 1

		count1 = 1
		for flag in ['te','ne']:
			for k in ['core','edge']:
				if g2v: self.note_ec[flag]['scale']['ts'][k].set(self.note_ec['e%i'%count1].get())
				else: self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec[flag]['scale']['ts'][k].get())
				count1 = count1+1

		for flag in ['ti','vt']:
			for k in ['core','edge']:
				if g2v: self.note_ec[flag]['scale']['ces'][k].set(self.note_ec['e%i'%count1].get())
				else: self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec[flag]['scale']['ces'][k].get())
				count1 = count1+1

		if g2v: 
			self.note_ec['te']['scale']['ece']['core'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1
			self.note_ec['ne']['scale']['refl']['core'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1
			self.note_ec['tscale']['minc'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1		
			self.note_ec['tscale']['maxc'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1
			self.note_ec['tscale']['mine'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1
			self.note_ec['tscale']['maxe'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1
			self.note_ec['tscale']['cn'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1
			self.note_ec['tscale']['en'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1				
			self.note_ec['shift']['ts'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1
			self.note_ec['shift']['ces'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1		
			self.note_ec['shift']['ece'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1
			self.note_ec['shift']['ref'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1						
			self.note_ec['ashift']['min'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1
			self.note_ec['ashift']['max'].set(self.note_ec['e%i'%count1].get())
			count1 = count1+1											
		else:
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['te']['scale']['ece']['core'].get())
			count1 = count1+1
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['ne']['scale']['refl']['core'].get())
			count1 = count1+1
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['tscale']['minc'].get())
			count1 = count1+1		
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['tscale']['maxc'].get())
			count1 = count1+1	
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['tscale']['mine'].get())
			count1 = count1+1	
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['tscale']['maxe'].get())
			count1 = count1+1		
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['tscale']['cn'].get())
			count1 = count1+1	
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['tscale']['en'].get())
			count1 = count1+1									
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['shift']['ts'].get())
			count1 = count1+1
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['shift']['ces'].get())
			count1 = count1+1
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['shift']['ece'].get())
			count1 = count1+1
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['shift']['ref'].get())
			count1 = count1+1			
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['ashift']['min'].get())
			count1 = count1+1
			self.put_vars_to_gui(['note_ec','e%i'%count1],self.note_ec['ashift']['max'].get())
			count1 = count1+1		
	
		if skip_mds: return

		if g2v:
			self.note_md['shot'].set(self.note_md['e1'].get())
			self.note_md['time'].set(self.note_md['e2'].get())
			self.note_md['dt']['ts'].set(self.note_md['e3'].get())
			self.note_md['dt']['ces'].set(self.note_md['e4'].get())
			self.note_md['dt']['tci'].set(self.note_md['e5'].get())
			self.note_md['dt']['ece'].set(self.note_md['e6'].get())
			self.note_md['dt']['ref'].set(self.note_md['e7'].get())
			self.note_md['dt']['ref2'].set(self.note_md['e9'].get())
			self.note_md['times'].set(self.note_md['e8'].get())
			self.note_md['dacrit'].set(self.note_md['e10'].get())
			self.note_md['duty'].set(self.note_md['e11'].get())
		else:
			self.put_vars_to_gui(['note_md','e1'],self.note_md['shot'].get())
			self.put_vars_to_gui(['note_md','e2'],self.note_md['time'].get())
			self.put_vars_to_gui(['note_md','e3'],self.note_md['dt']['ts'].get())
			self.put_vars_to_gui(['note_md','e4'],self.note_md['dt']['ces'].get())
			self.put_vars_to_gui(['note_md','e5'],self.note_md['dt']['tci'].get())
			self.put_vars_to_gui(['note_md','e6'],self.note_md['dt']['ece'].get())
			self.put_vars_to_gui(['note_md','e7'],self.note_md['dt']['ref'].get())
			self.put_vars_to_gui(['note_md','e9'],self.note_md['dt']['ref2'].get())
			self.put_vars_to_gui(['note_md','e8'],self.note_md['times'].get())
			self.put_vars_to_gui(['note_md','e10'],self.note_md['dacrit'].get())
			self.put_vars_to_gui(['note_md','e11'],self.note_md['duty'].get())

		return

	def put_vars_to_gui(self,var,inputv):
	
		if len(var) == 2:
			self.__dict__['%s'%var[0]]['%s'%var[1]].configure(state='normal')
			self.__dict__['%s'%var[0]]['%s'%var[1]].delete(0,'end')
			self.__dict__['%s'%var[0]]['%s'%var[1]].insert(10,inputv)
		elif len(var) == 3:
			self.__dict__['%s'%var[0]]['%s'%var[1]]['%s'%var[2]].configure(state='normal')
			self.__dict__['%s'%var[0]]['%s'%var[1]]['%s'%var[2]].delete(0,'end')
			self.__dict__['%s'%var[0]]['%s'%var[1]]['%s'%var[2]].insert(10,inputv)

		return

	def open_page(self):
		if self.curr_page == 3:
			for i in range(0,7):
				flag = self.fit.inter_list[i]
				vals = self.note_in[flag].get()
				self.note_ne['e%i'%(i+107)].configure(state='normal')
				self.note_ne['e%i'%(i+107)].delete(0,'end')
				self.note_ne['e%i'%(i+107)].insert(10,vals)
				self.note_ne['e%i'%(i+107)].configure(state='readonly')

		if self.curr_page > -1:
			ped_fix = dict()
			for flag in self.fit.prof_list:
				ped_fix[flag] =dict(); ped_fix[flag]['fix'] = False; ped_fix[flag]['fixe'] = False; ped_fix[flag]['ind'] = 0
				if not self.note_fn['func_type'][flag].get().lower() == self.fit.func['list'][0].lower(): 
					if not ((self.note_fn['func_type'][flag].get().lower() == self.fit.func['list'][4].lower()) or (self.note_fn['func_type'][flag].get().lower() == self.fit.func['list'][6].lower())): 
						ped_fix[flag]['ind'] = 3
						if self.note_fn['func_type'][flag].get().lower() == self.fit.func['list'][2].lower(): ped_fix[flag]['ind'] = 6
						if self.__dict__['note_%s'%flag]['param']['a%i'%(ped_fix[flag]['ind'])]['fix'].get() == 1: ped_fix[flag]['fix'] = True
						if self.note_fn['width_fix'][flag].get() == 1: ped_fix[flag]['fixe'] = True; ped_fix[flag]['fix'] = True;

		if ((self.curr_page == 1 and self.prev_page > 0)):		
			if not (ped_fix['te']['fixe'] and ped_fix['ne']['fixe']): 
				if self.ti_ped_width.get() == 1: self.note_fn['c3'].invoke()

			for i in range(4):
				flag = self.fit.prof_list[i]
				func_type = self.note_fn['func_type'][flag].get().lower()
				if not (func_type == 'core' or func_type == 'spline' or func_type == 'nspline'):
					if ped_fix[flag]['fixe']: 
						self.note_fn['width_fix'][flag].set(1)
						self.note_fn['e%i'%(i+5)].delete(0,'end')
						self.note_fn['e%i'%(i+5)].insert(10,self.__dict__['note_%s'%flag]['e%i'%(3*ped_fix[flag]['ind']-2)].get())
					else: self.note_fn['width_fix'][flag].set(0)

				if not (func_type == 'spline' or func_type == 'nspline'):
					if self.__dict__['note_%s'%flag]['param']['a1']['fix'].get() == 1:
						self.note_fn['sep_fix'][flag].set(1)
						self.note_fn['e%i'%(i+1)].delete(0,'end')
						self.note_fn['e%i'%(i+1)].insert(10,self.__dict__['note_%s'%flag]['e1'].get())
					else: self.note_fn['sep_fix'][flag].set(0)

		if ((self.curr_page> 1 and self.curr_page <6) or self.initialise_opt):
			
			if not self.initialise_opt: 
				tlist = [self.fit.prof_list[self.curr_page-2]]				
			else:	
				tlist = self.fit.prof_list
				self.initialise_opt = False
			
			for flag in tlist:
				func_type = self.note_fn['func_type'][flag].get().lower()		
				if not (func_type == 'core' or func_type == 'spline' or func_type == 'nspline'):
					self.__dict__['note_%s'%flag]['c%i'%(ped_fix[flag]['ind'])].configure(state='normal')
					self.__dict__['note_%s'%flag]['e%i'%(3*ped_fix[flag]['ind']-2)].configure(state='normal')
					self.__dict__['note_%s'%flag]['e%i'%(3*ped_fix[flag]['ind']-1)].configure(state='normal')
					self.__dict__['note_%s'%flag]['e%i'%(3*ped_fix[flag]['ind']+0)].configure(state='normal')
					if ped_fix[flag]['fixe'] == 1:
						self.__dict__['note_%s'%flag]['param']['a%i'%(ped_fix[flag]['ind'])]['fix'].set(1)
						self.__dict__['note_%s'%flag]['e%i'%(3*ped_fix[flag]['ind']-2)].delete(0,'end')
						self.__dict__['note_%s'%flag]['e%i'%(3*ped_fix[flag]['ind']-2)].insert(10,self.note_fn['width_val'][flag].get())

						self.__dict__['note_%s'%flag]['c%i'%(ped_fix[flag]['ind'])].configure(state='disabled')
						self.__dict__['note_%s'%flag]['e%i'%(3*ped_fix[flag]['ind']-2)].configure(state='readonly')
						self.__dict__['note_%s'%flag]['e%i'%(3*ped_fix[flag]['ind']-1)].configure(state='readonly')
						self.__dict__['note_%s'%flag]['e%i'%(3*ped_fix[flag]['ind']+0)].configure(state='readonly')					

					elif (ped_fix[flag]['fix'] == 1): self.__dict__['note_%s'%flag]['param']['a%i'%(ped_fix[flag]['ind'])]['fix'].set(1)
					else: self.__dict__['note_%s'%flag]['param']['a%i'%(ped_fix[flag]['ind'])]['fix'].set(0)

				if not (func_type == 'spline' or func_type == 'nspline'):
					self.__dict__['note_%s'%flag]['c1'].configure(state='normal')
					self.__dict__['note_%s'%flag]['e1'].configure(state='normal')
					self.__dict__['note_%s'%flag]['e2'].configure(state='normal')
					self.__dict__['note_%s'%flag]['e3'].configure(state='normal')
					if self.note_fn['sep_fix'][flag].get() == 1:
						self.__dict__['note_%s'%flag]['param']['a1']['fix'].set(1)
						self.__dict__['note_%s'%flag]['e1'].delete(0,'end')
						self.__dict__['note_%s'%flag]['e1'].insert(10,self.note_fn['sep_val'][flag].get())

						self.__dict__['note_%s'%flag]['c1'].configure(state='disabled')
						self.__dict__['note_%s'%flag]['e1'].configure(state='readonly')
						self.__dict__['note_%s'%flag]['e2'].configure(state='readonly')
						self.__dict__['note_%s'%flag]['e3'].configure(state='readonly')

					else: self.__dict__['note_%s'%flag]['param']['a1']['fix'].set(0)

		if (self.curr_page == 3):
			for i in range(0,7):
				vals = self.note_in['e%i'%(i+10)].get()
				self.note_ne['e%i'%(i+107)].configure(state='normal')
				self.note_ne['e%i'%(i+107)].delete(0,'end')
				self.note_ne['e%i'%(i+107)].insert(10,vals)
				self.note_ne['e%i'%(i+107)].configure(state='readonly')
				self.note_ne['e%i'%(i+114)].configure(state='normal')
				self.note_ne['e%i'%(i+114)].delete(0,'end')
				if float(self.note_ne['e%i'%(i+107)].get()) == 0.:
					self.note_ne['e%i'%(i+114)].insert(10,'-N/A-')
				else:	
					self.note_ne['e%i'%(i+114)].insert(10,'%5.3f'%self.fit.post['den_diff'][self.fit.inter_list[i]])
					if self.fit.post['den_diff'][self.fit.inter_list[i]] > 0.: self.note_ne['e%i'%(i+114)].configure(fg='red')
					else:self.note_ne['e%i'%(i+114)].configure(fg='blue')
				self.note_ne['e%i'%(i+114)].configure(state='readonly')

		if (self.curr_page == 6): self.make_mds_e2()

		self.prev_page = self.curr_page	
		return
	
	def run_fit(self):
		self.run_fit_page()
		if self.runmode == 'multi':
			if not (self.fit.fit_opt['mds']['times'] == self.note_md['e8'].get()):
				print('>>> Time slice list is changed. SYNC FIT first!')
				return

		self.sync_gui_opt()
		self.put_fit_opt()

		self.fit.main_run()
		self.fit.post['time'] = self.fit.fit_opt['mds']['time']
		self.draw_plot()

		for flag in self.fit.prof_list:
			func = self.note_fn['func_type'][flag].get().lower()
			if func in self.post_opt['param'][flag].keys():
				self.post_opt['param'][flag][func]['is'] = True
				self.post_opt['param'][flag][func]['val'] = copy.deepcopy(self.fit.param[flag])
			else:
				self.post_opt['param'][flag][func] = dict()
				self.post_opt['param'][flag][func]['is'] = True
				self.post_opt['param'][flag][func]['val'] = copy.deepcopy(self.fit.param[flag])

		if self.runmode=='multi': ttime = self.mfit_opt['times'][self.plotpage]
		else: ttime = self.fit.post['time']

		self.post_opt['param_m'][ttime] = copy.deepcopy(self.post_opt['param'])

		self.didfit = True
		self.run_post()
		if self.runmode == 'multi':
			self.mfit_opt[self.mfit_opt['times'][self.plotpage]] = copy.deepcopy([self.fit.fit_opt,self.fit.te_prof,self.fit.ne_prof,self.fit.ti_prof,self.fit.vt_prof,self.fit.post,self.fit.fit_eq])

			nel = np.zeros(7); bpped = np.zeros(2); width = np.zeros(8); ww = np.zeros(2); tci = np.zeros(7); tcie = np.copy(tci)
			nel = [self.fit.post['int01'],self.fit.post['int02'],self.fit.post['tci01'],self.fit.post['tci02'],self.fit.post['tci03'],self.fit.post['tci04'],self.fit.post['tci05']]
			bpped = [self.fit.post['bpped'],min(self.fit.post['bppede'],self.fit.post['bpped'])]
			for k in range(4):
				width[k]   = self.fit.post['width'][k]
				width[k+4] = min(self.fit.post['width'][k],self.fit.post['widthe'][k])
			ww = [self.fit.post['wmhd'],self.fit.post['wkin']]
			for kk in range(7):
				tci[kk] = self.fit.fit_opt[self.fit.inter_list[kk]]['val']
				tcie[kk] = self.fit.fit_opt[self.fit.inter_list[kk]]['sig']			
			self.mfit_opt['post'][self.mfit_opt['times'][self.plotpage]] = copy.deepcopy([nel,bpped,width,ww,tci,tcie])
			self.fit.write_kinprof()
			self.fit.read_kinprof_kprofile()
			try: rmtree('MFIT/%i'%self.mfit_opt['times'][self.plotpage])
			except: pass
			copytree('PROFILES','MFIT/%i'%self.mfit_opt['times'][self.plotpage])

		return

	def draw_plot(self):
		self.figure['name']['home'].canvas.draw_idle()
		if self.plot_type  == 1:  self.fit.draw_plot(fig_ex=self.figure['name']['home'],second=True)
		elif self.plot_type== 0:  self.fit.draw_plot(fig_ex=self.figure['name']['home'],second=False)
		elif self.plot_type== 2:  
			self.fit.draw_plot(fig_ex=self.figure['name']['home'],second=False)
			self.plot_type = 0
			self.make_plot_c()
		elif self.plot_type== 3: self.make_plot_d()
		elif self.plot_type== 4: self.make_plot_e()

		return

	def run_fit_page(self):

		if self.curr_page==0:
			for i in range(0,7):
				vals = self.note_in['e%i'%(i+10)].get()
				self.note_ne['e%i'%(i+107)].configure(state='normal')
				self.note_ne['e%i'%(i+107)].delete(0,'end')
				self.note_ne['e%i'%(i+107)].insert(10,vals)
				self.note_ne['e%i'%(i+107)].configure(state='readonly')

		if (self.curr_page > 1 and self.curr_page <6):
			curr_page = self.curr_page
			prev_page = self.prev_page
			self.curr_page = 1
			self.prev_page = curr_page
			self.open_page()
			self.curr_page = curr_page
			self.prev_page = prev_page
		if self.curr_page == 1:
			curr_page = self.curr_page
			prev_page = self.prev_page
			self.curr_page = 2
			self.prev_page = curr_page
			self.open_page()
			self.curr_page = curr_page
			self.prev_page = prev_page
		return

	def run_post(self):

		if self.curr_page==3:
			for i in range(0,7):
				vals = self.note_in['e%i'%(i+10)].get()
				self.note_ne['e%i'%(i+107)].configure(state='normal')
				self.note_ne['e%i'%(i+107)].delete(0,'end')
				self.note_ne['e%i'%(i+107)].insert(10,vals)
				self.note_ne['e%i'%(i+107)].configure(state='readonly')
				self.note_ne['e%i'%(i+114)].configure(state='normal')
				self.note_ne['e%i'%(i+114)].delete(0,'end')
				if float(self.note_ne['e%i'%(i+107)].get()) == 0.:
					self.note_ne['e%i'%(i+114)].insert(10,'-N/A-')
				else:	
					self.note_ne['e%i'%(i+114)].insert(10,'%5.3f'%self.fit.post['den_diff'][self.fit.inter_list[i]])
					if self.fit.post['den_diff'][self.fit.inter_list[i]] > 0.: self.note_ne['e%i'%(i+114)].configure(fg='red')
					else:self.note_ne['e%i'%(i+114)].configure(fg='blue')
				self.note_ne['e%i'%(i+114)].configure(state='readonly')

		count1 = 1
		for flag in self.fit.prof_list:
			for kk in ['xmin','xmax','ymin','ymax']:
				self.note_pl['e%i'%count1].delete(0,'end')
				self.note_pl['e%i'%count1].insert(10,'%4.2f'%self.fit.fit_opt['plot'][flag][kk])
			
				count1 = count1 + 1
	
		if self.fit.fit_opt['ascale']:
			self.note_ec['e3'].delete(0,'end')
			self.note_ec['e3'].insert(10,'%5.3f'%self.fit.fit_opt['scale']['ne']['ts']['core'])
			self.note_ec['e4'].delete(0,'end')
			self.note_ec['e4'].insert(10,'%5.3f'%self.fit.fit_opt['scale']['ne']['ts']['edge'])

		if self.fit.fit_opt['ashift']['is']:
			self.note_ec['e18'].delete(0,'end')
			self.note_ec['e18'].insert(10,'%5.3f'%self.fit.fit_opt['shift']['ti']['ces'])

		if self.plot_type ==2: self.plot_type ==0			
		return

	def run_save(self):

		if not (self.didfit or self.didmfit): return

		try: os.mkdir('PROFILES')
		except: pass		

		inputf = asksaveasfilename()
		if (len(inputf) == 0):      return

		self.fit.write_kinprof()

		for i in range(4):
			flag = self.fit.prof_list[i]
			filename = 'PROFILES/%s_fit.dat'%flag.upper()
			self.note_in['e%i'%(i+21)].delete(0,'end')			
			if not os.path.isfile(filename): continue

			self.note_in['e%i'%(i+21)].insert(10,filename)
			self.fit.fit_opt['file'][flag]['kfile'] = filename

		self.fit.read_kinprof_kprofile()
		self.draw_plot()

		f = open(inputf,'wb')
		pickle.dump(self.fit.fit_opt,f)
		f.close()
		f = open(inputf+'_param','wb')
		pickle.dump(self.post_opt['param'],f)
		f.close()
		f = open(inputf+'_param_m','wb')
		pickle.dump(self.post_opt['param_m'],f)
		f.close()

		copyfile(inputf,'fit_opt.save')
		copyfile(inputf+'_param','fit_opt.save_param')
		copyfile(inputf+'_param_m','fit_opt.save_param_m')

		print('>>> Fitting option is saved')

		if self.runmode == 'single':
			f = open('run_mode','w')
			f.write('single')
			f.close()
			print('>>> Run mode = SINGLE')

			f = open('result_single.save','wb')
			pickle.dump([self.fit.fit_opt,self.fit.te_prof,self.fit.ne_prof,self.fit.ti_prof,self.fit.vt_prof,self.fit.post,self.fit.fit_eq],f)
			f.close()
			try: rmtree('SPROFILES')
			except: pass
			copytree('PROFILES','SPROFILES')
			
		else:	
			f = open('run_mode','w')
			f.write('multi %s'%self.mfit_opt['times'])
			f.close()
			try: os.system('rm MFIT/fit*.save')
			except: pass
			self.times = copy.deepcopy(self.mfit_opt['times'])
			for i in self.mfit_opt['times']:
				try: rmtree('MPROFILES/%i'%i)
				except: pass
				copytree('MFIT/%i'%i,'MPROFILES/%i'%i)
				f = open('MFIT/fit_%05i.save'%i,'wb')
				pickle.dump(self.mfit_opt[i][0],f)
				f.close()
				f = open('MFIT/fit_%05i_param.save'%i,'wb')
				pickle.dump(self.mfit_opt[i][0],f)
				f.close()
			try: rmtree('MPROFILES/%i'%self.times[self.plotpage])
			except: pass				
			copytree('PROFILES','MPROFILES/%i'%self.times[self.plotpage])

			f = open('result_multi.save','wb')
			pickle.dump(self.mfit_opt,f)
			f.close()
			print('>>> Run mode = MULTI')

		return

	def read_gfitp(self):

		try:	f=open('gfitp_history.dat','r')
		except:	
			print('>>> NO FITP history')
			return

		while True:
			line = f.readline()
			if not line: break
			line = line.split()

			if line[0].lower()   == 'eq_file': self.fit.fit_opt['file']['gfile']= '../'+line[1]
			elif line[0].lower() == 'te_file': self.fit.fit_opt['file']['te']['ts']= '../'+line[1]
			elif line[0].lower() == 'ne_file': self.fit.fit_opt['file']['ne']['ts']= '../'+line[1]
			elif line[0].lower() == 'te_edge_file': self.fit.fit_opt['file']['te']['tse']= '../'+line[1]
			elif line[0].lower() == 'ne_edge_file': self.fit.fit_opt['file']['ne']['tse']= '../'+line[1]
			elif line[0].lower() == 'ti_file': self.fit.fit_opt['file']['ti']['ces']= '../'+line[1]
			elif line[0].lower() == 'vt_file': self.fit.fit_opt['file']['vt']['ces']= '../'+line[1]
			elif line[0].lower() == 'shot':    self.fit.fit_opt['mds']['shot']= int(float(line[1]))
			elif line[0].lower() == 'time':    self.fit.fit_opt['mds']['time']= int(float(line[1]))
			elif line[0].lower() == 'mds_dir': self.mds_dir = line[1]
			elif line[0].lower() == 'lineden':
				if round(float(line[1]),2) > 0.:
					self.fit.fit_opt['use_density_scale'] = True
					self.fit.fit_opt['target_density']    = float(line[1])
					self.fit.fit_opt['use_density_scale'] = True
					self.fit.fit_opt['line_type']         = 1

			elif line[0].lower() == 'zeff': self.fit.fit_opt['zeff']  = round(float(line[1]),2)
			elif line[0].lower() == 'zimp': self.fit.fit_opt['zimp']  = round(float(line[1]),2)

		f.close()

		return		

	def load_fitopt(self):

		inputf = askopenfilename()
		if (len(inputf) == 0):      return

		f = open(inputf,'rb')
		self.load_fitopt_dict(pickle.load(f),True)
		f.close()			

		if os.path.isfile(inputf+'_param'):
			print('>>> Load FIT params')
			f = open(inputf+'_param','rb')
			temp = pickle.load(f)
			f.close()			
			for flag in self.fit.prof_list:
				for key in temp[flag].keys():
					self.post_opt['param'][flag][key] = copy.deepcopy(temp[flag][key])
		if os.path.isfile(inputf+'_param_m'):
			f = open(inputf+'_param_m','rb')
			self.post_opt['param_m'] = pickle.load(f)
			f.close()
		self.restore_gui()
		return

	def load_fitopt_dict(self,fit_opt,input_skip=False):

		#if not input_skip:
		for i in self.fit.inter_list:
			self.fit.fit_opt[i]['val'] = fit_opt[i]['val']
			self.fit.fit_opt[i]['sig'] = fit_opt[i]['sig']

		self.fit.fit_opt['use_ti_width']      = fit_opt['use_ti_width']
		self.fit.fit_opt['use_ti_eped']       = fit_opt['use_ti_eped'] 
		self.fit.fit_opt['ped_scan_fit']      = fit_opt['ped_scan_fit']
		self.fit.fit_opt['line_type']         = fit_opt['line_type'] 
		self.fit.fit_opt['target_density']    = fit_opt['target_density'] 
		self.fit.fit_opt['use_density_scale'] = fit_opt['use_density_scale']
		
		if not input_skip:
			self.fit.fit_opt['file']['gfile']     = fit_opt['file']['gfile']

		for i in self.fit.prof_list:	
			self.fit.fit_opt['use_rho'][i]       = fit_opt['use_rho'][i]
			self.fit.fit_opt['sep_fix'][i]       = fit_opt['sep_fix'][i]
			self.fit.fit_opt['sep_val'][i]       = fit_opt['sep_val'][i] 
			self.fit.fit_opt['width_fix'][i]     = fit_opt['width_fix'][i]
			self.fit.fit_opt['width_val'][i]     = fit_opt['width_val'][i]
			self.fit.fit_opt['raw_fit'][i]       = fit_opt['raw_fit'][i]
			self.fit.fit_opt['sspline_order'][i] = fit_opt['sspline_order'][i]

			self.fit.fit_opt['oli'][i]['n']  = fit_opt['oli'][i]['n']
			if not input_skip:
				self.fit.fit_opt['file'][i]['kfile'] = fit_opt['file'][i]['kfile']

			for j in self.fit.__dict__['%s_list'%i]:
				if not j in fit_opt['file'][i].keys(): continue
				if not input_skip:
					self.fit.fit_opt['file'][i][j]   = fit_opt['file'][i][j]
				self.fit.fit_opt['avg'][i][j]		 = fit_opt['avg'][i][j]
				self.fit.fit_opt['shift'][i][j]	     = fit_opt['shift'][i][j]	
				self.fit.fit_opt['exclude'][i][j]	 = fit_opt['exclude'][i][j]
				self.fit.fit_opt['oli'][i][j]['use'] = fit_opt['oli'][i][j]['use']
				self.fit.fit_opt['oli'][i][j]['per'] = fit_opt['oli'][i][j]['per']
				self.fit.fit_opt['std'][i][j]['use'] = fit_opt['std'][i][j]['use']
				self.fit.fit_opt['std'][i][j]['raw'] = fit_opt['std'][i][j]['raw']
				self.fit.fit_opt['weight'][i][j]     = fit_opt['weight'][i][j]
				self.fit.fit_opt['psi_end'][i][j]    = fit_opt['psi_end'][i][j]		
		self.fit.fit_opt['func_type']['te']   = fit_opt['func_type']['te']
		self.fit.fit_opt['func_type']['ne']   = fit_opt['func_type']['ne']
		self.fit.fit_opt['func_type']['ti']   = fit_opt['func_type']['ti']
		self.fit.fit_opt['func_type']['vt']   = fit_opt['func_type']['vt']

		#if not input_skip:
		self.fit.fit_opt['amain'] = fit_opt['amain']
		self.fit.fit_opt['zeff']  = fit_opt['zeff']
		self.fit.fit_opt['zimp']  = fit_opt['zimp']
		self.fit.fit_opt['aimp']  = fit_opt['aimp']

		try:
			self.fit.fit_opt['dacrit']= fit_opt['dacrit']
			self.fit.fit_opt['duty']  = fit_opt['duty']
		except: pass

		try:
			for flag in self.fit.prof_list:
				self.fit.fit_opt['plot'][flag]['xmin'] = fit_opt['plot'][flag]['xmin']
				self.fit.fit_opt['plot'][flag]['xmax'] = fit_opt['plot'][flag]['xmax']
				self.fit.fit_opt['plot'][flag]['ymin'] = fit_opt['plot'][flag]['ymin']
				self.fit.fit_opt['plot'][flag]['ymax'] = fit_opt['plot'][flag]['ymax']
		except: pass

		try:
			for flag in self.fit.prof_list:
				for kk in self.fit.__dict__['%s_list'%flag]:
					self.fit.fit_opt['scale'][flag][kk]['core'] = fit_opt['scale'][flag][kk]['core']
					self.fit.fit_opt['scale'][flag][kk]['edge'] = fit_opt['scale'][flag][kk]['edge']
		except: pass

		try:
			self.fit.fit_opt['ashift']['is']  = fit_opt['ashift']['is']
			self.fit.fit_opt['ashift']['min'] = fit_opt['ashift']['min']
			self.fit.fit_opt['ashift']['max'] = fit_opt['ashift']['max']
		except: pass

		try:
			self.fit.fit_opt['mds']['shot']      = fit_opt['mds']['shot']
			self.fit.fit_opt['mds']['time']      = fit_opt['mds']['time']
			self.fit.fit_opt['mds']['times']     = fit_opt['mds']['times']
			self.fit.fit_opt['mds']['dt']['ts']  = fit_opt['mds']['dt']['ts']
			self.fit.fit_opt['mds']['dt']['ces'] = fit_opt['mds']['dt']['ces']
			self.fit.fit_opt['mds']['dt']['tci'] = fit_opt['mds']['dt']['tci']
			self.fit.fit_opt['mds']['dt']['ece'] = fit_opt['mds']['dt']['ece']
			self.fit.fit_opt['mds']['dt']['ref'] = fit_opt['mds']['dt']['ref']
			try: self.fit.fit_opt['mds']['dt']['ref2'] = fit_opt['mds']['dt']['ref2']
			except: pass
			try: self.fit.fit_opt['mds']['efit'] = fit_opt['mds']['efit']
			except: pass
		except: pass

		try:
			self.fit.fit_opt['scale']['ne']['ts']['minc']  = fit_opt['scale']['ne']['ts']['minc']
			self.fit.fit_opt['scale']['ne']['ts']['maxc']  = fit_opt['scale']['ne']['ts']['maxc']
			self.fit.fit_opt['scale']['ne']['ts']['mine']  = fit_opt['scale']['ne']['ts']['mine']
			self.fit.fit_opt['scale']['ne']['ts']['maxe']  = fit_opt['scale']['ne']['ts']['maxe']
			self.fit.fit_opt['scale']['ne']['ts']['cn']    = fit_opt['scale']['ne']['ts']['cn']
			self.fit.fit_opt['scale']['ne']['ts']['en']    = fit_opt['scale']['ne']['ts']['en']
			self.fit.fit_opt['ascale']                     = fit_opt['ascale']
		except: pass		
		try: self.fit.fit_opt['ascale1d'] = fit_opt['ascale1d']	
		except: pass

		return

	def load_fitopt_param(self):

		inputf = askopenfilename()
		if (len(inputf) == 0):      return

		f = open(inputf,'rb')
		temp = pickle.load(f)
		f.close()			
		for flag in self.fit.prof_list:
			for key in temp[flag].keys():
				self.post_opt['param'][flag][key] = copy.deepcopy(temp[flag][key])
		for flag in self.fit.prof_list:
			self.fitting_func_option_a(flag)

		return

	def reset_fitopt(self):

		temp = copy.deepcopy(self.fit.fit_opt['file'])
		self.fit.initialise_fitopt()
		for key in temp.keys(): self.fit.fit_opt['file'][key] = copy.deepcopy(temp[key])
		self.get_fit_opt()		
		self.sync_gui_opt(g2v=False)
		for flag in self.fit.prof_list:
			self.fitting_func_option_a(flag,False)
		self.fit_constraint_option_a()
		return

	def run_exit(self):

		self.root.destroy()
		return

	def read_multi_fitopt(self):

		file1 = 'mult_opt.dat'
		self.ishmode = True
		self.forcefit = False
		self.tavg = [0.,0.,0.,0.]		
		if not os.path.isfile(file1): return

		f = open(file1,'rb')
		force_opt = pickle.load(f)
		f.close()

		self.fit.fit_opt['mds']['shot'] = int(force_opt['shot'])
		self.fit.fit_opt['mds']['times']= force_opt['times']
		self.fit.fit_opt['mds']['time'] = force_opt['time']

		self.fit.fit_opt['mds']['dt']['ts']  = force_opt['tsdt']
		self.fit.fit_opt['mds']['dt']['ces'] = force_opt['cesdt']
		self.fit.fit_opt['mds']['dt']['tci'] = force_opt['tcidt']
		self.fit.fit_opt['mds']['dt']['ece'] = force_opt['ecedt']

		self.fit.fit_opt['exclude']['te']['ts'] = force_opt['ex_tst']
		self.fit.fit_opt['exclude']['ne']['ts'] = force_opt['ex_tsn']
		self.fit.fit_opt['exclude']['ti']['ces'] = force_opt['ex_cest']
		self.fit.fit_opt['exclude']['vt']['ces'] = force_opt['ex_cesv']
		self.fit.fit_opt['exclude']['te']['ece'] = force_opt['ex_ece']
	
		try:
			self.fit.fit_opt['exclude']['te']['tse'] = force_opt['ex_tset']	
			self.fit.fit_opt['exclude']['ne']['tse'] = force_opt['ex_tsen']
			self.fit.fit_opt['duty']  = force_opt['duty']
			self.fit.fit_opt['dacrit']= force_opt['dacrit']
		except: pass

		self.fit.fit_opt['sep_val']['te'] = float(force_opt['tsep'])
		self.fit.fit_opt['sep_val']['ti'] = float(force_opt['tsep'])
		self.fit.force_width = float(force_opt['twidmin'])

		self.fit.fit_opt['scale']['te']['ts']['edge'] = float(force_opt['tstem'])
		self.fit.fit_opt['scale']['ne']['ts']['edge'] = float(force_opt['tsnem'])
		self.fit.fit_opt['scale']['te']['ts']['core'] = float(force_opt['tstcm'])
		self.fit.fit_opt['scale']['ne']['ts']['core'] = float(force_opt['tsncm'])

		self.fit.fit_opt['scale']['ne']['ts']['minc'] = float(force_opt['minc'])
		self.fit.fit_opt['scale']['ne']['ts']['maxc'] = float(force_opt['maxc'])
		self.fit.fit_opt['scale']['ne']['ts']['mine'] = float(force_opt['mine'])
		self.fit.fit_opt['scale']['ne']['ts']['maxe'] = float(force_opt['maxe'])
		self.fit.fit_opt['ascale'] = force_opt['ascale']
		self.fit.fit_opt['ascale1d']=force_opt['ascale1d']
		self.fit.fit_opt['zeff'] = force_opt['zeff']

		self.fit.fit_opt['scale']['ne']['ts']['cn'] = int(force_opt['cn'])
		self.fit.fit_opt['scale']['ne']['ts']['en'] = int(force_opt['en'])

		for kk in self.fit.inter_list:
			self.fit.fit_opt[kk]['sig'] = float(force_opt[kk])

		flag = ['ts','ces','ece','tci']
		for k in range(4):
			self.tavg[k] = force_opt[flag[k]]

		self.fit.fit_opt['ashift']['is'] = force_opt['ashift']
		self.ishmode = force_opt['ishmode']
		self.forcefit= force_opt['forcefit']
	
		return

	def forced_fit(self):

		try: os.mkdir('MPROFILES')
		except: pass

		self.fit = fittool.fit_tool()
		self.fit.forcefit = True	
		self.fit.noprint = True
		self.didfit = True
		self.runmode = 'multi'
		file1 = 'mult_opt.dat'
		if not os.path.isfile(file1): return
		self.read_multi_fitopt()
		self.fit.forced_fit(self.ishmode)
		if not self.forcefit: return False
		
		shot   = self.fit.fit_opt['mds']['shot']
		times  = self.fit.fit_opt['mds']['times'].split(',')
		ts_ta  = self.tavg[0]
		ces_ta = self.tavg[1]
		tci_ta = self.tavg[3]
		ece_ta = self.tavg[2]

		self.filename = dict()
		for flag in self.fit.prof_list: self.filename[flag] = dict()
		if len(times) == 0: 
			print('>>> Times are not given...')
			return
		self.times = np.array(times,dtype='int')
		self.shot  = int(shot)
		efit_list  = get_efit_list2(self.shot)
		efit_no    = 1
		tlen = len(self.times)
		self.files = dict()
		try: os.mkdir('GFILES')
		except: pass
		try: os.mkdir('INPUTS')
		except: pass
		self.filename['gfile'] = dict()
		self.filename['tci']   = dict()
		print('>>> Load gfiles')
		for i in range(tlen):
			filename = 'GFILES/g%06i.%06i_%i'%(self.shot,self.times[i],efit_no)
			ttime = self.times[i]
			tind  = np.argmin(abs(efit_list['times'][efit_no]-ttime))
			ttime2= efit_list['times'][efit_no][tind]
			if not os.path.isfile(filename):
				filename0= efit_list['dirs'][efit_no]+'/g%06i.%06i'%(self.shot,ttime2)
				copyfile(filename0,filename)
			self.filename['gfile'][i] = filename
			self.filename['tci'][i] = None
			update_progress(float((i+1)/tlen))

		chan = ['TS','TSE','CES','TCI','ECE']
		for flag in self.fit.prof_list:
			self.filename[flag]['kfile'] = dict()
			for kk in self.fit.__dict__['%s_list'%flag]:
				self.filename[flag][kk] = dict()
				for ll in range(tlen):
					self.filename[flag]['kfile'][ll] = None
					self.filename[flag][kk][ll] = None	


		tas = [ts_ta,ts_ta,ces_ta,tci_ta,ece_ta]
		dirs = '%s/DATASAVE/%i'%(self.mds_dir,self.shot)
		self.get_elm_peak(dirs,self.times);	
		for k in range(5):
			flag = chan[k].lower()
			if float(tas[k]) == 0: continue
			print('>>> Load %s raw profiles'%flag.upper())
			if (flag == 'ts' or flag =='tse' or flag == 'ces'): 
				dirs = '%s/DATASAVE/%i'%(self.mds_dir,self.shot)
				ftemp,flist = self.make_ts_ces_raw_files(float(tas[k]),flag,dirs,'INPUTS')
				for kk in flist:
					if ftemp[kk]['isfile']:
						print('>>> There is %s - %s raw file'%(flag.upper(),kk.upper()))
						for ll in range(tlen):
							if (kk=='ne' or kk=='te'): 
								if flag == 'ts': self.filename[kk]['ts'][ll] = ftemp[kk][ll]
								else:            self.filename[kk]['tse'][ll]= ftemp[kk][ll]
							else:	self.filename[kk]['ces'][ll] = ftemp[kk][ll]

			if (flag == 'tci'):
				dirs = '%s/DATASAVE/%i'%(self.mds_dir,self.shot)
				ftemp,isfile = self.make_tci_raw_files(float(tas[k]),dirs,'INPUTS')
				if isfile: print('>>> There is %s - %s raw file'%(flag.upper(),'NE'))
				for ll in range(tlen):
					self.filename['tci'][ll] = ftemp[ll]

			if (flag == 'ece'):
				dirs = '%s/DATASAVE/%i'%(self.mds_dir,self.shot)
				ftemp = self.make_ece_raw_files(float(tas[k])/1.e3,dirs)
				if ftemp['isfile']: 
					print('>>> There is ECE raw filw')
					for ll in range(tlen):
						self.filename['te']['ece'][ll] = ftemp[ll]

		self.mfit_opt = dict()
		self.mfit_opt['times'] = copy.deepcopy(self.times)
		self.mfit_opt['post'] = dict()

		for flag in self.fit.prof_list: self.fit.__dict__['%s_prof'%flag]['fit_old'] = []
		tlen = len(self.times)
		nel = np.zeros(7); bpped = np.zeros(2); width = np.zeros(8); ww = np.zeros(2); tci = np.zeros(7); tcie = np.copy(tci)
		for i in range(tlen):
			self.reset_input_file()
			self.fit.fit_opt['file']['gfile'] = self.filename['gfile'][i]
			if not self.filename['tci'][i] == None:
				f = open(self.filename['tci'][i],'r')
				for tflag in self.fit.inter_list:
					line = f.readline()
					self.fit.fit_opt[tflag]['val'] = round(float(line.split()[1]),2)
			for flag in self.fit.prof_list:
				for kk in self.fit.__dict__['%s_list'%flag]:
					self.fit.fit_opt['file'][flag][kk] = self.filename[flag][kk][i]
			print('>>> Fitting for %05i [ms]'%self.times[i])
			for flag in self.fit.prof_list: self.fit.post['didfit'][flag] = False
			self.fit.first_run = True
			self.fit.fit_opt['mds']['time'] = self.times[i]
			self.fit.main_run()
			self.fit.write_kinprof()
			self.fit.read_kinprof_kprofile()
			self.fit.post['time'] = self.times[i]
			nel = [self.fit.post['int01'],self.fit.post['int02'],self.fit.post['tci01'],self.fit.post['tci02'],self.fit.post['tci03'],self.fit.post['tci04'],self.fit.post['tci05']]
			bpped = [self.fit.post['bpped'],min(self.fit.post['bppede'],self.fit.post['bpped'])]
			for k in range(4):
				width[k]   = self.fit.post['width'][k]
				width[k+4] = min(self.fit.post['width'][k],self.fit.post['widthe'][k])
			ww = [self.fit.post['wmhd'],self.fit.post['wkin']]
			for kk in range(7):
				tci[kk] = self.fit.fit_opt[self.fit.inter_list[kk]]['val']
				tcie[kk] = self.fit.fit_opt[self.fit.inter_list[kk]]['sig']
			copyfile('PROFILES/chease_kinprof_fit','MFIT/chease_kinprof_fit_%05i_temp'%self.times[i])
			try: copyfile('PROFILES/VT_fit.dat','MFIT/VT_fit.dat_%05i_temp'%self.times[i])
			except: pass
			self.mfit_opt[self.times[i]] = copy.deepcopy([self.fit.fit_opt,self.fit.te_prof,self.fit.ne_prof,self.fit.ti_prof,self.fit.vt_prof,self.fit.post,self.fit.fit_eq])
			self.mfit_opt['post'][self.times[i]] = copy.deepcopy([nel,bpped,width,ww,tci,tcie])

			try: rmtree('MFIT/%i'%self.times[i])
			except: pass
			copytree('PROFILES','MFIT/%i'%self.times[i])

		try: os.system('rm MFIT/fit*.save')
		except: pass
		self.times = copy.deepcopy(self.mfit_opt['times'])
		for i in self.mfit_opt['times']:
			try: rmtree('MPROFILES/%i'%i)
			except: pass
			move('MFIT/%i'%i,'MPROFILES/%i'%i)
			f = open('MFIT/fit_%05i.save'%i,'wb')
			pickle.dump(self.mfit_opt[i][0],f)
			f.close()
			f = open('MFIT/fit_%05i_param.save'%i,'wb')
			pickle.dump(self.mfit_opt[i][0],f)
			f.close()

		f = open('result_multi.save','wb')
		pickle.dump(self.mfit_opt,f)
		f.close()
		return True

	def __init__(self,single_only=False,efit_only=False,multi_only=False):

		self.lmfit_mod = False
		self.input_ready = False
		self.input_load = False
		self.single_only = single_only
		self.multi_only = multi_only
		self.efit_only = False
		if efit_only:
			self.efit_only = True
			self.single_only = True

		return

if __name__ == "__main__":

	import guifit

	try: os.mkdir('PROFILES')
	except: pass

	print(' ---------------------------------------------------------------')
	print('||               Function based Fitting tool Ver %3s           ||'%version['gfit'])
	print('||                   Kinetic Profile generator                 ||')
	print('%s'%author['gfit'])
	print('%s'%comment['gfit'])
	print(' ---------------------------------------------------------------\n')

	now = time.gmtime(time.time())
	years = now.tm_year, now.tm_mon, now.tm_mday
	hours = now.tm_hour, now.tm_min, now.tm_sec
	line = '#-- %i/%02i/%02i -- %02ih %02im %02is --'%(years[0],years[1],years[2],hours[0],hours[1],hours[2])
	os.system("echo '%s' >> log.gfit"%line)

	try: flag = sys.argv[1]
	except: flag = ''
	efit_only = False
	single_only = False
	multi_only = False
	if flag == '-fitp': single_only = True
	elif flag == '-efit': efit_only = True
	elif flag == '-mult': multi_only = True

	guifit = guifit.guifittool(single_only,efit_only,multi_only)

	if multi_only:
		guifit.mds_dir = './MDS'
		forcefit = guifit.forced_fit()
		if forcefit: exit()
		guifit.fit.forcefit2 = True
		guifit.fit.hmodefit = guifit.ishmode
		
	guifit.root = tk.Tk()
	guifit.root.title('KFIT')
	guifit.gfit()
