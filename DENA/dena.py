#!/usr/local/anaconda3/bin/python3
import os,sys, copy, time, pickle

import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilename,asksaveasfilename

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.widgets import MultiCursor
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

from get_efit import *
from progress import *
from exec_dirs import dummy_dir, rdena_db_dir, version
import eqdsk
import fittool

import rtefit

import numpy as np
from scipy.interpolate import interp1d
from shutil import move, copyfile, copytree, rmtree
from MDS import mds
import pickle

currdir = os.getcwd()

class kstar_density_tool:

	def _den_tool(self):

		if not self.read_mode: self._make_directories()
		self._load_fittool()
		self.root.protocol('WM_DELETE_WINDOW', lambda: self._exit())
		self._declare_variables()
		self._initialise_variables()
		self.note_in['gflag'].set(self.efitv)
		self._load_diag()
		self._shot_preset()
		if not self.nogui: self._make_home_canvas()
		self._make_note_frame()
		for k in range(4): self.root.grid_columnconfigure(k,weight=1)

		self._make_input_frame()
		self._make_post_frame()
		self._update_diagno_state()
		self._update_time_list()

		self._preset_tci_opt()
			
		if self.brutal:
			self._add()
			self._fit()
		else:
			if not self.nogui:
				self._draw_overview()
				self._draw_tci_fit()
			if self.read_mode: 
				self._load_prev()
				self._fit()
				if self.profiles['diag']['ts']: self._auto_scale()
				if self.nogui: 
					cm = self.note_in['CM'].get()
					em = self.note_in['EM'].get()
					print('>>> CM/EM: %s/%s'%(cm,em))
		if not self.nogui:
			self.root.resizable(0,0)
			self.root.mainloop()
			

		return

	def _load_prev(self):
		f = open(self.in_file,'rb')
		sav = pickle.load(f); 
		f.close()
		for key in sav.keys(): self.profiles[key] = sav[key]
		for tt in self.profiles['times']:
			line = '  %06i ms'%int(tt)
			self.note_in['l2'].insert('end',line)
			if not tt in self.profiles['scale'].keys():
				self.profiles['scale'][tt] = 1.
		return

	def _exit(self):
		exit()

	def _make_directories(self):
		for dirs in ['dena']:
			if not os.path.isdir(dirs.upper()): os.mkdir(dirs.upper())
		return

	def _make_home_canvas(self):

		self.figure['name']['home'] = plt.figure(1,figsize=self.figure['size']['home'])
		self.figure['name']['home'].tight_layout()

		self.figure['canvas']['home'] = FigureCanvasTkAgg(self.figure['name']['home'],master=self.root)

		self.figure['widget']['home'] = self.figure['canvas']['home'].get_tk_widget()
		self.figure['widget']['home'].grid(row=1,column=4,rowspan=2,sticky='w')

		toolbar_frame = tk.Frame(self.root)
		toolbar_frame.grid(row=0,column=4)
		
		self.figure['toolbar']['home'] = NavigationToolbar2Tk(self.figure['canvas']['home'],toolbar_frame)

		return

	def _clear_home_canvas(self,clear_list=True):
		self.figure['name']['home'].canvas.draw_idle()
		self.figure['name']['home'].clf()
		self.figure['pegend'] = dict(); self.figure['legend'] = dict();
		self.tciplot = []; self.figure['axes']['home'] = []; self.verplot = [];
		if clear_list: self.tfig_list = []
		self.pfig_list = [];
		return

	def _leftclick(self,event):
		page = self.home['nb'].tk.call(self.home['nb']._w,'identify','tab',event.x,event.y)
		if self.prev_page == page: return
		if page == 1: 
			self._sync_list()
			self._clear_home_canvas()
		elif page == 0: 
			self.tfig_list = []
			self._draw_overview()
			self._draw_tci_fit()
		self.prev_page = page;	
		return

	def _make_note_frame(self):

		titles = ['FIT_PROFILE','POST_PROCESS']
		self.titles = titles
		self.home['nb'] = ttk.Notebook(self.root,width=400,height=720)
		self.home['nb'].bind('<Button-1>',self._leftclick)
		for i in range(1,len(titles)+1):
			self.home['page%i'%i] = ttk.Frame(self.home['nb'])
			self.home['nb'].add(self.home['page%i'%i],text=titles[i-1])
			if i > 1: self.home['nb'].tab(i-1,state='disabled')
			for k in range(10): self.home['page%i'%i].grid_columnconfigure(k,weight=1)

		self.home['nb'].grid(row=1,column=0,columnspan=4,sticky='n')

		return

	def _make_input_frame(self):

		self.l1 = tk.Label(self.home['page1'], text="====================== EFIT INFO. =======================",justify='center')
		self.l1.grid(row=0,column=0,columnspan=10,pady=5)

		self.l1 = tk.Label(self.home['page1'], text="== AVAILABLE ==",justify='center')
		self.l1.grid(row=2,column=0,columnspan=5,pady=5,padx=45,sticky='w')
		self.l1 = tk.Label(self.home['page1'], text="== SELECTED ==",justify='center')
		self.l1.grid(row=2,column=5,columnspan=5,pady=5,padx=45,sticky='e')

		self.note_in['frame1'] = tk.Frame(self.home['page1'])
		self.note_in['frame1'].grid(row=3,column=0,columnspan=5,padx=20,sticky='w')
		self.note_in['s1'] = tk.Scrollbar(self.note_in['frame1'])
		self.note_in['s1'].pack(side='right',fill='y')
		self.note_in['l1'] = tk.Listbox(self.note_in['frame1'],yscrollcommand = self.note_in['s1'].set,height=7)
		self.note_in['l1'].pack(side='left',fill='x')
		self.note_in['s1']["command"] = self.note_in['l1'].yview
	
		self.note_in['frame2'] = tk.Frame(self.home['page1'])
		self.note_in['frame2'].grid(row=3,column=5,columnspan=5,padx=20,sticky='e')
		self.note_in['s2'] = tk.Scrollbar(self.note_in['frame2'])
		self.note_in['s2'].pack(side='right',fill='y')
		self.note_in['l2'] = tk.Listbox(self.note_in['frame2'],yscrollcommand = self.note_in['s2'].set,height=7)
		self.note_in['l2'].pack(side='left',fill='x')
		self.note_in['s2']["command"] = self.note_in['l2'].yview

		self.note_in['frame3'] = tk.Frame(self.home['page1'])
		self.note_in['frame3'].grid(row=4,column=0,columnspan=10)
		self.l1 = tk.Label(self.note_in['frame3'], text="Range[ms,ms]",justify='center')
		self.l1.pack(side='left')
		self.note_in['e1'] = tk.Entry(self.note_in['frame3'],width=13,justify='center')
		self.note_in['e1'].insert(10,'%i, %i'%(self.profiles['tmin'],self.profiles['tmax']))
		self.note_in['e1'].pack(side='left')

		self.l1 = tk.Label(self.note_in['frame3'], text="Tdel [ms]",justify='center')
		self.l1.pack(side='left',padx=10)	
		self.note_in['e2'] = tk.Entry(self.note_in['frame3'],width=9,justify='center')
		self.note_in['e2'].insert(10,self.profiles['delt'])
		self.note_in['e2'].pack(side='left',padx=3)
		self.note_in['e2'].bind('<Return>', lambda x: self._add())

		b1 = tk.Button(self.note_in['frame3'],  text="MAKE", width = 5,command=lambda: self._add())
		b1.pack(side='left',padx=3)

		self.l1 = tk.Label(self.home['page1'], text="==================== DIAGNO INFO. =====================",justify='center')
		self.l1.grid(row=5,column=0,columnspan=10,pady=5)

		flag = ['TS_CORE','TS_EDGE','TCI01','TCI02','TCI03','TCI04','TCI05','INT01','INT02','E_REFL']
		count = 2; rowv = 5
		for j in range(2):
			rowv = rowv + 1
			self.note_in['frame%i'%(j+4)] = tk.Frame(self.home['page1'])
			self.note_in['frame%i'%(j+4)].grid(row=rowv,column=0,columnspan=10)			
			for i in range(5):
				count = count + 1
				self.note_in['e%i'%count] = tk.Entry(self.note_in['frame%i'%(j+4)],width=10,justify='center')
				self.note_in['e%i'%count].insert(10,flag[count-3])
				self.note_in['e%i'%count].pack(side='left')
				self.note_in['e%i'%count].configure(state='readonly',fg='white',readonlybackground='tomato')
	
		if self.ts['core']['nch'] > 0:
			rowv = rowv + 1
			self.l1 = tk.Label(self.home['page1'], text="================ TS-CORE CH.SELECTION =================",justify='center')
			self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)	
			rowv = rowv + 1					
			j = 0;
			for k in range(1,self.ts['core']['nch']+1):
				self.note_in['cts_core%i'%k] = tk.Checkbutton(self.home['page1'],text='%02i'%(k),variable=self.note_in['ts_core%i'%k])
				self.note_in['cts_core%i'%k].grid(row=rowv,column=j)
				j = j + 1; 
				if (j==10): j=0; rowv = rowv + 1;

		if self.ts['edge']['nch'] > 0:
			rowv = rowv + 1
			self.l1 = tk.Label(self.home['page1'], text="================ TS-EDGE CH.SELECTION =================",justify='center')
			self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)				
			rowv = rowv + 1	
			j = 0;
			for k in range(1,self.ts['edge']['nch']+1):
				self.note_in['cts_edge%i'%k] = tk.Checkbutton(self.home['page1'],text='%02i'%(k),variable=self.note_in['ts_edge%i'%k])
				self.note_in['cts_edge%i'%k].grid(row=rowv,column=j)
				j = j + 1; 
				if (j==10): j=0; rowv = rowv + 1;

		rowv = rowv + 1
		self.l1 = tk.Label(self.home['page1'], text="================ TCI/INT/REFL SELECTION =================",justify='center')
		self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)
		rowv = rowv + 1; j= 0;
		for i in range(1,self.tci['nch']+1):
			self.note_in['ctci%i'%i] = tk.Checkbutton(self.home['page1'],text='TCI_%02i'%(i),variable=self.note_in['tci%i'%i])
			self.note_in['ctci%i'%i].grid(row=rowv,column=j,columnspan=2)
			j = j + 2; 
			if (j==10): j=0; rowv = rowv + 1;

		for i in range(1,self.int['nch']+1):
			self.note_in['cint%i'%i] = tk.Checkbutton(self.home['page1'],text='INT_%02i'%(i),variable=self.note_in['int%i'%i])
			self.note_in['cint%i'%i].grid(row=rowv,column=j,columnspan=2)
			j = j + 2; 
			if (j==10): j=0; rowv = rowv + 1;
		self.note_in['cref'] = tk.Checkbutton(self.home['page1'],text='REFL',variable=self.note_in['ref'])
		self.note_in['cref'].grid(row=rowv,column=j,columnspan=2)
	
		self.note_in['b1'] = tk.Button(self.home['page1'],  textvariable=self.checked, width = 8,command=lambda: self._dia_select())
		self.note_in['b1'].grid(row=rowv,column=7,columnspan=3)
		self.note_in['b1'].configure(state='disabled')

		rowv = rowv + 1
		self.l1 = tk.Label(self.home['page1'], text="================= FITTING DIAG. WEIGHT ==================",justify='center')		
		self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)

		labels = ['TS','REFL','INT01','INT02','TCI01','TCI02','TCI03','TCI04','TCI05']
		rowv = rowv + 1
		self.note_in['frame6'] = tk.Frame(self.home['page1'])
		self.note_in['frame6'].grid(row=rowv,column=0,columnspan=10)					
		for i in range(4):
			self.l1 = tk.Label(self.note_in['frame6'], text=labels[i])
			self.l1.pack(side='left')
			self.note_in['e%i'%(i+13)] = tk.Entry(self.note_in['frame6'],textvariable=self.note_in['W%s'%labels[i]],width=5,justify='center')
			self.note_in['e%i'%(i+13)].pack(side='left',padx=10)

		rowv = rowv + 1
		self.l1 = tk.Label(self.home['page1'], text="================== FITTING TCI SIGMA ===================",justify='center')		
		self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)			

		rowv = rowv + 1
		self.note_in['frame7'] = tk.Frame(self.home['page1'])
		self.note_in['frame7'].grid(row=rowv,column=0,columnspan=10)					
		for i in range(4,9):
			self.l1 = tk.Label(self.note_in['frame7'], text=labels[i])
			self.l1.pack(side='left')
			self.note_in['e%i'%(i+13)] = tk.Entry(self.note_in['frame7'],textvariable=self.note_in['W%s'%labels[i]],width=5,justify='center')
			self.note_in['e%i'%(i+13)].pack(side='left')

		labels = ['TS','TCI','INT','REFL']
		rowv = rowv + 1
		self.l1 = tk.Label(self.home['page1'], text="=================== DIAG AVG. DT [ms] ====================",justify='center')		
		self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)

		rowv = rowv + 1
		self.note_in['frame8'] = tk.Frame(self.home['page1'])
		self.note_in['frame8'].grid(row=rowv,column=0,columnspan=10)		
		for i in range(4):
			self.l1 = tk.Label(self.note_in['frame8'], text=labels[i])
			self.l1.pack(side='left')
			self.note_in['e%i'%(i+22)] = tk.Entry(self.note_in['frame8'],textvariable=self.note_in['A%s'%labels[i]],width=5,justify='center')
			self.note_in['e%i'%(i+22)].pack(side='left',padx=10)	


		labels = ['CM','EM']
		rowv = rowv + 1			
		self.l1 = tk.Label(self.home['page1'], text="===================== FITTING OPT. ======================",justify='center')		
		self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)

		rowv = rowv + 1
		self.note_in['frame10'] = tk.Frame(self.home['page1'])
		self.note_in['frame10'].grid(row=rowv,column=0,columnspan=10,sticky='w')	

		self.l1 = tk.Label(self.note_in['frame10'], text='[XMAP]')
		self.l1.pack(side='left',padx=3)

		self.note_in['r6']=tk.Radiobutton(self.note_in['frame10'], text="PSI", value=0, variable=self.note_in['rhofit'])
		self.note_in['r6'].pack(side='left')		

		self.note_in['r7']=tk.Radiobutton(self.note_in['frame10'], text="RHO", value=1, variable=self.note_in['rhofit'])
		self.note_in['r7'].pack(side='left')

		self.l1 = tk.Label(self.note_in['frame10'], text='[EFIT]')
		self.l1.pack(side='left',padx=6)

		self.note_in['r3']=tk.Radiobutton(self.note_in['frame10'], text="01", value=1, variable=self.note_in['gflag'])
		self.note_in['r3'].pack(side='left')		

		self.note_in['r4']=tk.Radiobutton(self.note_in['frame10'], text="02", value=2, variable=self.note_in['gflag'])
		self.note_in['r4'].pack(side='left')

		self.note_in['r5']=tk.Radiobutton(self.note_in['frame10'], text="RT1", value=0, variable=self.note_in['gflag'])
		self.note_in['r5'].pack(side='left')

		for i in range(3,6): self.note_in['r%i'%i].configure(state='disabled')	

		b1 = tk.Button(self.note_in['frame10'],  text="DO FIT", width = 6,command=lambda: self._fit())
		b1.pack(side='left',padx=3)					
		if self.read_mode: b1.config(state='disabled')

		rowv = rowv + 1
		self.note_in['frame9'] = tk.Frame(self.home['page1'])
		self.note_in['frame9'].grid(row=rowv,column=0,columnspan=10,sticky='w')
		for i in range(2):
			self.l1 = tk.Label(self.note_in['frame9'], text=labels[i])
			self.l1.pack(side='left')
			self.note_in['e%i'%(i+26)] = tk.Entry(self.note_in['frame9'],textvariable=self.note_in['%s'%labels[i]],width=4,justify='center')
			self.note_in['e%i'%(i+26)].pack(side='left',padx=2)

		self.l1 = tk.Label(self.note_in['frame9'], text='MINPH')
		self.l1.pack(side='left')
		self.note_in['e29'] = tk.Entry(self.note_in['frame9'],textvariable=self.note_in['MINPH'],width=4,justify='center')
		self.note_in['e29'].pack(side='left',padx=2)				

		self.l1 = tk.Label(self.note_in['frame9'], text='WIDTH')
		self.l1.pack(side='left')
		self.note_in['e30'] = tk.Entry(self.note_in['frame9'],textvariable=self.note_in['WIDTH'],width=5,justify='center')
		self.note_in['e30'].pack(side='left',padx=2)				

		self.note_in['chmode'] = tk.Checkbutton(self.note_in['frame9'],text='HMODE',variable=self.note_in['hmode'])
		self.note_in['chmode'].pack(side='left',padx=2)	

		rowv = rowv + 1
		self.note_in['frame11'] = tk.Frame(self.home['page1'])
		self.note_in['frame11'].grid(row=rowv,column=0,columnspan=10)

		self.note_in['cdocal'] = tk.Checkbutton(self.note_in['frame11'],text='CHCAL',variable=self.note_in['docal'])
		self.note_in['cdocal'].pack(side='left')		

		self.note_in['r1']=tk.Radiobutton(self.note_in['frame11'], text="CORE", value=1, variable=self.note_in['fcore'])
		self.note_in['r1'].pack(side='left')		
			
		self.note_in['m1'] = tk.OptionMenu(self.note_in['frame11'],self.note_in['CALCHC'],*self.ts_core_ch)
		self.note_in['m1'].pack(side='left',padx=2)
		self.note_in['m1'].config(width=1,state='disabled')		

		self.note_in['r2']=tk.Radiobutton(self.note_in['frame11'], text="EDGE", value=0, variable=self.note_in['fcore'])
		self.note_in['r2'].pack(side='left')	

		self.note_in['m2'] = tk.OptionMenu(self.note_in['frame11'],self.note_in['CALCHE'],*self.ts_edge_ch)
		self.note_in['m2'].pack(side='left',padx=2)
		self.note_in['m2'].config(width=1,state='disabled')

		self.l1 = tk.Label(self.note_in['frame11'], text='RATIO')
		self.l1.pack(side='left',padx=2)
		self.note_in['e28'] = tk.Entry(self.note_in['frame11'],textvariable=self.note_in['CALCH'],width=4,justify='center')
		self.note_in['e28'].pack(side='left',padx=2)				

		return

	def _make_post_frame(self):

		self.l1 = tk.Label(self.home['page2'], text="====================== PLOT OPT. =======================",justify='center')
		self.l1.grid(row=0,column=0,columnspan=10,pady=5)

		self.l1 = tk.Label(self.home['page2'], text="== SELECTED ==",justify='center')
		self.l1.grid(row=2,column=0,columnspan=5,pady=5,padx=5)

		self.note_ps['frame1'] = tk.Frame(self.home['page2'])
		self.note_ps['frame1'].grid(row=3,column=0,columnspan=5,padx=5,sticky='e')
		self.note_ps['s1'] = tk.Scrollbar(self.note_ps['frame1'])
		self.note_ps['s1'].pack(side='right',fill='y')
		self.note_ps['l1'] = tk.Listbox(self.note_ps['frame1'],yscrollcommand = self.note_ps['s1'].set,height=13)
		self.note_ps['l1'].pack(side='left',fill='x')
		self.note_ps['s1']["command"] = self.note_ps['l1'].yview
		self.note_ps['l1'].bind('<Double-1>',lambda x: self._click_list())

		self.note_ps['frame2'] = tk.Frame(self.home['page2'])
		self.note_ps['frame2'].grid(row=3,column=5,columnspan=5,padx=1)

		self.note_ps['frame211'] = tk.Frame(self.note_ps['frame2'])
		self.note_ps['frame211'].pack(anchor='w')

		self.l1 = tk.Label(self.note_ps['frame211'], text='[TYPE]')
		self.l1.pack(side='left',padx=3)

		self.note_ps['frame21'] = tk.Frame(self.note_ps['frame2'])
		self.note_ps['frame21'].pack(anchor='w')

		self.note_ps['r1']=tk.Radiobutton(self.note_ps['frame21'], text="1D", value=0, variable=self.note_ps['fig'], command=lambda: self._update_figtype())
		self.note_ps['r1'].pack(side='left')		

		self.note_ps['r2']=tk.Radiobutton(self.note_ps['frame21'], text="2D", value=1, variable=self.note_ps['fig'], command=lambda: self._update_figtype())
		self.note_ps['r2'].pack(side='left')

		self.note_ps['r3']=tk.Radiobutton(self.note_ps['frame21'], text="COM", value=2, variable=self.note_ps['fig'], command=lambda: self._update_figtype())
		self.note_ps['r3'].pack(side='left')	

		self.note_ps['r32']=tk.Radiobutton(self.note_ps['frame21'], text="DIFF", value=3, variable=self.note_ps['fig'], command=lambda: self._update_figtype())
		self.note_ps['r32'].pack(side='left')	

		self.note_ps['frame212'] = tk.Frame(self.note_ps['frame2'])
		self.note_ps['frame212'].pack(anchor='w')		

		self.note_ps['r7']=tk.Radiobutton(self.note_ps['frame212'], text="SCAN", value=4, variable=self.note_ps['fig'], command=lambda: self._update_figtype())
		self.note_ps['r7'].pack(side='left')		

		self.note_ps['r72']=tk.Radiobutton(self.note_ps['frame212'], text="TSCAN", value=5, variable=self.note_ps['fig'], command=lambda: self._update_figtype())
		self.note_ps['r72'].pack(side='left')		

		self.note_ps['frame22'] = tk.Frame(self.note_ps['frame2'])
		self.note_ps['frame22'].pack(anchor='w')

		self.l1 = tk.Label(self.note_ps['frame22'], text='[XMAP]')
		self.l1.pack(side='left',padx=3)

		self.note_ps['r4']=tk.Radiobutton(self.note_ps['frame22'], text="PSI", value=0, variable=self.note_ps['xmap'], command=lambda: self._update_xmap())
		self.note_ps['r4'].pack(side='left')		

		self.note_ps['r5']=tk.Radiobutton(self.note_ps['frame22'], text="RHO", value=1, variable=self.note_ps['xmap'], command=lambda: self._update_xmap())
		self.note_ps['r5'].pack(side='left')

		self.note_ps['r6']=tk.Radiobutton(self.note_ps['frame22'], text="R",   value=2, variable=self.note_ps['xmap'], command=lambda: self._update_xmap())
		self.note_ps['r6'].pack(side='left')

		self.note_ps['frame23'] = tk.Frame(self.note_ps['frame2'])
		self.note_ps['frame23'].pack(anchor='w')

		self.l1 = tk.Label(self.note_ps['frame23'], text='XMIN')
		self.l1.pack(side='left',padx=3)
		self.note_ps['e1'] = tk.Entry(self.note_ps['frame23'],textvariable=self.note_ps['xmin'],width=6,justify='center')
		self.note_ps['e1'].pack(side='left',padx=3)
		self.l1 = tk.Label(self.note_ps['frame23'], text='XMAX')
		self.l1.pack(side='left',padx=3)
		self.note_ps['e2'] = tk.Entry(self.note_ps['frame23'],textvariable=self.note_ps['xmax'],width=6,justify='center')
		self.note_ps['e2'].pack(side='left',padx=3)

		self.note_ps['frame24'] = tk.Frame(self.note_ps['frame2'])
		self.note_ps['frame24'].pack(anchor='w')
		self.l1 = tk.Label(self.note_ps['frame24'], text='YMIN')
		self.l1.pack(side='left',padx=3)
		self.note_ps['e3'] = tk.Entry(self.note_ps['frame24'],textvariable=self.note_ps['ymin'],width=6,justify='center')
		self.note_ps['e3'].pack(side='left',padx=3)
		self.l1 = tk.Label(self.note_ps['frame24'], text='YMAX')
		self.l1.pack(side='left',padx=3)
		self.note_ps['e4'] = tk.Entry(self.note_ps['frame24'],textvariable=self.note_ps['ymax'],width=6,justify='center')
		self.note_ps['e4'].pack(side='left',padx=3)		

		self.note_ps['frame25'] = tk.Frame(self.note_ps['frame2'])
		self.note_ps['frame25'].pack(anchor='e',pady=5)

		b1 = tk.Button(self.note_ps['frame25'],  text="APPLY", width = 6,command=lambda: self._apply_axis())
		b2 = tk.Button(self.note_ps['frame25'],  text="RESET", width = 6,command=lambda: self._reset_plot())
		b2.pack(side='left',padx=3)		
		b1.pack(side='left',padx=3)

		self.note_ps['frame26'] = tk.Frame(self.note_ps['frame2'])
		self.note_ps['frame26'].pack(anchor='e',pady=5)
		b3 = tk.Button(self.note_ps['frame26'],  text="DO PLOT", width = 7,command=lambda: self._draw_plots())
		b3.pack(side='left',padx=3)		

		#if self.ts['core']['nch'] == 0: return

		rowv = 3
		if self.ts['core']['nch'] > 0:
			rowv = rowv + 1
			self.l1 = tk.Label(self.home['page2'], text="================ TS-CORE CH.SELECTION =================",justify='center')
			self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)	
			rowv = rowv + 1					
			j = 0;
			for k in range(1,self.ts['core']['nch']+1):
				self.note_ps['cts_core%i'%k] = tk.Checkbutton(self.home['page2'],text='%02i'%(k),variable=self.note_in['ts_core%i'%k])
				self.note_ps['cts_core%i'%k].grid(row=rowv,column=j)
				j = j + 1; 
				if (j==10): j=0; rowv = rowv + 1;

		if self.ts['edge']['nch'] > 0:
			rowv = rowv + 1
			self.l1 = tk.Label(self.home['page2'], text="================ TS-EDGE CH.SELECTION =================",justify='center')
			self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)				
			rowv = rowv + 1	
			j = 0;
			for k in range(1,self.ts['edge']['nch']+1):
				self.note_ps['cts_edge%i'%k] = tk.Checkbutton(self.home['page2'],text='%02i'%(k),variable=self.note_in['ts_edge%i'%k])
				self.note_ps['cts_edge%i'%k].grid(row=rowv,column=j)
				j = j + 1; 
				if (j==10): j=0; rowv = rowv + 1;		

		rowv = rowv + 1
		self.l1 = tk.Label(self.home['page2'], text="============== TS CORE/EDGE VIEW OPTION ===============",justify='center')
		self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)		

		rowv = rowv + 1
		if self.ts['core']['nch'] > 0:
			rowv = rowv + 1					
			j = 0;
			for k in range(1,self.ts['core']['nch']+1):
				self.note_ps['cts_core%i_2'%k] = tk.Checkbutton(self.home['page2'],text='%02i'%(k),variable=self.note_ps['ts_core%i'%k],command=lambda: self._checkbutton_limit())
				self.note_ps['cts_core%i_2'%k].grid(row=rowv,column=j)
				j = j + 1; 
				if (j==10): j=0; rowv = rowv + 1;

		if self.ts['edge']['nch'] > 0:
			rowv = rowv + 1
			self.l1 = tk.Label(self.home['page2'], text="-- ",anchor='w')
			self.l1.grid(row=rowv,column=0,columnspan=2)				
			rowv = rowv + 1	
			j = 0;
			for k in range(1,self.ts['edge']['nch']+1):
				self.note_ps['cts_edge%i_2'%k] = tk.Checkbutton(self.home['page2'],text='%02i'%(k),variable=self.note_ps['ts_edge%i'%k],command=lambda: self._checkbutton_limit())
				self.note_ps['cts_edge%i_2'%k].grid(row=rowv,column=j)
				j = j + 1; 
				if (j==10): j=0; rowv = rowv + 1;	

		rowv = rowv + 1
		self.l1 = tk.Label(self.home['page2'], text="================= PLOT DIAGNO. OPTION ==================",justify='center')
		self.l1.grid(row=rowv,column=0,columnspan=10,pady=5)								

		rowv = rowv + 1
		self.note_ps['frame3'] = tk.Frame(self.home['page2'])
		self.note_ps['frame3'].grid(row=rowv,column=0,columnspan=10)	

		self.note_ps['cdocal'] = tk.Checkbutton(self.note_ps['frame3'],text='CHCAL',variable=self.note_in['docal'])
		self.note_ps['cdocal'].pack(side='left')		

		self.note_ps['r10']=tk.Radiobutton(self.note_ps['frame3'], text="CORE", value=1, variable=self.note_in['fcore'])
		self.note_ps['r10'].pack(side='left')					

		self.note_ps['m1'] = tk.OptionMenu(self.note_ps['frame3'],self.note_in['CALCHC'],*self.ts_core_ch)
		self.note_ps['m1'].pack(side='left',padx=3)
		self.note_ps['m1'].config(width=1,state='disabled')		

		self.note_ps['r20']=tk.Radiobutton(self.note_ps['frame3'], text="EDGE", value=0, variable=self.note_in['fcore'])
		self.note_ps['r20'].pack(side='left')	

		self.note_ps['m2'] = tk.OptionMenu(self.note_ps['frame3'],self.note_in['CALCHE'],*self.ts_edge_ch)
		self.note_ps['m2'].pack(side='left',padx=3)
		self.note_ps['m2'].config(width=1,state='disabled')

		self.l1 = tk.Label(self.note_ps['frame3'], text='RATIO')
		self.l1.pack(side='left')
		self.note_ps['e1'] = tk.Entry(self.note_ps['frame3'],textvariable=self.note_in['CALCH'],width=5,justify='center')
		self.note_ps['e1'].pack(side='left',padx=2)			

		rowv = rowv + 1
		self.note_ps['frame6'] = tk.Frame(self.home['page2'])
		self.note_ps['frame6'].grid(row=rowv,column=0,columnspan=10,sticky='w')	
		self.l1 = tk.Label(self.note_ps['frame6'], text='COREM')
		self.l1.pack(side='left',padx=5)
		self.note_ps['e4'] = tk.Entry(self.note_ps['frame6'],textvariable=self.note_in['CM'],width=6,justify='center')
		self.note_ps['e4'].pack(side='left',padx=5)				
		self.l1 = tk.Label(self.note_ps['frame6'], text='EDGEM')
		self.l1.pack(side='left',padx=5)
		self.note_ps['e5'] = tk.Entry(self.note_ps['frame6'],textvariable=self.note_in['EM'],width=6,justify='center')
		self.note_ps['e5'].pack(side='left',padx=5)			

		self.l1 = tk.Label(self.note_ps['frame6'], text='CN')
		self.l1.pack(side='left',padx=5)
		self.note_ps['e6'] = tk.Entry(self.note_ps['frame6'],textvariable=self.note_ps['CORE_N'],width=4,justify='center')
		self.note_ps['e6'].pack(side='left',padx=5)			
		self.l1 = tk.Label(self.note_ps['frame6'], text='EN')
		self.l1.pack(side='left',padx=5)
		self.note_ps['e7'] = tk.Entry(self.note_ps['frame6'],textvariable=self.note_ps['EDGE_N'],width=4,justify='center')
		self.note_ps['e7'].pack(side='left',padx=5)							

		rowv = rowv + 1; count = 8
		self.note_ps['frame7'] = tk.Frame(self.home['page2'])
		self.note_ps['frame7'].grid(row=rowv,column=0,columnspan=10,sticky='w')	
		flist = ['CMIN','CMAX','EMIN','EMAX']; vlist = ['CORE_MIN','CORE_MAX','EDGE_MIN','EDGE_MAX']
		for k in range(len(flist)):

			self.l1 = tk.Label(self.note_ps['frame7'], text=flist[k])
			self.l1.pack(side='left',padx=2)			
			self.note_ps['e%i'%count] = tk.Entry(self.note_ps['frame7'],textvariable=self.note_ps[vlist[k]],width=5,justify='center')
			self.note_ps['e%i'%count].pack(side='left',padx=2)		
			count = count + 1

		rowv = rowv + 1
		self.note_ps['frame8'] = tk.Frame(self.home['page2'])
		self.note_ps['frame8'].grid(row=rowv,column=0,columnspan=10,sticky='w')	
		self.l1 = tk.Label(self.note_ps['frame8'], text='COREW')
		self.l1.pack(side='left',padx=5)
		self.note_ps['e4'] = tk.Entry(self.note_ps['frame8'],textvariable=self.note_ps['COREM'],width=6,justify='center')
		self.note_ps['e4'].pack(side='left',padx=5)				
		self.l1 = tk.Label(self.note_ps['frame8'], text='EDGEW')
		self.l1.pack(side='left',padx=5)
		self.note_ps['e5'] = tk.Entry(self.note_ps['frame8'],textvariable=self.note_ps['EDGEM'],width=6,justify='center')
		self.note_ps['e5'].pack(side='left',padx=5)						

		rowv = rowv + 1
		self.note_ps['frame9'] = tk.Frame(self.home['page2'])
		self.note_ps['frame9'].grid(row=rowv,column=0,columnspan=10,sticky='e')

		b1 = tk.Button(self.note_ps['frame9'],  text="AUTO-SCALE", width = 10,command=lambda: self._auto_scale())
		b1.pack(side='left',padx=5)
						
		return

	def _checkbutton_limit(self):

		ch_index = []

		for i in range(self.ts['core']['nch']):
			if self.note_ps['ts_core%i'%(i+1)].get() == 1:
				ch_index.append(i)
		for i in range(self.ts['edge']['nch']):
			if self.note_ps['ts_edge%i'%(i+1)].get() == 1:
				ch_index.append(i+self.ts['core']['nch'])
		if len(ch_index) < 7:
			self.check_list = copy.deepcopy(ch_index)
		else:
			for ind in ch_index:
				try: index = self.check_list.index(ind)
				except: break
			if ind < self.ts['core']['nch']: self.note_ps['ts_core%i'%(ind+1)].set(0)
			else:
				jj = ind - self.ts['core']['nch']
				self.note_ps['ts_edge%i'%(jj+1)].set(0)
				print('>>> Maximum 6 selection !')

		return

	def _dia_select(self):

		if self.checked.get() == 'CHECK-A': val = 1; tag = 'UNCHECK-A'
		else: val = 0; tag = 'CHECK-A'

		for i in range(1,self.tci['nch']+1):
			if self.tci[i]['is']: self.note_in['tci%i'%(i)].set(val)
		for i in range(1,self.int['nch']+1):
			if self.int[i]['is']: self.note_in['int%i'%(i)].set(val)
		if self.ref['is']: self.note_in['ref'].set(val)
		self.checked.set(tag)
			
		return

	def _draw_overview(self):

		self._clear_home_canvas()
		self.figure['axes']['home'] = [self.figure['name']['home'].add_subplot(self.figure['gs']['home'][0,0])]
		self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(self.figure['gs']['home'][0,1],sharex=self.figure['axes']['home'][0]))
		for k in range(1,self.figure['nrow']):
			self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(self.figure['gs']['home'][k,0],sharex=self.figure['axes']['home'][0]))
			self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(self.figure['gs']['home'][k,1],sharex=self.figure['axes']['home'][0]))

		for k in range(self.figure['nrow']*2-2): self.figure['axes']['home'][k].axes.get_xaxis().set_visible(False)

		self.figure['name']['homecursor'] = MultiCursor(self.figure['name']['home'].canvas,self.figure['axes']['home'],horizOn=True,color='r',lw=1)
		line1, = self.figure['axes']['home'][0].plot(self.mds['ip'][0],abs(self.mds['ip'][1])/1.e6)
		line2, = self.figure['axes']['home'][0].plot(self.mds['da'][0],abs(self.mds['da'][1])/6.e21)
		slegend1 = [line1,line2]; plegend1=['$Ip$ [MA]','$D_\\alpha$ [a.u.]']

		if self.ts['core']['nch'] > 0: 
			line3, = self.figure['axes']['home'][0].plot(self.ts['core'][1][0],self.ts['core'][1][1]/1.e20)
			slegend1.append(line3); plegend1.append('TS_CORE02 [$10^{20}/m^3$]')
		self.figure['axes']['home'][0].legend(slegend1,plegend1,loc='upper right')

		if len(self.profiles['times']) > 0:
			self.verplot = [self.figure['axes']['home'][0].axvline(x=self.profiles['times'][0]/1000,color='magenta',linestyle='--',linewidth=1.0,zorder=5)]
			self.verplot.append(self.figure['axes']['home'][0].axvline(x=self.profiles['times'][-1]/1000,color='magenta',linestyle='--',linewidth=1.0,zorder=5))

		count = 1; 
		for ch in range(1,self.tci['nch']+1):
			if self.tci[ch]['is']: 
				lenx = len(self.tci[ch]['val'][0]); lenv = len(self.tci[ch]['val'][1]); lenm = min(lenx,lenv)
				self.tci[ch]['val'][0] = self.tci[ch]['val'][0][0:lenm]; self.tci[ch]['val'][1] = self.tci[ch]['val'][1][0:lenm]
				line1, = self.figure['axes']['home'][count].plot(self.tci[ch]['val'][0],self.tci[ch]['val'][1])
				self.figure['pegend'][count] = [line1]; self.figure['legend'][count] = ['TCI%02i [$10^{19}/m^3$]'%ch];
				self.figure['axes']['home'][count].legend(self.figure['pegend'][count],self.figure['legend'][count],loc='upper right')
				self.figure['axes']['home'][count].axhline(y=0.,color='gold',linestyle='--',linewidth=1.0,zorder=5)
				count = count + 1;

		for ch in range(1,self.int['nch']+1):
			if self.int[ch]['is']: 
				line1, = self.figure['axes']['home'][count].plot(self.int[ch]['val'][0],self.int[ch]['val'][1])
				self.figure['pegend'][count] = [line1]; self.figure['legend'][count] = ['INT%02i [$10^{19}/m^3$]'%ch];
				self.figure['axes']['home'][count].legend(self.figure['pegend'][count],self.figure['legend'][count],loc='upper right')
				self.figure['axes']['home'][count].axhline(y=0.,color='gold',linestyle='--',linewidth=1.0,zorder=5)
				count = count + 1;

		self.figure['axes']['home'][0].set_xlim(self.tmin/1000.,self.tmax/1000.+0.5)

		self.figure['axes']['home'][0].set_ylim([0,0.8])
		for k in range(1,self.figure['nrow']*2): self.figure['axes']['home'][k].set_ylim(-0.1)
		self.figure['axes']['home'][self.figure['nrow']*2-2].set_xlabel('time [s]')
		self.figure['axes']['home'][self.figure['nrow']*2-1].set_xlabel('time [s]')

		self.figure['name']['home'].tight_layout()

		return

	def _update_diagno_state(self):
		
		self.note_in['CALCHC'].set('1'); self.note_in['CALCHE'].set('1');

		if self.ts['core']['nch'] > 0: 
			self.note_in['e3'].configure(readonlybackground='steelblue')
			for i in range(1,self.ts['core']['nch']+1): self.note_in['ts_core%i'%i].set(1)
			for i in self.ts_core_ex_preset: self.note_in['ts_core%i'%i].set(0)
			self.note_in['CALCHC'].set('%02i'%(self.ts['core']['nch']-1)); self.note_in['m1'].config(state='normal'); self.note_ps['m1'].config(state='normal');

			for ch in self.ts_core_ex_preset:
				self.note_in['ts_core%i'%(ch+1)].set(0)
			self.note_in['b1'].configure(state='normal')
		else:
			self.note_in['e13'].configure(state='readonly')
			self.note_ps['r3'].configure(state='disabled')
			self.note_ps['r32'].configure(state='disabled')
			self.note_ps['r7'].configure(state='disabled')
			self.note_ps['r72'].configure(state='disabled')
			self.note_in['cdocal'].configure(state='disabled')
			self.note_ps['cdocal'].configure(state='disabled')
			self.note_in['r1'].configure(state='disabled')
			self.note_in['r2'].configure(state='disabled')
			self.note_ps['r10'].configure(state='disabled')
			self.note_ps['r20'].configure(state='disabled')

		if self.ts['edge']['nch'] > 0: 
			self.note_in['e4'].configure(readonlybackground='steelblue')
			for i in range(1,self.ts['edge']['nch']+1): self.note_in['ts_edge%i'%i].set(1)
			for i in self.ts_edge_ex_preset: self.note_in['ts_edge%i'%i].set(0)
			self.note_in['CALCHE'].set('%02i'%(1)); self.note_in['m2'].config(state='normal'); self.note_ps['m2'].config(state='normal');

			for ch in self.ts_edge_ex_preset:
				self.note_in['ts_edge%i'%(ch+1)].set(0)			

		count = 1; indcount = 17;
		for i in range(self.tci['nch']): 
			if self.tci[i+1]['is']: 
				self.note_in['e%i'%(5+i)].configure(readonlybackground='steelblue')
				self.note_in['tci%i'%(i+1)].set(1)
				self.note_in['WTCI%02i'%(i+1)].set('%2.1f'%self.tcis_preset[i])
				count = count + 1; 
				indcount = indcount + 1;
			else:
				self.note_in['ctci%i'%(i+1)].configure(state='disabled')
				self.note_in['e%i'%indcount].configure(state='readonly')
				indcount = indcount + 1;				


		indcount = 15;
		for i in range(self.int['nch']): 
			if self.int[i+1]['is']: 
				self.note_in['e%i'%(5+i)].configure(readonlybackground='steelblue')
				self.note_in['int%i'%(i+1)].set(1)
				self.note_in['WINT%02i'%(i+1)].set('%2.1f'%self.ints_preset[i])
				count = count + 1; 	
				indcount = indcount + 1;				
			else:
				self.note_in['cint%i'%(i+1)].configure(state='disabled')
				self.note_in['e%i'%indcount].configure(state='readonly')
				indcount = indcount + 1;							

		if self.ref['is']: 
				self.note_in['e10'].configure(readonlybackground='steelblue')
				self.note_in['ref'].set(1)
		else: 
			self.note_in['cref'].configure(state='disabled')
			self.note_in['e14'].configure(state='readonly')

		count = int(np.ceil(count/2.))
		self.figure['gs']['home'] = gridspec.GridSpec(count, 2,hspace=0.)
		self.figure['nrow'] = count

#		if self.out_mse == None: self.note_in['r4'].configure(state='disabled')
#		if self.out_rt1 == None: self.note_in['r5'].configure(state='disabled')

		return

	def _update_time_list(self):

		if self.note_in['gflag'].get() > 0: elist = self.efit_list['times'][self.note_in['gflag'].get()]
		else: elist = self.out_rt1;
		self.note_in['l1'].delete(0,'end')
		self.tlist = np.array([])
		tmin = 100000; tmax = -1;
		for k in elist: 
			tt = float(k)
			if (tt>=self.tmin and tt<=self.tmax):
				tmin = min(tt,tmin); tmax = max(tt,tmax)
				self.tlist = np.append(self.tlist,int(tt))
				line = '  %06i ms'%int(float(k))
				self.note_in['l1'].insert('end',line)
		self.profiles['tmin'] = tmin+1000; self.profiles['tmax'] = tmax; #-500;
		self.note_in['e1'].delete(0,'end'); self.note_in['e1'].insert(10,'%i, %i'%(tmin+1000,tmax)) #-500))

		return

	def _sync_list(self):

		tlist = list(self.note_in['l2'].get(0,'end'))
		self.note_ps['l1'].delete(0,'end')
		for tt in range(len(tlist)):
			self.note_ps['l1'].insert('end',tlist[tt])

		return

	def _add(self):
		self.note_in['l2'].delete(0,'end')
		self.profiles['times']		= np.array([],dtype='int')
		try: 
			self.profiles['tmin'] = int(float(self.note_in['e1'].get().split(',')[0]))
			self.profiles['tmax'] = int(float(self.note_in['e1'].get().split(',')[1]))
			self.profiles['delt'] = int(float(self.note_in['e2'].get()))
			if self.profiles['tmax']<=self.profiles['tmin']: return
			if self.profiles['delt'] <=0: return
			nrange = int((self.profiles['tmax'] - self.profiles['tmin'])/self.profiles['delt'])+1
			trange = np.linspace(self.profiles['tmin'],self.profiles['tmax'],nrange)
		except: return
		tt = -1;
		for time in trange:
			tt2 = self.tlist[np.argmin(abs(self.tlist-time))]
			if tt2 > tt:
				line = '  %06i ms'%int(tt2)
				self.note_in['l2'].insert('end',line)
				tt = int(tt2);
				self.profiles['times']	= np.append(self.profiles['times'],tt)

		if self.brutal: return

		self.figure['name']['home'].canvas.draw_idle()
		if len(self.verplot) > 0:
			self.verplot[0].remove(); self.verplot[1].remove()
		if self.brutal: return
		self.verplot = [self.figure['axes']['home'][0].axvline(x=self.profiles['times'][0]/1000,color='magenta',linestyle='--',linewidth=1.0,zorder=5)]
		self.verplot.append(self.figure['axes']['home'][0].axvline(x=self.profiles['times'][-1]/1000,color='magenta',linestyle='--',linewidth=1.0,zorder=5))	

		self.home['nb'].tab(1,state='disabled')		
		return

	def _fit(self):

		if self.read_mode:
			self.profiles['didfit'] = True
		else:
			tlen = len(self.profiles['times'])
			if tlen ==0: return
			self._set_fitopt()
			self._get_gfiles()
			self._make_ts_raw_files()
			self._make_tci_raw_files()
			for k in range(tlen): 
				self._dofit(k)
				update_progress(float((k+1)/tlen))
			self.profiles['didfit'] = True
			self._dump_profiles('dena_%05i.sav'%self.shotn)

		if not self.nogui:
			self._draw_tci_fit()
			self.home['nb'].tab(1,state='normal')

		return

	def _dump_profiles(self,filename):
		f = open(filename,'wb')
		pickle.dump(self.profiles,f)
		f.close()
		return

	def _draw_tci_fit(self):

		if not self.profiles['didfit']: return
		self.figure['name']['home'].canvas.draw_idle()
		lenp = len(self.tciplot)
		if lenp > 0:
			for i in range(lenp):
				plot = self.tciplot[i]; plot.remove(); 
				self.figure['pegend'][i+1].remove(plot); self.figure['legend'][i+1].remove('FIT');

		self.tciplot = []
		tlen = len(self.profiles['times'])
		count = 1;

		for ch in range(1,self.tci['nch']+1):
			if self.tci[ch]['is']:
				neline = np.zeros(tlen)
				for k in range(tlen): neline[k] = self.profiles['fit']['nel'][self.profiles['times'][k]][ch+1]
				line1, = self.figure['axes']['home'][count].plot(self.profiles['times']/1000,neline,linestyle='--',marker='x',color='r')
				self.figure['pegend'][count].append(line1); self.figure['legend'][count].append('FIT');
				self.figure['axes']['home'][count].legend(self.figure['pegend'][count],self.figure['legend'][count],loc='upper right')
				self.tciplot.append(line1)
				count = count + 1;

		for ch in range(1,self.int['nch']+1):
			if self.int[ch]['is']: 
				neline = np.zeros(tlen)
				for k in range(tlen): neline[k] = self.profiles['fit']['nel'][self.profiles['times'][k]][ch-1]
				line1, = self.figure['axes']['home'][count].plot(self.profiles['times']/1000,neline,linestyle='--',marker='x',color='r')
				self.figure['pegend'][count].append(line1); self.figure['legend'][count].append('FIT');
				self.figure['axes']['home'][count].legend(self.figure['pegend'][count],self.figure['legend'][count],loc='upper right')				
				self.tciplot.append(line1)			
				count = count + 1;		
		return

	def _apply_axis(self):

		if self.note_ps['fig'].get() > 1: return

		self.figure['name']['home'].canvas.draw_idle()		
		xmin = float(self.note_ps['xmin'].get());
		xmax = float(self.note_ps['xmax'].get());
		ymin = float(self.note_ps['ymin'].get());
		ymax = float(self.note_ps['ymax'].get());
		self.figure['axes']['home'][0].set_xlim(xmin,xmax)
		self.figure['axes']['home'][0].set_ylim(ymin,ymax)

		return

	def _reset_plot(self):

		self.figure['name']['home'].canvas.draw_idle()
		if len(self.tfig_list) == 0: return
		self._clear_home_canvas()
		self.tfig_list = []; self.pfig_list = [];

		return

	def _update_figtype(self):

		if self.prev_plot_type == self.note_ps['fig'].get(): return
		self.prev_plot_type = self.note_ps['fig'].get()
		self._clear_home_canvas()
		return

	def _update_xmap(self):

		if self.prev_xmap == self.note_ps['xmap'].get(): return
		self.prev_xmap = self.note_ps['xmap'].get()

		if len(self.figure['axes']['home']) == 0: return

		if self.note_ps['fig'].get() == 0: self._update_profiles_1d()
		if self.note_ps['fig'].get() == 1: self._draw_profiles_2d(update=True)

		return

	def _click_list(self):

		selection = self.note_ps['l1'].curselection()
		if len(selection) == 0: return
		line = self.note_ps['l1'].get(selection[0])
		time = int(float(line.split()[0]))

		try: self.prev_vline.remove()
		except: pass

		if self.note_ps['fig'].get() == 0: self._draw_profiles_1d()
		elif self.note_ps['fig'].get() == 1: 
			if len(self.figure['axes']['home']) == 0: return
			self.prev_vline = self.figure['axes']['home'][0].axvline(x=time/1000.,linestyle='--',color='magenta')
		return

	def _draw_plots(self,update=False):
		if self.note_ps['fig'].get() == 0:
			if not update: self._draw_profiles_1d()
			else: self._update_profiles_1d()
			
		elif self.note_ps['fig'].get() == 1:
			self._draw_profiles_2d()
		elif self.note_ps['fig'].get() == 2:
			self._draw_channel_gap()
		elif self.note_ps['fig'].get() == 3:
			self._draw_channel_gap(diff=True)
		elif self.note_ps['fig'].get() == 4:			
			self._draw_channel_gap2d()
		elif self.note_ps['fig'].get() == 5:
			self._draw_channel_gap2dt()
		return

	def _draw_profiles_1d(self):

		if not self.note_ps['fig'].get() == 0: return
		self.prev_plot_type = self.note_ps['fig'].get()

		selection = self.note_ps['l1'].curselection()
		if len(selection) == 0: return
		line = self.note_ps['l1'].get(selection[0])
		time = int(float(line.split()[0]))
		try: self.tfig_list.index(time); return
		except: pass 

		pnum = len(self.tfig_list);
		ists = False; 
		if self.ts['core']['nch']>0: ists = True;

		self.figure['name']['home'].canvas.draw_idle()
		xlabel = '$\\psi_N [a.u]$';
		if self.note_ps['xmap'].get() == 1: xlabel = '$\\rho_N [a.u]$';
		elif self.note_ps['xmap'].get() == 2: xlabel = '$R [m]$';
		if (pnum > 8): self.figure['axes']['home'][0].cla(); self.tfig_list = []; self.pfig_list= []; pnum = 0;
		if (pnum == 0):			
			self.figure['name']['homecursor'].visible = False	
			self.figure['axes']['home'] = [self.figure['name']['home'].add_subplot(1,1,1)]
			self.figure['axes']['home'][0].set_xlabel(xlabel)
			self.figure['axes']['home'][0].set_ylabel('Density [$10^{19}/m^3$]')
			self.figure['axes']['home'][0].axhline(y=0,linestyle='--',c='gold')

		nch = self.ts['core']['nch'] + self.ts['edge']['nch']
		vals = self.profiles['fit']['ne'][time]
		if self.note_ps['xmap'].get() < 2: 
			xfit = vals[0]; yfit = vals[1]
			if ists: xdat = vals[2]['xxr'][:nch]; ydat = vals[2]['raw'][:nch]; sdat = vals[2]['raws'][:nch]
			if self.note_ps['xmap'].get() == 1: 
				xfit = vals[3](xfit); 
				if ists: xdat = vals[3](xdat)
		else:
			if ists: 
				xdat = vals[2]['datR'][:nch]; ydat = vals[2]['raw'][:nch]; sdat = vals[2]['raws'][:nch]
				minR = min(xdat)-0.05; maxR = max(xdat); xfit = np.linspace(minR,maxR,50);
			else: xfit = np.linspace(1.75,2.3)
			xx = vals[4](xfit,vals[2]['datZ'][0]);
			for i in range(len(xx)): xx[i] = min(xx[i],1.1)
			yf = interp1d(vals[0],vals[1]); yfit = yf(xx)

		if self.ts['core']['nch']>0: ydat = self._do_scale(time,ydat)
		if pnum == 0: 
			for i in range(self.ts['core']['nch']): self.figure['axes']['home'][0].text(xdat[i],ydat[i],'C%02i'%(i+1))
			for i in range(self.ts['edge']['nch']):  
				j = self.ts['core']['nch'] + i;
				self.figure['axes']['home'][0].text(xdat[j],ydat[j],'E%02i'%(i+1))

		line1, = self.figure['axes']['home'][0].plot(xfit,yfit,linestyle='--',color='C%i'%(pnum+1))
		if ists: line2 = self.figure['axes']['home'][0].errorbar(xdat,ydat,sdat,fmt='x',c='C%i'%(pnum+1))
		self.pfig_list.append(line1)

		if pnum == 0:
			if self.note_ps['xmap'].get() < 2: self.figure['axes']['home'][0].set_xlim(0.,1.1); self.note_ps['xmin'].set('0.0'); self.note_ps['xmax'].set('1.1');
			else: self.figure['axes']['home'][0].set_xlim(1.7,2.3); self.note_ps['xmin'].set('1.7'); self.note_ps['xmax'].set('2.3');
			self.figure['axes']['home'][0].set_ylim(0.)
			self.figure['name']['home'].tight_layout()

		llegend = [];
		self.tfig_list.append(time)
		for tt in self.tfig_list: llegend.append('%g ms'%tt)
		self.figure['axes']['home'][0].legend(self.pfig_list,llegend)
		self.prev_xmap = self.note_ps['xmap'].get()
		return

	def _update_profiles_1d(self):

		pnum = len(self.tfig_list);
		if pnum == 0: return
		self._clear_home_canvas(False)
		self.figure['name']['home'].canvas.draw_idle()
		xlabel = '$\\psi_N [a.u]$';
		if self.note_ps['xmap'].get() == 1: xlabel = '$\\rho_N [a.u]$';
		elif self.note_ps['xmap'].get() == 2: xlabel = '$R [m]$';		

		ists = False; 
		if self.ts['core']['nch']>0: ists = True;
		
		self.figure['name']['homecursor'].visible = False	
		self.figure['axes']['home'] = [self.figure['name']['home'].add_subplot(1,1,1)]
		self.figure['axes']['home'][0].set_xlabel(xlabel)
		self.figure['axes']['home'][0].set_ylabel('Density [$10^{19}/m^3$]')
		self.figure['axes']['home'][0].axhline(y=0,linestyle='--',c='gold')		
		count = 0
		nch = self.ts['core']['nch'] + self.ts['edge']['nch']
		for time in self.tfig_list:
			count = count + 1
			vals = self.profiles['fit']['ne'][time]
			if self.note_ps['xmap'].get() < 2: 
				xfit = vals[0]; yfit = vals[1]
				if ists: xdat = vals[2]['xxr'][:nch]; ydat = vals[2]['raw'][:nch]; sdat = vals[2]['raws'][:nch]
				if self.note_ps['xmap'].get() == 1: 
					xfit = vals[3](xfit); 
					if ists: xdat = vals[3](xdat)
			else:
				if ists:
					xdat = vals[2]['datR'][:nch]; ydat = vals[2]['raw'][:nch]; sdat = vals[2]['raws'][:nch]
					minR = min(xdat)-0.05; maxR = max(xdat); xfit = np.linspace(minR,maxR,50);
				else: xfit = np.linspace(1.75,2.3)
				xx = vals[4](xfit,vals[2]['datZ'][0]);
				for i in range(len(xx)): xx[i] = min(xx[i],1.1)
				yf = interp1d(vals[0],vals[1]); yfit = yf(xx)

			if ists: ydat = self._do_scale(time,ydat)
			line1, = self.figure['axes']['home'][0].plot(xfit,yfit,linestyle='--',color='C%i'%count)
			if ists: line2 = self.figure['axes']['home'][0].errorbar(xdat,ydat,sdat,fmt='x',c='C%i'%count)
			self.pfig_list.append(line1)

		for i in range(self.ts['core']['nch']): self.figure['axes']['home'][0].text(xdat[i],ydat[i],'C%02i'%(i+1))
		for i in range(self.ts['edge']['nch']):  
			j = self.ts['core']['nch'] + i;
			self.figure['axes']['home'][0].text(xdat[j],ydat[j],'E%02i'%(i+1))

		if self.note_ps['xmap'].get() < 2: self.figure['axes']['home'][0].set_xlim(0.,1.1); self.note_ps['xmin'].set('0.0'); self.note_ps['xmax'].set('1.1');
		else: self.figure['axes']['home'][0].set_xlim(1.7,2.3); self.note_ps['xmin'].set('1.7'); self.note_ps['xmax'].set('2.3');
		self.figure['axes']['home'][0].set_ylim(0.)
		self.figure['name']['home'].tight_layout()			

		llegend = [];
		for tt in self.tfig_list: llegend.append('%g ms'%tt)
		self.figure['axes']['home'][0].legend(self.pfig_list,llegend)
		
		return

	def _draw_profiles_2d(self,update=False):

		if not self.note_ps['fig'].get() == 1: return
		self.figure['name']['homecursor'].visible = False	
		self._clear_home_canvas()

		self.figure['name']['home'].canvas.draw_idle()
		lent = len(self.profiles['times']); time = self.profiles['times'][0]; lenr = len(self.profiles['fit']['ne'][time][0]);
		xx = np.copy(self.profiles['fit']['ne'][time][0]);
		if self.note_ps['xmap'].get() == 2: lenr = 100; xx = np.linspace(1.7,2.3,lenr)
		ne_prof2 = np.zeros((lent,lenr))

		self.figure['name']['homecursor'].visible = False	
		self.figure['axes']['home'] = [self.figure['name']['home'].add_subplot(1,1,1)]


		for i in range(lent):
			time = self.profiles['times'][i]
			vals = self.profiles['fit']['ne'][time]
			if self.note_ps['xmap'].get() == 1:
				nef = interp1d(vals[3](vals[0]),vals[1])
				for j in range(len(vals[0])):
					if vals[0][j] < vals[3](vals[0])[-1]: ne_prof2[i,j] = nef(vals[0][j])
					else: ne_prof2[i,j] = vals[1][-1]
	
			elif self.note_ps['xmap'].get() == 2:
				xx2 = vals[4](xx,0.);
				for j in range(len(xx2)): xx2[j] = min(xx2[j],1.1)
				nef = interp1d(vals[0],vals[1]); ne_prof2[i,:] = nef(xx2)
			else: ne_prof2[i,:] = vals[1]

		xlabel = '$\\psi_N [a.u]$';
		if self.note_ps['xmap'].get() == 1: xlabel = '$\\rho_N [a.u]$';
		elif self.note_ps['xmap'].get() == 2: xlabel = '$R [m]$';

		A = self.figure['axes']['home'][0].contourf(self.profiles['times']/1.e3,xx,np.transpose(ne_prof2))
		self.figure['name']['home'].colorbar(A,orientation='vertical')
		self.figure['axes']['home'][0].set_xlabel('Time [s]')
		self.figure['axes']['home'][0].set_ylabel(xlabel)

		self.figure['name']['home'].tight_layout()

		return

	def _draw_channel_gap(self,diff=False):

		self._clear_home_canvas()
		self.figure['name']['homecursor'].visible = False
		core_ind = []; edge_ind = [];
		for i in range(self.ts['core']['nch']):
			if self.note_ps['ts_core%i'%(i+1)].get() == 1: core_ind.append(i)
		for i in range(self.ts['edge']['nch']):
			if self.note_ps['ts_edge%i'%(i+1)].get() == 1: edge_ind.append(i)			

		count = len(core_ind)+len(edge_ind)
		inds = []; llegend = [];
		for i in core_ind: 
			inds.append(i); 
			if diff: llegend.append('C%02i'%(i+1))
			else: llegend.append('C%02i-FIT'%(i+1)); llegend.append('C%02i-RAW'%(i+1))
		for i in edge_ind: 
			inds.append(i+self.ts['core']['nch']);
			if diff: llegend.append('E%02i'%(i+1))
			else: llegend.append('E%02i-FIT'%(i+1)); llegend.append('E%02i-RAW'%(i+1))

		if (count) > 6:
			print(">>> Max CH.# is 6"); return
		if count == 0: print('>>> Empty CH. selection'); return

		gs = gridspec.GridSpec(count, 1,hspace=0.)
		self.figure['axes']['home'] = [self.figure['name']['home'].add_subplot(gs[0,0])]
		time = self.profiles['times'][0]; lent = len(self.profiles['times'])

		vals = self.profiles['fit']['ne'][time]
		nch = self.ts['core']['nch'] + self.ts['edge']['nch']
		xxr = vals[2]['xxr'][0:nch]; raw = vals[2]['raw'][0:nch]; raws=vals[2]['raws'][0:nch]; vdiff = np.zeros((lent,nch)); vval = np.copy(vdiff); sig = np.copy(vdiff)
		for j in range(len(xxr)): xxr[j] = min(xxr[j],1.1) 
		nf = interp1d(vals[0],vals[1]); 
		raw = self._do_scale(time,raw)
		vdiff[0,:] = nf(xxr) - raw; vval[0,:] = raw; sig[0,:] = raws
		
		for i in range(1,lent):
			time = self.profiles['times'][i]
			vals = self.profiles['fit']['ne'][time];
			raw = vals[2]['raw'][0:nch]; raws=vals[2]['raws'][0:nch]; 
			xxr = np.copy(vals[2]['xxr'][0:nch]); nf = interp1d(vals[0],vals[1])
			for j in range(len(xxr)): xxr[j] = min(xxr[j],1.1) 
			raw = self._do_scale(time,raw)
			vdiff[i,:] = nf(xxr) - raw;  vval[i,:] = raw; sig[i,:] = raws

		if diff: 
			self.figure['axes']['home'][0].errorbar(self.profiles['times'],vdiff[:,inds[0]],yerr=sig[:,inds[0]],fmt='x--')
			self.figure['axes']['home'][0].legend([llegend[0]])
			self.figure['axes']['home'][0].axhline(y=0,linestyle='--',c='gold')				
		else:
			self.figure['axes']['home'][0].plot(self.profiles['times'],vval[:,inds[0]]+vdiff[:,inds[0]],linestyle='--',marker='x')
			self.figure['axes']['home'][0].errorbar(self.profiles['times'],vval[:,inds[0]],yerr=sig[:,inds[0]],fmt='x--')
			self.figure['axes']['home'][0].legend([llegend[0],llegend[1]])

		for k in range(1,count):
			self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(gs[k,0],sharex=self.figure['axes']['home'][0]))
			if diff: 
				self.figure['axes']['home'][k].errorbar(self.profiles['times'],vdiff[:,inds[k]],yerr=sig[:,inds[k]],fmt='x--')
				self.figure['axes']['home'][k].legend([llegend[k]])
				self.figure['axes']['home'][k].axhline(y=0,linestyle='--',c='gold')				
			else:
				self.figure['axes']['home'][k].plot(self.profiles['times'],vval[:,inds[k]]+vdiff[:,inds[k]],linestyle='--',marker='x')
				self.figure['axes']['home'][k].errorbar(self.profiles['times'],vval[:,inds[k]],yerr=sig[:,inds[0]],fmt='x--')
				self.figure['axes']['home'][k].legend([llegend[2*k],llegend[2*k+1]])

		self.figure['axes']['home'][-1].set_xlabel('time [ms]')

		for k in range(count-1): self.figure['axes']['home'][k].axes.get_xaxis().set_visible(False)

		return

	def _draw_channel_gap2d(self):

		self._clear_home_canvas()
		self.figure['name']['homecursor'].visible = False	

		gs = gridspec.GridSpec(1, 2); tcount = 0
		self.figure['axes']['home'] = []
		for i in range(2):
				self.figure['axes']['home'].append(self.figure['name']['home'].add_subplot(gs[0,i]))

		time = self.profiles['times'][0]; lent = len(self.profiles['times']);

		vals = self.profiles['fit']['ne'][time]
		nch = self.ts['core']['nch'] + self.ts['edge']['nch']
		xxr = vals[2]['xxr'][0:nch]; raw = vals[2]['raw'][0:nch]; raws=vals[2]['raws'][0:nch]; 
		vval = np.zeros((lent,nch)); vraw = np.copy(vval); vraws = np.copy(vval);
		for j in range(len(xxr)): xxr[j] = min(xxr[j],1.1) 
		nf = interp1d(vals[0],vals[1]); sigs = np.zeros(nch)
		for k in range(nch): sigs[k] = raws[k]/max(raw[k],0.3);
		raw = self._do_scale(time,raw)
		raws= self._do_scale(time,raws)
		vval[0,:] = nf(xxr); vraw[0,:] = raw; vraws[0,:] = raws
		for k in range(nch): vraws[0,k] = np.max([vraws[0,k],vraw[0,k]*0.02,0.02])
	
		for i in range(1,lent):
			time = self.profiles['times'][i]
			vals = self.profiles['fit']['ne'][time];
			raw = vals[2]['raw'][0:nch]; raws=vals[2]['raws'][0:nch]; 
			xxr = np.copy(vals[2]['xxr'][0:nch]); nf = interp1d(vals[0],vals[1])
			for k in range(nch): sigs[k] = sigs[k] + raws[k]/max(raw[k],0.3);
			for j in range(len(xxr)): xxr[j] = min(xxr[j],1.1)
			raw = self._do_scale(time,raw)
			raws= self._do_scale(time,raws)
			vval[i,:] = nf(xxr);  vraw[i,:] = raw; vraws[i,:] = raws;
			for k in range(nch): vraws[i,k] = np.max([vraws[i,k],vraw[i,k]*0.02,0.02])

		xmin = float(self.note_ps['CORE_MIN'].get()); xmax = float(self.note_ps['CORE_MAX'].get()); xn = int(self.note_ps['CORE_N'].get())
		ymin = float(self.note_ps['EDGE_MIN'].get()); ymax = float(self.note_ps['EDGE_MAX'].get()); yn = int(self.note_ps['EDGE_N'].get())
		xlim = np.linspace(xmin,xmax,xn); ylim = np.linspace(ymin,ymax,yn); ccore = np.linspace(1,self.ts['core']['nch']+1,self.ts['core']['nch']);
		ecore = np.linspace(1,self.ts['edge']['nch']+1,self.ts['edge']['nch']); 

		diffc = np.zeros((xn,self.ts['core']['nch']+1)); diffe = np.zeros((yn,self.ts['edge']['nch']+1))
		sigs = sigs / lent

		for k in range(self.ts['core']['nch']):
			for j in range(xn):
				factor = xlim[j]
				tempv = (vval[:,k] - vraw[:,k]*factor)/vraws[:,k]/factor
				diffc[j,k] = np.mean(tempv)
			diffc[:,k] = diffc[:,k] / max(abs(diffc[:,k]))

		for k in range(self.ts['edge']['nch']):
			for j in range(yn):
				factor = ylim[j]
				tempv = (vval[:,k+self.ts['core']['nch']] - vraw[:,k+self.ts['core']['nch']]*factor)/vraws[:,k+self.ts['core']['nch']]/factor
				diffe[j,k] = np.mean(tempv)
			diffe[:,k] = diffe[:,k] / max(abs(diffe[:,k]))

		ecore2 = np.linspace(1,self.ts['edge']['nch']+1,self.ts['edge']['nch']+1);
		ccore2 = np.linspace(1,self.ts['core']['nch']+1,self.ts['core']['nch']+1);

		self.figure['axes']['home'][0].pcolor(xlim,ccore2-0.5,-np.transpose(abs(diffc)));
		for k in range(1,self.ts['core']['nch']+1): 
			self.figure['axes']['home'][0].text(1.,k,'C%02i'%k,color='r')
			fc = sigs[k-1]; xx = [max(xmin,1-fc),min(xmax,1+fc)]; yy = [k,k];
			self.figure['axes']['home'][0].plot(xx,yy,color='magenta')

		self.figure['axes']['home'][0].axvline(x=1.,linestyle='--',c='gold')
		self.figure['axes']['home'][1].pcolor(ylim,ecore2-0.5,-np.transpose(abs(diffe)));
		for k in range(1,self.ts['edge']['nch']+1): 
			self.figure['axes']['home'][1].text(1.,k,'E%02i'%k,color='r')
			fc = sigs[k-1+self.ts['core']['nch']]; xx = [max(ymin,1-fc),min(ymax,1+fc)]; yy = [k,k];
			self.figure['axes']['home'][1].plot(xx,yy,color='magenta')

		self.figure['axes']['home'][1].axvline(x=1.,linestyle='--',c='gold')

		self.figure['axes']['home'][0].set_xlabel('MULTI. SCALE')
		self.figure['axes']['home'][1].set_xlabel('MULTI. SCALE')

		self.figure['axes']['home'][0].set_ylabel('TS CORE CH [#]')
		self.figure['axes']['home'][1].set_ylabel('TS EDGE CH [#]')

		self.figure['name']['home'].tight_layout()

		return

	def _do_scale(self,ttime,raw):
	
		core_ind2 = int(self.note_in['CALCHC'].get())-1; 
		edge_ind2 = int(self.note_in['CALCHE'].get())-1 + self.ts['core']['nch']; 
		ratio2 = float(self.note_in['CALCH'].get())
		factor1c = self.profiles['scale']['core'];  factor1e = self.profiles['scale']['edge']
		factor2c = float(self.note_in['CM'].get()); factor2e = float(self.note_in['EM'].get());

		if self.note_in['docal'].get() == 1:
			if self.note_in['fcore'].get() == 1: factor2e = factor2c
			else: factor2c = factor2e

		if self.profiles['scale']['fcore']: factor1e = factor1e * self.profiles['scale'][ttime]
		else: factor1c = factor1c / self.profiles['scale'][ttime]

		if self.note_in['docal'].get() == 1:
			factor2 = raw[core_ind2] / raw[edge_ind2]/ratio2 * (factor1e/factor1c)
			if self.note_in['fcore'].get() == 1: factor2e = factor2e * factor2
			else: factor2c = factor2c / factor2

		raws = np.copy(raw)
		for k in range(self.ts['core']['nch']):
			raws[k] = raws[k] * factor2c / factor1c;
		for k in range(self.ts['core']['nch'],self.ts['core']['nch']+self.ts['edge']['nch']):
			raws[k] = raws[k] * factor2e / factor1e;

		return (raws)

	def _draw_channel_gap2dt(self, auto=False):

		time = self.profiles['times'][0]; lent = len(self.profiles['times'])

		vals = self.profiles['fit']['ne'][time]
		nch = self.ts['core']['nch'] + self.ts['edge']['nch']
		xxr = vals[2]['xxr'][0:nch]; raw = vals[2]['raw'][0:nch]; raws=vals[2]['raws'][0:nch]; 
		vval = np.zeros((lent,nch)); vraw = np.copy(vval); vraws = np.copy(vval);
		for j in range(len(xxr)): xxr[j] = min(xxr[j],1.1) 
		nf = interp1d(vals[0],vals[1]); 
		raw = self._do_scale(time,raw)
		raws= self._do_scale(time,raws);
		vval[0,:] = nf(xxr); vraw[0,:] = raw; vraws[0,:] = raws;
		for k in range(nch): vraws[0,k] = np.max([vraws[0,k],vraw[0,k]*0.02,0.02])
		
		for i in range(1,lent):
			time = self.profiles['times'][i]
			vals = self.profiles['fit']['ne'][time];
			raw = vals[2]['raw'][0:nch]; raws=vals[2]['raws'][0:nch]; 
			xxr = np.copy(vals[2]['xxr'][0:nch]); nf = interp1d(vals[0],vals[1])
			for j in range(len(xxr)): xxr[j] = min(xxr[j],1.1) 
			raw = self._do_scale(time,raw)
			raws= self._do_scale(time,raws);
			vval[i,:] = nf(xxr);  vraw[i,:] = raw; vraws[i,:] = raws;
			for k in range(nch): vraws[i,k] = np.max([vraws[i,k],vraw[i,k]*0.02,0.02])

		xmin = float(self.note_ps['CORE_MIN'].get()); xmax = float(self.note_ps['CORE_MAX'].get()); xn = int(self.note_ps['CORE_N'].get())
		ymin = float(self.note_ps['EDGE_MIN'].get()); ymax = float(self.note_ps['EDGE_MAX'].get()); yn = int(self.note_ps['EDGE_N'].get())
		xlim = np.linspace(xmin,xmax,xn); ylim = np.linspace(ymin,ymax,yn); ccore = np.linspace(1,self.ts['core']['nch']+1,self.ts['core']['nch']);
		ecore = np.linspace(1,self.ts['edge']['nch']+1,self.ts['edge']['nch']);
		diffc = np.zeros((xn,self.ts['core']['nch'])); diffe = np.zeros((yn,self.ts['edge']['nch']))
		corem = float(self.note_ps['COREM'].get()); edgem = float(self.note_ps['EDGEM'].get());

		for k in range(self.ts['core']['nch']):
			for j in range(xn):
				factor = xlim[j]
				tempv = (vval[:,k] - vraw[:,k]*factor)/vraws[:,k]/factor
				diffc[j,k] = np.mean(abs(tempv))

		for k in range(self.ts['edge']['nch']):
			for j in range(yn):
				factor = ylim[j]
				tempv = (vval[:,k+self.ts['core']['nch']] - vraw[:,k+self.ts['core']['nch']]*factor)/vraws[:,k+self.ts['core']['nch']]/factor
				diffe[j,k] = np.mean(abs(tempv))
	
		difft = np.zeros((xn+1,yn+1))
		for i in range(xn):
			for j in range(yn):
				for k in range(1,self.ts['core']['nch']+1):
					if self.note_in['ts_core%i'%k].get() == 1:
						difft[i,j] = difft[i,j] + diffc[i,k-1]*corem
				for k in range(1,self.ts['edge']['nch']+1):				
					if self.note_in['ts_edge%i'%k].get() == 1:
						difft[i,j] = difft[i,j] + diffe[j,k-1]*edgem

		if auto: return (xlim,ylim,difft)

		xlim = np.append(xlim,xlim[-1]+0.01); xlim = xlim - 0.005
		ylim = np.append(ylim,ylim[-1]+0.01); ylim = ylim - 0.005

		
		if self.nogui: return

		self._clear_home_canvas((not auto))
		self.figure['name']['homecursor'].visible = False
		gs = gridspec.GridSpec(1, 1); tcount = 0
		self.figure['axes']['home'] = [self.figure['name']['home'].add_subplot(gs[0,0])]
		self.figure['axes']['home'][0].pcolor(xlim,ylim,-np.transpose(abs(difft)));
		self.figure['axes']['home'][0].axvline(x=1.,c='gold',linestyle='--')
		self.figure['axes']['home'][0].axhline(y=1.,c='gold',linestyle='--')
		self.figure['axes']['home'][0].set_xlabel('CORE MULTI')
		self.figure['axes']['home'][0].set_ylabel('EDGE MULTI')
		self.figure['name']['home'].tight_layout()

		return

	def _load_diag(self):

		tmin = 0.; tmax = 500000;
		print('>>> Load TS CORE...')
		g = mds('kstar',self.shotn)
		#Load TS_CORE
		radius = np.array([])
		self.ts['core']['rr'] = radius
		self.ts['edge']['rr'] = radius
		for ch in range(30):
			node_name1 = '\\TS_CORE%i.CORE%i_NE'%(ch+1,ch+1)
			node_name2 = '\\TS_CORE%i.CORE%i_NERRH'%(ch+1,ch+1)
			node_name3 = '\\TS_CORE%i.CORE%i_POS'%(ch+1,ch+1)
			self.ts['core'][ch] = g.get(node_name1)
			self.ts['core_err'][ch] = g.get(node_name2)
			if len(self.ts['core'][ch][0])==0: break			
			if len(self.ts['core_err'][ch][0])==0:self.ts['core_err'][ch] = [self.ts['core'][ch][0],self.ts['core'][ch][0]/10.]
			self.ts['core']['nch'] = ch + 1
			rr = g.get(node_name3)
			radius = np.append(radius,rr[1]/1000.)
		if self.ts['core']['nch'] > 0: 
			self.ts['core']['rr'] = np.copy(radius)
			tmin = max(tmin,self.ts['core'][0][0][0])
			tmax = min(tmax,self.ts['core'][0][0][-1])

		print('>>> Load TS EDGE...')
		radius = np.array([])
		for ch in range(30):
			node_name1 = '\\TS_EDGE%i.EDGE%i_NE'%(ch+1,ch+1)
			node_name2 = '\\TS_EDGE%i.EDGE%i_NERRH'%(ch+1,ch+1)
			node_name3 = '\\TS_EDGE%i.EDGE%i_POS'%(ch+1,ch+1)
			self.ts['edge'][ch] = g.get(node_name1)
			self.ts['edge_err'][ch] = g.get(node_name2)
			if len(self.ts['edge'][ch][0])==0: break			
			if len(self.ts['edge_err'][ch][0])==0:self.ts['edge_err'][ch] = [self.ts['edge'][ch][0],self.ts['edge'][ch][0]/10.]
			self.ts['edge']['nch'] = ch + 1
			rr = g.get(node_name3)
			radius = np.append(radius,rr[1]/1000.)
		if self.ts['edge']['nch'] > 0: 
			self.ts['edge']['rr'] = np.copy(radius)
			tmin = max(tmin,self.ts['edge'][0][0][0])
			tmax = min(tmax,self.ts['edge'][0][0][-1])
			self.ts['edge']['rr'][0] = self.ts['edge']['rr'][0] + 0.002

		if self.ts['edge']['nch'] == 0: print('>>> TS data unavail.');
		else: print('>>> TS CORE/EDGE CH. %i/%i [#] '%(self.ts['core']['nch'],self.ts['edge']['nch']))
		if self.ts['edge']['nch'] > 0: self.profiles['diag']['ts'] = True

		print('>>> Load INT...')
		line = '>>> INT CH '
		length = [1.9,2.75]
		for ch in range(self.int['nch']):
			node_name1 = '\\NE_INT%02i'%(ch+1)
			self.int[ch+1]['val'] = g.get(node_name1)
			if len(self.int[ch+1]['val'][0]) > 0: 
				self.int[ch+1]['is'] = True
				dt = (self.int[ch+1]['val'][0][10]-self.int[ch+1]['val'][0][9])*1000
				dscale = int(2./dt);
				self.int[ch+1]['val'] = [self.int[ch+1]['val'][0][0::dscale],self.int[ch+1]['val'][1][0::dscale]/length[ch]]
				line = line + '%i '%(ch+1)
				tmin = max(tmin,self.int[ch+1]['val'][0][0])
				tmax = min(tmax,self.int[ch+1]['val'][0][-1])
				self.profiles['diag']['int%i'%(ch+1)] = True

		if not line=='>>> INT CH ': print(line+ 'avail.')
		else: print(line+ 'unavail.')
		print('>>> Load TCI...')
		line = '>>> TCI CH '
		for ch in range(self.tci['nch']):
			node_name1 = '\\NE_TCI%02i'%(ch+1)
			self.tci[ch+1]['val'] = g.get(node_name1)
			if len(self.tci[ch+1]['val'][0]) > 0: 
				self.tci[ch+1]['is'] = True
				dt = (self.tci[ch+1]['val'][0][10]-self.tci[ch+1]['val'][0][9])*1000
				dscale = int(2./dt);
				self.tci[ch+1]['val'] = [self.tci[ch+1]['val'][0][0::dscale],self.tci[ch+1]['val'][1][0::dscale]]
				line = line + '%i '%(ch+1)
				tmin = max(tmin,self.tci[ch+1]['val'][0][0])
				tmax = min(tmax,self.tci[ch+1]['val'][0][-1])
				self.profiles['diag']['tci%i'%(ch+1)] = True

		if not line=='>>> TCI CH.': print(line+ 'avail.')
		else: print(line+ 'unavail.')
		print('>>> Load Ip/Da...')
		node_name1 = '\\pcrc03'
		node_name2 = '\\tor_ha10'
		self.mds['ip'] = g.get(node_name1)
		self.mds['da'] = g.get(node_name2)

		self.ts_core_ch = []; self.ts_edge_ch = [];
		for ch in range(5): self.ts_core_ch.append('%02i'%(self.ts['core']['nch']-ch))
		for ch in range(5): self.ts_edge_ch.append('%02i'%(ch+1))		

		lena = len(self.mds['da'][0]); lenv = len(self.mds['da'][1]); lend = lena- lenv
		if not lend == 0:
			if lend > 0: 
				for i in range(lend): self.mds['da'][1].append(0)
			else:
				for i in range(lend): self.mds['da'][0].append(0)

		self.tmin = tmin*1000.; self.tmax = tmax*1000.;
		print('>>> Load EFIT...')
		ok,year = get_year(self.shotn)
		if not ok==1: print('>>> Invalid shot & year...'); exit()
		efitok = False; self.out_rt1=None
		if self.note_in['gflag'].get() > 0:
			for i in range(1,6): 
				if (self.efit_list['isefit'][i] and (not efitok)): efitok = True
		else:
			currdir = os.getcwd()
			if not os.path.isdir('DENA/RT1'): os.mkdir('DENA/RT1')
			os.chdir('DENA/RT1')
			try: rtefit.run(self.shotn,0.,10000.); print('>>> RTEFIT avail.'); efitok=True
			except: print('>>> RTEFIT unavail.');
			self.out_rt1 = list()
			status, output = subprocess.getstatusoutput('ls')
			output = output.split()
			for g in output:
				if g.find('%s_'%self.shotn) > -1:
					tt = g.split('.')[0].split('_')[-1]
					self.out_rt1.append(tt)
			os.chdir(currdir)
		if not efitok: exit()

		return

	def _preset_tci_opt(self):
	
		if self.note_in['gflag'].get() == 0: shot_end = float(self.out_rt1[-1])/1000.; print('>>> Ends at %4.1fs'%shot_end);
		else: shot_end = float(self.efit_list['times'][self.note_in['gflag'].get()][-1])/1000.; print('>>> Ends at %4.1fs'%shot_end);
		if shot_end<1.: print('>>> No successful shot'); exit()
		print('>>> Shot #: %i'%self.shotn)
		for ch in range(self.tci['nch']):
			dat  = self.tci[ch+1]['val']
			tind = np.where(abs(dat[0]-shot_end-1.)<.5)
			try: davg = np.mean(dat[1][tind])/(np.max(dat[0][tind])-np.min(dat[0][tind]))
			except: davg = 0
			isdrift = False
			if abs(davg)>1.0: isdrift = True
			print('>>> TCI%02i end value %3.1f/Drift=%s'%(ch+1,davg,isdrift))
			if isdrift:
				self.note_in['tci%i'%(ch+1)].set(0)
		if self.brutal: 
			self.note_in['WTCI%02i'%(5)].set(0.05)
		
	def _make_ts_raw_files(self):

		scale = 1.e-18;
		print('>>> Generate TS inputs...')
		if self.note_in['docal'].get() == 1: self.profiles['scale']['docal'] = True
		else: self.profiles['scale']['docal'] = False

		if self.note_in['fcore'].get() == 1: self.profiles['scale']['fcore'] = True
		else: self.profiles['scale']['fcore'] = False		

		tot_nch = self.ts['core']['nch']+self.ts['edge']['nch']
		profv = dict(); profe = dict(); profr = np.zeros(tot_nch)
		profr[0:self.ts['core']['nch']] = self.ts['core']['rr']; profr[self.ts['core']['nch']:tot_nch] = self.ts['edge']['rr'];
		for k in range(len(self.profiles['times'])):
			time1 = self.profiles['times'][k]
			if not self.ts['core']['nch'] > 0:
				self.profiles['infile']['ts'][time1] = dummy_dir
				continue
			
			ltime = time1 - float(self.note_in['ATS'].get()) * 0.5
			utime = time1 + float(self.note_in['ATS'].get()) * 0.5
			core_ind1 = np.where((self.ts['core'][0][0]*1000 >= ltime))
			core_ind2 = np.where((self.ts['core'][0][0][core_ind1]*1000 <= utime))
			edge_ind1 = np.where((self.ts['edge'][0][0]*1000 >= ltime))
			edge_ind2 = np.where((self.ts['edge'][0][0][edge_ind1]*1000 <= utime))
			
			for i in range(tot_nch):
				if i < self.ts['core']['nch']: label = 'core'; j = i; ind1 = core_ind1; ind2 = core_ind2
				else: label = 'edge'; j = i - self.ts['core']['nch']; ind2 = edge_ind1; ind2 = edge_ind2

				profv[i] = np.nan_to_num(self.ts[label][j][1][ind1][ind2]) * scale
				profe[i] = np.nan_to_num(self.ts[label+'_err'][j][1][ind1][ind2]) * scale

			core_ind = int(self.note_in['CALCHC'].get())-1; edge_ind = int(self.note_in['CALCHE'].get())-1; ratio = float(self.note_in['CALCH'].get())
			factor = np.mean(profv[core_ind]) / np.mean(profv[edge_ind+self.ts['core']['nch']])/ratio
			if not self.note_in['docal'].get(): factor = 1.

			self.profiles['scale'][time1] = factor

			if self.note_in['fcore'].get() == 1:
				for i in range(self.ts['core']['nch'],tot_nch): profv[i] = profv[i]*factor; profe[i] = profe[i]*factor;
			else:
				for i in range(self.ts['core']['nch']): profv[i] = profv[i]/factor; profe[i] = profe[i]/factor;

			self.profiles['infile']['ts'][time1] = 'DENA/TS_NE_%ims.dat'%time1
			f = open(self.profiles['infile']['ts'][time1],'w')
			f.write('R[m]\tZ[m]\tVAL\tERR\n')
			for i in range(len(profv[0])):
				for j in range(tot_nch):
					f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(profr[j],0.,profv[j][i],profe[j][i]))
			f.close()

		return

	def _make_tci_raw_files(self):

		L_weight=[4,4,7.23,5.51,4.60,3.81,2.48]
		L_weight=[1,1,1,1,1,1,1]
		print('>>> Generate TCI/INT inputs...')

		for k in range(len(self.profiles['times'])):

			time1 = self.profiles['times'][k]
			iltime = time1 - float(self.note_in['AINT'].get()) * 0.5
			iutime = time1 + float(self.note_in['AINT'].get()) * 0.5			
			tltime = time1 - float(self.note_in['ATCI'].get()) * 0.5
			tutime = time1 + float(self.note_in['ATCI'].get()) * 0.5

			tot_nch = self.int['nch']+self.tci['nch']
			self.profiles['infile']['tciv'][time1] = np.zeros(tot_nch)
			self.profiles['infile']['tcis'][time1] = np.zeros(tot_nch)
			
			for i in range(tot_nch):
				if i<self.int['nch']:
					j = i+1; ltime = iltime; utime = iutime
					if self.int[j]['is']: dat = self.int[j]['val']
					else: dat = None
				else:
					j = i + 1 - self.int['nch']; ltime = tltime; utime = tutime;
					if self.tci[j]['is']: dat = self.tci[j]['val']
					else: dat = None

				if not dat == None:
					lenx = len(dat[0]); lenv = len(dat[1]); lenm = min(lenx,lenv);
					dat[0] = dat[0][0:lenm]; dat[1] = dat[1][0:lenm]
					ind1 = np.where(dat[0]*1000>ltime)
					ind2 = np.where(dat[0][ind1]*1000<utime)

					tciavg = np.sum(dat[1][ind1][ind2])/len(dat[1][ind1][ind2])
					ind3 = np.where(abs(dat[1][ind1][ind2]-tciavg) < 0.5)
					if len(dat[1][ind1][ind2][ind3])==0: continue

					tciavg = np.sum(dat[1][ind1][ind2][ind3])/len(dat[1][ind1][ind2][ind3])
				#	tcisig = (np.max(dat[1][ind1][ind2][ind3]) - np.min(dat[1][ind1][ind2][ind3])) / 3.
					tcisig = np.std(dat[1][ind1][ind2][ind3])* 1.5 / L_weight[i] #* 4.

					if tciavg<=0.: tciavg = 0.; tcisig = 0.
					self.profiles['infile']['tciv'][time1][i] = tciavg
					self.profiles['infile']['tcis'][time1][i] = tcisig


			self.profiles['infile']['tcif'][time1] = 'DENA/TCI_%ims.dat'%time1
			f = open(self.profiles['infile']['tcif'][time1],'w')
			for k in range(tot_nch):
				if self.profiles['infile']['tciv'][time1][k] > 0: isdata = 1
				else: isdata = 0
				val = self.profiles['infile']['tciv'][time1][k]; sig = self.profiles['infile']['tcis'][time1][k];
				f.write('%i %4.2f %4.3f \n'%(isdata,max(0.,val),max(0,sig)))
			f.close()
		return

	def _get_gfiles(self):

		print('>>> Download GFILES')
		gflag = ['RT1','EFIT01','EFIT02','EFIT03','EFIT04','EFIT05']
		tlen = len(self.profiles['times'])
		for flag in gflag:
			if not os.path.isdir('DENA/%s'%flag): os.mkdir('DENA/%s'%flag);

		efit_dir = 'DENA/%s'%gflag[self.note_in['gflag'].get()]
		shot = self.shotn

		if self.note_in['gflag'].get() > 0: self.profiles['gflag'] = 'EFIT%02i'%self.note_in['gflag'].get()
		else: self.profiles['gflag'] = 'EFITRT1'

		if self.note_in['gflag'].get() > 0:	
			efitdir = self.efit_list['dirs'][self.note_in['gflag'].get()]
			fileline = ''
			for i in range(tlen):
				ttime = self.profiles['times'][i]
				self.profiles['infile']['gfile'][ttime] = '%s/g%06i.%06i'%(efit_dir,shot,ttime)
				if not os.path.isfile(self.profiles['infile']['gfile'][ttime]):
					fileline = fileline + ' %sg%06i.%06i'%(efitdir,shot,ttime)
	
			pwd = os.getcwd(); os.chdir(efit_dir); 
			if not fileline == '':
#			        comm = 'scp %s:"%s" .'%(efit_address,fileline)
			        comm = 'cp %s .'%(fileline)
			        os.system(comm)
			os.chdir(pwd)
		else: 
			for i in range(tlen):
				ttime = self.profiles['times'][i]
				self.profiles['infile']['gfile'][ttime] = '%s/kstar_EFITRT1_%i_%06i.geqdsk'%(efit_dir,shot,ttime)

		return

	def _set_fitopt(self):

		use_rho = False; excn = '';
		if self.note_in['rhofit'].get() == 1: use_rho = True
		
		for k in range(self.ts['core']['nch']):
			if self.note_in['ts_core%i'%(k+1)].get() == 0:
				if excn == '': excn = '%i'%k
				else: excn = excn + ',%i'%k
		for k in range(self.ts['edge']['nch']):
			if self.note_in['ts_edge%i'%(k+1)].get() == 0:
				if excn == '': excn = '%i'%(k+self.ts['core']['nch'])
				else: excn = excn + ',%i'%(k+self.ts['core']['nch'])

		print('>>> Exclude %s'%excn)

		self.fit.fit_opt['use_rho']['ne'] = use_rho
		self.fit.fit_opt['exclude']['ne']['ts'] = excn
		self.fit.fit_opt['weight']['ne']['ts'] = float(self.note_in['WTS'].get())
		if not self.ts['core']['nch'] > 0: self.fit.fit_opt['weight']['ne']['ts'] = 0.
		self.fit.fit_opt['weight']['ne']['refl'] = float(self.note_in['WREFL'].get())
		if not self.ref['is']: self.fit.fit_opt['weight']['ne']['refl'] = 0.

		self.fit.noprint = True
		self.fit.fit_opt['oli']['ne']['ts']['use'] = True
		self.fit.fit_opt['oli']['ne']['ts']['per'] = 0.9
		self.fit.fit_opt['psi_end']['ne']['ts'] = 1.01
		self.fit.fit_opt['avg']['ne']['ts'] = False

		min_ped = float(self.note_in['MINPH'].get())
		wid_ped = float(self.note_in['WIDTH'].get())

		if self.note_in['hmode'].get() == 1:
			self.fit.fit_opt['sep_fix']['ne'] = True
#			self.fit.param['ne']['min'][0] = 0.1
#			self.fit.param['ne']['val'][0] = max(min_ped * 0.4,0.4)
#			self.fit.param['ne']['max'][0] = 1.2
			self.fit.fit_opt['sep_val']['ne'] = max(min_ped * 0.3,0.1)
			self.fit.fit_opt['width_fix']['ne'] = False
			self.fit.fit_opt['width_val']['ne'] = wid_ped
			self.fit.param['ne']['min'][2] = wid_ped*0.5
			self.fit.param['ne']['val'][2] = wid_ped
			self.fit.param['ne']['max'][2] = wid_ped*1.5
			self.fit.param['ne']['vary'][2] = True
			self.fit.param['ne']['min'][4] = 1.1
			self.fit.param['ne']['max'][4] = 2.7 #2.2
			self.fit.param['ne']['min'][5] = 1.1
			self.fit.param['ne']['max'][5] = 2.7 #2.2
			self.fit.param['ne']['min'][1] = max((min_ped- self.fit.fit_opt['sep_val']['ne']) / 2. / np.tanh(1),0.1)
		else:
#			self.fit.fit_opt['sep_fix']['ne'] = True
			self.fit.fit_opt['sep_val']['ne'] = max(min_ped * 0.3,0.1)
			self.fit.param['ne']['min'][2] = wid_ped*0.9
			self.fit.param['ne']['val'][2] = wid_ped
			self.fit.param['ne']['max'][2] = wid_ped*1.1
			self.fit.fit_opt['width_fix']['ne'] = True
			self.fit.fit_opt['width_val']['ne'] = max(wid_ped,0.1)
			self.fit.param['ne']['min'][4] = 1.05
			self.fit.param['ne']['max'][4] = 2.6 #3.0
			self.fit.param['ne']['min'][5] = 1.05
			self.fit.param['ne']['max'][5] = 2.6 #3.0
			self.fit.param['ne']['min'][1] = 0.0
#			self.fit.param['ne']['max'][1] = 2.0
#			self.fit.param['ne']['val'][1] = min_ped

		self.fit.fit_opt['scale']['ne']['ts']['core'] = float(self.note_in['CM'].get())
		self.fit.fit_opt['scale']['ne']['ts']['edge'] = float(self.note_in['EM'].get())

		if self.note_in['docal'].get() == 1:
			if self.note_in['fcore'].get() == 1: self.fit.fit_opt['scale']['ne']['ts']['edge'] = float(self.note_in['CM'].get())
			else: self.fit.fit_opt['scale']['ne']['ts']['core'] = float(self.note_in['EM'].get())

		self.profiles['scale']['core'] = self.fit.fit_opt['scale']['ne']['ts']['core']
		self.profiles['scale']['edge'] = self.fit.fit_opt['scale']['ne']['ts']['edge']

		return

	def _dofit(self,tind):

		time1 = self.profiles['times'][tind]
		self.fit.ne_prof['fit_old'] = [];
		self.fit.fit_opt['file']['gfile']     = None
		self.notci = True
		for k in range(len(self.fit.inter_list)):
			flag = self.fit.inter_list[k]
			self.fit.fit_opt[flag]['val'] =self.profiles['infile']['tciv'][time1][k]

			if k<2: insig = float(self.note_in['WINT%02i'%(k+1)].get()); usech = self.note_in['int%i'%(k+1)].get();
			else: insig = float(self.note_in['WTCI%02i'%(k-1)].get()); usech = self.note_in['tci%i'%(k-1)].get();
			if insig == 0.: insig = self.profiles['infile']['tcis'][time1][k]
			else: self.notci = False
			self.fit.fit_opt[flag]['sig'] = insig;
			if usech == 0: self.fit.fit_opt[flag]['sig'] = 0.
			else: self.notci = False

		if self.notci: self.fit.fit_opt['weight']['ne']['ts'] = 1.

		self.fit.fit_opt['file']['ne']['kfile'] = None
		for j in self.fit.__dict__['ne_list']: self.fit.fit_opt['file']['ne'][j] = None		

		self.fit.fit_opt['file']['gfile'] = self.profiles['infile']['gfile'][time1]

		self.fit.fit_opt['file']['ne']['ts'] = self.profiles['infile']['ts'][time1]
		self.fit.post['didfit']['ne'] = False
		self.fit.first_run = True
		self.fit.fit_opt['mds']['time'] = time1
		self.fit.main_run()

		self.fit.post['time'] = time1
		nel = [self.fit.post['int01'],self.fit.post['int02'],self.fit.post['tci01'],self.fit.post['tci02'],self.fit.post['tci03'],self.fit.post['tci04'],self.fit.post['tci05']]
		self.profiles['fit']['nel'][time1] = np.copy(np.array(nel))
		fit_eq = self.fit.fit_eq; fit_ne = self.fit.ne_prof;
		if self.note_in['rhofit'].get() == 1: fit_ne['ts']['xxr'] = fit_eq['rho_to_psi'](fit_ne['ts']['xxr'])
		tsx = np.copy(fit_ne['ts']['xxr']); 
		for i in range(len(tsx)): tsx[i] = min(tsx[i],1.1)
		tsf = interp1d(fit_eq['psin2'],fit_ne['fit2p']); tsv = tsf(tsx);
		self.profiles['fit']['ne'][time1] = copy.deepcopy([fit_eq['psin2'],fit_ne['fit2p'],fit_ne['ts'],fit_eq['psi_to_rho'],fit_eq['psif'],tsv,self.fit.post['chi']['ne']])

		return

	def _auto_scale(self):

		cmint = self.note_ps['CORE_MIN'].get();
		cmaxt = self.note_ps['CORE_MAX'].get()
		emint = self.note_ps['EDGE_MIN'].get()
		emaxt = self.note_ps['EDGE_MAX'].get()

		cw    = self.note_ps['COREM'].get()
		ew    = self.note_ps['EDGEM'].get()

		cn    = self.note_ps['CORE_N'].get()
		en    = self.note_ps['EDGE_N'].get()

		self.note_ps['CORE_MIN'].set('0.1');
		self.note_ps['CORE_MAX'].set('5.1');
		self.note_ps['EDGE_MIN'].set('0.1');
		self.note_ps['EDGE_MAX'].set('5.1');
		self.note_ps['CORE_N'].set('101');
		self.note_ps['EDGE_N'].set('101');
		self.note_in['CM'].set('1.0')
		self.note_in['EM'].set('1.0')

		if self.note_in['fcore'].get() == 1: self.note_ps['EDGEM'].set('0.0')
		else: self.note_ps['COREM'].set('0.0')
		
		xlim,ylim,difft = self._draw_channel_gap2dt(auto=True)
		lenx = len(xlim); leny = len(ylim);

		if self.note_in['fcore'].get() == 1:
			self.note_ps['EDGEM'].set(ew)
			diff = difft[:-2,0]
			ind = np.argmin(diff)
			self.note_in['CM'].set('%2.2f'%xlim[ind])
		else:
			self.note_ps['COREM'].set(cw)
			diff = difft[0,:-2]
			ind = np.argmin(diff)
			self.note_in['EM'].set('%2.2f'%ylim[ind])

		xlim,ylim,difft = self._draw_channel_gap2dt(auto=True)

		if not self.note_in['docal'].get() == 1:
			xlim,ylim,difft = self._draw_channel_gap2dt(auto=True)
			lenx = len(xlim); leny = len(ylim);

			if self.note_in['fcore'].get() == 1:
				ind = np.argmin(abs(xlim-1.))
				diff = difft[ind,:-2]
				ind = np.argmin(diff)
				self.note_in['EM'].set('%2.2f'%ylim[ind])
			else:
				ind = np.argmin(abs(ylim-1.))
				diff = difft[:-2,ind]
				ind = np.argmin(diff)
				self.note_in['CM'].set('%2.2f'%xlim[ind])
		else:
			ratios = np.linspace(0.5,2.0,41); vals = np.copy(ratios)
			if self.note_in['fcore'].get() == 1: self.note_ps['CORE_MIN'].set('0.9'); self.note_ps['CORE_MAX'].set('1.1'); self.note_ps['CORE_N'].set('3')
			else: self.note_ps['EDGE_MIN'].set('0.9'); self.note_ps['EDGE_MAX'].set('1.1'); self.note_ps['EDGE_N'].set('3')	

			for i in range(len(ratios)):
				self.note_in['CALCH'].set('%2.2f'%ratios[i])
				xlim,ylim,difft = self._draw_channel_gap2dt(auto=True)
				lenx = len(xlim); leny = len(ylim);

				if self.note_in['fcore'].get() == 1:
					ind = np.argmin(abs(xlim-1.))
					diff = difft[ind,:-2]
					ind = np.argmin(diff)
					vals[i] = ylim[ind]
				else:
					ind = np.argmin(abs(ylim-1.))
					diff = difft[:-2,ind]
					ind = np.argmin(diff)
					vals[i] = xlim[ind]

			ind = np.argmin(abs(vals-1.))
			self.note_in['CALCH'].set('%2.2f'%ratios[ind])

		self.note_ps['CORE_MIN'].set(cmint)
		self.note_ps['CORE_MAX'].set(cmaxt)
		self.note_ps['EDGE_MIN'].set(emint)
		self.note_ps['EDGE_MAX'].set(emaxt)
		self.note_ps['CORE_N'].set(cn)
		self.note_ps['EDGE_N'].set(en)
		if not self.nogui: self._draw_plots(True)		

		return

	def _declare_variables(self):

		self.figure    		   = dict()
		self.figure['name']    = dict()
		self.figure['size']    = dict()
		self.figure['canvas']  = dict()
		self.figure['widget']  = dict()
		self.figure['toolbar'] = dict()
		self.figure['legend']  = dict()
		self.figure['pegend']  = dict()
		self.figure['axes']    = dict()
		self.figure['gs']	   = dict()
		
		self.verplot           = []
		self.tciplot		   = []

		self.home       		= dict()
		self.note_in    		= dict()
		self.note_ps            = dict()

		self.profiles 			= dict()
		self.ts 				= dict()
		self.tci 			    = dict()
		self.int  		= dict()
		self.ref 		= dict()
		self.mds                = dict()
		return

	def _initialise_variables(self):

		self.ts_core_ex_preset      = []
		self.ts_edge_ex_preset      = []
		self.prev_page              = 0
		self.prev_plot_type         = -1
		self.prev_xmap              = 0
		self.prev_vline             = 0
		self.checked	            = tk.StringVar()
		self.checked.set('UNCHECK-A')

		self.tfig_list				= []
		self.pfig_list			    = []
		self.check_list				= []
		self.tci['nch'] = 5
		self.int['nch'] = 2

		self.profiles['shotn']		= 0
		self.profiles['didfit']     = False
		self.profiles['tmin']		= 0
		self.profiles['tmax']		= 0
		self.profiles['delt']		= 200
		self.profiles['times']		= np.array([])
		self.profiles['diag']       = dict()
		self.profiles['infile']     = dict()
		self.profiles['scale']      = dict()
		self.profiles['infile']['ts']     = dict()
		self.profiles['infile']['tciv']   = dict()
		self.profiles['infile']['tcis']   = dict()
		self.profiles['infile']['tcif']   = dict()
		self.profiles['infile']['gfile']  = dict()
		self.profiles['fit']			  = dict()
		self.profiles['fit']['nel']       = dict()
		self.profiles['fit']['ne']		  = dict()
		self.profiles['gflag'] 		      = 'EFIT01'

		self.ts['core'] = dict()
		self.ts['edge'] = dict()
		self.ts['core_err'] = dict()
		self.ts['edge_err'] = dict()		

		self.ts['core']['nch'] = 0
		self.ts['edge']['nch'] = 0

		for i in range(1,self.tci['nch']+1): 
			self.tci[i]       = dict()
			self.tci[i]['is'] = False

		for i in range(1,self.int['nch']+1): 
			self.int[i]       = dict()
			self.int[i]['is'] = False			
		self.ref['is'] = False

		self.profiles['diag']['ts'] = False
		self.profiles['diag']['ref']= False
		for i in range(1,6): self.profiles['diag']['tci%i'%i] = False
		for i in range(1,3): self.profiles['diag']['int%i'%i] = False

		for i in range(1,30):
			self.note_in['ts_core%i'%i] = tk.IntVar()
			self.note_in['ts_edge%i'%i] = tk.IntVar()
			self.note_ps['ts_core%i'%i] = tk.IntVar()
			self.note_ps['ts_edge%i'%i] = tk.IntVar()			
			self.note_in['tci%i'%i] = tk.IntVar()
			self.note_in['int%i'%i] = tk.IntVar()

		self.note_in['ref'] = tk.IntVar()

		self.note_in['WTS']    = tk.StringVar()
		self.note_in['WTCI01'] = tk.StringVar()
		self.note_in['WTCI02'] = tk.StringVar()
		self.note_in['WTCI03'] = tk.StringVar()
		self.note_in['WTCI04'] = tk.StringVar()
		self.note_in['WTCI05'] = tk.StringVar()
		self.note_in['WINT01'] = tk.StringVar()
		self.note_in['WINT02'] = tk.StringVar()
		self.note_in['WREFL']   = tk.StringVar()

		self.note_in['ATS']  = tk.StringVar()
		self.note_in['AINT'] = tk.StringVar()
		self.note_in['ATCI'] = tk.StringVar()
		self.note_in['AREFL'] = tk.StringVar()

		self.note_in['CM'] = tk.StringVar()
		self.note_in['EM'] = tk.StringVar()
		self.note_in['CALCH']    = tk.StringVar()
		self.note_in['CALCHC']   = tk.StringVar()
		self.note_in['CALCHE']   = tk.StringVar()
		self.note_in['MINPH'] = tk.StringVar()
		self.note_in['WIDTH'] = tk.StringVar()
		self.note_in['docal'] = tk.IntVar()
		self.note_in['fcore'] = tk.IntVar()
		self.note_in['rhofit'] = tk.IntVar()
		self.note_in['hmode'] = tk.IntVar()
		self.note_in['gflag'] = tk.IntVar()

		self.note_ps['fig']  = tk.IntVar()
		self.note_ps['xmap'] = tk.IntVar()
		self.note_ps['xmin'] = tk.StringVar()
		self.note_ps['xmax'] = tk.StringVar()
		self.note_ps['ymin'] = tk.StringVar()
		self.note_ps['ymax'] = tk.StringVar()

		self.note_ps['core_ch'] = tk.StringVar()
		self.note_ps['edge_ch'] = tk.StringVar()	
		self.note_ps['COREM'] = tk.StringVar()
		self.note_ps['EDGEM'] = tk.StringVar()		

		self.note_ps['CORE_MIN'] = tk.StringVar()
		self.note_ps['CORE_MAX'] = tk.StringVar()		
		self.note_ps['CORE_N'] = tk.StringVar()		

		self.note_ps['EDGE_MIN'] = tk.StringVar()
		self.note_ps['EDGE_MAX'] = tk.StringVar()		
		self.note_ps['EDGE_N'] = tk.StringVar()		

		self.figure['size']['home'] = (11,7.4)
		self.figure['size']['plot'] = (3.8,0.9)

		self.profiles['fitopt']		        = dict()
		self.profiles['fitopt']['weight']       = dict()
		self.profiles['fitopt']['weight']['ts'] = 1.
		self.profiles['fitopt']['weight']['ref']= 1.
		self.profiles['fitopt']['tci_sig']      = np.zeros(10)
		self.profiles['fitopt']['int_sig']      = np.zeros(10)
		
		return

	def _shot_preset(self):

		self.note_in['WTS'].set('0.0'); self.note_in['WREFL'].set('0.0'); self.note_in['WINT01'].set('0.0'); self.note_in['WINT02'].set('0.0');
		self.note_in['WTCI01'].set('0.0'); self.note_in['WTCI02'].set('0.0'); self.note_in['WTCI03'].set('0.0'); self.note_in['WTCI04'].set('0.0'); self.note_in['WTCI05'].set('0.0'); 
		self.note_in['ATS'].set('100.'); self.note_in['AREFL'].set('50.'); self.note_in['ATCI'].set('100.'); self.note_in['AINT'].set('100.')
		self.note_ps['xmin'].set('0.0'); self.note_ps['xmax'].set('1.1'); self.note_ps['ymin'].set('0.0'); self.note_ps['ymax'].set('5.0');
		self.note_ps['COREM'].set('0.35'); self.note_ps['EDGEM'].set('1.0');
		self.note_ps['CORE_MIN'].set('0.5'); self.note_ps['CORE_MAX'].set('1.5'); self.note_ps['CORE_N'].set('50')
		self.note_ps['EDGE_MIN'].set('0.5'); self.note_ps['EDGE_MAX'].set('1.5'); self.note_ps['EDGE_N'].set('50')

		self.note_in['hmode'].set(1); self.note_in['CM'].set('1.0'); self.note_in['EM'].set('1.0'); self.note_in['CALCH'].set('1.0');		
		#Rules for different years, can determin whether diagnostics are valid?.. if condits do not meet minimum requirement for diag. -> exit
		self.profiles['shotn'] = self.shotn
		self.tcis_preset = np.zeros(5)
		self.ints_preset = np.zeros(2)
		self.note_in['fcore'].set(1)
#		self.note_in['gflag'].set(1)
		self.note_in['MINPH'].set('0.4')
		self.note_in['WIDTH'].set('0.05')
		self.note_ps['xmap'].set(1)

		self.ts_core_ex_preset      = []
		self.ts_edge_ex_preset      = [7,8,9,10,11,12,13,14,15,16,17]
		return

	def _load_fittool(self):
		self.fit = fittool.fit_tool(nolog=True)
		if self.fast_mode:
			self.fit.fit_trunc = 11
			self.fit.line_trunc = 8

	def __init__(self,shotn):
		self.brutal    = False
		self.read_mode = False
		self.fast_mode = False
		self.nogui     = False
		self.efitv     = 0
		self.shotn     = shotn
		return

if __name__ == "__main__":

	import dena as dcal
	sta = time.time()
	now = time.gmtime(sta)
	inlen       = len(sys.argv)-1
	brutal_mode = False
	fast_mode   = False
	read_mode   = False
	nogui       = False
	efitv       = -1
	infile      = ''

	if inlen > 0:
		try:
			shotn = int(sys.argv[1])
			fast_mode = True
		except: print('>>> Invalid Shot# -> [#shot] [b/B/r/R] [infile]'); exit()
	else:   
		try:
			shotn = input('>>> Shot number [#]  ')
			shotn = int(shotn)
		except: print('>>> Invalid Shot number [#]  ')

	try:
		if   sys.argv[2]== 'b': brutal_mode = True; nogui = True;  efitv = 0;
		elif sys.argv[2]== 'B': brutal_mode = True; nogui = True;  efitv = 1;
		elif sys.argv[2]== 'r': read_mode   = True; nogui = False; efitv = 1;
		elif sys.argv[2]== 'R': read_mode   = True; nogui = True;  efitv = 1;
	except: pass

	efit_list = get_efit_list2(shotn)	
	line      = '>>> EFIT type: 0[RT1], '
	for i in range(1,6):
		if efit_list['isefit'][i]: line += '%i[%s], '%(i,efit_list['name'][i])

	if not (brutal_mode or read_mode):
		try: efitv = int(input(line))
		except: print('>>> Invalid EFIT'); exit()

		if not (efitv>=0 and efitv<6): print('>>> Invalid EFIT'); exit()
		if (efitv>0 and not efit_list['isefit'][efitv]): print('>>> Invalid EFIT'); exit()
		print(' ----------------------------------------------------------------')
		print('||            KSTAR Density diagnostics tool Ver %s            ||'%version['dena'])
		print('||                Profile & Channeal comparison                 ||')
		print('||      Developed by Plasma control group (PU) & PLARE (SNU)    ||')
		print('||        Report: sk42@princeton.edu / leeys1996@snu.ac.kr      ||')
		print('-----------------------------------------------------------------')
	else:   
		if brutal_mode: print('>>> Brutal mode')
		if read_mode:   print('>>> Read mode')
		if nogui:       print('>>> No GUI mode')

	if (read_mode):
		if inlen>2: infile = sys.argv[3];
		else: infile = rdena_db_dir+'/dena_%05i.sav'%shotn;
		if not os.path.isfile(infile): print('>>> No input for read'); exit();

	dcal = dcal.kstar_density_tool(shotn)
	dcal.brutal    = brutal_mode
	dcal.read_mode = read_mode
	dcal.fast_mode = fast_mode
	dcal.efit_list = efit_list
	dcal.nogui     = nogui
	dcal.efitv     = efitv
	dcal.in_file   = infile

	dcal.root = tk.Tk()
	dcal.root.title('KSTAR density diagnostic tool #%i'%shotn)
	dcal._den_tool()
	fin = time.time()
	print('>>> Elapsed time: %f s'%round((fin-sta),1))
