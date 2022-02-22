#!/usr/local/miniconda3/bin/python3
import os,sys, copy, time, pickle
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilename,asksaveasfilename
import numpy as np
from shutil import move, copyfile, copytree, rmtree
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from get_efit import *
from gfit_ece_multi import *
from progress import *
import knots_tool3
import eqdsk
from gefit_tool import make_batch_script, submit_batch_script, make_kfile_str, get_kfile_data, get_midRZ
from exec_dirs import version,author,comment,mds_dir,mds_tci,pythonc_exec,mds_over,mds_ref,mds_ts,mds_ces,mds_lit,mse_corr,mse_dir,python2_exec
from exec_dirs import efit_address,efit_source_dir,shotk, nubeam_dir2, efit_dir, mse_good_ch, gfit_dir, gzip_dir, mds_da, dena_dir
from MDS import mds
from netCDF4 import Dataset

currdir = os.getcwd()
class fgefittool:

	def fgefit(self):
		self._make_directories()
		
		self.root.protocol('WM_DELETE_WINDOW')
		self._declare_variables()
		self._initialise_variables()
		self._load_run()
		self._update_efit_list()

		self._make_home_canvas()
		self._make_note_frame()
		for k in range(4): self.root.grid_columnconfigure(k,weight=1)

		self._make_input_frame()
		self._make_profile_frame()
		self._make_mse_frame()
		self._make_beam_frame()
		self._make_efit_frame()

		b1 = tk.Button(self.root, text="RESET", height = 1, width = 5,command=lambda: self._reset_plot())	
		b1.grid(row=0, column=4,padx=5,sticky='e')	

		b1 = tk.Button(self.root,  text="SAVE", width = 5,command=lambda: self._save_run())	
		b1.grid(row=0,column=3,sticky='e')
		
		if len(self.note_in['times'])> 0.:
			print('>>> Load Previous setting')
			self._load(woplot=True)		
			for time in self.note_in['times']:
				self.note_in['l2'].insert('end',time)	
			self._draw_shot_plot()
			self._sync_gui_opt()
						
		self.root.resizable(0,0)
		self.root.mainloop()		

		return	
	
	def _get_years(self, shotn):

		for yy in shotk.keys():
			if shotn in shotk[yy]['shot']: return yy
		return yy

	def _make_directories(self):
		for dirs in ['profiles','efits','equil','output','datasave']:
			try: os.mkdir(dirs.upper())
			except: pass
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

	def _make_plot_canvas(self):
		self.figure['name']['plot'] = plt.figure(0,figsize=self.figure['size']['plot'])
		ax1 = self.figure['name']['plot'].add_subplot(1,1,1)

		self.figure['canvas']['plot'] = FigureCanvasTkAgg(self.figure['name']['plot'],master=self.home['page2'])
		self.figure['widget']['plot'] = self.figure['canvas']['plot'].get_tk_widget()
		self.figure['widget']['plot'].grid(row=16,column=0,columnspan=9)
		[ax] = self.figure['name']['plot'].axes
		ax.set_facecolor('b')
		ax.axis('off')
		return

	def _leftclick(self,event):
		page = self.home['nb'].tk.call(self.home['nb']._w,'identify','tab',event.x,event.y)
		if not page == '':
			if not self.curr_page == page:
				self.curr_page = page				
				self._sync_gui_opt(page)
				self._update_efit_list()
				self._update_option_menus()
				if page == 1: 
					self._make_profile_plot()
					self._load_fit()
				if page == 0: self._draw_shot_plot()

				if page == 2: self._make_mse_plot()

				if page == 3: self._make_beam_plot()

				if page == 4: self._make_efit_plot()

		return

	def _make_note_frame(self):

		titles = ['INPUT','PROFILE','MSE','H&CD','KEFIT','POST']
		self.titles = titles
		self.home['nb'] = ttk.Notebook(self.root,width=400,height=740)
		self.home['nb'].bind('<Button-1>',self._leftclick)
		for i in range(1,len(titles)+1):
			self.home['page%i'%i] = ttk.Frame(self.home['nb'])
			self.home['nb'].add(self.home['page%i'%i],text=titles[i-1])
			if i > 1: self.home['nb'].tab(i-1,state='disabled')
			for k in range(9): self.home['page%i'%i].grid_columnconfigure(k,weight=1)

		self.home['nb'].grid(row=1,column=0,columnspan=4,sticky='n')

		return

	def _make_input_frame(self):

		self.l1 = tk.Label(self.home['page1'], text="================== SHOT INFO. ===================",justify='center')
		self.l1.grid(row=0,column=0,columnspan=9,pady=5)

		self.l1 = tk.Label(self.home['page1'], text='SHOT [#]',justify='center')
		self.l1.grid(row=1, column=0,pady=3)

		self.note_in['e1'] = tk.Entry(self.home['page1'],width=9,justify='center')
		self.note_in['e1'].insert(10,self.note_in['shot'])
		self.note_in['e1'].grid(row=1,column=1,columnspan=2)
		self.note_in['e1'].bind('<Return>', lambda x: self._load())

		b1 = tk.Button(self.home['page1'],  text="LOAD", width = 4,command=lambda: self._load())
		b1.grid(row=1, column=3,columnspan=2,sticky='w')

		self.l1 = tk.Label(self.home['page1'], text='TIME [ms]',justify='center')
		self.l1.grid(row=1, column=5,pady=3)

		self.note_in['e2'] = tk.Entry(self.home['page1'],width=9,justify='center')
		self.note_in['e2'].grid(row=1,column=6,columnspan=2)		
		self.note_in['e2'].bind('<Return>', lambda x: self._add())

		b1 = tk.Button(self.home['page1'],  text="ADD", width = 3,command=lambda: self._add())
		b1.grid(row=1, column=8,columnspan=1,sticky='w')

		self.l1 = tk.Label(self.home['page1'], text="== AVAILABLE ==",justify='center')
		self.l1.grid(row=2,column=0,columnspan=5,pady=5,padx=45,sticky='w')
		self.l1 = tk.Label(self.home['page1'], text="== SELECTED ==",justify='center')
		self.l1.grid(row=2,column=4,columnspan=5,pady=5,padx=45,sticky='e')

		self.note_in['frame1'] = tk.Frame(self.home['page1'])
		self.note_in['frame1'].grid(row=3,column=0,columnspan=5,padx=20,sticky='w')
		self.note_in['s1'] = tk.Scrollbar(self.note_in['frame1'])
		self.note_in['s1'].pack(side='right',fill='y')
		self.note_in['l1'] = tk.Listbox(self.note_in['frame1'],yscrollcommand = self.note_in['s1'].set,height=7)
		self.note_in['l1'].pack(side='left',fill='x')
		self.note_in['s1']["command"] = self.note_in['l1'].yview
		self.note_in['l1'].bind('<ButtonRelease-1>',lambda x: self._click_list1())
		self.note_in['l1'].bind('<Double-1>',lambda x: self._add(True))
	
		self.note_in['frame2'] = tk.Frame(self.home['page1'])
		self.note_in['frame2'].grid(row=3,column=4,columnspan=5,padx=20,sticky='e')
		self.note_in['s2'] = tk.Scrollbar(self.note_in['frame2'])
		self.note_in['s2'].pack(side='right',fill='y')
		self.note_in['l2'] = tk.Listbox(self.note_in['frame2'],yscrollcommand = self.note_in['s2'].set,height=7)
		self.note_in['l2'].pack(side='left',fill='x')
		self.note_in['s2']["command"] = self.note_in['l2'].yview
		self.note_in['l2'].bind('<ButtonRelease-3>',lambda x: self._delete())

		self.l1 = tk.Label(self.home['page1'], text="Range[ms,ms]",justify='center')
		self.l1.grid(row=4,column=0,columnspan=2,pady=5,padx=3,sticky='w')
		self.l1 = tk.Label(self.home['page1'], text="Tdel [ms]",justify='center')
		self.l1.grid(row=4,column=5,pady=5)		

		self.note_in['e3'] = tk.Entry(self.home['page1'],width=13,justify='center')
		self.note_in['e3'].insert(10,self.note_in['trange'])
		self.note_in['e3'].grid(row=4,column=1,columnspan=4,padx=5)
		self.note_in['e4'] = tk.Entry(self.home['page1'],width=9,justify='center')
		self.note_in['e4'].insert(10,self.note_in['delt'])
		self.note_in['e4'].grid(row=4,column=6,columnspan=2)	
		self.note_in['e4'].bind('<Return>', lambda x: self._add2())

		b1 = tk.Button(self.home['page1'],  text="ADD", width = 3,command=lambda: self._add2())
		b1.grid(row=4, column=8,columnspan=1,sticky='w')

		self.l1 = tk.Label(self.home['page1'], text="================== DIAGNO INFO. ===================",justify='center')
		self.l1.grid(row=5,column=0,columnspan=9,pady=5)

		flag = ['TS','CES','ECE','INT1','INT2','TCI1','TCI2','TCI3','TCI4','TCI5','MSE','REFL']
		for i in range(4):
			for j in range(3):
				count = 5 + 3*i +j
				self.note_in['e%i'%count] = tk.Entry(self.home['page1'],width=13,justify='center')
				self.note_in['e%i'%count].insert(10,flag[count-5])
				self.note_in['e%i'%count].grid(row=6+i,column=3*j,columnspan=3,padx=3)
				self.note_in['e%i'%count].configure(state='readonly',fg='white',readonlybackground='tomato')

		self.l1 = tk.Label(self.home['page1'], text="=================== RUN SET-UP ====================",justify='center')
		self.l1.grid(row=10,column=0,columnspan=9,pady=5)				


#		b1 = tk.Button(self.home['page1'],  text="SAVE", width = 10,command=lambda: self._save_run())
#		b1.grid(row=11, column=0,columnspan=5)				
		return

	def _add(self,click=False,extnum=None):

		if not click:
			try: 
				if extnum == None: line = '  %06i ms'%int(float(self.note_in['e2'].get()))
				else: line = '  %06i ms'%int(float(extnum))
				tlist = list(self.note_in['l1'].get(0,'end'))
				ind = tlist.index(line)
			except: return

		else: 
			selection = self.note_in['l1'].curselection()
			if len(selection) == 0: return
			line = self.note_in['l1'].get(selection[0])

		ttime = float(line.split()[0])

		if ((ttime > self.note_in['tmax']) or (ttime < self.note_in['tmin'])): 
			print('>>> AVAIL TIME[ms] -> min:%5i,max:%5i'%(self.note_in['tmin'],self.note_in['tmax']))
			return

		tlist = self.note_in['l2'].get(0,'end')
		try: ind = tlist.index(line)
		except: ind = -1
		if ind > -1: return
		else:
			tlist = list(tlist)
			tlist.append(line)
			tlisti = []
			for k in tlist: tlisti.append('%i'%int(float(k.split()[0])))
			tlisti = np.array(tlisti,dtype='int')
			tlist = []
			tlisti = np.sort(tlisti)
			self.note_in['l2'].delete(0,'end')
			for k in tlisti:
				self.note_in['l2'].insert('end','  %06i ms'%int(k))

			self.figure['name']['home'].canvas.draw_idle()
			ax1 = self.figure['name']['home'].axes[0]
			ax2 = self.figure['name']['home'].axes[1]
			ax3 = self.figure['name']['home'].axes[2]
			ax4 = self.figure['name']['home'].axes[3]
			if not self.ver_plot == 0:
				self.ver_plot.remove()
				self.ver_plot = 0			
			lines1 = ax1.axvline(x=int(ttime)/1.e3,color='magenta',linestyle='--',linewidth=1.0,zorder=5)
			lines2 = ax2.axvline(x=int(ttime)/1.e3,color='magenta',linestyle='--',linewidth=1.0,zorder=5)
			lines3 = ax3.axvline(x=int(ttime)/1.e3,color='magenta',linestyle='--',linewidth=1.0,zorder=5)
			lines4 = ax4.axvline(x=int(ttime)/1.e3,color='magenta',linestyle='--',linewidth=1.0,zorder=5)
			self.ver_plott.append(line)
			self.ver_plots[1].append(lines1)
			self.ver_plots[2].append(lines2)
			self.ver_plots[3].append(lines3)
			self.ver_plots[4].append(lines4)
		return

	def _add2(self):

		try:
			line = self.note_in['e3'].get().split(',')
			delt = float(self.note_in['e4'].get())
			tmin = float(line[0])
			tmax = float(line[1])
			tt = tmin
			self._add(click=False,extnum=tt)
			while True:
				tt = tt + delt
				if tt > tmax: break
				self._add(click=False,extnum=tt)
		except: return	
		return

	def _delete(self):

		selection = self.note_in['l2'].curselection()
		if len(selection) == 0: return
		line = self.note_in['l2'].get(selection[0])
		self.note_in['l2'].delete(selection[0])
		ind = self.ver_plott.index(line)
		self.figure['name']['home'].canvas.draw_idle()
		self.ver_plott.remove(line)
		for k in range(1,5): 
			self.ver_plots[k][ind].remove()
			self.ver_plots[k].remove(self.ver_plots[k][ind])

		return

	def _click_list1(self):

		selection = self.note_in['l1'].curselection()
		if len(selection) == 0: return
		line = self.note_in['l1'].get(selection[0])
		time = float(line.split()[0])/1.e3

		self.figure['name']['home'].canvas.draw_idle()
		ax1 = self.figure['name']['home'].axes[0]
		if not self.ver_plot ==0: self.ver_plot.remove()
		self.ver_plot = ax1.axvline(x=time,color='goldenrod',linestyle='--',linewidth=1.0,zorder=5)
		return

	def _load(self,woplot=False):
		try: shotn = int(float(self.note_in['e1'].get()))
		except: return
		self.note_in['shot'] = shotn
		file1 = 'DATASAVE/%06i_list.dat'%shotn
		if not os.path.isfile(file1):
			ok, year = get_year(shotn)
			if ok==1: out_mag,out_mse = get_efit_list(year,shotn)
			else: return
			f = open(file1,'w')
			for k in out_mag:
				f.write('%i\n'%k)
			f.close()
		else:
			out_mag = []
			f = open(file1,'r')
			while True:
				line = f.readline()
				if not line: break
				out_mag.append(int(float(line)))
			f.close()

		self.note_in['l1'].delete(0,'end')
		for k in out_mag: 
			line = '  %06i ms'%int(float(k))
			self.note_in['l1'].insert('end',line)

		#LOAD MDS
		file2 = 'DATASAVE/%06i_ip.dat'%shotn
		file3 = 'DATASAVE/%06i_hcd.dat'%shotn
		file4 = 'DATASAVE/%06i_equ.dat'%shotn
		file5 = 'DATASAVE/%06i_etc.dat'%shotn
		file6 = 'DATASAVE/%06i_mse.dat'%shotn
		if not os.path.isfile(file2):
			g = mds('kstar',shotn)
			eq= mds('efit01',shotn)
			[t,ip] = g.get('\\pcrc03'); ip = -ip;
			self._dump_pickle(file2,[t,ip])
			print('>>> Load H&CD data...')
			self.hcd = self._load_hcd(g,file3)
			print('>>> Load EQUIL data...')
			self.equ = self._load_equil(eq,file4)
			self.etc = self._load_etc(g,file5)
			self.ipraw = [t,ip]
			print('>>> Load MSE data...')
			self.mse_data = self._load_mse(g,file6)

		else:
			self.ipraw = self._get_pickle(file2)
			print('>>> Load H&CD data...')
			self.hcd = self._get_pickle(file3)
			print('>>> Load EQUIL data...')
			self.equ = self._get_pickle(file4)
			self.etc = self._get_pickle(file5)	
			print('>>> Load MSE data...')
			self.mse_data = self._get_pickle(file6)

		try: self.mse_data['nch'] = len(self.mse_data['ch'])
		except: pass

		for years in shotk.keys():
			if ((shotn-shotk[years]['shot'][0])*(shotn-shotk[years]['shot'][-1])) < 0:
				tyear = years
				break
		sum_mse = 0
		for k in range(self.mse_data['nch']):
			sum_mse = sum_mse + self.note_ms['msec%02i'%(k+1)].get()
		
		if sum_mse == 0:
			for k in range(self.mse_data['nch']):
				self.note_ms['msec%02i'%(k+1)].set(mse_good_ch[tyear][k])

		self._load_diag(shotn)
		self.years = int(self._get_years(shotn))

		ind1 = np.where(self.ipraw[0]>2.)
		ind2 = np.where(self.ipraw[1][ind1]<1.e5)
		self.etime = self.ipraw[0][ind1][ind2][0]
		if not woplot: self._draw_shot_plot()
		
		self.home['nb'].tab(1,state='normal')
		self.home['nb'].tab(5,state='normal')
		if (self.diagno['mse']['is'] and self.diagno['TE']['is'] and self.diagno['NE']['is'] and self.diagno['TI']['is']): 
			self.home['nb'].tab(2,state='normal')
			self.home['nb'].tab(3,state='normal')
			self.home['nb'].tab(4,state='normal')

		if (self.diagno['mse']['is'] and (self.years>2019) and self.diagno['TI']['is']):
			self.home['nb'].tab(2,state='normal')
			self.home['nb'].tab(3,state='normal')
			self.home['nb'].tab(4,state='normal')

		return

	def _load_hcd(self,g,filename):

		t1,nbp1 =  g.get('\\nb11_pnb'); T1,nbe1 =  g.get('\\nb11_vg1');
		t2,nbp2 =  g.get('\\nb12_pnb'); T2,nbe2 =  g.get('\\nb12_vg1');
		t3,nbp3 =  g.get('\\nb13_pnb'); T3,nbe3 =  g.get('\\nb13_vg1');
		t4,nbp4 =  g.get('\\nb2a_pb'); T4,nbe4 =  g.get('\\nb2a_vg1');
		t5,nbp5 =  g.get('\\nb2b_pb'); T5,nbe5 =  g.get('\\nb2b_vg1');
		t6,nbp6 =  g.get('\\nb2c_pb'); T6,nbe6 =  g.get('\\nb2c_vg1');
		t7,ech1 =  g.get('\\ec1_pwr');  
		t8,ech2 =  g.get('\\ec2_pwr');
		t9,ech3 =  g.get('\\ec3_pwr');
		t10,ech4 =  g.get('\\ec4_pwr');
		hcd = dict()
		hcd['1A'] = [t1,nbp1,T1,nbe1]
		hcd['1B'] = [t2,nbp2,T2,nbe2]
		hcd['1C'] = [t3,nbp3,T3,nbe3]
		hcd['2A'] = [t4,nbp4,T4,nbe4]
		hcd['2B'] = [t5,nbp5,T5,nbe5]
		hcd['2C'] = [t6,nbp6,T6,nbe6]
		hcd['EC1'] = [t7,ech1]
		hcd['EC2'] = [t8,ech2]
		hcd['EC3'] = [t9,ech3]
		hcd['EC4'] = [t10,ech4]
		self._dump_pickle(filename,hcd)

		return hcd

	def _load_equil(self,eq,filename):
		equ = dict()
		equ['bp']   = eq.get('\\betap')
		equ['bn']   = eq.get('\\betan')
		equ['wmhd'] = eq.get('\\wmhd')
		equ['li']   = eq.get('\\li')
		equ['a0']   = eq.get('\\aminor')
		equ['r0']   = eq.get('\\r0')
		equ['f0']   = eq.get('\\fpol')
		equ['kappa']= eq.get('\\kappa')
		self._dump_pickle(filename,equ)

		return equ

	def _load_diag(self,shotn):

		currdir = os.getcwd()
		dir0   = currdir+'/PROFILES/MDS/'
		dirs   = currdir+'/PROFILES/MDS/DATASAVE'
		datdir = currdir+'/PROFILES/MDS/DATASAVE/%i'%shotn
		try: os.mkdir(dir0)
		except: pass		
		try: os.mkdir(dirs)
		except: pass
		try: os.mkdir(datdir)
		except: pass	
		os.chdir(dir0)

		#Load TS
		self.ts_ncore = 0
		if (not os.path.isfile(datdir+'/TE_size')):
			try: os.system('%s %i 0 0 0 1'%(mds_dir,shotn))
			except: pass
		if (os.path.isfile(datdir+'/TE_size')):
			f = open(datdir+'/TE_size','r')
			line = f.readline(); f.close()
			self.ts_ncore = int(line.split()[-1])
			print('>>> TS core CH#%i'%(self.ts_ncore))
		if (not os.path.isfile(datdir+'/NE_size')):
			try: os.system('%s %i 0 0 0 1'%(mds_dir,shotn))
			except: pass
		#Load CES
		if not os.path.isfile(datdir+'/TI_size'): os.system('%s %i 0 0 0 1'%(mds_dir,shotn))
		if not os.path.isfile(datdir+'/VT_size'): os.system('%s %i 0 0 0 1'%(mds_dir,shotn))
		#Load TCI
		if not os.path.isfile(datdir+'/TCI_size'): os.system('%s %i 0 0 0 1'%(mds_dir,shotn))
		#Load ECE
		if not os.path.isfile(datdir+'/ECE_size'): gefit_ece(shotn)

		self.diagno = self._load_diag_post(datdir,shotn)
		os.chdir(currdir)
		return

	def _load_diag_post(self,dirs,shotn):

		diagno = dict()
		for flag in ['TE','NE','TI','VT','ECE']:
			print('>>> Load %s data...'%flag)
			diagno[flag] = dict()
			file1 = dirs+'/%s_size'%flag
			if os.path.isfile(file1):
				f = open(file1,'r')
				line = f.readline().split()
				datn = int(float(line[-1]))
				if len(line)>1: datn2= int(float(line[-2]))
				f.close()
				if datn == 0: diagno[flag]['is'] = False
				else: 
					diagno[flag]['is'] = True
					if (flag=='TE' or flag=='NE'): 
						diagno[flag]['nch'] = datn2
					else:   diagno[flag]['nch'] = datn
			else: diagno[flag]['is'] = False
			if diagno[flag]['is']:
				diagno[flag]['dat'] = dict()
				diagno[flag]['R'] = np.zeros(diagno[flag]['nch'])
				for i in range(1,diagno[flag]['nch']+1):
					filen = dirs+'/%s%i.npz'%(flag.upper(),i)
					try:os.system('gzip -d %s.gz'%filen)
					except: pass
					npzfile = np.load(filen,mmap_mode='r')
					time = npzfile['arr_0']
					if not flag == 'ECE':
						self.note_in['tmin'] = max(self.note_in['tmin'],1000.*min(time))
						self.note_in['tmax'] = min(self.note_in['tmax'],1000.*max(time))
					if flag == 'ECE': 
						diagno[flag]['R'][i-1] = npzfile['arr_1'][i-1]
						diagno[flag]['dat'][i] = copy.deepcopy([time,npzfile['arr_2']])
					else:
						diagno[flag]['R'][i-1] = npzfile['arr_3']
						diagno[flag]['dat'][i] = copy.deepcopy([time,npzfile['arr_1'],npzfile['arr_2']])			
				
					os.system('gzip %s'%filen)

		file1 = dirs+'/TCI_size'	
		tci = ['int1','int2','tci1','tci2','tci3','tci4','tci5']
		print('>>> Load INTER data...')
		f = open(file1,'r')
		for k in range(7):
			diagno[tci[k]] = dict()
			line = f.readline().split()
			datn = int(float(line[0]))
			if datn > 0: diagno[tci[k]]['is'] = True
			else: diagno[tci[k]]['is'] = False

			if diagno[tci[k]]['is']:
				filen = dirs+'/'+'ne_%s.npz'%tci[k]
				comm='gzip -d '+filen+'.gz'
				os.system(comm)
				npzfile  = np.load(filen, mmap_mode='r')
				time=npzfile['arr_0']
				ne=npzfile['arr_1']
				comm='gzip '+filen
				os.system(comm)

				diagno[tci[k]]['dat'] = copy.deepcopy([time,ne])
		f.close()
		diagno['mse'] = dict(); diagno['mse']['is'] = False
		if len(self.mse_data['ch']) > 0: diagno['mse']['is'] = True
		diagno['refl'] = dict(); diagno['refl']['is'] = False

		flag1 = ['TE','TI','ECE','int1','int2','tci1','tci2','tci3','tci4','tci5','mse','refl']
		for i in range(4):
			for j in range(3):
				count = 5 + 3*i +j
				if diagno[flag1[count-5]]['is']:
					self.note_in['e%i'%count].configure(state='readonly',fg='white',readonlybackground='steelblue')		
				else:
					self.note_in['e%i'%count].configure(state='readonly',fg='white',readonlybackground='tomato')		
		return (diagno)

	def _load_etc(self,g,filename):

		etc = dict()
		etc['wdia'] = g.get('\\wtot_dlm03')
		[t,da] = g.get('\\pol_ha02')

		delt = t[101]-t[100];
		dmul = max(int(0.0001/delt),1)
		t = t[0::dmul]
		da=da[0::dmul]
		etc['dal'] = [t,da]
		self._dump_pickle(filename,etc)
		return etc

	def _load_mse(self,g,filename):

		mse_dat = dict()
		mse_dat['is'] = False
		R = g.get('\\RRRGAM')
		mse_dat['ch'] = R[0]
		if len(mse_dat['ch']) == 0: 
			self._dump_pickle(filename,mse_dat)
			return
		mse_dat['is'] = True
		mse_dat['RRRGAM'] = R[1]
		Z = g.get('\\ZZZGAM')
		mse_dat['ZZZGAM'] = Z[1]
		for k in range(6):
			A = g.get('\\AA%iGAM'%(k+1))
			mse_dat['AA%iGAM'%(k+1)] = A[1]

		T = g.get('\\TGAMMA01');
		mse_dat['time'] = T[0]
		mse_dat['TGAMMA01'] = T[1]
		T = g.get('\\SGAMMA01');
		mse_dat['SGAMMA01'] = T[1]

		for k in range(2,len(mse_dat['ch'])+1):
			T = g.get('\\TGAMMA%02i'%k);
			mse_dat['TGAMMA%02i'%k] = T[1]
			T = g.get('\\SGAMMA%02i'%k);
			mse_dat['SGAMMA%02i'%k] = T[1]

		self._dump_pickle(filename,mse_dat)

		return (mse_dat)

	def _dump_pickle(self,filename,dat):

		f = open(filename,'wb')
		pickle.dump(dat,f)
		f.close()
		return

	def _get_pickle(self,filename):
		f = open(filename,'rb')
		out = pickle.load(f)
		f.close()
		return out

	def _cut_init_end_index(self,arr,minv,maxv):

		arr2 = np.copy(arr)
		arr2[np.isnan(arr2)] = -1
		ind1 = np.where(arr2>=minv)
		ind2 = np.where(arr2[ind1]<=maxv)
		return ind1,ind2

	def _draw_shot_plot(self):

		self.figure['name']['home'].canvas.draw_idle()
		self.figure['name']['home'].clf()
		ax1 = self.figure['name']['home'].add_subplot(2,2,1)
		ax2 = self.figure['name']['home'].add_subplot(2,2,2,sharex=ax1)
		ax3 = self.figure['name']['home'].add_subplot(2,2,3,sharex=ax1)
		ax4 = self.figure['name']['home'].add_subplot(2,2,4,sharex=ax1)
		self.figure['name']['homecursor'] = MultiCursor(self.figure['name']['home'].canvas,[ax1,ax2,ax3,ax4],horizOn=True,color='r',lw=1)
		line, = ax1.plot(self.ipraw[0],abs(self.ipraw[1])/1.e6)
		slegend1 = [line]; plegend1=['Ip [MA]']
		ax1.set_xlim(-0.5,self.etime+1.)
		ax2.set_xlim(-0.5,self.etime+1.)
		ax1.set_xlabel('time [s]')
		ax2.set_xlabel('time [s]')
		ax3.set_xlabel('time [s]')
		ax4.set_xlabel('time [s]')
		ax1.set_ylabel('[a.u]')
		ax2.set_ylabel('[a.u]')
		ax3.set_ylabel('[a.u]')
		ax4.set_ylabel('[$10^{19}m^{-3}$]')

		for flag in ['1A','1B','1C','2A','2B','2C']:
			if (len(self.hcd[flag][0]))>0:
				line, = ax1.plot(self.hcd[flag][0],self.hcd[flag][1]*0.1)
				slegend1.append(line); plegend1.append('NB%s [10MW]'%flag)

		self.ver_plot = 0
		tlist = list(self.note_in['l2'].get(0,'end'))				

		self.ver_plott = []; 
		for k in range(1,5): self.ver_plots[k] = [];
		for k in tlist:
			line = k.split()[0]
			line1 = ax1.axvline(x=int(line)/1.e3,color='magenta',linestyle='--',linewidth=1.0,zorder=5)		
			line2 = ax2.axvline(x=int(line)/1.e3,color='magenta',linestyle='--',linewidth=1.0,zorder=5)		
			line3 = ax3.axvline(x=int(line)/1.e3,color='magenta',linestyle='--',linewidth=1.0,zorder=5)		
			line4 = ax4.axvline(x=int(line)/1.e3,color='magenta',linestyle='--',linewidth=1.0,zorder=5)		
			self.ver_plott.append(k)
			self.ver_plots[1].append(line1)
			self.ver_plots[2].append(line2)
			self.ver_plots[3].append(line3)
			self.ver_plots[4].append(line4)

		for flag in ['EC1','EC2','EC3','EC4']:
			if flag in self.hcd.keys():
				if (len(self.hcd[flag][0]))>0:
					line, = ax1.plot(self.hcd[flag][0],self.hcd[flag][1]/10000)	
					slegend1.append(line); plegend1.append('%s [10MW]'%flag)
		ax1.legend(slegend1,plegend1,loc='upper right')	

		slegend2 = []; plegend2 = ['$\\beta_p$','$\\beta_N$[%]','2$l_i$','$W_{MHD}$[200kJ]','$D_{\\alpha}$'];
		for flag in ['bp','bn']:
			line, = ax2.plot(self.equ[flag][0]/1.e3,self.equ[flag][1])
			slegend2.append(line); 
		line, = ax2.plot(self.equ['li'][0]/1.e3,2*self.equ['li'][1])
		slegend2.append(line)
		line, = ax2.plot(self.equ['wmhd'][0]/1.e3,self.equ['wmhd'][1]/2.e5)
		slegend2.append(line)
		len1 = len(self.etc['dal'][1])
		damax = max(self.etc['dal'][1][int(len1/4):int(len1/2)])
		damin = min(self.etc['dal'][1][0:int(len1/2)])
		if damin < 0.: 
			self.etc['dal'][1] = self.etc['dal'][1] - 1.2*damin
			damax = damax - 1.2*damin

		ind1,ind2 = self._cut_init_end_index(self.etc['dal'][1],-damax,2.*damax)

		line, = ax2.plot(self.etc['dal'][0][ind1][ind2],self.etc['dal'][1][ind1][ind2]/damax*0.6)
		slegend2.append(line)
		ax2.legend(slegend2,plegend2,loc='upper right')

		plegend3 = []; slegend3 = []; plegend4 = []; slegend4 =[];
		scale = [1.e3,  1.e19, 1.e3,  1.e2, 1.e3]
		maxv  = [8.e3,  1.e20, 1.5e4, 6.e2, 8.e3] 
		flags = ['TE','NE','TI','VT','ECE']
		titles = ['TS-TE[keV]','TS-NE[19/m3]','CES-TI[keV]','CES-VT[100km/s]','ECE[keV]']
		for k in range(5):
			flag = flags[k]
			if self.diagno[flag]['is']:
				ind1,ind2 = self._cut_init_end_index(self.diagno[flag]['dat'][5][1],0.,maxv[k])
				line, = ax3.plot(self.diagno[flag]['dat'][5][0][ind1][ind2],self.diagno[flag]['dat'][5][1][ind1][ind2]/scale[k])
				slegend3.append(line)
				plegend3.append('CH5.%s'%titles[k])

		for flag in ['int1','int2','tci1','tci2','tci3','tci4','tci5']:
			if self.diagno[flag]['is']:
				
				ind1,ind2 = self._cut_init_end_index(self.diagno[flag]['dat'][1],0.,7.0) 
				if flag == 'int1': line, = ax4.plot(self.diagno[flag]['dat'][0][ind1][ind2],self.diagno[flag]['dat'][1][ind1][ind2]/1.9)
				elif flag == 'int2': line, = ax4.plot(self.diagno[flag]['dat'][0][ind1][ind2],self.diagno[flag]['dat'][1][ind1][ind2]/2.75)
				else: line, = ax4.plot(self.diagno[flag]['dat'][0][ind1][ind2],self.diagno[flag]['dat'][1][ind1][ind2]/1.e0)
				slegend4.append(line)
				plegend4.append(flag.upper())

		ax3.legend(slegend3,plegend3,loc='upper right')
		ax4.legend(slegend4,plegend4,loc='upper right')

		ax1.set_ylim(0.)
		ax2.set_ylim(0.)
		ax3.set_ylim(0.)
		ax4.set_ylim(0.)

		self.figure['name']['home'].tight_layout()

		return

	def _make_profile_frame(self):

		self.l1 = tk.Label(self.home['page2'], text="================= TIME SLICE INFO. ==================",justify='center')
		self.l1.grid(row=0,column=0,columnspan=9,pady=5)

		self.l1 = tk.Label(self.home['page2'], text="== RAW-DATA ==",justify='center')
		self.l1.grid(row=2,column=0,columnspan=5,pady=5,padx=45,sticky='w')
		self.l1 = tk.Label(self.home['page2'], text="== POST-FIT ==",justify='center')
		self.l1.grid(row=2,column=4,columnspan=5,pady=5,padx=50,sticky='e')

		self.note_pr['frame1'] = tk.Frame(self.home['page2'])
		self.note_pr['frame1'].grid(row=3,column=0,columnspan=5,padx=20,sticky='w')
		self.note_pr['s1'] = tk.Scrollbar(self.note_pr['frame1'])
		self.note_pr['s1'].pack(side='right',fill='y')
		self.note_pr['l1'] = tk.Listbox(self.note_pr['frame1'],yscrollcommand = self.note_pr['s1'].set,height=7)
		self.note_pr['l1'].pack(side='left',fill='x')
		self.note_pr['s1']["command"] = self.note_pr['l1'].yview
		self.note_pr['l1'].bind('<Double-1>',lambda x: self._raw_draw())
	
		self.note_pr['frame2'] = tk.Frame(self.home['page2'])
		self.note_pr['frame2'].grid(row=3,column=4,columnspan=5,padx=20,sticky='e')
		self.note_pr['s2'] = tk.Scrollbar(self.note_pr['frame2'])
		self.note_pr['s2'].pack(side='right',fill='y')
		self.note_pr['l2'] = tk.Listbox(self.note_pr['frame2'],yscrollcommand = self.note_pr['s2'].set,height=7)
		self.note_pr['l2'].pack(side='left',fill='x')
		self.note_pr['s2']["command"] = self.note_pr['l2'].yview
		self.note_pr['l2'].bind('<Double-1>',lambda x: self._fit_draw())


		flag1 = ['TS-TE','TS-NE','CES-TI','CES-VT','ECE']
		flag3 = ['TS-TE','TSE-TE','TS-NE','TSE-NE','CES-TI','CES-VT','ECE']
		flag2 = ['TE','NE','TI','VT']
		self.note_pr['m1'] = tk.OptionMenu(self.home['page2'],self.note_pr['rawtype'],*flag1)
		self.note_pr['m1'].grid(row=4,column=0,columnspan=3,pady=5)
		self.note_pr['m1'].config(width=6)
		self.note_pr['m2'] = tk.OptionMenu(self.home['page2'],self.note_pr['fittype'],*flag2)
		self.note_pr['m2'].grid(row=4,column=5,columnspan=2,pady=5)
		self.note_pr['m2'].config(width=4)

		b1 = tk.Button(self.home['page2'],  text="CLEAR", width = 4,command=lambda: self._clf_plot(1))
		b1.grid(row=4, column=3,columnspan=2,sticky='w')
		b1 = tk.Button(self.home['page2'],  text="CLEAR", width = 4,command=lambda: self._clf_plot(2))
		b1.grid(row=4, column=7,columnspan=2,sticky='w')				
		b1 = tk.Button(self.home['page2'],  text="UPDATE", width = 4,command=lambda: self._update_plot(1))
		b1.grid(row=6, column=3,columnspan=2,sticky='w')
		b1 = tk.Button(self.home['page2'],  text="UPDATE", width = 4,command=lambda: self._update_plot(2))
		b1.grid(row=6, column=7,columnspan=2,sticky='w')

		self.l1 = tk.Label(self.home['page2'], text="YMAX",justify='center')
		self.l1.grid(row=5,column=0,pady=3)
		self.l1 = tk.Label(self.home['page2'], text='YMIN',justify='center')
		self.l1.grid(row=6,column=0)
		self.l1 = tk.Label(self.home['page2'], text="YMAX",justify='center')
		self.l1.grid(row=5,column=5)
		self.l1 = tk.Label(self.home['page2'], text="YMIN",justify='center')
		self.l1.grid(row=6,column=5)

		for k in range(4):
			self.note_pr['e%i'%(k+1)] = tk.Entry(self.home['page2'],width=6,justify='center')
		self.note_pr['e1'].insert(10,'4.0')
		self.note_pr['e2'].insert(10,'0.0')
		self.note_pr['e3'].insert(10,'4.0')
		self.note_pr['e4'].insert(10,'0.0')
		self.note_pr['e1'].grid(row=5,column=1,columnspan=2)
		self.note_pr['e2'].grid(row=6,column=1,columnspan=2)
		self.note_pr['e3'].grid(row=5,column=6,columnspan=1)
		self.note_pr['e4'].grid(row=6,column=6,columnspan=1)

		self.l1 = tk.Label(self.home['page2'], text="================ PROF. FITTING CONST. =================",justify='center')
		self.l1.grid(row=7,column=0,columnspan=9,pady=5)		

		title = ['TST','TSET','TSN','TSEN','CEST','CESV','ECE']
		for k in range(7):
			flag = flag3[k]
			self.l1 = tk.Label(self.home['page2'], text='Ch.Del [%s]'%title[k], justify='center')
			self.l1.grid(row=8+k,column=0,columnspan=2,pady=3,padx=3,sticky='w')
			self.note_pr['e%i'%(k+5)] = tk.Entry(self.home['page2'],width=26,justify='center')
			self.note_pr['e%i'%(k+5)].insert(10,self.note_pr['exclude'][flag])
			self.note_pr['e%i'%(k+5)].grid(row=8+k,column=1,columnspan=5,sticky='e',padx=2)
		title = ['TSEP[keV]','TWIDMIN','ZEFF','TS-TCM','TS-TEM','TS-NCM','TS-NEM']
		flag0 = ['TSEP','TWIDTH','ZEFF','TS_TCM','TS_TEM','TS_NCM','TS_NEM']
		for k in range(7):
			self.l1 = tk.Label(self.home['page2'], text=title[k], justify='center')
			self.l1.grid(row=8+k,column=6)
			self.note_pr['e%i'%(k+12)] = tk.Entry(self.home['page2'],width=6,justify='center')
			self.note_pr['e%i'%(k+12)].insert(10,self.note_pr[flag0[k]])
			self.note_pr['e%i'%(k+12)].grid(row=8+k,column=7,columnspan=2)

		self.l1 = tk.Label(self.home['page2'], text="================ DENSITY PROF. SCALE =================",justify='center')
		self.l1.grid(row=15,column=0,columnspan=9,pady=5)					

		title = ['CMAX','CMIN','EMAX','EMIN']
		flag0 = ['NCMAX','NCMIN','NEMAX','NEMIN']
		for k in range(4):
			self.l1 = tk.Label(self.home['page2'], text=title[k], justify='center')
			self.l1.grid(row=16,column=2*k)
			self.note_pr['e%i'%(k+19)] = tk.Entry(self.home['page2'],width=5,justify='center')
			self.note_pr['e%i'%(k+19)].insert(10,self.note_pr[flag0[k]])
			self.note_pr['e%i'%(k+19)].grid(row=16,column=2*k+1)

		title = ['CN','EN','DACRIT','DUTY']
		flag0 = ['NCN','NEN','DACRIT','DUTY']
		for k in range(4):
			self.l1 = tk.Label(self.home['page2'], text=title[k], justify='center')
			self.l1.grid(row=17,column=2*k)
			self.note_pr['e%i'%(k+23)] = tk.Entry(self.home['page2'],width=5,justify='center')
			self.note_pr['e%i'%(k+23)].insert(10,self.note_pr[flag0[k]])
			self.note_pr['e%i'%(k+23)].grid(row=17,column=2*k+1)			

		self.l1 = tk.Label(self.home['page2'], text="================ PROF TIME AVG.[ms] =================",justify='center')
		self.l1.grid(row=18,column=0,columnspan=9,pady=5)	

		flag0 = ['TS','CES','ECE','TCI']
		for k in range(4):
			self.l1 = tk.Label(self.home['page2'], text=flag0[k], justify='center')
			self.l1.grid(row=19,column=2*k)
			self.note_pr['e%i'%(k+27)] = tk.Entry(self.home['page2'],width=5,justify='center')
			self.note_pr['e%i'%(k+27)].insert(10,self.note_pr['TAVG'][flag0[k]])
			self.note_pr['e%i'%(k+27)].grid(row=19,column=2*k+1)		

		self.l1 = tk.Label(self.home['page2'], text="================= PROF INTER WEIG. ==================",justify='center')
		self.l1.grid(row=20,column=0,columnspan=9,pady=5)					
		flag0 = ['INT1','INT2','TCI1','TCI2']
		for k in range(4):
			self.l1 = tk.Label(self.home['page2'], text=flag0[k], justify='center')
			self.l1.grid(row=21,column=2*k)
			self.note_pr['e%i'%(k+31)] = tk.Entry(self.home['page2'],width=5,justify='center')
			self.note_pr['e%i'%(k+31)].insert(10,self.note_pr[flag0[k]])
			self.note_pr['e%i'%(k+31)].grid(row=21,column=2*k+1)		

		flag0 = ['TCI3','TCI4','TCI5']
		for k in range(3):
			self.l1 = tk.Label(self.home['page2'], text=flag0[k], justify='center')
			self.l1.grid(row=22,column=2*k)
			self.note_pr['e%i'%(k+35)] = tk.Entry(self.home['page2'],width=5,justify='center')
			self.note_pr['e%i'%(k+35)].insert(10,self.note_pr[flag0[k]])
			self.note_pr['e%i'%(k+35)].grid(row=22,column=2*k+1)	

		self.l1 = tk.Label(self.home['page2'], text="================ PROF FITTING OPT. =================",justify='center')
		self.l1.grid(row=23,column=0,columnspan=9,pady=5)				

		self.l1 = tk.Label(self.home['page2'], text='SHIFT', justify='center')
		self.l1.grid(row=24,column=0)
		self.note_pr['c2'] = tk.Checkbutton(self.home['page2'],variable=self.note_pr['ASHIFT'])
		self.note_pr['c2'].grid(row=24,column=1)

		self.l1 = tk.Label(self.home['page2'], text='HMODE', justify='center')
		self.l1.grid(row=24,column=2)
		self.note_pr['c3'] = tk.Checkbutton(self.home['page2'],variable=self.note_pr['HMODE'])
		self.note_pr['c3'].grid(row=24,column=3)	

		self.l1 = tk.Label(self.home['page2'], text='NSCALE', justify='center')
		self.l1.grid(row=24,column=4)
		self.note_pr['c1'] = tk.Checkbutton(self.home['page2'],variable=self.note_pr['ASCALE'])
		self.note_pr['c1'].grid(row=24,column=5)	

		self.l1 = tk.Label(self.home['page2'], text='FORCE', justify='center')
		self.l1.grid(row=25,column=0)
		self.note_pr['c4'] = tk.Checkbutton(self.home['page2'],variable=self.note_pr['FORCE_FIT'])
		self.note_pr['c4'].grid(row=25,column=1)	

		self.l1 = tk.Label(self.home['page2'], text='1DSCALE', justify='center')
		self.l1.grid(row=24,column=6)
		self.note_pr['c5'] = tk.Checkbutton(self.home['page2'],variable=self.note_pr['ASCALE1D'])
		self.note_pr['c5'].grid(row=24,column=7)

		self.note_pr['m3'] = tk.OptionMenu(self.home['page2'],self.note_pr['runeq'],*self.note_nb['gflag'],command=lambda value: self._sync_gui_opt(1))
		self.note_pr['m3'].grid(row=26,column=0,columnspan=2,pady=5)		
		self.note_pr['m3'].config(width=5)

		b1 = tk.Button(self.home['page2'],  text="DENA", width = 7,command=lambda: self._run_dena())
		b1.grid(row=26, column=2,columnspan=2,pady=5)

		b1 = tk.Button(self.home['page2'],  text="ELM CHK", width = 10,command=lambda: self._get_elm_peak())
		b1.grid(row=26, column=4,columnspan=2,pady=5)
		b1 = tk.Button(self.home['page2'],  text="RUN FIT", width = 10,command=lambda: self._run_fit())
		b1.grid(row=26, column=6,columnspan=2,pady=5)
		return	

	def _make_profile_plot(self):

		self.figure['name']['home'].canvas.draw_idle()
		self.figure['name']['home'].clf()
		self.figure['name']['homecursor'].visible = False
	
		ax1 = self.figure['name']['home'].add_subplot(2,2,1)
		ax2 = self.figure['name']['home'].add_subplot(2,2,2)
		ax3 = self.figure['name']['home'].add_subplot(2,2,3)
		ax4 = self.figure['name']['home'].add_subplot(2,2,4)
		ax1.cla(); ax2.cla(); ax3.cla(); ax4.cla()
		ax1.set_xlabel('R[m]')
		ax2.set_xlabel('$\\psi_N$')
		ax3.set_xlabel('times [s]')
		ax4.set_xlabel('times [s]')
		ax1.set_xlim(1.75,2.3)
		ax2.set_xlim(0.,1.1)
		ax1.set_ylabel('RAW [a.u]')
		ax2.set_ylabel('FIT [a.u]')
		ax3.set_xlabel('times [s]')
		ax3.set_ylabel('W [kJ]')
		ax4.set_ylabel('[$10^{19}m^{-3}$]')

		for k in range(4): 
			self.figure['legend'][k] = []
			self.figure['pegend'][k] = []

		return

	def _reduce_index(self,instr,ind):

		arr = instr.split(',')
		if len(arr) == 0: return instr
		if len(arr[0].split()) == 0: return instr
		arr = np.array(arr,dtype='int');
		line = '%i'%(arr[0]-ind)
		for i in range(1,len(arr)):
			line = line + ',%i'%(arr[i]-ind)
		return line

	def _make_fitopt_input(self,filename):

		f = open(filename,'wb')
		force_opt = dict()
		force_opt['shot'] = float(self.note_in['e1'].get())
		tlist = list(self.note_pr['l1'].get(0,'end'))
		force_opt['times'] = ''
		count = 0
		for time in tlist:
			tt = float(time.split()[0])
			if count ==0: force_opt['times'] = force_opt['times'] + '%i'%int(tt)
			else: force_opt['times'] = force_opt['times'] + ',%i'%int(tt)
			count = count + 1

		force_opt['time']  = int(float(tlist[0].split()[0]))
		force_opt['tsdt']  = int(float(self.note_pr['e27'].get()))
		force_opt['cesdt'] = int(float(self.note_pr['e28'].get()))
		force_opt['ecedt'] = int(float(self.note_pr['e29'].get()))
		force_opt['tcidt'] = int(float(self.note_pr['e30'].get()))

		force_opt['ex_tst'] = self._reduce_index(self.note_pr['e5'].get(),1)
		force_opt['ex_tset']= self._reduce_index(self.note_pr['e6'].get(),1)
		force_opt['ex_tsn'] = self._reduce_index(self.note_pr['e7'].get(),1)
		force_opt['ex_tsen']= self._reduce_index(self.note_pr['e8'].get(),1)
		force_opt['ex_cest']= self._reduce_index(self.note_pr['e9'].get(),1)
		force_opt['ex_cesv']= self._reduce_index(self.note_pr['e10'].get(),1)
		force_opt['ex_ece'] = self._reduce_index(self.note_pr['e11'].get(),1)

		force_opt['tsep']     = float(self.note_pr['e12'].get())
		force_opt['twidmin']  = float(self.note_pr['e13'].get())
		force_opt['zeff']     = float(self.note_pr['e14'].get())
		force_opt['tstcm']    = float(self.note_pr['e15'].get())
		force_opt['tstem']    = float(self.note_pr['e16'].get())
		force_opt['tsncm']    = float(self.note_pr['e17'].get())
		force_opt['tsnem']    = float(self.note_pr['e18'].get())

		force_opt['maxc']    = float(self.note_pr['e19'].get())
		force_opt['minc']    = float(self.note_pr['e20'].get())
		force_opt['maxe']    = float(self.note_pr['e21'].get())
		force_opt['mine']    = float(self.note_pr['e22'].get())
		force_opt['cn']      = float(self.note_pr['e23'].get())
		force_opt['en']      = float(self.note_pr['e24'].get())
		force_opt['dacrit']  = float(self.note_pr['e25'].get())
		force_opt['duty']    = float(self.note_pr['e26'].get())

		if self.note_pr['ASCALE'].get() == 1: force_opt['ascale'] = True
		else: force_opt['ascale'] = False

		if self.note_pr['ASCALE1D'].get() == 1: force_opt['ascale1d'] = True		
		else: force_opt['ascale1d'] = False

		force_opt['ts']      = float(self.note_pr['e27'].get())
		force_opt['ces']     = float(self.note_pr['e28'].get())
		force_opt['ece']     = float(self.note_pr['e29'].get())
		force_opt['tci']     = float(self.note_pr['e30'].get())

		force_opt['int01']   = float(self.note_pr['e31'].get())
		force_opt['int02']   = float(self.note_pr['e32'].get())
		force_opt['tci01']   = float(self.note_pr['e33'].get())
		force_opt['tci02']   = float(self.note_pr['e34'].get())
		force_opt['tci03']   = float(self.note_pr['e35'].get())
		force_opt['tci04']   = float(self.note_pr['e36'].get())
		force_opt['tci05']   = float(self.note_pr['e37'].get())

		if self.note_pr['ASHIFT'].get() == 1: force_opt['ashift'] = True
		else: force_opt['ashift'] = False
		if self.note_pr['HMODE'].get() == 1: force_opt['ishmode'] = True
		else: force_opt['ishmode'] = False
		if self.note_pr['FORCE_FIT'].get() == 1: force_opt['forcefit'] = True
		else: force_opt['forcefit'] = False

		pickle.dump(force_opt,f)
		f.close()

		return

	def _raw_draw(self):

		rscale= [1.e3,1.e3,1.e3,1.e3,1.e0]
		scale = [1.e3,1.e19,1.e3,1.e2,1.e3]
		tsch  = self.ts_ncore
		tstcm = float(self.note_pr['e15'].get())
		tstem = float(self.note_pr['e16'].get())
		tsncm = float(self.note_pr['e17'].get())
		tsnem = float(self.note_pr['e18'].get())

		selection = self.note_pr['l1'].curselection()
		if len(selection) == 0: return
		line = self.note_pr['l1'].get(selection[0])
		time = float(line.split()[0])/1.e3

		self.figure['name']['home'].canvas.draw_idle()
		ax1 = self.figure['name']['home'].axes[0]

		flag1 = ['TS-TE','TS-NE','CES-TI','CES-VT','ECE']
		flag2 = ['TE','NE','TI','VT','ECE']
		indexc= [5,7,9,10,11]
		ind0  = flag1.index(self.note_pr['rawtype'].get())
		flag  = flag2[ind0]
		if not self.diagno[flag]['is']: return

		t = self.diagno[flag]['dat'][1][0]
		ind = np.argmin(abs(t-time))
		val = np.copy(self.diagno[flag]['R'])
		err = np.copy(self.diagno[flag]['R'])
		for k in range(self.diagno[flag]['nch']):
			val[k] = self.diagno[flag]['dat'][k+1][1][ind]
			if not flag == 'ECE': err[k] = self.diagno[flag]['dat'][k+1][2][ind]

		factors = np.ones(len(val))
		if (flag == 'TE'): 
			for k in range(len(val)): 
				if (k>=tsch): factors[k] = tstem	
				else: factors[k] = tstcm
		if (flag == 'NE'): 		
			for k in range(len(val)): 
				if (k>=tsch): factors[k] = tsnem
				else: factors[k] = tsncm

		if flag == 'ECE': line, = ax1.plot(self.diagno[flag]['R'],val*factors/scale[ind0],'-x')
		else:
			xx = self.diagno[flag]['R']/rscale[ind0]
			yy = val*factors/scale[ind0]
			ye = err*factors/scale[ind0]
			if (flag=='TE' and flag=='NE'):
				line = ax1.errorbar(xx[:self.ts_ncore],yy[:self.ts_ncore],yerr=ye[:self.ts_ncore])
				line = ax1.errorbar(xx[self.ts_ncore:],yy[self.ts_ncore:],yerr=ye[self.ts_ncore:])
			else:
				line = ax1.errorbar(xx,yy,ye)
		for k in range(self.diagno[flag]['nch']):
			if flag == 'ECE': 
				txt = '%i'%(self.diagno[flag]['nch']-k)
				ax1.text(self.diagno[flag]['R'][k]/rscale[ind0],val[k]*factors[k]/scale[ind0],txt)
			else:
				k2 = k+1; txt = '%i'%k2 
				if (flag=='TE' or flag=='NE'):
					if k>=self.ts_ncore: k2 = k-self.ts_ncore+1; txt='E%i'%k2
					else: txt='C%i'%k2
				ax1.text(self.diagno[flag]['R'][k]/rscale[ind0],val[k]*factors[k]/scale[ind0],txt)

		self.figure['legend'][0].append('%s - %i ms'%(flag1[ind0],t[ind]*1.e3))
		self.figure['pegend'][0].append(line)
	
		indexs = self.note_pr['e%i'%indexc[ind0]].get().split(',');
		indexe = []
		if ind0 < 2: indexe = self.note_pr['e%i'%(indexc[ind0]+1)].get().split(',');
		if len(indexs) > 0:
			if indexs[0] == '': indexs = []; 
		if len(indexe) > 0:
			if indexe[0] == '': indexe = [];
		if len(indexs) > 0:
			for i in indexs:
				if flag == 'ECE': 
					ii = self.diagno[flag]['nch']-int(float(i))-2
					ax1.scatter(self.diagno[flag]['R'][ii],val[ii]*factors[ii]/scale[ind0],s=10,c='gray')
				else:
					ii = int(float(i))-1;
					ax1.errorbar(self.diagno[flag]['R'][ii]/rscale[ind0],val[ii]*factors[ii]/scale[ind0],yerr=err[ii]*factors[ii]/scale[ind0],c='gray')
		if len(indexe) > 0:
			for i in indexe:
				ii = int(float(i)) + self.ts_ncore - 1
				ax1.errorbar(self.diagno[flag]['R'][ii]/rscale[ind0],val[ii]*factors[ii]/scale[ind0],yerr=err[ii]*factors[ii]/scale[ind0],c='gray')

		ymax2 = 1.2*max((val*factors)[0:15]/scale[ind0])
		ymax = float(self.note_pr['e1'].get())
		if ymax < ymax2:
			self.note_pr['e1'].delete(0,'end')
			self.note_pr['e1'].insert(10,'%4.2f'%ymax2)
			ax1.set_ylim(0,ymax2)					
	
		ax1.legend(self.figure['pegend'][0],self.figure['legend'][0])
		self._update_plot(1)
		self.figure['name']['home'].tight_layout()

		return

	def _fit_draw(self):

		selection = self.note_pr['l2'].curselection()
		if len(selection) == 0: return
		line = self.note_pr['l2'].get(selection[0])
		time = int(line.split()[0])

		data = self.mfit_opt[time]

		count = len(self.figure['pegend'][1])
		if count == 10: self._clf_plot(2)

		self.figure['name']['home'].canvas.draw_idle()
		ax1 = self.figure['name']['home'].axes[1]

		flag1 = ['TE','NE','TI','VT']
		flag2 = ['ts','ts','ces','ces']
		ind0 = flag1.index(self.note_pr['fittype'].get())
		flag = flag1[ind0].lower()
		flags= flag2[ind0]
		nts  = not self.diagno['TE']['is']
		if (nts and (flag1[ind0].lower() == 'te')): flags = 'ece'
		if not data[5]['didfit'][flag]: return

		count = len(self.figure['pegend'][1])

		ymax1 = max(data[ind0+1]['fit2'])
		ymax2 = max(data[ind0+1][flags]['raw'][0:15])
		ymax2 = ymax1 * 1.1

		line, = ax1.plot(data[6]['psin2'],data[ind0+1]['fit2'],'--',color='C%i'%count,zorder=2)
		self.figure['legend'][1].append('%s - %i ms %s'%(flag,time,self.note_pr['runeq'].get()))
		self.figure['pegend'][1].append(line)
		len1 = len(data[ind0+1][flags]['xxr'])
		rind  = np.array(np.linspace(0,len1-1,len1),dtype='int')
		rindo = np.delete(rind,data[5]['fit_ind'][flag][flags])
		rindi = np.copy(data[5]['fit_ind'][flag][flags])
		len3 = len(data[ind0+1][flags]['avg'])

		if (data[0]['avg'][flag][flags] and len3 > 0):
			len2 = int(len1 / len3)
			for i in range(1,len2): rindi = np.hstack([rindi,data[5]['fit_ind'][flag][flags]+len3*i])
			rindo = np.delete(rind,rindi)
		if not data[0]['file'][flag][flags]==None:
			line = ax1.errorbar(data[ind0+1][flags]['xxr'][rindi],data[ind0+1][flags]['raw'][rindi],yerr=data[ind0+1][flags]['raws'][rindi],fmt='x',c='C%i'%count,zorder=1)
			ax1.errorbar(data[ind0+1][flags]['xxr'][rindo],data[ind0+1][flags]['raw'][rindo],yerr=data[ind0+1][flags]['raws'][rindo],fmt='x',c='gray',zorder=1)

		if (flags=='ts' and not data[0]['file'][flag][flags]==None):
			flage = 'tse'; len1e = len(data[ind0+1][flage]['xxr'])
			rinde = np.array(np.linspace(0,len1e-1,len1e),dtype='int')
			rindeo= np.delete(rinde,data[5]['fit_ind'][flag][flage])
			rindei= np.copy(data[5]['fit_ind'][flag][flage])
			len3e = len(data[ind0+1][flage]['avg']) 
			if (data[0]['avg'][flag][flage] and len3e > 0):			
				len2e = int(len1e / len3e)
				for i in range(1,len2e): rindei = np.hstack([rindei,data[5]['fit_ind'][flag][flage]+len3e*i])
				rindeo = np.delete(rinde,rindei)
					
			ax1.errorbar(data[ind0+1][flage]['xxr'][rindei],data[ind0+1][flage]['raw'][rindei],yerr=data[ind0+1][flage]['raws'][rindei],fmt='x',c='C%i'%count,zorder=1)
			ax1.errorbar(data[ind0+1][flage]['xxr'][rindeo],data[ind0+1][flage]['raw'][rindeo],yerr=data[ind0+1][flage]['raws'][rindeo],fmt='x',c='gray',zorder=1)	

		if flag == 'te':
			if not data[0]['file'][flag]['ece']==None: 
				ymax1 = max(data[ind0+1]['ece']['raw'][0:15])
				len1 = len(data[ind0+1]['ece']['xxr'])
				rind = np.array(np.linspace(0,len1-1,len1),dtype='int')
				rindo = np.delete(rind,data[5]['fit_ind'][flag]['ece'])
				rindi = data[5]['fit_ind'][flag]['ece']
				line = ax1.scatter(data[ind0+1]['ece']['xxr'][rindi],data[ind0+1]['ece']['avg'][rindi],s=20,marker='+',c='C%i'%count,zorder=1)
				ax1.scatter(data[ind0+1]['ece']['xxr'][rindo],data[ind0+1]['ece']['avg'][rindo],s=20,marker='+',c='gray',zorder=1)

		ymax2= max(ymax1,ymax2)
		ymax = float(self.note_pr['e3'].get())
		if ymax < ymax2:
			self.note_pr['e3'].delete(0,'end')
			self.note_pr['e3'].insert(10,'%4.2f'%ymax2)
			ax1.set_ylim(0,ymax2)		

		ax1.legend(self.figure['pegend'][1],self.figure['legend'][1])
		self._update_plot(2)
		self.figure['name']['home'].tight_layout()	
		return

	def _clf_plot(self,type):

		self.figure['name']['home'].canvas.draw_idle()
		if type==1: ax = self.figure['name']['home'].axes[0]
		elif type==2: ax = self.figure['name']['home'].axes[1]
		elif type==3: ax = self.figure['name']['home'].axes[0]
		elif type>=4: ax = self.figure['name']['home'].axes[0]
		ax.cla()
		if type==1:ax.set_xlabel('R[m]')
		elif type==2: ax.set_ylabel('$\\psi_N$')
		elif type==3: ax.set_ylabel('R[m]')
		if type==1:ax.set_xlim(1.75,2.3)	
		elif type==2: ax.set_xlim(0,1.2)
		elif type==3: ax.set_xlim(1.7,2.4)
		if type==1: ax.set_ylabel('RAW [a.u]')	
		elif type==2: ax.set_ylabel('RAW [a.u]')	
		elif type==3: 
			ax.set_ylabel('T-GAMMA [a.u]')	
			ax.axhline(y=0,linestyle='--',c='gold')

		if type < 3:
			self.figure['legend'][type-1] = []
			self.figure['pegend'][type-1] = []
		elif type >= 3:
			self.figure['legend'][0] = []
			self.figure['pegend'][0] = []

		self._update_plot(type)
		return

	def _update_plot(self,type):

		self.figure['name']['home'].canvas.draw_idle()
		if type == 1: 
			ax = self.figure['name']['home'].axes[0]
			ymax = float(self.note_pr['e1'].get())
			ymin = float(self.note_pr['e2'].get())
			ax.set_xlim(1.75,2.3)
		elif type == 2: 
			ax = self.figure['name']['home'].axes[1]
			ymax = float(self.note_pr['e3'].get())
			ymin = float(self.note_pr['e4'].get())
			ax.set_xlim(0.,1.1)
		elif type == 3: 
			ax = self.figure['name']['home'].axes[0]
			ymax = float(self.note_ms['e2'].get())
			ymin = float(self.note_ms['e3'].get())
			ax.set_xlim(1.7,2.4)

		elif type == 4:
			ax = self.figure['name']['home'].axes[0]
			ymax = float(self.note_nb['e1'].get())
			ymin = float(self.note_nb['e2'].get())

		elif type == 5:
			ax = self.figure['name']['home'].axes[0]
			ymax = float(self.note_ef['e1'].get())
			ymin = float(self.note_ef['e2'].get())

		ax.set_ylim(ymin,ymax)
		self.figure['name']['home'].tight_layout()
		return

	def _reset_plot(self):

		if self.curr_page == 0: self._draw_shot_plot()
		if self.curr_page == 1:
			self.figure['name']['home'].canvas.draw_idle()
			ymax1 = float(self.note_pr['e1'].get())
			ymin1 = float(self.note_pr['e2'].get())			
			ymax2 = float(self.note_pr['e3'].get())
			ymin2 = float(self.note_pr['e4'].get())			
			self.figure['name']['home'].axes[0].set_ylim(ymin1,ymax1)
			self.figure['name']['home'].axes[0].set_xlim(1.75,2.3)
			self.figure['name']['home'].axes[1].set_ylim(ymin2,ymax2)
			self.figure['name']['home'].axes[1].set_xlim(0.,1.1)

		nlist = ['ms','nb','ef']

		if self.curr_page > 1: 
			self.figure['name']['home'].canvas.draw_idle()
			ax = self.figure['name']['home'].axes[0]
			ax.set_xlim(self.__dict__['note_%s'%nlist[self.curr_page-2]]['xmin'],self.__dict__['note_%s'%nlist[self.curr_page-2]]['xmax'])
			ax.set_ylim(self.__dict__['note_%s'%nlist[self.curr_page-2]]['ymin'],self.__dict__['note_%s'%nlist[self.curr_page-2]]['ymax'])

		return

	def _check_gfiles(self,eqver=''):

		dirs = 'EFITS/%s'%eqver
		shot = int(float(self.note_in['e1'].get()))
		tlist = list(self.note_pr['l1'].get(0,'end'))

		for time in tlist:
			tt = int(float(time.split()[0]))
			filen = dirs+'/g%06i.%06i'%(shot,tt)
			if not os.path.isfile(filen): 
				print('>>> No %s'%filen)
				return False

		return True

	def _run_dena(self):

		currdir = os.getcwd()
		dirs    = currdir +'/PROFILES'
		os.chdir(dirs)
		if os.path.isfile('run_state'): os.remove('run_state')
		os.system('%s %i &'%(dena_dir,self.note_in['shot']))
		count = 0
		os.chdir(currdir)
		while count <=6000:
			if os.path.isfile(dirs+'/run_state'): break
			time.sleep(2)
			count +=1
		return

	def _get_elm_peak(self):

		currdir = os.getcwd()
		dirs = currdir+'/PROFILES/MDS/DATASAVE/%i'%self.note_in['shot']

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

		tlist = list(self.note_pr['l1'].get(0,'end'))	
		tmin  = float(tlist[0].split()[0]); tmax = float(tlist[-1].split()[0])
		tmin  = tmin/1.e3 - 0.3
		tmax  = tmax/1.e3 + 0.3

		ind1 = np.where(self.da[0]>tmin);
		ind2 = np.where(self.da[0][ind1]<tmax);
		ind3 = int(round(0.003/(self.da[0][4]-self.da[0][3])))
		ind4 = int(round(0.00003/(self.da[0][4]-self.da[0][3])))

		tt = self.da[0][ind1][ind2]
		yy = self.da[1][ind1][ind2]
		lent = len(tt)

		mean_sig = np.mean(yy);
		if mean_sig>0: base = mean_sig*0.9
		else: base = -mean_sig*1.1

		dacrit = float(self.note_pr['e25'].get())
		ii=0
		while ii<lent:
			if (yy[ii]-base) > dacrit*base:
				self.peaks.append(tt[ii]*1.e3);
				ii += ind3
			else: ii += ind4
		
		dirs2  = dirs+'/elm_peaks.npz'
		np.savez(dirs+'/elm_peaks',self.peaks)

		tmin = tmin*1000; tmax = tmax*1000
		os.system('%s %s %s %i %i %s'%(pythonc_exec,mds_da,dirs+'/DA.npz',tmin,tmax,dirs2))

		return

	def _run_fit(self):

		currdir = os.getcwd()
		isgfile = self._check_gfiles(self.note_pr['runeq'].get())
		if not isgfile: 
			if self.note_pr['runeq'].get() == 'EFIT01': 
				print('>>> Get EFIT01...')
				self._get_gfile()

		try:rmtree('PROFILES/GFILES')
		except: pass
		try:os.mkdir('PROFILES/MFIT')
		except: pass
		copytree('EFITS/%s'%self.note_pr['runeq'].get(),'PROFILES/GFILES')

		os.chdir('PROFILES')
		self._make_fitopt_input('mult_opt.dat')
		if os.path.isdir('MPROFILES_%s'%self.note_pr['runeq'].get()):
			move('MPROFILES_%s'%self.note_pr['runeq'].get(),'MPROFILES')
		if os.path.isfile('result_multi.save_%s'%self.note_pr['runeq'].get()):
			move('result_multi.save_%s'%self.note_pr['runeq'].get(),'result_multi.save')
			with open('run_mode','w') as f: f.write('multi')

		try: os.system('%s -mult'%gfit_dir)
		except: pass
		os.chdir(currdir)
		file1 = 'PROFILES/result_multi.save'
		if not os.path.isfile(file1): return
		f = open(file1,'rb')
		self.mfit_opt = pickle.load(f)
		f.close()
		#self._make_profile_plot()
		self._save_run()
		try: rmtree('PROFILES/MPROFILES_%s'%self.note_pr['runeq'].get())
		except: pass
		move('PROFILES/MPROFILES','PROFILES/MPROFILES_%s'%self.note_pr['runeq'].get())
		move('PROFILES/result_multi.save','PROFILES/result_multi.save_%s'%self.note_pr['runeq'].get())
		self._sync_gui_opt(1)
		self._make_profile_plot()
		self._load_fit()
		return

	def _load_fit(self,woplot=False):

		file1 = 'PROFILES/result_multi.save_%s'%self.note_pr['runeq'].get()
		if not os.path.isfile(file1): return

		tlist = list(self.note_pr['l1'].get(0,'end'))
		tlist2 = copy.deepcopy(tlist)
		for k in range(len(tlist)):
			tt = float(tlist[k].split()[0])
			tlist2[k] = int(tt)

		self.note_pr['l2'].delete(0,'end')
		for time in self.mfit_opt['times']:
			try: 
				ind = tlist2.index(time)
				self.note_pr['l2'].insert('end',tlist[ind])
			except: pass

		self._draw_0d_profiles()

		return

	def _draw_0d_profiles(self):

		lscale  = [1.0,0.5,1.5,2.0,2.5,3.0,3.5]
		lscalet = ['x1.0','x0.5','x1.5','x2.0','x2.5','x3.0','x3.5']

		self.figure['name']['home'].canvas.draw_idle()
		if (len(self.figure['name']['home'].axes) <3): self._make_profile_plot()

		ax3 = self.figure['name']['home'].axes[2]
		ax4 = self.figure['name']['home'].axes[3]

		ax3.cla(); ax4.cla();

		wth = np.copy(self.mfit_opt['times'])
		wfa = np.copy(wth)
		tci = np.zeros((len(wth),7))
		for k in range(len(wth)):
			wth[k] = self.mfit_opt['post'][self.mfit_opt['times'][k]][3][1]
			wfa[k] = self.mfit_opt['post'][self.mfit_opt['times'][k]][3][0] - wth[k]
			for j in range(7):
				tci[k,j] = self.mfit_opt['post'][self.mfit_opt['times'][k]][0][j]

		plegend3 = []; slegend3 = [];
		line, = ax3.plot(self.equ['wmhd'][0]/1.e3,self.equ['wmhd'][1]/1.e3)
		plegend3.append(line); slegend3.append('$W_{MHD}$')
		line, = ax3.plot(self.mfit_opt['times']/1.e3,wth,'-x')
		plegend3.append(line); slegend3.append('$W_{TH}$')
		line, = ax3.plot(self.mfit_opt['times']/1.e3,wfa,'-x')
		plegend3.append(line); slegend3.append('$W_{GAP}$')		
		ax3.legend(plegend3,slegend3)
		ax3.set_xlabel('time [s]')
		ax3.set_ylabel('W [kJ]')
		ax3.set_ylim(0.)

		slegend4 = []; plegend4 = [];
		flag1 = ['int1','int2','tci1','tci2','tci3','tci4','tci5']
		flag2 = ['int01','int02','tci01','tci02','tci03','tci04','tci05']
		scale = [1./1.9,1./2.75,1.,1.,1.,1.,1.]
		lcount = 0
		for k in range(len(flag1)):
			flag = flag1[k]
			if (self.diagno[flag]['is'] and self.mfit_opt['post'][self.mfit_opt['times'][0]][5][k]>0.): 
				line, = ax4.plot(self.diagno[flag]['dat'][0],self.diagno[flag]['dat'][1]/1.e0*lscale[lcount] * scale[k],zorder=0)
				plegend4.append(line); slegend4.append('%s(%s)'%(flag.upper(),lscalet[lcount]))
				line, = ax4.plot(self.mfit_opt['times']/1.e3,tci[:,k]*lscale[lcount],'-x',zorder=1)
				plegend4.append(line); slegend4.append('%s-F'%flag.upper())
				lcount = lcount + 1

		ax4.legend(plegend4,slegend4)
		ax4.set_ylim(0.)	
		ax4.set_xlabel('time [s]')
		ax4.set_ylabel('[$10^{19}m^{-3}$]')

		return

	def _get_gfile(self,efit02=False):

		currdir = os.getcwd()
		dirs = 'EFITS/%s'%self.note_pr['runeq'].get()
		if not self.note_pr['runeq'].get() == 'EFIT01': return

		try: os.mkdir(dirs)
		except: pass
		os.chdir(dirs)
		shotn = int(float(self.note_in['e1'].get()))
		for years in shotk.keys():
			if ((shotn-shotk[years]['shot'][0])*(shotn-shotk[years]['shot'][-1])) < 0:
				tyear = years
				break
		efitdir = efit_source_dir + shotk[tyear]['mag'] + 'EXP%06i/'%shotn
		if efit02: efitdir = efit_source_dir + shotk[tyear]['mse'] + 'EXP%06i/'%shotn

		tlist = list(self.note_pr['l1'].get(0,'end'))
		if (len(tlist) ==0 ): 
			os.chdir(currdir)
			return
		fileline = ''
		for time in tlist:
			tt = int(float(time.split()[0]))
			gfilen = 'g%06i.%06i'%(shotn,tt)
			kfilen = 'k%06i.%06i'%(shotn,tt)
			if not os.path.isfile(gfilen): fileline = fileline + ' %sg%06i.%06i'%(efitdir,shotn,tt)
			if not os.path.isfile(kfilen): fileline = fileline + ' %sk%06i.%06i'%(efitdir,shotn,tt)
		
		if not fileline == '':
			comm = 'scp %s:"%s" .'%(efit_address,fileline)
			comm = 'cp %s .'%(fileline)
			os.system(comm)

		os.chdir(currdir)
		return

	def _make_mse_frame(self):

		self.l1 = tk.Label(self.home['page3'], text="================= TIME SLICE INFO. ==================",justify='center')
		self.l1.grid(row=0,column=0,columnspan=9,pady=5)

		self.l1 = tk.Label(self.home['page3'], text="== POST_FIT ==",justify='center')
		self.l1.grid(row=2,column=0,columnspan=5,pady=5,padx=45,sticky='w')

		self.note_ms['frame1'] = tk.Frame(self.home['page3'])
		self.note_ms['frame1'].grid(row=3,column=0,columnspan=5,padx=20,rowspan=5,sticky='w')
		self.note_ms['s1'] = tk.Scrollbar(self.note_ms['frame1'])
		self.note_ms['s1'].pack(side='right',fill='y')
		self.note_ms['l1'] = tk.Listbox(self.note_ms['frame1'],yscrollcommand = self.note_ms['s1'].set,height=7)
		self.note_ms['l1'].pack(side='left',fill='x')
		self.note_ms['s1']["command"] = self.note_ms['l1'].yview
		self.note_ms['l1'].bind('<Double-1>',lambda x: self._raw_mse_draw())

		self.l1 = tk.Label(self.home['page3'], text="TAVG [ms]",justify='center')
		self.l1.grid(row=5,column=5,pady=3,padx=5,sticky='w')
		self.l1 = tk.Label(self.home['page3'], text='YMAX',justify='center')
		self.l1.grid(row=6,column=5,padx=5,sticky='w')
		self.l1 = tk.Label(self.home['page3'], text="YMIN",justify='center')
		self.l1.grid(row=7,column=5,padx=5,sticky='w')

		for k in range(3):
			self.note_ms['e%i'%(k+1)] = tk.Entry(self.home['page3'],width=6,justify='center')
		self.note_ms['e1'].insert(10,self.note_ms['TAVG'])
		self.note_ms['e2'].insert(10,'0.5')
		self.note_ms['e3'].insert(10,'-0.5')
		self.note_ms['e1'].grid(row=5,column=6)
		self.note_ms['e2'].grid(row=6,column=6)
		self.note_ms['e3'].grid(row=7,column=6)

		b1 = tk.Button(self.home['page3'],  text="UPDATE", width = 4,command=lambda: self._update_plot(3))
		b1.grid(row=7, column=7,columnspan=2,sticky='w')			
		b1 = tk.Button(self.home['page3'],  text="CLEAR", width = 4,command=lambda: self._clf_plot(3))
		b1.grid(row=6, column=7,columnspan=2,sticky='w')	

		self.l1 = tk.Label(self.home['page3'], text="================= MSE CHANNEL OPT. ==================",justify='center')
		self.l1.grid(row=8,column=0,columnspan=9,pady=5)

		self.note_ms['frame2'] = tk.Frame(self.home['page3'],width=380)
		self.note_ms['frame2'].grid(row=9,column=0,columnspan=9)
		for k in range(10): self.note_ms['frame2'].grid_columnconfigure(k,weight=1)

		count = 0
		for i in range(3):
			for k in range(10):
				count = count + 1
				if count > self.note_in['nch']: break
				self.note_ms['c%i'%count] = tk.Checkbutton(self.note_ms['frame2'],variable=self.note_ms['msec%02i'%count],command=lambda: self._mse_check())
				self.note_ms['c%i'%count].grid(row=i,column=k,padx=7)				

		self.l1 = tk.Label(self.home['page3'], text="================== MSE WEIGHT OPT. ===================",justify='center')
		self.l1.grid(row=10,column=0,columnspan=9,pady=5)		

		self.note_ms['frame3'] = tk.Frame(self.home['page3'],width=380)
		self.note_ms['frame3'].grid(row=11,column=0,columnspan=9)
		for k in range(5): self.note_ms['frame3'].grid_columnconfigure(k,weight=1)

		count = 0
		for i in range(5):
			for k in range(5):
				count = count + 1
				if count > self.note_in['nch']: break
				self.note_ms['e%i'%(count+3)] = tk.Entry(self.note_ms['frame3'],width=9,justify='center')
				self.note_ms['e%i'%(count+3)].insert(10,self.note_ms['mses%02i'%count])
				self.note_ms['e%i'%(count+3)].grid(row=i,column=k,padx=3)

		self.l1 = tk.Label(self.home['page3'], text="USE EXP.MSE",justify='center')
		self.l1.grid(row=12,column=0,columnspan=2,pady=10)						
		self.note_ms['c1'] = tk.Checkbutton(self.home['page3'],variable=self.note_ms['use_exp_sgam'])
		self.note_ms['c1'].grid(row=12,column=2)

		self.l1 = tk.Label(self.home['page3'], text="SGAMMA",justify='center')
		self.l1.grid(row=12,column=4)	
		self.note_ms['e29'] = tk.Entry(self.home['page3'],width=7,justify='center')
		self.note_ms['e29'].insert(10,'0.02')
		self.note_ms['e29'].grid(row=12,column=5)

		b1 = tk.Button(self.home['page3'],  text="SET", width = 4,command=lambda: self._put_mse())
		b1.grid(row=12, column=6,sticky='w')			
		b1 = tk.Button(self.home['page3'],  text="RESET", width = 4,command=lambda: self._reset_mse())
		b1.grid(row=12, column=7,columnspan=2,sticky='w')			

		return

	def _mse_check(self):

		self._clf_plot(3)
		self._raw_mse_draw()
		return

	def _make_mse_plot(self):

		self.figure['name']['home'].canvas.draw_idle()
		self.figure['name']['home'].clf()
		self.figure['name']['homecursor'].visible = False
	
		ax1 = self.figure['name']['home'].add_subplot(1,1,1)

		ax1.cla();
		ax1.set_xlabel('R[m]')
		ax1.set_xlim(1.7,2.4)
		ax1.set_ylabel('T-GAMMA [a.u]')
		ax1.axhline(y=0,linestyle='--',c='gold')
		self.figure['name']['home'].tight_layout()

		for k in range(4): 
			self.figure['legend'][k] = []
			self.figure['pegend'][k] = []

		return

	def _raw_mse_draw(self):

		selection = self.note_ms['l1'].curselection()
		if len(selection) == 0: return
		line = self.note_ms['l1'].get(selection[0])
		time = float(line.split()[0])/1.e3

		self.figure['name']['home'].canvas.draw_idle()
		ax1 = self.figure['name']['home'].axes[0]

		if not self.mse_data['is']: return

		t = self.mse_data['time']
		dt = float(self.note_ms['e1'].get())
		tmin = time - 0.5*dt/1.e3
		tmax = time + 0.5*dt/1.e3
		ind1 = np.where(t>=tmin)
		ind2 = np.where(t[ind1]<=tmax)

		len2 = len(self.mse_data['ch'])
		RRR = np.zeros(len2)
		TGA = np.zeros(len2)
		SGA = np.zeros(len2)

		for k in range(len2):
			RRR[k] = self.mse_data['RRRGAM'][k]
			TGA[k] = np.mean(self.mse_data['TGAMMA%02i'%(k+1)][ind1][ind2])
			SGA[k] = np.mean(self.mse_data['SGAMMA%02i'%(k+1)][ind1][ind2])

		count = len(self.figure['pegend'][0]) + 1
		for k in range(len2):
			if self.note_ms['msec%02i'%(k+1)].get() == 1:
				line = ax1.errorbar(RRR[k],TGA[k],yerr=SGA[k],fmt='x',c='C%i'%count)
			else:
				ax1.errorbar(RRR[k],TGA[k],yerr=SGA[k],fmt='x',c='gray')
			ax1.text(RRR[k],TGA[k],str(k+1))

		self.figure['pegend'][0].append(line)
		self.figure['legend'][0].append('%i ms'%(time*1.e3))

		ax1.legend(self.figure['pegend'][0],self.figure['legend'][0])
		self._update_plot(3)
		self.figure['name']['home'].tight_layout()	
		return

	def _put_mse(self):

		line = self.note_ms['e29'].get()
		for k in range(self.note_in['nch']):
			self.note_ms['e%i'%(k+4)].delete(0,'end')
			self.note_ms['e%i'%(k+4)].insert(10,line)
		return

	def _reset_mse(self):

		selection = self.note_ms['l1'].curselection()
		if len(selection) == 0: return
		line = self.note_ms['l1'].get(selection[0])
		time = float(line.split()[0])/1.e3

		self.figure['name']['home'].canvas.draw_idle()
		ax1 = self.figure['name']['home'].axes[0]

		if not self.mse_data['is']: return

		t = self.mse_data['time']
		dt = float(self.note_ms['e1'].get())
		tmin = time - 0.5*dt/1.e3
		tmax = time + 0.5*dt/1.e3
		ind1 = np.where(t>=tmin)
		ind2 = np.where(t[ind1]<=tmax)

		len2 = len(self.mse_data['ch'])
		for k in range(len2):
			val = np.mean(self.mse_data['SGAMMA%02i'%(k+1)][ind1][ind2])
			self.note_ms['e%i'%(k+4)].delete(0,'end')
			self.note_ms['e%i'%(k+4)].insert(10,str(val))			
			
		return

	def _make_beam_frame(self):

		self.l1 = tk.Label(self.home['page4'], text="================= TIME SLICE INFO. ==================",justify='center')
		self.l1.grid(row=0,column=0,columnspan=9,pady=5)

		self.l1 = tk.Label(self.home['page4'], text="== FIT-DATA ==",justify='center')
		self.l1.grid(row=2,column=0,columnspan=5,pady=5,padx=45,sticky='w')
		self.l1 = tk.Label(self.home['page4'], text="== POST-NB ==",justify='center')
		self.l1.grid(row=2,column=4,columnspan=5,pady=5,padx=50,sticky='e')

		self.note_nb['frame1'] = tk.Frame(self.home['page4'])
		self.note_nb['frame1'].grid(row=3,column=0,columnspan=5,padx=20,sticky='w')
		self.note_nb['s1'] = tk.Scrollbar(self.note_nb['frame1'])
		self.note_nb['s1'].pack(side='right',fill='y')
		self.note_nb['l1'] = tk.Listbox(self.note_nb['frame1'],yscrollcommand = self.note_nb['s1'].set,height=7)
		self.note_nb['l1'].pack(side='left',fill='x')
		self.note_nb['s1']["command"] = self.note_nb['l1'].yview

		self.note_nb['frame2'] = tk.Frame(self.home['page4'])
		self.note_nb['frame2'].grid(row=3,column=4,columnspan=5,padx=20,sticky='e')
		self.note_nb['s2'] = tk.Scrollbar(self.note_nb['frame2'])
		self.note_nb['s2'].pack(side='right',fill='y')
		self.note_nb['l2'] = tk.Listbox(self.note_nb['frame2'],yscrollcommand = self.note_nb['s2'].set,height=7)
		self.note_nb['l2'].pack(side='left',fill='x')
		self.note_nb['s2']["command"] = self.note_nb['l2'].yview
		self.note_nb['l2'].bind('<Double-1>',lambda x: self._draw_beam())

		flag1 = ['PBEAM','INCUR','BECUR','BSCUR','TOTCUR','PHEAT','WFAST','INFRAC']
		self.note_nb['m1'] = tk.OptionMenu(self.home['page4'],self.note_nb['plottype'],*flag1)
		self.note_nb['m1'].grid(row=4,column=4,columnspan=3,pady=5,sticky='e')		
		self.note_nb['m1'].config(width=7)
		b1 = tk.Button(self.home['page4'],  text="CLEAR", width = 4,command=lambda: self._clf_plot(4))
		b1.grid(row=4, column=7,columnspan=2,sticky='w')			

		b1 = tk.Button(self.home['page4'],  text="UPDATE", width = 4,command=lambda: self._update_plot(4))
		b1.grid(row=5, column=7,columnspan=2,sticky='w')		
		
		self.l1 = tk.Label(self.home['page4'], text="YMAX",justify='center')
		self.l1.grid(row=5,column=3,padx=5,sticky='e')		
		self.note_nb['e1'] = tk.Entry(self.home['page4'],width=6,justify='center')
		self.note_nb['e1'].insert(10,'2.0')
		self.note_nb['e1'].grid(row=5, column=4)

		self.l1 = tk.Label(self.home['page4'], text="YMIN",justify='center')
		self.l1.grid(row=5,column=5,padx=5,sticky='e')		
		self.note_nb['e2'] = tk.Entry(self.home['page4'],width=6,justify='center')
		self.note_nb['e2'].insert(10,'0.0')
		self.note_nb['e2'].grid(row=5, column=6)		

		self.l1 = tk.Label(self.home['page4'], text="================= EQUIL RUN OPT. ==================",justify='center')
		self.l1.grid(row=6,column=0,columnspan=9,pady=5)

		self.l1 = tk.Label(self.home['page4'], text="GFILE",justify='center')
		self.l1.grid(row=7,column=0,padx=5,sticky='e')				
		self.note_nb['m2'] = tk.OptionMenu(self.home['page4'],self.note_nb['runeq'],*self.note_nb['gflag'])
		self.note_nb['m2'].grid(row=7,column=1,columnspan=2,pady=5,sticky='e')		
		self.note_nb['m2'].config(width=5)

		self.l1 = tk.Label(self.home['page4'], text="BSMODEL",justify='center')
		self.l1.grid(row=7,column=3,padx=5,sticky='e')				
		self.note_nb['m3'] = tk.OptionMenu(self.home['page4'],self.note_nb['BSMODEL'],*self.note_nb['bslist'])
		self.note_nb['m3'].grid(row=7,column=4,columnspan=2,pady=5,sticky='e')				
		self.note_nb['m3'].config(width=7)

		self.l1 = tk.Label(self.home['page4'], text="HAGCORE",justify='center')	
		self.l1.grid(row=7,column=6,sticky='e')		
		self.note_nb['c1'] = tk.Checkbutton(self.home['page4'],variable=self.note_nb['HAGCORE'])
		self.note_nb['c1'].grid(row=7,column=7)		

		self.note_nb['frame3'] = tk.Frame(self.home['page4'])
		self.note_nb['frame3'].grid(row=8,column=0,columnspan=9)
		flag1 = ['HAGCOREPSI','CORENEO','NSCALE','BSMULTI','NS','NT','MAPNS','MAPNT','IPCRIT','BSCRIT']
		flag2 = ['HAGPSI','CENEO','NSCALE','BSMULTI','NS [#]','NT [#]','MAPNS','MAPNT','IPCRIT','BSCRIT']

		count = 0
		for i in range(4):
			for j in range(3):
				count = count + 1
				if count > 10: break
				self.l1 = tk.Label(self.note_nb['frame3'], text=flag2[count-1],justify='center')
				self.l1.grid(row=i,column=2*j,padx=5)				
				self.note_nb['e%i'%(count+2)] = tk.Entry(self.note_nb['frame3'],width=6,justify='center')
				self.note_nb['e%i'%(count+2)].insert(10,self.note_nb[flag1[count-1]])
				self.note_nb['e%i'%(count+2)].grid(row=i,column=2*j+1,padx=5)			

		self.l1 = tk.Label(self.home['page4'], text="================= NUBEAM RUN OPT. ==================",justify='center')
		self.l1.grid(row=9,column=0,columnspan=9,pady=5)				

		self.note_nb['frame4'] = tk.Frame(self.home['page4'])
		self.note_nb['frame4'].grid(row=10,column=0,columnspan=9)
		flag1 = ['NPROC','RUNSTEP','RUNAVG','RUNDT','AVGDT','B-DIFF','MFILE','SFILE','IFILE']
		flag2 = ['NPROC','RUNSTEP','RUNAVG','RUNDT','AVGDT','B-DIFF','M-FILE','S-FILE','I-FILE']

		count = 0
		for i in range(2):
			for j in range(3):
				count = count + 1
				if count > 6: break
				self.l1 = tk.Label(self.note_nb['frame4'], text=flag2[count-1],justify='center')
				self.l1.grid(row=i,column=2*j,padx=5)				
				self.note_nb['e%i'%(count+12)] = tk.Entry(self.note_nb['frame4'],width=6,justify='center')
				self.note_nb['e%i'%(count+12)].insert(10,self.note_nb[flag1[count-1]])
				self.note_nb['e%i'%(count+12)].grid(row=i,column=2*j+1,padx=5)	

		count = 6
		for i in range(3):
			count = count + 1
			self.l1 = tk.Label(self.note_nb['frame4'], text=flag2[count-1],justify='center')
			self.l1.grid(row=i+3,column=0,padx=5,sticky='w')				
			self.note_nb['e%i'%(count+12)] = tk.Entry(self.note_nb['frame4'],width=35,justify='center')
			self.note_nb['e%i'%(count+12)].insert(10,self.note_nb[flag1[count-1]])
			self.note_nb['e%i'%(count+12)].grid(row=i+3,column=1,columnspan=4,padx=5,sticky='w')
			b1 = tk.Button(self.note_nb['frame4'],  text="OPEN", width = 4,command=lambda: self._ask_file_name(count+13))
			b1.grid(row=i+3, column=5,sticky='w')			

		b1 = tk.Button(self.home['page4'],  text="RUN BEAM", width = 10,command=lambda: self._run_beam_select())
		b1.grid(row=11, column=0,columnspan=3,pady=5)
		b1 = tk.Button(self.home['page4'],  text="RUN BEAM A", width = 15,command=lambda: self._run_beams())
		b1.grid(row=11, column=3,columnspan=3,pady=5)			
		b1 = tk.Button(self.home['page4'],  text="UPATE RUN", width = 10,command=lambda: self._load_beam())
		b1.grid(row=11, column=6,columnspan=3,pady=5)	
		return

	def _make_beam_plot(self):

		self.figure['name']['home'].canvas.draw_idle()
		self.figure['name']['home'].clf()
		self.figure['name']['homecursor'].visible = False
	
		ax1 = self.figure['name']['home'].add_subplot(1,1,1)

		ax1.cla();
		ax1.set_xlabel('$\\psi_N$')
		ax1.set_xlim(0.,1.1)
		ax1.set_ylabel('[a.u]')
		self.figure['name']['home'].tight_layout()
		for k in range(4): 
			self.figure['legend'][k] = []
			self.figure['pegend'][k] = []

		return

	def _find_beam_run(self,runeq=''):

		tlist = list(self.note_nb['l1'].get(0,'end'))
		self.note_nb['l2'].delete(0,'end')
		self.bdata = dict()
		for time in tlist:
			tt = int(float(time.split()[0]))
			dirs = 'EQUIL/%s/%i/'%(runeq,tt)
			zdfile = dirs+'nubeam_iter_result'
			if os.path.isfile(zdfile):
				self.bdata[tt] = dict()
				self.note_nb['l2'].insert('end',time)
		return

	def _load_beam(self,runeq=''):

		self._find_beam_run(self.note_nb['runeq'].get())	
		tlist = list(self.note_nb['l2'].get(0,'end'))
		if len(tlist) == 0: return
		for key in self.bdata.keys():
			if runeq == '': dirs = 'EQUIL/%s/%i/'%(self.note_nb['runeq'].get(),key)
			else:	dirs = 'EQUIL/%s/%i/'%(runeq,key)
			f = open(dirs+'NUBEAM/nubeam_out0d','r')
			self.bdata[key]['pe'] = float(f.readline().split()[-1])
			self.bdata[key]['pi'] = float(f.readline().split()[-1])
			self.bdata[key]['ptot'] = float(f.readline().split()[-1])
			self.bdata[key]['wfast'] = float(f.readline().split()[-1])
			self.bdata[key]['bbn'] = float(f.readline().split()[-1])
			self.bdata[key]['btn'] = float(f.readline().split()[-1])
			self.bdata[key]['totn'] = float(f.readline().split()[-1])
			f.close()
			f = open(dirs+'pre_prof','r')
			datn = int(float(f.readline()))
			self.bdata[key]['prest'] = np.zeros((datn,3))
			for i in range(datn):
				line = f.readline().split()
				self.bdata[key]['prest'][i,0] = float(line[0])
				self.bdata[key]['prest'][i,1] = float(line[1])
				self.bdata[key]['prest'][i,2] = float(line[2])
			f.close()

			f = open(dirs+'OUTPUT/EFIT_JCONST','r')
			line = f.readline()
			datn = int(float(f.readline()))
			self.bdata[key]['jconst'] = np.zeros((datn,3))
			for i in range(datn):
				line = f.readline().split()
				self.bdata[key]['jconst'][i,0] = float(line[0])
				self.bdata[key]['jconst'][i,1] = float(line[2])
				self.bdata[key]['jconst'][i,2] = float(line[3])
			f.close()			

			f = open(dirs+'nubeam_iter_result','r')
			for i in range(3): line = f.readline().split()
			f.close()
			
			self.bdata[key]['wth'] = float(line[1])
			self.bdata[key]['wfast'] = float(line[2])
			self.bdata[key]['wmhd'] = float(line[3])
			self.bdata[key]['wdia'] = float(line[4])
			self.bdata[key]['ifast'] = float(line[5])*1.e6
			self.bdata[key]['vloop']= float(line[6])

			self.bdata[key]['currt'] = np.zeros((datn,5))
			f = open(dirs+'CHEASE/curr_prof','r')
			line = f.readline().split()
			for i in range(datn): 
				line = f.readline().split()
				self.bdata[key]['currt'][i,0] = float(line[0])
				self.bdata[key]['currt'][i,1] = float(line[1])
				self.bdata[key]['currt'][i,2] = float(line[3])
				self.bdata[key]['currt'][i,3] = float(line[4])
				self.bdata[key]['currt'][i,4] = float(line[5])
			f.close()

			f = open(dirs+'OUTPUT/BS_Profile','r')
			for i in range(datn+2): line = f.readline().split()
			self.bdata[key]['ipbs']   =  float(line[-2])
			line = f.readline().split()
			self.bdata[key]['iptot']  =  float(line[-2])
			f.close()			

		return

	def _run_beams(self):

		isgfile = self._check_gfiles(self.note_nb['runeq'].get())
		if not isgfile: return

		tlist = list(self.note_nb['l1'].get(0,'end'))
		for time in tlist:
			tt = float(time.split()[0])
			self._run_beam(tt)
		self._load_beam()
		return

	def _run_beam_select(self):

		isgfile = self._check_gfiles(self.note_nb['runeq'].get())
		if not isgfile: return
		time = self.note_nb['l1'].curselection()
		if len(time) == 0: return
		time = self.note_nb['l1'].get(time[0])
		
		tt = float(time.split()[0])
		self._run_beam(tt)
		self._load_beam()

		return

	def _run_beam(self,ttime=0):

		currdir = os.getcwd()
		flags = ['1A','1B','1C','2A','2B','2C']
		powers   = np.zeros(6)
		energies = np.zeros(6)

		delt = 100
		tmax = (ttime+0.5*delt)/1.e3
		tmin = (ttime-0.5*delt)/1.e3	

		for k in range(len(flags)):
			flag = flags[k]
			t = self.hcd[flag][0]
			if len(t) >0.:
				ind1 = np.where(t>=tmin)
				ind2 = np.where(t[ind1]<=tmax)
				pw = np.mean(self.hcd[flag][1][ind1][ind2])
				en = np.mean(self.hcd[flag][3][ind1][ind2])
			else:
				pw = 0.
				en = 0.

			if pw < 0.1: pw= 0.;
			powers[k] = pw; energies[k] = en;

		t = self.etc['wdia'][0]
		if len(t) >0.:
			ind1 = np.where(t>=tmin)
			ind2 = np.where(t[ind1]<=tmax)
			wdia = np.mean(self.etc['wdia'][1][ind1][ind2])
		else:
			wdia = 0.

		self.nbeam = 0.;
		powerl = '%13.7e'%(powers[0]*1.e6)
		energyl = '%13.7e'%energies[0]
		for k in range(len(flags)):
			if powers[k] > 0.: self.nbeam = self.nbeam + 1
			if k>0:
				powerl = powerl +',%13.7e'%(powers[k]*1.e6)
				energyl = energyl +',%13.7e'%energies[k]

		if self.nbeam == 0.:	
			print('>>> No beam...')
			return
		self.nbeam = 3
		runid = self.note_nb['runeq'].get()
		shotn = int(float(self.note_in['e1'].get()))

		run_dir = currdir + '/EQUIL/%s/%i'%(runid,ttime)
		nubeam_dir = run_dir + '/NUBEAM'
		profile_dir = run_dir + '/PROFILES'
		vt_dir = currdir+'/PROFILES/MPROFILES_%s/%i/VT_fit.dat'%(runid,ttime)
		eq_dir = currdir+'/EFITS/%s/g%06i.%06i'%(runid,shotn,ttime)
		p_dir  = currdir+'/PROFILES/MPROFILES_%s/%i/chease_kinprof_fit'%(runid,ttime)

		out_dir1 = run_dir + '/OUTPUT'
		out_dir2 = run_dir + '/PROFILES'
		out_file1= run_dir + '/nubeam_iter_result'
		try:	os.mkdir(currdir + '/EQUIL/%s'%runid)
		except:	pass	
		try:	os.mkdir(run_dir)
		except:	pass	
		try: os.remove(out_file1)
		except: pass

		os.chdir(run_dir)
		self._make_chease_opt('chease_opt',vt_file=vt_dir,eq_file=eq_dir,p_file=p_dir,wdia=wdia)
		self._make_nubeam_opt('nubeam_opt',eqtype=1,ttime=ttime,bpower=powerl,benergy=energyl)

		try:	rmtree(nubeam_dir)
		except:	pass
		os.mkdir(nubeam_dir)
		try:	os.mkdir(save_dir)
		except:	pass
		try:	os.mkdir(profile_dir)
		except:	pass

		filename = run_dir+'/chease_batch'
		log_e = run_dir+'/chease_batch.e'
		log_o = run_dir+'/chease_batch.o'
		command = 'cd ' + run_dir + '\n'	
		command = command + nubeam_dir2 +' %s %s %s \n'%(eq_dir,'f',self.note_nb['e18'].get())

		make_batch_script(filename,None,log_e,log_o,command,'CHEASE_EMSE_%s_%i'%(self.note_in['e1'].get(),ttime))	
		runid = submit_batch_script(filename)
		os.chdir(currdir)

		return

	def _make_chease_opt(self,filename='',vt_file='',eq_file='',p_file='',wdia=0.):

		f = open(filename,'w')

		f.write('!-- run type \n')

		if self.nbeam > 0:	f.write('RUN_MODE = NUBEAM \n')
		else:	f.write('RUN_MODE = NORMAL \n')

		f.write('!-- input type \n')
		f.write('EQDSK = %s\n'%eq_file)
		f.write('kinetic_profile_type = 1 \n')
		f.write('chease_kinetic_file = %s \n'%p_file)
		f.write('ne_file =  \n')
		f.write('te_file =  \n')
		f.write('ti_file =  \n')
		f.write('vt_file = %s \n'%vt_file)
		if self.nbeam > 0:
			f.write('USE_EXT_P = True \n')
			f.write('USE_EXT_J = True \n')
		else:
			f.write('USE_EXT_P = False \n')
			f.write('USE_EXT_J = False \n')		

		f.write('ZEFF  = %s \n'%self.note_pr['e14'].get())
		f.write('ZIMP  = 6.0 \n')
		f.write('WDIA = %s \n'%wdia)
		f.write('APF = 0.0 \n')
		f.write('BND_PSIN = 0.995 \n')
		f.write('Beta_criterion_type = 0 \n')
		f.write('!-- run option (EXT_VLOOP = 0.0 then automatic run)\n')
		f.write('VLOOP_MOD = 1.0 \n')
		f.write('EXT_VLOOP = 0.0 \n')
		f.write('Current_ITERN = 25 \n')
		f.write('NUBEAM_ITERN = 5 \n')
		f.write('RELAX = 0.6 \n')
		if (self.note_nb['BSMODEL'].get().lower() == 'nhager' or self.note_nb['BSMODEL'].get().lower() == 'neo'):
			f.write('USE_NEO = True \n')
			f.write('USE_HAGER = False \n')
			f.write('USE_CHANG = False \n')
		elif self.note_nb['BSMODEL'].get().lower() == 'hager':
			f.write('USE_NEO = False \n')
			f.write('USE_HAGER = True \n')
			f.write('USE_CHANG = False \n')
		elif self.note_nb['BSMODEL'].get().lower() == 'csauter':
			f.write('USE_NEO = False \n')
			f.write('USE_HAGER = False \n')
			f.write('USE_CHANG = True \n')
		else:	
			f.write('USE_NEO = False \n')
			f.write('USE_HAGER = False \n')
			f.write('USE_CHANG = False \n')
		if self.note_nb['HAGCORE'].get() == 1:	f.write('HAG_CORE_MOD=True \n')
		else:	f.write('HAG_CORE_MOD=False \n')
		f.write('HAG_CORE_MOD_PSIN = %s \n'%self.note_nb['e3'].get())
		f.write('Core_neo = %s \n'%self.note_nb['e4'].get())
		f.write('DENSITY_SCALE = %s \n'%self.note_nb['e5'].get())
		f.write('BSMULTI = %s\n'%self.note_nb['e6'].get())
		f.write('NS = %s\n'%self.note_nb['e7'].get())
		f.write('NT = %s\n'%self.note_nb['e8'].get())
		f.write('MAP_NS = %s\n'%self.note_nb['e9'].get())
		f.write('MAP_NT = %s\n'%self.note_nb['e10'].get())
		f.write('IP_CRIT = %s\n'%self.note_nb['e11'].get())
		f.write('BS_CRIT = %s\n'%self.note_nb['e12'].get())
		f.write('!-- Output option \n')
		f.write('EFIT_CONST = True \n')
		f.write('ADJUST_PROF = False\n')

		f.close()
		return	

	def _make_nubeam_opt(self,filename='',eqtype=1,ttime=0,bpower='',benergy=''):

		f = open(filename,'w')
		f.write('ZEFF = %s\n'%self.note_pr['e14'].get())
		f.write('NBEAM = %i\n'%self.nbeam)
		f.write('BPOWER = %s\n'%bpower)
		f.write('BENERGY = %s\n'%benergy)
		f.write('DIFFUSIVITY = %s\n'%self.note_nb['e18'].get())
		f.write('SHOT = %s\n'%self.note_in['e1'].get())
		f.write('TIME = %s\n'%int(ttime))
		f.write('eqdsk = geqdsk\n')
		f.write('NPROC = %s\n'%self.note_nb['e13'].get())
		f.write('RUN_STEP = %s\n'%self.note_nb['e14'].get())
		f.write('RUN_AVG = %s\n'%self.note_nb['e15'].get())
		f.write('RUN_DT = %s\n'%self.note_nb['e16'].get())
		f.write('AVG_DT = %s\n'%self.note_nb['e17'].get())
		if not self.note_nb['e19'].get() == '':	f.write('MFILE = %s\n'%self.note_nb['e19'].get())
		if not self.note_nb['e20'].get() == '':	f.write('SFILE = %s\n'%self.note_nb['e20'].get())
		if not self.note_nb['e21'].get() == '':	f.write('IFILE = %s\n'%self.note_nb['e21'].get())

		f.close()	

		return

	def _draw_beam(self):

		selection = self.note_nb['l2'].curselection()
		if len(selection) == 0: return
		line = self.note_nb['l2'].get(selection[0])
		time = int(line.split()[0])

		tlist = list(self.note_nb['l2'].get(0,'end'))
		tmin = float(tlist[0].split()[0])
		tmax = float(tlist[-1].split()[0])

		self.figure['name']['home'].canvas.draw_idle()
		ax = self.figure['name']['home'].axes[0]
		plottype = self.note_nb['plottype'].get()

		if   plottype == 'PBEAM' : plotspace = 'r'
		elif plottype == 'INCUR' : plotspace = 'r'
		elif plottype == 'BECUR' : plotspace = 'r'
		elif plottype == 'BSCUR' : plotspace = 'r'
		elif plottype == 'TOTCUR': plotspace = 'r'
		elif plottype == 'PHEAT' : plotspace = 't'
		elif plottype == 'WFAST' : plotspace = 't'
		elif plottype == 'INFRAC': plotspace = 't'

		self.note_nb['plotspace'] = plotspace
		if not plottype == self.note_nb['plottype_old']:
			self._clf_plot(4)
			self.figure['pegend'][0] = []
			self.figure['legend'][0] = []
			self.note_nb['plottype_old'] = plottype
			ymax = 0.
			self.note_nb['e1'].delete(0,'end')
			self.note_nb['e1'].insert(10,'0.')
			if plotspace =='r': ax.set_xlabel('$\\psi_N$')
			else: ax.set_xlabel('times [ms]')

		count = len(self.figure['legend'][0]) + 1

		times = []
		for key in self.bdata.keys(): times.append(key)
		times = np.array(times,dtype='float')
		times = np.sort(times)

		if plotspace == 'r': self.note_nb['xmin'] = 0.; self.note_nb['xmax'] = 1.
		else: self.note_nb['xmin'] = 0.7*min(times); self.note_nb['xmax'] = 1.3*max(times);
		
		ax.set_xlim(self.note_nb['xmin'],self.note_nb['xmax'])

		vals = np.copy(times)
		data = self.bdata[time]

		if   plottype == 'PBEAM' : 
			line, = ax.plot(data['prest'][:,0],data['prest'][:,2]/1.e3,c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $P_{BEAM}$-%s'%(time,self.note_nb['runeq'].get()))
			self.note_nb['ymax'] = 1.1 * max(data['prest'][:,2]) / 1.e3
			line2 = '[kPa]'
		elif plottype == 'INCUR' :
			line, = ax.plot(data['currt'][:,0],data['currt'][:,1]/1.e6,c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $J_{IND}$-%s'%(time,self.note_nb['runeq'].get()))
			self.note_nb['ymax'] = 1.1 * max(data['currt'][:,1]) / 1.e6
			line2 = '[MA/$m^2$]'
		elif plottype == 'BECUR' : 
			line, = ax.plot(data['currt'][:,0],data['currt'][:,3]/1.e6,c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $J_{BEAM}$-%s'%(time,self.note_nb['runeq'].get()))		
			self.note_nb['ymax'] = 1.1 * max(data['currt'][:,2]) / 1.e6
			line2 = '[MA/$m^2$]'
		elif plottype == 'BSCUR' : 
			line, = ax.plot(data['currt'][:,0],data['currt'][:,2]/1.e6,c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $J_{BS}$-%s'%(time,self.note_nb['runeq'].get()))		
			self.note_nb['ymax'] = 1.1 * max(data['currt'][:,3]) / 1.e6
			line2 = '[MA/$m^2$]'
		elif plottype == 'TOTCUR':
			line, = ax.plot(data['currt'][:,0],data['currt'][:,4]/1.e6,c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $J_{TOT}$-%s'%(time,self.note_nb['runeq'].get()))		
			self.note_nb['ymax'] = 1.1 * max(data['currt'][:,4]) / 1.e6
			line2 = '[MA/$m^2$]'
		elif plottype == 'PHEAT' : 
			for k in range(len(times)): vals[k] = self.bdata[int(times[k])]['ptot']/1.e6
			line, = ax.plot(times,vals,'-x',c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $P_{HEAT}$-%s'%(time,self.note_nb['runeq'].get()))		
			self.note_nb['ymax'] = 1.2 * max(vals)
			line2 = '[MW]'
		elif plottype == 'WFAST' : 
			for k in range(len(times)): vals[k] = self.bdata[int(times[k])]['wfast']
			line, = ax.plot(times,vals,'-x',c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $W_{FAST}$-%s'%(time,self.note_nb['runeq'].get()))		
			self.note_nb['ymax'] = 1.2 * max(vals)
			line2 = '[kJ]'
		elif plottype == 'INFRAC': 
			for k in range(len(times)): vals[k] = (self.bdata[int(times[k])]['ipbs']+self.bdata[int(times[k])]['ifast'])/self.bdata[int(times[k])]['iptot']
			line, = ax.plot(times,vals,'-x',c='C%i'%count)		
			self.figure['legend'][0].append('%i ms Frac -%s'%(time,self.note_nb['runeq'].get()))
			self.note_nb['ymax'] = 1.2 * max(vals)
			line2 = '[%]'
		self.figure['pegend'][0].append(line)

		ymax = float(self.note_nb['e1'].get())
		if ymax < self.note_nb['ymax']:
			self.note_nb['e1'].delete(0,'end')
			self.note_nb['e1'].insert(10,'%4.2f'%self.note_nb['ymax'])

		ax.set_ylim(self.note_nb['ymin'],self.note_nb['ymax'])
		ax.set_ylabel(line2)
		ax.legend(self.figure['pegend'][0],self.figure['legend'][0])
		self.figure['name']['home'].tight_layout()			
		return

	def _make_efit_frame(self):

		self.l1 = tk.Label(self.home['page5'], text="================= TIME SLICE INFO. ==================",justify='center')
		self.l1.grid(row=0,column=0,columnspan=9,pady=5)

		self.l1 = tk.Label(self.home['page5'], text="== PREPARED ==",justify='center')
		self.l1.grid(row=2,column=0,columnspan=5,pady=5,padx=45,sticky='w')
		self.l1 = tk.Label(self.home['page5'], text="== POST-EFIT ==",justify='center')
		self.l1.grid(row=2,column=4,columnspan=5,pady=5,padx=50,sticky='e')

		self.note_ef['frame1'] = tk.Frame(self.home['page5'])
		self.note_ef['frame1'].grid(row=3,column=0,columnspan=5,padx=20,sticky='w')
		self.note_ef['s1'] = tk.Scrollbar(self.note_ef['frame1'])
		self.note_ef['s1'].pack(side='right',fill='y')
		self.note_ef['l1'] = tk.Listbox(self.note_ef['frame1'],yscrollcommand = self.note_ef['s1'].set,height=7)
		self.note_ef['l1'].pack(side='left',fill='x')
		self.note_ef['s1']["command"] = self.note_ef['l1'].yview
		self.note_ef['l1'].bind('<Double-1>',lambda x: self._get_kfile())

		self.note_ef['frame2'] = tk.Frame(self.home['page5'])
		self.note_ef['frame2'].grid(row=3,column=4,columnspan=5,padx=20,sticky='e')
		self.note_ef['s2'] = tk.Scrollbar(self.note_ef['frame2'])
		self.note_ef['s2'].pack(side='right',fill='y')
		self.note_ef['l2'] = tk.Listbox(self.note_ef['frame2'],yscrollcommand = self.note_ef['s2'].set,height=7)
		self.note_ef['l2'].pack(side='left',fill='x')
		self.note_ef['s2']["command"] = self.note_ef['l2'].yview
		self.note_ef['l2'].bind('<Double-1>',lambda x: self._draw_efit())

		flag1 = ['PRESS','JTOR','QPROF','MSE','MCOIL','FLOOP','PSCALE','CHISQR','WMHD','BPOL','LI']
		self.note_ef['m1'] = tk.OptionMenu(self.home['page5'],self.note_ef['plottype'],*flag1)
		self.note_ef['m1'].grid(row=4,column=4,columnspan=3,pady=5,sticky='e')		
		self.note_ef['m1'].config(width=7)
		b1 = tk.Button(self.home['page5'],  text="CLEAR", width = 4,command=lambda: self._clf_plot(5))
		b1.grid(row=4, column=7,columnspan=2,sticky='w')			

		b1 = tk.Button(self.home['page5'],  text="UPDATE", width = 4,command=lambda: self._update_plot(5))
		b1.grid(row=5, column=7,columnspan=2,sticky='w')		
		
		self.l1 = tk.Label(self.home['page5'], text="YMAX",justify='center')
		self.l1.grid(row=5,column=3,padx=5,sticky='e')		
		self.note_ef['e1'] = tk.Entry(self.home['page5'],width=6,justify='center')
		self.note_ef['e1'].insert(10,'2.0')
		self.note_ef['e1'].grid(row=5, column=4)

		self.l1 = tk.Label(self.home['page5'], text="YMIN",justify='center')
		self.l1.grid(row=5,column=5,padx=5,sticky='e')		
		self.note_ef['e2'] = tk.Entry(self.home['page5'],width=6,justify='center')
		self.note_ef['e2'].insert(10,'0.0')
		self.note_ef['e2'].grid(row=5, column=6)		

		self.l1 = tk.Label(self.home['page5'], text="================= EFIT TYPE OPT. ==================",justify='center')
		self.l1.grid(row=6,column=0,columnspan=9,pady=5)

		self.l1 = tk.Label(self.home['page5'], text="GFILE",justify='center')
		self.l1.grid(row=7,column=0,padx=5,sticky='e')				
		self.note_ef['m2'] = tk.OptionMenu(self.home['page5'],self.note_ef['runeq'],*self.note_nb['gflag'])
		self.note_ef['m2'].grid(row=7,column=1,columnspan=3,pady=5,sticky='w')		
		self.note_ef['m2'].config(width=5)

		b1 = tk.Button(self.home['page5'],  text="COIL CONFIGURE", width = 25,command=lambda: self._custom_field_config())
		b1.grid(row=7, column=4,columnspan=4)	

		self.l1 = tk.Label(self.home['page5'], text="================= EFIT KNOTS OPT. ==================",justify='center')
		self.l1.grid(row=8,column=0,columnspan=9,pady=5)

		self.note_ef['frame3'] = tk.Frame(self.home['page5'])
		self.note_ef['frame3'].grid(row=9,column=0,columnspan=9)
		flag1 = ['PCNN','PENN','FCNN','FENN','JCNN','JENN','JSTART']
		flag2 = ['PCNN','PENN','FCNN','FENN','JCNN','JENN','JPSI']

		count = 0
		for i in range(2):
			for j in range(4):
				count = count + 1
				if count > 7: break
				self.l1 = tk.Label(self.note_ef['frame3'], text=flag2[count-1],justify='center')
				self.l1.grid(row=i,column=2*j,padx=2)				
				self.note_ef['e%i'%(count+2)] = tk.Entry(self.note_ef['frame3'],width=5,justify='center')
				self.note_ef['e%i'%(count+2)].insert(10,self.note_ef[flag1[count-1]])
				self.note_ef['e%i'%(count+2)].grid(row=i,column=2*j+1,padx=1)	

		self.l1 = tk.Label(self.note_ef['frame3'], text="HMODE",justify='center')
		self.l1.grid(row=1,column=6,padx=3,sticky='e')	
		self.note_ef['c1'] = tk.Checkbutton(self.note_ef['frame3'],variable=self.note_ef['HMODE'])
		self.note_ef['c1'].grid(row=1,column=7)	
		#self.note_ef['c1'].config(state='disabled')

		self.l1 = tk.Label(self.home['page5'], text="================= EFIT CONSTRAINTS ==================",justify='center')
		self.l1.grid(row=10,column=0,columnspan=9,pady=5)		

		self.note_ef['frame4'] = tk.Frame(self.home['page5'])
		self.note_ef['frame4'].grid(row=11,column=0,columnspan=9)		
		flag1 = ['RITER','MAXITER','RELAX','CONVERG','FWTCUR','SUPP0','QCONST']
		flag2 = ['RITER','MITER','RELAX','CONV','FWTCUR','SUPP0','Q0']

		count = 0
		for i in range(2):
			for j in range(4):
				count = count + 1
				if count > 7: break
				self.l1 = tk.Label(self.note_ef['frame4'], text=flag2[count-1],justify='center')
				self.l1.grid(row=i,column=2*j,padx=2)				
				self.note_ef['e%i'%(count+9)] = tk.Entry(self.note_ef['frame4'],width=5,justify='center')
				self.note_ef['e%i'%(count+9)].insert(10,self.note_ef[flag1[count-1]])
				self.note_ef['e%i'%(count+9)].grid(row=i,column=2*j+1,padx=1)	

		self.l1 = tk.Label(self.note_ef['frame4'], text="Q0FIX",justify='center')
		self.l1.grid(row=1,column=6,padx=3,sticky='e')	
		self.note_ef['c2'] = tk.Checkbutton(self.note_ef['frame4'],variable=self.note_ef['QFIX'])
		self.note_ef['c2'].grid(row=1,column=7)			

		self.l1 = tk.Label(self.note_ef['frame4'], text="JCONST",justify='center')
		self.l1.grid(row=2,column=0,padx=2,sticky='e')	
		self.note_ef['c3'] = tk.Checkbutton(self.note_ef['frame4'],variable=self.note_ef['JCONST'])
		self.note_ef['c3'].grid(row=2,column=1)		
		self.note_ef['c3'].config(state='disabled')

		self.l1 = tk.Label(self.note_ef['frame4'], text="PSCALE",justify='center')
		self.l1.grid(row=2,column=2,padx=2,sticky='e')	
		self.note_ef['c4'] = tk.Checkbutton(self.note_ef['frame4'],variable=self.note_ef['PSCALE'])
		self.note_ef['c4'].grid(row=2,column=3)		

		self.l1 = tk.Label(self.note_ef['frame4'], text="PSCALE2",justify='center')
		self.l1.grid(row=2,column=4,padx=2,sticky='e')	
		self.note_ef['c4'] = tk.Checkbutton(self.note_ef['frame4'],variable=self.note_ef['PSCALE2'])
		self.note_ef['c4'].grid(row=2,column=5)				

		self.l1 = tk.Label(self.note_ef['frame4'], text="ERMSE",justify='center')
		self.l1.grid(row=2,column=6,padx=2,sticky='e')	
		self.note_ef['c6'] = tk.Checkbutton(self.note_ef['frame4'],variable=self.note_ef['ERMSE'])
		self.note_ef['c6'].grid(row=2,column=7)	

		b1 = tk.Button(self.home['page5'],  text="RUN", width = 5,command=lambda: self._run_efits())
		b1.grid(row=12, column=2,columnspan=2,pady=5)

		b1 = tk.Button(self.home['page5'],  text="RUN ALL", width = 10,command=lambda: self._run_efitm())
		b1.grid(row=12, column=4,columnspan=3,pady=5)			
		b1 = tk.Button(self.home['page5'],  text="UPATE RUN", width = 10,command=lambda: self._find_efit_run(load_data=True))
		b1.grid(row=12, column=7,columnspan=3,pady=5)			

		self.l1 = tk.Label(self.home['page5'], text="NW",justify='center')
		self.l1.grid(row=12,column=0,sticky='w')	
		self.note_ef['m3'] = tk.OptionMenu(self.home['page5'],self.note_ef['efitgrid'],'65','129','257')
		self.note_ef['m3'].grid(row=12,column=0,columnspan=3,pady=5,padx=25,sticky='w')		
		self.note_ef['m3'].config(width=2)		

		return

	def _make_efit_plot(self):

		self.figure['name']['home'].canvas.draw_idle()
		self.figure['name']['home'].clf()
		self.figure['name']['homecursor'].visible = False
	
		ax1 = self.figure['name']['home'].add_subplot(1,1,1)

		ax1.cla();
		ax1.set_xlabel('$\\psi_N$')
		ax1.set_xlim(0.,1.1)
		ax1.set_ylabel('[a.u]')
		self.figure['name']['home'].tight_layout()
		for k in range(4): 
			self.figure['legend'][k] = []
			self.figure['pegend'][k] = []

		return

	def _check_efit_input(self):

		if (self.note_ef['runeq'].get() == self.note_nb['runeq_old']):
			self.note_nb['runeq_old'] = self.note_ef['runeq'].get()
			self._load_beam(self.note_ef['runeq'].get())
		tlist = list(self.note_nb['l2'].get(0,'end'))
		self.note_ef['l1'].delete(0,'end')
		for time in tlist: self.note_ef['l1'].insert('end',time)
		return

	def _make_efit_mse(self,ttime=0):

		time = ttime/1.e3
		t = self.mse_data['time']
		dt = float(self.note_ms['e1'].get())
		tmin = time - 0.5*dt/1.e3
		tmax = time + 0.5*dt/1.e3
		ind1 = np.where(t>=tmin)
		ind2 = np.where(t[ind1]<=tmax)

		len2 = len(self.mse_data['ch'])
		mse_dat = dict()
		mse_dat['rrrgam']  = np.zeros(len2)
		mse_dat['tgamma0'] = np.zeros(len2)
		mse_dat['tgamma']  = np.zeros(len2)
		mse_dat['sgamma']  = np.zeros(len2)
		mse_dat['fwtgam']  = np.zeros(len2)
		mse_dat['rrrgam']  = np.zeros(len2)
		mse_dat['zzzgam']  = np.zeros(len2)
		mse_dat['aa1gam']  = np.zeros(len2)
		mse_dat['aa2gam']  = np.zeros(len2)
		mse_dat['aa3gam']  = np.zeros(len2)
		mse_dat['aa4gam']  = np.zeros(len2)
		mse_dat['aa5gam']  = np.zeros(len2)
		mse_dat['aa6gam']  = np.zeros(len2)

		for k in range(len2):
			mse_dat['rrrgam'][k]  = self.mse_data['RRRGAM'][k]
			mse_dat['zzzgam'][k]  = self.mse_data['ZZZGAM'][k]
			mse_dat['tgamma0'][k] = np.mean(self.mse_data['TGAMMA%02i'%(k+1)][ind1][ind2])
			mse_dat['tgamma'][k]  = np.mean(self.mse_data['TGAMMA%02i'%(k+1)][ind1][ind2])
			if self.note_ms['use_exp_sgam'].get() == 1: mse_dat['sgamma'][k]  = np.mean(self.mse_data['SGAMMA%02i'%(k+1)][ind1][ind2])
			else: mse_dat['sgamma'][k] = float(self.note_ms['e%i'%(k+4)].get())
			mse_dat['fwtgam'][k]  = self.note_ms['msec%02i'%(k+1)].get()

			for j in range(1,7):
				mse_dat['aa%igam'%j][k] = self.mse_data['AA%iGAM'%j][k]

		flags = ['1A','1B','1C','2A','2B','2C']
		powers   = np.zeros(6)
		energies = np.zeros(6)

		delt = 100
		tmax = (ttime+0.5*delt)/1.e3
		tmin = (ttime-0.5*delt)/1.e3	

		for k in range(len(flags)):
			flag = flags[k]
			t = self.hcd[flag][0]
			if len(t) >0.:
				ind1 = np.where(t>=tmin)
				ind2 = np.where(t[ind1]<=tmax)
				pw = np.mean(self.hcd[flag][1][ind1][ind2])
				en = np.mean(self.hcd[flag][3][ind1][ind2])
			else:
				pw = 0.
				en = 0.

			if pw < 0.1: pw= 0.;
			powers[k] = pw; energies[k] = en;				

		return mse_dat, powers, energies

	def _make_pf_knots(self,ttime=0):
		
		kn = knots_tool3.knots()

		kn.spline_start = 0.2		#PRES
		kn.end_knot     = 0.985
		kn.start_knot   = 0.3
		kn.knots_shift  = 0.05
		kn.coren        = int(float(self.note_ef['e3'].get()))
		kn.edgen        = int(float(self.note_ef['e4'].get()))				
		kn.input_datan  = 501
		kn.delmin_core  = 0.1
		kn.delmin_edge  = 0.01
		kn.minloc = 0.
	
		psin = np.copy(self.bdata[ttime]['prest'][:,0])
		pres = np.copy(self.bdata[ttime]['prest'][:,1])
		presf = interp1d(psin,pres)
		pp   = np.copy(pres)
		deps = 1.e-5
		for k in range(1,len(psin)-1):
			psin2 = psin[k]+deps; psin1 = psin[k]-deps;
			pp[k] = (presf(psin2)-presf(psin1)) /2. / deps

		if self.note_ef['HMODE'].get() == 1:	pknot = kn.lmfit_fit(psin,pp/min(pp),2,True)
		else:	pknot = kn.lmfit_fit(psin,pp/min(pp),1,True)
		
		if self.note_ef['HMODE'].get() == 1:
			pp = pp * np.sign(pp[-5])
			ind = np.where(psin>0.8)
			loc = np.argmax(pp[ind])
			ind1 = np.where(psin[ind]<psin[ind][loc])
			loc = np.argmin(pp[ind][ind1])
			minloc = min(psin[ind][ind1][loc]+0.03,0.92)
		else:
			minloc = 1.
		
		kn.spline_start = 0.2			#FFp
		kn.end_knot     = 0.985
		kn.start_knot   = 0.3
		kn.knots_shift  = 0.05
		kn.coren        = int(float(self.note_ef['e5'].get()))
		kn.edgen        = int(float(self.note_ef['e6'].get()))
		kn.input_datan  = 501
		kn.delmin_core  = 0.1
		kn.delmin_edge  = 0.01
		kn.minloc       = minloc

		if self.note_ef['HMODE'].get() == 1: fknot = kn.lmfit_fit(self.bdata[ttime]['jconst'][:,0],self.bdata[ttime]['jconst'][:,0],2,False)
		else:   fknot = kn.lmfit_fit(self.bdata[ttime]['jconst'][:,0],self.bdata[ttime]['jconst'][:,0],1,False)			

		return (pknot,fknot)

	def _make_current_knots(self,ttime):

		coren = int(float(self.note_ef['e7'].get()))
		edgen = int(float(self.note_ef['e8'].get()))
		knots = float(self.note_ef['e9'].get())
		knote = 1.

		xx = np.copy(self.bdata[ttime]['jconst'][:,0]);	yy = np.copy(self.bdata[ttime]['jconst'][:,1]);

		ind = np.where(xx>0.8);	    maxloc = xx[ind][np.argmax(yy[ind])];
		ind = np.where(xx<maxloc);	minloc = xx[ind][np.argmin(yy[ind])];

		if edgen < 6:
			print('Edge [#] > 5 is needed')
			edgen = 6;

		x = [maxloc, knote-0.006, knote-0.003, knote]

		pedn1 = (edgen-5)/2
		if not pedn1 == int(pedn1):
			pedn1 = int(pedn1-0.5);	pedn2 = pedn1 + 1;
		else:
			pedn1 = int(pedn1);	pedn2 = pedn1;

		dx = np.linspace(maxloc,0.994,pedn1+2)
		for i in range(1,pedn1+1):	x.append(dx[i])
		dx = np.linspace(min(minloc+0.03,maxloc-0.01),maxloc,pedn2+2)
		for i in range(0,pedn2+1):	x.append(dx[i])
		dx = np.linspace(minloc-0.05,min(minloc+0.03,maxloc-0.01),int(pedn2/2)+2)
		for i in range(0,int(pedn2/2)+1):	x.append(dx[i])

		if coren < 2:
			print('Edge [#] > 1 is needed')
			coren = 2;	

		x.append(knots)
		dx = np.linspace(knots,minloc-0.05,coren+1)
		for i in range(1,coren):	x.append(dx[i])

		x = np.sort(x)
		for i in range(len(x)):
			x[i] = round(x[i],3)
		return x

	def _make_efit_bnd_constraint(self,ttime):

		shotn = int(float(self.note_in['e1'].get()))
		efit_dir = os.getcwd() + '/EFITS/%s/g%06i.%06i'%(self.note_ef['runeq'].get(),shotn,ttime)

		eq = eqdsk.eqdsk(efit_dir)
		eq.read_eqdsk(efit_dir)
		rbdyc = get_midRZ(eq.rzbdy,1.9)
		return rbdyc

	def _run_efits(self):

		selection = self.note_ef['l1'].curselection()
		if len(selection) == 0: return
		line = self.note_ef['l1'].get(selection[0])
		time = int(float(line.split()[0]))
	
		currdir = os.getcwd()
		if not self.t_close1: return
		isgfile = self._check_gfiles(self.note_ef['runeq'].get())
		if not isgfile: return
	
		runeq = self.note_ef['runeq'].get().split('0')
		if runeq[0] == 'EFIT': runind = 1
		else: runind = int(float(runeq[1])) + 1

		shotn = int(float(self.note_in['e1'].get()))
		tlist = list(self.note_ef['l2'].get(0,'end'))
		self.note_ef['l2'].delete(0,'end')

		filename = 'EFITS/%s/k%06i.%06i'%(self.note_ef['runeq'].get(),shotn,time)
		self._read_kfile(filename)
		isrun = True
		try:    self._run_efit(runind,time)
		except: isrun = False
		os.chdir(currdir)
		if not os.path.isfile('EFITS/KIN%02i/g%06i.%06i'%(runind,shotn,time)): isrun = False
		if isrun: 
			if not line in tlist: tlist.append(line)

		for item in tlist: self.note_ef['l2'].insert('end',item)

		self._find_efit_run(load_data=True)

		return			

	def _run_efitm(self):
		
		currdir = os.getcwd()
		if not self.t_close1: return

		isgfile = self._check_gfiles(self.note_ef['runeq'].get())
		if not isgfile: return		

		runeq = self.note_ef['runeq'].get().split('0')
		if runeq[0] == 'EFIT': runind = 1
		else: runind = int(float(runeq[1])) + 1

		tlist = list(self.note_ef['l1'].get(0,'end'))
		shotn = int(float(self.note_in['e1'].get()))
		self.note_ef['l2'].delete(0,'end')
		for time in tlist:
			tt = int(float(time.split()[0]))
			filename = 'EFITS/%s/k%06i.%06i'%(self.note_ef['runeq'].get(),shotn,tt)

			self._read_kfile(filename)				
			isrun = True
			try:	self._run_efit(runind,tt)
			except: isrun = False
			os.chdir(currdir)
			if not os.path.isfile('EFITS/KIN%02i/g%06i.%06i'%(runind,shotn,tt)): isrun = False

			if isrun: self.note_ef['l2'].insert('end',time)

		os.chdir(currdir)
		self._find_efit_run(load_data=True)
		return

	def _run_efit(self,runind=0,ttime=0):

		riter = int(float(self.note_ef['e10'].get()))
		if riter > 100:
			self._run_efit_riter()
			return

		currdir = os.getcwd()
		mse_dat, pbeam, ebeam = self._make_efit_mse(ttime)
		if self.note_ef['HMODE'].get() == 1:
			pknot, fknot = self._make_pf_knots(ttime)
			jknot = self._make_current_knots(ttime)
		else:
			pknot = [0.0,0.4,1.0]
			fknot = [0.0,0.4,1.0]
			jknot = [0.93,0.95,0.98,1.0]

		rbdyc = self._make_efit_bnd_constraint(ttime)

		kfile_dir0= currdir +'/EFITS/%s/'%(self.note_ef['runeq'].get())
		efitdir   = currdir +'/EFITS/'
		save_dir  = currdir +'/EFITS/KIN%02i/'%(runind)
		try:	os.mkdir(save_dir)
		except:	pass

		shot = int(float(self.note_in['e1'].get()));	time = ttime;

		kfile_dir  = efitdir  + 'kfile_run'
		kfile_dir2 = save_dir + 'k%06i.%06i'%(shot,time)
		gfile_dir  = efitdir  + 'g%06i.%06i'%(shot,time)
		gfile_dir2 = save_dir + 'g%06i.%06i'%(shot,time)
		mfile_dir  = efitdir  + 'm%06i.%06i'%(shot,time)
		mfile_dir2 = save_dir + 'm%06i.%06i'%(shot,time)
		afile_dir  = efitdir  + 'a%06i.%06i'%(shot,time)
		afile_dir2 = save_dir + 'a%06i.%06i'%(shot,time)
		pfile_dir  = currdir  + '/PROFILES/MPROFILES_%s/%i/chease_kinprof_fit'%(self.note_ef['runeq'].get(),time)
		pfile_dir2 = save_dir + 'p%06i.%06i'%(shot,time)
		vfile_dir  = currdir  + '/PROFILES/MPROFILES_%s/%i/VT_fit.dat'%(self.note_ef['runeq'].get(),time)

		fit_dir    = efitdir  + 'fitout.dat'
		fit_dir2   = save_dir + 'fitout.dat_%06i.%06i'%(shot,time)
		w_dir2     = save_dir + 'wmhd.dat_%06i.%06i'%(shot,time)
		mse_dir2   = save_dir + 'mse.dat_%06i.%06i'%(shot,time)
		er_mse_dir2= save_dir + 'er_mse.dat_%06i.%06i'%(shot,time)
		pres_dir2  = save_dir + 'pres.dat_%06i.%06i'%(shot,time)
		j_dir2     = save_dir + 'jconst.dat_%06i.%06i'%(shot,time)
		map_dir2   = save_dir + 'map.dat_%06i.%06i'%(shot,time)
		hfile_dir2 = save_dir + 'history.dat_%06i.%06i'%(shot,time)

		run_history = dict()
		run_history['isermse'] = False
		run_history['isrit'] = False

		if self.note_ef['efitgrid'].get() == '257':
			efit_exec = efitdir +'efit257'
			efit_exec2 = efit_dir+'/efit257'
		elif self.note_ef['efitgrid'].get() == '129':
			efit_exec = efitdir +'efit129'
			efit_exec2 = efit_dir+'/efit129'
		else:
			efit_exec = efitdir +'efit65'
			efit_exec2 = efit_dir+'/efit65'

		if not os.path.isfile(efit_exec):	
			try:	copyfile(efit_exec2,efit_exec)
			except:	
				print('Error-')
				pass
		os.system('chmod 777 %s'%efit_exec)
		
		pres_mult = 1.; err=0.5;
		wmhd = self.mfit_opt['post'][time][3][0]
		if self.note_ef['PSCALE'].get() == 1: 
			pres_mult = wmhd/(self.bdata[time]['wth']+self.bdata[time]['wfast'])			
		else: err = -1

		if self.note_ef['ERMSE'].get() == 0:
			while (((self.note_ef['PSCALE'].get() == 1) and err > 0.05) or (err == -1)):
				print('>>> Pscale value %4.2f'%pres_mult)
				self._make_kfile(pknot=pknot,fknot=fknot,jknot=jknot,rbdc=rbdyc,ttime=ttime,pres_mult=pres_mult,mse_data=mse_dat)
				self._write_efit_input(kfile_dir)
				os.chdir(efitdir)
				os.system(efit_exec)
				efit_out = self._read_efit_result('fitout.dat',afile_dir,mfile_dir)
				err = abs(wmhd/efit_out['wmhd']-1.)
				if err > 0.05: pres_mult = wmhd/efit_out['wmhd'] * pres_mult

				if not os.path.isfile(gfile_dir):	
					print('>>> RUN FAIL..')
					os.chdir(currdir)
					return
				else:
					print('>>> RUN FINISHED..')
				if self.note_ef['PSCALE'].get() == 0: err = 0.


		elif(self.note_ef['ERMSE'].get() == 1):
			print('>>> Er correction (by Y.H.LEE) is started')
			if float(pbeam[0]) > 0.:   mseb = ebeam[0];
			elif float(pbeam[1]) > 0.: mseb = ebeam[1];
			else: mseb = ebeam[2];

			a1 = self.note_ef['e11'].get()
			a2 = self.note_ef['e12'].get()
			a3 = self.note_ef['e13'].get()
			self.note_ef['e11'].delete(0,'end')
			self.note_ef['e12'].delete(0,'end')
			self.note_ef['e13'].delete(0,'end')
			self.note_ef['e11'].insert(10,'201')
			self.note_ef['e12'].insert(10,'0.5')
			self.note_ef['e13'].insert(10,'1.e-4')

			self._make_kfile(pknot=pknot,fknot=fknot,jknot=jknot,rbdc=rbdyc,ttime=ttime,pres_mult=pres_mult,mse_data=mse_dat)
			self._write_efit_input(kfile_dir)
			os.chdir(efitdir)
			os.system(efit_exec)
			if not os.path.isfile(gfile_dir):
				print('>>> RUN FAIL..')
				return

			run_history['isermse'] = True
			efit_out = self._read_efit_result('fitout.dat',afile_dir,mfile_dir)
			run_history[0] = copy.deepcopy([mse_dat,efit_out])
			print('>>> MSE 0/1 iteration')

			mse_dat['tgamma'] = Er_mse_corr(gfile=gfile_dir,pfile=pfile_dir,vfile=vfile_dir,tgamma=mse_dat['tgamma0'],rgamma=mse_dat['rrrgam'],ebeam=mseb)
			
			self.note_ef['e11'].delete(0,'end')
			self.note_ef['e12'].delete(0,'end')
			self.note_ef['e13'].delete(0,'end')
			self.note_ef['e11'].insert(10,a1)
			self.note_ef['e12'].insert(10,a2)
			self.note_ef['e13'].insert(10,a3)

			self._make_kfile(pknot=pknot,fknot=fknot,jknot=jknot,rbdc=rbdyc,ttime=ttime,pres_mult=pres_mult,mse_data=mse_dat)
			self._write_efit_input(kfile_dir)
			os.chdir(efitdir)
			os.system(efit_exec)
			if not os.path.isfile(gfile_dir):
				print('>>> RUN FAIL..')
				os.chdir(currdir)
				return
			else:
				print('>>> MSE 1/1 iteration')

			efit_out = self._read_efit_result('fitout.dat',afile_dir,mfile_dir)
			run_history[1] = copy.deepcopy([mse_dat,efit_out])

			print('>>> RUN FINISHED..')	

		if (self.note_ef['PSCALE2'].get() == 1):
			delp = -0.03
			chi_old = efit_out['chi']; pres_mult = pres_mult + delp
			print('>>> Pscale value %f, XI = %f'%(pres_mult-delp,chi_old))
			self._make_kfile(pknot=pknot,fknot=fknot,jknot=jknot,rbdc=rbdyc,ttime=ttime,pres_mult=pres_mult,mse_data=mse_dat)
			self._write_efit_input(kfile_dir)
			os.chdir(efitdir)
			os.system(efit_exec)
			efit_out = self._read_efit_result('fitout.dat',afile_dir,mfile_dir)
			chi = efit_out['chi']
			print('>>> Pscale value %f, XI = %f'%(pres_mult,chi))
			if chi > chi_old: delp =  - delp
			pres_mult_old = pres_mult
			pres_mult = pres_mult + delp
			scale_p = True
			chi_old = chi

			while (scale_p):
				self._make_kfile(pknot=pknot,fknot=fknot,jknot=jknot,rbdc=rbdyc,ttime=ttime,pres_mult=pres_mult,mse_data=mse_dat)
				self._write_efit_input(kfile_dir)
				os.chdir(efitdir)
				os.system(efit_exec)
				efit_out = self._read_efit_result('fitout.dat',afile_dir,mfile_dir)
				chi = efit_out['chi']

				print('>>> Pscale value %f, XI = %f'%(pres_mult,chi))

				if chi_old > chi: 
					pres_mult = 2.*pres_mult - pres_mult_old
					print('>>> XI[O/N] %s %s'%(chi_old,chi))
				else: 
					pres_mult = pres_mult_old
					self._make_kfile(pknot=pknot,fknot=fknot,jknot=jknot,rbdc=rbdyc,ttime=ttime,pres_mult=pres_mult,mse_data=mse_dat)
					self._write_efit_input(kfile_dir)
					os.chdir(efitdir)
					os.system(efit_exec)
					efit_out = self._read_efit_result('fitout.dat',afile_dir,mfile_dir)
					print('>>> Scale Ended')
					scale_p = False

				pres_mult_old = (pres_mult+pres_mult_old)/2.
				chi_old = chi
				if not os.path.isfile(gfile_dir):	
					print('>>> RUN FAIL..')
					os.chdir(currdir)
					return
				else:
					print('>>> RUN FINISHED..')

		efit_out['pscale'] = pres_mult
		run_history[0] = copy.deepcopy([mse_dat,efit_out])

		move(kfile_dir,    kfile_dir2)
		move(gfile_dir,    gfile_dir2)
		move(mfile_dir,    mfile_dir2)
		move(fit_dir,        fit_dir2)
		move(afile_dir,    afile_dir2)
		copyfile(pfile_dir,pfile_dir2)
	
		f = open(hfile_dir2,'wb')
		pickle.dump(run_history,f)
		f.close()

		os.chdir(currdir)
		return

	def _run_efit_riter(self):

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
				print('Error-')
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

		psin_temp =copy.deepcopy(self.psin)
		pres_temp =copy.deepcopy(self.pres)
		pp_temp   =copy.deepcopy(self.pp)
		ffp_temp  =copy.deepcopy(self.ffp)
		jconst_temp =copy.deepcopy(self.jconst)
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


				make_kfile(self,True)
				write_efit_input(self,kfile_dir)
				os.chdir(efitdir)
				os.system(efit_exec)
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
				self.ts_run = False

				make_kfile(self,True)
				write_efit_input(self,kfile_dir)
				os.chdir(efitdir)
				os.system(efit_exec)
				if not os.path.isfile(gfile_dir):
					print('>>> RUN FAIL..')
					return
				self.ts_run = True
				run_history = dict()
				nn = int(float(self.MenuVar15.get()))
				run_history[0] = np.copy(self.tgamma)
				run_history['q'] = np.zeros(nn+1)
				run_history['xi']= np.zeros(nn+1)
				read_efit_result(self,'fitout.dat')
				run_history['q'][0] = self.efit_q[0]
				run_history['xi'][0]= self.efit_chi
				print('>>> ITER #%i/%i > MSE #%i/%i iteration'%(riteri+1,riter+1,0.,nn))
				for i in range(nn):
					self.tgamma3 = Er_mse_corr(gfile=gfile_dir,pfile='../../'+self.e4.get(),vfile='../../'+self.e45.get(),tgamma=self.tgamma,rgamma=self.rrrgam,ebeam=float(self.__dict__['StrVar%d'%ind].get()))
					if ((i== nn-1) and (riteri == (riter))):
						self.StrVar38.set(a1)
						self.StrVar39.set(a2)
						self.StrVar40.set(a3)
					make_kfile(self,True)
					write_efit_input(self,kfile_dir)
					run_history[i+1] = np.copy(self.tgamma3)
					os.chdir(efitdir)
					os.system(efit_exec)
					if not os.path.isfile(gfile_dir):
						print('>>> RUN FAIL..')
						os.chdir(self.currdir)
						return
					else:
						print('>>> ITER #%i/%i > MSE #%i/%i iteration'%(riteri+1,riter+1,i+1,nn))
	
					read_efit_result(self,'fitout.dat')
					run_history['q'][i+1] = self.efit_q[0]
					run_history['xi'][i+1]= self.efit_chi
				
				f = open(temp_dir+'/mse_er.dat_%i'%riteri,'w')
				line = '%9s'%'R[m]'
				for i in range(nn+1):line = line +'\t%9s'%('tgam#%1i'%i)
				print(line);    f.write(line+'\n');
				for i in range(len(self.rrrgam)):
					line = '%9.6f'%self.rrrgam[i]
					for j in range(nn+1):line = line +'\t%9.6f'%(run_history[j][i])
					print(line);    f.write(line+'\n');
			
				line = '%9s'%'xisq'
				for i in range(nn+1):line = line +'\t%9.6f'%run_history['xi'][i]
				print(line);    f.write(line+'\n');
				line = '%9s'%'q0'
				for i in range(nn+1):line = line +'\t%9.6f'%run_history['q'][i]
				print(line);    f.write(line+'\n');
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

			self.psin =copy.deepcopy(psin_temp)
			self.pres =copy.deepcopy(pres_temp)
			self.jconst =copy.deepcopy(jconst_temp)
			self.pp   =copy.deepcopy(pp_temp)
			self.ffp  =copy.deepcopy(ffp_temp)

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


		for string in self.efit_index:
			menu.add_command(label=string, command=lambda value=string: option_menu_update(self.MenuVar10,value,self.button_func9b,True))
	
		self.efit_list3 = []
		os.chdir(self.currdir)

		return

	def _read_kfile(self,filename):

		self.kfile['in1'] = []
		self.kfile['in2'] = []
		self.kfile['in3'] = []
		self.kfile['endline'] = ''
		
		count = 1
		
		f4 = open(filename,'r')
		
		while True:
			line = f4.readline()
			if not line: break
			
			if (line.find('&IN1') > -1):
				count = 1
			if (line.find('&INWANT') > -1):
				count = 2
			if (line.find('&INS') > -1):
				count = 3
			if not (line.find('/') > -1):

				if (line.find('MAG') == -1 and line.find('MSE') == -1):
					if (count == 1):
						self.kfile['in1'].append(line)
					elif (count == 2 and line.find('&INWANT') == -1):
						self.kfile['in2'].append(line)
					elif (count == 3 ):
						self.kfile['in3'].append(line)
				else:
					self.kfile['endline'] = line.split('M')[0]
		f4.close()		
		for k in range(1,4): self.kfile['in%02i'%k] =copy.deepcopy(self.kfile['in%i'%k])
		return

	def _load_kfile(self):

		self.kfile['fwtsi']  = get_kfile_data(self.kfile['in1'],'FWTSI')
		self.kfile['fwtmp2'] = get_kfile_data(self.kfile['in1'],'FWTMP2')

		self.note_ef['FLOOPN'] = len(self.kfile['fwtsi'])
		self.note_ef['MCOILN'] = len(self.kfile['fwtmp2'])

		print('>>> Probe signals ->','Magnetic [#]',self.note_ef['MCOILN'],'Flux-loop [#]',self.note_ef['FLOOPN'])

		for k in range(self.note_ef['MCOILN']): self.note_ef['MCOIL%03i'%(k+1)].set(int(self.kfile['fwtmp2'][k]))
		for k in range(self.note_ef['FLOOPN']): self.note_ef['FLOOP%03i'%(k+1)].set(int(self.kfile['fwtsi'][k]))
		self.note_ef['KLOAD'] = True
		return

	def _get_kfile(self):

		selection = self.note_ef['l1'].curselection()
		if len(selection) == 0: return
		line = self.note_ef['l1'].get(selection[0])
		time = int(float(line.split()[0]))

		filen = 'EFITS/%s/k%06i.%06i'%(self.note_ef['runeq'].get(),self.note_in['shot'],time)

		self._read_kfile(filen)
		self._load_kfile()
		self.note_ef['KLOAD']  = True

		return

	def _make_kfile(self,pknot='',fknot='',jknot='',rbdc=2.3,ttime=0,pres_mult=1.,mse_data=[]):

		prest = self.bdata[ttime]['prest']
		jconst= self.bdata[ttime]['jconst']
		rpress = np.zeros(51)
		if (self.note_ef['HMODE'].get() == 1):
			rpress0 = -1.0*np.linspace(0,0.9,25)
			rpress1 = -1.0*np.linspace(0.9,1.0,51-25)
			rpress[0:25]  = rpress0
			rpress[25:51] = rpress1

		else:
			rpress = -1.0*np.linspace(0,1.,51)
		
		self.kfile['rpress'] = np.copy(rpress)
		for k in range(1,4): self.kfile['in%i'%k] = copy.deepcopy(self.kfile['in%02i'%k])

		self.kfile['in1'].append('\n  mxiter = -%s\n'%self.note_ef['e11'].get())
		self.kfile['in1'].append('  relax = %s \n'%self.note_ef['e12'].get())
		self.kfile['in1'].append('  pcurbd = 0.0 \n')
		self.kfile['in1'].append('  fcurbd = 0.0 \n\n')

		if self.note_ef['JCONST'].get() == 1: jknot = np.array(jknot,dtype='double')
		else:	jknot = np.array([0,1.])
		self.kfile['in1'].append('  kppfnc  = 6, kppknt = %i, pptens = 1\n'%len(pknot))

		line = '  ppknt   = '
		for i in range(len(pknot)): line = line + '%4.3f '%float(pknot[i])
		self.kfile['in1'].append(line+'\n\n')

		self.kfile['in1'].append('  kfffnc  = 6, kffknt = %i, fftens = 1\n'%len(fknot))
		line = '  ffknt   = '
		for i in range(len(fknot)): line = line + '%4.3f '%float(fknot[i])
		self.kfile['in1'].append(line+'\n\n')
		
		self.kfile['in1'].append('  error   =  %s, \n'%self.note_ef['e13'].get())
		self.kfile['in1'].append('  errmin  =  %s, \n'%self.note_ef['e13'].get())
		self.kfile['in1'].append(' \n')
		self.kfile['in1'].append('  fwtcur  = %s \n'%self.note_ef['e14'].get())
		self.kfile['in1'].append('  fwtqa   = %i \n'%self.note_ef['QFIX'].get())
		self.kfile['in1'].append('  qvfit   = %s \n\n'%self.note_ef['e16'].get())		

		if (self.note_ef['JCONST'].get() == 1):
			self.kfile['in1'].append('  RZEROJ  = %i*0 \n'%len(jknot))
			self.kfile['in1'].append('  KZEROJ  = %i \n\n'%len(jknot))

		self.kfile['in1'].append('  NBDRY = 1 \n')
		self.kfile['in1'].append('  RBDRY = %s\n'%rbdc)
		self.kfile['in1'].append('  ZBDRY = 0.\n')

		if self.note_ef['KLOAD']: 
			self.kfile['in1'].append('  fwtsi = \n')
			line = '  '
			count = 0
			for i in range(self.note_ef['FLOOPN']):
				count = count + 1
				line = line + '%i '%self.note_ef['FLOOP%03i'%(i+1)].get()
				if (count == 10):
					line = line + '\n  '
					count = 0;
			self.kfile['in1'].append(line)
			self.kfile['in1'].append('\n\n\n')

			self.kfile['in1'].append('  fwtmp2 = \n')
			line = '  '
			count = 0
			for i in range(self.note_ef['MCOILN']):
				count = count + 1
				line = line + '%i '%self.note_ef['MCOIL%03i'%(i+1)].get()
				if (count == 20):
					line = line + '\n  '
					count = 0;
			self.kfile['in1'].append(line)
			self.kfile['in1'].append('\n\n\n')
		
		line = make_kfile_str(self.kfile['rpress'],'RPRESS')
		self.kfile['in1'].append(line)
		self.kfile['in1'].append('\n')
		
		presf = interp1d(prest[:,0],prest[:,1],'cubic')
		presff= interp1d(prest[:,0],prest[:,2],'cubic')
		self.kfile['prest']      = presf(-self.kfile['rpress'])
		self.kfile['prestf']     = presff(-self.kfile['rpress'])
		self.kfile['sigpres']    = np.copy(self.kfile['prest']*0.1)
		self.kfile['sigpres'][0] = self.kfile['prest'][0]*0.03

		core_sup = float(self.note_ef['e15'].get())
		alpha = 0.
		w1 = pres_mult*alpha + (1-alpha)
		w2 = (pres_mult-1)*(1-alpha)

		w1 = np.copy(self.kfile['prest']); w2 = np.copy(w1)
		for k in range(len(w1)):
			w1[k] = 1 - (1-pres_mult) * max(np.tanh((0.85+self.kfile['rpress'][k])/0.2),0.)
			w2[k] = w1[k] #pres_mult

		if float(core_sup) > 0.:
			self.kfile['prest'][0] = core_sup * (w1[1]*self.kfile['prest'][1] + (w2[1]-w1[1])*self.kfile['prestf'][1])
		
		line = make_kfile_str(w1*self.kfile['prest'] + (w2-w1)*self.kfile['prestf'],'PRESSR')
		self.kfile['in1'].append(line)
		self.kfile['in1'].append('\n')

		line = make_kfile_str(self.kfile['sigpres'],'SIGPRE')
		self.kfile['in1'].append(line)
		self.kfile['in1'].append('\n\n')
		
		self.kfile['in1'].append(' NPRESS  =        %3i,\n'%len(self.kfile['prest']))
		self.kfile['in1'].append(' NBEAM   =        %3i,\n'%len(self.kfile['prest']))
		self.kfile['in1'].append('\n')

		#Dummy fast ion info (can be updated later)
		line = make_kfile_str(-self.kfile['rpress'],'SIBEAM')
		self.kfile['in1'].append(line)
		self.kfile['in1'].append('\n')

		line = make_kfile_str(self.kfile['prest']*0.3,'PBEAM')
		self.kfile['in1'].append(line)
		self.kfile['in1'].append('\n')

		line = make_kfile_str(self.kfile['prest']/self.kfile['prest'][0]*4.e18,'DNBEAM')
		self.kfile['in1'].append(line)
		self.kfile['in1'].append('\n')
		
		line = make_kfile_str(self.kfile['prest']/self.kfile['prest'][0]*4.e18*1.627*1.e-27,'DMASS')
		self.kfile['in1'].append(line)
		self.kfile['in1'].append('\n\n')
		
		self.kfile['in1'].append(' NMASS  =        %3i,\n'%len(self.kfile['prest']))
		self.kfile['in1'].append(' KPRFIT   =      1,\n')
		self.kfile['in1'].append('\n')
		
		jconstf = interp1d(jconst[:,0],jconst[:,1],'cubic')
		if self.note_ef['JCONST'].get()==1:
			jc = jconstf(jknot)
			line = make_kfile_str(jknot,'SIZEROJ')
			self.kfile['jc'] = line + '\n'
			line = make_kfile_str(jc,'VZEROJ')
			self.kfile['jc'] = self.kfile['jc'] + line + '\n'

		self.kfile['in3t'] = []
		self.kfile['in3t'].append('  &INS\n')
		for key in ['tgamma','sgamma','fwtgam','rrrgam','zzzgam','aa1gam','aa2gam','aa3gam','aa4gam','aa5gam','aa6gam']:
			line = make_kfile_str(mse_data[key],str(key).upper())
			self.kfile['in3t'].append(line)
			self.kfile['in3t'].append('\n')

		self.kfile['in3t'].append(line)
		self.kfile['in3t'].append('\n\n')
		self.kfile['in3t'].append(' MSEBKP     =       0,\n')
		self.kfile['in3t'].append(' MSEFITFUN  =       1 \n')
		self.kfile['in3t'].append('\n')
			
		return

	def _custom_field_config(self):

		if not self.note_ef['KLOAD']: 
			print('>>> Load K-file first...')
			return

		if not self.t_close1:
			print('>>> Coil custom Window is already opened...')
			return			

		self.t1 = tk.Toplevel(self.root)
		self.t1.wm_title("Coil Config")
		self.t_close1 = False		
		self.t1.protocol('WM_DELETE_WINDOW', lambda: self._detect_close())

		self.t1.resizable(0,0)

		self.l1 = tk.Label(self.t1, text="======= M-COIL =======",justify='center')
		self.l1.grid(row=0,column=0,columnspan=10,pady=5)
		self.l1 = tk.Label(self.t1, text="======= F-LOOP =======",justify='center')
		self.l1.grid(row=0,column=10,columnspan=10,pady=5)

		count = 0
		for i in range(20):
			for j in range(10):
				count = count + 1
				if (count > self.note_ef['MCOILN']): break
				self.note_ef['cm_%03i'%count] = tk.Checkbutton(self.t1,variable=self.note_ef['MCOIL%03i'%count],bg='lime')
				self.note_ef['cm_%03i'%count].grid(row=i+1,column=j)	

		count = 0
		for i in range(20):
			for j in range(10):
				count = count + 1
				if (count > self.note_ef['FLOOPN']): break
				self.note_ef['cf_%03i'%count] = tk.Checkbutton(self.t1,variable=self.note_ef['FLOOP%03i'%count],bg='magenta')
				self.note_ef['cf_%03i'%count].grid(row=i+1,column=j+10)	
		return

	def _detect_close(self):

		try:	self.t1.destroy()
		except:	print('>>> Error-1');pass
		self.t_close1 = True
		return

	def _write_efit_input(self,filename,type=1):

		f = open(filename,'w')
		
		if type == 1:	
			for item in self.kfile['in1']:
				f.write(item)
		else:		
			for item in self.kfile['in01']:
				f.write(item)
		f.write('/\n')
		f.write('&INWANT\n')
		if type == 1:	
			for item in self.kfile['in2']:
				f.write(item)
		else:		
			for item in self.kfile['in02']:
				f.write(item)
		if type == 1:
			if (self.note_ef['JCONST'].get()==1):
				f.write(self.kfile['jc'])
			f.write('/\n')

		if type == 1:	
			for item in self.kfile['in3t']:
				f.write(item)
		else:
			for item in self.kfile['in03']:
				f.write(item)
			
		f.write('/\n')
		if type == 1:	f.write(self.kfile['endline']+' KIN')
		else:	f.write(self.kfile['endline']+' MSE')
		f.close()
		
		return

	def _read_efit_result(self,fitout,afile,mfile):

		efit_out = dict()

		f = open(afile,'r')
		for k in range(5): line = f.readline().split()
		efit_out['chi'] = float(line[0])
		for k in range(12): line = f.readline().split()
		efit_out['wmhd'] = float(line[-1])/1.e3
		f.close()
		f=Dataset(mfile,'r',format='NETCDF4')
		efit_out['tangam'] = np.array(f['tangam'][0])
		efit_out['rrgam']  = np.array(f['rrgam'][0])
		efit_out['cmgam']  = np.array(f['cmgam'][0])
		f.close()

		f = open(fitout,'r')

		efit_out['chi_ploop'] = np.zeros(self.note_ef['FLOOPN'])
		efit_out['chi_mloop'] = np.zeros(self.note_ef['MCOILN']) 

		efit_out['m_ploop'] = np.zeros(self.note_ef['FLOOPN'])
		efit_out['m_mloop'] = np.zeros(self.note_ef['MCOILN'])
		efit_out['c_ploop'] = np.zeros(self.note_ef['FLOOPN'])
		efit_out['c_mloop'] = np.zeros(self.note_ef['MCOILN'])

		clen1 = int(np.ceil(self.note_ef['MCOILN']/8))
		clen2 = int(np.ceil(self.note_ef['FLOOPN']/8))
		clen3 = int(np.ceil(self.note_in['nch']/8))

		efit_out['num'] = 1
		
		efit_out['conve'] = 0.
		efit_out['bp'] = 0.

		count = 0
		while True:
			line = f.readline()
			if not line:	break
			if line.find('chi psi loops') > -1:
				for i in range(self.note_ef['FLOOPN']):
					if (int(i/8)==(i/8)):	
						line = f.readline().split()
						j = 0
					efit_out['chi_ploop'][i] = float(line[j])
					j = j + 1;

			elif line.find('inner magnetic probes') > -1:
				for i in range(self.note_ef['MCOILN']):
					if (int(i/8)==(i/8)):	
						line = f.readline().split()
						j = 0
					efit_out['chi_mloop'][i] = float(line[j])
					j = j + 1;
			elif line.find('betap') > -1:
				efit_out['bp'] = float(line.split('=')[2].split('li')[0])

			elif line.find('EFITD') > -1:
				efit_out['num']  = int(line.split('EFITD')[1].split('dx2')[0])
				efit_out['psin'] = np.zeros(efit_out['num'])
				efit_out['pp']   = np.zeros(efit_out['num'])
				efit_out['ffp']  = np.zeros(efit_out['num'])
				efit_out['q']    = np.zeros(efit_out['num'])
				efit_out['jav']  = np.zeros(efit_out['num'])
				for k in range(8): line = f.readline()
				line = f.readline().split()
				line = f.readline().split()
				efit_out['betat'] = float(line[2])
				efit_out['betap'] = float(line[5])
				efit_out['bp']    = float(line[5])
				efit_out['li']    = float(line[8])


			elif line.find('plasma summary') > -1:
				count = 1
				line = f.readline()
				for i in range(efit_out['num']):
					line = f.readline().split()
					efit_out['psin'][i] = float(line[1])
					efit_out['pp'][i]   = float(line[3])
					efit_out['ffp'][i]  = float(line[5])
					efit_out['q'][i]    = float(line[8])

			elif (line.find('pflux') > -1 and count == 1):
				for i in range(efit_out['num']):
					line = f.readline().split()
					efit_out['jav'][i] = float(line[2])

			elif line.find('calculated psi-loops') > -1:
				for i in range(self.note_ef['FLOOPN']):
					if (int(i/4)==(i/4)):	
						line = f.readline().split()
						j = 0
					efit_out['c_ploop'][i] = float(line[j])
					j = j + 1;
			elif line.find('measured psi-loops') > -1:
				for i in range(self.note_ef['FLOOPN']):
					if (int(i/4)==(i/4)):	
						line = f.readline().split()
						j = 0
					efit_out['m_ploop'][i] = float(line[j])
					j = j + 1;

			elif line.find('calculated total plasma current') > -1:
				efit_out['ip'] = float(f.readline())

			elif line.find('calculated magnetic probes') > -1:
				for i in range(self.note_ef['MCOILN']):
					if (int(i/4)==(i/4)):	
						line = f.readline().split()
						j = 0
					efit_out['c_mloop'][i] = float(line[j])
					j = j + 1;
			elif line.find('measured magnetic probes') > -1:
				for i in range(self.note_ef['MCOILN']):
					if (int(i/4)==(i/4)):	
						line = f.readline().split()
						j = 0
					efit_out['m_mloop'][i] = float(line[j])
					j = j + 1;				

			elif line.find('iteration') > -1:
				line = f.readline().split()
				for i in range(999):
					try:	line = f.readline().split()
					except:	line = []
					if len(line) == 0.: break
					efit_out['conve'] = float(line[1])
		
		f.close()

		return efit_out

	def _find_efit_run(self,load_data = False):

		if self.note_ef['runeq'].get() == 'EFIT01': dirs = 'EFITS/KIN01/'
		else: dirs = 'EFITS/KIN%02i/'%(int(float(self.note_ef['runeq'].get().split('N')[1]))+1)
	
		tlist = list(self.note_ef['l1'].get(0,'end'))
		self.note_ef['l2'].delete(0,'end')
		for time in tlist:
			tt = int(float(time.split()[0]))
			filen = dirs + 'g%06i.%06i'%(self.note_in['shot'],tt)
			if os.path.isfile(filen): self.note_ef['l2'].insert('end',time)

		if load_data: self._load_efit_runs()

		return

	def _update_efit_list(self):

		self.note_pr['gflag'] = ['EFIT01']
		self.note_nb['gflag'] = ['EFIT01']
		self.note_ef['gflag'] = ['EFIT01']

		for k in range(1,100):
			if os.path.isdir('EFITS/KIN%02i'%k):
				self.note_pr['gflag'].append('KIN%02i'%k)
				self.note_nb['gflag'].append('KIN%02i'%k)
			if os.path.isdir('EQUIL/KIN%02i'%k):
				self.note_ef['gflag'].append('KIN%02i'%k)

		return

	def _load_efit_runs(self):

		self.efdata = dict()
		tlist = list(self.note_ef['l2'].get(0,'end'))
		if self.note_ef['runeq'].get() == 'EFIT01': dirs = 'EFITS/KIN01'
		else:
			dirs = 'EFITS/KIN%02i'%(int(float(self.note_ef['runeq'].get().split('N')[-1]))+1)

		for time in tlist:
			tt = int(float(time.split()[0]))
			filen = '%s/history.dat_%06i.%06i'%(dirs,self.note_in['shot'],tt)
			mfile = '%s/m%06i.%06i'%(dirs,self.note_in['shot'],tt)
			if os.path.isfile(filen):
				f = open(filen,'rb')
				self.efdata[tt] = dict()
				self.efdata[tt]['f'] = copy.deepcopy(pickle.load(f))
				f.close()
				self.efdata[tt]['m'] = Dataset(mfile,'r')
		return

	def _draw_efit(self):

		selection = self.note_ef['l2'].curselection()
		if len(selection) == 0: return
		line = self.note_ef['l2'].get(selection[0])
		time = int(line.split()[0])

		tlist = list(self.note_ef['l2'].get(0,'end'))
		tmin = float(tlist[0].split()[0])
		tmax = float(tlist[-1].split()[0])

		self.figure['name']['home'].canvas.draw_idle()
		ax = self.figure['name']['home'].axes[0]
		plottype = self.note_ef['plottype'].get()
		
		if   plottype == 'PRESS' : plotspace = 'r'
		elif plottype == 'JTOR'  : plotspace = 'r'
		elif plottype == 'QPROF' : plotspace = 'r'
		elif plottype == 'MSE'   : plotspace = 'rr'
		elif plottype == 'MCOIL' : plotspace = 'i'
		elif plottype == 'FLOOP' : plotspace = 'i'
		elif plottype == 'PSCALE': plotspace = 't'
		elif plottype == 'CHISQR': plotspace = 't'
		elif plottype == 'WMHD':   plotspace = 't'
		elif plottype == 'BPOL':   plotspace = 't'
		elif plottype == 'LI':     plotspace = 't'

		self.note_ef['plotspace'] = plotspace
		ymax = float(self.note_ef['e1'].get())
		if not plottype == self.note_ef['plottype_old']:
			self._clf_plot(5)
			self.figure['pegend'][0] = []
			self.figure['legend'][0] = []
			self.note_ef['plottype_old'] = plottype
			ymax = 0.
			self.note_ef['e1'].delete(0,'end')
			self.note_ef['e1'].insert(10,'0.')

		if    plotspace =='r': ax.set_xlabel('$\\psi_N$')
		elif  plotspace =='rr':ax.set_xlabel('$R [m]$')
		elif  plotspace =='i': ax.set_xlabel('Ch.[#]')
		else: ax.set_xlabel('times [ms]')

		count = len(self.figure['legend'][0]) + 1

		times = []
		for key in self.efdata.keys(): times.append(key)
		times = np.array(times,dtype='float')
		times = np.sort(times)

		if plotspace == 'r': self.note_ef['xmin'] = 0.; self.note_ef['xmax'] = 1.;
		elif plotspace == 'rr': self.note_ef['xmin'] = 1.7; self.note_ef['xmax'] = 2.3;
		elif plotspace == 'i': self.note_ef['xmin'] = 0.; self.note_ef['xmax'] = None;
		else: self.note_ef['xmin'] = 0.7*min(times); self.note_ef['xmax'] = 1.3*max(times);

		ax.set_xlim(self.note_ef['xmin'],self.note_ef['xmax'])
		
		vals = np.copy(times)
		fdata1 = self.efdata[time]['f'][0][1]
		fdata2 = self.efdata[time]['f'][0][0]
		mdata = self.efdata[time]['m']

		ymin = 0.

		if self.note_ef['runeq'].get() == 'EFIT01': efitname = 'KIN01'
		else: efitname = 'KIN%02i'%(int(float(self.note_ef['runeq'].get().split('N')[1]))+1)

		if   plottype == 'PRESS' : 
			line, = ax.plot(-mdata['rpress'][0],mdata['cpress'][0]/1.e3,c='C%i'%count)
			line1 = ax.errorbar(-mdata['rpress'][0],mdata['pressr'][0]/1.e3,yerr=mdata['sigpre'][0]/1.e3,fmt='x',c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $P_{BEAM}$-%s'%(time,efitname))
			self.note_ef['ymax'] = 1.1 * max(mdata['cpress'][0]) / 1.e3
			line2 = '[kPa]'
		elif plottype == 'JTOR' :
			jsign = np.sign(fdata1['jav'][0])
			line, = ax.plot(fdata1['psin'],jsign*fdata1['jav']/1.e6,c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $J_{IND}$-%s'%(time,efitname))
			self.note_ef['ymax'] = 1.1 * max(jsign*fdata1['jav']) / 1.e6
			line2 = '[MA/$m^2$]'
		elif plottype == 'QPROF' : 
			line, = ax.plot(fdata1['psin'],fdata1['q'],c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $q$-%s'%(time,efitname))
			self.note_ef['ymax'] = 1.1 * max(fdata1['q'])
			line2 = '[a.u]'
		elif plottype == 'MSE' : 
			line, = ax.plot(fdata1['rrgam'],fdata1['cmgam'],c='C%i'%count)
			line2 = ax.errorbar(fdata1['rrgam'],fdata1['tangam'],yerr=mdata['siggam'][0],fmt='x',c='C%i'%count)
			line3 = ax.axhline(y=0.,c='gold',linestyle='--')
			self.figure['legend'][0].append('%i ms $TGAM$-%s'%(time,efitname))
			self.note_ef['ymax'] = 1.1 * max(fdata1['cmgam']) ; ymin = 1.1 * min(fdata1['cmgam']) 
			line2 = '[a.u]'
		elif plottype == 'MCOIL':
			len1  = len(mdata['expmpi'][0])
			line, = ax.plot(np.linspace(1,len1,len1),mdata['cmpr2'][0]+2.,'--',c='C%i'%count)	
			line2 = ax.scatter(np.linspace(1,len1,len1),mdata['expmpi'][0]+2.,marker='x',c='C%i'%count)	
			line2 = ax.plot(np.linspace(1,len1,len1),mdata['saimpi'][0],'--d',c='C%i'%count)	
			self.figure['legend'][0].append('%i ms $MCOIL$-%s'%(time,self.note_ef['runeq'].get()))		
			self.note_ef['ymax'] = 1.1 * max(mdata['cmpr2'][0]+2.)
			ax.set_xlim(0,len1+1)
			line2 = '[a.u]'
		elif plottype == 'FLOOP' : 
			len1  = len(mdata['silopt'][0])
			line, = ax.plot(np.linspace(1,len1,len1),mdata['csilop'][0]+2.,'--',c='C%i'%count)	
			line2 = ax.scatter(np.linspace(1,len1,len1),mdata['silopt'][0]+2.,marker='x',c='C%i'%count)	
			line2 = ax.plot(np.linspace(1,len1,len1),mdata['saisil'][0],'--d',c='C%i'%count)	
			self.figure['legend'][0].append('%i ms $FLOOP$-%s'%(time,self.note_ef['runeq'].get()))		
			self.note_ef['ymax'] = 1.1 * max(mdata['silopt'][0]+2.)
			ax.set_xlim(0,len1+1)
			line2 = '[a.u]'
		elif plottype == 'PSCALE' : 
			for k in range(len(times)): vals[k] = self.efdata[int(times[k])]['f'][0][1]['pscale']
			line, = ax.plot(times,vals,'-x',c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $P_{SCALE}$-%s'%(time,self.note_ef['runeq'].get()))		
			self.note_ef['ymax'] = 1.2 * max(vals)
			line2 = '[a.u]'
		elif plottype == 'CHISQR': 
			for k in range(len(times)): vals[k] = self.efdata[int(times[k])]['f'][0][1]['chi']
			line, = ax.plot(times,vals,'-x',c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $\\chi$-%s'%(time,self.note_ef['runeq'].get()))		
			self.note_ef['ymax'] = 1.2 * max(vals)
			line2 = '[a.u]'			
		elif plottype == 'WMHD': 
			for k in range(len(times)): vals[k] = self.efdata[int(times[k])]['f'][0][1]['wmhd']
			line, = ax.plot(times,vals,'-x',c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $W_{MHD}$-%s'%(time,self.note_ef['runeq'].get()))		
			self.note_ef['ymax'] = 1.2 * max(vals)
			line2 = '[kJ]'
		elif plottype == 'BPOL': 
			for k in range(len(times)): vals[k] = self.efdata[int(times[k])]['f'][0][1]['bp']
			line, = ax.plot(times,vals,'-x',c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $\\beta_{P}$-%s'%(time,self.note_ef['runeq'].get()))		
			self.note_ef['ymax'] = 1.2 * max(vals)
			line2 = '[a.u]'
		elif plottype == 'LI': 
			for k in range(len(times)): vals[k] = self.efdata[int(times[k])]['f'][0][1]['li']
			line, = ax.plot(times,vals,'-x',c='C%i'%count)		
			self.figure['legend'][0].append('%i ms $l_{i}$-%s'%(time,self.note_ef['runeq'].get()))		
			self.note_ef['ymax'] = 1.2 * max(vals)
			line2 = '[a.u]'			


		self.figure['pegend'][0].append(line)

		
		if ymax < self.note_ef['ymax']:
			self.note_ef['e1'].delete(0,'end')
			self.note_ef['e1'].insert(10,'%4.2f'%self.note_ef['ymax'])
			ax.set_ylim(ymin,self.note_ef['ymax'])

		self.note_ef['e2'].delete(0,'end')
		self.note_ef['e2'].insert(10,'%4.2f'%ymin)

		ax.set_ylabel(line2)
		ax.legend(self.figure['pegend'][0],self.figure['legend'][0])
		self.figure['name']['home'].tight_layout()			
		return

	def _update_option_menus(self):

		#GFILE
		menu1 = self.note_pr['m3']["menu"]
		menu1.delete(0,'end')
		menu2 = self.note_nb['m2']["menu"]
		menu2.delete(0,'end')		
		menu3 = self.note_ef['m2']["menu"]
		menu3.delete(0,'end')				

		selec1 = self.note_pr['runeq'].get()
		selec2 = self.note_nb['runeq'].get()
		selec3 = self.note_ef['runeq'].get()

		for string in self.note_pr['gflag']: 
			menu1.add_command(label=string,command=lambda value=string: self._option_menu_update(self.note_pr['runeq'],value,self._sync_gui_opt,1))
		for string in self.note_nb['gflag']: 
			menu2.add_command(label=string,command=lambda value=string: self._option_menu_update(self.note_nb['runeq'],value,self._sync_gui_opt,3))
		for string in self.note_ef['gflag']: 
			menu3.add_command(label=string,command=lambda value=string: self._option_menu_update(self.note_ef['runeq'],value,self._sync_gui_opt,4))
		if selec1 in self.note_pr['gflag']: self.note_pr['runeq'].set(selec1)
		if selec2 in self.note_nb['gflag']: self.note_nb['runeq'].set(selec2)
		if selec3 in self.note_ef['gflag']: self.note_ef['runeq'].set(selec3)

		return

	def _option_menu_update(self,var1,value,func,funcvar):

		var1.set(value)
		if not funcvar == None:
			func(funcvar)
		else:
			func()
		return

	def _ask_file_name(self,inds):

		inputf = askopenfilename()
		if (len(inputf) == 0):      return

		self.note_nb['e%i'%inds].delete(0,'end')
		self.note_nb['e%i'%inds].insert(10,inputf)
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

		self.home       		= dict()
		self.note_in    		= dict()
		self.note_pr  	        = dict()
		self.note_pr['exclude'] = dict()
		self.note_ms 	 	    = dict()
		self.note_nb  		    = dict()
		self.note_ef  		    = dict()
		self.note_ec  		    = dict()

		self.kfile = dict()

		self.ver_plot  = 0
		self.ver_plots = dict()
		self.ver_plots[1] = []
		self.ver_plots[2] = []
		self.ver_plots[3] = []
		self.ver_plots[4] = []
		self.ver_plott = []

		self.post_opt  = dict()

		self.curr_page = 0
		self.bdata = dict()
		self.efdata = dict()

		self.func_list = ['CORE','MTANH','PTANH','EPED','SPLINE','EPED2','NSPLINE']
		self.initialise_opt = True
		self.nbeam = 0

		for i in range(4): self.__dict__['t_close%i'%i] = True


		return

	def _initialise_variables(self):

		self.note_pr['gflag'] = ['EFIT01']
		self.note_nb['gflag'] = ['EFIT01']
		self.note_ef['gflag'] = ['EFIT01']

		self.note_pr['runeq_old'] = 'None'
		self.note_nb['runeq_old'] = 'None'
		self.note_ef['runeq_old'] = 'None'

		self.figure['size']['home'] = (11,7.6)
		self.figure['size']['plot'] = (3.8,0.9)

		self.note_in['shot'] = ''
		self.note_in['times']= []
		self.note_in['trange'] = ''
		self.note_in['delt'] = ''
		self.note_in['tmin'] = 0.
		self.note_in['tmax'] = 1.e5
		self.note_in['nch']  = 25

		self.note_pr['rawtype'] = tk.StringVar()
		self.note_pr['fittype'] = tk.StringVar()
		self.note_pr['rawtype'].set('TS-TE')
		self.note_pr['fittype'].set('TE')

		flag1 = ['TS-TE','TS-NE','TSE-TE','TSE-NE','CES-TI','CES-VT','ECE']
		for flag in flag1: self.note_pr['exclude'][flag] = ''

		self.note_pr['TS_TCM'] = '1.0'
		self.note_pr['TS_TEM'] = '1.0'
		self.note_pr['TS_NCM'] = '1.0'
		self.note_pr['TS_NEM'] = '1.0'
		self.note_pr['TSEP']   = '0.15'
		self.note_pr['TWIDTH'] = '0.04'
		self.note_pr['NCMAX']  = '1.5'
		self.note_pr['NCMIN']  = '0.5'
		self.note_pr['NEMAX']  = '1.5'
		self.note_pr['NEMIN']  = '0.5'
		self.note_pr['NCN']    = '6'
		self.note_pr['NEN']    = '6'
		self.note_pr['TAVG']   = dict()
		self.note_pr['ZEFF']   = '2.0'
		self.note_pr['DACRIT'] = '1.0'
		self.note_pr['DUTY']   = '30'
		for flag in ['TS','CES','ECE','TCI']: self.note_pr['TAVG'][flag] = '100'
		self.note_pr['TAVG']['TCI'] = '100'

		flag1 = ['INT1','INT2','TCI1','TCI2','TCI3','TCI4','TCI5']
		for flag in flag1: self.note_pr[flag] = '0.0'
		self.note_pr['TAVG']['ECE'] = '0.1'
		self.note_pr['TCI4'] = '0.3'
		self.note_pr['TCI2'] = '0.5'
		self.note_pr['TCI3'] = '0.3'
			
		self.note_pr['FORCE_FIT'] = tk.IntVar()
		self.note_pr['HMODE']     = tk.IntVar()
		self.note_pr['HMODE'].set(1)
		self.note_pr['ASHIFT']    = tk.IntVar()
		self.note_pr['ASCALE']    = tk.IntVar()
		self.note_pr['ASCALE1D']  = tk.IntVar()
		self.note_pr['runeq']     = tk.StringVar()
		self.note_pr['runeq'].set('EFIT01')

		self.note_pr['ASHIFT'].set(1)
		self.note_pr['FORCE_FIT'].set(1) ##- PRESET

		self.note_ms['TAVG'] = '50'
		for k in range(1,51):
			self.note_ms['msec%02i'%k] = tk.IntVar()
			self.note_ms['mses%02i'%k] = '0.02'
			self.note_ms['msec%02i'%k].set(0)
			
		self.note_ms['use_exp_sgam'] = tk.IntVar()
		self.note_ms['xmin'] = 1.7
		self.note_ms['xmax'] = 2.4

		self.note_nb['runeq'] = tk.StringVar()
		self.note_nb['plottype'] = tk.StringVar()
		self.note_nb['BSMODEL'] = tk.StringVar()

		self.note_nb['runeq'].set('EFIT01')
		self.note_nb['plottype'].set('PBEAM')
		self.note_nb['BSMODEL'].set('csauter')

		self.note_nb['HAGCORE'] = tk.IntVar()
		self.note_nb['HAGCORE'].set(1)

		self.note_nb['bslist'] = ['csauter','sauter','neo','hager']
		self.note_nb['HAGCOREPSI'] = '0.3'
		self.note_nb['CORENEO'] = '0.1'
		self.note_nb['NSCALE'] = '0.0'
		self.note_nb['BSMULTI'] = '1.0'
		self.note_nb['NS'] = '70'
		self.note_nb['NT'] = '70'
		self.note_nb['MAPNS'] = '200'
		self.note_nb['MAPNT'] = '200'
		self.note_nb['IPCRIT'] = '1.e-4'
		self.note_nb['BSCRIT'] = '1.e-4'

		self.note_nb['NPROC'] = '28'
		self.note_nb['RUNSTEP']= '20'
		self.note_nb['RUNAVG'] = '0'
		self.note_nb['RUNDT'] = '0.01'
		self.note_nb['AVGDT'] = '0.01'
		self.note_nb['MFILE'] = ''
		self.note_nb['SFILE'] = ''
		self.note_nb['IFILE'] = ''
		self.note_nb['B-DIFF'] = '0.0'

		self.note_nb['plottype_old'] = ''
		self.note_nb['ymax'] = 3.0
		self.note_nb['ymin'] = 0.
		self.note_nb['xmin'] = 0.
		self.note_nb['xmax'] = 1.

		self.note_ef['runeq'] = tk.StringVar()
		self.note_ef['plottype'] = tk.StringVar()
		self.note_ef['efitgrid'] = tk.StringVar()

		self.note_ef['runeq'].set('EFIT01')
		self.note_ef['plottype'].set('PRESS')
		self.note_ef['efitgrid'].set('257')

		self.note_ef['HMODE']  = tk.IntVar()
		self.note_ef['JCONST'] = tk.IntVar()
		self.note_ef['QFIX']   = tk.IntVar()
		self.note_ef['ERMSE']  = tk.IntVar()
		self.note_ef['PSCALE'] = tk.IntVar()
		self.note_ef['PSCALE2'] = tk.IntVar()

		self.note_ef['HMODE'].set(1)
		self.note_ef['JCONST'].set(1)
		self.note_ef['PSCALE'].set(1)

		self.note_ef['MCOILN'] = 1
		self.note_ef['FLOOPN'] = 1
		self.note_ef['KLOAD']  = False

		for k in range(150): self.note_ef['MCOIL%03i'%(k+1)] = tk.IntVar()
		for k in range(100): self.note_ef['FLOOP%03i'%(k+1)] = tk.IntVar()

		self.note_ef['RITER']   = '0'
		self.note_ef['MAXITER'] = '101'
		self.note_ef['RELAX']   = '0.5'
		self.note_ef['CONVERG'] = '2.e-4'
		self.note_ef['FWTCUR']  = '2.0'
		self.note_ef['QCONST']  = '1.'
		self.note_ef['SUPP0']   = '0.'

		self.note_ef['SNAPFILE'] = ''
		self.note_ef['PCNN']  = '2'
		self.note_ef['PENN']  = '5'
		self.note_ef['FCNN']  = '1'
		self.note_ef['FENN']  = '5'
		self.note_ef['JCNN']  = '2'
		self.note_ef['JENN']  = '12'
		self.note_ef['JSTART']= '0.6'

		self.note_ef['plottype_old'] = ''
		self.note_ef['ymax'] = 3.0
		self.note_ef['ymin'] = 0.0		
		self.note_ef['xmax'] = 1.
		self.note_ef['xmin'] = 0.
		
		return

	def _sync_gui_opt(self,sync_type=-1):

		#PROFILE
		if (sync_type ==1 or sync_type <0):
			if (self.note_pr['runeq'].get() == self.note_pr['runeq_old'] and not self.curr_page == 1): return
			print('>>> Sync Profile...')
			self.note_pr['l1'].delete(0,'end')
			self.note_pr['l2'].delete(0,'end')
			tlist = list(self.note_in['l2'].get(0,'end'))
			if len(tlist)== 0: return
			for item in tlist:
				if self.note_pr['runeq'].get() == 'EFIT01':
					self.note_pr['l1'].insert('end',item)
				else:
					tt = int(float(item.split()[0]))
					filen = 'EFITS/%s/g%06i.%06i'%(self.note_pr['runeq'].get(),self.note_in['shot'],tt)
					if os.path.isfile(filen): self.note_pr['l1'].insert('end',item)

			if not (self.note_pr['runeq'].get() == self.note_pr['runeq_old']):
				self.note_pr['runeq_old'] = self.note_pr['runeq'].get()
				file1 = 'PROFILES/result_multi.save_%s'%self.note_pr['runeq'].get()
				if os.path.isfile(file1):
					f = open(file1,'rb'); self.mfit_opt = pickle.load(f); f.close()
				else: return

				tlist = list(self.note_pr['l1'].get(0,'end'))
				tlist2 = copy.deepcopy(tlist)
				for k in range(len(tlist)):
					tt = float(tlist[k].split()[0])
					tlist2[k] = int(tt)							

				self.note_pr['l2'].delete(0,'end')
				for time in self.mfit_opt['times']:
					try: 
						ind = tlist2.index(time)
						self.note_pr['l2'].insert('end',tlist[ind])
					except: pass	

				if sync_type == 1: self._draw_0d_profiles()

		#MSE
		if (sync_type == 2 or sync_type < 0):
			self.note_ms['l1'].delete(0,'end')
			tlist = list(self.note_in['l2'].get(0,'end'))
			for time in tlist:
				self.note_ms['l1'].insert('end',time)

		#BEAM
		if (sync_type == 3 or sync_type < 0):
			if (self.note_nb['runeq'].get() == self.note_nb['runeq_old']): return
			print('>>> Sync H&CD ...')
			self.note_nb['runeq_old'] = self.note_nb['runeq'].get()
			self.note_nb['l1'].delete(0,'end')
			if not (self.note_nb['runeq'].get() == self.note_pr['runeq_old']):
				self.note_pr['runeq_old'] = self.note_nb['runeq'].get()
				file1 = 'PROFILES/result_multi.save_%s'%self.note_nb['runeq'].get()
				if os.path.isfile(file1):
					f = open(file1,'rb'); self.mfit_opt = pickle.load(f); f.close()

			tlist = list(self.note_in['l2'].get(0,'end'))
			tlist2 = copy.deepcopy(tlist)
			for k in range(len(tlist)):
				tt = float(tlist[k].split()[0])
				tlist2[k] = int(tt)							

			self.note_nb['l1'].delete(0,'end')
			for time in self.mfit_opt['times']:
				try: 
					ind = tlist2.index(time)
					self.note_nb['l1'].insert('end',tlist[ind])
				except: pass
			self._load_beam()

		#EFITS
		if (sync_type == 4 or sync_type < 0):
			self._check_efit_input()
			self._find_efit_run()
			if (self.note_ef['runeq'].get() == self.note_ef['runeq_old']): return
			print('>>> Sync EFIT ...')		
			self.note_ef['runeq_old'] = self.note_ef['runeq'].get()	
			self._load_efit_runs()

		return

	def _run_exit(self):

		self.root.destroy()
		return

	def _save_run(self):

		note_in = dict()
		note_in['shot'] = self.note_in['e1'].get()
		note_in['times'] = list(self.note_in['l2'].get(0,'end'))

		note_pr = dict()
		note_pr['exclude'] = dict()
		flag1 = ['TS-TE','TSE-TE','TS-NE','TSE-NE','CES-TI','CES-VT','ECE']
		for k in range(7): note_pr['exclude'][flag1[k]] = self.note_pr['e%i'%(k+5)].get()
		flag1 = ['TSEP','TWIDTH','ZEFF','TS_TCM','TS_TEM','TS_NCM','TS_NEM','NCMAX','NCMIN','NEMAX','NEMIN','NCN','NEN','DACRIT','DUTY']
		for k in range(15): note_pr[flag1[k]] = self.note_pr['e%i'%(k+12)].get()
		note_pr['TAVG'] = dict()
		flag1 = ['TS','CES','ECE','TCI']
		for k in range(4): note_pr['TAVG'][flag1[k]] = self.note_pr['e%i'%(k+27)].get()
		flag1 = ['INT1','INT2','TCI1','TCI2','TCI3','TCI4','TCI5']
		for k in range(7): note_pr[flag1[k]] = self.note_pr['e%i'%(k+31)].get()
		
		flag1 = ['FORCE_FIT','HMODE','ASHIFT','ASCALE','ASCALE1D']
		for flag in flag1: note_pr[flag] = self.note_pr[flag].get()

		note_ms = dict()
		note_ms['TAVG']	= self.note_ms['e1'].get()
		note_ms['NCH']  = len(self.mse_data['ch'])
		for k in range(len(self.mse_data['ch'])):
			note_ms['msec%02i'%(k+1)] = self.note_ms['msec%02i'%(k+1)].get()
			note_ms['mses%02i'%(k+1)] = self.note_ms['e%i'%(k+4)].get()		
		note_ms['use_exp_sgam'] = self.note_ms['use_exp_sgam'].get()

		note_nb = dict()
		flag1 = ['HAGCOREPSI','CORENEO','NSCALE','BSMULTI','NS','NT','MAPNS','MAPNT','IPCRIT','BSCRIT']
		for k in range(10): note_nb[flag1[k]] = self.note_nb['e%i'%(k+3)].get()
		flag1 = ['NPROC','RUNSTEP','RUNAVG','RUNDT','AVGDT','B-DIFF','MFILE','SFILE','IFILE']
		for k in range(9): note_nb[flag1[k]] = self.note_nb['e%i'%(k+13)].get()

		note_nb['HAGCORE'] = self.note_nb['HAGCORE'].get()
		note_nb['BSMODEL'] = self.note_nb['BSMODEL'].get()

		note_ef = dict()
		flag1 = ['PCNN','PENN','FCNN','FENN','JCNN','JENN','JSTART']
		for k in range(7): note_ef[flag1[k]] = self.note_ef['e%i'%(k+3)].get()
		flag1 = ['RITER','MAXITER','RELAX','CONVERG','FWTCUR','SUPP0','QCONST']
		for k in range(7): note_ef[flag1[k]] = self.note_ef['e%i'%(k+10)].get()

		flag1 = ['HMODE','JCONST','QFIX','ERMSE','PSCALE','PSCALE2','efitgrid']
		for flag in flag1: note_ef[flag] = self.note_ef[flag].get()

		for k in range(150): note_ef['MCOIL%03i'%(k+1)] = self.note_ef['MCOIL%03i'%(k+1)].get()
		for k in range(100): note_ef['FLOOP%03i'%(k+1)] = self.note_ef['FLOOP%03i'%(k+1)].get()

		note_pr['runeq'] = self.note_pr['runeq'].get()
		note_nb['runeq'] = self.note_nb['runeq'].get()
		note_ef['runeq'] = self.note_ef['runeq'].get()

		f = open('save.dat','wb')
		pickle.dump([note_in,note_pr,note_ms,note_nb,note_ef],f)
		f.close()

		return

	def _load_run(self):

		if not os.path.isfile('save.dat'): return
		f = open('save.dat','rb')
		[note_in,note_pr,note_ms,note_nb,note_ef] = pickle.load(f)
		f.close()

		for key in note_in.keys(): self.note_in[key] = copy.deepcopy(note_in[key])
		for key in note_pr.keys(): self.note_pr[key] = copy.deepcopy(note_pr[key])
		for key in note_ms.keys(): self.note_ms[key] = copy.deepcopy(note_ms[key])
		for key in note_nb.keys(): self.note_nb[key] = copy.deepcopy(note_nb[key])
		for key in note_ef.keys(): self.note_ef[key] = copy.deepcopy(note_ef[key])

		flag1 = ['FORCE_FIT','HMODE','ASHIFT','ASCALE','ASCALE1D']
		for flag in flag1: 
			try:
				self.note_pr[flag] = tk.IntVar()
				self.note_pr[flag].set(note_pr[flag])
			except: pass
		

		for k in range(self.note_in['nch']):
			self.note_ms['msec%02i'%(k+1)] = tk.IntVar()
			self.note_ms['msec%02i'%(k+1)].set(note_ms['msec%02i'%(k+1)])	
		self.note_ms['use_exp_sgam'] = tk.IntVar()
		self.note_ms['use_exp_sgam'].set(note_ms['use_exp_sgam'])

		try:
			self.note_nb['HAGCORE'] = tk.IntVar()
			self.note_nb['HAGCORE'].set(note_nb['HAGCORE'])
			self.note_nb['BSMODEL'] = tk.StringVar()
			self.note_nb['BSMODEL'].set(note_nb['BSMODEL'])
		except: pass

		try:
			flag1 = ['HMODE','JCONST','QFIX','ERMSE','PSCALE','PSCALE2']
			for flag in flag1:
				self.note_ef[flag] = tk.IntVar()
				self.note_ef[flag].set(note_ef[flag])
		except: pass

		try:	
			self.note_ef['efitgrid'] = tk.StringVar()
			self.note_ef['efitgrid'].set(note_ef['efitgrid'])
		except: pass

#		self.note_ef['HMODE'].set(0)

		try:
			for k in range(150): 
				self.note_ef['MCOIL%03i'%(k+1)] = tk.IntVar()
				self.note_ef['MCOIL%03i'%(k+1)].set(note_ef['MCOIL%03i'%(k+1)])
			for k in range(100): 
				self.note_ef['FLOOP%03i'%(k+1)] = tk.IntVar()
				self.note_ef['FLOOP%03i'%(k+1)].set(note_ef['FLOOP%03i'%(k+1)])			
		except: pass

		self.note_pr['runeq'] = tk.StringVar()
		self.note_nb['runeq'] = tk.StringVar()
		self.note_ef['runeq'] = tk.StringVar()

		try: self.note_pr['runeq'].set(note_pr['runeq'])
		except: pass
		try: self.note_nb['runeq'].set(note_nb['runeq'])
		except: pass
		try: self.note_ef['runeq'].set(note_ef['runeq'])
		except: pass		

		if self.note_pr['runeq'].get() == '': self.note_pr['runeq'].set('EFIT01')
		if self.note_nb['runeq'].get() == '': self.note_nb['runeq'].set('EFIT01')
		if self.note_ef['runeq'].get() == '': self.note_ef['runeq'].set('EFIT01')

		#self.note_pr['runeq_old'] = self.note_pr['runeq'].get()
		#self.note_nb['runeq_old'] = self.note_nb['runeq'].get()
		#self.note_ef['runeq_old'] = self.note_ef['runeq'].get()

		return

	def __init__(self,single_only=False):

		self.lmfit_mod = False
		self.input_ready = False
		self.input_load = False
		self.single_only = single_only

		return

if __name__ == "__main__":

	import fgefit

	try:
		os.mkdir('PROFILES')
	except: pass
	try: os.mkdir('MDS')
	except: pass

	print(' ---------------------------------------------------------------')
	print('||               Fast profile Fitting tool Ver %s             ||'%version['fgefit'])
	print('||               Kinetic Profile & EFIT generator              ||')
	print('%s'%author['fgefit'])
	print('%s'%comment['fgefit'])
	print(' ---------------------------------------------------------------\n')

	try: flag = sys.argv[1]
	except: flag = ''

	fgefit = fgefit.fgefittool()

	fgefit.root = tk.Tk()
	fgefit.root.title('FGEFIT')
	fgefit.fgefit()
