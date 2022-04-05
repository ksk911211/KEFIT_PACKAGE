#!/usr/local/anaconda3/bin/python3
import os,sys
import tkinter as tk
import numpy as np
import eped_fit
from shutil import move, copyfile
from tkinter.filedialog import askopenfilename,asksaveasfilename

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

import matplotlib.pyplot as plt
import copy

currdir = os.getcwd()

class eped_gtool:

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
			
	def trans_vars2(self,var1,type):
	
		if (type == 4):
			if (var1 == 1):
				return  True
			else:
				return False
		elif (type == 2):
			return (float(var1))
		elif (type == 1):
			return (int(float(var1)))
		elif (type == 3):
			return var1
		elif (type == 5):
		
			for i in range(len(self.var_list)):
			
				if (var1.lower() == self.var_list[i].lower()):
				
					break
			return (i+1)
	
	def button_func1(self):
		
		self.write_fit_opt()
		self.fit.read_namelist('fit_opt')
		if (self.lmfit_mod):
			self.set_lmfit_param()

		self.fit.main_run()
		self.initialise_lmfit_param()

		if (self.fit.use_ti_width):
			self.e25.delete(0, 'end')
			self.e26.delete(0, 'end')
			self.StrVar13.set(round(self.fit.popt_ti[2],5))
			self.StrVar14.set(round(self.fit.popt_ti[2],5))
			self.e25.insert(0,self.StrVar13.get())
			self.e26.insert(0,self.StrVar14.get())
			self.CheckVar23.set(1)
			self.CheckVar24.set(1)
		else:
			self.CheckVar23.set(0)
			self.CheckVar24.set(0)			


#		self.fit.draw_plot(self.fig)
		self.fig.canvas.draw_idle()
		self.fit.draw_plot(self.fig)
		self.did_fit = True

		try:
			if not (self.t3_close):
				self.fig2.canvas.draw_idle()
				self.fit.draw_plot(self.fig2,True)
		except:
			self.t3_close = True
			pass
		
		return

	def button_func2(self):
		try:
			os.mkdir('PROFILES')
		except:
			pass

		if not (self.fit.no_ne):	self.StrVar46.set('PROFILES/NE_non_scaled.dat')
		if not (self.fit.no_te):	self.StrVar45.set('PROFILES/TE_fit.dat')
		if not (self.fit.no_ti):	self.StrVar47.set('PROFILES/TI_fit.dat')
		self.fit.write_kinprof()
		print('>>> Fitting results saved...')
		self.write_fit_opt()
		self.save = 'fit_opt_save'

		options = {}
		options['initialdir'] = currdir
		try:	self.save = asksaveasfilename(**options)
		except:		pass
		if (self.save == ''):
			return

		if (self.lmfit_mod):
			self.fit.lmfit_write_param(self.fit.param)

		copyfile('fit_opt','fit_opt_save')
		try:
			copyfile('fit_opt',self.save)
			print('>>> Save file to %s...'%self.save)
		except:
			print('>>> Save file to fit_opt...')
		try:
			copyfile('fit_opt_param','fit_opt_save'+'_param')
			copyfile('fit_opt_param',self.save+'_param')
		except:
			print('>>> Default lmfit values are used, do not save lmfit_params...')
		self.input = 'fit_opt_save'

		self.fit.read_namelist('fit_opt')
		self.fit.read_kinprof_kprofile()
		self.fit.ne_fit2 = np.copy(self.fit.ne_fit2 / self.fit.scaled_density)
		if not (self.fit.no_ne):	self.fit.den_scale(2)
		self.fig.canvas.draw_idle()
		self.fit.draw_plot(self.fig)

		try:
			if not (self.t3_close):
				self.fig2.canvas.draw_idle()
				self.fit.draw_plot(self.fig2,True)
		except:
			self.t3_close = True
			pass		

		return
		
	def button_func3(self):

		self.input = askopenfilename()
		if not (len(self.input) == 0):
			self.reopen = True
			self.file_fixed = True

			self.te_file = self.fit.te_file
			self.ne_file = self.fit.ne_file
			self.ti_file = self.fit.ti_file
			self.vt_file = self.fit.vt_file

			self.te_dat_file = self.fit.te_dat_file
			self.ne_dat_file = self.fit.ne_dat_file
			self.ti_dat_file = self.fit.ti_dat_file
			self.vt_dat_file = self.fit.vt_dat_file

			self.target_density = self.fit.target_density
			self.eqdsk_name = self.fit.eqdsk_name

			self.root.quit()
		return

	def button_func4(self):
		self.reopen = False
		print('>>> Exiting....')
		self.root.quit()
		return

	def button_func5(self):

		self.input = 'fit_opt_save'
		self.reopen = True
		self.root.quit()
		return

	def button_func6(self):
		
		self.input = 'fit_opt~'
		self.reopen = True
		self.root.quit()

		return

	def button_func7(self):

		if not self.t_close:
			print('>>> Window is already opened...')
			return
		
		self.t = tk.Toplevel(self.root)
		self.t.wm_title("Input file dirs")
		self.t_close = False
		self.t.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(1))
		l = tk.Label(self.t, text="This is window #")

		self.l2 = tk.Label(self.t, text = '====== EQDSK FILE ======',justify='center')
		self.l2.grid(row=1,column=1)
		self.l2 = tk.Label(self.t, text=" EQDSK FILE")
		self.l2.grid(row=2, column=0)
		self.e200 = tk.Entry(self.t,width=40,justify='center')
		self.e200.insert(0,self.StrVar44.get())
		self.e200.grid(row=2, column=1, padx=5)

		b1 = tk.Button(self.t, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar44,self.e200),height = 1,width = 3)
		b1.grid(row=2, column=2)		

		self.l2 = tk.Label(self.t, text = '======= FIT FILE =======',justify='center')
		self.l2.grid(row=3,column=1)		

		self.l2 = tk.Label(self.t, text=" TE FILE [eV]")
		self.l2.grid(row=4, column=0)
		self.l2 = tk.Label(self.t, text=" NE FILE [(18)/m3]")
		self.l2.grid(row=5, column=0)
		self.l2 = tk.Label(self.t, text=" TI FILE [eV]")
		self.l2.grid(row=6, column=0)		

		self.e201 = tk.Entry(self.t,width=40,justify='center')
		self.e201.insert(0,self.StrVar45.get())
		self.e201.grid(row=4, column=1, padx=5)
		self.e202 = tk.Entry(self.t,width=40,justify='center')
		self.e202.insert(0,self.StrVar46.get())
		self.e202.grid(row=5, column=1, padx=5)
		self.e203 = tk.Entry(self.t,width=40,justify='center')
		self.e203.insert(0,self.StrVar47.get())
		self.e203.grid(row=6, column=1, padx=5)

		b2 = tk.Button(self.t, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar45,self.e201),height = 1,width = 3)
		b2.grid(row=4, column=2)	
		b3 = tk.Button(self.t, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar46,self.e202),height = 1,width = 3)
		b3.grid(row=5, column=2)			
		b4 = tk.Button(self.t, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar47,self.e203),height = 1,width = 3)
		b4.grid(row=6, column=2)		

		self.l2 = tk.Label(self.t, text = '======= RAW FILE =======',justify='center')
		self.l2.grid(row=8,column=1)		

		self.l2 = tk.Label(self.t, text=" TE FILE [eV]")
		self.l2.grid(row=9, column=0)
		self.l2 = tk.Label(self.t, text=" NE FILE [(18)/m3")
		self.l2.grid(row=10, column=0)
		self.l2 = tk.Label(self.t, text=" TI FILE [eV]")
		self.l2.grid(row=11, column=0)		

		self.e205 = tk.Entry(self.t,width=40,justify='center')
		self.e205.insert(0,self.StrVar49.get())
		self.e205.grid(row=9, column=1, padx=5)
		self.e206 = tk.Entry(self.t,width=40,justify='center')
		self.e206.insert(0,self.StrVar50.get())
		self.e206.grid(row=10, column=1, padx=5)
		self.e207 = tk.Entry(self.t,width=40,justify='center')
		self.e207.insert(0,self.StrVar51.get())
		self.e207.grid(row=11, column=1, padx=5)

		b6 = tk.Button(self.t, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar49,self.e205),height = 1,width = 3)
		b6.grid(row=9, column=2)	
		b7 = tk.Button(self.t, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar50,self.e206),height = 1,width = 3)
		b7.grid(row=10, column=2)		
		b8 = tk.Button(self.t, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar51,self.e207),height = 1,width = 3)
		b8.grid(row=11, column=2)	

		b10 = tk.Button(self.t, text="DONE", bg = "gray",command=lambda: self.button_func7a(),height = 1,width = 3)
		b10.grid(row=13, column=0,columnspan=3, pady=10)		

		return

	def button_func7a(self):

		self.t.destroy()
		self.t_close = True
		return

	def button_func8(self,var,var2):

		self.input_file = askopenfilename()
		if not(self.input_file == ''):
			var.set(self.input_file)
			var2.delete(0, 'end')
			var2.insert(0,self.input_file)

		return

	def button_func9(self,var):
	
			self.write_fit_opt('fit_opt',True)
			var.destroy()
			return

	def button_func10(self):

		self.lmfit_mod = True

		if not self.t2_close:
			print('>>> LMfit Param Window is already opened...')
			return
	
		self.t2 = tk.Toplevel(self.root)
		self.t2.wm_title("LMFIT Params")
		self.t2_close = False		
		self.t2.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(2))		

		self.l2 = tk.Label(self.t2, text = '======== TE PARAMS ========',justify='center')
		self.l2.grid(row=1,column=1,columnspan=4,pady=15)
		self.l2 = tk.Label(self.t2, text=" Var")
		self.l2.grid(row=2, column=0)
		self.l2 = tk.Label(self.t2, text=" Value")
		self.l2.grid(row=2, column=1)
		self.l2 = tk.Label(self.t2, text=" min")
		self.l2.grid(row=2, column=2)
		self.l2 = tk.Label(self.t2, text=" Max")
		self.l2.grid(row=2, column=3)
		self.l2 = tk.Label(self.t2, text=" Fix")
		self.l2.grid(row=2, column=4)

		self.l2 = tk.Label(self.t2, text = '======== NE PARAMS ========',justify='center')
		self.l2.grid(row=1,column=5,columnspan=4)
		self.l2 = tk.Label(self.t2, text=" Value")
		self.l2.grid(row=2, column=5)
		self.l2 = tk.Label(self.t2, text=" min")
		self.l2.grid(row=2, column=6)
		self.l2 = tk.Label(self.t2, text=" Max")
		self.l2.grid(row=2, column=7)
		self.l2 = tk.Label(self.t2, text=" Fix")
		self.l2.grid(row=2, column=8)

		self.l2 = tk.Label(self.t2, text=" a1")
		self.l2.grid(row=3, column=0)	
		self.l2 = tk.Label(self.t2, text=" a2")
		self.l2.grid(row=4, column=0)
		self.l2 = tk.Label(self.t2, text=" a3")
		self.l2.grid(row=5, column=0)
		self.l2 = tk.Label(self.t2, text=" a4")
		self.l2.grid(row=6, column=0)	
		self.l2 = tk.Label(self.t2, text=" a5")
		self.l2.grid(row=7, column=0)	
		self.l2 = tk.Label(self.t2, text=" a6")
		self.l2.grid(row=8, column=0)		
	
		self.e301 = self.make_box(self.t2,6,self.StrVar101,3,1,2,2)
		self.e302 = self.make_box(self.t2,6,self.StrVar102,4,1,2,2)
		self.e303 = self.make_box(self.t2,6,self.StrVar103,5,1,2,2)
		self.e304 = self.make_box(self.t2,6,self.StrVar104,6,1,2,2)
		self.e305 = self.make_box(self.t2,6,self.StrVar105,7,1,2,2)
		self.e306 = self.make_box(self.t2,6,self.StrVar106,8,1,2,2)
		self.e307 = self.make_box(self.t2,6,self.StrVar107,3,2,2,2)
		self.e308 = self.make_box(self.t2,6,self.StrVar108,4,2,2,2)
		self.e309 = self.make_box(self.t2,6,self.StrVar109,5,2,2,2)
		self.e310 = self.make_box(self.t2,6,self.StrVar110,6,2,2,2)
		self.e311 = self.make_box(self.t2,6,self.StrVar111,7,2,2,2)
		self.e312 = self.make_box(self.t2,6,self.StrVar112,8,2,2,2)
		self.e313 = self.make_box(self.t2,6,self.StrVar113,3,3,2,2)
		self.e314 = self.make_box(self.t2,6,self.StrVar114,4,3,2,2)
		self.e315 = self.make_box(self.t2,6,self.StrVar115,5,3,2,2)
		self.e316 = self.make_box(self.t2,6,self.StrVar116,6,3,2,2)
		self.e317 = self.make_box(self.t2,6,self.StrVar117,7,3,2,2)
		self.e318 = self.make_box(self.t2,6,self.StrVar118,8,3,2,2)

		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar31)
		self.e300.grid(row=3, column=4)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar32)
		self.e300.grid(row=4, column=4)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar33)
		self.e300.grid(row=5, column=4)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar34)
		self.e300.grid(row=6, column=4)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar35)
		self.e300.grid(row=7, column=4)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar36)
		self.e300.grid(row=8, column=4)

		self.e319 = self.make_box(self.t2,6,self.StrVar119,3,5,2,2)
		self.e320 = self.make_box(self.t2,6,self.StrVar120,4,5,2,2)
		self.e321 = self.make_box(self.t2,6,self.StrVar121,5,5,2,2)
		self.e322 = self.make_box(self.t2,6,self.StrVar122,6,5,2,2)
		self.e323 = self.make_box(self.t2,6,self.StrVar123,7,5,2,2)
		self.e324 = self.make_box(self.t2,6,self.StrVar124,8,5,2,2)
		self.e325 = self.make_box(self.t2,6,self.StrVar125,3,6,2,2)
		self.e326 = self.make_box(self.t2,6,self.StrVar126,4,6,2,2)
		self.e327 = self.make_box(self.t2,6,self.StrVar127,5,6,2,2)
		self.e328 = self.make_box(self.t2,6,self.StrVar128,6,6,2,2)
		self.e329 = self.make_box(self.t2,6,self.StrVar129,7,6,2,2)
		self.e330 = self.make_box(self.t2,6,self.StrVar130,8,6,2,2)
		self.e331 = self.make_box(self.t2,6,self.StrVar131,3,7,2,2)
		self.e332 = self.make_box(self.t2,6,self.StrVar132,4,7,2,2)
		self.e333 = self.make_box(self.t2,6,self.StrVar133,5,7,2,2)
		self.e334 = self.make_box(self.t2,6,self.StrVar134,6,7,2,2)
		self.e335 = self.make_box(self.t2,6,self.StrVar135,7,7,2,2)
		self.e336 = self.make_box(self.t2,6,self.StrVar136,8,7,2,2)

		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar37)
		self.e300.grid(row=3, column=8)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar38)
		self.e300.grid(row=4, column=8)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar39)
		self.e300.grid(row=5, column=8)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar40)
		self.e300.grid(row=6, column=8)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar41)
		self.e300.grid(row=7, column=8)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar42)
		self.e300.grid(row=8, column=8)

		self.l2 = tk.Label(self.t2, text = '======== TI PARAMS ========',justify='center')
		self.l2.grid(row=11,column=1,columnspan=4,pady=15)
		self.l2 = tk.Label(self.t2, text=" Var")
		self.l2.grid(row=12, column=0)
		self.l2 = tk.Label(self.t2, text=" Value")
		self.l2.grid(row=12, column=1)
		self.l2 = tk.Label(self.t2, text=" min")
		self.l2.grid(row=12, column=2)
		self.l2 = tk.Label(self.t2, text=" Max")
		self.l2.grid(row=12, column=3)
		self.l2 = tk.Label(self.t2, text=" Fix")
		self.l2.grid(row=12, column=4)

		self.l2 = tk.Label(self.t2, text=" a1")
		self.l2.grid(row=13, column=0)	
		self.l2 = tk.Label(self.t2, text=" a2")
		self.l2.grid(row=14, column=0)
		self.l2 = tk.Label(self.t2, text=" a3")
		self.l2.grid(row=15, column=0)
		self.l2 = tk.Label(self.t2, text=" a4")
		self.l2.grid(row=16, column=0)	
		self.l2 = tk.Label(self.t2, text=" a5")
		self.l2.grid(row=17, column=0)	
		self.l2 = tk.Label(self.t2, text=" a6")
		self.l2.grid(row=18, column=0)			

	
		self.e337 = self.make_box(self.t2,6,self.StrVar137,13,1,2,2)
		self.e338 = self.make_box(self.t2,6,self.StrVar138,14,1,2,2)
		self.e339 = self.make_box(self.t2,6,self.StrVar139,15,1,2,2)
		self.e340 = self.make_box(self.t2,6,self.StrVar140,16,1,2,2)
		self.e341 = self.make_box(self.t2,6,self.StrVar141,17,1,2,2)
		self.e342 = self.make_box(self.t2,6,self.StrVar142,18,1,2,2)
		self.e343 = self.make_box(self.t2,6,self.StrVar143,13,2,2,2)
		self.e344 = self.make_box(self.t2,6,self.StrVar144,14,2,2,2)
		self.e345 = self.make_box(self.t2,6,self.StrVar145,15,2,2,2)
		self.e346 = self.make_box(self.t2,6,self.StrVar146,16,2,2,2)
		self.e347 = self.make_box(self.t2,6,self.StrVar147,17,2,2,2)
		self.e348 = self.make_box(self.t2,6,self.StrVar148,18,2,2,2)
		self.e349 = self.make_box(self.t2,6,self.StrVar149,13,3,2,2)
		self.e350 = self.make_box(self.t2,6,self.StrVar150,14,3,2,2)
		self.e351 = self.make_box(self.t2,6,self.StrVar151,15,3,2,2)
		self.e352 = self.make_box(self.t2,6,self.StrVar152,16,3,2,2)
		self.e353 = self.make_box(self.t2,6,self.StrVar153,17,3,2,2)
		self.e354 = self.make_box(self.t2,6,self.StrVar154,18,3,2,2)

		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar43)
		self.e300.grid(row=13, column=4)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar44)
		self.e300.grid(row=14, column=4)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar45)
		self.e300.grid(row=15, column=4)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar46)
		self.e300.grid(row=16, column=4)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar47)
		self.e300.grid(row=17, column=4)
		self.e300 = tk.Checkbutton(self.t2,variable=self.CheckVar48)
		self.e300.grid(row=18, column=4)

		b11 = tk.Button(self.t2, text="RESET", bg = "violet",command=lambda: self.button_func12(),height = 1,width = 3)
		b11.grid(row=21, column=3,columnspan=2, pady=10)

		b12 = tk.Button(self.t2, text="LOAD", bg = "yellow",command=lambda: self.button_func13(),height = 1,width = 3)
		b12.grid(row=21, column=5,columnspan=2, pady=10)		

		b13 = tk.Button(self.t2, text="SAVE", bg = "lime",command=lambda: self.button_func11(),height = 1,width = 3)
		b13.grid(row=21, column=7,columnspan=2, pady=10)

		if not self.fit.param_change:
			self.fit.param2_old = copy.deepcopy(self.fit.param)
		self.interact_main_to_lmfit()
	
		return

	def button_func11(self):

		self.update_lmfit_param()
		self.set_lmfit_param()
		self.interact_lmfit_to_main()
		self.fit.param_change = True
		self.fit.param2 = copy.deepcopy(self.fit.param)
		self.t2.destroy()
		self.t2_close = True
		return

	def button_func12(self):

		self.fit.lmfit_init_params()
		self.initialise_lmfit_param()
		self.update_lmfit_param(False)

		return

	def button_func13(self):

		lmfit_input = askopenfilename()
		if (lmfit_input == ''):
			return

		fit_types = [self.fit.fit_type_te,self.fit.fit_type_ne,self.fit.fit_type_ti,self.fit.fit_type_vt]
		self.fit.lmfit_load_param(self.fit.param,fit_types,lmfit_input)
		self.initialise_lmfit_param()
		self.update_lmfit_param(False)

		return

	def button_func14(self,typev=1):

		if typev == 1:
			line = 'y=a1 +  a2 * (1 - x^a3)^a4'
		elif typev == 2:
			line = 'mtanh(x,b) = ((1.0 + b*x)*np.exp(x) - np.exp(-x)) / (np.exp(x) + np.exp(-x))\n'
			line = line + 'F(x) = (a2 - a1)/ 2.0 * (mtanh((a4-x)/2/a3,a5) + 1.0) + a1\n'
			line = line + 'y = F(x) + (a6-F(x))* np.exp((-x/(a7+0.1))^a7)'
		elif typev == 3:
			line = 'y = (a1 -  a2) * (1.0 + a3*x + a4 * x**2 + a5 * x**3) * 0.5 * (1.0 - np.tanh((x-a7)/a6)) + a2'
		elif typev == 4:
			line = 'y = a1 + a2*(np.tanh(1) - np.tanh((x - 1.0 + 0.5*(a3))*2/a3)) + a4 * (1 - (x/(1-a3))^a5)^a6'
		elif typev == 5:
			line = 'y = a1 + a2*(np.tanh((1-a7)*2.0/a3) - np.tanh((x - a7)*2/a3)) + a4 * (1 - (x/(a7-a3/2))^a5)^a6'

		print('=----------------------------------------------------------------------------=')
		print('FIT FUNC TYPE = %i'%typev)
		print(line)
		print('=----------------------------------------------------------------------------=')

		return

	def button_func15(self):

		if not (self.t3_close):
			print('>>> 2nd Prof Window is already opened...')
			return			

		self.t3 = tk.Toplevel(self.root) #self.root
		self.t3.title("2nd_derivative profiles")
		self.t3_close = False		
		self.t3.protocol('WM_DELETE_WINDOW', lambda: self.detect_close(3))

		self.fig2, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2,2,figsize=(12,7))
		ax1.set_title('TE_prime')
		ax2.set_title('NE_prime')
		ax3.set_title('TI_prime')
		ax4.set_title('VT_prime')
		if not (self.fit.use_rho):
			ax1.set_xlabel('$\psi_n$ [a.u]')
			ax2.set_xlabel('$\psi_n$ [a.u]')
			ax3.set_xlabel('$\psi_n$ [a.u]')
			ax4.set_xlabel('$\psi_n$ [a.u]')
		else:
			ax1.set_xlabel('$\\rho_t$ [a.u]')
			ax2.set_xlabel('$\\rho_t$ [a.u]')
			ax3.set_xlabel('$\\rho_t$ [a.u]')
			ax4.set_xlabel('$\\rho_t$ [a.u]')			
		self.fig2.tight_layout()

		#self.fig.subplots_adjust(left=0.05, right=0.993, top=0.97, bottom=0.07)

		self.canvas2 = FigureCanvasTkAgg(self.fig2,master=self.t3)
		self.plot_widget2 = self.canvas2.get_tk_widget()
		self.plot_widget2.grid(rowspan=20,row=1,column=20)

		toolbar_frame = tk.Frame(self.t3)
		toolbar_frame.grid(column=20,row=0)
		
		self.toolbar2 = NavigationToolbar2TkAgg(self.canvas2,toolbar_frame)

		self.fig2.canvas.draw_idle()
		if (self.did_fit):
			self.fit.draw_plot(self.fig2,True)

#		self.t3.mainloop()

		return

	def detect_close(self,win_num):

		#if tk.tkMessageBox.askokcancel("Quit", "Do you really wish to quit?"):
		if (win_num == 3):
			self.t3.destroy()
			self.t3_close = True
		elif (win_num == 2):
			self.t2.destroy()
			self.t2_close = True
		elif (win_num == 1):
			self.t.destroy()
			self.t_close = True
		return

	def make_box(self,tk_win,width,var,row,column,padx,pady):

		tk_obj = tk.Entry(tk_win,width=width,justify='center')
		tk_obj.insert(0,var.get())
		tk_obj.grid(row=row, column=column, padx=padx, pady=pady)

		return tk_obj

	def renew_box(self,var1,var2,copyf=True):

		if (copyf):
			var2.set(var1.get())
		var1.delete(0, 'end')
		var1.insert(0,var2.get())
		return

	def gui_fitopt(self,var):

		var.wm_title("Input file dirs")
		l = tk.Label(var, text="This is window #")

		self.l2 = tk.Label(var, text = '====== EQDSK FILE ======',justify='center')
		self.l2.grid(row=1,column=1)
		self.l2 = tk.Label(var, text=" EQDSK FILE")
		self.l2.grid(row=2, column=0)
		self.e210 = tk.Entry(var,width=40,justify='center')
		self.e210.insert(0,self.StrVar44.get())
		self.e210.grid(row=2, column=1, padx=5)

		b1 = tk.Button(var, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar44,self.e210),height = 1,width = 3)
		b1.grid(row=2, column=2)		

		self.l2 = tk.Label(var, text = '======= FIT FILE =======',justify='center')
		self.l2.grid(row=3,column=1)		

		self.l2 = tk.Label(var, text=" TE FILE [eV]")
		self.l2.grid(row=4, column=0)
		self.l2 = tk.Label(var, text=" NE FILE [(18)/m3]")
		self.l2.grid(row=5, column=0)
		self.l2 = tk.Label(var, text=" TI FILE [eV]")
		self.l2.grid(row=6, column=0)		

		self.e211 = tk.Entry(var,width=40,justify='center')
		self.e211.insert(0,self.StrVar45.get())
		self.e211.grid(row=4, column=1, padx=5)
		self.e212 = tk.Entry(var,width=40,justify='center')
		self.e212.insert(0,self.StrVar46.get())
		self.e212.grid(row=5, column=1, padx=5)
		self.e213 = tk.Entry(var,width=40,justify='center')
		self.e213.insert(0,self.StrVar47.get())
		self.e213.grid(row=6, column=1, padx=5)

		b2 = tk.Button(var, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar45,self.e211),height = 1,width = 3)
		b2.grid(row=4, column=2)	
		b3 = tk.Button(var, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar46,self.e212),height = 1,width = 3)
		b3.grid(row=5, column=2)			
		b4 = tk.Button(var, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar47,self.e213),height = 1,width = 3)
		b4.grid(row=6, column=2)	

		self.l2 = tk.Label(var, text = '======= RAW FILE =======',justify='center')
		self.l2.grid(row=8,column=1)		

		self.l2 = tk.Label(var, text=" TE FILE [eV]")
		self.l2.grid(row=9, column=0)
		self.l2 = tk.Label(var, text=" NE FILE [(18)/m3]")
		self.l2.grid(row=10, column=0)
		self.l2 = tk.Label(var, text=" TI FILE [eV]")
		self.l2.grid(row=11, column=0)
		self.l2 = tk.Label(var, text=" LINE_N(19)")
		self.l2.grid(row=13, column=0)	

		self.e215 = tk.Entry(var,width=40,justify='center')
		self.e215.insert(0,self.StrVar49.get())
		self.e215.grid(row=9, column=1, padx=5)
		self.e216 = tk.Entry(var,width=40,justify='center')
		self.e216.insert(0,self.StrVar50.get())
		self.e216.grid(row=10, column=1, padx=5)
		self.e217 = tk.Entry(var,width=40,justify='center')
		self.e217.insert(0,self.StrVar51.get())
		self.e217.grid(row=11, column=1, padx=5)
		self.e113 = tk.Entry(var,width=40,justify='center')
		self.e113.insert(0,self.StrVar53.get())
		self.e113.grid(row=13, column=1, padx=5)		

		b6 = tk.Button(var, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar49,self.e215),height = 1,width = 3)
		b6.grid(row=9, column=2)	
		b7 = tk.Button(var, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar50,self.e216),height = 1,width = 3)
		b7.grid(row=10, column=2)		
		b8 = tk.Button(var, text="OPEN", bg = "gray",command=lambda: self.button_func8(self.StrVar51,self.e217),height = 1,width = 3)
		b8.grid(row=11, column=2)	

		b10 = tk.Button(var, text="DONE", bg = "gray",command=lambda: self.button_func9(var),height = 1,width = 3)
		b10.grid(row=14, column=0,columnspan=3, pady=10)

		var.mainloop()

		return

	def gui_fit(self):

		self.fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2,2,figsize=(13.5,9))
		ax1.set_title('TE [keV]')
		ax2.set_title('NE [10(19)/m3]')
		ax3.set_title('TI [keV]')
		ax4.set_title('[a.u]')
		if not (self.fit.use_rho):
			ax1.set_xlabel('$\psi_n$ [a.u]')
			ax2.set_xlabel('$\psi_n$ [a.u]')
			ax3.set_xlabel('$\psi_n$ [a.u]')
			ax4.set_xlabel('$\psi_n$ [a.u]')
		else:
			ax1.set_xlabel('$\\rho_t$ [a.u]')
			ax2.set_xlabel('$\\rho_t$ [a.u]')
			ax3.set_xlabel('$\\rho_t$ [a.u]')
			ax4.set_xlabel('$\\rho_t$ [a.u]')

		ax1.set_ylabel('TE [keV]')
		ax2.set_ylabel('NE [10(19)/m3]')
		ax3.set_ylabel('TI [keV]')
		ax4.set_ylabel('[a.u]')

		self.fig.tight_layout()

		#self.fig.subplots_adjust(left=0.05, right=0.993, top=0.97, bottom=0.07)
		self.fit.read_kinprof_kprofile()
		self.fit.read_raw_kinprof()	
		#self.fit.main_run()
		self.set_axis(ax1,ax2,ax3,ax4)

		self.declare_vars()
		self.initialise_vars()			


		self.canvas = FigureCanvasTkAgg(self.fig,master=self.root)
		self.plot_widget = self.canvas.get_tk_widget()
		self.plot_widget.grid(rowspan=20,row=1,column=20)

		toolbar_frame = tk.Frame(self.root)
		toolbar_frame.grid(column=20,row=0)
		
		self.toolbar = NavigationToolbar2TkAgg(self.canvas,toolbar_frame)

		#Fitting option
		self.l1 = tk.Label(self.root, text="=============== Pedestal scan option ===============",justify='center')
		self.l1.grid(row=0, column=0,columnspan=9,padx=5,pady=15)
		 
		self.l2 = tk.Label(self.root, text="Use_TI_ped_width")
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar27)
		self.l2.grid(row=1, column=0, columnspan = 2)
		self.e2.grid(row=1, column=2)
		
		self.l1 = tk.Label(self.root, text="================== Fitting option ==================",justify='center')
		self.l1.grid(row=2, column=0,columnspan=9,padx=5)		

		self.l2 = tk.Label(self.root, text=" TE [keV]")
		self.l2.grid(row=4, column=0)
		self.l2 = tk.Label(self.root, text=" NE [(19)]")
		self.l2.grid(row=5, column=0)
		self.l2 = tk.Label(self.root, text=" TI [keV]")
		self.l2.grid(row=6, column=0)
		
		self.l2 = tk.Label(self.root, text="RAW")
		self.l2.grid(row=3, column=1)
		self.l2 = tk.Label(self.root, text="AVG")
		self.l2.grid(row=3, column=2)
		self.l2 = tk.Label(self.root, text="OUT")
		self.l2.grid(row=3, column=3)	
		self.l2 = tk.Label(self.root, text="SEP")
		self.l2.grid(row=3, column=4)	
		self.l2 = tk.Label(self.root, text="WIDTH")
		self.l2.grid(row=3, column=5)			
		self.l3 = tk.Label(self.root, text=" USE_STD ")
		self.l3.grid(row=3, column=6)				

		#Fitting RAW
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar15)
		self.e2.grid(row=4, column=1)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar16)
		self.e2.grid(row=5, column=1)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar17)
		self.e2.grid(row=6, column=1)

		#Fitting AVG
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar2)
		self.e2.grid(row=4, column=2)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar3)
		self.e2.grid(row=5, column=2)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar4)
		self.e2.grid(row=6, column=2)

		#Fitting Outlier
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar19)
		self.e2.grid(row=4, column=3)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar20)
		self.e2.grid(row=5, column=3)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar21)
		self.e2.grid(row=6, column=3)

		#Fitting SEP_CONST
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar6)
		self.e2.grid(row=4, column=4)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar7)
		self.e2.grid(row=5, column=4)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar8)
		self.e2.grid(row=6, column=4)

		#Fitting WIDTH_CONST
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar23)
		self.e2.grid(row=4, column=5)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar24)
		self.e2.grid(row=5, column=5)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar25)
		self.e2.grid(row=6, column=5)
			
		#Fitting STD
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar11)
		self.e2.grid(row=4, column=6)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar12)
		self.e2.grid(row=5, column=6)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar13)
		self.e2.grid(row=6, column=6)
	
		self.l9 = tk.Label(self.root, text = 'RAW_STD')
		self.l9.grid(row=3,column=7)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar64)
		self.e2.grid(row=4, column=7)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar65)
		self.e2.grid(row=5, column=7)
		self.e2 = tk.Checkbutton(self.root,variable=self.CheckVar66)
		self.e2.grid(row=6, column=7)

		self.l1 = tk.Label(self.root, text="================= Fitting variables =================",justify='center')
		self.l1.grid(row=7, column=0,columnspan=9,padx=5)

		b8 = tk.Button(self.root, text="PARAMS",bg = "lightsteelblue",command=lambda: self.button_func10())

		b8.grid(row=9, column=6,columnspan=2)
		
		self.l2 = tk.Label(self.root, text=" TE [keV]")
		self.l2.grid(row=9, column=0)
		self.l2 = tk.Label(self.root, text=" NE [(19)]")
		self.l2.grid(row=10, column=0)
		self.l2 = tk.Label(self.root, text=" TI [keV]")
		self.l2.grid(row=11, column=0)	
		
		self.l9 = tk.Label(self.root, text = 'Range')
		self.l9.grid(row=8,column=1)
		#Fitting range
		self.e10 = tk.Entry(self.root,width=6,justify='center')
		self.e10.insert(10,self.StrVar5.get())
		self.e10.grid(row=9, column=1, padx=5)
		self.e11 = tk.Entry(self.root,width=6,justify='center')
		self.e11.insert(10,self.StrVar6.get())
		self.e11.grid(row=10, column=1, padx=5)
		self.e12 = tk.Entry(self.root,width=6,justify='center')
		self.e12.insert(10,self.StrVar7.get())
		self.e12.grid(row=11, column=1, padx=5)

		self.l9 = tk.Label(self.root, text = 'SEP')
		self.l9.grid(row=8,column=2)
		#Separatrix value
		self.e16 = tk.Entry(self.root,width=6,justify='center')
		self.e16.insert(0,self.StrVar9.get())
		self.e16.grid(row=9, column=2, padx=5)
		self.e18 = tk.Entry(self.root,width=6,justify='center')
		self.e18.insert(0,self.StrVar10.get())
		self.e18.grid(row=10, column=2)
		self.e20 = tk.Entry(self.root,width=6,justify='center')
		self.e20.insert(0,self.StrVar11.get())
		self.e20.grid(row=11, column=2,padx=5, pady=5)


		self.l9 = tk.Label(self.root, text = 'WIDTH')
		self.l9.grid(row=8,column=3)
		#Width value	
		self.e25 = tk.Entry(self.root,width=6,justify='center')
		self.e25.insert(0,self.StrVar13.get())
		self.e25.grid(row=9, column=3, padx=5)
		self.e26 = tk.Entry(self.root,width=6,justify='center')
		self.e26.insert(0,self.StrVar14.get())
		self.e26.grid(row=10, column=3)
		self.e27 = tk.Entry(self.root,width=6,justify='center')
		self.e27.insert(0,self.StrVar15.get())
		self.e27.grid(row=11, column=3)

		self.l9 = tk.Label(self.root, text = 'OLC')
		self.l9.grid(row=8,column=4)
		#OLI_CUT	
		self.e39 = tk.Entry(self.root,width=6,justify='center')
		self.e39.insert(10,self.StrVar28.get())
		self.e39.grid(row=9, column=4, padx=5)
		self.e40 = tk.Entry(self.root,width=6,justify='center')
		self.e40.insert(10,self.StrVar29.get())
		self.e40.grid(row=10, column=4)
		self.e41 = tk.Entry(self.root,width=6,justify='center')
		self.e41.insert(10,self.StrVar30.get())
		self.e41.grid(row=11, column=4)

		self.l9 = tk.Label(self.root, text = 'OLCN')
		self.l9.grid(row=8,column=5)
		#OLI_CUT	
		self.e43 = tk.Entry(self.root,width=6,justify='center')
		self.e43.insert(10,self.StrVar32.get())
		self.e43.grid(row=9, column=5, padx=5)
		self.e44 = tk.Entry(self.root,width=6,justify='center')
		self.e44.insert(10,self.StrVar33.get())
		self.e44.grid(row=10, column=5)
		self.e45 = tk.Entry(self.root,width=6,justify='center')
		self.e45.insert(10,self.StrVar34.get())
		self.e45.grid(row=11, column=5)
			
		self.l29=tk.Label(self.root, text = '====== Removal data point ======',justify='center')
		self.l29.grid(row=12,column=1, columnspan=5)

		self.l30 = tk.Label(self.root, text="TE ",anchor='e')
		self.e30 = tk.Entry(self.root,width=30,justify='center')
		self.e30.insert(10,self.StrVar17.get())
		self.l30.grid(row=13, column=0)
		self.e30.grid(row=13, column=1,columnspan=5)

		self.l31 = tk.Label(self.root, text="NE",anchor='e')
		self.e31 = tk.Entry(self.root,width=30,justify='center')
		self.e31.insert(10,self.StrVar18.get())
		self.l31.grid(row=14, column=0)
		self.e31.grid(row=14, column=1,columnspan=5)

		self.l32 = tk.Label(self.root, text="TI",anchor='e')
		self.e32 = tk.Entry(self.root,width=30,justify='center')
		self.e32.insert(10,self.StrVar19.get())
		self.l32.grid(row=15, column=0)
		self.e32.grid(row=15, column=1,columnspan=5)

		self.l100 = tk.Label(self.root, text = '======== ETC ========',justify='center')
		self.l100.grid(row=16,column=1, columnspan=4)
			
		self.l101 = tk.Label(self.root, text="ZEFF ",anchor='e')
		self.e101 = tk.Entry(self.root,width=6,justify='center')
		self.e101.insert(10,self.StrVar25.get())
		self.l101.grid(row=17, column=0)
		self.e101.grid(row=17, column=1)

		self.l102 = tk.Label(self.root, text="ZIMP ",anchor='e')
		self.e102 = tk.Entry(self.root,width=6,justify='center')
		self.e102.insert(10,self.StrVar26.get())
		self.l102.grid(row=17, column=2)
		self.e102.grid(row=17, column=3)

		self.l103 = tk.Label(self.root, text="LINEN(19) ",anchor='e')
		self.e103 = tk.Entry(self.root,width=6,justify='center')
		self.e103.insert(10,self.StrVar53.get())
		self.l103.grid(row=17, column=4)
		self.e103.grid(row=17, column=5)

		self.l104 = tk.Label(self.root, text="AMAIN",anchor='e')
		self.e104 = tk.Entry(self.root,width=6,justify='center')
		self.e104.insert(10,self.StrVar27.get())
		self.l104.grid(row=18, column=0)
		self.e104.grid(row=18, column=1)

		self.l105 = tk.Label(self.root, text="AIMP ",anchor='e')
		self.e105 = tk.Entry(self.root,width=6,justify='center')
		self.e105.insert(10,self.StrVar54.get())
		self.l105.grid(row=18, column=2,pady=5)
		self.e105.grid(row=18, column=3)

		self.l39=tk.Label(self.root, text = '==== Done ====',justify='center')
		self.l39.grid(row=12, column=6,columnspan=4)

		b1 = tk.Button(self.root, text="Do Fitting", bg = "gold",command=lambda: self.button_func1(),height = 1,width = 15)
		b2 = tk.Button(self.root, text="Save Fitting",bg = "lime",command=lambda: self.button_func2(),height = 1,width = 15)
		b3 = tk.Button(self.root, text="Load Setting",bg = "cadetblue",command=lambda: self.button_func3(),height = 1,width = 15)
		b4 = tk.Button(self.root, text="End Fitting",bg = "pink",command=lambda: self.button_func4(),height = 1,width = 15)
		b5 = tk.Button(self.root, text="Reset Fitting",bg = "violet",command=lambda: self.button_func5(),height = 1,width = 15)
		#b6 = tk.Button(self.root, text="Load default",bg = "peru",command=lambda: self.button_func6(),height = 1,width = 15)
		b7 = tk.Button(self.root, text="Select files",bg = "peru",command=lambda: self.button_func7(),height = 1,width = 15)
		b11 = tk.Button(self.root, text="2ND PROF",bg = "lightsteelblue",command=lambda: self.button_func15(), height = 1,width = 8)

		b1.grid(row=13, column=6,columnspan=4)
		b2.grid(row=14, column=6,columnspan=4)
		b3.grid(row=16, column=6,columnspan=4)
		b4.grid(row=18, column=6,columnspan=4)
		b5.grid(row=15, column=6,columnspan=4)
		#b6.grid(row=19, column=6,columnspan=4)
		b7.grid(row=17, column=6,columnspan=4)
		b11.grid(row=10, column=6,columnspan=2)

		self.root.mainloop()
		return

	def initialise_vars(self):
	
		self.CheckVar1.set(self.trans_vars(self.fit.ped_scan_fit,4))
		
		self.CheckVar2.set(self.trans_vars(self.fit.te_use_avg_dat,4))
		self.CheckVar3.set(self.trans_vars(self.fit.ne_use_avg_dat,4))
		self.CheckVar4.set(self.trans_vars(self.fit.ti_use_avg_dat,4))
		self.CheckVar5.set(self.trans_vars(self.fit.vt_use_avg_dat,4))

		self.CheckVar15.set(self.trans_vars(self.fit.te_raw_fit,4))
		self.CheckVar16.set(self.trans_vars(self.fit.ne_raw_fit,4))
		self.CheckVar17.set(self.trans_vars(self.fit.ti_raw_fit,4))
		self.CheckVar18.set(self.trans_vars(self.fit.vt_raw_fit,4))

		self.CheckVar19.set(self.trans_vars(self.fit.te_use_oli,4))
		self.CheckVar20.set(self.trans_vars(self.fit.ne_use_oli,4))
		self.CheckVar21.set(self.trans_vars(self.fit.ti_use_oli,4))
		self.CheckVar22.set(self.trans_vars(self.fit.vt_use_oli,4))
		
		self.CheckVar6.set(self.trans_vars(self.fit.fixed_sep_te,4))
		self.CheckVar7.set(self.trans_vars(self.fit.fixed_sep_ne,4))
		self.CheckVar8.set(self.trans_vars(self.fit.fixed_sep_ti,4))
		self.CheckVar9.set(self.trans_vars(self.fit.fixed_sep_vt,4))
		self.CheckVar10.set(self.trans_vars(self.fit.use_ti_eped,4))
		self.CheckVar27.set(self.trans_vars(self.fit.use_ti_width,4))
		
		self.CheckVar23.set(self.trans_vars(self.fit.fixed_width_te,4))
		self.CheckVar24.set(self.trans_vars(self.fit.fixed_width_ne,4))
		self.CheckVar25.set(self.trans_vars(self.fit.fixed_width_ti,4))
		self.CheckVar26.set(self.trans_vars(self.fit.fixed_width_vt,4))
		
		self.CheckVar11.set(self.trans_vars(self.fit.std_te,4))
		self.CheckVar12.set(self.trans_vars(self.fit.std_ne,4))
		self.CheckVar13.set(self.trans_vars(self.fit.std_ti,4))
		self.CheckVar14.set(self.trans_vars(self.fit.std_vt,4))

		self.CheckVar63.set(self.trans_vars(self.fit.use_rho,4))

		self.CheckVar64.set(self.trans_vars(self.fit.use_raws_te,4))
		self.CheckVar65.set(self.trans_vars(self.fit.use_raws_ne,4))
		self.CheckVar66.set(self.trans_vars(self.fit.use_raws_ti,4))
		self.CheckVar67.set(self.trans_vars(self.fit.use_raws_vt,4))
		
		self.StrVar1.set(self.trans_vars(self.fit.fit_type_te,5))
		self.StrVar2.set(self.trans_vars(self.fit.fit_type_ne,5))
		self.StrVar3.set(self.trans_vars(self.fit.fit_type_ti,5))
		self.StrVar4.set(self.trans_vars(self.fit.fit_type_vt,5))

		self.StrVar5.set(self.trans_vars(self.fit.psi_end_te,2))
		self.StrVar6.set(self.trans_vars(self.fit.psi_end_ne,2))
		self.StrVar7.set(self.trans_vars(self.fit.psi_end_ti,2))
		self.StrVar8.set(self.trans_vars(self.fit.psi_end_vt,2))
		
		self.StrVar9.set(self.trans_vars(self.fit.te_sep,2))
		self.StrVar10.set(self.trans_vars(self.fit.ne_sep,2))
		self.StrVar11.set(self.trans_vars(self.fit.ti_sep,2))
		self.StrVar12.set(self.trans_vars(self.fit.vt_sep,2))
		
		self.StrVar13.set(self.trans_vars(self.fit.te_width,2))
		self.StrVar14.set(self.trans_vars(self.fit.ne_width,2))
		self.StrVar15.set(self.trans_vars(self.fit.ti_width,2))
		self.StrVar16.set(self.trans_vars(self.fit.vt_width,2))
		
		self.StrVar17.set(self.trans_vars(self.fit.te_excn,6))
		self.StrVar18.set(self.trans_vars(self.fit.ne_excn,6))
		self.StrVar19.set(self.trans_vars(self.fit.ti_excn,6))
		self.StrVar20.set(self.trans_vars(self.fit.vt_excn,6))
		
		self.StrVar21.set(self.trans_vars(self.fit.s_te,2))
		self.StrVar22.set(self.trans_vars(self.fit.s_ne,2))
		self.StrVar23.set(self.trans_vars(self.fit.s_ti,2))
		self.StrVar24.set(self.trans_vars(self.fit.s_vt,2))

		self.StrVar25.set(self.trans_vars(self.fit.zeff,2))
		self.StrVar26.set(self.trans_vars(self.fit.zimp,2))
		self.StrVar27.set(self.trans_vars(self.fit.amain,2))
		self.StrVar54.set(self.trans_vars(self.fit.aimp,2))

		self.StrVar28.set(self.trans_vars(self.fit.te_oli_cut,2))
		self.StrVar29.set(self.trans_vars(self.fit.ne_oli_cut,2))
		self.StrVar30.set(self.trans_vars(self.fit.ti_oli_cut,2))
		self.StrVar31.set(self.trans_vars(self.fit.vt_oli_cut,2))

		self.StrVar32.set(self.trans_vars(self.fit.te_oli_cutn,1))
		self.StrVar33.set(self.trans_vars(self.fit.ne_oli_cutn,1))
		self.StrVar34.set(self.trans_vars(self.fit.ti_oli_cutn,1))
		self.StrVar35.set(self.trans_vars(self.fit.vt_oli_cutn,1))

		self.StrVar36.set(self.trans_vars(self.fit.te_height,2))
		self.StrVar37.set(self.trans_vars(self.fit.ne_height,2))
		self.StrVar38.set(self.trans_vars(self.fit.ti_height,2))
		self.StrVar39.set(self.trans_vars(self.fit.vt_height,2))

		self.StrVar40.set(self.trans_vars(self.fit.te_core,2))
		self.StrVar41.set(self.trans_vars(self.fit.ne_core,2))
		self.StrVar42.set(self.trans_vars(self.fit.ti_core,2))
		self.StrVar43.set(self.trans_vars(self.fit.vt_core,2))

		self.StrVar44.set(self.trans_vars(self.fit.eqdsk_name,3))

		self.StrVar45.set(self.trans_vars(self.fit.te_file,3))
		self.StrVar46.set(self.trans_vars(self.fit.ne_file,3))
		self.StrVar47.set(self.trans_vars(self.fit.ti_file,3))
		self.StrVar48.set(self.trans_vars(self.fit.vt_file,3))

		self.StrVar49.set(self.trans_vars(self.fit.te_dat_file,3))
		self.StrVar50.set(self.trans_vars(self.fit.ne_dat_file,3))
		self.StrVar51.set(self.trans_vars(self.fit.ti_dat_file,3))
		self.StrVar52.set(self.trans_vars(self.fit.vt_dat_file,3))

		self.StrVar53.set(self.trans_vars(self.fit.target_density,2))

		self.initialise_lmfit_param()

		self.t3_close = True
		self.t2_close = True
		self.t_close = True
		self.did_fit = False

		return

	def initialise_lmfit_param(self):

		count = 101
		for i in ['te','ne','ti']:
			for j in ['val','min','max']:
				for k in range(6):
					self.__dict__['StrVar%d'%count].set(self.trans_vars(self.fit.param[i][j][k],2))
					count += 1

		count = 31
		for i in ['te','ne','ti']:
				for j in range(6):
					self.__dict__['CheckVar%d'%count].set(self.trans_vars((self.fit.param[i]['vary'][j]==False),4))
					count += 1	

		return

	def set_lmfit_param(self):

		count = 101
		for i in ['te','ne','ti']:
			for j in ['val','min','max']:
				for k in range(6):
					self.fit.param[i][j][k] = self.trans_vars2(self.__dict__['StrVar%d'%count].get(),2)
					count += 1

		count = 31
		for i in ['te','ne','ti']:
				for j in range(6):
					self.fit.param[i]['vary'][j] = (self.trans_vars2(self.__dict__['CheckVar%d'%count].get(),4)==False)
					count += 1

		return						

	def update_lmfit_param(self,copyf=True):

		for i in range(101,155):
			self.renew_box(self.__dict__['e%d'%(i+200)],self.__dict__['StrVar%d'%(i)],copyf)

		return

	def make_same_val(self,var1,var2,e1,e2,typev=1):

		var2.set(e2.get())
		if not (typev == 1):
			if (typev == 2 ):	
				temp = abs(float(var2.get()))
			elif (typev == 3 ):	
				temp = -abs(float(var2.get()))

			temp2 = str(temp)
			var1.set(temp2)
		else:
			var1.set(var2.get())

		e1.delete(0, 'end')
		e1.insert(0,var1.get())

		return

	def interact_main_lmfit(self,var_ind,var_ind1,var_ind2,var_ind3,e_ind1,e_ind2):

		
		var1 = self.__dict__['StrVar%d'%var_ind1]
		var2 = self.__dict__['StrVar%d'%var_ind2]
		var3 = self.__dict__['StrVar%d'%(var_ind1+6)]
		var4 = self.__dict__['StrVar%d'%(var_ind1+12)]
	
		e1 = self.__dict__['e%d'%e_ind1]
		e2 = self.__dict__['e%d'%e_ind2]
		e3 = self.__dict__['e%d'%(e_ind1+6)]
		e4 = self.__dict__['e%d'%(e_ind1+12)]

		if not (var_ind==None):
			var  = self.__dict__['CheckVar%d'%var_ind]
			var5 = self.__dict__['CheckVar%d'%var_ind3]
			if (var.get() == 1):
				self.make_same_val(var1,var2,e1,e2)		
				var5.set(var.get())
			else:
				var5.set(var.get())
				if (float(e2.get()) < 0.0):
					self.make_same_val(var3,var2,e3,e2,2)
				if (float(e2.get()) > 0.0):
					self.make_same_val(var4,var2,e4,e2)
		else:
				if (float(e2.get()) < 0.0):
					self.make_same_val(var3,var2,e3,e2,2)
				if (float(e2.get()) > 0.0):
					self.make_same_val(var4,var2,e4,e2)
		return

	def interact_main_to_lmfit(self):
		#sep
		self.interact_main_lmfit(6,101,9,31,301,16) 	#te
		self.interact_main_lmfit(7,119,10,37,319,18)	#ne
		self.interact_main_lmfit(8,137,11,43,337,20)	#ti

		#width
		self.interact_main_lmfit(23,103,13,33,303,25) 		#te
		self.interact_main_lmfit(24,121,14,39,321,26)		#ne
		self.interact_main_lmfit(25,139,15,45,339,27)		#ti
		return

	def interact_lmfit_main(self,var_ind,var_ind1,var_ind2,var_ind3,e_ind1,e_ind2):

		var1 = self.__dict__['StrVar%d'%var_ind1]
		var2 = self.__dict__['StrVar%d'%var_ind2]
		var3 = self.__dict__['StrVar%d'%(var_ind1+6)]
		var4 = self.__dict__['StrVar%d'%(var_ind1+12)]
	
		e1 = self.__dict__['e%d'%e_ind1]
		e2 = self.__dict__['e%d'%e_ind2]
		e3 = self.__dict__['e%d'%(e_ind1+6)]
		e4 = self.__dict__['e%d'%(e_ind1+12)]

		if not (var_ind==None):
			var  = self.__dict__['CheckVar%d'%var_ind]
			var5 = self.__dict__['CheckVar%d'%var_ind3]
			if (var5.get() == 1):
				self.make_same_val(var2,var1,e2,e1)		
				var.set(var5.get())
			else:
				var.set(var5.get())
				if (float(e2.get()) < 0.0):
					self.make_same_val(var2,var3,e2,e3,3)
				if (float(e2.get()) > 0.0):
					self.make_same_val(var2,var4,e2,e4)
		else:
				if (float(e2.get()) < 0.0):
					self.make_same_val(var2,var3,e2,e3,3)
				if (float(e2.get()) > 0.0):
					self.make_same_val(var2,var4,e2,e4)
		return

	def interact_lmfit_to_main(self):
		#sep
		self.interact_lmfit_main(6,101,9,31,301,16) 	#te
		self.interact_lmfit_main(7,119,10,37,319,18)	#ne
		self.interact_lmfit_main(8,137,11,43,337,20)	#ti

		#width
		self.interact_lmfit_main(23,103,13,33,303,25) 		#te
		self.interact_lmfit_main(24,121,14,39,321,26)		#ne
		self.interact_lmfit_main(25,139,15,45,339,27)		#ti
		return

	def write_fit_opt(self,filename='fit_opt',wtype=False):
	
		f = open(filename,'w')
		f.write('!--- EQDSK option \n')

		if not (wtype):

			f.write('EQDSK = %s \n'%self.trans_vars2(self.StrVar44.get(),3))

			if (self.fit.use_chease):
				f.write('USE_CHEASE = True \n')
			else:
				f.write('USE_CHEASE = False \n')

			f.write('USE_RHO = %s \n'%self.trans_vars2(self.CheckVar63.get(),4))
			f.write(' \n')
			f.write('!--- Fitted profile \n')
			f.write('NE_file = %s  \n'%self.trans_vars2(self.StrVar46.get(),3))
			f.write('TE_file = %s  \n'%self.trans_vars2(self.StrVar45.get(),3))		
			f.write('TI_file = %s  \n'%self.trans_vars2(self.StrVar47.get(),3))		
			f.write(' \n')
			f.write('!--- Raw data  \n')
			f.write('NE_dat_file = %s  \n'%self.trans_vars2(self.StrVar50.get(),3))		
			f.write('TE_dat_file = %s  \n'%self.trans_vars2(self.StrVar49.get(),3))		
			f.write('TI_dat_file = %s  \n'%self.trans_vars2(self.StrVar51.get(),3))			

			f.write('Density_scale = %f \n'%self.trans_vars2(self.e103.get(),2))

		else:	

			try:
				f.write('EQDSK = %s \n'%self.trans_vars2(self.e210.get(),3))
			except:
				f.write('EQDSK = %s \n'%self.fit.eqdsk_name)
			try:
				if (self.fit.use_chease):
					f.write('USE_CHEASE = True \n')
				else:
					f.write('USE_CHEASE = False \n')
			except:
				f.write('USE_CHEASE = False \n')
			
			f.write('USE_RHO = %s \n'%self.trans_vars2(self.CheckVar63.get(),4))
			f.write(' \n')
			f.write('!--- Fitted profile \n')

			try:
				f.write('NE_file = %s  \n'%self.trans_vars2(self.e212.get(),3))
			except:
				f.write('NE_file = %s  \n'%self.fit.ne_file)
			try:
				f.write('TE_file = %s  \n'%self.trans_vars2(self.e211.get(),3))		
			except:
				f.write('TE_file = %s  \n'%self.fit.te_file)
			try:
				f.write('TI_file = %s  \n'%self.trans_vars2(self.e213.get(),3))		
			except:
				f.write('TI_file = %s  \n'%self.fit.ti_file)

			f.write(' \n')
			f.write('!--- Raw data  \n')

			try:
				f.write('NE_dat_file = %s  \n'%self.trans_vars2(self.e216.get(),3))		
			except:
				f.write('NE_dat_file = %s  \n'%self.fit.ne_dat_file)	
			try:
				f.write('TE_dat_file = %s  \n'%self.trans_vars2(self.e215.get(),3))		
			except:
				f.write('TE_dat_file = %s  \n'%self.fit.te_dat_file)
			try:
				f.write('TI_dat_file = %s  \n'%self.trans_vars2(self.e217.get(),3))		
			except:
				f.write('TI_dat_file = %s  \n'%self.fit.ti_dat_file)

			if not (self.e113.get() == ''):
				f.write('Density_scale = %f \n'%self.trans_vars2(self.e113.get(),2))
			else:
				f.write('Density_scale = %f \n'%1.0)

		if (wtype):
			f.close()
			return
		
		f.write(' \n')
		f.write('!--- Fitting option  \n')
		f.write('USE_TI_WIDTH = %s \n'%self.trans_vars2(self.CheckVar27.get(),4))
		
		f.write('TE_AVG_DAT = %s \n'%self.trans_vars2(self.CheckVar2.get(),4))
		f.write('NE_AVG_DAT = %s \n'%self.trans_vars2(self.CheckVar3.get(),4))
		f.write('TI_AVG_DAT = %s \n'%self.trans_vars2(self.CheckVar4.get(),4))

		
		f.write('TE_RAW_FIT = %s \n'%self.trans_vars2(self.CheckVar15.get(),4))
		f.write('NE_RAW_FIT = %s \n'%self.trans_vars2(self.CheckVar16.get(),4))
		f.write('TI_RAW_FIT = %s \n'%self.trans_vars2(self.CheckVar17.get(),4))

		
		f.write('TE_USE_OUT = %s \n'%self.trans_vars2(self.CheckVar19.get(),4))
		f.write('NE_USE_OUT = %s \n'%self.trans_vars2(self.CheckVar20.get(),4))
		f.write('TI_USE_OUT = %s \n'%self.trans_vars2(self.CheckVar21.get(),4))


		f.write('TE_OLI_CUT = %f \n'%self.trans_vars2(self.e39.get(),2))
		f.write('NE_OLI_CUT = %f \n'%self.trans_vars2(self.e40.get(),2))
		f.write('TI_OLI_CUT = %f \n'%self.trans_vars2(self.e41.get(),2))


		f.write('TE_OLI_CUTN = %i \n'%self.trans_vars2(self.e43.get(),2))
		f.write('NE_OLI_CUTN = %i \n'%self.trans_vars2(self.e44.get(),2))
		f.write('TI_OLI_CUTN = %i \n'%self.trans_vars2(self.e45.get(),2))
	

		f.write(' \n')
		f.write('!--- Fitting fun type \n')
		f.write('FIT_type_TE = %i \n'%self.trans_vars2(self.StrVar1.get(),5))
		f.write('FIT_type_NE = %i \n'%self.trans_vars2(self.StrVar2.get(),5))
		f.write('FIT_type_TI = %i \n'%self.trans_vars2(self.StrVar3.get(),5))
		f.write(' \n')
		f.write('!--- Fitting range \n')
		f.write('psi_end_te = %f \n'%self.trans_vars2(self.e10.get(),2))
		f.write('psi_end_ne = %f \n'%self.trans_vars2(self.e11.get(),2))
		f.write('psi_end_ti = %f \n'%self.trans_vars2(self.e12.get(),2))
		f.write(' \n')
		f.write('!--- Sep value const. \n')
		f.write('fix_te_sep = %s \n'%self.trans_vars2(self.CheckVar6.get(),4))
		f.write('fix_ne_sep = %s \n'%self.trans_vars2(self.CheckVar7.get(),4))
		f.write('fix_ti_sep = %s \n'%self.trans_vars2(self.CheckVar8.get(),4))
		f.write('ne_sep = %f \n'%self.trans_vars2(self.e18.get(),2))
		f.write('te_sep = %f \n'%self.trans_vars2(self.e16.get(),2))
		f.write('ti_sep = %f  \n'%self.trans_vars2(self.e20.get(),2))
		f.write(' \n')
		f.write('!--- Pedestal width const. \n')
		f.write('FIX_TE_width = %s \n'%self.trans_vars2(self.CheckVar23.get(),4))
		f.write('FIX_NE_width = %s \n'%self.trans_vars2(self.CheckVar24.get(),4))
		f.write('FIX_TI_width = %s \n'%self.trans_vars2(self.CheckVar25.get(),4))
		
		f.write('TE_width = %f \n'%self.trans_vars2(self.e25.get(),2))
		f.write('NE_width = %f \n'%self.trans_vars2(self.e26.get(),2))
		f.write('TI_width = %f \n'%self.trans_vars2(self.e27.get(),2))
		f.write(' \n')

		f.write('!--- Removal point \n')
		f.write('ne_exc = %s \n'%self.trans_vars2(self.e31.get(),3))
		f.write('te_exc = %s \n'%self.trans_vars2(self.e30.get(),3))
		f.write('ti_exc = %s \n'%self.trans_vars2(self.e32.get(),3))
		f.write(' \n')
		f.write('!--- Smooth spline weight \n')
		f.write('STD_TE = %s \n'%self.trans_vars2(self.CheckVar11.get(),4))
		f.write('STD_NE = %s \n'%self.trans_vars2(self.CheckVar12.get(),4))
		f.write('STD_TI = %s \n'%self.trans_vars2(self.CheckVar13.get(),4))

		f.write('RAW_STD_TE = %s \n'%self.trans_vars2(self.CheckVar64.get(),4))
		f.write('RAW_STD_NE = %s \n'%self.trans_vars2(self.CheckVar65.get(),4))
		f.write('RAW_STD_TI = %s \n'%self.trans_vars2(self.CheckVar66.get(),4))

		f.write('\n')
		f.write('!--- Plasma params \n')
		f.write('AMAIN = %f \n'%self.fit.amain)
		f.write('ZEFF = %f \n'%self.trans_vars2(self.e101.get(),2))
		f.write('AIMP = %f \n'%self.trans_vars2(self.e105.get(),2))
		f.write('ZIMP = %f \n'%self.trans_vars2(self.e102.get(),2))
		f.write('AMAIN = %f \n'%self.trans_vars2(self.e104.get(),2))
	
		f.write('\n')
		
		f.close()
		
		return
	
	def set_axis(self, ax1, ax2, ax3, ax4):
		
		if not (self.fit.te_dat_file == None):
			ax1.axis([-0.05,1.25,min(-0.2,1.1*min(self.fit.te_datx)),1.1*max(self.fit.te_datx)])
		elif not (self.fit.te_file == None):
			ax1.axis([-0.05,1.25,min(-0.2,1.1*min(self.fit.tek)),1.1*max(self.fit.tek)])
		else:
			ax1.axis([-0.05,1.25,0.,1])

		if not (self.fit.ne_dat_file == None):
			ax2.axis([-0.05,1.25,min(-0.2,1.1*min(self.fit.ne_datx)),1.1*max(self.fit.ne_datx)])
		elif not (self.fit.ne_file == None):
			ax2.axis([-0.05,1.25,min(-0.2,1.1*min(self.fit.nek)),1.1*max(self.fit.nek)])
		else:
			ax2.axis([-0.05,1.25,0.,1])
				
		if not(self.fit.ti_dat_file == None):
			ax3.axis([-0.05,1.25,min(-0.2,1.1*min(self.fit.ti_datx)),1.1*max(self.fit.ti_datx)])
		elif not (self.fit.ti_file == None):
			ax3.axis([-0.05,1.25,min(-0.2,1.1*min(self.fit.tik)),1.1*max(self.fit.tik)])
		else:
			ax3.axis([-0.05,1.25,0.,1])
				
		return
		
	def declare_vars(self):
		#Vars
		for i in range(1,70):
			self.__dict__['CheckVar%d'%i] = tk.IntVar()

		for i in range(1,201):
			self.__dict__['StrVar%d'%i] = tk.StringVar()

		return

	def __init__(self):

		self.lmfit_mod = False

		return
		
if __name__ == "__main__":

	import eped_gfit

	#print(' ---------------------------------------------------------------')
	#print('||               Function based Fitting tool Ver 2.0           ||')
	#print('||                          by S.K.Kim                         ||')
	#print('||                   Kinetic Profile generator                 ||')
	#print(' ---------------------------------------------------------------\n')

	try:
		os.mkdir('PROFILES')
	except:
		pass
	

	guifit = eped_gfit.eped_gtool()	
	
	count = 0		

	if not (os.path.isfile('fit_opt')):
		guifit.root = tk.Tk()
		guifit.declare_vars()
		guifit.gui_fitopt(guifit.root)

	if not (os.path.isfile('fit_opt~')):
		copyfile('fit_opt','fit_opt~')

	guifit.fit = eped_fit.eped_ftool('chease_opt')				
	guifit.var_list = ['CORE','MTANH','PTANH','EPED','SPLINE','EPED2']

	guifit.reopen = True
	guifit.file_fixed = False
#	guifit.input = 'fit_opt_save'
	guifit.input = 'fit_opt'

	guifit.root = tk.Tk()
	guifit.root.title('GFIT tool')
	guifit.fit.read_namelist(guifit.input)
	guifit.gui_fit()

	while guifit.reopen:

		print(">>> Load input file from '%s'"%guifit.input)

		try:
			guifit.fit.read_namelist(guifit.input)

			if (guifit.file_fixed):
				guifit.fit.te_file = guifit.te_file
				guifit.fit.ne_file = guifit.ne_file
				guifit.fit.ti_file = guifit.ti_file
				guifit.fit.vt_file = guifit.vt_file

				guifit.fit.te_dat_file = guifit.te_dat_file
				guifit.fit.ne_dat_file = guifit.ne_dat_file
				guifit.fit.ti_dat_file = guifit.ti_dat_file
				guifit.fit.vt_dat_file = guifit.vt_dat_file

				guifit.fit.target_density = guifit.fit.target_density
				guifit.fit.eqdsk_name = guifit.eqdsk_name

				guifit.file_fixed = False

			guifit.gui_fit()
		except:
			if (count == 0):
				print('>>> No saved fitting configuration')
				guifit.input = 'fit_opt~'
				count = count + 1
			else:
				guifit.reopen = False

