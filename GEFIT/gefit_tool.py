import os,sys
import numpy as np
import matplotlib.pyplot as plt
from shutil import move,copytree,rmtree,copyfile
from scipy.interpolate import interp1d,interp2d
import eqdsk
import knots_tool3
from batch_run import *
from get_efit import *
from copy import deepcopy
from exec_dirs import kindata_dir,mds_dir, mse_good_ch
from gefit_mse import mse_rad_err
import ch_tool
import pickle

def find_run_index(run_index, run_index2,use_last_index=True):

	i = 0;
	while True:

		dirs = os.getcwd() + '/%i'%i
		if os.path.isdir(dirs):	run_index.append(str(i))
		else:	break
		i = i + 1

	if run_index[-1] == '--':	run_index2 = ['--']
	else:
		len2 = len(run_index)
		run_index2 = []
		for i in range(len2-1):
			run_index2.append(run_index[i])
		if len(run_index2) == 0:
			run_index2 = ['--']
		
	dirs = os.getcwd() +'/0/last_index'
	if os.path.isfile(dirs):
		f = open(dirs,'r')
		last_index = f.readline().split()[0]
		f.close()
		if use_last_index: print('>>> Last run index -> %s'%last_index)
	else:	last_index = run_index[-1]
		
	return	(run_index, run_index2, last_index)

def draw_psirz(iseqdsk,fig,eq,sim):

	if not iseqdsk:
		return
	fig.canvas.draw_idle()
	[ax1] = fig.axes
	ax1.cla()
	ax1.contourf(eq.RR,eq.ZZ,eq.psirz,50)
	ax1.plot(eq.rzbdy[:,0],eq.rzbdy[:,1],'r')
	ax1.plot(eq.rzlim[:,0],eq.rzlim[:,1],'black')
	ax1.axis('scaled')
	ax1.set_title('Equilibrium')
	ax1.set_xlabel('R [m]')
	ax1.set_ylabel('Z [m]')
	ax1.set_xlim([min(eq.RR[0,:]),max(eq.RR[0,:])])
	ax1.set_ylim([min(eq.ZZ[:,0]),max(eq.ZZ[:,0])])

	save_dir2 = sim.currdir + '/%s/CSOLVE/eq.info'%(sim.MenuVar1.get())
	f = open(save_dir2,'w')
	f.write('%f\t%f'%(eq.rmag,eq.wmhd))
	f.close()

	return

def update_psirz(iseqdsk,fig,eq,rzbdy,sim):

	draw_psirz(iseqdsk,fig,eq,sim)

	if not iseqdsk:
		fig.canvas.draw_idle()
	[ax1] = fig.axes
	ax1.plot(rzbdy[:,0],rzbdy[:,1],'--',color='pink')
	ax1.axis('scaled')

	R_max = (max(rzbdy[:,0])+min(rzbdy[:,0]))*0.5
	Z_max = 0.
	R_min = (max(rzbdy[:,0])+min(rzbdy[:,0]))*0.5
	Z_min = 0.
	ind = np.zeros(4)

	for i in range(len(rzbdy[:,0])):

		if (rzbdy[i,0] > R_max):
			R_max = rzbdy[i,0]
			ind[0] = i
		if (rzbdy[i,0] < R_min):
			R_min = rzbdy[i,0]
			ind[1] = i
		if (rzbdy[i,1] > Z_max):
			Z_max = rzbdy[i,1]
			ind[2] = i
		if (rzbdy[i,1] < Z_min):
			Z_min = rzbdy[i,1]
			ind[3] = i

	r0 = (R_min + R_max) * 0.5
	a0 = (R_max - R_min) * 0.5
	z0 = (Z_min + Z_max) * 0.5

	tri_u  = (r0 - rzbdy[int(ind[2]),0])/a0
	tri_l  = (r0 - rzbdy[int(ind[3]),0])/a0
	elong  = (Z_max - Z_min) / a0 / 2.

	inform = 'tri_u = %4.3f\n'%tri_u
	inform = inform + 'tri_l = %4.3f\n'%tri_l
	inform = inform + 'elong = %4.3f\n'%elong

	ax1.text(eq.rmag,eq.zmag,inform,color='lime')

	plt.draw()	

	return	

def draw_eqdsk(ax1,ax2,eq):

	ax1.plot(eq.psin,-eq.pp/max(abs(eq.pp)),eq.psin,eq.ffp/max(abs(eq.ffp)))
	ax1.set_title("P'($\psi$), FF'($\psi$)")
	ax1.set_xlabel('Normalized radius ($\psi_N$)')
	ax1.set_ylabel('Normalized value [a.u]')
	ax1.legend(["P'($\psi$)","FF'($\psi$)"])
	ax1.set_xlim((-0.05,1.05))

	jplowf = np.copy(eq.psin)
	jphighf = np.copy(eq.psin)

	Rlow = interp1d(eq.prhoR[:,0],eq.prhoR[:,2],'slinear')
	Rhigh  = interp1d(eq.prhoR[:,0],eq.prhoR[:,3],'slinear')

	mu0 = np.pi * 4.0 * 1.e-7

	for i in range(len(eq.psin)):
		jplowf[i] = - Rlow(eq.psin[i]) * eq.pp[i]  - eq.ffp[i]/mu0/Rlow(eq.psin[i])
		jphighf[i] = -Rhigh(eq.psin[i]) * eq.pp[i] - eq.ffp[i]/mu0/Rhigh(eq.psin[i])

	ax2.plot(eq.psin,jplowf/1.e6,eq.psin,jphighf/1.e6)
	ax2.set_title("$J_{\phi}$($\psi$)")
	ax2.set_xlabel('Normalized radius ($\psi_N$)')
	ax2.set_ylabel('Current density [$MA/m^2$]')
	ax2.legend(["LFS","HFS","AVERAGE"])
	ax2.set_xlim((-0.05,1.05))

	return

def draw_kin_file(dirs,fig,ax1,ax2,ax3,eq=None,skip=False):

	f = open(dirs,'r')
	num = int(float(f.readline()))
	result = np.zeros((num,4))
	line = f.readline().split()
	zeff = float(line[0])
	zimp = float(line[1])
	amain = float(line[2])
	aimp = float(line[3])

	for i in range(num):
		line = f.readline().split()
		result[i,0] = float(line[0])
		result[i,1] = float(line[1])
		result[i,2] = float(line[2])
		result[i,3] = float(line[3])
	f.close()

	linen = 0;	wkin = 0;
	if not eq==None:
		linen,wkin = pos_kin_prof(result[:,0],result[:,2],result[:,1],result[:,3],zeff,zimp,eq)		

	if skip:	return (zeff,zimp,amain,zimp,linen,wkin)

	ax1.plot(result[:,0],result[:,2])
	ax1.set_title("NE[10(19)/m3]")
	ax1.set_xlabel('Normalized radius ($\psi_N$)')
	ax1.set_ylabel('NE[10(19)/m3]')
	ax1.set_xlim((-0.05,1.05))

	ax2.plot(result[:,0],result[:,1])
	ax2.set_title("TE[keV]")
	ax2.set_xlabel('Normalized radius ($\psi_N$)')
	ax2.set_ylabel('TE[keV]')
	ax2.set_xlim((-0.05,1.05))	

	ax3.plot(result[:,0],result[:,3])
	ax3.set_title("TI[keV]")
	ax3.set_xlabel('Normalized radius ($\psi_N$)')
	ax3.set_ylabel('TI[keV]')
	ax3.set_xlim((-0.05,1.05))				

	fig.tight_layout()
	
	return (zeff,zimp,amain,zimp,linen,wkin)

def draw_rot_file(dirs,fig,ax1,skip=False):
	if skip:	return
	f = open(dirs,'r')
	result = np.zeros((101,2))
	for i in range(2):	line = f.readline()

	for i in range(101):
		line = f.readline().split()
		result[i,0] = float(line[1])
		result[i,1] = float(line[2])

	ax1.plot(result[:,0],result[:,1])
	ax1.set_title("VT[km/s]")
	ax1.set_xlabel('Normalized radius ($\psi_N$)')
	ax1.set_ylabel('VT[km/s]')
	ax1.set_xlim((-0.05,1.05))
	fig.tight_layout()
	return	

def pos_kin_prof(psin,ne,te,ti,zeff,zimp,eq):

	linen = 0.;	wkin = 0.;
	pres = ne*(te + (1.0 - (zeff-1.)/zimp)*ti)
	pres = pres * 1.602*1.e-19*1.e19*1.e3
	Vf = interp1d(eq.avolp[:,0],eq.avolp[:,2],'cubic')
	vol = Vf(psin)
	nef = interp1d(psin,ne,'cubic')
	RR = np.zeros(201)
	NN = np.zeros(201)

	psin2 = np.linspace(0,2.0,301)
	R2 =np.copy(psin2);  R3 =np.copy(psin2);
	num = len(eq.prhoR[:,0])

	dRdp = (eq.prhoR[-1,2]-eq.prhoR[num-2,2])/(eq.prhoR[-1,0]-eq.prhoR[num-2,0])
	dRdp2 = (eq.prhoR[-1,3]-eq.prhoR[num-2,3])/(eq.prhoR[-1,0]-eq.prhoR[num-2,0])

	Rf = interp1d(eq.prhoR[:,0],eq.prhoR[:,2],'cubic')
	Rf2 = interp1d(eq.prhoR[:,0],eq.prhoR[:,3],'cubic')

	for i in range(301):
		if psin2[i] < 1.0:
			R2[i] = Rf(psin2[i])
			R3[i] = Rf2(psin2[i])
		else:
			R2[i] = dRdp*(psin2[i]-1.) + eq.prhoR[-1,2]
			R3[i] = dRdp2*(psin2[i]-1.) + eq.prhoR[-1,3]

	Rf = interp1d(psin2,R2,'cubic')
	Rf2 = interp1d(psin2,R3 ,'cubic')

	for i in range(101):
		psin = float(i)/100.
		RR[100-i] = Rf2(psin)
		RR[100+i] = Rf(psin)
		NN[100-i] = nef(psin)
		NN[100+i] = NN[100-i]

	linen = np.trapz(NN,x=RR) / 1.9*2.0
	wkin = np.trapz(pres,x=vol) * 1.5

	return linen, wkin

def read_comment(sim,filename):

	try:	textfile = open(filename,'r')
	except:	
		sim.textpad.delete('1.0',sim.END)
		return
	contents = textfile.read()
	textfile.close()
	sim.textpad.delete(1.0,sim.END)
	sim.textpad.insert('1.0',contents)
	return

def save_comment(sim,filename):

	textfile = open(filename,'w')
	contents = sim.textpad.get('1.0',sim.END).split('\n')
	nline = len(contents)
	for i in range(nline-2):
		textfile.write(contents[i]+'\n')
	for i in range(2):
		if not len(contents[nline-2+i])==0: textfile.write(contents[i]+'\n')
	textfile.close()
	
	return

def modify_fit_opt(filename,type,eq=None,te=None,ne=None,te_edge=None,ne_edge=None,ti=None,vt=None,refl='',ece=''):

	f = open(filename,'rb')
	fit_opt = pickle.load(f)
	f.close()

	if not eq==None: fit_opt['file']['gfile']     = '../../../%s'%eq
	if not te==None: fit_opt['file']['te']['ts']  = '../../../%s'%te
	if not ne==None: fit_opt['file']['ne']['ts']  = '../../../%s'%ne
	if not te_edge==None: fit_opt['file']['te']['tse'] = '../../../%s'%te_edge
	if not ne_edge==None: fit_opt['file']['ne']['tse'] = '../../../%s'%ne_edge
	if not ti==None: fit_opt['file']['ti']['ces'] = '../../../%s'%ti
	if not vt==None: fit_opt['file']['vt']['ces'] = '../../../%s'%vt

	if not refl=='': fit_opt['file']['ne']['refl'] = '../../../%s'%refl
	if not ece =='': fit_opt['file']['te']['ece']  = '../../../%s'%ece

	f = open(filename,'wb')
	pickle.dump(fit_opt,f)
	f.close()
	return

def modify_gfitp_opt(sim,filename,te='',ne='',te_edge='',ne_edge='',ti='',vt='',refl='',ece=''):

	f = open(filename,'r')
	f2 = open(filename+'_temp','w')
	while True:
		line = f.readline()
		if not line: break

		if line.lower().find('eq_file')>-1:
			line2 = 'eq_file ../../%s\n'%sim.e3.get()
		elif (line.lower().find('te_file')>-1 and not te == ''):
			line2 = 'te_file ../../%s\n'%te
		elif (line.lower().find('ne_file')>-1 and not ne == ''):
			line2 = 'ne_file ../../%s\n'%ne
		elif (line.lower().find('te_edge_file')>-1 and not te_edge == ''):
			line2 = 'te_edge_file ../../%s\n'%te_edge
		elif (line.lower().find('ne_edge_file')>-1 and not ne_edge == ''):
			line2 = 'ne_edge_file ../../%s\n'%ne_edge
		elif (line.lower().find('ti_file')>-1 and not ti == ''):
			if not line.lower().find('ti_file2')>-1:
				line2 = 'ti_file ../../%s\n'%ti
		elif (line.lower().find('vt_file')>-1 and not vt == ''):
			line2 = 'vt_file ../../%s\n'%vt

		elif (line.lower().find('refl_file')>-1 and not refl == ''):
			line2 = 'refl_file ../../%s\n'%refl
		elif (line.lower().find('ece_file')>-1 and not ece == ''):
			line2 = 'ece_file ../../%s\n'%ece

		elif line.lower().find('zeff')>-1:
			line2 = 'zeff %s\n'%sim.e7.get()
		elif line.lower().find('zimp')>-1:
			line2 = 'zimp %s\n'%sim.e14.get()			
		elif line.lower().find('lineden')>-1:
			line2 = 'lineden %s\n'%sim.e8.get()
		elif line.lower().find('bs_model')>-1:
			line2 = 'bs_model %s\n'%sim.MenuVar3.get().lower()
		elif line.lower().find('bs_model')>-1:
			line2 = 'bsmulti %s\n'%sim.e13.get()		

		else:
			line2 = line

		f2.write(line2)	

	f.close()
	f2.close()
	move(filename+'_temp',filename)

	return

def read_gfit_opt(sim,filename):

	f = open(filename,'r')

	while True:
		line = f.readline()
		if not line:	break
		
		if line.lower().find('ne_file') > -1:
			sim.ne_file = line.split()[1].strip()
		elif line.lower().find('te_file') > -1:
			sim.te_file = line.split()[1].strip()
		elif (line.lower().find('ti_file') > -1 and line.lower().find('ti_file2') == -1):
			sim.ti_file = line.split()[1].strip()
		elif line.lower().find('vt_file') > -1:
			sim.vt_file = line.split()[1].strip()
		elif line.lower().find('eq_file') > -1:
			sim.eq_file = line.split()[1].strip()

	f.close()

	return

def make_nubeam_config(sim):

	power = []; energy = [];
	
	power.append(sim.StrVar151.get());	energy.append(sim.StrVar152.get());
	power.append(sim.StrVar153.get());	energy.append(sim.StrVar154.get());
	power.append(sim.StrVar155.get());	energy.append(sim.StrVar156.get());
	power.append(sim.StrVar157.get());	energy.append(sim.StrVar158.get());
	power.append(sim.StrVar159.get());	energy.append(sim.StrVar160.get());
	power.append(sim.StrVar161.get());	energy.append(sim.StrVar162.get());

	for i in range(6):
		if float(power[5-i]) > 0.:	break

	if i == 5:
		if float(power[0]) == 0:	nbeam = 0;
		else:	nbeam = 1;
	else:
		nbeam = 6-i

	line = str(float(power[0])*1.e6)
	line2 = energy[0]

	for i in range(5):
		line = line + ',%s'%str(float(power[i+1])*1.e6)
		line2 = line2 + ',%s'%energy[i+1]

	if (nbeam > 0 and nbeam <=3): nbeam = 3
	if (nbeam > 0 and nbeam > 3): nbeam = 6

	return (line,line2,nbeam)

def make_chease_opt(sim,filename='chease_opt',vt=None):

	f = open(filename,'w')

	f.write('!-- run type \n')
	if sim.CheckVar22.get() == 1:	
		if sim.nbeam > 0:	f.write('RUN_MODE = NUBEAM \n')
		else:	f.write('RUN_MODE = NORMAL \n')
	else:	f.write('RUN_MODE = NORMAL \n')
	f.write('!-- input type \n')
	f.write('EQDSK = ../../../%s\n'%sim.e3.get())
	f.write('kinetic_profile_type = 1 \n')
	f.write('chease_kinetic_file = ../../../%s \n'%sim.e4.get())
	f.write('ne_file =  \n')
	f.write('te_file =  \n')
	f.write('ti_file =  \n')
	if not vt==None: f.write('vt_file = ../../../%s \n'%vt)
	if sim.nbeam > 0:
		f.write('USE_EXT_P = True \n')
		f.write('USE_EXT_J = True \n')
	else:
		f.write('USE_EXT_P = False \n')
		f.write('USE_EXT_J = False \n')		

	f.write('ZEFF  = %s \n'%sim.e7.get())
	f.write('ZIMP  = %s \n'%sim.e14.get())
	f.write('WDIA = %s \n'%sim.e6.get())
	f.write('APF = 0.0 \n')
	f.write('VLOOP_MOD = %s \n'%sim.StrVar102.get())
	f.write('BND_PSIN = %f \n'%sim.DoubleVar1.get())
	f.write('Beta_criterion_type = 0 \n')
	f.write('!-- run option (EXT_VLOOP = 0.0 then automatic run)\n')
	f.write('EXT_VLOOP = %s \n'%sim.StrVar103.get())
	f.write('Current_ITERN = %s \n'%sim.StrVar105.get())
	f.write('NUBEAM_ITERN = %s \n'%sim.StrVar106.get())
	f.write('RELAX = %s \n'%sim.StrVar107.get())
	if (sim.MenuVar3.get().lower() == 'nhager' or sim.MenuVar3.get().lower() == 'neo'): 
		f.write('USE_NEO = True \n')
		f.write('USE_HAGER = False \n')
		f.write('USE_CHANG = False \n')
	elif sim.MenuVar3.get().lower() == 'hager':
		f.write('USE_NEO = False \n')
		f.write('USE_HAGER = True \n')
		f.write('USE_CHANG = False \n')
	elif sim.MenuVar3.get().lower() == 'csauter':
		f.write('USE_NEO = False \n')
		f.write('USE_HAGER = False \n')
		f.write('USE_CHANG = True \n')
	else:	
		f.write('USE_NEO = False \n')
		f.write('USE_HAGER = False \n')
		f.write('USE_CHANG = False \n')
	if sim.CheckVar21.get() == 1:	f.write('HAG_CORE_MOD=True \n')
	else:	f.write('HAG_CORE_MOD=False \n')
	f.write('HAG_CORE_MOD_PSIN = %s \n'%sim.StrVar104.get())
	f.write('Core_neo = %s \n'%sim.StrVar101.get())
	f.write('DENSITY_SCALE = %s \n'%sim.e8.get())
	f.write('BSMULTI = %s\n'%sim.e13.get())
	f.write('NS = %s\n'%sim.StrVar108.get())
	f.write('NT = %s\n'%sim.StrVar109.get())
	f.write('MAP_NS = %s\n'%sim.StrVar110.get())
	f.write('MAP_NT = %s\n'%sim.StrVar111.get())
	f.write('IP_CRIT = %s\n'%sim.StrVar112.get())
	f.write('BS_CRIT = %s\n'%sim.StrVar113.get())
	f.write('!-- Output option \n')
	f.write('EFIT_CONST = True \n')
	if sim.CheckVar2.get() == 1:	f.write('ADJUST_PROF = True\n')
	else:	f.write('ADJUST_PROF = False\n')

	f.close()
	return

def make_nubeam_opt(sim,filename='nubeam_opt',eqtype=1):

	f = open(filename,'w')
	f.write('ZEFF = %s\n'%sim.e7.get())
	f.write('NBEAM = %i\n'%sim.nbeam)
	f.write('BPOWER = %s\n'%sim.bpower)
	f.write('BENERGY = %s\n'%sim.benergy)
	f.write('DIFFUSIVITY = %f\n'%round(float(sim.e12.get()),2))
	f.write('SHOT = %s\n'%sim.e1.get())
	f.write('TIME = %s\n'%sim.e2.get())
	if eqtype == 1:	f.write('eqdsk = geqdsk\n')
	else:	f.write('eqdsk = %s\n'%sim.e3.get())
	f.write('NPROC = %s\n'%sim.StrVar118.get())
	f.write('RUN_STEP = %s\n'%sim.StrVar114.get())
	f.write('RUN_AVG = %s\n'%sim.StrVar115.get())
	f.write('RUN_DT = %s\n'%sim.StrVar116.get())
	f.write('AVG_DT = %s\n'%sim.StrVar117.get())
	if not sim.StrVar119.get() == '':	f.write('MFILE = %s\n'%sim.StrVar119.get())
	if not sim.StrVar120.get() == '':	f.write('SFILE = %s\n'%sim.StrVar120.get())
	if not sim.StrVar121.get() == '':	f.write('IFILE = %s\n'%sim.StrVar121.get())

	f.close()

	return

def find_bdiff_index(sim,rtype=1,reset=False):

	if sim.MenuVar1.get() == '--':	
		index =[]
		return index

	if rtype == 1:
		if sim.MenuVar5.get().lower() == 'smse':	run_dir = sim.currdir + '/%s/CSOLVE'%sim.MenuVar1.get()
		elif sim.MenuVar5.get().lower() == 'emse':	run_dir = sim.currdir + '/%s/NUBEAM'%sim.MenuVar1.get()
	else:
		if sim.MenuVar7.get().lower() == 'smse':	run_dir = sim.currdir + '/%s/CSOLVE'%sim.MenuVar1.get()
		elif sim.MenuVar7.get().lower() == 'emse':	run_dir = sim.currdir + '/%s/NUBEAM'%sim.MenuVar1.get()

	os.chdir(run_dir)
	index = []
	for i in range(1000):
		dirs = run_dir + '/Bdiff_%03i/nubeam_iter_result'%i
		if os.path.isfile(dirs):	
			index = np.append(index,float(i)/100)

	os.chdir(sim.currdir)

	if sim.MenuVar5.get().lower() == sim.MenuVar7.get().lower():	rtype = 2

	if rtype == 2:

		ind = int(float(sim.StrVar106.get()))
		sim.iter_index = np.linspace(1,ind,ind)
		if int(sim.MenuVar9.get()) == 0:	sim.MenuVar9.set(str(ind))

		menu = sim.m9["menu"]
		menu.delete(0, "end")
		init = sim.MenuVar9.get()
		for string in sim.iter_index:
			string = str(int(string))
			menu.add_command(label=string, command=lambda value=int(string): sim.MenuVar9.set(value))
		sim.MenuVar9.set(init)

		menu = sim.m6["menu"]
		menu.delete(0, "end")
		if len(index) == 0.:	
			sim.MenuVar6.set('')
			return index

		menu.delete(0, "end")
		init = sim.MenuVar6.get()
		for string in index:
			menu.add_command(label=string, command=lambda value=str(round(string,2)): sim.MenuVar6.set(value))
		try:	
			if(float(init) in index):	sim.MenuVar6.set(init)
			else:	sim.MenuVar6.set('')
		except:	sim.MenuVar6.set(init)

		if reset:
			sim.MenuVar9.set(sim.StrVar106.get())					
			sim.MenuVar6.set(str(round(index[-1],2)))
	return index

def run_smse(sim,bdiff,chease_exec):

	run_dir = sim.currdir + '/%s/CHEASE/RUN_%03i'%(sim.MenuVar1.get(),round(bdiff*100.,0))
	vt_dir = sim.e45.get()
	out_dir1 = run_dir + '/OUTPUT'
	out_dir2 = run_dir + '/PROFILES'
	save_dir = sim.currdir + '/%s/CSOLVE/Bdiff_%03i'%(sim.MenuVar1.get(),round(bdiff*100.,0))
	save_dir2 = sim.currdir + '/%s/CSOLVE/eq.info'%(sim.MenuVar1.get())
	try:	os.mkdir(run_dir)
	except:	pass
	os.chdir(run_dir)
	make_chease_opt(sim,'chease_opt',vt_dir)
	make_nubeam_opt(sim,'nubeam_opt')

	if sim.nbeam == 0.:	sim.StrVar106.set(1)

	run_n = int(float(sim.StrVar106.get()))

	filename = run_dir+'/chease_batch'
	log_e = run_dir+'/chease_batch.e'
	log_o = run_dir+'/chease_batch.o'
	command = 'cd ' + run_dir + '\n'
	command = command + chease_exec+' chease_opt nomap \n'
	command = command + 'mkdir ' + save_dir + '\n'

	command = command + 'mv ' + out_dir1+'/nubeam_iter_result '+ save_dir+'/nubeam_iter_result'+'\n'

	#if sim.nbeam > 0:
	command = command + 'rmdir ' + save_dir+'/PROFILES'+'\n'
	command = command + 'mv ' + out_dir2+' '+save_dir+'/PROFILES'+'\n'
	
	for i in range(run_n+1):
		command = command + 'mv ' +out_dir1+'/geqdsk_%i'%i+' '+save_dir+'/geqdsk_%i'%i+'\n'
		command = command + 'mv ' +out_dir1+'/chease_namelist_%i'%i+' '+save_dir+'/chease_namelist_%i'%i+'\n'
		command = command + 'mv ' +out_dir1+'/EXPEQ_%i'%i+' '+save_dir+'/EXPEQ_%i'%i+'\n'
		command = command + 'mv ' +out_dir1+'/EFIT_JCONST_%i'%i+' '+save_dir+'/EFIT_JCONST_%i'%i+'\n'
		command = command + 'mv ' +out_dir1+'/Vneo.dat_%i'%i+' '+save_dir+'/Vneo.dat_%i'%i+'\n'
		command = command + 'mv ' +out_dir1+'/curr_prof_%i'%i+' '+save_dir+'/curr_prof_%i'%i+'\n'
		command = command + 'mv ' +out_dir2+'/neo_coefs_%i'%i+' '+save_dir+'/neo_coefs_%i'%i+'\n'
		command = command + 'mv ' +out_dir2+'/nubeam_out0d_%i'%i+' '+save_dir+'/nubeam_out0d_%i'%i+'\n'
		command = command + 'mv ' +out_dir2+'/nubeam_out1d_%i'%i+' '+save_dir+'/nubeam_out1d_%i'%i+'\n'
		

	make_batch_script(filename,None,log_e,log_o,command,'CHEASE_SMSE%03i'%round(bdiff*100.,0))
	runid = submit_batch_script(filename)

	os.chdir(sim.currdir)
	return

def run_emse(sim,bdiff,nubeamrun_exec,python2_exec):

	if sim.nbeam == 0.:	
		print('>>> No beam...')
		return

	run_dir = sim.currdir + '/%s/CHEASE/RUN_%03i'%(sim.MenuVar1.get(),round(bdiff*100.,0))
	nubeam_dir = run_dir + '/NUBEAM'
	profile_dir = run_dir + '/PROFILES'
	vt_dir = sim.e45.get()
	out_dir1 = run_dir + '/OUTPUT'
	out_dir2 = run_dir + '/PROFILES'
	save_dir = sim.currdir + '/%s/NUBEAM/Bdiff_%03i'%(sim.MenuVar1.get(),round(bdiff*100.,0))
	try:	os.mkdir(run_dir)
	except:	pass	
	os.chdir(run_dir)
	make_chease_opt(sim,'chease_opt',vt_dir)
	make_nubeam_opt(sim,'nubeam_opt')

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
	command = command + nubeamrun_exec+' ../../../%s %s %s \n'%(sim.e3.get(),sim.MenuVar1.get(),str(bdiff))

	make_batch_script(filename,None,log_e,log_o,command,'CHEASE_EMSE%03i'%round(bdiff*100.,0))	
	runid = submit_batch_script(filename)

	os.chdir(sim.currdir)

	return

def read_smse_result(sim,bdiff_index,fig,ind=2,type=1):

	len1 = len(bdiff_index)
	len2 = int(float(sim.StrVar106.get()))+1
	if len1 == 0:
		print('>>> No available result')
		return

	yy	= np.zeros((len1,len2,6))
	xx = np.zeros(len1)

	for i in range(len1):
		bdiff = round(float(bdiff_index[i]),2)
		xx [i] = bdiff
		save_dir = sim.currdir + '/%s/CSOLVE/Bdiff_%03i'%(sim.MenuVar1.get(),round(bdiff*100.,0))
		file_dir = save_dir + '/nubeam_iter_result'
		f = open(file_dir,'r')
		line = f.readline()
		j = 0
		while True:
			line = f.readline()
			if not line: break
			line = line.split()
			if j == len2:	break
			for k in range(6):
				yy[i,j,k] = float(line[k+1])
			j = j + 1

		if j < (len2-1):
			for t in range(j,len2+1):
				for k in range(6):
					yy[i,t,k] = yy[i,j-1,k]
		f.close()

	bp = yy[0,0,0]
	wmhd = yy[0,0,2]
	save_dir2 = sim.currdir + '/%s/CSOLVE/eq.info'%(sim.MenuVar1.get())
	f = open(save_dir2,'r')
	line = f.readline().split()
	f.close()
	rmag = float(line[0])
	if int(sim.MenuVar1.get()) > 0:
		f = open(sim.currdir + '/0/CSOLVE/eq.info','r')
		line = f.readline().split()
		f.close()
		rmag = float(line[0])
		try:    wmhd = float(line[1])/1.e3
		except:	pass
		try:	bp = bp * wmhd / yy[0,0,2]
		except:	pass
		
	fig.canvas.draw_idle()
	xdat = np.linspace(0,max(0.6,max(xx)),51)

	if type == 1:
		[ax1,ax2,ax3,ax4] = fig.axes
		ax1.cla();	ax2.cla();	ax3.cla();	ax4.cla();

		ax1.scatter(xx,yy[:,ind,0])
		ax1.set_title(r'$\beta_p$'+' [a.u]')
		ax1.set_xlabel('Diffusivity [a.u]')
		ax1.set_ylabel(r'$\beta_p$'+' [a.u]')
		ax1.axhline(y=bp,color='goldenrod',linestyle='--')
		if (len1 > 1):
			z = np.polyfit(xx,yy[:,ind,0],1)
			p = np.poly1d(z)
			ax1.plot(xdat,p(xdat),color='blueviolet',linestyle='--',linewidth = 1.0)

			val = (bp-p[0])/p[1]
			val = max(val,0.)

			sim.e201.configure(state='normal')
			sim.e201.delete(0,'end')
			sim.e201.insert(10,str(round(val,2)))
			sim.e201.configure(state='readonly')

		ax2.scatter(xx,yy[:,ind,2])
		ax2.set_title("$W_{MHD}$ [kJ]")
		ax2.set_xlabel('Diffusivity [a.u]')
		ax2.set_ylabel("$W_{MHD}$ [kJ]")
		ax2.axhline(y=wmhd,color='goldenrod',linestyle='--')
		if (len1 > 1):
			z = np.polyfit(xx,yy[:,ind,2],1)
			p = np.poly1d(z)
			ax2.plot(xdat,p(xdat),color='blueviolet',linestyle='--',linewidth = 1.0)		

			val = (wmhd-p[0])/p[1]
			val = max(val,0.)

			sim.e202.configure(state='normal')
			sim.e202.delete(0,'end')
			sim.e202.insert(10,str(round(val,2)))	
			sim.e202.configure(state='readonly')

		ax3.scatter(xx,yy[:,ind,3])
		ax3.set_title("$W_{DIA}$ [kJ]")
		ax3.set_xlabel('Diffusivity [a.u]')
		ax3.set_ylabel("$W_{DIA}$ [kJ]")
		ax3.axhline(y=float(sim.e6.get()),color='goldenrod',linestyle='--')
		if (len1 > 1):
			z = np.polyfit(xx,yy[:,ind,3],1)
			p = np.poly1d(z)
			ax3.plot(xdat,p(xdat),color='blueviolet',linestyle='--',linewidth = 1.0)	

			val = (float(sim.e6.get())-p[0])/p[1]
			val = max(val,0.)

			sim.e203.configure(state='normal')
			sim.e203.delete(0,'end')
			sim.e203.insert(10,str(round(val,2)))			
			sim.e203.configure(state='readonly')		

		ax4.scatter(xx,yy[:,ind,5])
		ax4.set_title("$R_{axis}$ [m]")
		ax4.set_xlabel('Diffusivity [a.u]')
		ax4.set_ylabel("$R_{axis}$ [m]")
		ax4.axhline(y=rmag,color='goldenrod',linestyle='--')
		if (len1 > 1):
			z = np.polyfit(xx,yy[:,ind,5],1)
			p = np.poly1d(z)
			ax4.plot(xdat,p(xdat),color='blueviolet',linestyle='--',linewidth = 1.0)

			val = (rmag-p[0])/p[1]
			val = max(val,0.)

			sim.e204.configure(state='normal')
			sim.e204.delete(0,'end')
			sim.e204.insert(10,str(round(val,2)))			
			sim.e204.configure(state='readonly')
	
	elif type == 2:
		xx2 = np.linspace(1,len2-1,len2-1)
		[ax1,ax2,ax3,ax4,ax5,ax6] = fig.axes
		ax1.cla();	ax2.cla();	ax3.cla();	ax4.cla();	ax5.cla();	ax6.cla();
		ax1.plot(xx2,yy[ind,1:len2+1,0],'--o')
		ax1.set_title(r'$\beta_p$'+' [a.u]')
		ax1.set_xlabel('Iteration [#]')
		ax1.set_ylabel(r'$\beta_p$'+' [a.u]')

		ax2.plot(xx2,yy[ind,1:len2+1,1],'--o')
		ax2.set_title(r'$l_i$'+' [a.u]')
		ax2.set_xlabel('Iteration [#]')
		ax2.set_ylabel(r'$l_i$'+' [a.u]')	
		
		ax3.plot(xx2,yy[ind,1:len2+1,2],'--o')
		ax3.set_title("$W_{MHD}$ [kJ]")
		ax3.set_xlabel('Iteration [#]')
		ax3.set_ylabel("$W_{MHD}$ [kJ]")
		
		ax4.plot(xx2,yy[ind,1:len2+1,3],'--o')
		ax4.set_title("$W_{DIA}$ [kJ]")
		ax4.set_xlabel('Iteration [#]')
		ax4.set_ylabel("$W_{DIA}$ [kJ]")	

		ax5.plot(xx2,yy[ind,1:len2+1,4],'--o')
		ax5.set_title("$W_{fast}$ [kJ]")
		ax5.set_xlabel('Iteration [#]')
		ax5.set_ylabel("$W_{fast}$ [kJ]")			

		ax6.plot(xx2,yy[ind,1:len2+1,5],'--o')
		ax6.set_title("$R_{axis}$ [m]")
		ax6.set_xlabel('Iteration [#]')
		ax6.set_ylabel("$R_{axis}$ [m]")
		

	fig.tight_layout()
	
	return

def read_emse_result(sim,bdiff_index,fig):

	len1 = len(bdiff_index)
	if len1 == 0:
		print('>>> No available result')
		return

	yy = np.zeros((len1,2,7))
	xx = np.zeros(len1)

	for i in range(len1):
		bdiff = round(float(bdiff_index[i]),2)
		xx [i] = bdiff
		save_dir = sim.currdir + '/%s/NUBEAM/Bdiff_%03i'%(sim.MenuVar1.get(),round(bdiff*100.,0))
		file_dir = save_dir + '/nubeam_iter_result'
		f = open(file_dir,'r')
		line = f.readline()
		for j in range(2):
			line = f.readline()
			if not line: break
			line = line.split()
			for k in range(7):
				yy[i,j,k] = float(line[k])

		f.close()

	wmhd = yy[0,0,3]
	save_dir2 = sim.currdir + '/%s/CSOLVE/eq.info'%(sim.MenuVar1.get())
	f = open(save_dir2,'r')
	line = f.readline().split()
	f.close()
	rmag = float(line[0])
	if int(sim.MenuVar1.get()) > 0:
		f = open(sim.currdir + '/0/CSOLVE/eq.info','r')
		line = f.readline().split()
		f.close()
		rmag = float(line[0])
		try:    wmhd = float(line[1])/1.e3
		except:	pass

	fig.canvas.draw_idle()
	xdat = np.linspace(0,max(0.6,max(xx)),51)
	[ax1,ax2,ax3,ax4] = fig.axes
	ax1.cla();	ax2.cla();	ax3.cla();	ax4.cla();

	ax1.scatter(xx,yy[:,1,6])
	ax1.set_title("$V_{Loop}$ [V]")
	ax1.set_xlabel('Diffusivity [a.u]')
	ax1.set_ylabel("$V_{Loop}$ [V]")
	if (len1 > 1):
		z = np.polyfit(xx,yy[:,1,6],1)
		p = np.poly1d(z)
		ax1.plot(xdat,p(xdat),color='blueviolet',linestyle='--',linewidth = 1.0)		

	ax2.scatter(xx,yy[:,1,3])
	ax2.set_title("$W_{MHD}$ [kJ]")
	ax2.set_xlabel('Diffusivity [a.u]')
	ax2.set_ylabel("$W_{MHD}$ [kJ]")
	ax2.axhline(y=wmhd,color='goldenrod',linestyle='--')
	if (len1 > 1):
		z = np.polyfit(xx,yy[:,1,3],1)
		p = np.poly1d(z)
		ax2.plot(xdat,p(xdat),color='blueviolet',linestyle='--',linewidth = 1.0)		

		val = (wmhd-p[0])/p[1]
		val = max(val,0.)

		sim.e202.configure(state='normal')
		sim.e202.delete(0,'end')
		sim.e202.insert(10,str(round(val,2)))	
		sim.e202.configure(state='readonly')

	ax3.scatter(xx,yy[:,1,4])
	ax3.set_title("$W_{DIA}$ [kJ]")
	ax3.set_xlabel('Diffusivity [a.u]')
	ax3.set_ylabel("$W_{DIA}$ [kJ]")
	ax3.axhline(y=float(sim.e6.get()),color='goldenrod',linestyle='--')
	if (len1 > 1):
		z = np.polyfit(xx,yy[:,1,4],1)
		p = np.poly1d(z)
		ax3.plot(xdat,p(xdat),color='blueviolet',linestyle='--',linewidth = 1.0)	

		val = (float(sim.e6.get())-p[0])/p[1]
		val = max(val,0.)

		sim.e203.configure(state='normal')
		sim.e203.delete(0,'end')
		sim.e203.insert(10,str(round(val,2)))			
		sim.e203.configure(state='readonly')		

	ax4.scatter(xx,yy[:,1,2])
	ax4.set_title("$W_{fast}$ [kJ]")
	ax4.set_xlabel('Diffusivity [a.u]')
	ax4.set_ylabel("$W_{fast}$ [kJ]")
	if (len1 > 1):
		z = np.polyfit(xx,yy[:,1,2],1)
		p = np.poly1d(z)
		ax4.plot(xdat,p(xdat),color='blueviolet',linestyle='--',linewidth = 1.0)		

	fig.tight_layout()

	return

def var_split(line):

	line2 = line.split()
	line3 = []
	for item in line2:
		item2 = item.strip().split(',')
		for it in item2: 
			if not it=='': line3.append(it)
	return line3

def get_kfile_data(data,var,firstonly=False,debug=False):

	var = var.lower()
	outdata = []
	istart = -1; invar = False
	for i in range(len(data)):
		if (data[i].find('=') > -1 and invar ): iend = i; invar = False;
		elif (data[i].lower().find(var) > -1 and istart==-1): istart = i; invar = True;
		elif (data[i].lower().find(var) > -1 and not firstonly): istart = i; invar = True;

	if invar: iend = i+1;
	iline = iend - istart
	for j in range(iline):
		if j==0: line = var_split(data[istart+j].split('=')[1])
		else: line = var_split(data[istart+j])
		if debug: print(line,data[istart+j])
		for item in line:
			if item.find('*')>-1:
				item2 = item.split('*')
				ii = int(item2[0])
				iv = float(item2[1])
				for i in range(ii): outdata.append(iv)
			else: outdata.append(float(item))
		
	out = np.array(outdata)
	return (out)
	
def read_kfile(sim,filename,extlib=False):

	filename2 = kindata_dir
	sim.kfile_in1 = []
	sim.kfile_in2 = []
	sim.kfile_in3 = []
	sim.kfile_in4 = []
	
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
		sim.endline = ''
		if not (line.find('/') > -1):

			if (line.find('MAG') == -1 and line.find('MSE') == -1 and line.find('KIN') == -1):
				if (count == 1):
					sim.kfile_in1.append(line)
				elif (count == 2 and line.find('&INWANT') == -1):
					sim.kfile_in2.append(line)
				elif (count == 3 ):
					sim.kfile_in3.append(line)
			else:
				sim.endline = line.split('M')[0]

	f4.close()
	
	sim.ismse = True
	if (count < 3):
		sim.fake_mse = True
		sim.ismse = False
	if not extlib: sim.CheckVar24.set(sim.trans_vars(sim.ismse,4))	
	f4 = open(filename2,'r')
	
	while True:
		line = f4.readline()
		if not line: break
		
		sim.kfile_in4.append(line)
	f4.close()
	
	sim.kfile_in01 = deepcopy(sim.kfile_in1)
	sim.kfile_in02 = deepcopy(sim.kfile_in2)
	sim.kfile_in03 = deepcopy(sim.kfile_in3)
	sim.kfile_in04 = deepcopy(sim.kfile_in4)

	return

def load_kfile(sim,coil_skip=False,extlib=False):

	if not coil_skip:
		sim.fwtsi  = get_kfile_data(sim.kfile_in1,'FWTSI')
		sim.fwtmp2 = get_kfile_data(sim.kfile_in1,'FWTMP2')
	if (sim.ismse):
		sim.tgamma = get_kfile_data(sim.kfile_in3,'TGAMMA')
		sim.sgamma = get_kfile_data(sim.kfile_in3,'SGAMMA')
		sim.fwtgam = get_kfile_data(sim.kfile_in3,'FWTGAM')
		sim.sgamma = get_kfile_data(sim.kfile_in3,'SGAMMA')
		sim.rrrgam = get_kfile_data(sim.kfile_in3,'RRRGAM')
		sim.zzzgam = get_kfile_data(sim.kfile_in3,'ZZZGAM')
		sim.aa1gam = get_kfile_data(sim.kfile_in3,'AA1GAM')
		sim.aa2gam = get_kfile_data(sim.kfile_in3,'AA2GAM')
		sim.aa3gam = get_kfile_data(sim.kfile_in3,'AA3GAM')
		sim.aa4gam = get_kfile_data(sim.kfile_in3,'AA4GAM')
		sim.aa5gam = get_kfile_data(sim.kfile_in3,'AA5GAM')
		sim.aa6gam = get_kfile_data(sim.kfile_in3,'AA6GAM')
	else:
		sim.sgamma = get_kfile_data(sim.kfile_in4,'SGAMMA')
		sim.fwtgam = get_kfile_data(sim.kfile_in4,'FWTGAM')
		sim.sgamma = get_kfile_data(sim.kfile_in4,'SGAMMA')
		sim.rrrgam = get_kfile_data(sim.kfile_in4,'RRRGAM')
		sim.zzzgam = get_kfile_data(sim.kfile_in4,'ZZZGAM')
		sim.aa1gam = get_kfile_data(sim.kfile_in4,'AA1GAM')
		sim.aa2gam = get_kfile_data(sim.kfile_in4,'AA2GAM')
		sim.aa3gam = get_kfile_data(sim.kfile_in4,'AA3GAM')
		sim.aa4gam = get_kfile_data(sim.kfile_in4,'AA4GAM')
		sim.aa5gam = get_kfile_data(sim.kfile_in4,'AA5GAM')
		sim.aa6gam = get_kfile_data(sim.kfile_in4,'AA6GAM')
		sim.tgamma = np.copy(sim.sgamma)

	sim.sgamma2 = np.copy(sim.sgamma)
	sim.dtgamma = np.zeros(len(sim.sgamma))

	sim.rpress = get_kfile_data(sim.kfile_in4,'RPRESS')
	sim.clen1 = len(sim.fwtmp2);	sim.clen2 = len(sim.fwtsi);	sim.clen3 = len(sim.fwtgam);
	print('>>> Probe signals ->','Magnetic [#]',sim.clen1,'PSI-loop [#]',sim.clen2,'MSE [#]',sim.clen3)

	if not extlib:
		if sim.CheckVar23.get() == 1:	return

	if coil_skip:	return
	
	if not extlib:
		for i in range(0,sim.clen1):
			sim.__dict__['CoilVar%d'%(i+1)].set(int(sim.fwtmp2[i]))
		for i in range(150,150+sim.clen2):
			sim.__dict__['CoilVar%d'%(i+1)].set(int(sim.fwtsi[i-150]))
		for i in range(200,200+sim.clen3):
			sim.__dict__['CoilVar%d'%(i+1)].set(int(sim.fwtgam[i-200]))

	return

def get_mse_mod(sim,filename):
	data = []
	if not os.path.isfile(filename):
		f = open(filename,'w')
		f.write('R[m] dtgam \n')
		print(sim.rrrgam)
		for i in range(len(sim.rrrgam)):
			f.write('%5.3f %5.3f \n'%(sim.rrrgam[i],0.))
		f.close()
		data = np.copy(sim.rrrgam)
		data = data* 0.;
	else:
		f = open(filename,'r')
		line = f.readline()
		while True:
			line = f.readline()
			if not line: break
			dtgam = float(line.split()[1])
			data.append(dtgam)
	data = np.array(data,dtype='float')
	return data

def make_kfile_str(data,var):

	out = ' ' + var+' = '
	len2 = len(data)
	count = 0
	for i in range(len2):
	
		count = count + 1
		out = out + '%19.13e,'%data[i]
		if (count == 7):
			count = 0
			out = out + '\n     '
			
		if (np.sum(abs(data)) == 0.0):
			out = ' ' + var+' = %i*0.000000000000000E+000  ,'%len(data)
			
	return (out)

def make_kfile(sim,mse_shift=True):

	rpress = np.zeros(51)
	if (sim.MenuVar8.get().lower() == 'hmode'):
		rpress0 = -1.0*np.linspace(0,0.9,25)
		rpress1 = -1.0*np.linspace(0.9,1.0,51-25)
		rpress[0:25]  = rpress0
		rpress[25:51] = rpress1

	else:
		rpress = -1.0*np.linspace(0,1.,51)
	
	sim.rpress = rpress	

	sim.kfile_in1 = deepcopy(sim.kfile_in01)
	sim.kfile_in2 = deepcopy(sim.kfile_in02)
	sim.kfile_in3 = deepcopy(sim.kfile_in03)
	sim.kfile_in4 = deepcopy(sim.kfile_in04)	

	sim.kfile_in1.append('\n  mxiter = -%s\n'%sim.StrVar38.get())
	sim.kfile_in1.append('  relax = %s \n'%sim.StrVar39.get())	
	sim.kfile_in1.append('  pcurbd = 0.0 \n')
	sim.kfile_in1.append('  fcurbd = 0.0 \n\n')

	pknot = sim.StrVar17.get().split(',');	fknot = sim.StrVar18.get().split(',');	jknot = sim.StrVar19.get().split(',')
	if sim.CheckVar4.get() == 1: jknot = np.array(jknot,dtype='double')
	else:	jknot = np.array([0,1.])
	sim.kfile_in1.append('  kppfnc  = 6, kppknt = %i, pptens = 1\n'%len(pknot))

	line = '  ppknt   = '
	for i in range(len(pknot)):
		line = line + '%4.3f '%float(pknot[i])
	sim.kfile_in1.append(line+'\n\n')

	sim.kfile_in1.append('  kfffnc  = 6, kffknt = %i, fftens = 1\n'%len(fknot))
	line = '  ffknt   = '
	for i in range(len(fknot)):
		line = line + '%4.3f '%float(fknot[i])
	sim.kfile_in1.append(line+'\n\n')
	
	sim.kfile_in1.append('  error   =  %s, \n'%sim.StrVar40.get())
	sim.kfile_in1.append('  errmin  =  %s, \n'%sim.StrVar40.get())
	sim.kfile_in1.append(' \n')
	sim.kfile_in1.append('  fwtcur  = %s \n'%sim.StrVar41.get())
	if sim.CheckVar5.get() == 1:	sim.kfile_in1.append('  fwtqa   = 1  \n')
	else:	sim.kfile_in1.append('  fwtqa   = 0  \n')

	if len(sim.StrVar42.get().split(',')) > 1:
		qvfit = float(sim.StrVar42.get().split(',')[0])
		mse_rshift = float(sim.StrVar42.get().split(',')[1])
	else:
		qvfit = float(sim.StrVar42.get().split(',')[0])
		mse_rshift = 0.

	sim.kfile_in1.append('  qvfit   = %s \n\n'%qvfit)		

	if (sim.CheckVar4.get() == 1):
		sim.kfile_in1.append('  RZEROJ  = %i*0 \n'%len(jknot))
		sim.kfile_in1.append('  KZEROJ  = %i \n\n'%len(jknot))

	len1 = len(sim.StrVar51.get().split(','))
	len2 = len(sim.StrVar52.get().split(','))

	if (len1 == len2 and len1 > 0. and sim.CheckVar6.get() == 1):	
		sim.kfile_in1.append('  NBDRY = %i\n'%(len1))
		sim.kfile_in1.append('  RBDRY = %s\n'%sim.StrVar51.get())
		sim.kfile_in1.append('  ZBDRY = %s\n'%sim.StrVar52.get())

	sim.kfile_in1.append('  fwtsi = \n')
	line = '  '
	count = 0
	for i in range(sim.clen2):
		count = count + 1
		line = line + '%i '%sim.__dict__['CoilVar%d'%(i+151)].get()
		if (count == 10):
			line = line + '\n  '
			count = 0;
	sim.kfile_in1.append(line)
	sim.kfile_in1.append('\n\n\n')

	sim.kfile_in1.append('  fwtmp2 = \n')
	line = '  '
	count = 0
	for i in range(sim.clen1):
		count = count + 1
		line = line + '%i '%sim.__dict__['CoilVar%d'%(i+1)].get()
		if (count == 20):
			line = line + '\n  '
			count = 0;
	sim.kfile_in1.append(line)
	sim.kfile_in1.append('\n\n\n')
	
	line = make_kfile_str(sim.rpress,'RPRESS')
	sim.kfile_in1.append(line)
	sim.kfile_in1.append('\n')
	
	presf = interp1d(sim.psin,sim.pres,'cubic')
	sim.prest = presf(-rpress)
	sim.sigpres = np.copy(sim.prest*0.1)
	sim.sigpres[0] = sim.prest[0]*0.03
	if len(sim.StrVar43.get().split(',')) > 1:
		core_sup = float(sim.StrVar43.get().split(',')[0])
		pres_mult = float(sim.StrVar43.get().split(',')[1])
	else:
		core_sup = float(sim.StrVar43.get().split(',')[0])
		pres_mult = 1.

	if float(core_sup) > 0.:
		sim.prest[0] = core_sup * sim.prest[1]

	pres_mult2 = np.copy(sim.prest)
	for i in range(len(pres_mult2)):
		pres_mult2[i] = 1. + (pres_mult-1) * 0.5 * (1.+np.tanh((0.8+rpress[i])/0.05))
	
	line = make_kfile_str(sim.prest*pres_mult2,'PRESSR')
	sim.kfile_in1.append(line)
	sim.kfile_in1.append('\n')

	line = make_kfile_str(sim.sigpres,'SIGPRE')
	sim.kfile_in1.append(line)
	sim.kfile_in1.append('\n\n')
	
	sim.kfile_in1.append(' NPRESS  =        %3i,\n'%len(sim.prest))
	sim.kfile_in1.append(' NBEAM   =        %3i,\n'%len(sim.prest))
	sim.kfile_in1.append('\n')

	#Dummy fast ion info (can be updated later)
	line = make_kfile_str(-rpress,'SIBEAM')
	sim.kfile_in1.append(line)
	sim.kfile_in1.append('\n')

	line = make_kfile_str(sim.prest*0.3,'PBEAM')
	sim.kfile_in1.append(line)
	sim.kfile_in1.append('\n')

	line = make_kfile_str(sim.prest/sim.prest[0]*4.e18,'DNBEAM')
	sim.kfile_in1.append(line)
	sim.kfile_in1.append('\n')
	
	line = make_kfile_str(sim.prest/sim.prest[0]*4.e18*1.627*1.e-27,'DMASS')
	sim.kfile_in1.append(line)
	sim.kfile_in1.append('\n\n')
	
	sim.kfile_in1.append(' NMASS  =        %3i,\n'%len(sim.prest))
	sim.kfile_in1.append(' KPRFIT   =      1,\n')
	sim.kfile_in1.append('\n')
	
	jconstf = interp1d(sim.jconst[:,0],sim.jconst[:,1],'cubic')
	if sim.CheckVar4.get()==1:
		jc = jconstf(jknot)
		line = make_kfile_str(jknot,'SIZEROJ')
		sim.line_jc = line + '\n'
		line = make_kfile_str(jc,'VZEROJ')
		sim.line_jc = sim.line_jc + line + '\n'

	if (sim.CheckVar3.get() == 1):	
		tgamma = np.copy(sim.tgamma2)
		print('>>> Use SMSE')
	else:	
		if sim.ts_run: tgamma = np.copy(sim.tgamma3)
		else: tgamma = np.copy(sim.tgamma)
		print('>>> USE EMSE')

	for i in range(171,171+len(sim.sgamma)):	
		sim.sgamma2[i-171] = float(sim.__dict__['StrVar%d'%i].get())
	for i in range(61,  61+len(sim.sgamma)):
		sim.dtgamma[i-61] = float(sim.__dict__['StrVar%d'%i].get())

	dtgamma = np.copy(sim.dtgamma)
	if (sim.CheckVar3.get() == 1 and sim.CheckVar7.get() == 0): dtgamma = np.zeros(len(sim.dtgamma))
	if (sim.CheckVar3.get() == 0 and sim.CheckVar7.get() == 1): dtgamma = np.zeros(len(sim.dtgamma))
	sim.kfile_in3t = []
	sim.kfile_in3t.append('  &INS\n')
	line=make_kfile_str(tgamma+dtgamma,'TGAMMA')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n')
	line=make_kfile_str(sim.sgamma2,'SGAMMA')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n')
	line=make_kfile_str(sim.fwtgam,'FWTGAM')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n')
	if (sim.CheckVar3.get() == 1):	
		line=make_kfile_str(sim.rrrgam+mse_rshift,'RRRGAM')
	else:	
		if mse_shift:	line=make_kfile_str(sim.rrrgam+mse_rshift+float(sim.StrVar44.get())/1.e2,'RRRGAM')
		else:	line=make_kfile_str(sim.rrrgam+mse_rshift,'RRRGAM')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n')
	line=make_kfile_str(sim.zzzgam,'ZZZGAM')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n')
	line=make_kfile_str(sim.aa1gam,'AA1GAM')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n')
	line=make_kfile_str(sim.aa2gam,'AA2GAM')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n')
	line=make_kfile_str(sim.aa3gam,'AA3GAM')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n')
	line=make_kfile_str(sim.aa4gam,'AA4GAM')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n')
	line=make_kfile_str(sim.aa5gam,'AA5GAM')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n')
	line=make_kfile_str(sim.aa6gam,'AA6GAM')
	sim.kfile_in3t.append(line)
	sim.kfile_in3t.append('\n\n')
	sim.kfile_in3t.append(' MSEBKP  =           0,\n')
	sim.kfile_in3t.append(' MSEFITFUN       =           1 \n')

	sim.kfile_in3t.append('\n')
	sim.kfile_in3t.append(' FWTGAM = \n')
	line = '  '
	count = 0
	for i in range(sim.clen3):
		count = count + 1
		line = line + '%i '%sim.__dict__['CoilVar%d'%(i+201)].get()
		if (count == 5):
			line = line + '\n  '
			count = 0;
	sim.kfile_in3t.append(line)
		
	return

def write_efit_input(sim,filename,type=1):

	f = open(filename,'w')
	
	if type == 1:	
		for item in sim.kfile_in1:
			f.write(item)
	else:		
		for item in sim.kfile_in01:
			f.write(item)
	f.write('/\n')
	f.write('&INWANT\n')
	if type == 1:	
		for item in sim.kfile_in2:
			f.write(item)
	else:		
		for item in sim.kfile_in02:
			f.write(item)
	if type == 1:
		if (sim.CheckVar4.get()):
			f.write(sim.line_jc)
		f.write('/\n')

	if type == 1:	
		for item in sim.kfile_in3t:
			f.write(item)
	else:
		for item in sim.kfile_in03:
			f.write(item)
		
	f.write('/\n')
	if type == 1:	f.write(sim.endline+' KIN')
	else:	f.write(sim.endline+' MSE')
	f.close()
	
	return

def get_efit_constraint_eq(sim,skip=False):

	if not sim.ismse:
		print('>>> No MSE data in KFILE, use synthetic MSE')
		sim.MenuVar7.set('sMSE')

	bdiff = round(float(sim.MenuVar6.get()),2)
	if not skip: print('>>> EFIT constraint with run mode = %s, iter = %s with bdiff = %s'%(sim.MenuVar7.get(),sim.MenuVar9.get(),sim.MenuVar6.get()))

	if sim.MenuVar7.get().lower() == 'smse':
		eq_file = sim.currdir + '/%s/CSOLVE/Bdiff_%03i/geqdsk_%s'%(sim.MenuVar1.get(),round(bdiff*100.,0),sim.MenuVar9.get())
		jc_file = sim.currdir + '/%s/CSOLVE/Bdiff_%03i/EFIT_JCONST_%s'%(sim.MenuVar1.get(),round(bdiff*100.,0),sim.MenuVar9.get())
		pc_file = ''
	else:
		eq_file = sim.e3.get()
		jc_file = sim.currdir + '/%s/NUBEAM/Bdiff_%03i/EFIT_JCONST'%(sim.MenuVar1.get(),round(bdiff*100.,0))
		pc_file = sim.currdir + '/%s/NUBEAM/Bdiff_%03i/pre_prof'%(sim.MenuVar1.get(),round(bdiff*100.,0))

	load_eq_data(sim,eq_file,pc_file,skip)
	get_efit_constraint_cu(sim,jc_file)

	return

def renew_efit_constraint_pj(eqfile,kinfile,beam_pf,beam_jf):

	run_dir = os.getcwd()

	ch2 = ch_tool.chease('chease_opt')
	ch2.eqdsk_name = eqfile
	ch2.chkin_file = kinfile
	ch2.load_eqdsk()
	ch2.read_hagermap('CHEASE/hager_map')
	ch2.adjust_prof = False

	ch2.pres_file = beam_pf
	ch2.curr_file = beam_jf

	if beam_pf == None: ch2.use_ext_pressure = False
	else: ch2.use_ext_pressure = True

	if beam_jf == None: ch2.use_ext_current = False
	else: ch2.use_ext_current = True

	ch2.kinprof_type = 1
	ch2.load_kin_profile()

	ch2.interpol()
	ch2.make_efit_constraint()

	fastpv = np.trapz(ch2.pres_ex, x=ch2.vol)
	fastppv = np.trapz(ch2.pres_exp,x=ch2.vol)
	thermp = np.trapz(ch2.pt,      x=ch2.vol)

	wmhd = 1.5 * (thermp + fastpv)  * 1.e-3
	wdia = 1.5 * (thermp + fastppv) * 1.e-3

	f = open('EFIT_PCONST','w')
	f.write('%i\n'%len(ch2.psin))
	for i in range(len(ch2.psin)): f.write('%9.6f\t%9.6f\n'%(ch2.psin[i],ch2.pt[i]+ch2.pres_ex[i]))
	f.close()
		
	return

def renew_efit_constraint_kin(sim,gfile0,gfile1):

	eq1 = eqdsk.eqdsk(gfile0)
	eq2 = eqdsk.eqdsk(gfile1)
	eq1.read_eqdsk(gfile0)
	eq2.read_eqdsk(gfile1)
	eq1.make_grid()
	eq2.make_grid()
	eq1.get_flux_contour()
	eq2.get_flux_contour()
	eq1.make_rho_R_psin()
	eq2.make_rho_R_psin()

	psi_n = np.linspace(0.,1.,401)
	rhof1 = interp1d(eq1.prhoR[:,0],eq1.prhoR[:,1],'quadratic')
	psif2 = interp1d(eq2.prhoR[:,1],eq2.prhoR[:,0],'quadratic')
#	for kk in range(len(eq1.prhoR[:,0])): print(eq1.prhoR[kk,0],eq1.prhoR[kk,1],eq2.prhoR[kk,0],eq2.prhoR[kk,1])

	rho_temp = rhof1(psi_n)
	psi_temp = psif2(rho_temp)
	mapf  = interp1d(psi_n,psi_temp,'quadratic')

	efitdir = sim.currdir +'/%s/EFIT/'%(sim.MenuVar1.get())
	shot = int(float(sim.e1.get()));       time = int(float(sim.e2.get()));
	ofile = 'chease_kinprof_new'
	f = open('map_init.dat','w')
	f.write('%i\n'%len(eq1.prhoR[:,0]))
	for kk in range(len(eq1.prhoR[:,0])): f.write('%9.6f\t%9.6f\n'%(eq1.prhoR[kk,0],eq1.prhoR[kk,1]))
	f.close()

	ifile =  efitdir+'p%06i.%06i'%(shot,time)

	f1 = open(ifile,'r')
	f2 = open(ofile,'w')
	for i in range(2):
		line = f1.readline()
		f2.write(line)

	while True:
		line = f1.readline()
		if not line: break
		line = line.split()
		psi_n = float(line[0])
		line[0] = mapf(psi_n)
		line2 = '%9.6f'%float(line[0])
		for i in range(len(line)-1): line2  = line2+'\t%9.6f'%float(line[i+1])
		line2 = line2 + '\n'
		f2.write(line2)
		
	f1.close()
	f2.close()

	return

def get_efit_constraint_cu(sim,filename,riter=False):

	f = open(filename,'r')
	line = f.readline()
	line = f.readline()
	num = int(float(line))


	if (sim.MenuVar7.get().lower() == 'smse' and not riter):	sim.jconst = np.zeros((num,2))
	else:	sim.jconst = np.zeros((num,3))

	if (sim.MenuVar7.get().lower() == 'smse' and not riter):

		for i in range(num):			
			line = f.readline().split()
			sim.jconst[i,0] = float(line[1])
			sim.jconst[i,1] = float(line[2])

	else:

		for i in range(num):

			line = f.readline().split()
			sim.jconst[i,0] = float(line[0])
			sim.jconst[i,1] = float(line[2])
			sim.jconst[i,2] = float(line[3])

	f.close()

	return

def save_eq_data(f_file,dat):

	f = open(f_file,'w')
	f.write('%i\n'%len(dat))
	for i in range(len(dat)):	f.write('%f\n'%dat[i])
	f.close()
	return

def load_eq_data(sim,filename,filename2=None,skip=False):

	eq = eqdsk.eqdsk(filename)
	eq.read_eqdsk_file(True)

	if not skip:
			
		f_file = sim.currdir +'/%s/CSOLVE/tgam.info'%sim.MenuVar1.get()
		
		if not os.path.isfile(f_file):	
			eq.read_eqdsk_file()
			eq.construct_2d_field()
		
			sim.tgamma2 = np.copy(sim.tgamma)
			brf = interp2d(eq.r,eq.z,eq.br)
			bzf = interp2d(eq.r,eq.z,eq.bz)
			btf = interp2d(eq.r,eq.z,eq.bt)		
	
			for i in range(len(sim.rrrgam)):
			
				R = sim.rrrgam[i]
				Z = sim.zzzgam[i]
				
				br = brf(R,Z)
				bz = bzf(R,Z)
				bt = btf(R,Z) * +1.0
				
				sim.tgamma2[i] = sim.aa1gam[i]*bz/(sim.aa2gam[i]*bt+sim.aa3gam[i]*br+sim.aa4gam[i]*bz)
			
			save_eq_data(f_file,sim.tgamma2)

		else:
			sim.tgamma2 = np.copy(sim.tgamma)
			f = open(f_file,'r')
			nn = int(f.readline())
			for i in range(nn):	sim.tgamma2[i] = float(f.readline())
			f.close()
		
	sim.pres = np.copy(eq.pres)
	sim.pp = np.copy(eq.pp)
	sim.ffp = np.copy(eq.ffp)
	sim.psin = np.copy(eq.psin)
	sim.bnd = np.copy(eq.rzbdy)
	sim.q2 = eq.q[0]

	if sim.MenuVar7.get().lower() == 'emse':

		sim.psin = np.linspace(0,1.,301)
		sim.pp = np.copy(sim.psin)
		sim.res = np.copy(sim.pp)
		num2 = 301

		f = open(filename2,'r')
		num = int(float(f.readline()))
		dat = np.zeros((num,2))
		for i in range(num):
			line = f.readline().split()
			dat[i,0] = float(line[0])
			dat[i,1] = float(line[1])
		f.close()

		pf = interp1d(dat[:,0],dat[:,1],'cubic')
		deps = 1.e-5
		for i in range(1,len(sim.psin)-1):
			val1 = pf(sim.psin[i]-deps)
			val2 = pf(sim.psin[i]+deps)

			sim.pp[i] = (val2-val1)/2./deps

		sim.pp[0] = 2.0*sim.pp[1] - sim.pp[2]
		sim.pp[-1] = 2.0*sim.pp[num2-2] - sim.pp[num2-3]

		sim.pres = pf(sim.psin)

	return

def draw_efit_constraint(sim,fig):

	[ax1,ax2,ax3,ax4] = fig.axes

	sim.ax1 = ax1
	sim.ax2 = ax2
	sim.ax3 = ax3
	sim.ax4 = ax4

	draw_ff_constraint(sim)
	draw_pp_constraint(sim)
	draw_j_constraint(sim)
	draw_mse_constraint(sim)

	fig.tight_layout()

	return

def draw_ff_constraint(sim):
	sim.fig6.canvas.draw_idle(); sim.ax1.cla();	
	if sim.MenuVar7.get().lower() == 'smse':	
		sim.ax1.plot(sim.psin,sim.ffp)
		draw_knots(sim.psin,sim.ffp,sim.e18.get(),sim.ax1)
	else:	
		sim.ax1.plot(sim.jconst[:,0],sim.jconst[:,2])
		draw_knots(sim.jconst[:,0],sim.jconst[:,2],sim.e18.get(),sim.ax1)

	sim.ax1.set_title("FF'")
	sim.ax1.set_xlabel('$\psi_N$ [a.u]')
	sim.ax1.set_ylabel("FF' [a.u]")
	

	if sim.MenuVar7.get().lower() == 'smse':	
		ind = np.where(sim.psin>0.2); minv = min(sim.ffp[ind]); maxv = max(sim.ffp[ind])
	else:	
		ind = np.where(sim.jconst[:,0]>0.2); minv = np.min(sim.jconst[ind,2]); maxv = np.max(sim.jconst[ind,2])

	maxv = maxv*1.5
	if minv < 0.:	minv = 1.5*minv
	else:	minv = 0.5*minv

	sim.ax1.set_xlim([0.2,1.05])
	sim.ax1.set_ylim([minv, maxv])
	sim.StrVar18.set(sim.e18.get())

	return

def draw_j_constraint(sim):
	sim.fig6.canvas.draw_idle(); sim.ax2.cla();	
	sim.ax2.plot(sim.jconst[:,0],sim.jconst[:,1])
	draw_knots(sim.jconst[:,0],sim.jconst[:,1],sim.e19.get(),sim.ax2)
	sim.ax2.set_title("Current constraint")
	sim.ax2.set_xlabel('$\psi_N$ [a.u]')
	sim.ax2.set_ylabel("Current constraint [a.u]")
	sim.ax2.set_xlim([0.2,1.05])
	sim.StrVar19.set(sim.e19.get())
	return

def draw_pp_constraint(sim):
	sim.fig6.canvas.draw_idle(); sim.ax3.cla();
	sim.ax3.plot(sim.psin,sim.pp/min(sim.pp))
	draw_knots(sim.psin,sim.pp/min(sim.pp),sim.e17.get(),sim.ax3)
	sim.ax3.set_title("P'")
	sim.ax3.set_xlabel('$\psi_N$ [a.u]')
	sim.ax3.set_ylabel("P' [a.u]")
	sim.ax3.set_xlim([0.2,1.05])
	sim.StrVar17.set(sim.e17.get())
	return

def draw_mse_constraint(sim):
	sim.fig6.canvas.draw_idle(); sim.ax4.cla();	
	sim.mse_shift = float(sim.e44.get())/1.e2

	if float(sim.StrVar151.get()) > 0.:   mseb = 'NB1A'; ind = 152
	elif float(sim.StrVar153.get()) > 0.: mseb = 'NB1B'; ind = 154
	else: mseb = 'NB1C'; ind=156
	sshot = int(float(sim.StrVar1.get()))
	if (sshot >= 20756 and sshot < 21759): mseb = 'NB1B'; ind = 154
	if (sshot >= 23057 and sshot < 23134): mseb = 'NB1C'; ind = 156
	if (sshot >= 25283):
		if sim.aa1gam[0] < 0.83: mseb = 'NB1A'; ind = 152
		else: mseb = 'NB1B'; ind = 154
	print('>>> Beam used in MSE', mseb)

	xerr = mse_rad_err[mseb]/1.e3

	dtgamma = np.copy(sim.dtgamma); dtgamma2 = np.copy(sim.dtgamma);
	if (sim.CheckVar7.get()==0): dtgamma2 = np.zeros(len(sim.dtgamma))
	else: dtgamma = np.zeros(len(sim.dtgamma))

	xmag1 = get_mse_center(sim.rrrgam,sim.tgamma2+dtgamma2)
	indp  = np.where(sim.rrrgam>1.0)
	for i in range(sim.clen3):
		if (sim.rrrgam[i] < 1.0 and sim.__dict__['CoilVar%d'%(i+201)].get() == 1):
			sim.__dict__['CoilVar%d'%(i+201)].set(0)

	if sim.MenuVar7.get().lower() == 'emse':	
		line1, = sim.ax4.plot(sim.rrrgam[indp],sim.tgamma2[indp],linestyle='--',color='lime')	
	else:
		for i in range(sim.clen3):
			if sim.rrrgam[i] < 1.0: continue
			if sim.__dict__['CoilVar%d'%(i+201)].get() == 1:
				line1 = sim.ax4.scatter(sim.rrrgam[i],sim.tgamma2[i]+dtgamma2[i],c='lime')
			else:
				line1 = sim.ax4.scatter(sim.rrrgam[i],sim.tgamma2[i]+dtgamma2[i],c='gray')

	datx = [];	daty = [];
	if sim.ismse:	
		for i in range(sim.clen3):
			if sim.rrrgam[i] < 1.0: continue
			if sim.__dict__['CoilVar%d'%(i+201)].get() == 1:
				line2 = sim.ax4.errorbar(sim.rrrgam[i]+sim.mse_shift,sim.tgamma[i]+dtgamma[i],xerr= xerr[i],yerr = 1.5*sim.sgamma2[i],fmt='*',markersize='7',c='magenta',ecolor='r',capthick=1)
				datx.append(sim.rrrgam[i]+sim.mse_shift);	daty.append(sim.tgamma[i]+dtgamma[i]);
			else:
				line2 = sim.ax4.errorbar(sim.rrrgam[i]+sim.mse_shift,sim.tgamma[i]+dtgamma[i],xerr= xerr[i],yerr = 1.5*sim.sgamma2[i],fmt='*',markersize='7',c='gray',ecolor='gray',capthick=1)


	datx = np.array(datx,dtype='double');	daty = np.array(daty,dtype='double');
	#xmag2 = get_mse_center(datx,daty)
	try:	xmag2 = get_mse_center(datx,daty)
	except:	xmag2 = 0.

	sim.ax4.set_title("TGAMMA-MSE")
	sim.ax4.set_xlabel('R [m]')
	sim.ax4.set_ylabel("TGAMMA-MSE")
	max1 = max(sim.tgamma2)	* 2.0
	min1 = min(sim.tgamma2)	* 2.0
	sim.ax4.set_ylim([min1,max1])
	if sim.ismse:	sim.ax4.legend([line1,line2],['Synthetic-MSE ($R_{MAG}$=%4.3f [m])'%xmag1,'Experimental-MSE ($R_{MAG}$=%4.3f [m])'%xmag2])
	else:	sim.ax4.legend([line1],['Synthetic-MSE ($R_{MAG}$=%4.3f)'%xmag1])
	return

def draw_knots(x,y,knots,ax):

	if (len(knots.split()) == 0.):	return

	knots = knots.split(',')
	yf = interp1d(x,y,'cubic')
	xx = len(knots)
	xx = np.zeros(xx)
	for i in range(len(xx)):
		xx[i] = float(knots[i])

	ax.scatter(xx,yf(xx),color='magenta',marker='+',s=50)
	return

def draw_efit_riter(sim,ax1,ax2,ax3,ax4):
	
	cmap = ['#ff7f0e','#2ca02c','m','#1f77b4','#d62728','#8c564b','#bcbd22']
	index = int(sim.MenuVar10.get())

	sim.fig18.canvas.draw_idle();
	ax1.cla(); ax2.cla(); ax3.cla(); ax4.cla()

	efitdir  = sim.currdir +'/%s/EFIT/'%(sim.MenuVar14.get())
	temp_dir = efitdir+ '/RESULT/EFIT_ITER_%i/'%index

	if not os.path.isdir(temp_dir):
		print('>>> No RIteration history...')
		return

	shot = int(float(sim.e1.get()));       time = int(float(sim.e2.get()));

	os.chdir(temp_dir)
	f = open(temp_dir+'/niter','r')
	riter = int(f.readline())
	f.close()

	slegend = []
	sline = dict()
	for k in range(1,13): 
		sline[k] = dict()
		for j in range(1,6): sline[k][j] = '{:12s}'.format('')

	for k in range(0,6): sline[2][k] = '\n{:12s}'.format('')

	sline[1][0] = '{:12s}'.format('# ITERATION ')
	sline[2][0]='\n{:12s}'.format('XI2     [#]')
	sline[3][0] = '{:12s}'.format('Conv.   [#]')
	sline[4][0] = '{:12s}'.format('Q0      [#]')
	sline[5][0] = '{:12s}'.format('Q95     [#]')
	sline[6][0] = '{:12s}'.format('IP     [MA]')
	sline[7][0] = '{:12s}'.format('WMHD   [kJ]')
	sline[8][0] = '{:12s}'.format('BETAP   [#]')
	sline[9][0] = '{:12s}'.format('BETA    [%]')
	sline[10][0]= '{:12s}'.format('LI      [#]')
	sline[11][0]= '{:12s}'.format('RMAG    [m]')
	sline[12][0]= '{:12s}'.format('ROUT    [m]')
	
	for k in range(riter+1):

		f3 = open('wmhd.dat_%i'%k,'r')		
		f2 = open('qpres_%i.dat'%k,'r')
		f1 = open('prho_%i.dat'%k,'r')
		
		f4 = open('pres.dat_%i'%k,'r')
		f5 = open('jconst.dat_%i'%k,'r')

		datn = int(f1.readline())
		xy = np.zeros((datn,2))
		for kk in range(datn):
			line = f1.readline().split()
			xy[kk,0] = float(line[0])
			xy[kk,1] = float(line[1])
		ax1.plot(xy[:,0],xy[:,1],color=cmap[k],linestyle='--')
		datn = int(f2.readline())
		xy = np.zeros((datn,3))
		for kk in range(datn):
			line = f2.readline().split()
			xy[kk,0] = float(line[0])
			xy[kk,1] = float(line[1])*1.e-3
			xy[kk,2] = float(line[2])
		ax2.plot(xy[:,0],xy[:,1],color=cmap[k],linestyle='--')
	#	ax4.plot(xy[:,0],xy[:,2],color=cmap[k],linestyle='--')
		read_efit_result(sim,'fitout.dat_%i'%k)

		wline = f3.readline().split()
		ip = float(wline[3])*1.e6
		area = float(wline[4])		

		ax3.plot(sim.efit_psin,sim.efit_jav*area/ip,color=cmap[k],linestyle='--')	
		datn = int(f4.readline())
		xy = np.zeros((datn,3))
		for kk in range(datn):
			line = f4.readline().split()
			xy[kk,0] = -float(line[0])
			xy[kk,1] = float(line[1])*1.e-3
			xy[kk,2] = float(line[2])*1.e-3

		ax2.errorbar(xy[:,0],xy[:,1],xy[:,2],fmt='D',markersize=2,c=cmap[k],capsize=3,capthick=0.5)
		datn = int(f5.readline())
		xy = np.zeros((datn,2))
		for kk in range(datn):
			line = f5.readline().split()
			xy[kk,0] = float(line[0])
			xy[kk,1] = float(line[1])
		ax3.scatter(xy[:,0],xy[:,1],s=50,marker='+',color=cmap[k])
		f1.close()
		f2.close()
		f3.close()
		f4.close()
		f5.close()
		qf = interp1d(sim.efit_psin,sim.efit_q)

		for j in range(1,6):			
			if j == (k+1):
				sline[1][j]  = sline[1][j]  + ' {:7d}'.format(k+1)
				sline[2][j]  = sline[2][j]  + ' {:7.2f}'.format(sim.efit_chi)
				sline[3][j]  = sline[3][j]  + ' {:7.0e}'.format(sim.efit_conve)
				sline[4][j]  = sline[4][j]  + ' {:7.3f}'.format(qf(0))
				sline[5][j]  = sline[5][j]  + ' {:7.3f}'.format(qf(0.95))
				sline[6][j]  = sline[6][j]  + ' {:7.3f}'.format(-sim.efit_ip*1.e-6)
				sline[7][j]  = sline[7][j]  + ' {:7.1f}'.format(float(wline[2]))
				sline[8][j]  = sline[8][j]  + ' {:7.3f}'.format(sim.efit_betap)
				sline[9][j]  = sline[9][j]  + ' {:7.3f}'.format(sim.efit_betat)
				sline[10][j] = sline[10][j] + ' {:7.3f}'.format(sim.efit_li)
				sline[11][j] = sline[11][j] + ' {:7.3f}'.format(float(wline[1]))
				sline[12][j] = sline[12][j] + ' {:7.3f}'.format(float(wline[0]))
			else:
				for kk in range(1,13): sline[kk][j] = sline[kk][j] + ' {:7s}'.format('')

		slegend.append('Iter #%i'%(k+1))

	for j in range(riter+2,6):
		for k in range(1,6):
			if j == (k+1):
				sline[1][j]  = sline[1][j]  + ' {:7d}'.format(k+1)
				for kk in range(2,13): sline[kk][j] = sline[kk][j] + ' {:7s}'.format('      -')
			else:
				for kk in range(1,13): sline[kk][k] = sline[kk][k] + ' {:7s}'.format('')

	ax4.set_xlim([0,1])
	ax4.set_ylim([0,1])
	ax4.axis('off')
		
	for j in range(6):
		tline = ''
		for k in range(1,13): tline = tline + sline[k][j] + '\n'
		if j == 0:
			ax4.text(0,1,tline,verticalalignment='top',horizontalalignment='left',transform = ax4.transAxes,linespacing=1.2,family='monospace')
		else: 
			ax4.text(0,1,tline,verticalalignment='top',horizontalalignment='left',transform = ax4.transAxes,linespacing=1.2,family='monospace',color=cmap[j-1])
		
	f = open('map_init.dat','r')
	datn = int(f.readline())
	xy = np.zeros((datn,2))
	for kk in range(datn):
		line = f.readline().split()
		xy[kk,0] = float(line[0])
		xy[kk,1] = float(line[1])

	xx= np.linspace(0,1.,41)
	temp = interp1d(xy[:,0],xy[:,1])
	yy = temp(xx)
	ax1.scatter(xx,yy,color='gray',marker='x',s=20)


	ax1.set_title('Radial mapping')
	ax1.set_ylabel('$\\rho_N$ [a.u]')
	ax2.set_title('Pressure constraint')
	ax2.set_ylabel('Pressure [kPa]')
	ax3.set_title('Current constraint')
	ax3.set_ylabel('Nomalized density [a.u]')
	#ax4.set_title('q-profile')
	#ax4.set_ylabel('q [a.u]')

	ax1.set_xlabel('$\psi_N$ [a.u]')
	ax2.set_xlabel('$\psi_N$ [a.u]')
	ax3.set_xlabel('$\psi_N$ [a.u]')
	#ax4.set_xlabel('$\psi_N$ [a.u]')
	ax2.legend(slegend)
	ax3.legend(slegend)
	slegend.append('Reference')
	ax1.legend(slegend)
	#ax4.legend(slegend)

	sim.fig18.tight_layout()
	os.chdir(sim.currdir)

	return

def make_pf_knots(sim,type):
	
	kn = knots_tool3.knots()
	if type == 1:
		kn.spline_start = float(sim.StrVar20.get())		
		kn.end_knot = float(sim.StrVar21.get())			
		kn.start_knot = float(sim.StrVar22.get())			
		kn.knots_shift = float(sim.StrVar23.get())	
		kn.coren = int(float(sim.StrVar24.get()))
		kn.edgen = int(float(sim.StrVar25.get()))				
		kn.input_datan = int(float(sim.StrVar26.get()))
		kn.delmin_core = float(sim.StrVar46.get())
		kn.delmin_edge = float(sim.StrVar47.get())
		kn.minloc = 0.
		if sim.MenuVar8.get().lower() == 'hmode':	knot = kn.lmfit_fit(sim.psin,sim.pp/min(sim.pp),2,True)
		else:	knot = kn.lmfit_fit(sim.psin,sim.pp/min(sim.pp),1,True)
		
		if sim.MenuVar8.get().lower() == 'hmode':
			pp = sim.pp * np.sign(sim.pp[-5])
			ind = np.where(sim.psin>0.8)
			loc = np.argmax(pp[ind])
			ind1 = np.where(sim.psin[ind]<sim.psin[ind][loc])
			loc = np.argmin(pp[ind][ind1])
			sim.minloc = min(sim.psin[ind][ind1][loc]+0.03,0.92)
		else:
			sim.minloc = 1.
	else:
		kn.spline_start = float(sim.StrVar27.get())		
		kn.end_knot = float(sim.StrVar28.get())			
		kn.start_knot = float(sim.StrVar29.get())			
		kn.knots_shift = float(sim.StrVar30.get())	
		kn.coren = int(float(sim.StrVar31.get()))
		kn.edgen = int(float(sim.StrVar32.get()))				
		kn.input_datan = int(float(sim.StrVar33.get()))
		kn.delmin_core = float(sim.StrVar48.get())
		kn.delmin_edge = float(sim.StrVar49.get())
		try:	kn.minloc = sim.minloc
		except:
			print('>>> Do pressure first!')
			return
		if sim.MenuVar8.get().lower() == 'hmode':	
			if sim.MenuVar7.get().lower() == 'smse':	knot = kn.lmfit_fit(sim.psin,sim.ffp,2,False)
			else:	knot = kn.lmfit_fit(sim.jconst[:,0],sim.jconst[:,2],2,False)
		else:
			if sim.MenuVar7.get().lower() == 'smse':	knot = kn.lmfit_fit(sim.psin,sim.ffp,1,False)
			else:	knot = kn.lmfit_fit(sim.jconst[:,0],sim.jconst[:,2],1,False)			

	sknot = '%s'%round(knot[0],3)
	for i in range(1,len(knot)):
		sknot = sknot + ',%s'%round(knot[i],3)

	return sknot

def find_max_current(xx,yy):

	yf = interp1d(xx,yy)
	scan_xs = np.linspace(0.8,1.0,6,'linear')
	maxy= np.zeros(scan_xs.shape[0]-1)
	maxx= np.ones(scan_xs.shape[0]-1) * 0.95
	
	for i in range(scan_xs.shape[0]-1):
	
		xxs = np.linspace(scan_xs[i],scan_xs[i+1],501)
		yys = yf(xxs)
		maxind = np.argmax(yys) 

		if (maxind>0 and maxind<(xxs.shape[0]-1)):
			maxy[i] = yys[maxind]
			maxx[i] = xxs[maxind]

	maxind = np.argmax(maxy)
	return maxx[maxind]

def make_current_knots(sim,coren,edgen,knots,knote):

	xx = np.copy(sim.jconst[:,0]);	yy = np.copy(sim.jconst[:,1]);

	ind = np.where(xx>0.8);	maxloc = xx[ind][np.argmax(yy[ind])];
	maxloc = find_max_current(xx,yy)
	ind = np.where(xx<maxloc);	minloc = xx[ind][np.argmin(yy[ind])];
	print('>>> Min current location = %4.3f'%minloc)
	print('>>> Max current location = %4.3f'%maxloc)
	if edgen < 6:
		print('>>> Edge [#] > 5 is needed')
		edgen = 6;

	x = [maxloc, knote-0.006, knote-0.003, knote]

	pedn1 = (edgen-5)/2
	if not pedn1 == int(pedn1):
		pedn1 = int(pedn1-0.5);	pedn2 = pedn1 + 1;
	else:
		pedn1 = int(pedn1);	pedn2 = pedn1;

	dx = np.linspace(maxloc,knote-0.006,pedn1+2)
	for i in range(1,pedn1+1):	x.append(dx[i])
	
	dx = np.linspace(min(minloc,maxloc),maxloc,pedn2+2)
	for i in range(0,pedn2+1):	x.append(dx[i])

	dx = np.linspace(minloc-0.05,min(minloc+0.03,maxloc-0.01),int(pedn2/2)+2)
	for i in range(0,int(pedn2/2)+1):	x.append(dx[i])

	if coren < 2:
		print('>>> Edge [#] > 1 is needed')
		coren = 2;	

	x.append(knots)
	dx = np.linspace(knots,minloc-0.05,coren+1)
	for i in range(1,coren):	x.append(dx[i])

	x = np.sort(x)
	knot = '%s'%round(x[0],3)
	for i in range(1,len(x)):
		knot = knot + ',%s'%round(x[i],3)

	return knot

def get_mse_center(datx,daty):

	ind = np.where(abs(daty)<0.5)
	
	xx = datx[ind];	yy = daty[ind]
	ind2 = np.where(xx>1.);
	xx = xx[ind2]; yy = yy[ind2];

	datf = interp1d(xx,yy,'cubic')
	xmin = min(xx)
	xmax = max(xx)

	x0 = 0.5*(xmin+xmax)
	err0 = datf(x0)

	if err0>0.:	x1 = get_mse_center2(xmin,xmax,x0+0.01)
	else:	x1 =  get_mse_center2(xmin,xmax,x0-0.01)
	err1 = datf(x1)
	count = 0
	while (abs(err1) > 1.e-5 and count < 100):
		x2 = (x1-x0)/(err1-err0+1.e-5)*(0.-err0) + x0
		x2 = get_mse_center2(xmin,xmax,x2)
		x0 = x1;
		x1 = x2;
		err0 = err1
		err1 = datf(x1)
		count = count + 1

	return x1

def get_mse_center2(minx,maxx,x):

	x = max(minx,x)
	x = min(maxx,x)
	return x

def find_efit_runs(sim):

	efitdir = sim.currdir +'/%s/EFIT/RESULT/'%(sim.MenuVar14.get())
	shot = int(float(sim.e1.get()));	time = int(float(sim.e2.get()));
	kfile_dir = efitdir + 'k%06i.%06i_kin'%(shot,time)
	gfile_dir = efitdir + 'g%06i.%06i_kin'%(shot,time)
	list1 = []
	for i in range(0,100):
		gfile_dir2 = gfile_dir + '_%i'%i
		if os.path.isfile(gfile_dir2):	list1.append(i)

	if len(list1) == 0.:	list1 = ['--']

	return list1

def read_efit_result(sim,filename):

	f = open(filename,'r')

	sim.chi_ploop = np.zeros(sim.clen2)
	sim.chi_mloop = np.zeros(sim.clen1) 

	sim.m_ploop = np.zeros(sim.clen2)
	sim.m_mloop = np.zeros(sim.clen1)
	sim.c_ploop = np.zeros(sim.clen2)
	sim.c_mloop = np.zeros(sim.clen1)

	clen1 = int(np.ceil(sim.clen1/8))
	clen2 = int(np.ceil(sim.clen2/8))
	clen3 = int(np.ceil(sim.clen3/8))
	sim.efit_chi = 0.
	sim.efit_num = 1
	
	sim.efit_conve = 0.
	sim.efit_bp = 0.

	count = 0
	while True:
		line = f.readline()
		if not line:	break
		if line.find('chi psi loops') > -1:
			for i in range(sim.clen2):
				if (int(i/8)==(i/8)):	
					line = f.readline().split()
					j = 0
				sim.chi_ploop[i] = float(line[j])
				j = j + 1;

		elif line.find('inner magnetic probes') > -1:
			for i in range(sim.clen1):
				if (int(i/8)==(i/8)):	
					line = f.readline().split()
					j = 0
				sim.chi_mloop[i] = float(line[j])
				j = j + 1;
		elif line.find('shot #') > -1:
			sim.efit_chi = float(line.split('=')[-1])
		elif line.find('betap') > -1:
			sim.efit_bp = float(line.split('=')[2].split('li')[0])

		elif line.find('EFITD') > -1:
			sim.efit_num = int(line.split('EFITD')[1].split('dx2')[0])
			sim.efit_psin = np.zeros(sim.efit_num)
			sim.efit_pp = np.zeros(sim.efit_num)
			sim.efit_ffp = np.zeros(sim.efit_num)
			sim.efit_q = np.zeros(sim.efit_num)
			sim.efit_jav = np.zeros(sim.efit_num)
			for k in range(8): line = f.readline()
			line = f.readline().split()
			sim.efit_chi   = float(line[9])
			line = f.readline().split()
			sim.efit_betat = float(line[2])
			sim.efit_betap = float(line[5])
			sim.efit_bp = float(line[5])
			sim.efit_li    = float(line[8])


		elif line.find('plasma summary') > -1:
			count = 1
			line = f.readline()
			for i in range(sim.efit_num):
				line = f.readline().split()
				sim.efit_psin[i] = float(line[1])
				sim.efit_pp[i] = float(line[3])
				sim.efit_ffp[i] = float(line[5])
				sim.efit_q[i] = float(line[8])

		elif (line.find('pflux') > -1 and count == 1):
			for i in range(sim.efit_num):
				line = f.readline().split()
				sim.efit_jav[i] = float(line[2])

		elif line.find('calculated psi-loops') > -1:
			for i in range(sim.clen2):
				if (int(i/4)==(i/4)):	
					line = f.readline().split()
					j = 0
				sim.c_ploop[i] = float(line[j])
				j = j + 1;
		elif line.find('measured psi-loops') > -1:
			for i in range(sim.clen2):
				if (int(i/4)==(i/4)):	
					line = f.readline().split()
					j = 0
				sim.m_ploop[i] = float(line[j])
				j = j + 1;

		elif line.find('calculated total plasma current') > -1:
			sim.efit_ip = float(f.readline())

		elif line.find('calculated magnetic probes') > -1:
			for i in range(sim.clen1):
				if (int(i/4)==(i/4)):	
					line = f.readline().split()
					j = 0
				sim.c_mloop[i] = float(line[j])
				j = j + 1;
		elif line.find('measured magnetic probes') > -1:
			for i in range(sim.clen1):
				if (int(i/4)==(i/4)):	
					line = f.readline().split()
					j = 0
				sim.m_mloop[i] = float(line[j])
				j = j + 1;				

		elif line.find('iteration') > -1:
			line = f.readline().split()
			for i in range(999):
				try:	line = f.readline().split()
				except:	line = []
				if len(line) == 0.: break
				sim.efit_conve = float(line[1])
	
	f.close()

	return

def draw_efit_runs(sim,filename,map_name,ax1,ax2,ax3,ax4,index,skip=False):
	read_efit_result(sim,filename)
	if skip:	return
	sim.fig8.canvas.draw_idle()
	sim.efit_pp = sim.efit_pp * np.sign(sim.efit_pp[20])
	sim.efit_jav = sim.efit_jav * np.sign(sim.efit_jav[0])

	line1, = ax1.plot(sim.efit_psin,sim.efit_pp/1.e3)
	ax2.plot(sim.efit_psin,sim.efit_ffp)
	ax3.plot(sim.efit_psin,sim.efit_jav/1.e6)
	ax4.plot(sim.efit_psin,sim.efit_q)
	ax4.axhline(y=1.,color='gold',linestyle='--',linewidth = 1.0);

	radius,ismap = find_q1_surface(sim.efit_psin,sim.efit_q,map_name)

	if not ismap: line = ', q=1 at $\psi_N$ = '
	else: line = ', q=1 at R[m] = '
	if radius[0]==-1:	line = ''
	else: 
		for i in range(len(radius)): line = line + '%3.2f '%radius[i]

	sim.efit_plg1.append(line1)
	sim.efit_slg1.append('#%s-%i($\chi^2$=%3.1f)'%(sim.MenuVar14.get(),index,sim.efit_chi))
	sim.efit_slg12.append('#%s-%i($\chi^2$=%3.1f)'%(sim.MenuVar14.get(),index,sim.efit_chi)+line)

	ax1.set_title("P'")
	ax1.set_xlabel('Normalized radius ($\psi_N$)')
	ax1.set_ylabel("P' [kPa/Wb]")
	ax2.set_title("FF'")
	ax2.set_xlabel('Normalized radius ($\psi_N$)')
	ax2.set_ylabel("FF' [a.u]")	
	ax3.set_title("Averaged current density")
	ax3.set_xlabel('Normalized radius ($\psi_N$)')
	ax3.set_ylabel("<j> [MA/m2]")	
	ax4.set_title("q")
	ax4.set_xlabel('Normalized radius ($\psi_N$)')
	ax4.set_ylabel("q [a.u]")
	

	ax1.legend(sim.efit_plg1,sim.efit_slg1)
	ax2.legend(sim.efit_plg1,sim.efit_slg1)
	ax3.legend(sim.efit_plg1,sim.efit_slg1)
	ax4.legend(sim.efit_plg1,sim.efit_slg12)

	sim.fig8.tight_layout()
	
	return

def find_q1_surface(psin,qq,map_dir):

	ismap = False
	if os.path.isfile(map_dir):
		f = open(map_dir,'r')
		datlen = int(f.readline())
		psin2 = np.zeros(datlen)
		RP = np.zeros(datlen)
		for i in range(datlen):
			line = f.readline().split()
			psin2[i] = float(line[0])
			RP[i]   = float(line[2])
		f.close()
		pRf = interp1d(psin2,RP)
		ismap = True

	qtarget = 1.
	if min(qq) > qtarget:	return [-1], ismap
	radius = np.array([])
	q1r_old = -1
	if not (len(qq)==len(psin)): print('len qq/psin',len(qq),len(psin))
	for i in range(len(qq)-1):
		if ((qq[i]-qtarget)*(qq[i+1]-qtarget) <= 0.):	
			q1r = (psin[i+1]-psin[i])/(qq[i+1]-qq[i])*(qtarget-qq[i]) + psin[i]
			if ismap: q1r = pRf(q1r)
			if not (q1r == q1r_old):
				q1r_old = q1r
				radius = np.append(radius,q1r)
	return radius, ismap

def run_efit(sim):

	#sim.StrVar16.set(sim.e16.get())
	sim.StrVar17.set(sim.e17.get())
	sim.StrVar18.set(sim.e18.get())
	sim.StrVar38.set(sim.e38.get())
	sim.StrVar39.set(sim.e39.get())
	sim.StrVar40.set(sim.e40.get())
	sim.StrVar41.set(sim.e41.get())
	sim.StrVar42.set(sim.e42.get())
	sim.StrVar43.set(sim.e43.get())

	sim.StrVar51.set(sim.e51.get())
	sim.StrVar52.set(sim.e52.get())

	return

def draw_efit_boundary(sim,ax1,ax2,ax3,ax4,ax5,ax6,ax7):

	index = sim.MenuVar10.get()
	if index == '--':	return
	index = int(float(index))	

	if sim.efit_first_plot:	sim.efit_plg2 = [];	sim.efit_slg2 = [];	sim.efit_plg3 =[];	sim.efit_slg3 = [];	sim.efit_list3 = []
	if len(sim.efit_list3) == 0.:	sim.efit_plg3 =[];      sim.efit_slg3 = [];

	efitdir = sim.currdir +'/%s/EFIT/'%(sim.MenuVar14.get())
	shot = int(float(sim.e1.get()));	time = int(float(sim.e2.get()));
	sim.efit_index = find_efit_runs(sim)
	save_dir = efitdir + 'RESULT/'
	kfile_dir2 = save_dir + 'k%06i.%06i_kin_%i'%(shot,time,index)
	gfile_dir2 = save_dir + 'g%06i.%06i_kin_%i'%(shot,time,index)
	fit_dir2 = save_dir + 'fitout.dat_%i'%index			

	sim.fig9.canvas.draw_idle()
	eq2 = eqdsk.eqdsk(sim.e3.get())
	gfile0 = sim.currdir +'/%s/EFIT/RESULT/g%06i.%06i_kin_00'%(0,shot,time)
	gfile0m= sim.currdir +'/%s/EFIT/RESULT/g%06i.%06i_kin_0'%(0,shot,time)
	if not os.path.isfile(gfile0): gfile0 = gfile0m
	if (sim.CheckVar24.get() == 1):	eq3 = eqdsk.eqdsk(gfile0)
	eq3m= eqdsk.eqdsk(gfile0m)
	if sim.efit_first_plot:
		eq2.read_eqdsk(sim.e3.get())
		eq3m.read_eqdsk(gfile0m)
		if (sim.CheckVar24.get() == 1):	
			eq3.read_eqdsk(gfile0)
			sim.efit_rmag1 = eq3.rmag
		else:
			sim.efit_rmag1 = 0.

		eq2.make_grid()
		line1, = ax1.plot(eq2.rzbdy[:,0],eq2.rzbdy[:,1])
		ax1.scatter(eq2.rmag,eq2.zmag,marker='+',s=100)
		sim.efit_plg2.append(line1)
		sim.efit_slg2.append('Reference')
		sim.efit_rmax0 = get_midRZ(eq2.rzbdy,eq2.rmag) #max(eq2.rzbdy[:,0])
		sim.efit_rmax1 = get_midRZ(eq3m.rzbdy,eq2.rmag) #max(eq3m.rzbdy[:,0])

		f = open(sim.currdir+'/0/CSOLVE/eq.info','r')
		line = f.readline().split()
		f.close()
		sim.efit_wmhd0 = float(line[1])/1.e3

		f = open(sim.currdir+'/%s/CSOLVE/eq.info'%sim.MenuVar1.get(),'r')
		line = f.readline().split()
		sim.efit_rmag0 = float(line[0])
		if (sim.CheckVar3.get() == 1 and sim.MenuVar7.get().lower() == 'smse'):	
			fname = sim.currdir+'/%s/CSOLVE/Bdiff_%03i/geqdsk_%s'%(sim.MenuVar1.get(),int(round(100.*float(sim.MenuVar6.get()),0)),sim.MenuVar9.get())
			eq2.read_eqdsk(fname)
			sim.efit_rmag2 = sim.efit_rmag0
			sim.efit_rmag0 = eq2.rmag

		sim.efit_first_plot = False
	
	ind = []; rmax = []; chi = [];	rmag = [];	wmhd = [];	cerr = [];	bp = [];

	if sim.efit_first_plot:	ax2.cla();      ax3.cla();      ax4.cla();      ax5.cla();	ax6.cla();	ax7.cla()

	if not sim.MenuVar14.get() in sim.efit_list3:
		sim.efit_list3.append(sim.MenuVar14.get())
		sim.efit_slg3.append('# %s-'%sim.MenuVar14.get())

		fit_dir0 = sim.currdir +'/0/EFIT/RESULT/fitout.dat_0'
		try:	
			read_efit_result(sim,fit_dir0)
			bp_ref = sim.efit_bp
		except:	bp_ref = 0.
		for i in sim.efit_index:		
			gfile_dir2 = save_dir + 'g%06i.%06i_kin_%i'%(shot,time,i)
			fit_dir2 = save_dir + 'fitout.dat_%i'%i
			w_dir2 = save_dir +   'wmhd.dat_%i'%i
			map_dir2 = save_dir + 'map.dat_%i'%i
			eq2.read_eqdsk(gfile_dir2)
			read_efit_result(sim,fit_dir2)
			f = open(w_dir2,'r')
			line = f.readline().split()
			f.close()
			ind.append(i);
			rmax.append(float(line[0]))
			rmag.append(eq2.rmag)
			chi.append(sim.efit_chi)
			cerr.append(sim.efit_conve)
			bp.append(sim.efit_bp)
			wmhd.append(float(line[2]))
		line1, = ax2.plot(ind,rmax,'o--')#,color='red')
		ax3.plot(ind,chi,'o--')#,color='red')
		ax4.plot(ind,rmag,'o--')#,color='red')
		ax5.plot(ind,wmhd,'o--')#,color='red')
		ax6.plot(ind,cerr,'o--')
		ax7.plot(ind,bp,'o--')
		ax2.axhline(y=sim.efit_rmax1,color='magenta',linestyle='--',linewidth = 1.0)
		ax2.text(0.0,sim.efit_rmax1,'MAG-EFIT',color='magenta')
		if (float(sim.MenuVar1.get()) > 0):	
			ax2.axhline(y=sim.efit_rmax0,color='gold',linestyle='--',linewidth = 1.0)
			ax2.text(0.0,sim.efit_rmax0,'INP-REF',color='gold')

		if sim.CheckVar3.get() == 1:	
			if (sim.MenuVar7.get().lower() == 'smse'):
				ax4.axhline(y=sim.efit_rmag0,color='gold',linestyle='--',linewidth = 1.0);	ax4.text(0.0,sim.efit_rmag0,'S-MSE',color='gold')
				ax4.axhline(y=sim.efit_rmag2,color='gray',linestyle='--',linewidth = 1.0);      ax4.text(0.0,sim.efit_rmag2,'INP-REF',color='gray')
			else:
				ax4.axhline(y=sim.efit_rmag0,color='gold',linestyle='--',linewidth = 1.0);      ax4.text(0.0,sim.efit_rmag0,'INP-REF',color='gold')
		else:
			ax4.axhline(y=sim.efit_rmag0,color='gold',linestyle='--',linewidth = 1.0);      ax4.text(0.0,sim.efit_rmag0,'INP-REF',color='gold')
		
		if sim.efit_rmag1 >0.: ax4.axhline(y=sim.efit_rmag1,color='magenta',linestyle='--',linewidth = 1.0);	ax4.text(0.0,sim.efit_rmag1,'MSE-EFIT',color='magenta')
						
		ax5.axhline(y=sim.efit_wmhd0,color='gold',linestyle='--',linewidth = 1.0);	ax5.text(0.0,sim.efit_wmhd0,'MAG-EFIT',color='gold')
		if bp_ref > 0.:	ax7.axhline(y=bp_ref,color='gold',linestyle='--',linewidth = 1.0);	ax7.text(0.0,bp_ref,'MAG-EFIT',color='gold')

		sim.efit_plg3.append(line1)

	ax1.set_title("Plasma boundary")
	ax1.set_xlabel('R [m]')
	ax1.set_ylabel("Z [m]")
	ax2.set_title("$R_{mid,LFS}$")
	ax2.set_xlabel('Iteration [#]')
	ax2.set_ylabel("$R_{mid,LFS}$ [m]")	
	ax3.set_title("$\chi^2$")
	ax3.set_xlabel('Iteration [#]')
	ax3.set_ylabel("$\chi^2$[a.u]")		
	ax4.set_title("$R_{mag}$")
	ax4.set_xlabel('Iteration [#]')
	ax4.set_ylabel("$R_{mag}$ [m]")			
	ax5.set_title("$W_{MHD}$")
	ax5.set_xlabel('Iteration [#]')
	ax5.set_ylabel("$W_{MHD}$ [kJ]")	
	ax7.set_title('$\\beta_{p}$')
	ax7.set_xlabel('Iteration [#]')
	ax7.set_ylabel('$\\beta_{p}$ [a.u]')	
	ax6.set_title('Convergence Error')
	ax6.set_xlabel('Iteration [#]')
	ax6.set_ylabel('Error [a.u]')
	ax6.set_yscale('log')
	
	gfile_dir2 = save_dir + 'g%06i.%06i_kin_%i'%(shot,time,index)
	eq2.read_eqdsk(gfile_dir2)
	line1, = ax1.plot(eq2.rzbdy[:,0],eq2.rzbdy[:,1],linestyle='--')
	ax1.scatter(eq2.rmag,eq2.zmag,marker='+',s=100)
	sim.efit_plg2.append(line1)
	sim.efit_slg2.append('#%s-%i'%(sim.MenuVar14.get(),index))
	ax1.legend(sim.efit_plg2,sim.efit_slg2)
	ax2.legend(sim.efit_plg3,sim.efit_slg3)
	ax3.legend(sim.efit_plg3,sim.efit_slg3)
	ax4.legend(sim.efit_plg3,sim.efit_slg3)
	ax5.legend(sim.efit_plg3,sim.efit_slg3)
	ax6.legend(sim.efit_plg3,sim.efit_slg3)
	ax7.legend(sim.efit_plg3,sim.efit_slg3)

	sim.fig9.tight_layout()
	return

def get_midRZ(rzbdy,rmag):
	
	rr = rzbdy[:,0]
	zz = rzbdy[:,1]
	ind = np.where(rr > rmag)
	rr = rr[ind]
	zz = zz[ind]
	ind = np.argsort(zz)
	rr = rr[ind]
	zz = zz[ind]
	rzf = interp1d(zz,rr,'quadratic')
	
	return rzf(0)

def draw_efit_coils(sim,ax1,ax2,ax3,ax4):

	sim.fig10.canvas.draw_idle()
	ax1.cla();	ax2.cla();	ax3.cla();	ax4.cla();
	index = sim.MenuVar10.get()
	if index == '--':	return
	index = int(float(index))	

	efitdir = sim.currdir +'/%s/EFIT/'%(sim.MenuVar14.get())
	shot = int(float(sim.e1.get()));	time = int(float(sim.e2.get()));
	sim.efit_index = find_efit_runs(sim)
	save_dir = efitdir + 'RESULT/'
	kfile_dir2 = save_dir + 'k%06i.%06i_kin_%i'%(shot,time,index)
	gfile_dir2 = save_dir + 'g%06i.%06i_kin_%i'%(shot,time,index)
	afile_dir2 = save_dir + 'a%06i.%06i_kin_%i'%(shot,time,index)
	fit_dir2 = save_dir + 'fitout.dat_%i'%index
	mse_dir2 = save_dir + 'mse.dat_%i'%index
	mer_dir2 = save_dir + 'mse_er.dat_%i'%index

	read_efit_result(sim,fit_dir2)

	xx1 = np.linspace(1,sim.clen1,sim.clen1);	xx2 = np.linspace(1,sim.clen2,sim.clen2);

	ax1.plot(xx1,sim.chi_mloop,'--',color='lime',linewidth=0.7)
	for i in range(sim.clen1):
		if sim.__dict__['CoilVar%i'%(1+i)].get() == 1:	ax1.scatter(xx1[i],sim.chi_mloop[i],color='lime')
		else:		ax1.scatter(xx1[i],sim.chi_mloop[i],color='gray')
	dely = max(sim.chi_mloop)+0.4
	ax1.scatter(xx1,sim.m_mloop+dely,color='orange',s=4)
	ax1.plot(xx1,sim.c_mloop+dely,linestyle='--',color='blue',linewidth=0.7)

	ax2.plot(xx2,sim.chi_ploop,'--',color='lime')
	for i in range(sim.clen2):
		if sim.__dict__['CoilVar%i'%(151+i)].get() == 1:	ax2.scatter(xx2[i],sim.chi_ploop[i],color='lime')
		else:		ax2.scatter(xx2[i],sim.chi_ploop[i],color='gray')

	ax2.scatter(xx2,sim.m_ploop+2.,color='orange')
	ax2.plot(xx2,sim.c_ploop+2.,linestyle='--',color='blue')

	xx3,yy = read_mse_result(mse_dir2)
	yyt = []
	plegend = []
	slegend = []
	er_corr = False
	color_ind = ['dodgerblue','orange','gold','black']
	if os.path.isfile(mer_dir2):
		er_corr = True
		f = open(mer_dir2,'r')
		line = f.readline()
		for i in range(sim.clen3):
			line = f.readline()
			line2 = line.split()
			yy[i,1] = float(line2[-1])
			yyt.append(line.split('\n')[0])

		er_itern = len(yyt[0].split())-2
		yy2 = np.zeros((sim.clen3,er_itern))
		for i in range(sim.clen3):
			line = yyt[i].split()
			for j in range(er_itern):
				yy2[i,j] = float(line[j+1])
		for i in range(er_itern): 
			plegend.append('')
			slegend.append('Er corr #%i'%(i))
		
	if float(sim.StrVar151.get()) > 0.:   mseb = 'NB1A'; ind = 152
	elif float(sim.StrVar153.get()) > 0.: mseb = 'NB1B'; ind = 154
	else: mseb = 'NB1C'; ind=156

	sshot = int(float(sim.StrVar1.get()))
	if (sshot >= 20756 and sshot < 21759): mseb = 'NB1B'; ind = 154
	if (sshot >= 23057 and sshot < 23134): mseb = 'NB1C'; ind = 156
	if (sshot >= 25283):
		if sim.aa1gam[0] < 0.83: mseb = 'NB1A'; ind = 152
		else: mseb = 'NB1B'; ind = 154
	print('>>> Beam used in MSE', mseb)

	xerr = mse_rad_err[mseb]/1.e3
	xerr2 = np.copy(xerr)
	indp = np.where(xx3[:,0]>1.);
	xx3t = xx3[indp[0],:]; yyt = yy[indp[0],:];
	xf = interp1d(xx3t[:,0],xx3t[:,1])
	for i in range(sim.clen3-1):
		if xx3[i,0] < 1.: continue
		try:	a1=abs(xf(xx3[i,0]+xerr[i]) - xx3[i,1])
		except: a1 =0.
		try:    a2=abs(xf(xx3[i,0]-xerr[i]) - xx3[i,1])
		except: a2 = 0.
		xerr2[i] = max(a1,a2)
	try: xerr2[-1] = abs(xx3[sim.clen3-1,1] - xf(xx3[sim.clen3-1,0]-xerr[-1]))
	except: xerr2[-1] = 0.;
	try: xerr2[0]  = abs(xf(xx3[0,0]+xerr[0]) - xx3[0,1])
	except: xerr2[0] = 0.;
	ax3.plot(xx3t[:,0],yyt[:,0],'--',color='lime')
	for i in range(sim.clen3):
		if xx3[i,0] < 1.: continue
		if sim.__dict__['CoilVar%i'%(201+i)].get() == 1:
			ax3.errorbar(xx3[i,0],yy[i,1], xerr = xerr[i], yerr = 1.5*yy[i,2],fmt='*',markersize='7',c='magenta',ecolor='r',capthick=1)
			if er_corr: 
				for j in range(er_itern): 
					line = ax3.errorbar(xx3[i,0],yy2[i,j], yerr = 0.,fmt='x',markersize='4',c=color_ind[j],ecolor=color_ind[j],capthick=1)
					plegend[j] = line
		else:
			ax3.errorbar(xx3[i,0],yy[i,1], xerr=xerr[i], yerr = 1.5*yy[i,2],fmt='*',markersize='7',c='gray',ecolor='gray',capthick=1)
			if er_corr: 
				for j in range(er_itern): ax3.errorbar(xx3[i,0],yy2[i,j], yerr = 0.,fmt='x',markersize='4',c='gray',ecolor='gray',capthick=1)

	
	ax4.plot(xx3t[:,1],yyt[:,0],'--',color='lime')
	for i in range(sim.clen3):
		if xx3[i,0] < 1.: continue
		if sim.__dict__['CoilVar%i'%(201+i)].get() == 1:
			ax4.errorbar(xx3[i,1],yy[i,1], xerr= xerr2[i],yerr = 1.5*yy[i,2],fmt='*',markersize='7',c='magenta',ecolor='r',capthick=1)
			if er_corr: 
				for j in range(er_itern): ax4.errorbar(xx3[i,1],yy2[i,j], yerr = 0.,fmt='x',markersize='4',c=color_ind[j],ecolor=color_ind[j],capthick=1)
		else:
			ax4.errorbar(xx3[i,1],yy[i,1], xerr =xerr2[i],yerr = 1.5*yy[i,2],fmt='*',markersize='7',c='gray',ecolor='gray',capthick=1)
			if er_corr: 
				for j in range(er_itern): ax4.errorbar(xx3[i,1],yy2[i,j], yerr = 0.,fmt='x',markersize='4',c='gray',ecolor='gray',capthick=1)

	ax1.set_title('Magnetic-Probe [$\chi^2$=%3.1f]'%sum(sim.chi_mloop))
	ax1.set_xlabel('Channel [#]')
	ax1.set_ylabel('[a.u]')
	ax2.set_title('PSI-Loop [$\chi^2$=%3.1f]'%sum(sim.chi_ploop))
	ax2.set_xlabel('Channel [#]')
	ax2.set_ylabel('[a.u]')
	ax3.set_title('MSE-TGAMMA-R')
	ax3.set_xlabel('R [m]')
	ax3.set_ylabel('[a.u]')
	ax4.set_title('MSE-TGAMMA-$\psi_N$')
	ax4.set_xlabel('$\psi_N$ [a.u]')
	ax4.set_ylabel('[a.u]')
	if er_corr: 
		try:
	 		ax3.legend(plegend,slegend)
 			ax4.legend(plegend,slegend)
		except: pass

	ax3.set_ylim([-0.3,0.3])
	ax4.set_ylim([-0.3,0.3])

	sim.fig10.tight_layout()

	return

def draw_efit_mse(sim,ax1,mse_rad_err,mseb):

	sim.fig14.canvas.draw_idle()
	
	if sim.ismse:
		for i in range(sim.clen3):
			if sim.rrrgam[i] <1.0: continue
			if sim.__dict__['CoilVar%i'%(201+i)].get() == 1:
				ax1.errorbar(sim.rrrgam[i],sim.tgamma[i],xerr = mse_rad_err[i], yerr = 1.5*sim.sgamma[i],fmt='*',markersize='7',c='magenta',ecolor='r',capthick=1)
			else:
				ax1.errorbar(sim.rrrgam[i],sim.tgamma[i],xerr = mse_rad_err[i], yerr = 1.5*sim.sgamma[i],fmt='*',markersize='7',c='gray',ecolor='gray',capthick=1)
	else:
		ax1.text(0.5,0.5,'No MSE DATA',color='r')
	
	ax1.set_title('MSE-TGAMMA-R')
	ax1.set_xlabel('R [m]')
	ax1.set_ylabel('[a.u]')
	sim.fig14.tight_layout()
	ax1.set_ylim([-0.3,0.5])
	
	if sim.ismse: ax1.legend(['MSE calibrated with %s'%mseb])

	return

def draw_efit_press(sim,ax1):

	sim.fig12.canvas.draw_idle()
	ax1.cla();
	index = sim.MenuVar10.get()
	if index == '--':	return
	index = int(float(index))	

	efitdir = sim.currdir +'/%s/EFIT/'%(sim.MenuVar14.get())
	shot = int(float(sim.e1.get()));	time = int(float(sim.e2.get()));
	sim.efit_index = find_efit_runs(sim)
	save_dir = efitdir + 'RESULT/'
	gfile_dir2 = save_dir + 'g%06i.%06i_kin_%i'%(shot,time,index)
	pres_dir2 = save_dir + 'pres.dat_%i'%index
	eq2 = eqdsk.eqdsk(gfile_dir2)
	eq2.read_eqdsk(gfile_dir2)
	eq2.make_grid()
	xx, yy = read_pres_result(pres_dir2)

	ax1.plot(eq2.psin,eq2.pres/1.e3,'--',color='green')
	ax1.errorbar(-xx,yy[:,0]/1.e3, yerr = yy[:,1]/1.e3,fmt='o',markersize='4',c='magenta',ecolor='r',capthick=1)
	ax1.set_title('PRESSURE')
	ax1.set_xlabel('$\psi_N$ [a.u]')
	ax1.set_ylabel('Pressure [kPa]')

	sim.fig12.tight_layout();

	return

def draw_efit_j(sim,ax1):

	sim.fig13.canvas.draw_idle()
	ax1.cla();
	index = sim.MenuVar10.get()
	if index == '--':	return
	index = int(float(index))	

	efitdir = sim.currdir +'/%s/EFIT/'%(sim.MenuVar14.get())
	shot = int(float(sim.e1.get()));	time = int(float(sim.e2.get()));
	sim.efit_index = find_efit_runs(sim)
	save_dir = efitdir + 'RESULT/'
	gfile_dir2 = save_dir + 'g%06i.%06i_kin_%i'%(shot,time,index)
	j_dir2 = save_dir + 'jconst.dat_%i'%index
	eq2 = eqdsk.eqdsk(gfile_dir2)
	eq2.read_eqdsk(gfile_dir2)
	eq2.make_grid()
	eq2.construct_volume()
	isjconst = True
	try:	xx, yy = read_j_result(j_dir2)
	except:	
		print('>>> No current constraint')
		isjconst = False

	try:	ax1.plot(sim.efit_psin,abs(sim.efit_jav/eq2.ip*eq2.area),'--',color='green')
	except:	print('>>> Press Load run!')
	if isjconst:
		ax1.scatter(xx,yy,marker='+',s=40,color='magenta')
	ax1.set_title('$<j_{\phi}>_A$')
	ax1.set_xlabel('$\psi_N$ [a.u]')
	ax1.set_ylabel('Normalised current density [a.u]')

	sim.fig13.tight_layout();	
	return

def read_mse_result(filename):

	f = open(filename,'r')
	num = int(f.readline())
	xx = np.zeros((num,2));	yy = np.zeros((num,3))
	for i in range(num):
		line = f.readline().split()
		xx[i,0] = float(line[0])
		xx[i,1] = float(line[1])
		yy[i,0] = float(line[2])
		yy[i,1] = float(line[3])
		yy[i,2] = float(line[4])

	f.close()
	return (xx,yy)

def read_pres_result(filename):
	f = open(filename,'r')
	num = int(f.readline())
	xx = np.zeros(num);	yy = np.zeros((num,2))
	for i in range(num):
		line = f.readline().split()
		xx[i] = float(line[0])
		yy[i,0] = float(line[1])
		yy[i,1] = float(line[2])

	f.close()
	return (xx,yy)	

def read_j_result(filename):
	f = open(filename,'r')
	num = int(f.readline())
	xx = np.zeros(num);	yy = np.zeros(num)
	for i in range(num):
		line = f.readline().split()
		xx[i] = float(line[0])
		yy[i] = float(line[1])

	f.close()
	return (xx,yy)	

def efit_post_process(sim,index,gfile_dir,w_dir,mse_dir,pres_dir,j_dir,map_dir,skip=False):

	sim.MenuVar10.set(str(index))
	print('>>> Post Processing...')
	eq2 = eqdsk.eqdsk(gfile_dir)
	eq2.read_eqdsk(gfile_dir)
	eq2.make_grid()
	eq2.get_flux_contour()
	eq2.make_rho_R_psin()
	eq2.construct_volume()
	eq2.construct_2d_field()

	brf = interp2d(eq2.r,eq2.z,eq2.br)
	bzf = interp2d(eq2.r,eq2.z,eq2.bz)
	btf = interp2d(eq2.r,eq2.z,eq2.bt)		

	len1 =len(sim.rrrgam)
	tgamma = np.zeros(len1)
	psin = np.zeros(len1)
	for i in range(len1):
	
		if sim.CheckVar3.get() == 0:	R = sim.rrrgam[i]+float(sim.StrVar44.get())/1.e2
		else:	R = sim.rrrgam[i]
		Z = sim.zzzgam[i]
		
		sign_ip = np.sign(-eq2.ip)
		br = brf(R,Z) * sign_ip
		bz = bzf(R,Z) * sign_ip
		bt = btf(R,Z) * +1.0
		tgamma[i] = sim.aa1gam[i]*bz/(sim.aa2gam[i]*bt+sim.aa3gam[i]*br+sim.aa4gam[i]*bz) #* np.sign(-eq2.ip) ## SIGN CORR
		psin[i] = eq2.psif(R,Z)*np.sign(R-eq2.rmag)
	
	f = open(w_dir,'w')
	f.write('%f\t%f\t%f\t%f\t%f'%(get_midRZ(eq2.rzbdy,eq2.rmag),eq2.rmag,eq2.wmhd/1.e3,eq2.ip/1.e6,eq2.area))
	f.close()

	f = open(map_dir,'w')
	f.write('%i\n'%len(eq2.prhoR[:,0]))
	for i in range(len(eq2.prhoR[:,0])):	f.write('%9.8f\t%9.8f\t%9.8f\t%9.8f\n'%(eq2.prhoR[i,0],eq2.prhoR[i,1],eq2.prhoR[i,2],eq2.prhoR[i,3]))
	f.close()

	dtgamma = np.copy(sim.dtgamma); dtgamma2 = np.copy(sim.dtgamma);
	if (sim.CheckVar7.get()==0): dtgamma2 = np.zeros(len(sim.dtgamma))
	else: dtgamma = np.zeros(len(sim.dtgamma))

	f = open(mse_dir,'w')
	f.write('%i\n'%len1)
	for i in range(len1):	
		if sim.CheckVar3.get() == 0:	R = sim.rrrgam[i]+float(sim.StrVar44.get())/1.e2
		else:	R = sim.rrrgam[i]
		if sim.CheckVar3.get() == 0:
			#print('>>> SMSE CHECK..')
			f.write('%f\t%f\t%f\t%f\t%f\n'%(R,psin[i],tgamma[i],sim.tgamma[i]+dtgamma[i],sim.sgamma2[i]))
		else:
			#print('>>> EMSE CHECK..',sim.CheckVar3.get())
			f.write('%f\t%f\t%f\t%f\t%f\n'%(R,psin[i],tgamma[i],sim.tgamma2[i]+dtgamma2[i],sim.sgamma2[i]))
	f.close()

	f = open(pres_dir,'w')
	f.write('%i\n'%len(sim.rpress))
	for i in range(len(sim.rpress)):	
		if not skip:	f.write('%f\t%f\t%f\n'%(sim.rpress[i],sim.prest[i],sim.sigpres[i]))
		else:	f.write('%f\t%f\t%f\n'%(sim.rpress[i],0.,0.))
	f.close()	

	if not skip:
		if (sim.CheckVar4.get() == 1):
			jknot = sim.StrVar19.get().split(',')
			jknot = np.array(jknot,dtype='double')	
			jconstf = interp1d(sim.jconst[:,0],sim.jconst[:,1],'cubic')	
			jc = jconstf(jknot)

			f = open(j_dir,'w')
			f.write('%i\n'%len(jknot))
			for i in range(len(jknot)):	
				f.write('%f\t%f\n'%(jknot[i],jc[i]))
			f.close()	

	else:
		if (sim.CheckVar4.get() == 1):
			f = open(j_dir,'w')
			f.write('%i\n'%3)
			for i in range(3):
				f.write('%f\t%f\n'%(1.,0.))
			f.close()

	print('>>> Done...!')

	return

def option_menu_update(var1,value,func,funcvar=None):

	var1.set(value)
	if not funcvar == None:
		func(funcvar)
	else:
		func()
	return

def delete_efit_runs(sim):
	sim.fig8.canvas.draw_idle()
	sim.eax1.cla()
	sim.eax2.cla()
	sim.eax3.cla()
	sim.eax4.cla()
	sim.efit_plg1 = []
	sim.efit_slg1 = []
	sim.efit_plg2 = []
	sim.efit_slg2 = []	
	sim.efit_slg12 = []

	sim.eax1.set_title("P'")
	sim.eax1.set_xlabel('Normalized radius ($\psi_N$)')
	sim.eax1.set_ylabel("P' [kPa/Wb]")
	sim.eax2.set_title("FF'")
	sim.eax2.set_xlabel('Normalized radius ($\psi_N$)')
	sim.eax2.set_ylabel("F' [a.u]")	
	sim.eax3.set_title("<j>")
	sim.eax3.set_xlabel('Normalized radius ($\psi_N$)')
	sim.eax3.set_ylabel("<j> [MA/m2]")	
	sim.eax4.set_title("q")
	sim.eax4.set_xlabel('Normalized radius ($\psi_N$)')
	sim.eax4.set_ylabel("q [a.u]")		

	sim.efit_first_plot = True
	if sim.window_bnd:
		sim.fig9.canvas.draw_idle()
		sim.eax5.cla()
		sim.eax6.cla()
		sim.eax7.cla()
		sim.eax8.cla()
		sim.eax9.cla()
		sim.eax16.cla()
		sim.eax17.cla()
	if sim.window_coil:
		sim.eax10.cla()	
		sim.eax11.cla()
		sim.eax12.cla()
		sim.eax13.cla()	
	if sim.window_pres:
		sim.eax14.cla()	
	if sim.window_j:
		sim.eax15.cla()
	if sim.window_bnd:
		sim.eax5.set_title("Plasma boundary")
		sim.eax5.set_xlabel('R [m]')
		sim.eax5.set_ylabel("Z [m]")
		sim.eax6.set_title("$R_{max}$")
		sim.eax6.set_xlabel('Iteration [#]')
		sim.eax6.set_ylabel("$R_{max}$ LFS [m]")	
		sim.eax7.set_title("$\chi^2$")
		sim.eax7.set_xlabel('Iteration [#]')
		sim.eax7.set_ylabel("$\chi^2 [a.u]")		
		sim.eax8.set_title("$R_{mag}$")
		sim.eax8.set_xlabel('Iteration [#]')
		sim.eax8.set_ylabel("$R_{mag}$ [m]")		
		sim.eax9.set_title("$W_{MHD}$")
		sim.eax9.set_xlabel('Iteration [#]')
		sim.eax9.set_ylabel("$W_{MHD}$ [kJ]")
		sim.eax17.set_title('$\\beta_{p}$')
		sim.eax17.set_xlabel('Iteration [#]')
		sim.eax17.set_ylabel('$\\beta_{p}$ [a.u]')
		sim.eax16.set_title('Convergence Error')
		sim.eax16.set_xlabel('Iteration [#]')
		sim.eax16.set_ylabel('Error [a.u]')	
	if sim.window_riter:
		sim.eax19.cla()
		sim.eax20.cla()
		sim.eax21.cla()
		sim.eax22.cla()

	return

def call_mse_data(dirs,shotn,time,sim=None):
	
	currdir = os.getcwd()
	os.chdir(dirs)
	if shotn==None:	
		shotn=0
		os.system('%s '%(mds_dir))
	else:
		os.system('%s %i %i'%(mds_dir,shotn,time))
	#gtime,ktime,ismse,isefit = get_efit_time_files(shotn,time,dirs)
	gtime = time; ktime = time; ismse = False; isefit = False;
	if os.path.isfile('result.dat'):
		f = open('result.dat','r')
		while True:
			line = f.readline()
			if not line: break
			if line.find('TIME') > -1: 
				gtime = int(line.split()[1]);
				ktime = gtime;
				isefit= True;
				break
		f.close()

	if not sim==None:	sim.CheckVar24.set(sim.trans_vars(ismse,4))
	os.chdir(currdir)

	return gtime,ktime,isefit

def read_mse_data(filename):

	g_file =''
	k_file =''
	te_file =''
	ne_file =''
	ti_file =''
	vt_file =''
	wdia   = 0.
	NBI1AP = 0.
	NBI1AE = 0.
	NBI1BP = 0.
	NBI1BE = 0.
	NBI1CP = 0.
	NBI1CE = 0.
	NBI2AP = 0.
	NBI2AE = 0.
	NBI2BP = 0.
	NBI2BE = 0.
	NBI2CP = 0.
	NBI2CE = 0.
	ECHPW  = 0.

	f = open(filename,'r')
	dat = [];	count = 0;
	while True:
		count = count + 1
		line = f.readline()
		if not line: break
		if len(line.split())>1:	dat.append(line.split()[1])
		else:	
			if count < 7:	dat.append('')
			else:	dat.append(0.)
	f.close()


	return dat


def rescale_kinprof(ifile,ofile,scale=1.):

	with open(ifile,'r') as f1:
		ndat = int(float(f1.readline()))
		params = np.array(f1.readline().split(),dtype='float')
		line = f1.readline().split()
		nline = len(line)
		if nline == 5: dats = np.zeros((ndat,5))
		else: dats = np.zeros((ndat,6))
		dats[0,:] = np.array(line,dtype='float')
		for i in range(1,ndat):
			dats[i,:] = np.array(f1.readline().split(),dtype='float')
	dats[:,2] = dats[:,2] * scale;
	dats[:,4] = dats[:,4] * scale;
	print('>>> Rescale density profile with %f'%scale)
	with open(ofile,'w') as f1:
		f1.write('%i\n'%ndat)
		f1.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(params[0],params[1],params[2],params[3]))
		for i in range(ndat):
			if nline == 5: f1.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(dats[i,0],dats[i,1],dats[i,2],dats[i,3],dats[i,4]))
			else: 	       f1.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(dats[i,0],dats[i,1],dats[i,2],dats[i,3],dats[i,4],dats[i,5]))

	return
