#!/usr/bin/env python

from __future__ import print_function
import os

rundir=os.environ['HOME']+'/NUBEAM_RUNDIR'
debug=False
sge_command='qsub -sync y -cwd -V -q "ukstar-node*" -b y'
pbs_command='qsub -Wblock=true -V -q compute -N NUBEAM'
parallel_env='make'
#geqdsk='g016498.002500'
#plasma_state='s016498.002500'
#zeff=2.0
#dbeam=0.0
#nbeam=3
#beam_power=[1697000,1087000,1751000]
#beam_energy=[100.14,66.815,84.064]
#hostfile='host11'
#nproc=18

# Hard coded Params
xplasma_data_list=['rho','vol','area']
nubeam_data_list=['pbe','pbi','curbeam','tqbi','tqbe','tqbth','tqbjxb','sbedep','nbeami','eperp_beami','epll_beami']
ffulla=[0.60,0.40,0.40,0.60,0.60,0.60]
fhalfa=[0.22,0.35,0.35,0.22,0.22,0.22]
###################

import numpy as np
import datetime
import shutil
from os.path import basename
from os import system
from os.path import isfile
from scipy.io.netcdf import netcdf_file
from os import environ
from shutil import move
import matplotlib.pyplot as plt
import glob
import errno
from exec_dirs import mpirun, plasma_state_test_exec, nubeam_comp_exec, adasdir, preactdir

def _jobs_execute(*args,**kwargs):
	if system('which qsub | grep sge > /dev/null')==0:
		return sge_execute(*args,**kwargs)
	if system('which qsub | grep pbs > /dev/null')==0:
		return pbs_execute(*args,**kwargs)
	return

def jobs_execute(executable,nproc,tag=None):
	if tag is not None:
		log='-e log.'+basename(executable)+'_'+tag+'_err -o log.'+basename(executable)+'_'+tag
	else:
		log='-e log.'+basename(executable)+'_err -o log.'+basename(executable)

	if system('which qsub | grep sge > /dev/null')==0:
		comm=sge_command+' '+log+' -pe make %d'%nproc
	if system('which qsub | grep pbs > /dev/null')==0:
		comm=pbs_command+' '+log+' -l select=1:ncpus=%d:mpiprocs=%d '%(nproc,nproc)+'-v NUBEAM_WORKPATH=`pwd` --'

	command=comm+' '+mpirun+' -np %d'%(nproc)+' '+executable
	print('Running',tag)
#	print(command)
	system(command)

def dec(func,*args,**kwargs):
	def _dec(*args,**kwargs):
		if debug:
			print('call '+func.__name__)
			print(args,kwargs)
		return func(*args,**kwargs)
	return _dec

def mkdir(dir):
	try:
		os.makedirs(dir)
	except OSError as ose:
		if ose.errno==errno.EEXIST and os.path.isdir(dir):
			pass
		else:
			raise

def cat(filename):
	with open(filename,'r') as f:
		return f.read()

def execute(executable,tag=None,quiet=True):
	if tag is None:
		command=executable+' &> log.'+basename(executable) #+' 2> log.'+basename(executable)+'_err'
	else:
		command=executable+' &> log.'+basename(executable)+'_'+tag #+' 2> log.'+basename(executable)+'_'+tag+'_err'
	if not quiet:
		print(command)
	system(command)

def sge_execute(executable,nproc,tag=None):
	if tag is not None:
		sge_command=qsub_command+' -pe %s %d'%(parallel_env,nproc)+' -e log.'+basename(executable)+'_'+tag+'_err -o log.'+basename(executable)+'_'+tag
	else:
		sge_command=qsub_command+' -pe %s %d'%(parallel_env,nproc)+' -e log.'+basename(executable)+'_err -o log.'+basename(executable)
	command=sge_command+' '+mpirun+' '+executable
	print(command)
	system(command)

def pbs_execute(executable,nproc,tag=None):
	from batch_run import run_batch_script_mpi

	if tag is not None:
		log_e='log.'+basename(executable)+'_'+tag+'_err'
		log_o='log.'+basename(executable)+'_'+tag
	else:
		log_e='log.'+basename(executable)+'_err'
		log_o='log.'+basename(executable)

#	command=mpirun+' -np %d %s'%(nproc,executable)
	command=executable
	run_batch_script_mpi('NUBEAM',log_e,log_o,1,nproc,nproc,command)
	#input()

def listize(inst):
	if not isinstance(inst,list):
		if isinstance(inst,np.ndarray):
			return list(inst)
		else:
			return [inst]
	return inst

def string6(integer):
	if integer<10:
		return '00000%d'%integer
	elif integer<100:
		return '0000%d'%integer
	elif integer<1000:
		return '000%d'%integer
	elif integer<10000:
		return '00%d'%integer
	elif integer<100000:
		return '0%d'%integer
	elif integer<1000000:
		return '%d'%integer
	else:
		return '000000'

def gaussian_filter(x,y,xx,sig):
	yy=np.zeros(len(xx))
	for i in range(len(xx)):
		yy[i]=0.
		weight=np.zeros(len(x))
		for j in range(len(x)):
			weight[j]=np.exp(-0.5*((x[j]-xx[i])/sig)**2)/sig
		weight=weight/np.sum(weight)
		for j in range(len(x)):
			yy[i]=yy[i]+weight[j]*y[j]
	return yy

		
def check_existence(name):
	if not isfile(name):
		raiserror(name+' doesn\'t exist')

def get_nc(file,variable):
#	from netCDF4 import Dataset
#	f=Dataset(file,'r')
#	data=f.variables[variable].getValue()
#	f.close()
#	return data
	with netcdf_file(file,'r',mmap=False) as f:
		return f.variables[variable].data

def set_nc(file,variable,value):
#	from netCDF4 import Dataset
#	f=Dataset(file,'r+')
#	f.variables[variable].assignValue(value)
#	f.close()
	with netcdf_file(file,'a') as f:
		f.variables[variable].assignValue(value)

def raiserror(log):
	print(' !!!!!!!!!!!!!!!!!! '+log+' !!!!!!!!!!!!!!!!!! ')
	exit()

def stage_input_file(files):
	for file in listize(files):
		check_existence(file)
		try:
			shutil.copy(file,'.')
		except:
			pass

@dec
def make_inputf(zeff,dbeam,nbeam,beam_power,beam_energy,geqdsk,NE,TE,TI,VT,plasma_state,mdescr,sconfig):
	with open('inputf','w') as f:
		f.write('%s\n'%basename(geqdsk))
		f.write('%s\n'%basename(mdescr))
		f.write('%s\n'%basename(sconfig))
		f.write('%s\n'%plasma_state)
		f.write('CKSTAR\n')
		f.write('2\n')
		f.write('%f\n'%zeff)
		f.write('%s\n'%basename(NE))
		f.write('%s\n'%basename(TE))
		f.write('%s\n'%basename(TI))
		f.write('%s\n'%basename(VT))
		f.write('%f\n'%dbeam)
		f.write('3\n')
		f.write('81. 91. 101.\n')
		f.write('1. 1. 1.\n')
		f.write('65\n')
		f.write('%d\n'%nbeam)
		for beam_index in range(nbeam):
			f.write('%f %f\n'%(listize(beam_power)[beam_index],listize(beam_energy)[beam_index]))
		for beam_index in range(nbeam):
			f.write('%f %f\n'%(ffulla[beam_index],fhalfa[beam_index]))
		f.write('1.e18\n')
		f.write('2\n')
		f.write('0.015\n0.0045\n')

@dec
def make_plasma_state(zeff,dbeam,nbeam,beam_power,beam_energy,geqdsk,NE,TE,TI,VT,plasma_state,
					mdescr='/home/ksk911211/NTCC/mdescr_18483.dat',sconfig='/home/ksk911211/NTCC/sconfig_18483.dat'):
	make_inputf(zeff=zeff,dbeam=dbeam,nbeam=nbeam,beam_power=beam_power,beam_energy=beam_energy,geqdsk=geqdsk,NE=NE,TE=TE,TI=TI,VT=VT,plasma_state=plasma_state,mdescr=mdescr,sconfig=sconfig)
	stage_input_file([geqdsk,NE,TE,TI,VT,mdescr,sconfig])
	execute(plasma_state_test_exec)

@dec
def set_nubeam_env_init():
	environ['NUBEAM_ACTION']='INIT'
	environ['NUBEAM_POSTPROC']=''
	environ['ADASDIR']=adasdir
	environ['PREACTDIR']=preactdir

@dec
def set_nubeam_env_step(run_step,time_step):
	environ['NUBEAM_ACTION']='STEP'
	environ['NUBEAM_POSTPROC']='SUMMARY_TEST'
	environ['ADASDIR']=adasdir
	environ['PREACTDIR']=preactdir
	environ['NUBEAM_REPEAT_COUNT']='%dx%f'%(run_step,time_step)

@dec
def set_nubeam_env_avg(time_avg):
	environ['NUBEAM_ACTION']='STEP'
	environ['NUBEAM_POSTPROC']='SUMMARY_TEST'
	environ['ADASDIR']=adasdir
	environ['PREACTDIR']=preactdir
	environ['NUBEAM_REPEAT_COUNT']='%dx%f'%(1,time_avg)

@dec
def make_nubeam_files(plasma_state,init_file,step_file):
	with open('nubeam_init_files.dat','w') as f:
		f.write('&NUBEAM_FILES\n')
		f.write(' input_plasma_state = "%s"\n'%plasma_state)
		f.write(' plasma_state_update = "state_changes.cdf"\n')
		f.write(' init_namelist = "%s"\n'%basename(init_file))
		f.write(' /')

	with open('nubeam_step_files.dat','w') as f:
		f.write('&NUBEAM_FILES\n')
		f.write(' input_plasma_state = "%s"\n'%plasma_state)
		f.write(' plasma_state_update = "state_changes.cdf"\n')
		f.write(' step_namelist = "%s"\n'%basename(step_file))
		f.write(' /')

@dec
def clean_directory():
# Erase plasma state test input files
	os.system('rm inputf cur_state.cdf *.txt &> /dev/null')
# Erase log files
	os.system('rm log.* *.log *.log~ *.msgs &> /dev/null')
# Erase nubeam related output files
	os.system('rm *.cdf *.cdf~ &> /dev/null')

	return 

@dec
def nubeam_init():
	
	set_nubeam_env_init()
#	if hostfile is None:
#		execute('%s -np %d %s'%(mpirun,nproc,nubeam_comp_exec),'init')
#	else:
#		check_existence(hostfile)
#		execute('%s -hostfile %s -np %d %s'%(mpirun,hostfile,nproc,nubeam_comp_exec),'init')

	execute(nubeam_comp_exec,'init')

@dec
def nubeam_step(run_step,time_step,plasma_state,nproc=1):
	if run_step==0:
		return

	set_nubeam_env_step(run_step,time_step)
#	if hostfile is None:
#		execute('%s -np %d %s'%(mpirun,nproc,nubeam_comp_exec),'step')
#	else:
#		check_existence(hostfile)
#		execute('%s -hostfile %s -np %d %s'%(mpirun,hostfile,nproc,nubeam_comp_exec),'step')
	jobs_execute(nubeam_comp_exec,nproc,'step')
	set_nc(plasma_state,'t0',get_nc(plasma_state,'t0')+run_step*time_step)
	check_existence('state_changes.cdf')
	move('state_changes.cdf','state_changes_0.cdf')

@dec
def nubeam_avg(run_avg,time_avg,plasma_state,nproc=1):
	if run_avg==0:
		return

	set_nubeam_env_avg(time_avg)
	for step in range(run_avg):
#		if hostfile is None:
#			execute('%s -np %d %s'%(mpirun,nproc,nubeam_comp_exec),'avg%d'%(step+1))
#		else:
#			check_existence(hostfile)
#			execute('%s -hostfile %s -np %d %s'%(mpirun,hostfile,nproc,nubeam_comp_exec),'avg%d'%(step+1))
		jobs_execute(nubeam_comp_exec,nproc,'avg%d'%(step+1))
		set_nc(plasma_state,'t0',get_nc(plasma_state,'t0')+time_avg)
		check_existence('state_changes.cdf')
		move('state_changes.cdf','state_changes_%d.cdf'%(step+1))

def get_avg_data(variable,run_avg):
	try:
#		print('Read data from state_changes_0~%d.cdf'%(run_avg))
		return np.mean([get_nc('state_changes_%d.cdf'%step,variable) for step in range(run_avg+1)],axis=0)
	except IOError:
#		print('Error occurred, Read data from state_changes_1~%d.cdf'%(run_avg))
		return np.mean([get_nc('state_changes_%d.cdf'%(step+1),variable) for step in range(run_avg)],axis=0)

@dec
def get_nubeam_data(run_avg):
	nubeam_data=dict()
	for data in nubeam_data_list:
		nubeam_data[data]=get_avg_data(data,run_avg)
	return nubeam_data

@dec
def get_xplasma_data(plasma_state):
	xplasma_data=dict()
	for data in xplasma_data_list:
		xplasma_data[data]=get_nc(plasma_state,data)
	return xplasma_data

def postprocess_nubeam_data(nubeam_data,xplasma_data):
	pp_nubeam_data=dict()
	rho=(xplasma_data['rho'][1:]+xplasma_data['rho'][:-1])*0.5
	dvol=xplasma_data['vol'][1:]-xplasma_data['vol'][:-1]
	darea=xplasma_data['area'][1:]-xplasma_data['area'][:-1]
	pp_nubeam_data['pbe']=nubeam_data['pbe']/dvol
	pp_nubeam_data['pbi']=nubeam_data['pbi']/dvol
	pp_nubeam_data['petot']=np.sum(nubeam_data['pbe'])
	pp_nubeam_data['pitot']=np.sum(nubeam_data['pbi'])
	pp_nubeam_data['ptot']=pp_nubeam_data['petot']+pp_nubeam_data['pitot']
	pp_nubeam_data['curbeam']=nubeam_data['curbeam']/darea
	pp_nubeam_data['nbeami']=nubeam_data['nbeami'][0]*1e-19
	pp_nubeam_data['pperp']=nubeam_data['nbeami'][0]*nubeam_data['eperp_beami'][0]*1e-19
	pp_nubeam_data['ppll']=nubeam_data['nbeami'][0]*nubeam_data['epll_beami'][0]*2e-19
	pp_nubeam_data['pfast']=(2.*pp_nubeam_data['pperp']+pp_nubeam_data['ppll'])/3.
	pp_nubeam_data['wfast']=1.5*1602.*np.sum(pp_nubeam_data['pfast']*dvol)
	pp_nubeam_data['tqb']=(nubeam_data['tqbi']+nubeam_data['tqbe']+nubeam_data['tqbth']+nubeam_data['tqbjxb'])/dvol
	pp_nubeam_data['tqbt']=np.sum(pp_nubeam_data['tqb']*dvol)
	pp_nubeam_data['sbedep']=nubeam_data['sbedep']/dvol*1e-19
	pp_nubeam_data['rho']=rho

	for key in pp_nubeam_data.keys():
		if key in ['rho','petot','pitot','ptot','wfast','tqbt']:
			continue
		pp_nubeam_data[key]=gaussian_filter(pp_nubeam_data['rho'],pp_nubeam_data[key],np.linspace(0,1,101),0.05)

	pp_nubeam_data['rho']=np.linspace(0,1,101)

	return pp_nubeam_data

def plot_dict(dic):
	if 'rho' not in dic.keys():
		print('no rho in dict, return')
		return
	fig=plt.figure()

	nop=len(dic.keys())-1
	nop2=nop/2 if nop%2==0 else nop/2+1
	imax=1
	for key in dic.keys():
		if key is not 'rho':
			ax=plt.subplot(2,nop2,imax)
			ax.plot(dic['rho'],dic[key],label=key)
			ax.legend()
			imax+=1
	plt.tight_layout()
	plt.show()

@dec
def run_nubeam_from_file(file):
	from InputObject import InputObject
	a=InputObject(file)
	run_nubeam(**a)

def appropriate_dir():
	mkdir(rundir)
	index=1
	savedir=rundir+'/%d'%index
	runs=glob.glob(rundir+'/*')
	while (savedir in runs) and (os.path.isfile(savedir+'/success')):
		index+=1
		savedir=rundir+'/%d'%index
	return savedir

def describerun(simname,run_step,run_avg,time_step,time_avg,zeff,dbeam,nbeam,beam_power,beam_energy,geqdsk,NE,TE,TI,VT,nproc,shot,time,plasma_state,mdescr,sconfig,init_file,step_file):
#	with open('runoption','w') as f:
	runoption=''
	runoption+='Simulation name: %s\n'%simname
	runoption+='Run steps: %d\n'%run_step
	runoption+='Run average steps: %d\n'%run_avg
	runoption+='Time steps: %f\n'%time_step
	runoption+='Time average steps: %f\n'%time_avg
	runoption+='Zeff: %f\n'%zeff
	runoption+='D_beam: %f\n'%dbeam
	runoption+='N_beam: %d\n'%nbeam
	for i in range(nbeam):
		runoption+='Beam power #%d: %f\n'%(i+1,beam_power[i])
		runoption+='Beam energy #%d: %f\n'%(i+1,beam_energy[i])
	runoption+='Number of processors: %d\n'%nproc
	runoption+='Shot number: %d\n'%shot
	runoption+='Time: %f\n'%time
	runoption+='Name of plasma state file: %s\n'%plasma_state
	runoption+='NE profile: \n'+cat(NE)+'\n'
	runoption+='TE profile: \n'+cat(TE)+'\n'
	runoption+='TI profile: \n'+cat(TI)+'\n'
	runoption+='VT profile: \n'+cat(VT)+'\n'
	runoption+='GEQDSK file contains: \n'+cat(geqdsk)+'\n'
	runoption+='MDESCR file contains: \n'+cat(mdescr)+'\n'
	runoption+='SCONFIG file contains: \n'+cat(sconfig)+'\n'
	runoption+='NUBEAM INIT file contains: \n'+cat(init_file)+'\n'
	runoption+='NUBEAM STEP file contains: \n'+cat(step_file)+'\n'
	return runoption

@dec
def run_nubeam(simname,run_step,run_avg,time_step,time_avg,zeff,dbeam,nbeam,beam_power,beam_energy,geqdsk,NE='NE.dat',TE='TE.dat',TI='TI.dat',VT='VT.dat',nproc=1,shot=0,time=0,plasma_state=None,
				mdescr='/home/ksk911211/NTCC/mdescr_18483.dat',sconfig='/home/ksk911211/NTCC/sconfig_18483.dat',
				init_file='/home/ksk911211/NTCC/nubeam_init.dat',step_file='/home/ksk911211/NTCC/nubeam_step.dat',force=False):
	cwd=os.getcwd()
	walltime=datetime.datetime.now()
	print('Start at: ',walltime)
	clean_directory()
#	runoption=describerun(simname,run_step,run_avg,time_step,time_avg,zeff,dbeam,nbeam,beam_power,beam_energy,geqdsk,NE,TE,TI,VT,nproc,shot,time,plasma_state,mdescr,sconfig,init_file,step_file)
#	rundir_case=find_run(runoption)
#	if rundir_case and (not force):
#		print('Run already exists, path: '+rundir_case)
#		os.chdir(rundir_case)
#	else:
#		savedir=appropriate_dir()
#		os.makedirs(savedir)
#		mkdir(savedir)
		
#		os.chdir(savedir)
#		with open('runoption','w') as f:
#			f.write(runoption)

	if plasma_state is None:
		plasma_state='s%s.%s'%(string6(shot),string6(time))
	make_plasma_state(zeff=zeff,dbeam=dbeam,nbeam=nbeam,beam_power=beam_power,beam_energy=beam_energy,geqdsk=geqdsk,NE=NE,TE=TE,TI=TI,VT=VT,plasma_state=plasma_state,mdescr=mdescr,sconfig=sconfig)
	make_nubeam_files(plasma_state,init_file,step_file)
	stage_input_file([init_file,step_file])

	nubeam_init()
	nubeam_step(run_step=run_step,time_step=time_step,plasma_state=plasma_state,nproc=nproc)
	nubeam_avg(run_avg=run_avg,time_avg=time_avg,plasma_state=plasma_state,nproc=nproc)

#	for i in range(run_avg+1):
#		shutil.copy('state_changes_%d.cdf'%i,cwd)
#	shutil.copy(plasma_state,cwd)

#	with open('success','w') as f:
#		f.write('%s\n'%datetime.datetime.now())

#	os.chdir(cwd)

	nubeam_data=get_nubeam_data(run_avg=run_avg)
	xplasma_data=get_xplasma_data(plasma_state=plasma_state)

#	dir_index=0
#	while True:
#		savedir=rundir+'/'+os.environ['USER']+'_'+simname+'_%d'%dir_index
#		if os.path.isdir(savedir):
#			dir_index+=1
#		else:
#			break

#	os.mkdir(savedir)
#	for step in range(run_avg+1):
#		shutil.copy('state_changes_%d.cdf'%step,savedir)
#	for file in [plasma_state,geqdsk,init_file,step_file]:
#		shutil.copy(file,savedir)
#	with open(savedir+'/nubeam_config','w') as f:
#		f.write('run_step: %d\nrun_avg: %d\ntime_step: %f\ntime_avg: %f\n'%(run_step,run_avg,time_step,time_avg))

	print('End of NUBEAM run, wall time: ',datetime.datetime.now()-walltime)

	return postprocess_nubeam_data(nubeam_data,xplasma_data)

def find_run(runoption):
	runs=glob.glob(rundir+'/*')
	for run in runs:
		if (cat('%s/runoption'%run)==runoption) and (os.path.isfile('%s/success'%run)):
			return run
	return False
			

if __name__=='__main__':
	while True:
		mode=raw_input('0: Run NUBEAM, 1: Run NUBEAM from file, 2: Plot result ')
		if mode=='0':
			simname=raw_input('Simname? ')
			run_step=int(raw_input('Run step? '))
			run_avg=int(raw_input('Run avg? '))
			time_step=float(raw_input('Time step? '))
			time_avg=float(raw_input('Time average? '))
			zeff=float(raw_input('Zeff? '))
			dbeam=float(raw_input('Dbeam? '))
			nbeam=int(raw_input('Nbeam? '))
			beam_power=[]; beam_energy=[];
			for i in range(nbeam):
				beam_power.append(float(raw_input('Beam power #%d? [W] '%(i+1))))
				beam_energy.append(float(raw_input('Beam energy #%d? [keV] '%(i+1))))
			geqdsk=raw_input('GEQDSK file? ')
			nproc=int(raw_input('Number of processors? '))
			plasma_state=raw_input('Plasma state file? ')
			mdescr=raw_input('MDESCR? (/home/ksk911211/NTCC/mdescr_18483.dat) ') or '/home/ksk911211/NTCC/mdescr_18483.dat'
			sconfig=raw_input('SCONFIG? (/home/ksk911211/NTCC/sconfig_18483.dat) ') or '/home/ksk911211/NTCC/sconfig_18483.dat'
			init_file=raw_input('NUBEAM init file? (/home/ksk911211/NTCC/nubeam_init.dat) ') or '/home/ksk911211/NTCC/nubeam_init.dat'
			step_file=raw_input('NUBEAM step file? (/home/ksk911211/NTCC/nubeam_step.dat) ') or '/home/ksk911211/NTCC/nubeam_step.dat'
#			pp_data=run_nubeam('ptl5000_whgt1',run_step=20,run_avg=10,time_step=0.01,time_avg=0.01,zeff=2,dbeam=0,nbeam=3,
#								beam_power=[1697000,1087000,1751000],beam_energy=[100.14,66.815,84.064],geqdsk='g016498.002500',
#								nproc=40,hostfile='host11',plasma_state='s016498.002500')
			pp_data=run_nubeam(simname,run_step=run_step,run_avg=run_avg,time_step=time_step,time_avg=time_avg,zeff=zeff,dbeam=dbeam,
								nbeam=nbeam,beam_power=beam_power,beam_energy=beam_energy,geqdsk=geqdsk,nproc=nproc,plasma_state=plasma_state,
								mdescr=mdescr,sconfig=sconfig,init_file=init_file,step_file=step_file)
		elif mode=='1':
			file=raw_input('NUBEAM input file? ')
			pp_data=run_nubeam_from_file(file)
		elif mode=='2':
			plasma_state=raw_input('Name of plasma state cdf file? ')
			run_avg=int(raw_input('Run averages? '))
			nubeam_data=get_nubeam_data(run_avg=run_avg)
			xplasma_data=get_xplasma_data(plasma_state=plasma_state)
			pp_data=postprocess_nubeam_data(nubeam_data,xplasma_data)
			plot_dict(pp_data)



