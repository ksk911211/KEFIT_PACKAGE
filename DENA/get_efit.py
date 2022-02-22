import os,sys
import numpy as np
import subprocess
import getpass
from exec_dirs import efit_source_dir, efit_address, shotk, years
import time as ttt

user_id = getpass.getuser()

def get_year(shotn):
	ok = 0
	for i in years:
		for j in range(len(shotk[i]['shot'])):
			if int(shotn) == int(shotk[i]['shot'][j]):	
				ok = 1
				return (ok,i)
	return(ok,i)

def get_efit_list(year,shotn,efit_no1=1,efit_no2=2):
#	command_mag = 'ssh -X -l %s %s ls -l %s%sEXP%06i/g0*'%(user_id,efit_address,efit_source_dir,shotk[year]['mag'],shotn)
	command_mag = 'ls -l %s%sEXP%06i/g0*'%(efit_source_dir,shotk[year][efit_no1],shotn)
#	command_mse = 'ssh -X -l %s %s ls -l %s%sEXP%06i/g0*'%(user_id,efit_address,efit_source_dir,shotk[year]['mse'],shotn)
	command_mse = 'ls -l %s%sEXP%06i/g0*'%(efit_source_dir,shotk[year][efit_no2],shotn)
	
	stat,out_mag = subprocess.getstatusoutput(command_mag)
	stat,out_mse = subprocess.getstatusoutput(command_mse)

	try:
		if out_mag.find('No such file or directory') > -1:	out_mag = None
	except:	out_mag = None
	try:
		if out_mse.find('No such file or directory') > -1:	out_mse = None
	except:
		out_mse = None

	if not out_mag == None:
		out_mag = out_mag.split('\n')
		out_mag = out_mag[1:]
		for i in range(len(out_mag)):	
			try:	out_mag[i] = int(out_mag[i].split('.')[-1])
			except:	
				print(out_mag[i])
				exit()
	if not out_mse == None:
		out_mse = out_mse.split('\n')
		out_mse = out_mse[1:]
		for i in range(len(out_mse)):	out_mse[i] = int(out_mse[i].split('.')[-1])

	return out_mag, out_mse

def get_efit_list2(shotn):
	
	ok, year = get_year(shotn)
	if not ok==1: print('>>> No EFIT for %i'%shot)
	else: print('>>> Year %s'%year)
	efit_list = dict()
	efit_list['isefit'] = dict()
	efit_list['times']  = dict()
	efit_list['dirs']   = dict()
	efit_list['name']   = dict()
	efit_name = ['MAG','MSE','-','MAG+','MSE+']
	for efit_no in range(1,6):
		efit_list[efit_no] = dict()
		command = 'ls -l %s/%s/EXP%06i/g0*'%(efit_source_dir,shotk[year][efit_no],shotn)
		stat,out = subprocess.getstatusoutput(command)
		isefit = True
		try:
			if out.find('No such file or directory') > -1: isefit = False
		except: isefit = False

		if isefit:
			out = out.split('\n')
			out = out[1:]
			for i in range(len(out)):
				try: out[i] = int(out[i].split('.')[-1])
				except: print(out[i]); isefit = False
		else: out = [];
		efit_list['name'][efit_no]  = efit_name[efit_no-1]
		efit_list['isefit'][efit_no]= isefit
		efit_list['times'][efit_no] = np.array(out,dtype='float')
		efit_list['dirs'][efit_no]  = '%s/%s/EXP%06i/'%(efit_source_dir,shotk[year][efit_no],shotn)
		print('>>> EFIT%02i status: %s'%(efit_no,isefit))
	return efit_list

def get_efit_time_files(shotn,time,dirs='',maglist=[],mselist=[],printlist = False,mag_only=False):

	if dirs=='':	dirs='.'
	isefit = True
	ok,year = get_year(shotn)
	if maglist == []:
		out_mag, out_mse = get_efit_list(str(year),shotn)
	else:
		out_mag = maglist; out_mse = mselist

	magnone = out_mag == None
	msenone = out_mse == None
	try: magnone = magnone.any()
	except: pass
	try: msenone = msenone.any()
	except: pass

	if (magnone and msenone):	
		print('>>> There is no efit files!')
		isefit = False	
		return 0,0,False,isefit

	if msenone:
		ismse = False
		print('>>> There is no mse efit!')
	else:	ismse = True	

	out_mag = np.array(out_mag,dtype='int')
	if ismse:	out_mse = np.array(out_mse,dtype='int')
	
	gind = np.argmin(abs(out_mag-time))

	if ismse:	kind = np.argmin(abs(out_mse-time))
	else:		kind = gind

	gtime = out_mag[gind]
	if ismse:	ktime = out_mse[kind]
	else:		ktime = gtime

#	gfile = 'scp %s:%s%sEXP%06i/g%06i.%06i %s'%(efit_address,efit_source_dir,shotk[year]['mag'],shotn,shotn,gtime,dirs)
	gfile = 'cp %s%sEXP%06i/g%06i.%06i %s'%(efit_source_dir,shotk[year][1],shotn,shotn,gtime,dirs)
	if ismse:	
#		kfile = 'scp %s:%s%sEXP%06i/k%06i.%06i %s'%(efit_address,efit_source_dir,shotk[year]['mse'],shotn,shotn,ktime,dirs)
		kfile = 'cp %s%sEXP%06i/k%06i.%06i %s'%(efit_source_dir,shotk[year][2],shotn,shotn,ktime,dirs)

	else:		
#		kfile = 'scp %s:%s%sEXP%06i/k%06i.%06i %s'%(efit_address,efit_source_dir,shotk[year]['mag'],shotn,shotn,ktime,dirs)
		kfile = 'cp %s%sEXP%06i/k%06i.%06i %s'%(efit_source_dir,shotk[year][1],shotn,shotn,ktime,dirs)
	if not mag_only: ttt.sleep(0.4)	
	if not os.path.isfile(dirs+'/g%06i.%06i'%(shotn,gtime)):	os.system(gfile)
	
	if not mag_only:
		ttt.sleep(0.4)
		if not os.path.isfile(dirs+'/k%06i.%06i'%(shotn,ktime)):	os.system(kfile)

	if not printlist: return gtime,ktime,ismse,isefit
	else: return gtime,ktime,ismse,isefit,out_mag,out_mse

