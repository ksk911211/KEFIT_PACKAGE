#!/usr/local/anaconda3/bin/python3
import os,sys
from gefit_tool import read_kfile, load_kfile, get_kfile_data, make_kfile_str
import numpy as np
import matplotlib.pyplot as plt
from MDS import mds
from shutil import move
from exec_dirs import shotk, years, efit_source_dir,mpraw_file,mpraw_dat,efit_dir
import pickle

class mp_coil:
	def __init__(self):

		self.mpraw_file = mpraw_file
		self.mpraw_dat  = mpraw_dat

		self.mp_coils = dict()
		self.coils_ip = dict()
		self.mp_coils_comp = dict()

		self.col_index = ['j','n','b','f']
		self.row_index = ['t','m','b']

		for ind in range(1,53): 
			self.mp_coils['MP4P%02iR'%ind] = 0.
			self.mp_coils['MP4P%02iZ'%ind] = 0.
			self.mp_coils_comp['MP4P%02iR'%ind] = dict()
			self.mp_coils_comp['MP4P%02iZ'%ind] = dict()
			for col in self.col_index:
				for row in self.row_index:
					self.mp_coils_comp['MP4P%02iR'%ind][row+col] = 0.
					self.mp_coils_comp['MP4P%02iZ'%ind][row+col] = 0.

		for ind in range(1,46):
			self.mp_coils['FL%02i'%ind]    = 0.
			self.mp_coils_comp['FL%02i'%ind] = dict()
			for col in self.col_index:
				for row in self.row_index:
					self.mp_coils_comp['FL%02i'%ind][row+col] = 0.

		for col in self.col_index:
			for row in self.row_index:
				self.coils_ip[row+col]=0.
		return

def load_coil(x):
	g = mds('kstar',int(x.shot))
	for col in x.col_index:
		for row in x.row_index:
			x.coils_ip[row+col] = g.get('\\pcrmp%s%suli'%(row,col))
	return

def read_rmpcoil(x):
	
	if os.path.isfile(x.mpraw_dat):
		print('>>> Read MP pickle...')
		f = open(x.mpraw_dat,'rb')
		x.mp_coils_comp = pickle.load(f)
		f.close()
		return
	print('>>> Read MP raw...')
	f = open(x.mpraw_file,'r')
	line = f.readline()
	for i in range(45+42*2):
		line = f.readline().split()
		flag = line[0].split("'")[-1]
		for row in range(3):
			for  col in range(4):
				rowf = x.row_index[row]
				colf = x.col_index[col]
				x.mp_coils_comp[flag][rowf+colf] = float(line[col+4*row+2])
	f.close()
	f = open(x.mpraw_dat,'wb')
	pickle.dump(x.mp_coils_comp,f)
	f.close()
	return

def compensate_coil(x):

	ind_sum   = 0
	for row in x.row_index:
		for col in x.col_index:
			flag = row + col
			coil_ip = x.coils_ip[flag][1][x.tind]
			for i in range(45): 
				coil_flag     = 'FL%02i'%(i+1)
				x.coils[i]    = x.coils[i]   - coil_ip * x.mp_coils_comp[coil_flag][flag]
			for i in range(41):
				i2 = i;
				if i>8: i2 = i + 1
				coil_flagr    = 'MP4P%02iR'%(i2+1)
				coil_flagz    = 'MP4P%02iZ'%(i2+1)
				x.expmp2[i]   = x.expmp2[i]    - coil_ip * x.mp_coils_comp[coil_flagz][flag]
				x.expmp2[i+41]= x.expmp2[i+41] - coil_ip * x.mp_coils_comp[coil_flagr][flag]
				
			ind_sum = ind_sum + 1
	return

def check_kfile(x):

	if not os.path.isfile(x.kfile): print('>>> No available k-file...'); exit()
	read_kfile(x,x.kfile,True)
	load_kfile(x,False,True)
	for item in x.kfile_in1:
		if item.find('corr_3d')>-1: print('>>> 3D corrected version... exit'); exit()
	print('>>> 3D coil compensation...')
	return

def get_coils(x):
	x.shot  = get_kfile_data(x.kfile_in1,'ISHOT')
	load_coil(x)
	read_rmpcoil(x)

def adjust_kfile(x):
	x.time   = get_kfile_data(x.kfile_in1,'ITIME',True)
	x.expmp2 = get_kfile_data(x.kfile_in1,'EXPMP2')
	x.coils  = get_kfile_data(x.kfile_in1,'COILS')
	x.tind = np.argmin(abs(x.coils_ip['tj'][0]*1000.-x.time))
	compensate_coil(x)
	line=make_kfile_str(x.expmp2,'EXPMP2')
	
	x.kfile_in1.append(line+'\n')
	line=make_kfile_str(x.coils,'COILS')
	x.kfile_in1.append(line)
	x.kfile_in1.append('\n !corr_3d\n')
	return

def write_efit(x,filename):

	f = open(filename,'w')
	for item in x.kfile_in1: f.write(item)
	f.write('/\n'); f.write('&INWANT\n')
	for item in x.kfile_in2: f.write(item)
	if len(x.kfile_in3)>1:
		f.write('/\n'); f.write('/\n'); f.write('&INS\n')
		for item in x.kfile_in3: f.write(item)
	f.write('/\n')		
	f.close()
	return

x = mp_coil()

if len(sys.argv)==2:
	x.kfile = sys.argv[1]
	check_kfile(x)
	get_coils(x)
	adjust_kfile(x)
	move(x.kfile,x.kfile+'_old')
	write_efit(x,x.kfile)
	exit()

if len(sys.argv)==2: print('>>> No Input...'); exit()
shotn = int(sys.argv[1]); itime = int(sys.argv[2]); idel = int(sys.argv[3]);
iscoil = False
for yy in years:
	if shotn in shotk[yy]['shot']: break
dirs = efit_source_dir+shotk[yy][1]+'/EXP%06i/'%shotn
print(dirs)
if os.path.isdir('EXP%06i'%shotn): exit();
os.mkdir('EXP%06i'%shotn);
os.chdir('EXP%06i'%shotn);
while itime < 1.e6:
	x.kfile = dirs+'k%06i.%06i'%(shotn,itime)
	if not os.path.isfile(x.kfile): itime = itime + idel; continue;
	read_kfile(x,x.kfile,True)
	load_kfile(x,False,True)
	if not iscoil: get_coils(x); iscoil = True;
	adjust_kfile(x)
	print('%06i ms'%itime)
	write_efit(x,'kfile_run')
	command = efit_dir + '/efit65p'
	os.system(command)
	if not os.path.isfile('g%06i.%06i'%(x.shot,itime)): itime = itime + idel; continue;
#	move('g%06i.%06i'%(x.shot,itime),'OUT/g%06i.%06i'%(x.shot,itime))
#	move('m%06i.%06i'%(x.shot,itime),'OUT/m%06i.%06i'%(x.shot,itime))
#	move('a%06i.%06i'%(x.shot,itime),'OUT/a%06i.%06i'%(x.shot,itime))
	itime = itime + idel
