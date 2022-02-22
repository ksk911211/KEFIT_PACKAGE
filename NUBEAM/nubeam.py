import os,sys
import numpy as np
from nubeamtool import *
from exec_dirs import Mfile,Sfile,Ifile,stepfile

try:
	f = open('nubeam_opt','r')
except:
	print('no nubeam option file')
	exit()
	
Zeff = 2
Nbeam = 2
Power = [1400000,1363000]
Energy = [90.155,76.312]
Bdiff = 0.
Shot = 10000
Time = 5000
Nproc = 28
Runstep = 20
Runavg = 0
Rundt = 0.01
Avgdt = 0.01

Mfile2 = Mfile
Sfile2 = Sfile
Ifile2 = Ifile
	
while True:

	line = f.readline()
	if not line: break
	
	line2 = line.split('=')
	
	if (line2[0].lower().find('zeff') > -1):
		Zeff = float(line2[1].split()[0])
	if (line2[0].lower().find('nbeam') > -1):
		Nbeam = int(float(line2[1].split()[0]))
	if (line2[0].lower().find('bpower') > -1):
		Power = line2[1].split(',')
	if (line2[0].lower().find('benergy') > -1):
		Energy = line2[1].split(',')
	if (line2[0].lower().find('diffusivity') > -1):
		Bdiff = float(line2[1].split()[0])
	if (line2[0].lower().find('shot') > -1):
		Shot = int(line2[1].split()[0])
	if (line2[0].lower().find('time') > -1):
		Time = int(line2[1].split()[0])	
	if (line2[0].lower().find('eqdsk') > -1):
		Eqdsk = line2[1].split()[0]
	if (line2[0].lower().find('nproc') > -1):
		Nproc = int(line2[1].split()[0])
	if (line2[0].lower().find('run_step') > -1):
		Runstep = int(line2[1].split()[0])
	if (line2[0].lower().find('run_avg') > -1):
		Runavg = int(line2[1].split()[0])
	if (line2[0].lower().find('run_dt') > -1):
		Rundt = float(line2[1].split()[0])
	if (line2[0].lower().find('avg_dt') > -1):
		Avgdt = float(line2[1].split()[0])
	if (line2[0].lower().find('mfile') > -1):
		Mfile2 = line2[1]
	if (line2[0].lower().find('sfile') > -1):
		Sfile2 = line2[1]
	if (line2[0].lower().find('ifile') > -1):
		Ifile2 = line2[1]				
		
f.close()		


currdir = os.getcwd()+'/'
Eqdsk = currdir+Eqdsk
	
power = np.zeros(Nbeam)
energy = np.zeros(Nbeam)

for i in range(Nbeam):
	power[i] = float(Power[i])
	energy[i] = float(Energy[i])

sname = '00000'+str(Shot)
if (Shot >= 10):
	sname = '0000'+str(Shot)
elif (Shot >= 1.e2):
	sname = '000'+str(Shot)
elif (Shot >= 1.e3):
	sname = '00'+str(Shot)
elif (Shot >= 1.e4):
	sname = '0'+str(Shot)
elif (Shot >= 1.e5):
	sname = str(Shot)
	
tname = '%06i'%Time
	
print('Powers',power)
print('Energy',energy)

plasma_name = 's'+sname+'.'+tname

if (not Ifile == ''):
	pp_data=run_nubeam('ch_tool',run_step=Runstep,run_avg=Runavg,time_step=Rundt,time_avg=Avgdt,
			zeff=Zeff,dbeam=Bdiff,nbeam=Nbeam,beam_power=power,beam_energy=energy,
			geqdsk=Eqdsk,NE=currdir+'NE.dat',TE=currdir+'TE.dat',TI=currdir+'TI.dat',VT=currdir+'VT.dat',nproc=Nproc,plasma_state=plasma_name,
			init_file=Ifile2,mdescr=Mfile2,sconfig=Sfile2,step_file=stepfile)
else:	
	pp_data=run_nubeam('ch_tool',run_step=Runstep,run_avg=Runavg,time_step=Rundt,time_avg=Avgdt,
			zeff=Zeff,dbeam=Bdiff,nbeam=Nbeam,beam_power=power,beam_energy=energy,
			geqdsk=Eqdsk,NE=currdir+'NE.dat',TE=currdir+'TE.dat',TI=currdir+'TI.dat',VT=currdir+'VT.dat',nproc=Nproc,plasma_state=plasma_name,mdescr=Mfile2,sconfig=Sfile2,step_file=stepfile)


nubeam_data=get_nubeam_data(run_avg=0)

xplasma_data=get_xplasma_data(plasma_state=plasma_name)

pp_data=postprocess_nubeam_data(nubeam_data,xplasma_data)

pbeam = np.zeros(101)
pperp = np.zeros(101)
ppll = np.zeros(101)
prho = np.zeros(101)
pcur = np.zeros(101)


for i in range(101):
	pbeam[i] = 2.0 / 3.0 * (pp_data['pperp'][i] + pp_data['ppll'][i]*0.5)
	pperp[i] = pp_data['pperp'][i] * 1.602 * 1.e3
	ppll[i] = pp_data['ppll'][i] * 1.602 * 1.e3

	pbeam[i] = pbeam[i] * 1.602 * 1.e3
	prho[i] = pp_data['rho'][i]
	pcur[i] = pp_data['curbeam'][i]

f1 = open('chease_pres','w');
f2 = open('chease_curr','w')
f1.write('101\n')
f2.write('101\n')

for i in range(101):

	f1.write('%9.6f %9.6f %9.6f \n'%(prho[i],ppll[i],pperp[i]))
	f2.write('%9.6f %9.6f \n'%(prho[i],pcur[i]))

f1.close()
f2.close()

with open('nubeam_out1d','w') as f:
	nubeam_data_list=['rho','pbe','pbi','curbeam','nbeami','pperp','ppll','pfast','tqb','sbedep']
	nubeam_unit_list=['[-]','[W/m^3]','[W/m^3]','[A/m^2]','[10^19/m^3]','[keV/m^3]','[keV/m^3]','[keV/m^3]','[N/m^2]','[10^19/m^3/s]']
	ndata=len(nubeam_data_list)
	f.write((ndata*'{:23s}').format(*[' '.join([x,y]) for x,y in zip(nubeam_data_list,nubeam_unit_list)])+'\n')
	for i in range(len(pp_data['rho'])):
		f.write((ndata*'{:23.7e}').format(*[pp_data[x][i] for x in nubeam_data_list])+'\n')

from scipy.io.netcdf import netcdf_file
import glob

with open('nubeam_out0d','w') as f:
	f.write('petot [W] %f\n'%pp_data['petot'])
	f.write('pitot [W] %f\n'%pp_data['pitot'])
	f.write('ptot [W] %f\n'%pp_data['ptot'])
	f.write('wfast [J] %f\n'%pp_data['wfast'])
	f.write('torque [Nm] %f\n'%pp_data['tqbt'])
	with netcdf_file(glob.glob('*_scalars_out.cdf')[0],'r') as scalars_out:
		f.write('beam-beam neutron rate [#/s] %e\n'%(scalars_out.variables['bbntot'].data))
		f.write('beam-thermal neutron rate [#/s] %e\n'%(scalars_out.variables['btneut'].data))
		f.write('total neutron rate [#/s] %e\n'%(scalars_out.variables['bbntot'].data+scalars_out.variables['btneut'].data))

print('Normal end of nubeam.py')

