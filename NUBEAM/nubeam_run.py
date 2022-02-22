#!/usr/local/anaconda3/bin/python3
import os,sys
import ch_tool
from shutil import copyfile, move
import numpy as np
from exec_dirs import python3_exec
from exec_dirs import nubeam_dir as nubeam_exec

eq_file = sys.argv[1]
run_num = sys.argv[2]
bdiff = float(sys.argv[3])

run_dir = os.getcwd()
nubeam_dir = run_dir + '/NUBEAM'
profile_dir = run_dir + '/PROFILES'
out_dir1 = run_dir + '/OUTPUT'
out_dir2 = run_dir + '/PROFILES'
if run_num == 'f': save_dir = run_dir
else: save_dir = run_dir + '/../../NUBEAM/Bdiff_%03i'%(round(bdiff*100.,0))

ch = ch_tool.chease('chease_opt')

ch.load_eqdsk()
ch.read_hagermap('CHEASE/hager_map')
ch.use_ext_pressure = False
ch.use_ext_current = False

ch.load_kin_profile()
ch.kinprof_type = 2
ch.write_kinprof()
ch.kinprof_type = 1

ch.load_kin_profile()
ch.interpol()
ch.read_rho_psi_R('CHEASE/RHO_PSI_R')

f4 = open('CHEASE/log.chease','r')
linec = 0
while True:
	line = f4.readline()
	if not line: break

	if (line.find('MKSA') > -1):
		linec = 1

	if ((line.find('POLOIDAL BETA') > -1) and (linec == 1 )):
		bpc = float(line.split()[0])
	if ((line.find('WMHD')>-1) and (linec == 1)):
		wmhd2 = float(line.split()[0])
	if ((line.find('LI')>-1) and (linec == 1)):
		li = float(line.split()[0])
f4.close()
	
copyfile('TI.dat_new','NUBEAM/TI.dat')
copyfile('NE.dat_new','NUBEAM/NE.dat')
copyfile('TE.dat_new','NUBEAM/TE.dat')
copyfile('VT.dat_new','NUBEAM/VT.dat')

move('TI.dat_new','PROFILES/TI.dat_0')
move('NE.dat_new','PROFILES/NE.dat_0')
move('TE.dat_new','PROFILES/TE.dat_0')
move('VT.dat_new','PROFILES/VT.dat_0')

copyfile('nubeam_opt','NUBEAM/nubeam_opt')
copyfile(eq_file,'NUBEAM/geqdsk')
os.chdir(nubeam_dir)
os.system(python3_exec + ' ' + nubeam_exec)
os.chdir(run_dir)
ch.use_ext_pressure = True
ch.use_ext_current = True
ch.pres_file = nubeam_dir+'/chease_pres'
ch.curr_file = nubeam_dir+'/chease_curr'

ch.make_efit_constraint()

fastpv  = np.trapz(ch.pres_ex,x=ch.vol)
fastppv = np.trapz(ch.pres_exp,x=ch.vol)
thermp  = np.trapz(ch.pt,x=ch.vol)

wmhd = 1.5 * (thermp +  fastpv) * 1.e-3
wdia = 1.5 * (thermp +  fastppv) * 1.e-3

f = open(save_dir+'/nubeam_iter_result','w')
f.write('Beam diff\tW_th [kJ]\tW_fast [kJ]\tW_mhd [kJ]\tW_dia [kJ]\tI_FAST[MA]\tVLOOP[V]\n')
f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(bpc,  0.0             ,0.0             ,wmhd2,0.0,0.0            ,0.0))
f.write('%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n'%(bdiff,1.5*thermp*1.e-3,1.5*fastpv*1.e-3,wmhd,wdia,ch.excurt*1.e-6,ch.vloop2))
f.close()

f = open(save_dir+'/pre_prof','w')
f.write('%i\n'%len(ch.psin))
for i in range(len(ch.psin)):
	f.write('%9.6f\t%9.6f\t%9.6f\n'%(ch.psin[i],ch.pt[i]+ch.pres_ex[i],ch.pres_ex[i]))
f.close()

copyfile(out_dir1+'/EFIT_JCONST',  save_dir+'/EFIT_JCONST')
copyfile(out_dir1+'/Vneo.dat',     save_dir+'/Vneo.dat')
copyfile(nubeam_dir+'/chease_pres',save_dir+'/chease_pres')	
copyfile(nubeam_dir+'/chease_curr',save_dir+'/chease_curr')
copyfile(nubeam_dir+'/nubeam_out1d',save_dir+'/nubeam_out1d')
copyfile(nubeam_dir+'/nubeam_out0d',save_dir+'/nubeam_out0d')
copyfile(run_dir+'/neo_coefs',save_dir+'/neo_coefs')
print('>>> RUN finished')
