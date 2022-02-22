#!/usr/local/anaconda3/bin/python3

import numpy as np
import os, sys
import subprocess
from shutil import copyfile
import ja_tool as jt

namelist = sys.argv[1]
nw = int(sys.argv[2])
nh = int(sys.argv[3])

sim = jt.simulation(namelist)

nsize1 = np.size(sim.modens.split(','))
nsize2 = np.size(sim.modenes.split(','))

moden = np.zeros(shape=(nsize1,1))
modene = np.zeros(shape=(nsize2,1))

for i in range(nsize1):
	moden[i] = int(sim.modens.split(',')[i])
for i in range(nsize2):
	modene[i] = int(sim.modenes.split(',')[i])

# --- directories
hel_temp = sim.run_dir +'/HEL_TEMP'
mis_temp = sim.run_dir +'/MIS_TEMP'
elite_temp = sim.run_dir + '/ELI_TEMP'
equ_save_dir = sim.run_dir + '/equ_inp'
fname = 'equ_scan_'+str(nw)+'_'+str(nh)
work_dir = hel_temp +'/'+str(nw)+'_'+str(nh)
stab_dir = mis_temp +'/'+str(nw)+'_'+str(nh)
stabe_dir = elite_temp + '/'+str(nw)+'_'+str(nh)

savename='/gr_'+sim.RUN_ID

sim.make_directory(hel_temp,'HELENA')
sim.make_directory(work_dir,'WORK')
if (sim.run_stab == 'mishka'):
	sim.make_directory(mis_temp,'MISHKA')
	sim.make_directory(stab_dir,'STAB MISHKA')
	savename= stab_dir + savename
elif (sim.run_stab == 'elite'):
	sim.make_directory(elite_temp,'ELITE')
	sim.make_directory(stabe_dir,'STAB ELITE')
	savename = stabe_dir + savename
	
f3 = open(savename,'w')
	
os.chdir(work_dir)
namelist = equ_save_dir + '/' + fname
os.system('cp ' + namelist + ' fort.10')

namelist_prof = equ_save_dir + '/chease_kinprof'
namelist_rot = equ_save_dir + '/chease_vtor'

if(sim.use_kin_prof):
	copyfile(namelist_prof,'chease_kinprof')
if(sim.use_rot):
	copyfile(namelist_rot,'chease_vtor')

helinp = 'fort.10'
hel_shape,hel_prof,hel_phys,hel_num,xiab,bvac,rvac,eps = sim.get_HELENA(helinp)
NPTS,NPTP = sim.get_num_opt(hel_num)

print('mapping r #',NPTS,'mapping t #',NPTP)
if(NPTS < 1 or NPTP < 1):
	print('failed to read grid #')
	exit()
	
os.chdir(work_dir)
namelist = equ_save_dir + '/' + fname
os.system('cp ' + namelist + ' fort.10')
outname = work_dir +'/fort.20'

os.system('mv fort.10 tempin')
xiab0, btor0 ,initb0 = sim.new_helinp('tempin',1.0,1.0,1.0)
xiab1 = xiab0/1.00
btor1 = btor0*1.00
initb1,isrun1 = sim.RUN_HEL2('tempin',xiab1,btor1,initb0)
am0,am20, jm0, qa0, ball0, jme = sim.read_hel2(outname,NPTS)

os.system('mv eliteinp eliteinp2')
elitename = work_dir + '/eliteinp2'

xiab1 = xiab0/1.01
btor1 = btor0*1.01

initb1,isrun1 = sim.RUN_HEL2('tempin',xiab1,btor1,initb1)
am1, am21, jm1, qa1, ball1, jme = sim.read_hel2(outname,NPTS)
w1 = 1.01

if (sim.run_stab == 'mishka'):
	nn = np.size(moden)
	print('SELECTED STABILITY CODE = MISHKA')
elif (sim.run_stab == 'elite'):
	nn = np.size(modene)
	print('SELECTED STABILITY CODE = RUN ELITE')
	if(sim.use_rot):
		sim.read_rot_prof(sim.rot_prof_name)
		sim.elite_rot_prof(sim.PSI_ROT,sim.FROT)

for i in range(nn):
	if (sim.run_stab == 'mishka'):
		os.chdir(work_dir)
		n = moden[i]	
		qdel = np.ceil(qa0*n) - qa0*n
		mm = np.ceil(qa0*n)

		qdel1 = mm - qa1*n
		w = 1.0
		print(qdel1,qdel,w1,w)

		if(qdel1 == qdel):
			w1 = 1.01
			xiab1 = xiab0/w1
			btor1 = btor0*w1
			initb1,isrun1 = sim.RUN_HEL2('tempin',xiab1,btor1,initb1)
			am1, am21, jm1, qa1, ball1, jme = sim.read_hel2(outname,NPTS)
			qdel1 = mm - qa1*n
			print(qdel1,qdel,w1,w)

		if (qdel > 0.8):
			mm = mm - 1
			qdel = qdel -1
			qdel1 = qdel1 - 1
			print(qdel1,qdel,w1,w)	

		iterc = 0
		while ( abs(sim.qdelfix-qdel1) > sim.qdel_crit and iterc < 3):
			if (qdel1 == qdel):
				print('q saddle point')
				if (iterc >= 3):
					w2 = w1
				else:
					w2 = w1 + (sim.qdelfix-qdel1)/(qdel1-qdel+1.e-4)*(w1-w)
					iterc = iterc + 1
			else:
				w2 = w1 + (sim.qdelfix-qdel1)/(qdel1-qdel)*(w1-w)

			print(qdel1,qdel,w1,w)
			xiab1 = xiab0/w2
			btor1 = btor0*w2
			initb1,isrun1 = sim.RUN_HEL2('tempin',xiab1,btor1,initb1)		
			am1, am21, jm1, qa1, ball1, jme = sim.read_hel2(outname,NPTS)

			qdel = qdel1
			w = w1
			qdel1 = mm - qa1*n
			w1= w2

		print('qdelfix,qdel1,qa0,qa1,n,w1')
		print(sim.qdelfix,qdel1,qa0,qa1,n,w1)

		misinp2 = work_dir + '/fort.12'
		os.chdir(stab_dir)
		os.system('cp ' + misinp2 + ' ./')

		if not(sim.highq):
			if (qa1 > 5.0):
				sim.highq = True

		sim.write_mishka_inp('fort.10',n,sim.gridn,sim.psis,sim.xr1,sim.sig1,sim.xr2,sim.sig2,sim.ias)

		if (sim.highq):
			mn = 41
			if (n > 2):
				mn = 41
			if (n > 5):
				mn = 51
			if (n > 10):
				mn = 71
			if (n > 16):
				mn = 91
		else:
			mn = 31
			if (n > 2):
				mn = 31
			if (n > 5):
				mn = 41
			if (n > 10):
				mn = 51
			if (n > 16):
				mn = 71	

		mis_dir1 = sim.mis_dir + str(mn)+'_501'
		print(n,mn,mis_dir1)
		os.system(mis_dir1)
		try:
			status, output = subprocess.getstatusoutput('cat fort.20 | grep INST')
			gr_temp = output.split(':')[1].split()
			gr_nn = -1*int(gr_temp[0])
			gr_chk = int(gr_temp[1])
			if (gr_chk < 20):
				gr_n1 = float(gr_temp[2])
			else:
				gr_n1 = 0.0;
			npsi,psi,dpdpsi,ffp,te,ne,ti,zeff,zimp,alpw = sim.read_elite(elitename)
			gr_n1 = np.sqrt(gr_n1)
			if not (gr_n1 >=0.):	gr_n1 = 0.
			gr_n2 = gr_n1 * alpw / float(gr_nn)
		except:
			gr_nn = moden[i]
			gr_n1 = 0.0
			gr_n2 = 0.0
		print('growth rate',gr_n1)
		f3.write('%i %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n'%(nw,nh,am0,am20,jm0/1.e6,am1,am21,jm1/1.e6,gr_nn,gr_n1,gr_n2,qdel1,qa1,ball0,ball1,jme/1.e6))
		
	elif(sim.run_stab == 'elite'):
	
		os.chdir(stabe_dir)
		os.system('cp ' + elitename +' ./eqin')
		for ii in range(1):
			n = int(modene[i])
			
			if (n < 5):
				print('mode n is too low for ELITE, try to increase it')

			if (n > 0):
				inp_name = 'work_' + str(n)	
				inp_name2 = 'work_' + str(n) + '.in'	
				elite_eq = sim.elite_dir + '/eliteeq5.7 ' + inp_name
				elite_vac = sim.elite_dir + '/elitevac5.7 ' + inp_name
				elite_run = sim.elite_dir + '/elite5.7 ' + inp_name
				elite_out = stabe_dir + '/' + inp_name + '.gamma'
				elite_eqout = stabe_dir + '/' + inp_name + '.eqout'
				elite_jedge = stabe_dir + '/' + inp_name + '.jedge'

				sim.write_elite_inp(inp_name2,n,sim.qdelfix,sim.ngride,NPTP,sim.psise,sim.ndist,sim.ias)
				
				status, output = subprocess.getstatusoutput(elite_eq)
				status, output = subprocess.getstatusoutput(elite_vac)
				status, output = subprocess.getstatusoutput(elite_run)
				
				isqmon = sim.read_elite_eq(elite_eqout,elite_jedge)
				if (isqmon):
					qmon = 0.0
				else:
					qmon = -1.0	
				
				file = open(elite_out)
				line = file.readline()
				line = file.readline()
				line = file.readline()
				file.close()
#				print(line)	
#				npsi,psi,dpdpsi,ffp,te,ne,ti,zeff,zimp,alpw = read_elite(elitename)
				gr_nn = int(line.split()[0])
				gr_n1 = float(line.split()[1])
				gr_n2 = float(line.split()[2])/4.0  #### change later
				qdelf0 = float(line.split()[4])
				print(n,'growth rate',gr_n1)
				f3.write('%i %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n'%(nw,nh,am0,am20,jm0/1.e6,am0,am20,jm0/1.e6,gr_nn,gr_n1,gr_n2,qdelf0,qa0,ball0,qmon,jme/1.e6))

f3.close()	
