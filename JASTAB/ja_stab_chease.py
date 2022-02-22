#!/usr/local/anaconda3/bin/python3

import numpy as np
import os, sys
import subprocess
from shutil import copyfile
import ja_tool as jt
from scipy.interpolate import interp1d

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
chease_temp = sim.run_dir +'/CHE_TEMP'
mis_temp = sim.run_dir +'/MIS_TEMP'
elite_temp = sim.run_dir + '/ELI_TEMP'
equ_save_dir = sim.run_dir + '/equ_inp'
fname = 'equ_scan_'+str(nw)+'_'+str(nh)
work_dir = chease_temp +'/'+str(nw)+'_'+str(nh)
stab_dir = mis_temp +'/'+str(nw)+'_'+str(nh)
stabe_dir = elite_temp + '/'+str(nw)+'_'+str(nh)

savename='/gr_'+sim.RUN_ID

sim.make_directory(chease_temp,'CHEASE')
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
	
namelist_expeq = equ_save_dir + '/' + fname + '_expeq'
namelist_name = equ_save_dir + '/' + fname + '_name'
namelist_prof = equ_save_dir + '/chease_kinprof'
namelist_rot = equ_save_dir + '/chease_vtor'

os.chdir(work_dir)
copyfile(namelist_expeq, 'EXPEQ')
copyfile(namelist_name,'chease_namelist')

if(sim.use_kin_prof):
	copyfile(namelist_prof,'chease_kinprof')
if(sim.use_rot):
	copyfile(namelist_rot,'chease_vtor')

sim.read_chease_bnd('EXPEQ')
if (sim.run_stab == 'mishka'):
	sim.modify_chease_namelist('chease_namelist',sim.CNS,sim.CNT,sim.CNPSI,sim.CNCHI)
elif (sim.run_stab == 'elite'):
	sim.modify_chease_namelist('chease_namelist',sim.CNSE,sim.CNTE,sim.CNPSIE,sim.CNCHIE)
	
sim.run_chease()

sim.read_chease_bnd('EXPEQ.OUT')

CHEout = work_dir + '/NJA'

copyfile('EXPEQ','EXPEQ_temp_init')
copyfile('chease_namelist','chease_namelist_temp')
copyfile('EXPEQ.OUT','EXPEQ')
copyfile('EXPEQ.OUT','EXPEQ_temp')
sim.run_chease()
sim.read_chease_out(CHEout)

qa0 = sim.q[-1]
qa2 = qa0
qa0 = abs(qa0)
qf = interp1d(sim.psin,sim.q,'cubic')
q95 = qf(0.95)
psin = np.copy(sim.psin)
jav = np.copy(sim.jav2)
pprime = np.copy(sim.pprime)
IPEXP = abs(sim.IPEXP)
B0EXP = abs(sim.B0EXP)

print('init qa0',qa0)

jm0 = sim.find_max_ja(0.8,sim.psin,sim.jav1,1)
am0 = sim.find_max_ja(0.8,sim.psin,sim.alp1,1)
am20 = sim.find_max_ja(0.8,sim.psin,sim.alp2,1)
ball0 = sim.find_max_ja(0.8,sim.psin,sim.ball,0)

	
copyfile('NELITE','NELITE2')
elitename = work_dir + '/NELITE2'

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
	
		if (qdel > 0.8):
			mm = mm - 1
			qdel1 = mm - qa0*n
			print('qdel0, new qdel, m0, m1',qdel,qdel1,mm+1,mm)
		else:
			print('qdel, m',qdel,mm)

		qa_target = float(mm - sim.qdelfix) / float(n)
		print('qa-target',qa_target)
	
		sim.write_chease_namelist(None,sim.CNSE,sim.CNTE,sim.CNPSIE,sim.CNCHIE,B0EXP,IPEXP)	
		sim.write_chease_expeq(None,psin,pprime,jav,B0EXP,IPEXP)
		sim.qval_change_chease_namelist('chease_namelist',1.0,np.sign(qa2)*qa_target)
		sim.run_chease()
		
		sim.read_chease_out(CHEout)
		
		jm1 = sim.find_max_ja(0.8,sim.psin,sim.jav1,1)
		jme = sim.jav1[-1]
		am1 = sim.find_max_ja(0.8,sim.psin,sim.alp1,1)
		am21 = sim.find_max_ja(0.8,sim.psin,sim.alp2,1)
		ball1 = sim.find_max_ja(0.8,sim.psin,sim.ball,0)
		qa1 = sim.q[-1]
		qa2 = qa1
		qa1 = abs(qa1)
		qf = interp1d(sim.psin,sim.q,'cubic')
		q95 = qf(0.95)
		qdel1 = np.ceil(qa1*n) - qa1*n

		print('qdelfix,qdel1,qa0,qa1,n,mm')
		print(sim.qdelfix,qdel1,qa0,qa1,n,mm)

		if(abs(sim.qdelfix-qdel1)>sim.qdel_crit):
			print('q converge fail')
			mm = mm + 1
			qa_target = float(mm - sim.qdelfix) / float(n)
			print('qa-target',qa_target)
			sim.write_chease_namelist(None,sim.CNSE,sim.CNTE,sim.CNPSIE,sim.CNCHIE,B0EXP,IPEXP)
			sim.write_chease_expeq(None,psin,pprime,jav,B0EXP,IPEXP)
			sim.qval_change_chease_namelist('chease_namelist',1.0,np.sign(qa2)*qa_target)
			sim.run_chease()
	
			sim.read_chease_out(CHEout)

			jm1 = sim.find_max_ja(0.8,sim.psin,sim.jav1,1)
			jme = sim.jav1[-1]
			am1 = sim.find_max_ja(0.8,sim.psin,sim.alp1,1)
			am21 = sim.find_max_ja(0.8,sim.psin,sim.alp2,1)
			ball1 = sim.find_max_ja(0.8,sim.psin,sim.ball,0)
			qa1 = sim.q[-1]
			qa2 = qa1
			qa1 = abs(qa1)
			qdel1 = np.ceil(qa1*n) - qa1*n

			print('qdelfix,qdel1,qa0,qa1,n,mm')
			print(sim.qdelfix,qdel1,qa0,qa1,n,mm)
		
		misinp2 = work_dir + '/NMISHKA'
		os.chdir(stab_dir)
		copyfile(misinp2,'fort.12')

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
			if not (gr_n1 >= 0.):	gr_n1 = 0.
			gr_n2 = gr_n1 * alpw / float(gr_nn)
		except:
			gr_nn = moden[i]
			gr_n1 = 0.0
			gr_n2 = 0.0
		print('growth rate',gr_n1)
		f3.write('%i %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n'%(nw,nh,am0,am20,jm0/1.e6,am1,am21,jm1/1.e6,gr_nn,gr_n1,gr_n2,qdel1,abs(q95),ball0,ball1,jme/1.e6))
		
	elif(sim.run_stab == 'elite'):
	
		os.chdir(stabe_dir)
		copyfile(elitename,'eqin')
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

				sim.write_elite_inp(inp_name2,n,sim.qdelfix,sim.ngride,sim.CNCHIE+1,sim.psise,sim.ndist,sim.ias)
				
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
				f3.write('%i %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %i %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n'%(nw,nh,am0,am20,jm0/1.e6,am0,am20,jm0/1.e6,gr_nn,gr_n1,gr_n2,qdelf0,abs(q95),ball0,qmon,jme/1.e6))
f3.close()	
