#!/usr/local/anaconda3/bin/python3

import numpy as np
import os, sys
import time
import ja_tool as jt
from shutil import copyfile, move
from progress import update_progress
from nodelist import make_nodelist
from batch_run import *
from exec_dirs import chease_dir,jastab_dir,jastabc_dir,jastabh_dir,node_machine,node_force,node_default

namelist0 = None
try:
	namelist0 = sys.argv[1]
except:
	pass

if (namelist0 == None):
	print('Command should include [input_file] [run_option]')
	print('Available option: collect, clear, clearall, run(default)')
	exit()

try:
	if (sys.argv[2] == 'collect'):
		collect = True
	else:
		collect = False
	if (sys.argv[2] == 'clear'):
		remove = True
	else:
		remove = False
	if (sys.argv[2] == 'clearall'):
		removeall = True
	else:
		removeall = False

except:
	collect = False
	remove = False
	removeall = False

currdir = os.getcwd()
namelist = currdir + '/' + namelist0

# --- Initialise
sim = jt.simulation(namelist)
run_dir = sim.run_dir
sim.make_directory(run_dir,'run')
copyfile(namelist,run_dir+'/'+namelist0)

if (collect):
        print('Finalizing process --> Collect data')
        sim.batch_run = False
        if (sim.run_equ == 'helena'):
                sim.generate_output_helena()
        elif (sim.run_equ == 'chease'):
                sim.generate_output_chease()
        exit()

if (remove or removeall):
	print('Remove temp data')
	sim.batch_run = False
	if (remove):
		sim.clear_equ_dat()
	else:
		sim.clear_equ_dat(True)
	exit()

use_batch = False

try:
	use_batch = sys.argv[2]
	use_batch = sim.str_to_logic(use_batch)
except:
	pass

# --- Main script directories
ja_gen = jastab_dir

if (sim.run_equ == 'helena'):
	stab_gen = jastabh_dir
elif (sim.run_equ == 'chease'):
	stab_gen = jastabc_dir

# --- Input file check

if not (os.path.isfile(sim.init_HELinp)):
	print('No equilibrium file')
	print('File dir = %s'%sim.init_HELinp)
	exit()
if (sim.use_kin_prof):
	if not (os.path.isfile(sim.kin_prof_name)):
		print('No kinetic profile')
		print('File dir = %s'%sim.kin_prof_name)
		exit()

# --- Main
if(sim.use_rot):
	if not (sim.run_stab == 'elite'):
		print('Our mishka do not support rotation effect ')
		exit()

if not (sim.batch_run and not use_batch):

	start_time = time.time()
	if (sim.equ_type == 1 and not sim.use_kin_prof):
		print('Kinetic profile is needed for equ_type 1')
		exit()

	if(sim.equ_type == 0): #conver eqdsk to 

		eqdsk_dir = run_dir + '/EQDSK'
		sim.make_directory(eqdsk_dir,'eqdsk')
		os.chdir(eqdsk_dir)
		import ch_tool
		ch = ch_tool.chease()
		ch.kinprof_type = 1
		ch.chkin_file = sim.kin_prof_name

		copyfile(sim.init_HELinp,'geqdsk')
		ch.eqdsk_name = 'geqdsk'
		ch.target_psin = sim.target_bnd
		try:
			ch.load_eqdsk(8,sim.use_prev)
			if (sim.use_prev):	print('>>> use previous CHEASE namelist')
		except:
			ch.load_eqdsk(8)

		if (sim.run_equ == 'helena'):
			copyfile('CHEASE/NHELENA',currdir+'/input/fort_'+sim.RUN_ID)
			sim.init_HELinp = currdir +'/input/fort_'+sim.RUN_ID
		elif(sim.run_equ == 'chease'):
			copyfile('CHEASE/chease_namelist',currdir+'/input/chease_namelist_'+sim.RUN_ID)
			copyfile('CHEASE/EXPEQ',currdir+'/input/EXPEQ_'+sim.RUN_ID)

			sim.init_HELinp = currdir +'/input/EXPEQ_'+sim.RUN_ID
			sim.init_CHEinp = currdir+'/input/chease_namelist_' + sim.RUN_ID
		os.chdir(currdir)

	elif(sim.equ_type == 1):

		eqdsk_dir = run_dir + '/EQDSK'
		sim.make_directory(eqdsk_dir,'eqdsk')
		os.chdir(eqdsk_dir)
		copyfile(sim.init_HELinp,'geqdsk')
		copyfile(sim.kin_prof_name,'chease_kinprof')
		f = open('eqdsk_opt','w')
		f.write('CHEASE_NR = 60\n')
		f.write('CHEASE_NT = 60\n')
		f.write('CHEASE_NRMAP = 200\n')
		f.write('CHEASE_NTMAP = 150\n')
		f.write('ZEROEDGEP = True\n')
		f.close()

		sim.make_chease_opt()

		os.system(chease_dir + ' chease_opt nomap')

		if (sim.adjust_prof and sim.equ_type == 1):
			#copyfile(sim.kin_prof_name,sim.kin_prof_name+'_old')
			move('OUTPUT/chease_kinprof_new',sim.kin_prof_name+'_new')

		copyfile('CHEASE/chease_namelist','CHEASE/temp')
		f = open('CHEASE/temp','r')
		f2 = open('CHEASE/chease_namelist','w')
		while True:
			line = f.readline()
			if not line: break
			if (line.find('!') == -1 and line.find('NIDEAL') > -1):
				line = 'NIDEAL = 8 \n'
			f2.write(line)

		f.close()
		f2.close()

		if (sim.run_equ == 'helena'):
			copyfile('CHEASE/NHELENA',currdir+'/input/fort_'+sim.RUN_ID)
			sim.init_HELinp = currdir +'/input/fort_'+sim.RUN_ID
		elif(sim.run_equ == 'chease'):
			copyfile('CHEASE/chease_namelist',currdir+'/input/chease_namelist_'+sim.RUN_ID)
			copyfile('CHEASE/EXPEQ',currdir+'/input/EXPEQ_'+sim.RUN_ID)

			sim.init_HELinp = currdir +'/input/EXPEQ_'+sim.RUN_ID
			sim.init_CHEinp = currdir+'/input/chease_namelist_' + sim.RUN_ID

		os.chdir(currdir)

	if (sim.run_equ == 'chease'):
		sim.modify_chease_namelist2(sim.init_CHEinp)	
	
	#print('dir')
#	nodelist_t = sim.nodelist.split(',')
#	node_num = np.size(nodelist_t)
	#nodelist1 = nodelist_t[0]
	batch_num = []

#	for i in range(node_num-1):
#		nodelist1 = nodelist1 + ',' + nodelist_t[i+1].split()[0]

	if node_machine.lower() == 'fusma':	nodelist1 = make_nodelist(sim.nodelist)
	else:	nodelist1 = ''
		
	nw = int(sim.grid_n1)
	nh = int(sim.grid_n2)

	print('>>> Run dir',run_dir)

	script_dir = run_dir + '/script'
	log_dir = run_dir + '/batch_log'

	sim.make_directory(run_dir,'run')
	sim.make_directory(log_dir,'batch log')
	sim.make_directory(script_dir,'script dir')
		
	scriptname = 'jascan_' + sim.RUN_ID
	
	if (sim.use_kin_prof):
		sim.read_kin_prof(sim.kin_prof_name)

	if (sim.run_equ == 'helena'):
		sim.make_scan_equilibrium_helena()
	elif (sim.run_equ == 'chease'):
		sim.make_scan_equilibrium_chease()

	for i in range(nw):
		for j in range(nh):
			f2name = scriptname + '_batch.' + str(i+1) + '_' + str(j+1)
			f2name_e = log_dir + '/' + f2name + '.e'
			f2name_o = log_dir + '/' + f2name + '.o'
			f2name = script_dir + '/' + f2name
			
			command = 'cd ' + currdir + ' \n'
			command = command + sim.python_dir + ' '+ stab_gen + ' ' + namelist +' '+str(i+1)+' '+str(j+1)
			command = command +'\n'
			make_batch_script(f2name,nodelist1,f2name_e,f2name_o,command,'JASTAB_SCAN')

			runid = submit_batch_script(f2name)
			batch_num.append(str(runid))
			
	print('>>> Scripts are submitted...')
	print('>>> Collect stability data...')
	batch_chk = 1
	bat_count0 = 1
	while ( bat_count0 > 0):
		output = batch_script_status()
		bat_count = 0
		for i in range(nw):
			for j in range(nh):
				index1 = (nh)*i + j
				batch_chk = output.find(batch_num[index1])
				if (batch_chk > -1):
					bat_count = bat_count + 1;

		if not (bat_count == bat_count0):
			bat_count0 = bat_count
			if not sim.batch_run:
				update_progress(1.0 - bat_count0/nw/nh)	
			else:
				prog_line = '>>> Progress -> ' + str(round(100-(100*bat_count0/nw/nh),0)) + ' %'
				print(prog_line)
		if (bat_count0 > 0): 
			time.sleep(30)

	print('>>> Job scripts are finished...')
	print('>>> Finalizing process --> Collect data...')
	
	if (sim.run_equ == 'helena'):
		sim.generate_output_helena()
	elif (sim.run_equ == 'chease'):
		sim.generate_output_chease()
		
	copyfile(namelist,run_dir+'/'+namelist0)

	print('>>> Elapsed time = %9.6f (s)'%(time.time()-start_time))

	print('>>> RUN FINISHED!')
	
else:
	log_dir = run_dir + '/batch_log'
	script_dir = run_dir + '/script'
	sim.make_directory(log_dir,'batch log')
	sim.make_directory(script_dir,'script dir')
	fbatch = script_dir + '/ja_batch_' + str(sim.RUN_ID)

	fbatch_e = log_dir + '/ja_batch_' + str(sim.RUN_ID) + '.e'
	fbatch_o = log_dir + '/ja_batch_' + str(sim.RUN_ID) + '.o'

	command = 'cd ' + currdir + ' \n'
	command = command + sim.python_dir + ' ' + ja_gen + ' ' + namelist0 + ' True \n'
#	file.write('mv ' + fbatch + '* ' + run_dir +  '/ \n')
	command = command + '\n'

	make_batch_script(fbatch,node_default,fbatch_e,fbatch_o,command,'JASTAB')

	runid = submit_batch_script(fbatch)
	print('Job %s is submitted'%runid)
	print('>>> Main script is submitted!')
