import os,sys
import subprocess
import numpy as np
import time
from exec_dirs import scheduler, node_default, node_force, qsub_exec, qstat_exec

def make_batch_script(filename,nodelist,log_e,log_o,command,run_name=None):

	if nodelist == None:
		nodelist2 = node_default
	else:
		nodelist2 = nodelist

	if node_force:	nodelist2 = node_default

	f = open(filename,'w')
	f.write('#!/bin/bash -l\n')
	
	if scheduler.lower() == 'pbs':	
		if not run_name == None:	f.write('#PBS -N '+ run_name + ' \n')
		else:	f.write('#PBS -N test \n')
		f.write('#PBS -q '+ nodelist2 + ' \n')
		f.write('#PBS -V \n')
		f.write('#PBS -e '+ log_e + ' \n')
		f.write('#PBS -o '+ log_o + ' \n')
		f.write('cd $PBS_O_WORKDIR\n')
		f.write('module load intel\n')
		f.write('module load gefit\n')
	elif scheduler.lower() == 'sge':
		if not run_name == None:        f.write('#$ -N '+ run_name + ' \n')
		else:   f.write('#$ -N test \n')
		f.write('#$ -q '+ nodelist2 + ' \n')
		f.write('#$ -V \n')
		f.write('#$ -e '+ log_e + ' \n')
		f.write('#$ -o '+ log_o + ' \n')

	f.write(command)
	f.close()
	return

def submit_batch_script(filename):

	os.system('chmod 777 '+filename)
	status, output = subprocess.getstatusoutput(qsub_exec + ' ' + filename)
	if status == 0:	print(output)
	else:	print('Submitting the script is failed')

	if scheduler.lower() == 'pbs':	runid = int(output.split('.')[0])
	elif scheduler.lower() == 'sge':	runid = int(output.split()[2])

	return runid

def batch_script_status():

	status, output = subprocess.getstatusoutput(qstat_exec)

	return output

def read_node_status():

	stat,out = subprocess.getstatusoutput('qstat -f > qstat.log')

	nodename = []
	used = []
	load = []
	nastat = []

	f = open('qstat.log','r')
	while True:
		line = f.readline()
		if not line: break

		if (line.find('lx24-amd64') > -1):
			line = line.split()
			nodename.append(line[0])
			used.append(line[2])
			load.append(line[3])
			if len (line) > 5:
				print('>>> ',line[0],'has problem...')
				nastat.append('na')
			else:
				nastat.append('ok')

	f.close()
	os.remove('qstat.log')

	len2 = len(nodename)
	load2 = np.zeros(len2)
	used2 = np.zeros(shape=(len2,3))
	name2 = [0] * len2
	nastat2 = [0] * len2

	for i in range(len2):

		name = nodename[i].split('node')
		if (name[1].find('master')>-1):
			name = name[1].split('master')
		ind = int(name[1])-1

		try:
			load2[ind] = float(load[i])
		except:
			nastat[i] == 'na'
		used2[ind,0] = float(used[i].split('/')[0])
		used2[ind,1] = float(used[i].split('/')[1])
		used2[ind,2] = float(used[i].split('/')[2])
		
		name2[ind] = nodename[i].split('@')[1] 
		if nastat[i] == 'na':
			nastat2[ind] = 'na'

	return (load2,used2,name2,nastat2)

