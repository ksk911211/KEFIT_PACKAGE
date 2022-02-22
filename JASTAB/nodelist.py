#!/usr/bin/sh

import os,sys

def machine_file(node,filename='machine1'):

	node_n = int(float(node.split('node')[1]))
	
	if (node_n == 1):
		print('Do not run on nodemaster..')
		exit()
	elif (node_n < 6):
		proc_n = 8
	elif (node_n < 9):
		proc_n = 20
	elif (node_n < 11):
		proc_n = 12
	elif (node_n < 14):
		proc_n = 28
	else:
		print('Wrong node number...')
		exit()
		
	self.proc_n = proc_n

	f4 = open(filename,'w')
	
	for i in range(proc_n):
		f4.write('%s\n'%node)
	f4.close()

	return
	
def make_nodelist(node_list):
	node_list = node_list.lower()
	node_list = node_list.replace(' ','')
	node_list = node_list.split(',')
	for i in range(len(node_list)):
	
		node_n = int(float(node_list[i].split('node')[1]))
	
		if (node_n == 1):
			print('Do not run on nodemaster..')
			exit()
		elif (node_n < 6):
			node_list[i] = 'all.q@'+node_list[i]
		elif (node_n < 9 and node_n > 5):
			node_list[i] = 'E5-2680@'+node_list[i]
		elif (node_n < 11 and node_n > 8):
			node_list[i] = 'old_group.q@'+node_list[i]
		elif (node_n < 14 and node_n > 10):
			node_list[i] = 'E5-2690v4.q@'+node_list[i]
		else:
			print('Wrong node number...')
			exit()
	
	nodelist = node_list[0]
	for i in range(len(node_list)-1):
		nodelist = nodelist+','+node_list[i+1]
		
	return nodelist
