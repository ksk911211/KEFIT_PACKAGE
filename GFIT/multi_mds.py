from multiprocessing import Process, Queue
from MDS import mds 
import numpy as np 

class multi_mds():
	
	def __init__(self,shot,device,CH,t0,t1,sampling):
		self._shot = shot
		self._device = device
		self._ch = CH 
		self._t0 = t0
		self._t1 = t1 
		self._sampling = sampling
		# device is 'MCT' or 'ECE' 

	def list_mds(self):
		self._lst =[]
		ch = list(self._ch)
		for i in ch:
			num = '%i'%(i)
			num = num.zfill(2)
			if self._device == 'MCT':
				name = 'MC1T' + num
			elif self._device == 'ECE':
				name = 'ECE' + num
			else: pass; 
			self._lst.append(name)

	def mds_load(self,node,output):
		g = mds('kstar',self._shot)
		res_name = 'resample(\\' + node + ',' +str(self._t0) + ',' + str(self._t1) + ',' + str(self._sampling) + ')'
		tt, yy = g.get(res_name)
		ch_num = int(node[-2:]) 
		output.put([ch_num,tt,yy]);
		
		msg = 'Parallel MDS loading for ...' + node 
		print(msg) 

	def multi_load(self):
		self.list_mds()
		procs = [] 
		outputs = []
		self.t = [] 
		self.y = [] 
		self.ch = []

		for name in self._lst:
			output = Queue()
			outputs.append(output);
			proc = Process(target=self.mds_load,args=(name,output))
			procs.append(proc)

		for p in procs:
			p.start()	
			
		for i in range(0,len(procs)):
			ch_num, tt, yy = outputs[i].get();
			outputs[i].close(); 

			self.ch.append(ch_num)
			self.t.append(tt)
			self.y.append(yy) 
