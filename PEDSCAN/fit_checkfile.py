import os,sys
import numpy as np
from shutil import copyfile, move

class chk_file:

	def readfiles(self,filename,flag):

		if   (flag.lower() == 'te'):	  critval = 100
		elif (flag.lower() == 'ti'):      critval = 100
		elif (flag.lower() == 'ne'):      critval = 0.8
		else: critval = 1

		f = open(filename,'r')
		datx = []
		datz = []
		datv = []
		dats = []

		line = f.readline()

		while True:
			line = f.readline()
			if not line: break
			dat = line.split()
			if (len(dat) > 2):
				ind = 2
			else:
				ind = 1

			if (dat[ind].find('***') == -1 and dat[ind].lower().find('nan')==-1):
				if ((float(dat[0])>0. and float(dat[0])<2.1 and float(dat[ind])>critval) or (float(dat[0])>=2.1)):
					datx = np.append(datx,float(dat[0]))
					datv = np.append(datv,float(dat[ind]))
					if (ind == 3):
						datz = np.append(datz,float(dat[1]))
						
					else:
						datz = np.append(datz,0.)

					try:
						if (float(dat[ind+1]) > 0.):
							dats = np.append(dats,float(dat[ind+1]))
						else:
							dats = np.append(dats,float(dat[ind])/10.)
					except:
						dats = np.append(dats,0.5)#float(dat[ind])/10.)

		f.close()
		return (datx,datz,datv,dats)


	def renewfiles(self,filename,flag='vt'):

		datx,datz,datv,dats = self.readfiles(filename,flag)
		ind = []

		for i in range(len(datx)-1):
				if (datx[i] > datx[i+1]):
					ind = np.append(ind,i)

		if len(ind)==0:
			Flag = False
			return datx,datz,datv,dats,Flag

		ind = np.append(ind,len(datx)-1)
		ind2 = np.copy(ind)
		ind2[0] = ind[0]+1
		for i in range(len(ind)-1):
			ind2[i+1] = ind[i+1]-ind[i]

		min_ind = int(min(ind2))
		max_ind = int(max(ind2))

		if (max_ind == min_ind):
			Flag = False
			return datx,datz,datv,dats, Flag

		datx_old = np.copy(datx)
		datz_old = np.copy(datz)
		datv_old = np.copy(datv)
		dats_old = np.copy(dats)

		datx = np.zeros(max_ind*len(ind))
		datz = np.zeros(max_ind*len(ind))
		datv = np.zeros(max_ind*len(ind))
		dats = np.zeros(max_ind*len(ind))

		for i in range(len(ind2)):
			if (ind2[i] == max_ind): break

		if (i==0):
			datx_ref = datx_old[0:int(ind2[0])]
			datz_ref = datz_old[0:int(ind2[0])]
			datv_ref = datv_old[0:int(ind2[0])]
			dats_ref = dats_old[0:int(ind2[0])]
		else:
			datx_ref = datx_old[int(ind[i-1])+1:int(ind[i]+1)]
			datz_ref = datz_old[int(ind[i-1])+1:int(ind[i]+1)]
			datv_ref = datv_old[int(ind[i-1])+1:int(ind[i]+1)]
			dats_ref = dats_old[int(ind[i-1])+1:int(ind[i]+1)]
		count  = 0
		counti = 0

		for i in range(len(datx)):

			if (datx_old[count] == datx_ref[counti]):
				datx[i] = datx_old[count]
				datz[i] = datz_old[count]
				datv[i] = datv_old[count]
				dats[i] = dats_old[count]
				count = count + 1
			else:
				datx[i] = datx_ref[counti]
				datz[i] = datz_ref[counti]
				datv[i] = datv_ref[counti]
				dats[i] = dats_ref[counti]	

			counti = counti + 1

			if (counti == len(datx_ref)):	counti = 0
			if (count == len(datx_old)):	count = count - 1

			if np.isnan(datv[i]):
				datv[i] = 0.
			if np.isnan(dats[i]):
				dats[i] = 0.1

			Flag = True

		return datx,datz,datv,dats,Flag

	def writefiles(self,filename,flag):

		datx,datz,datv,dats,Flag = self.renewfiles(filename,flag)

		if (Flag):
			copyfile(filename,filename+'_source')

		f = open(filename,'w')
		if (flag.lower() == 'te'):
			line = ' R[m]   Z[m]   TE[eV]        Error[eV] \n'
		elif (flag.lower() == 'ne'):
			line = ' R[m]   Z[m]   NE[10^18/m3]  Error[10^18/m3] \n'
		elif (flag.lower() == 'vt'):
			line = ' R[m]   Z[m]   VT[km/s]      Error[km/s] \n'
		elif (flag.lower() == 'ti'):
			line = ' R[m]   Z[m]   TI[eV]        Error[eV] \n'

		f.write(line)
		for i in range(len(datx)):
			line = ' %9.6f %9.6f %9.6f %9.6f \n'%(datx[i],datz[i],datv[i],dats[i])
			f.write(line)
		f.close()

		if (Flag):
			print('>>> %s data file is corrected'%flag)

		return Flag
