import os,sys
import numpy as np
import matplotlib.pyplot as plt


file_n = '/Users/sk42//Works/2022_KSTAR_EXP/n12_RMP/PROFILES/p031509.006100_kin_4'
variables = ['omeg','omgeb','omepp','omegp']

def read_pfile(file_n):
	dat = dict();
	with open(file_n,'r') as f:
		while True:
			line = f.readline()
			if not line: break
			line = line.split()
			var_n= line[2].split('(')[0];
			try: units= line[2].split('(')[1][:-1]
			except: units=''
			varn = int(line[0])
			dat[var_n] = dict();
			dat[var_n]['val'] = np.zeros((varn,3));
			dat[var_n]['unit']= units
			for i in range(int(line[0])): 
				line = f.readline().split()
				dat[var_n]['val'][i,:] = np.array(line,dtype='float')
	return dat

if not os.path.isfile(file_n): print('No such file'); exit()
dat = read_pfile(file_n);
llegend = []; slegend = [];
print(dat.keys())
for vname in variables:
	line1, = plt.plot(dat[vname]['val'][:,0],dat[vname]['val'][:,1])
	llegend.append(line1); slegend.append('%s [%s]'%(vname,dat[vname]['unit']))
plt.legend(llegend,slegend)
plt.xlabel('$\\psi_N$')
plt.show()
