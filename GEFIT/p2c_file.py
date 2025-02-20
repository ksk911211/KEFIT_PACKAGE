import os,sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

pfile_n = 'p154383.03230_w8099'
cfile_n = 'c154383.03230_w8099'
npsi    = 401

def read_pfile(pfile_n):
        dat = dict();
        with open(pfile_n,'r') as f:
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

if not os.path.isfile(pfile_n): print('No such file'); exit()
dat = read_pfile(pfile_n);

psi_n = np.linspace(0,1,npsi)
tefactor = 1.e3;
tifactor = 1.e3;
vtfactor = 1.e2;

if dat['te']['unit'].lower()=='kev':     tefactor = 1;
if dat['ti']['unit'].lower()=='kev':     tifactor = 1;
if dat['vtor1']['unit'].lower()=='km/s': vtfactor = 1;

nefactor = float('1.e+%s'%dat['ne']['unit']. lower().split('^')[1].split('/')[0])*1.e-19
nifactor = float('1.e+%s'%dat['ni']['unit']. lower().split('^')[1].split('/')[0])*1.e-19
nIfactor = float('1.e+%s'%dat['nz1']['unit'].lower().split('^')[1].split('/')[0])*1.e-19

tef = interp1d(dat['te']['val'][:,0],   dat['te']['val'][:,1]   *tefactor)
nef = interp1d(dat['ne']['val'][:,0],   dat['ne']['val'][:,1]   *nefactor)
tif = interp1d(dat['ti']['val'][:,0],   dat['ti']['val'][:,1]   *tifactor)
nif = interp1d(dat['ni']['val'][:,0],   dat['ni']['val'][:,1]   *nifactor)
nIf = interp1d(dat['nz1']['val'][:,0],  dat['nz1']['val'][:,1]  *nIfactor)
vtf = interp1d(dat['vtor1']['val'][:,0],dat['vtor1']['val'][:,1]*vtfactor)

te  = tef(psi_n)
ne  = nef(psi_n)
ti  = tif(psi_n)
ni  = nif(psi_n)
nI  = nIf(psi_n)
vt  = vtf(psi_n)

ZI  = dat['Z']['val'][0][0]
AI  = dat['Z']['val'][0][2]
Zi  = dat['Z']['val'][1][0]
Ai  = dat['Z']['val'][1][2]

ni_ov_ne = (ni+nI) / ne;
Zeff = 1+(1-ni_ov_ne)*ZI;

zeff  = np.mean(Zeff[np.where(abs(psi_n-0.8)<0.1)])
with open(cfile_n,'w') as f:
        f.write('%4i\n'%npsi)
        f.write('%8.6f\t%8.6f\t%8.6f\t%8.5f\n'%(zeff,ZI,Ai,AI))
        for i in range(npsi):
                f.write('%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.5f\n'%(psi_n[i],te[i],ne[i],ti[i],ni[i]+nI[i],vt[i]))
