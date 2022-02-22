import numpy as np
from MDS import mds

def _get_mse(shot):
	mse_dat = dict()
	mse_dat['is'] = False
	g = mds('kstar',int(shot))
	R = g.get('\\RRRGAM')
	mse_dat['ch'] = R[0]
	if len(mse_dat['ch']) == 0: return mse_dat
	mse_dat['is'] = True
	mse_dat['RRRGAM'] = R[1]
	Z = g.get('\\ZZZGAM')
	mse_dat['ZZZGAM'] = Z[1]
	for k in range(6):
		A = g.get('\\AA%iGAM'%(k+1))
		mse_dat['AA%iGAM'%(k+1)] = A[1]
	
	T = g.get('\\TGAMMA01');
	mse_dat['time'] = T[0]
	mse_dat['TGAMMA01'] = T[1]
	T = g.get('\\SGAMMA01');
	mse_dat['SGAMMA01'] = T[1]
	
	for k in range(2,len(mse_dat['ch'])+1):
		T = g.get('\\TGAMMA%02i'%k);
		mse_dat['TGAMMA%02i'%k] = T[1]
		T = g.get('\\SGAMMA%02i'%k);
		mse_dat['SGAMMA%02i'%k] = T[1]

	g.close()
	
	return mse_dat

def _make_mse(shot,time,dt):

	mse_dat = _get_mse(shot)
	if not mse_dat['is']: return mse_dat
	mse_dat['nch'] = len(mse_dat['ch'])
	mse_dat['tgam'] = np.zeros(mse_dat['nch'])
	mse_dat['sgam'] = np.zeros(mse_dat['nch'])
	mse_dat['fwt']  = np.ones(mse_dat['nch'])

	ind1    = np.where(mse_dat['time']>(time-0.5*dt))
	ind2    = np.where(mse_dat['time'][ind1]<(time+0.5*dt))
	for i in range(mse_dat['nch']):
		ind = i + 1
		mse_dat['tgam'][i] = np.mean(mse_dat['TGAMMA%02i'%ind][ind1][ind2])
		mse_dat['sgam'][i] = np.mean(mse_dat['SGAMMA%02i'%ind][ind1][ind2])
		mse_dat['fwt'][i]  = 1.
		if (mse_dat['RRRGAM'][i] < 1.or mse_dat['RRRGAM'][i]>2.1): mse_dat['fwt'][i] = 0.
	print('>>> Load MSE from MDS...')
	return mse_dat

def _read_mse_ext(n_file):

	f = open(n_file,'r')
	f.close()

	return
