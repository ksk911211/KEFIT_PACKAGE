## Concept is developed by Y.H.Lee
import os, sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import eqdsk

mse_chn = 25
mse_rad_err = dict()
mse_rad_err['NB1A'] = np.ones(mse_chn) #mm
mse_rad_err['NB1B'] = np.ones(mse_chn)
mse_rad_err['NB1C'] = np.ones(mse_chn)

mse_rad_err['NB1A'][0:10]  = [32.1913, 27.7642, 23.6813, 20.2123, 17.0620, 14.5040, 12.2766, 10.4646, 8.85010, 7.63049] 
mse_rad_err['NB1A'][10:20] = [6.62952, 5.89612, 5.40161, 5.11340, 5.02441, 5.11572, 5.39307, 5.83228, 6.41235, 7.13867]     
mse_rad_err['NB1A'][20:25] = [8.00610, 8.94043, 10.0667, 11.2344, 12.5549]

mse_rad_err['NB1B'][0:10]  = [124.427, 114.534, 105.092, 96.7568, 88.8531, 82.1125, 75.9073, 70.5181, 65.3152, 60.9685]
mse_rad_err['NB1B'][10:20] = [56.9159, 53.3745, 50.2976, 47.5640, 45.1602, 43.1921, 41.3933, 39.9307, 38.7666, 37.8420]
mse_rad_err['NB1B'][20:25] = [37.1472, 36.6963, 36.4189, 36.3428, 36.4434]

mse_rad_err['NB1C'][0:10]  = [32.1259, 32.8662, 33.3674, 33.6091, 33.6267, 33.4414, 33.0685, 32.5482, 31.8264, 31.0122 ]
mse_rad_err['NB1C'][10:20] = [30.0271, 28.9324, 27.7430, 26.4329, 25.0051, 23.5601, 21.9141, 20.2148, 18.4749, 16.6453 ]
mse_rad_err['NB1C'][20:25] = [14.7280, 12.8577, 10.7783, 8.75977, 6.59985]

def mse_coef(Eb,massn=2.):

	kev2j = 1.60218*1e-16
	mi    = massn * 1.6735575*1e-27 
	vb    = np.sqrt(2.* Eb * kev2j / mi)

	mse_omega = np.ones(mse_chn)
	mse_alpha = np.ones(mse_chn)

	mse_omega[0:13] = [-4.20953,-3.68296,-3.18795,-2.75943,-2.36308,-2.03564,-1.74603,-1.50704,-1.29152,-1.12724,-0.992012,-0.893677,-0.829453]
	mse_omega[13:25] = [-0.796097,-0.793655,-0.819618,-0.877419,-0.963249,-1.07442,-1.21326,-1.37982,-1.56078,-1.78113,-2.01253,-2.27760]
	
	mse_alpha[0:15] = [148.411,147.306,146.200,145.174,144.148,143.221,142.317,141.481,140.619,139.844,139.065,138.325,137.623,136.936,136.264]
	mse_alpha[15:25]= [135.647,135.003,134.393,133.814,133.247,132.692,132.182,131.649,131.160,130.664]
		
	mse_omega = np.array(mse_omega) * np.pi / 180.
	mse_alpha = np.array(mse_alpha) * np.pi / 180.

	AA1 = np.zeros(mse_chn)
	AA2 = np.zeros(mse_chn)
	AA3 = np.zeros(mse_chn)
	AA4 = np.zeros(mse_chn)
	AA5 = np.zeros(mse_chn)
	AA6 = np.zeros(mse_chn)
	
	for i in range(mse_chn):
		AA1[i] = -np.cos(mse_omega[i]+mse_alpha[i])
		AA2[i] = +np.sin(mse_alpha[i])
		AA3[i] = +np.cos(mse_alpha[i])
		AA4[i] = 0.
		AA5[i] = -np.cos(mse_omega[i]) / vb
		AA6[i] = -1./vb

	return (AA1,AA2,AA3,AA4,AA5,AA6)

def make_vtheta(eq,psin,ne,ni,te,ti,rbp,rbm,delta,bp,zeff,zimp,model=1):
	
	ndat   = len(ne)
	ki     = np.zeros(ndat)
	ft     = np.zeros(ndat)
	vpol   = np.zeros(ndat)
	me_hag = 9.109 * 1.e-31
	e0     = 1.602 * 1.e-19
	qf     = interp1d(eq.psin,eq.q)
	
	for i in range(ndat):
		r = rbp[i] - eq.rmag
		rc =  (rbp[i] + rbm[i])*0.5
		q = qf(psin[i])
		eps_hag = r/rc
		eps_hag2 = +0.67*(1.-1.4*(delta[i]**2))*eps_hag

		ne_hag  = ne[i] * 1.e19
		ni_hag  = ni[i] * 1.e19
		te_hag  = te[i] * 1.e3
		ti_hag  = ti[i] * 1.e3 #eV
		if (eps_hag <= 0.0): eps_hag = 1.e-8
		Z1_hag = ne_hag/ni_hag
		Z2_hag = zeff
		Z_hag  = np.power((1.**2.0)*Z1_hag*Z2_hag,0.25)
	
		zle_hag = 31.3 - np.log(np.sqrt(ne_hag)/te_hag)
		zli_hag = 30.0 - np.log((Z_hag **3.0)*np.sqrt(ni_hag)/(ti_hag**1.5))
	
		ft[i] = 1. - np.sqrt((1.-eps_hag)/(1.+eps_hag)) * (1.-eps_hag2)/(1.+2.*np.sqrt(eps_hag2))
	
		nues_hag = 6.921*1.e-18 * abs(q) * rc * ne_hag * Z2_hag      * zle_hag / (te_hag**2) / (eps_hag**1.5)
		nuis_hag = 4.900*1.e-18 * abs(q) * rc * ni_hag * (Z_hag**4.0) * zli_hag / (ti_hag**2) / (eps_hag**1.5) / 2.
	
		deltapsi = rc * np.sqrt(eps_hag)*np.sqrt(2*te_hag*me_hag/e0)
	
		psid_hag = abs(eq.smag-eq.sbdy) * (1-psin[i])
		H= (0.6/Z2_hag**4)*np.exp(-abs(psid_hag/(3.3*np.log(eps_hag**1.5* nues_hag + 2. )*deltapsi)))
		ft[i] = ft[i] * (1-H)

		alp0_hag = -1.17*(1.-ft[i])/(1.-0.22*ft[i]-0.19*(ft[i]**2))
		ki[i] = ((alp0_hag + 0.25*(1-(ft[i]**2))*np.sqrt(nuis_hag))/ (1.+0.5*np.sqrt(nuis_hag)) + 0.315*(nuis_hag**2)*(ft[i]**6)) / (1.+0.15*(nuis_hag**2)*(ft[i]**6))

		if model > 1:

			if i == 0:	dti = (ti[1]-ti[0])/(rbp[1]-rbp[0])
			elif i == ndat-1:	dti = (ti[-1]-ti[-2])/(rbp[-1]-rbp[-2])
			else:	dti = (ti[i+1]-ti[i-1])/(rbp[i+1]-rbp[i-1])
			dti = dti * 1.e3
			R_mag = eq.rmag # magnetic axis 
			fc = 1 - 1.46*np.sqrt(r/rbp[i]) + 0.46*(r/rbp[i])*np.sqrt(r/rbp[i])
			K1 = 0.8839 * fc / (0.3477 + 0.4058 * fc)
			btor = abs(eq.fpol[0]) / rbp[i]
			vpol[i] = K1 * dti * btor / (bp[i]**2 + btor**2)
		else:
			if i==0:	dti = (ti[1]-ti[0])/(psin[1]-psin[0])
			elif i==ndat-1: dti = (ti[-1]-ti[-2])/(psin[-1]-psin[-2])
			else :	dti = (ti[i+1]-ti[i-1])/(psin[i+1]-psin[i-1])
			btor = abs(eq.bcentr)
			dti = dti / abs(eq.sbdy - eq.smag) * 1.e3
			ki[0] = 0.
			vpol[i] = ki[i] * eq.fpol[0] * bp[i] / (btor**2) * dti
			vpol[i] = vpol[i]
	return (vpol)

def make_komega(RR,PP,VT,VP,BT,BR,BZ):

	k = np.copy(VT)
	omega = np.copy(VT)

	for i in range(len(VT)):
		k[i] = VP[i] / np.sqrt(BR[i]**2+BZ[i]**2)
		omega[i] = (VT[i] + k[i]*BT[i])/RR[i]

	return k, omega

def get_vals(gfile,pfile,vfile='VT_fit.dat',mseR=[1.8]):

	f = open(pfile,'r')
	ndat = int(f.readline())
	psin2 = np.zeros(ndat);	ni1 = np.zeros(ndat); ti1 = np.zeros(ndat); ne1 = np.zeros(ndat); te1 = np.zeros(ndat)
	line = f.readline().split()
	zimp = float(line[1]);	zeff = float(line[0]);
	for i in range(ndat):
		line = f.readline().split()
		psin2[i] = float(line[0])
		ni1[i]   = float(line[4])
		ti1[i]   = float(line[3])
		ne1[i]   = float(line[2])
		te1[i]   = float(line[1])
	f.close()
	pif = interp1d(psin2,ni1*ti1*1.e19*1.e3*1.602*1.e-19)
	pi1 = np.zeros(ndat)
	for i in range(1,ndat-1): 
		deps1 = psin2[i] - 1.e-4
		deps2 = psin2[i] + 1.e-4
		pi1[i] = (pif(deps2)-pif(deps1))/2.e-4
	pi1[0] = 2.*pi1[1] - pi1[2]
	pi1[-1]= 2.*pi1[-2]- pi1[-3]	

	f = open(vfile,'r')
	psin1 = np.zeros(101);	vt1 = np.zeros(101);
	line = f.readline();	line = f.readline()
	for i in range(101):
		line = f.readline().split()
		psin1[i] = float(line[1])
		vt1[i]    = float(line[2])*1.e3
	f.close()

	eq = eqdsk.eqdsk(gfile)
	eq.read_eqdsk(gfile)
	eq.make_grid()
	eq.make_rho_R_psin()
	eq.get_flux_contour()
	eq.construct_2d_field()
	eq.shape_postprocess(True)
	ind = int(eq.nh/3) 
	for i in range(eq.nh-3):
		if (eq.z[i]*eq.z[i+1]<=0.):
			ind = i+1
			break

	brf = interp1d(eq.r,eq.br[ind,:]*np.sign(-eq.ip))
	bzf = interp1d(eq.r,eq.bz[ind,:]*np.sign(-eq.ip))
	btf = interp1d(eq.r,eq.bt[ind,:])
	ppf = interp1d(psin2,pi1)
	pf  = interp1d(eq.psin,eq.pres)
	r1f = interp1d(eq.prhoR[:,0],eq.prhoR[:,2])
	nif = interp1d(psin2,ni1)

	ind2 = np.where(mseR<1.0);
	mseR2 = np.copy(mseR)
	mseR2[ind2] = 2.0
	
	br = brf(mseR2)
	bz = bzf(mseR2)
	bt = btf(mseR2)

	r1 = r1f(psin1)
	br1 = brf(r1)
	bz1 = bzf(r1)
	bt1 = btf(r1)

	rbpf = interp1d(eq.prhoR[:,0],eq.prhoR[:,2])
	rbmf = interp1d(eq.prhoR[:,0],eq.prhoR[:,3])

	rbp = rbpf(psin2)
	rbm = rbmf(psin2)
	br2 = brf(rbp);	bz2 = bzf(rbp);
	bp2 = np.sqrt(br2**2+bz2**2)
	delta = 0.5*(eq.triu+eq.tril) * (psin2**2)

	vpol = make_vtheta(eq,psin2,ne1,ni1,te1,ti1,rbp,rbm,delta,bp2,zeff,zimp,model=1)
	vpf  = interp1d(psin2,vpol)
	vp1  = vpf(psin1)
	k,omega = make_komega(r1,psin1,vt1,vp1,bt1,br1,bz1)
	kf = interp1d(psin1,k)
	omegaf = interp1d(psin1,omega)
	vt = np.copy(mseR)
	vp = np.copy(mseR)
	pn = np.copy(mseR)
	pnf1 = interp1d(eq.prhoR[:,2],eq.prhoR[:,0])
	pnf2 = interp1d(eq.prhoR[:,3],eq.prhoR[:,0])
	for i in range(len(mseR)):
		if mseR2[i] < eq.rmag: pn[i] = pnf2(mseR2[i])
		elif mseR2[i] < eq.prhoR[-1,2]:	pn[i] = pnf1(mseR2[i])
		else:	pn[i] = 1. #+ (pnf1(eq.prhoR[-1,2])-pnf1(eq.prhoR[-1,2]-0.1))/0.1*(mseR[i]-eq.prhoR[-1,2])
		k1 = kf(pn[i])
		w1 = omegaf(pn[i])
		vp[i] = k1 * np.sqrt(br[i]**2+bz[i]**2) * bz[i] / abs(bz[i])
		vt[i] = -k1*bt[i] + w1*mseR2[i]

	pp = ppf(pn)
	p  = pf(pn)
	ni = nif(pn) * 1.e19
	
	e0 = 1.602*1.e-19
	er = np.zeros(len(ni))
	ppr= np.copy(p)
	for i in range(len(er)):
		er[i] = pp[i]*mseR2[i]*bz[i]/e0/ni[i] + vt[i]*bz[i] + vp[i]*bt[i]#/zimp
#		er[i] = pp[i]*mseR[i]*bz[i]/zimp/e0/ni[i] + vt[i]*bz[i] + vp[i]*bt[i]/zimp

	return (br,bz,bt,er)

def Er_mse_corr(gfile='',pfile='',vfile='VT_fit.dat',tgamma=[0.],rgamma=[1.8],ebeam=90):

	br,bz,bt,er = get_vals(gfile,pfile,vfile,rgamma)
	AA1,AA2,AA3,AA4,AA5,AA6 = mse_coef(ebeam)

	tgamma2 = np.copy(tgamma)
	for i in range(len(tgamma)):
		BZ = (tgamma[i]*(AA2[i]*bt[i] + AA3[i]*br[i])- AA5[i]*er[i]) / (AA1[i]-tgamma[i]*AA4[i])
		tgamma2[i] = (AA1[i]*BZ)/ (AA2[i]*bt[i] + AA3[i]*br[i] + AA4[i]*BZ)

	return tgamma2
