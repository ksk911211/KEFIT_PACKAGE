import numpy as np
from scipy.interpolate import interp1d
import os,sys
import eqdsk
## - https://omfit.io/_modules/omfit_classes/omfit_osborne.html#OMFITpFile
class GENpFile:
    """
    OMFIT class used to interface with Osborne pfiles
    cocosio = 1  # pFile's have CCW phi and CW theta

    :param filename: filename passed to OMFITobject class

    :param \**kw: keyword dictionary passed to OMFITobject class
    """

    if __name__ == "__main__":

        import pfile_convert
        pfile = pfile_convert.GENpFile()
        pfile._make_pfile()

    def __init__(self):

        self.legacy = False
        self.V      = dict()
        self.electron_charge       = 1.602e-19
        self.filename              = 'a.out'

        self.descriptions = dict()
        self.descriptions['ne']    = 'Electron density'
        self.descriptions['te']    = 'Electron temperature'
        self.descriptions['ni']    = 'Ion density'
        self.descriptions['ti']    = 'Ion temperature'
        self.descriptions['nb']    = 'Fast ion density'
        self.descriptions['pb']    = 'Fast ion pressure'
        self.descriptions['ptot']  = 'Total pressure'
        self.descriptions['omeg']  = 'Toroidal rotation: VTOR/R'
        self.descriptions['omegp'] = 'Poloidal rotation: Bt * VPOL / (RBp)'
        self.descriptions['omgvb'] = 'VxB rotation term in the ExB rotation frequency: OMEG + OMEGP'
        self.descriptions['omgpp'] = 'Diamagnetic term in the ExB rotation frequency: (P_Carbon)/dpsi / (6*n_Carbon)'
        self.descriptions['omgeb'] = 'ExB rotation frequency: OMGPP + OMGVB = Er/(RBp)'
        self.descriptions['er']    = 'Radial electric field from force balance: OMGEB * RBp'
        self.descriptions['ommvb'] = 'Main ion VXB term of Er/RBp, considered a flux function'
        self.descriptions['ommpp'] = 'Main ion pressure term of Er/RBp, considered a flux function'
        self.descriptions['omevb'] = 'Electron VXB term of Er/RBp, considered a flux function'
        self.descriptions['omepp'] = 'Electron pressure term of Er/RBp, considered a flux function'
        self.descriptions['kpol']  = 'KPOL=VPOL/Bp : V_vector = KPOL*B_vector + OMGEB * PHI_Vector'
        self.descriptions['N Z A'] = 'N Z A of ION SPECIES'
        self.descriptions['omghb'] = 'Hahm-Burrell form for the ExB velocity shearing rate: OMGHB = (RBp)**2/Bt * d (Er/RBp)/dpsi'
        self.descriptions['nz1']   = 'Density of the 1st impurity species'
        self.descriptions['vtor1'] = 'Toroidal velocity of the 1st impurity species'
        self.descriptions['vpol1'] = 'Poloidal velocity of the 1st impurity species'
        self.descriptions['nz2']   = 'Density of the 2nd impurity species'
        self.descriptions['vtor2'] = 'Toroidal velocity of the 2nd impurity species'
        self.descriptions['vpol2'] = 'Poloidal velocity of the 2nd impurity species'
        # There may be more impurity species, but let's stop here for now.

        self.units = dict()
        self.units['ne']    = '10^20/m^3'
        self.units['te']    = 'KeV'
        self.units['ni']    = '10^20/m^3'
        self.units['ti']    = 'KeV'
        self.units['nb']    = '10^20/m^3'
        self.units['pb']    = 'KPa'
        self.units['ptot']  = 'KPa'
        self.units['omeg']  = 'kRad/s'
        self.units['omegp'] = 'kRad/s'
        self.units['omgvb'] = 'kRad/s'
        self.units['omgpp'] = 'kRad/s'
        self.units['omgeb'] = 'kRad/s'
        self.units['ommvb'] = ''
        self.units['ommpp'] = ''
        self.units['omevb'] = ''
        self.units['omepp'] = ''
        self.units['er']    = 'kV/m'
        self.units['kpol']  = 'km/s/T'
        self.units['N Z A'] = ''
        self.units['omghb'] = ''
        self.units['nz1']   = '10^20/m^3'
        self.units['vtor1'] = 'km/s'
        self.units['vpol1'] = 'km/s'
        self.units['nz2']   = '10^20/m^3'
        self.units['vtor2'] = 'km/s'
        self.units['vpol2'] = 'km/s'

        for key in list(self.descriptions.keys()):
            if key in ('N Z A',):
                continue
            self.V[key] = dict()
            self.V[key]['data'] = np.array([0])
            self.V[key]['description'] = self.descriptions[key]
            self.V[key]['psinorm'] = np.array([0])
            self.V[key]['units'] = self.units[key]
            self.V[key]['derivative'] = np.array([0])

    def _read_eqfile(self,gfile_name):

        eq = eqdsk.eqdsk(gfile_name);
        eq.read_eqdsk_file()

        self.gfile = dict()
        self.gfile['psi_norm'] = eq.psin
        self.gfile['pres']     = eq.pres * 1.e-3 #kPa
        self.gfile['fpol']     = eq.fpol

        self.gfile['psi_v']    = eq.psirz
        self.gfile['ip']       = eq.ip
        self.gfile['bcentr']   = eq.bcentr
        self.gfile['psi_r']    = eq.R
        self.gfile['psi_z']    = eq.Z
        self.gfile['prhoR']    = eq.prhoR
        self.gfile['psia']     = eq.sbdy - eq.smag;

    def _read_profile(self,pfile_name):

        self.pfile = dict()
        if not os.path.isfile(pfile_name): exit()
        f = open(pfile_name,'r')
        ngrid = int(float(f.readline()))
        self.pfile['ngrid'] = ngrid;
        dat   = f.readline().split()
        self.pfile['zeff']  = float(dat[0])
        self.pfile['zimp']  = float(dat[1])
        self.pfile['amain'] = float(dat[2])
        self.pfile['aimp']  = float(dat[3])

        self.pfile['psi_norm'] = np.zeros(ngrid)
        self.pfile['te']    = np.zeros(ngrid) #keV
        self.pfile['ne']    = np.zeros(ngrid) #10(20)/m3
        self.pfile['ti']    = np.zeros(ngrid) #keV
        self.pfile['ni']    = np.zeros(ngrid) #10(20)/m3
        self.pfile['nit']   = np.zeros(ngrid) #10(20)/m3
        self.pfile['timp']  = np.zeros(ngrid) #keV
        self.pfile['nimp']  = np.zeros(ngrid) #10(20)/m3
        self.pfile['pt']    = np.zeros(ngrid) #kPa
        self.pfile['ptot']  = np.zeros(ngrid) #kPa
        self.pfile['pb']    = np.zeros(ngrid) #kPa
        self.pfile['nb']    = np.zeros(ngrid) #10(20)/m3
        self.pfile['vtor']  = np.zeros(ngrid) #km/s

        for i in range(ngrid):
            dat   = f.readline().split()
            self.pfile['psi_norm'][i] = float(dat[0])
            self.pfile['te'][i]       = float(dat[1])
            self.pfile['ne'][i]       = float(dat[2]) * 1.e-1; #19->20
            self.pfile['ti'][i]       = float(dat[3])
            self.pfile['nit'][i]      = float(dat[4]) * 1.e-1; #19->20
            self.pfile['vtor'][i]     = float(dat[5]) * -1.    #clockwise

        f.close()

        self.ni_over_nit   = self.pfile['zimp']*(self.pfile['zimp']-self.pfile['zeff'])
        self.ni_over_nit   = self.ni_over_nit / (self.pfile['zimp']-1.) / (self.pfile['zimp']+1.-self.pfile['zeff'])

        self.nimp_over_ni  = (self.pfile['zeff']-1.)/self.pfile['zimp']/(self.pfile['zimp']-self.pfile['zeff'])

        self.pfile['timp'] = np.copy(self.pfile['ti'])
        self.pfile['ni']   = self.pfile['nit'] * self.ni_over_nit
        self.pfile['nimp'] = self.pfile['ni']  * self.nimp_over_ni

        self.pfile['pte']   = self.electron_charge * 1.e20 * self.pfile['ne']*self.pfile['te']
        self.pfile['pti']   = self.electron_charge * 1.e20 * self.pfile['ni']*self.pfile['ti']
        self.pfile['ptimp'] = self.electron_charge * 1.e20 * self.pfile['nimp']*self.pfile['ti']
        self.pfile['pt']    = self.pfile['pte'] + self.pfile['pti'] + self.pfile['ptimp']

        # 1st impurity: Carbon 
        self.pfile['nz1']  = self.pfile['nimp']
        self.pfile['vtor1']= self.pfile['vtor'] * -1. #(CCW is +, KSTAR CW)

        # neglecting 2nd impurity
        self.pfile['nz2']  = np.zeros(ngrid)
        self.pfile['vtor2']= np.zeros(ngrid)
        self.pfile['vpol2']= np.zeros(ngrid)
       
    def _read_neofile(self,neo_file):

        if not os.path.isfile(neo_file):
            self.neofile = dict()
            self.neofile['ngrid'] = 101
            self.neofile['psi_norm'] = np.linspace(0.,1.,101)
            self.neofile['ki'] = np.zeros(101)
            self.neofile['b02']= np.ones(101)
            return

        self.neofile = dict()
        f = open(neo_file,'r')
        line = f.readline()
        ngrid = int(float(f.readline()))
        self.neofile['ngrid'] = ngrid
        self.neofile['psi_norm'] = np.zeros(ngrid)
        self.neofile['ki'] = np.zeros(ngrid)        
        self.neofile['b02']= np.zeros(ngrid)
        
        for i in range(ngrid):
            dat = f.readline().split()
            self.neofile['psi_norm'][i] = float(dat[0])
            self.neofile['ki'][i]  = float(dat[1])
            self.neofile['b02'][i] = float(dat[5])
        f.close()

    def _read_beamfile(self,beam_file):

        if not os.path.isfile(beam_file):
            self.beamfile = dict()
            ngrid = 101
            self.beamfile['ngrid'] = ngrid
            self.beamfile['rho_norm'] = np.linspace(0.,1.,ngrid)
            self.beamfile['p_para']   = np.zeros(ngrid) #kPa
            self.beamfile['p_perp']   = np.zeros(ngrid) #kPa
            self.beamfile['p_tot']    = np.zeros(ngrid) #kPa
            self.beamfile['density']  = np.ones(ngrid)*0.05 #10(20)/m3
            self.beamfile['torque']   = np.zeros(ngrid) #Nm
            return

        f = open(beam_file,'r')
        ngrid = -1;
        while True:
            line = f.readline()
            if not line: break
            ngrid += 1
        f.close()
        f = open(beam_file,'r')
        line = f.readline()    
        self.beamfile = dict()
        self.beamfile['ngrid'] = ngrid
        self.beamfile['rho_norm'] = np.zeros(ngrid)
        self.beamfile['p_para']   = np.zeros(ngrid) #kPa
        self.beamfile['p_perp']   = np.zeros(ngrid) #kPa
        self.beamfile['p_tot']    = np.zeros(ngrid) #kPa
        self.beamfile['density']  = np.zeros(ngrid) #10(20)/m3
        self.beamfile['torque']   = np.zeros(ngrid) #Nm

        for i in range(ngrid):
            dat = f.readline().split()
            self.beamfile['rho_norm'][i] = float(dat[0])
            self.beamfile['p_para'][i]   = float(dat[6]) * 1.602
            self.beamfile['p_perp'][i]   = float(dat[5]) * 1.602
            self.beamfile['density'][i]  = float(dat[4]) * 1.e-1
            self.beamfile['torque'][i]   = float(dat[8])

        self.beamfile['p_tot'] = (2.*self.beamfile['p_perp']+self.beamfile['p_para'])/3.

        f.close()

    def _map(self):

        self.map=dict()
        self.map['rho_to_psi'] = interp1d(self.gfile['prhoR'][:,1],self.gfile['prhoR'][:,0])
        self.map['psi_to_Rlow']= interp1d(self.gfile['prhoR'][:,0],self.gfile['prhoR'][:,2])
        self.map['psi_to_Rhig']= interp1d(self.gfile['prhoR'][:,0],self.gfile['prhoR'][:,3])
        self.map['Rlow_to_psi']= interp1d(self.gfile['prhoR'][:,2],self.gfile['prhoR'][:,0])

    def _post_field(self):

        ngrid    = self.pfile['ngrid']
        self.pfile['Rlow'] = self.map['psi_to_Rlow'](self.pfile['psi_norm'])
        self.pfile['Btheta'] = np.zeros(ngrid)
        fpf   = interp1d(self.gfile['psi_norm'],self.gfile['fpol'])
        self.pfile['fpol'] = fpf(self.pfile['psi_norm'])
        for i in range(1,ngrid-1):
            self.pfile['Btheta'][i] = (self.pfile['psi_norm'][i+1]-self.pfile['psi_norm'][i-1])/ (self.pfile['Rlow'] [i+1]-self.pfile['Rlow'] [i-1])/ self.pfile['Rlow'] [i] * self.gfile['psia']
            if self.pfile['psi_norm'][i] > 0.993: self.pfile['Btheta'][i] = self.pfile['Btheta'][i-1]
        self.pfile['Btheta'][-1] = self.pfile['Btheta'][-2];

    def _post_beam_profile(self):

        ngrid    = self.pfile['ngrid']
        psi_norm = self.map['rho_to_psi'](self.beamfile['rho_norm'])
        pbf      = interp1d(psi_norm,self.beamfile['p_tot'])
        nbf      = interp1d(psi_norm,self.beamfile['density'])
        ptf      = interp1d(self.gfile['psi_norm'],self.gfile['pres'])

        for i in range(ngrid):
            self.pfile['ptot'][i] = ptf(self.pfile['psi_norm'][i])
            self.pfile['pb'][i]   = self.pfile['ptot'][i] - self.pfile['pt'][i]
            if self.pfile['pb'][i]<0.: self.pfile['pb'][i] = 0.
            if not self.beamfile['p_tot'][0]==0:
               self.pfile['nb'][i]   = self.pfile['pb'][i]/(pbf(self.pfile['psi_norm'][i])+1.e-3)*nbf(self.pfile['psi_norm'][i])
            else:
               self.pfile['nb'][i]   = nbf(self.pfile['psi_norm'][i])

    def _post_neo_profile(self):

        ngrid = self.pfile['ngrid']
        self.pfile['vpol']    = np.zeros(ngrid)
        self.pfile['vpoli']   = np.zeros(ngrid)
        self.pfile['vpole']   = np.zeros(ngrid)
        self.pfile['vpolimp'] = np.zeros(ngrid)
        self.pfile['vpolE']   = np.zeros(ngrid)
        kif   = interp1d(self.neofile['psi_norm'],self.neofile['ki'])       
        dti   = self._derivative(self.pfile['psi_norm'],self.pfile['ti'] )  / self.gfile['psia']
        dpi   = self._derivative(self.pfile['psi_norm'],self.pfile['pti'])  / self.gfile['psia']
        dpe   = self._derivative(self.pfile['psi_norm'],self.pfile['pte'])  / self.gfile['psia']
        dpimp = self._derivative(self.pfile['psi_norm'],self.pfile['ptimp'])/ self.gfile['psia']
        b02f  = interp1d(self.neofile['psi_norm'],self.neofile['b02'])
        self.pfile['b02']    = b02f(self.pfile['psi_norm'])

        for i in range(ngrid):
                psi_norm     = self.pfile['psi_norm'][i]
                RR           = self.pfile['Rlow'][i]
                fpol         = self.pfile['fpol'][i]
                btheta       = self.pfile['Btheta'][i]
                b02          = self.pfile['b02'][i]
                ne           = self.pfile['ne'][i]
                ni           = self.pfile['ni'][i]
                nimp         = self.pfile['nimp'][i]

                self.pfile['vpol'][i]    = -dti[i]   * kif(psi_norm) * fpol * btheta / b02
                self.pfile['vpoli'][i]   = -dpi[i]   * (RR**2) * btheta/ni/1.e20 /   fpol /self.electron_charge
                self.pfile['vpole'][i]   = +dpe[i]   * (RR**2) * btheta/ne/1.e20 /   fpol /self.electron_charge
                self.pfile['vpolimp'][i] = -dpimp[i] * (RR**2) * btheta/nimp/1.e20 / fpol /self.electron_charge/self.pfile['zimp']

    def _post_rotation_profile(self):

        ngrid = self.pfile['ngrid']

        self.pfile['omeg']  = np.zeros(ngrid) #kRad/s
        self.pfile['omegp'] = np.zeros(ngrid)
        self.pfile['omgvb'] = np.zeros(ngrid)
        self.pfile['omgpp'] = np.zeros(ngrid)
        self.pfile['omgeb'] = np.zeros(ngrid)

        self.pfile['er']    = np.zeros(ngrid)
        self.pfile['ommgp'] = np.zeros(ngrid)
        self.pfile['ommvb'] = np.zeros(ngrid)
        self.pfile['ommpp'] = np.zeros(ngrid)

        self.pfile['omevb'] = np.zeros(ngrid)
        self.pfile['omepp'] = np.zeros(ngrid)

        self.pfile['kpol']  = np.zeros(ngrid)
        self.pfile['omghb'] = np.zeros(ngrid) 


        #Less accurate in that it relies on main-ion force balance.
        for i in range(ngrid):
            btheta  = self.pfile['Btheta'][i] + 1.e-4
            self.pfile['omeg'][i] = +self.pfile['vtor'][i]   /self.pfile['Rlow'][i]
            self.pfile['ommgp'][i]= +self.pfile['vpol'][i]   *self.pfile['fpol'][i]/(self.pfile['Rlow'][i]**2)/btheta
            self.pfile['ommvb'][i]= +self.pfile['omeg'][i]   +self.pfile['ommgp'][i]
            self.pfile['ommpp'][i]= +self.pfile['vpoli'][i]  *self.pfile['fpol'][i]/(self.pfile['Rlow'][i]**2)/btheta
            self.pfile['omgeb'][i]= +self.pfile['ommpp'][i]  +self.pfile['ommvb'][i]

            self.pfile['er'][i]   = +self.pfile['omgeb'][i]  *self.pfile['Rlow'][i]* btheta

            self.pfile['omegp'][i]= +self.pfile['vpol'][i]   *self.pfile['fpol'][i]/(self.pfile['Rlow'][i]**2)/btheta
            self.pfile['omgvb'][i]= +self.pfile['omeg'][i]   +self.pfile['omegp'][i]
            self.pfile['omgpp'][i]= +self.pfile['vpolimp'][i]*self.pfile['fpol'][i]/(self.pfile['Rlow'][i]**2)/btheta
            self.pfile['omgvb'][i]= +self.pfile['omgeb'][i]  -self.pfile['omgpp'][i]

            self.pfile['omepp'][i]= +self.pfile['vpole'][i]  *self.pfile['fpol'][i]/(self.pfile['Rlow'][i]**2)/btheta
            self.pfile['omevb'][i]= +self.pfile['omgeb'][i]  -self.pfile['omepp'][i]
            
            self.pfile['kpol'][i] = +self.pfile['vpol'][i]   /btheta

        dw_dpsi = self._derivative(self.pfile['psi_norm'],self.pfile['omgeb'])/self.gfile['psia']
        for i in range(ngrid):
            btor = np.sqrt(self.pfile['b02'][i])
            self.pfile['omghb'][i]= +(self.pfile['Rlow'][i]*self.pfile['Btheta'][i])**2 / btor * dw_dpsi[i]

        self.pfile['vpol1'] = self.pfile['vpol']

    def _derivative(self,x,y):

        f = interp1d(x,y,'slinear')
        deps = 1.e-4;
        dy= np.copy(x)

        for i in range(1,len(y)-1):
            xx   = x[i]
            dy[i]= (f(xx+deps)-f(xx-deps))/2./deps;

        xx = x[0]
        dy[0] = (f(xx+deps)-f(xx))/deps
        xx = x[-1]
        dy[-1]= (f(xx)-f(xx-deps))/deps

        return dy

    def _put_variables(self):

        self.V['N Z A'] = dict()
        self.V['N Z A']['description'] = 'N Z A of ION SPECIES'
        self.V['N Z A']['N'] = np.array([self.pfile['zimp'], 1, 1])
        self.V['N Z A']['Z'] = np.array([self.pfile['zimp'], 1, 1])
        self.V['N Z A']['A'] = np.array([self.pfile['aimp'], self.pfile['amain'], self.pfile['amain']])

        thermal_list = ['ne','te','ni','ti']
        fast_list    = ['nb','pb','ptot',]
        flow_list    = ['omeg','omegp','omgvb','omgpp','er','ommvb','ommpp','omevb','omepp','kpol','omghb']
        imp_list     = ['nz1','vtor1','vpol1','nz2','vtor2','vpol2']
        for key in thermal_list:
            self.V[key]['psinorm']    = self.pfile['psi_norm']
            self.V[key]['data']       = self.pfile[key]
            self.V[key]['derivative'] = self._derivative(self.pfile['psi_norm'],self.pfile[key])

        for key in fast_list:
            self.V[key]['psinorm']    = self.pfile['psi_norm']
            self.V[key]['data']       = self.pfile[key]
            self.V[key]['derivative'] = self._derivative(self.pfile['psi_norm'],self.pfile[key])
        for key in flow_list:
            self.V[key]['psinorm']    = self.pfile['psi_norm']
            self.V[key]['data']       = self.pfile[key]
            self.V[key]['derivative'] = self._derivative(self.pfile['psi_norm'],self.pfile[key])
        for key in imp_list:
            self.V[key]['psinorm']    = self.pfile['psi_norm']
            self.V[key]['data']       = self.pfile[key]
            self.V[key]['derivative'] = self._derivative(self.pfile['psi_norm'],self.pfile[key])

    def _save(self):

        lines = []
        for key in self.V.keys():
            if key == 'N Z A':
                if self.legacy:
                    break
                lines.append('%i N Z A of ION SPECIES\n' % (len(self.V[key]['A']),))
                for i in range(len(self.V[key]['A'])):
                    lines.append(" %f   %f   %f\n" % (self.V[key]['N'][i], self.V[key]['Z'][i], self.V[key]['A'][i]))
            else:
                if len(self.V[key]['data']) == 1:
                    continue
                lines.append("%i psinorm %s(%s) d%s/dpsiN\n" % (len(self.V[key]['data']), key, self.V[key]['units'], key))
                for i in range(len(self.V[key]['data'])):
                    lines.append(" %f   %f   %f\n" % (self.V[key]['psinorm'][i], self.V[key]['data'][i], self.V[key]['derivative'][i]))
        with open(self.filename, 'w') as f:
            f.writelines(lines)

    def _make_pfile(self,gfile,pfile,bfile,nfile):
        self._read_eqfile(gfile)
        self._read_profile(pfile)
        self._read_beamfile(bfile)
        self._read_neofile(nfile)

        self._map()
        self._post_field()
        self._post_beam_profile()
        self._post_neo_profile()
        self._post_rotation_profile()
        self._put_variables()
        self._save()
