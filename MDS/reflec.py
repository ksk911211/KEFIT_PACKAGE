import ctypes as _C
from MDSplus import Connection
#from MDSplus._mdsshr import _load_library, MdsException 
from MDSplus._mdsshr import MdsshrException as MdsException
import tkinter
#from Tkinter import *

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

#ConnectToMds=_load_library('MdsIpShr').ConnectToMds
#ConnectToMds.argtypes=[_C.c_char_p]
#DisconnectFromMds = _load_library('MdsIpShr').DisconnectFromMds
#DisconnectFromMds.argtypes = [_C.c_int]

class _Connection( Connection):
    """
    Updating 'Connection' class in 'MDSplus' to manange the connection to the server
    (1) hanging off the connection when termination
    (2) retry the connection by 'reconnect' method
    Written by D. K. Oh
    Last Modification : Aug 2012
    """
    def __del__(self):
        self.closeConnection()
        
    def closeConnection(self):
        if self.socket != -1:
            if False: #DisconnectFromMds(self.socket) == 0: 
                raise Exception("Error in disconnection")
            else:
                self.socket = -1
                
    def reconnect(self):
        if self.hostspec == None:
             raise MdsException("Error: no host specified")
        else:
             if self.socket != -1:
                print(self.socket)
                try:
                    self.closeConnection()
                except:
                    raise Exception("Error in resetting connection to %s" %(self.hostspec,))
             self.socket = ConnectToMds(self.hostspec)
             if self.socket == -1:
                raise Exception("Error connecting to %s" %(self.hostspec,))   
 #      self.socket = -1

class MDS(object):
    """
    Implementation of a connection to the MDSplus tree based on the class mds by Y. M. Jeon
    Written by D. K. Oh
    Last modification : Aug 2012
    """
    __DefaultTree = "KSTAR"
#    __DefaultServer = "ssh://ksk911211@203.230.126.229"
    __DefaultServer = "172.17.100.200:8005"

    def __init__(self, shot=None, tree =__DefaultTree, server=__DefaultServer):
        try:                    
            self.alist = {"tree":None, "shot":None, "server":None}
            self.__mds__ = _Connection(server)
            if shot is not None:
                self.open(shot, tree)
        except MdsException:
            raise MdsException("Error in the connection %s" %(server))
        except:
            raise Exception(" Unknown error in the connection %s" %(server))
        else:
            self.alist = {"tree":tree, "shot":shot, "server":server}
            
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        self.disconnect()
        
    def reset(self, shot=None, tree=None, server=__DefaultServer):
        if tree is None:
            tree = self.alist["tree"]
        if shot is None:
            shot = self.alist["shot"]
        
        self.alist = {"tree":tree, "shot":shot, "server":server}
        alist = self.alist
        __mds = self.__mds__ 
        __mds.hostspec = alist["server"]
        __mds.reconnect()
        if (shot is not None) and (tree is not None):
            self.open(tree, shot)

    def open(self, shot, tree =__DefaultTree):
        if shot is None:
            shot = self.alist["shot"]
        if shot is None:
            self.alist["tree"] = None
            self.alist["shot"] = None
            raise MdsException("Error in open : shot number is not specified")
        else:    
            self.close( self.alist["shot"], self.alist["tree"])
            try:
                self.__mds__.openTree(tree, shot)
            except:
                self.alist["tree"] = None
                self.alist["shot"] = None
                raise MdsException("Error in open : unknown error")    
            else:
                self.alist["tree"] = tree
                self.alist["shot"] = shot
                return self.__mds__

    def close(self, shot=None, tree=None): 
        if tree is None:
            tree = self.alist["tree"]
        if shot is None:
            shot = self.alist["shot"]
        if (shot is not None) and (tree is not None):
            try:
                self.__mds__.closeTree(tree, shot)
            except:
                raise MdsException("Error in close : unknown error")
            else:
                self.alist["shot"] = None
                self.alist["tree"] = None

    def disconnect(self):
        self.__mds__.closeConnection()
        
    def get_T0(self):
        try:
            ret_str = self.__mds__.get('\T0_STR').data()
        except:
            ret_str = None            
            raise MdsException("Error in get")
        return ret_str

    def get_sig(self, sigstr):
        try:
            t = self.__mds__.get('dim_of(%s)' %(sigstr)).data();
            v = self.__mds__.get(sigstr).data();
        except:
            t = numpy.ndarray([])
            v = numpy.ndarray([])             
            raise MdsException("Error in get")
        return t,v;

#-------------------------------------------------------------------------------#
def interp1(xp,x,y):
    import numpy as np
    
    x_=np.array([])
    y_=np.array([])
    minxp=np.min(xp)
    if minxp<x[0]:
        x_=np.append(x_,minxp)
        y_=np.append(y_,y[1]+(y[1]-y[0])/(x[1]-x[0])*(minxp-x[1]))
        
    x_=np.append(x_,x)
    y_=np.append(y_,y)

    maxxp=np.max(xp)
    if maxxp>x[-1]:
        x_=np.append(x_,maxxp)
        y_=np.append(y_,y[-2]+(y[-1]-y[-2])/(x[-1]-x[-2])*(maxxp-x[-2]))

    return np.interp(xp,x_,y_)

def gprofdat(shot,treename):
    import numpy as np
    import pickle
    REFS1 = "dim_of(\if_qx:foo"
    REFS2 = ")[0:*:100]"
    prof = dict(); prof['isref'] = False
    prof['time'] = np.zeros(8)
    prof['prof'] = dict()
#    with MDS(server="ssh://ksk911211@203.230.126.229") as mds:
    with MDS(server="172.17.100.200:8005") as mds:
        try:
            eq=mds.open(shot=shot, tree=treename)
        except: 
            print("Error #2")
        else:
            print(mds.alist)
            try: t = eq.get('%s%i%s'%(REFS1,0,REFS2))
            except: 
                print("No Reflectometry avail")
                filename = 'DATASAVE/%i'%shot+'/REF.npz'
                f = open(filename,'wb')
                pickle.dump(prof,f)
                f.close()
                comm='gzip '+filename
                os.system(comm)
                return
            prof['isref'] = True
            prof['time'][0] = t[0]
            prof['prof']['is'] = np.zeros(8)
            for k in range(7):
                try: t = eq.get('%s%i%s'%(REFS1,k+1,REFS2))
                except:
                  print("No Reflectometry avail")
                  prof['isref'] = False
                  filename = 'DATASAVE/%i'%shot+'/REF.npz'
                  f = open(filename,'wb')
                  pickle.dump(prof,f)
                  f.close()
                  comm='gzip '+filename
                  os.system(comm)
                  return
                  
                prof['time'][k+1] = t[0]
            istart = -1
            for k in range(8): 
                try:
                    s = eq.get('SetTimeContext($,$)',prof['time'][k], prof['time'][k]+0.02).data()
                    t = eq.get('dim_of(\\r_ref)').data()
                    r = eq.get('\\r_ref').data()
                    n = eq.get('\\ne_ref').data()
                    if istart == -1:
                        for i in range(len(r[0,:])-2):
                            if not (r[0,i] == r[0,i+1]):
                                istart = i+1; break;
                    r = r[:,istart:-1]; n = n[:,istart:-1]
                    rscale = int(len(r[0,:])/40)
                    
                    n = n[:,0::rscale]
                    r = r[:,0::rscale]
                    t = t[0::4]
                    n = n[0::4,:]
                    r = r[0::4,:]
                    prof['prof']['is'][k] = 1
                    prof['prof'][k] = dict()
                    prof['prof'][k]['t'] = t
                    prof['prof'][k]['r'] = r
                    prof['prof'][k]['n'] = n 
                except: print('Ref Time slice #%i/8 not avail'%(k+1))
    #            except: pass
            filename = 'DATASAVE/%i'%shot+'/REF.npz'
            f = open(filename,'wb')
            pickle.dump(prof,f)
            f.close()
            comm='gzip '+filename
            os.system(comm)
        mds.close()
    return 

def post_data(shot,ttime=0.,dt=0.1,dirs=''):
    import numpy as np
    import pickle
    from scipy.interpolate import interp1d
    from matplotlib import gridspec
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    filename = 'DATASAVE/%i'%shot+'/REF.npz'
    comm='gzip -d '+ filename+'.gz'
    os.system(comm)
    f=open(filename,'rb')
    prof = pickle.load(f)
    f.close()
    comm='gzip '+ filename
    os.system(comm)
    fig = plt.figure("Edge Reflec Time",figsize=(6,6))
    gs = gridspec.GridSpec(4,1)
    axes = [0,0]
    axes[0] = plt.subplot(gs[0])
    axes[0].set_title('VCO-Q band reflectometer')
    axes[0].set_xlabel('time [s]')
    axes[0].set_ylabel('V [a.u]')
    axes[0].set_ylim(0,1)
    if (not prof['isref']):
       axes[0].text(0.5,0.5,'No data',color='r')
       axes[0].set_xlim(0,1)
       plt.show()
       return
    nomds = np.sum(prof['prof']['is']) == 0
    axes[0].plot([prof['time'][0]-0.1,prof['time'][0]],[0,0],color='gray')
    tinclude = False; tind = -1;
    for k in range(8):
        if prof['prof']['is'][k] == 1: 
           if (ttime >= prof['prof'][k]['t'][0] and ttime <= prof['prof'][k]['t'][-1]): 
               tinclude = True; tind = k;
           axes[0].plot([prof['time'][k],prof['time'][k],prof['time'][k]+0.02,prof['time'][k]+0.02,],[0,0.5,0.5,0.],color='lime')
        else: axes[0].plot([prof['time'][k],prof['time'][k],prof['time'][k]+0.02,prof['time'][k]+0.02,],[0,0.5,0.5,0.],color='gray')
        if k<7.: axes[0].plot([prof['time'][k],prof['time'][k+1]],[0,0],color='gray')
    axes[0].plot([prof['time'][7],prof['time'][7]+0.1],[0,0],color='gray') 
    axes[0].plot([ttime,ttime],[0.3,0.7],color='magenta',linestyle='--')
    axes[0].set_xlim(prof['time'][0]-0.1,prof['time'][7]+0.1)
    if (not tinclude and ttime >0.): axes[0].text(prof['time'][0],0.8,'Given time slice has no reflectometer data',color='r') 

    tmin = ttime - 0.5*dt
    tmax = ttime + 0.5*dt
 
    if tinclude:
      tmin = max(tmin,min(prof['prof'][tind]['t']))
      tmax = min(tmax,max(prof['prof'][tind]['t']))
      axes[1] = plt.subplot(gs[1:4])
      rmin = max(prof['prof'][tind]['r'][:,-1])+0.0002; rmax = min(prof['prof'][tind]['r'][:,0])-0.0002
      rr = np.linspace(rmin,rmax,50); tlen = len(prof['prof'][tind]['t'])
      nn = np.zeros((50,tlen))
      for k in range(tlen):
        nf = interp1d(prof['prof'][tind]['r'][k,:],prof['prof'][tind]['n'][k,:])
        nn[:,k] = nf(rr)
      A=axes[1].contourf(prof['prof'][tind]['t'],rr,nn/1.e19)
      axes[1].axvline(x=ttime,color='magenta',linestyle='--') 
      axes[1].axvline(x=tmin,color='orange',linestyle='--')
      axes[1].axvline(x=tmax,color='orange',linestyle='--')
      axes[1].set_ylabel('R[m]')
      axes[1].set_xlabel('time [s]')
      axes[1].set_title('Reflectometer $n_e$ [$10^{19}m^{-3}$]')
      axes[1].set_xlim(min(prof['prof'][tind]['t'])-0.001,max(prof['prof'][tind]['t'])+0.001)

      divider = make_axes_locatable(axes[1])
      cax = divider.append_axes('right', size='5%', pad=0.05)
      fig.colorbar(A,cax=cax)
    plt.tight_layout()
    plt.show()

    if tinclude:
     ind1 = np.where(prof['prof'][tind]['t']>=tmin)
     ind2 = np.where(prof['prof'][tind]['t'][ind1]<=tmax)
     f = open(dirs,'w')
     f.write('R[m]    Z[m]    VAL     ERR\n')
     for k in range(len(rr)):
        avg = np.sum(nn[k,:][ind1][ind2])/len(nn[k,:][ind1][ind2])/1.e18
        std = np.std(nn[k,:][ind1][ind2])/1.e18
        f.write('%9.6f %9.6f %9.6f %9.6f\n'%(rr[k],0.,avg,std))
     f.close()

    return

# python ts5.py shot time dt multiplier percentage DAselect TE/NE
if __name__=='__main__':
    import os,sys
    import numpy as np

    try: os.mkdir('DATASAVE')
    except: pass
    try: os.mkdir('DATASAVE/'+sys.argv[1])
    except: pass
    only_mds = 'n'
    if not os.path.isfile('DATASAVE/'+sys.argv[1]+'/REF.npz.gz'):
        gprofdat(int(float(sys.argv[1])),'kstar')
    try: only_mds = sys.argv[5]
    except: pass
    if only_mds == 'y': exit()
 #   for i in range(22916,23135): 
 #      try: os.mkdir('DATASAVE/%i'%i)
 #      except: pass
 #      gprofdat(i,'kstar')
        
    post_data(int(float(sys.argv[1])),float(sys.argv[2]),float(sys.argv[3]),sys.argv[4])
