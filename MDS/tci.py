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
            if False:
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

    NE_INT1='\\ne_inter01'
    NE_INT2='\\ne_inter02'
    NE_TCI1='\\ne_tci01'
    NE_TCI2='\\ne_tci02'
    NE_TCI3='\\ne_tci03'
    NE_TCI4='\\ne_tci04'
    NE_TCI5='\\ne_tci05'
    names = [NE_INT1,NE_INT2,NE_TCI1,NE_TCI2,NE_TCI3,NE_TCI4,NE_TCI5]
    names2= ['ne_int1','ne_int2','ne_tci1','ne_tci2','ne_tci3','ne_tci4','ne_tci5']
    ISTCI = [False,False,False,False,False,False,False]
    LETCI = [0,0,0,0,0,0,0]
# connection to MDSPlus data
    with MDS(server="172.17.100.200:8005") as mds:
        try:
            eq=mds.open(shot=shot, tree=treename)
        except: 
            print("Error #2")
        else:
            print(mds.alist)
            for k in range(7):
                dat = names[k]
                try:
                    ne=eq.get(dat).data()
                    time=eq.get('dim_of(%s)'%names[k]).data()
                    dt = (time[10] - time[9])*1000
                    dscale = int(2./dt)
                    scale = 1.
                    #if shot < 22222: 
                    #  if k==0: scale= 1./1.9
                    #  if k==1: scale= 1./2.75
                    ne = ne[0::dscale] * scale
                    time = time[0::dscale]
                    len_ne = len(ne); len_t = len(time);
                    if len_ne > len_t: ne = ne[0:len_t]
                    if len_ne < len_t: time  = time[0:len_ne] 
                    ISTCI[k] = True
                    LETCI[k] = len(time)
                    filename = 'DATASAVE/'+sys.argv[1]+'/'+names2[k]+'.npz'
                    np.savez(filename,time,ne)
                    comm='gzip '+filename
                    os.system(comm)
                except:print('%s not avail'%names2[k])

        f = open('DATASAVE/'+sys.argv[1]+'/TCI_size','w')
        for i in range(7): 
            if ISTCI[i]: f.write('1 %i\n'%LETCI[i])
            else: f.write('0 %i\n'%LETCI[i])
        f.close()
        mds.close()
    return

def post_data(shot,time,dt,dirs,noplot):
    import numpy as np

    NE_conv=[1.9,2.75,1.,1.,1.,1.,1.]
    isdata =[0,0,0,0,0,0,0] 
    size   =[0,0,0,0,0,0,0]
    name   =['INT1','INT2','TCI1','TCI2','TCI3','TCI4','TCI5']
#    L_weight=[4,4,7.23,5.51,4.60,3.81,2.48]
    L_weight=[1,1,1,1,1,1,1]
    names2 =['ne_int1','ne_int2','ne_tci1','ne_tci2','ne_tci3','ne_tci4','ne_tci5']
    tciavg =[0,0,0,0,0,0,0]
    tcisig =[1,1,1,1,1,1,1]
    axes   =[0,0,0,0,0,0,0]

    utime = (time+0.5*dt)/1000
    ltime = (time-0.5*dt)/1000

    f=open('DATASAVE/'+sys.argv[1]+'/TCI_size','r')
    for k in range(7):
        stmp=f.readline()
        isdata[k]=int(stmp.split()[0])
        size[k]=int(stmp.split()[1])
    f.close()  
    if not noplot:
        fig = plt.figure("TCI Time",figsize=(6,6))
        gs = gridspec.GridSpec(7, 1)
    tmin = 0; tmax = 0;
    for k in range(7):
        if not noplot:
            if k==0: axes[k] = plt.subplot(gs[k])
            else: axes[k] = plt.subplot(gs[k],sharex=axes[0])
        if isdata[k] == 0:continue

        filename = 'DATASAVE/'+sys.argv[1]+'/'+names2[k]+'.npz'
        comm='gzip -d '+filename+'.gz'
        os.system(comm)
        npzfile  = np.load(filename, mmap_mode='r')
        time=npzfile['arr_0']
        ne=npzfile['arr_1']
        comm='gzip '+filename
        os.system(comm)

        tmin = min(tmin,min(abs(time)))
        tmax = max(tmax,max(abs(time)))
        lent = len(time); lenn = len(ne); lenm = min(lent,lenn)
        time = time[:lenm]; ne = ne[:lenm]
        if not noplot: axes[k].plot(time,ne/NE_conv[k])

        ind1 = np.where(time>ltime)
        ind2 = np.where(time[ind1]<utime)
        tciavg[k] = np.sum(ne[ind1][ind2])/len(ne[ind1][ind2])/NE_conv[k]
        ind3 = np.where(abs(ne[ind1][ind2]/NE_conv[k]-tciavg[k]) < 0.5)
        if (len(ne[ind1][ind2][ind3])==0): continue
        tciavg[k] = np.sum(ne[ind1][ind2][ind3])/len(ne[ind1][ind2][ind3])/NE_conv[k]
        tcisig[k] = (np.max(ne[ind1][ind2][ind3]) - np.min(ne[ind1][ind2][ind3])) / NE_conv[k] / 4.
        tcisig[k] = np.std(ne[ind1][ind2][ind3]) / NE_conv[k] * 1.5 / L_weight[k]  #*4

        if tciavg[k]<=0.: tciavg[k]= 0.; tcisig[k] = 0.
    
    if not noplot:
        axes[0].set_title('INTERFEROMETER [$10^{19} m^{-3}$]')
        axes[-1].set_xlabel('time [s]')
        for k in range(6): 
            axes[k].axes.get_xaxis().set_visible(False)
        for k in range(7):
            if isdata[k] == 0: 
                axes[k].plot([tmin,tmax],[0.,0.])
                axes[k].text(0.5*(tmin+tmax),0.5,'No data',color='red')
                axes[k].set_ylim(0,1.0)
            else:axes[k].set_ylim(0)

            axes[k].set_xlim(tmin,tmax)
            axes[k].legend([name[k]],loc='upper right')
            axes[k].axvline(x=ltime,linestyle='--',color='gold')
            axes[k].axvline(x=utime,linestyle='--',color='gold')

        fig.tight_layout()
        plt.subplots_adjust(hspace=0.)
        plt.show(block=True)
    f = open(dirs,'w')
    for k in range(7):
        f.write('%i %4.2f %4.3f \n'%(isdata[k],tciavg[k],tcisig[k]))
    f.close()

    return

# python ts5.py shot time dt multiplier percentage DAselect TE/NE
if __name__=='__main__':
    import os,sys
    try: os.mkdir('DATASAVE')
    except: pass
    try: os.mkdir('DATASAVE/'+sys.argv[1])
    except: pass
    only_mds = 'n'
    try: only_mds = sys.argv[5]
    except: pass
    noplot = False
    if not os.path.isfile('DATASAVE/'+sys.argv[1]+'/TCI_size'):
        gprofdat(int(float(sys.argv[1])),'kstar')
    if only_mds == 'y': exit()
    if only_mds == 'px': noplot = True
    post_data(int(float(sys.argv[1])),int(float(sys.argv[2])),int(float(sys.argv[3])),sys.argv[4],noplot)
