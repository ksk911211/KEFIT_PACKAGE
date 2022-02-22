import ctypes as _C
from MDSplus import Connection
from MDSplus._mdsshr import _load_library, MdsException 
#from MDSplus._mdsshr import MdsshrException as MdsException

import Tkinter
#from Tkinter import *

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import numpy as np

#ConnectToMds=_load_library('MdsIpShr').ConnectToMds
#ConnectToMds.argtypes=[_C.c_char_p]
#DisconnectFromMds = _load_library('MdsIpShr').DisconnectFromMds
#DisconnectFromMds.argtypes = [_C.c_int]

offset_factor = 1.001

CES_RR = dict()
CES_RR[1] = [1.795,1.800,1.850,1.900,1.950,2.000,2.050,2.100,2.140,2.160,2.170,2.180,2.190,2.200,2.205,2.210,2.215,2.220,2.225,2.230,2.235,2.240,2.245,2.250,2.255,2.265,2.275,2.280,2.285,2.290,2.295,2.300] # 2011
CES_RR[2] = [1.800,1.850,1.900,1.950,2.000,2.050,2.100,2.150,2.170,2.180,2.190,2.200,2.205,2.210,2.215,2.220,2.225,2.230,2.235,2.240,2.245,2.250,2.255,2.260,2.265,2.270,2.275,2.280,2.285,2.290,2.295,2.300] # 2012
CES_RR[3] = [1.795,1.850,1.900,1.950,2.000,2.050,2.100,2.150,2.170,2.180,2.190,2.200,2.205,2.210,2.215,2.220,2.225,2.230,2.235,2.240,2.245,2.250,2.255,2.260,2.265,2.270,2.275,2.280,2.285,2.290,2.295,2.300] # 2013/2014
CES_RR[4] = [1.801,1.822,1.843,1.874,1.895,1.945,1.995,2.016,2.047,2.078,2.099,2.125,2.150,2.171,2.192,2.203,2.213,2.223,2.228,2.233,2.238,2.243,2.248,2.253,2.259,2.264,2.269,2.273,2.280,2.286,2.291,2.296] # 2015

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
                raise Exception, "Error in disconnection"
             else:
                self.socket = -1
                
    def reconnect(self):
        if self.hostspec == None:
             raise MdsException, "Error: no host specified"
        else:
             if self.socket != -1:
                print self.socket
                try:
                    self.closeConnection()
                except:
                    raise Exception, "Error in resetting connection to %s" %(self.hostspec,)
             self.socket = ConnectToMds(self.hostspec)
             if self.socket == -1:
                raise Exception, "Error connecting to %s" %(self.hostspec,)   

class MDS(object):
    """
    Implementation of a connection to the MDSplus tree based on the class mds by Y. M. Jeon
    Written by D. K. Oh
    Last modification : Aug 2012
    """
    __DefaultTree = "KSTAR"
#    __DefaultServer = "172.17.250.21:8005"
    __DefaultServer = "172.17.100.200:8005"

    def __init__(self, shot=None, tree =__DefaultTree, server=__DefaultServer):
        try:                    
            self.alist = {"tree":None, "shot":None, "server":None}
            self.__mds__ = _Connection(server)
            if shot is not None:
                self.open(shot, tree)
        except MdsException:
            raise MdsException, "Error in the connection %s" %(server)
        except:
            raise Exception, " Unknown error in the connection %s" %(server)
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
            raise MdsException, "Error in open : shot number is not specified"
        else:    
            self.close( self.alist["shot"], self.alist["tree"])
            try:
                self.__mds__.openTree(tree, shot)
            except:
                self.alist["tree"] = None
                self.alist["shot"] = None
                raise MdsException, "Error in open : unknown error"    
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
                raise MdsException, "Error in close : unknown error"
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
            raise MdsException, "Error in get"
        return ret_str

    def get_sig(self, sigstr):
        try:
            t = self.__mds__.get('dim_of(%s)' %(sigstr)).data();
            v = self.__mds__.get(sigstr).data();
        except:
            t = numpy.ndarray([])
            v = numpy.ndarray([])             
            raise MdsException, "Error in get"
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

# select times, only if not already in the 2nd list
def _select(event):
    tmp='%.5f'%tabtime[int(lb1.curselection()[0])]
    liste=lb2.get(0,Tkinter.END)
    if tmp not  in liste:
        lb2.insert(Tkinter.END,float(tmp))
    _update2(-1.0)
    _replotDA(-1.0)

def _selectall():
    liste1=lb1.get(0,Tkinter.END)
    lb2.delete(0,Tkinter.END)
    for i in liste1:
        lb2.insert(Tkinter.END,i)
    _update2(-1.0)
    _replotDA(-1.0)

def _preselect():
    liste1=lb1.get(0,Tkinter.END)
    lb2.delete(0,Tkinter.END)

#    print nbpick
    for i in liste1:
        ok=1
        for j in range(nbpick):
            if pdeb[j] < float(i) and float(i) < pfin[j]:
                ok=0

        if ok == 1:
            lb2.insert(Tkinter.END,i)



    _update2(-1.0)
    _replotDA(-1.0)    

# plot profile contains in 2nd list
def _update2(temps):
#    fig = Figure(figsize=(5, 4), dpi=100)
    fig.clf()
    liste=lb2.get(0,Tkinter.END)
    rmin=float(ermint.get())
    rmax=float(ermaxt.get())
    ymin=float(eymint.get())
    ymax=float(eymaxt.get())
    a=fig.add_subplot(111,xlim=(rmin,rmax),ylim=(ymin,ymax))
    a.cla()
    a.set_ylim(ymin=ymin,ymax=ymax)
    a.set_xlim(xmin=rmin,xmax=rmax)
    cptj=-1
    for i in liste:
        for j in range(taille):
            value='%.5f'%tabtime[j]
            if i == float(value):
                couleur='blue'
                if i==temps:
                    cptj=j
                else:
#                    a.scatter(rCES[j][0:taille3-1], tCES[j][0:taille3-1],c=couleur)
                    a.errorbar(rCES[j][0:taille3-1], tCES[j][0:taille3-1],yerr=tCESerr[j][0:taille3-1],c=couleur,fmt='.')

    if cptj >=0: print tCES[cptj][0]
    if cptj >=0:
        couleur='red'
#        a.scatter(rCES[cptj][0:taille3-1], tCES[cptj][0:taille3-1],c=couleur)
        a.errorbar(rCES[cptj][0:taille3-1], tCES[cptj][0:taille3-1],yerr=tCESerr[cptj][0:taille3-1],c=couleur,fmt='o')
    if choix==0:
        fig.suptitle('TI [keV]')
    else:
        fig.suptitle('VT [km/s]')

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().place(y=350,x=250)
    canvas.draw()

def _plotred(event):
    liste=lb2.get(0,Tkinter.END)
    temps=liste[lb2.curselection()[0]]
#    print temps
    _update2(temps)
    _replotDA(temps)

# replot Dalpha
def _replotDA(temps):
    tmin=float(etmint.get())
    tmax=float(etmaxt.get())
#    etmint.set(str(max(tmin,tALPHA[0])))
#    etmaxt.set(str(min(tmax,tALPHA[taille2-1])))
#    tmin=float(etmint.get())
#    tmax=float(etmaxt.get())
#    fig2 = Figure(figsize=(5, 3), dpi=100)
    fig2.clf()
    if dasel ==1:
        fig2.add_subplot(111,xlim=(tmin,tmax)).plot(tALPHA, DALPHA,color='black')

    xvar=[]
    yvar=[]
    liste=lb2.get(0,Tkinter.END)
    for i in liste:
        xvar.append(i)
        yvar.append(0.)

    fig2.add_subplot(111,xlim=(tmin,tmax)).scatter(xvar, yvar,color='blue',marker='.')

    xvarr=temps
    yvarr=0.
    fig2.add_subplot(111,xlim=(tmin,tmax)).scatter(xvarr, yvarr,color='red',marker='+')

    if dasel ==1:
        yvarp=np.zeros(nbpick)
#        fig2.add_subplot(111,xlim=(tmin,tmax)).scatter(timepeak, yvarp,color='green',marker='+')
        for i in range(len(timepeak)):
           fig2.add_subplot(111,xlim=(tmin,tmax)).axvline(x=timepeak[i],color='green',linestyle='--',linewidth=0.8)
    fig2.suptitle('$D_{\\alpha}$')
    canvas2 = FigureCanvasTkAgg(fig2, master=root)  # A tk.DrawingArea.
    canvas2.draw()
    canvas2.get_tk_widget().place(y=30,x=250)
    canvas2.draw()

def _replotR():
#    rmin=float(ermint.get())
#    rmax=float(ermaxt.get())
#    ermint.set(str(max(rmin,rCES[0][0])))
#    ermaxt.set(str(min(rmax,rCES[0][taille3-1])))

#    ymin=float(eymint.get())
#    ymax=float(eymaxt.get())
#    eymint.set(str(max(ymin,0.)))
#    eymaxt.set(str(min(ymax,np.amax(tCES))))

    _update2(-1.)
#    fig2 = Figure(figsize=(5, 3), dpi=100)
#    fig2.add_subplot(111,xlim=(tmin,tmax)).plot(tALPHA[0:taille2-1], DALPHA[0:taille2-1],color='black')
#    canvas2 = FigureCanvasTkAgg(fig2, master=root)  # A tk.DrawingArea.
#    canvas2.draw()
#    canvas2.get_tk_widget().place(y=30,x=250)
#    canvas2.draw()

# remove profiles from the 2nd list
def _remove():
    lb2.delete(Tkinter.ANCHOR)
    _update2(-1.)
    _replotDA(-1)

def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL state

def _f1():
    fig.add_subplot(111).plot(t2, 2 * np.cos(2 * np.pi * t2))
    canvas.draw()

def _f2():
    fig.add_subplot(111).plot(t2, 2 * np.cos(np.pi * t2))
    canvas.draw()

def _todo():
    print 'To do...'

# print profiles contains in the 2nd list only
def _generate2():
    if choix==0:
        fname='ti_%d_%d.dat'%(shot,mtime*1000)
        f=open(fname,'w')
        f.write(' R[m]	Z[m]	TI[eV]	Error[eV]\n')
        liste=lb2.get(0,Tkinter.END)

        for i in liste:
            for j in range(taille):
                if i == float('%.5f'%tabtime[j]):
                    for k in range(taille3):
                        f.write(' %.3f  0.000    %.3f   %.3f\n'%(rCES[j][k],1000.*tCES[j][k],1000.*tCESerr[j][k]))
        f.close()
    else:
        fname='vt_%d_%d.dat'%(shot,mtime*1000)
        f=open(fname,'w')
        f.write(' R[m]	Z[m]	VT[km/s]	Error[km/s]\n')
        liste=lb2.get(0,Tkinter.END)

        for i in liste:
            for j in range(taille):
                if i == float('%.5f'%tabtime[j]):
                    for k in range(taille3):
                        f.write(' %.3f  0.000    %.3f   %.3f\n'%(rCES[j][k],tCES[j][k],tCESerr[j][k]))
        f.close()

    print fname,' generated ',taille,taille3

def gprofdat(shot,treename,temps,dt):
    import numpy as np

    tail=-1
    tail2=-1
    DA='\\TOR_HA10' # Dalpha
    CES_RT='\\CES_RT' # CES radius
    CES_C=['\\CES_TI','\\CES_VT'] # CES TI
    CES_conv=[1000.,1.]

# connection to MDSPlus data
#    with MDS(server="172.17.250.21:8005") as mds:
    with MDS(server="172.17.100.200:8005") as mds:
        try:
            eq=mds.open(shot=shot, tree=treename)
        except: 
            print "Error #2"
        else:
            try:
                print mds.alist
                err=0
# read CES data
                time=eq.get('dim_of(%s)'%'\\CES_TI01').data() # time
                size=len(time)

                let=-1
                for i in range(size-1):
                    if (let==-1) and (time[i] <= temps) and (temps < time[i+1]):
                        let=i

                if size <2:
                    err=1
                else:

# loop on probe (radial)
                    for j in range(nb_ces):
                        dat=CES_RT+str(j+1).zfill(2)
                        if shot>12388: radius=eq.get(dat).data()
                        elif shot > 11726: radius = CES_RR[4][j] * np.ones(size) * 1.e3
                        elif shot > 8355: radius  = CES_RR[3][j] * np.ones(size) * 1.e3
                        elif shot > 6471: radius  = CES_RR[2][j] * np.ones(size) * 1.e3
                        else: radius = CES_RR[1][j] * np.ones(size) * 1.e3
                        dat=CES_C[choix]+str(j+1).zfill(2)
                        ti=eq.get(dat).data()
#                        dat=CES_VT+str(j+1).zfill(2)
#                        vt=eq.get(dat).data()

                        dat=CES_C[choix]+str(j+1).zfill(2)+':err_bar'
                        try:
                            tierr=eq.get(dat).data()
                        except:
                            tierr=[x / 10. for x in ti]
                            print 'No CES err_bar => 10% ',CES_C[choix]+str(j+1)


                        tailj=0
# loop on time
                        for i in range(size):
                            if (temps-dt/2. <= time[i]) and (time[i] < temps+dt/2.):
                                tabtime[tailj]=time[i]
                                rCES[tailj][j]=radius[i]/1000.
                                tCES[tailj][j]=ti[i]/CES_conv[choix] #keV
                                tCESerr[tailj][j]=tierr[i]/CES_conv[choix] #keV
                                tailj=tailj+1

# remove NaN 
                    for i in range(tailj):  
                        tCES[i]=np.nan_to_num(tCES[i])
                        tCESerr[i]=np.nan_to_num(tCESerr[i])

            except:
                print "Error CES"
            else:
                print "End of CES reading"

            print 'verif', nb_ces

            err=0
                # read Dalpha data
            time=eq.get('dim_of(%s)'%DA).data()               
            size=len(time)
#                print size,tabtime[0],tabtime[tail-1],time
            if size <2:
                err=1
            else:
                tail2=0

                dat=DA
                Dalpha=eq.get(dat).data()

                alpha_offset = 0.
                if (Dalpha[0] < 0.):alpha_offset = -offset_factor * Dalpha[0]
		Dalpha = Dalpha + alpha_offset
                for i in range(size):
                    if (tabtime[0] <= time[i] ) and (time[i]<= tabtime[tailj-1]):
                        tALPHA[tail2]=time[i]
                        DALPHA[tail2]=Dalpha[i]
                        tail2=tail2+1
                        

#    print tail
    return (tailj,tail2,nb_ces)

    mds.close()


def gprofdatdat(shot,treename,temps,dt):
    import numpy as np

    tail=-1
    tail2=-1
    CES_conv=[1000.,1.]

    f=open(savfile+'_size','r')
    stmp=f.readline()
    size=int(stmp.split()[0])
    size2=int(stmp.split()[1])
    nbprobe=int(stmp.split()[2])
    f.close()    

    time=np.zeros(size)
    ti=np.zeros(size)
    tierr=np.zeros(size)
    for k in range(nbprobe):
        fichier=savfile+str(k+1)+'.npz'
# unzip
        comm='gzip -d '+fichier+'.gz'
        os.system(comm)

        npzfile  = np.load(fichier, mmap_mode='r')
        time=npzfile['arr_0']
        ti=npzfile['arr_1']
        tierr=npzfile['arr_2']
        radius=npzfile['arr_3']

# zip again
        comm='gzip '+fichier
        os.system(comm)

        tailj=0 # time
# loop on time
        for j in range(size):
            if (temps-dt/2. <= time[j]) and (time[j] < temps+dt/2.):
                tabtime[tailj]=time[j]
                rCES[tailj][k]=radius/1000.
                tCES[tailj][k]=ti[j]/CES_conv[choix] #keV
                tCESerr[tailj][k]=tierr[j]/CES_conv[choix] #keV
                tailj=tailj+1

        for i in range(tailj):
            tCES[i]=np.nan_to_num(tCES[i])
            tCESerr[i]=np.nan_to_num(tCESerr[i])  

#DA                            
    time=np.zeros(size2)
    Dalpha=np.zeros(size2)
    istart=0
# unzip
    if dasel ==1:
        comm='gzip -d '+dafile+'.gz'
        os.system(comm)

#    f=open(dafile,'r')
#    for i in range(size2):
#        stmp=f.readline()
#        time[i]=float(stmp.split()[0])
#        Dalpha[i]=float(stmp.split()[1])
#    f.close()
        npzfile  = np.load(dafile, mmap_mode='r')
        time=npzfile['arr_0']
        Dalpha=npzfile['arr_1']
        alpha_offset = 0.
        if (Dalpha[0] < 0.):alpha_offset = -offset_factor * Dalpha[0]
        Dalpha = Dalpha + alpha_offset


# zip again
        comm='gzip '+dafile
        os.system(comm)

        tail2=0
        istart=-1
        for i in range(size2):
            if (tabtime[0] <= time[i] ) and (time[i]<= tabtime[tailj-1]):
                if istart ==-1:
                    istart=i
#            tALPHA[tail2]=time[i]
#            DALPHA[tail2]=Dalpha[i]
                tail2=tail2+1
    else:
        tail2=0

    return (tailj,tail2,nbprobe,istart,time,Dalpha)

def gprofdatmds(shot,treename,temps,dt):
    import numpy as np

    tail=-1
    tail2=-1
    DA='\\TOR_HA10' # Dalpha
    CES_RT='\\CES_RT' # CES radius
    CES_C=['\\CES_TI','\\CES_VT'] # CES TI
    CES_conv=[1000.,1.]


# connection to MDSPlus data
#    with MDS(server="172.17.250.21:8005") as mds:
    with MDS(server="172.17.100.200:8005") as mds:
        try:
            eq=mds.open(shot=shot, tree=treename)
        except: 
            print "Error #2"
        else:
            try:
                print mds.alist

                err=0
# read CES data
                time=eq.get('dim_of(%s)'%'\\CES_TI01').data() # time
                size=len(time)
                
                if size <2:
                    err=1
                else:
# loop on probe (radial)
                    for j in range(nb_ces):
                        dat=CES_RT+str(j+1).zfill(2)
#                        print dat
                        if  shot  > 12388: radius=eq.get(dat).data()
                        elif shot > 11726: radius = CES_RR[4][j] * np.ones(size) * 1.e3
                        elif shot >  8355: radius = CES_RR[3][j] * np.ones(size) * 1.e3
                        elif shot >  6471: radius = CES_RR[2][j] * np.ones(size) * 1.e3
                        else:              radius = CES_RR[1][j] * np.ones(size) * 1.e3

                        dat=CES_C[choix]+str(j+1).zfill(2)
                        ti=eq.get(dat).data()
#                        dat=CES_VT+str(j+1).zfill(2)
#                        vt=eq.get(dat).data()

                        dat=CES_C[choix]+str(j+1).zfill(2)+':err_bar'
                        try:
                            tierr=eq.get(dat).data()
                        except:
                            tierr=[x / 10. for x in ti]
                            print 'No CES err_bar => 10% ',CES_C[choix]+str(j+1)

                        tailj=0
# loop on time
                        for i in range(size):
                            if (temps-dt/2. <= time[i]) and (time[i] < temps+dt/2.):
                                tabtime[tailj]=time[i]
                                rCES[tailj][j]=radius[i]/1000.
                                tCES[tailj][j]=ti[i]/CES_conv[choix] #keV
                                tCESerr[tailj][j]=tierr[i]/CES_conv[choix] #keV
                                tailj=tailj+1

                        fichier=savfile+str(j+1)+'.npz'
                        np.savez(fichier,time,ti,tierr,radius[0])

                        comm='gzip '+fichier
                        os.system(comm)

# remove NaN 
                    for i in range(tailj):  
                        tCES[i]=np.nan_to_num(tCES[i])
                        tCESerr[i]=np.nan_to_num(tCESerr[i])

            except:
                print "Error CES"
            else:
                print "End of CES reading"

            err=0
                # read Dalpha data
            time=eq.get('dim_of(%s)'%DA).data()               
            size2=len(time)
#                print size,tabtime[0],tabtime[tail-1],time
            if size2 <2:
                err=1
            else:
                tail2=0

                dat=DA
                Dalpha=eq.get(dat).data()
                alpha_offset = 0.
                if (Dalpha[0] < 0.):alpha_offset = -offset_factor * Dalpha[0]
                Dalpha = Dalpha + alpha_offset

                tstart=-1
                for i in range(size2):
                    if (tabtime[0] <= time[i] ) and (time[i]<= tabtime[tailj-1]):
                        if tstart==-1:
                            tstart=i
#                        tALPHA[tail2]=time[i]
#                        DALPHA[tail2]=Dalpha[i]
                        tail2=tail2+1
                        
                if not os.path.isfile(dafile+'.gz'):
#                    f=open(dafile,'w')
#                    for i in range(size2):
#                        f.write('%f %f\n'%(time[i],Dalpha[i]))
#                    f.close()
                    np.savez(dafile,time,Dalpha)
                    comm='gzip '+dafile
                    os.system(comm)

            
            f=open(savfile+'_size','w')
            f.write('%d %d %d'%(size,size2,nb_ces))
            f.close()
#    print tail
    mds.close()

    return (tailj,tail2,nb_ces,tstart,time,Dalpha)


# python ti2.py shot time dt
if __name__=='__main__':
   """ """
   import os, sys

   exename=os.path.basename(__file__);
   nargs=len(sys.argv);
   
#   treename = 'HEATING'
   treename='kstar'
   shot=int(sys.argv[1])
   mtime=float(sys.argv[2])
   dt=float(sys.argv[3])
   meanmult=float(sys.argv[4])
   pourcent=float(sys.argv[5])
   dasel=int(sys.argv[6]) # 1 =DA, 0 = no DA
   choix=int(sys.argv[7]) # 0=ti, 1=vt
   try: mds_only = sys.argv[8]
   except: mds_only = 'n'

   root = Tkinter.Tk()
   if choix==0:
       root.wm_title("TI")
   else:
       root.wm_title("VT")
   root.geometry("900x820")
#   root.pack(fill=BOTH, expand=1)

   fig = Figure(figsize=(5, 4), dpi=100)
#   t = np.arange(0, 3, .01)
#fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))

   nb_ces=32
   sizmax=1000
   shape=(sizmax,nb_ces)
   tabtime=np.zeros(sizmax)
   rCES=np.zeros(shape)
   tCES=np.zeros(shape)
   tCESerr=np.zeros(shape)
#   tALPHA=np.zeros(sizmax)
#   DALPHA=np.zeros(sizmax)

   if choix==0:
       savfile='./DATASAVE/'+str(shot)+'/TI'
   else:
       savfile='./DATASAVE/'+str(shot)+'/VT'
   dafile='./DATASAVE/'+str(shot)+'/DA.npz'
# data already saved
   if os.path.isfile(savfile+'1.npz.gz'):
       print 'Data read from DATASAVE'
       (taille,taille2,taille3,istart,tALPHA,DALPHA)=gprofdatdat(shot,treename,mtime,dt)
   else:
       print 'Data read from KSTAR MDSPlus tree'
       (taille,taille2,taille3,istart,tALPHA,DALPHA)=gprofdatmds(shot,treename,mtime,dt)
       if mds_only == 'y': exit()
# 1st time...need data from MDSPlus Tree
#       (taille,taille2,taille3,istart,tALPHA,DALPHA)=gprofdatmds(shot,treename,mtime,dt)
# taille = size time CES, taille2= size Dalpha, taille3 = size radial CES
   print taille,taille2,taille3,istart
   if os.path.isfile('./DATASAVE/CES_err'): os.system('rm ./DATASAVE/CES_err')
   if not os.path.isfile(savfile+'1.npz.gz'): 
     f = open('./DATASAVE/CES_err','w')
     f.write('sad')
     f.close()
   
   fig2 = Figure(figsize=(5, 3), dpi=100)
   if dasel ==1:
       fig2.add_subplot(111,xlim=(tALPHA[istart],tALPHA[istart+taille2-1])).plot(tALPHA, DALPHA,color='black')
   

# for Pre select.
# select only time points that are outside 30% zone after a DA peak
# looking for peaks
   inred=0
   nbpick=0
   pdeb=[]
   pfin=[]
   timepeak=[]
   if dasel ==1:
       moyenne=np.mean(DALPHA[istart:istart+taille2-1])
       for i in range(taille2):
           if DALPHA[istart+i] > meanmult*moyenne and inred==0:
               if len(timepeak) == 0 or (tALPHA[istart+i]-timepeak[-1]) > 0.002:
                   timepeak.append(tALPHA[istart+i])
                   nbpick=nbpick+1
                   inred=1
           if DALPHA[istart+i] < meanmult*moyenne and inred==1:
               inred=0
           
       print nbpick,'peaks (in green)'
       if nbpick == 0:
           print 'no peak found'
           dasel=0
       else:
           for i in range(nbpick-1):
               pdeb.append(timepeak[i]-0.001)
               pfin.append(timepeak[i]+(timepeak[i+1]-timepeak[i])*pourcent/100.)
           pdeb.append(timepeak[-1]-0.001)
           if nbpick > 1:
               pfin.append(timepeak[-1]+(timepeak[-2]-timepeak[-1])*pourcent/100.)
           else:
               pfin.append(timepeak[-1]-0.05)

   l1=Tkinter.Label(root, text="Shot:%d"%shot)
   l1.place(y=10,x=10)

   if dasel ==1:
       l2=Tkinter.Label(root, text="Dalpha")
   else: 
       l2=Tkinter.Label(root, text="No Dalpha peak")
   l2.place(y=10, x=450)

# upper listbox
   frame1=Tkinter.Frame(root)
   frame1.place(y=40,x=10)
   
   lb1=Tkinter.Listbox(frame1)
   lb1.bind("<Double-Button-1>",_select)
   lb1.pack(side = 'left',fill = 'y' )

   sc1=Tkinter.Scrollbar(frame1,orient=Tkinter.VERTICAL,command=lb1.yview)
   sc1.pack(side="right", fill="y")
   lb1.config(yscrollcommand=sc1.set)

   for i in range(taille):
       value='%.5f'%tabtime[i]
#       value=str(tabtime[i])
       lb1.insert(Tkinter.END,float(value))


   buttonS = Tkinter.Button(master=root, text="Select 1", command=lambda:_select(0))
   buttonS.place(y=210,x=10)

   buttonSA = Tkinter.Button(master=root, text="Select All", command=_selectall)
   buttonSA.place(y=210,x=90)
   buttonSp = Tkinter.Button(master=root, text="Pre Select", command=_preselect)
   buttonSp.place(y=240,x=90)
   if dasel ==0:
       buttonSp.config(state='disabled')

#lower listbox
   frame2=Tkinter.Frame(root)
   frame2.place(y=290,x=10)

   lb2=Tkinter.Listbox(frame2)
   lb2.pack(side = 'left',fill = 'y' )
   sc2=Tkinter.Scrollbar(frame2,orient=Tkinter.VERTICAL,command=lb2.yview)
   sc2.pack(side="right", fill="y")
   lb2.config(yscrollcommand=sc2.set)
   lb2.bind("<Double-Button-1>",_plotred)

   buttonR = Tkinter.Button(master=root, text="Remove", command=_remove)
   buttonR.place(y=460,x=60)

# upper figure (Dalpha)
   canvas2 = FigureCanvasTkAgg(fig2, master=root)  # A tk.DrawingArea.
   canvas2.draw()
   canvas2.get_tk_widget().place(y=30,x=250)

# lower figure (profiles)   
   canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
   canvas.draw()
   canvas.get_tk_widget().place(y=350,x=250)
   
   ltmin=Tkinter.Label(root, text="Tmin")
   ltmin.place(y=250, x=770)
   etmint=Tkinter.StringVar()
   etmin = Tkinter.Entry(root,width=10,textvariable=etmint)
   etmin.place(y=250, x=810)
   if dasel==1:
       value='%.5f'%tALPHA[istart]
   else:
       value='%.5f'%tabtime[0]
   etmint.set(value)

   ltmax=Tkinter.Label(root, text="Tmax")
   ltmax.place(y=270, x=770)
   etmaxt=Tkinter.StringVar()
   etmax = Tkinter.Entry(root,width=10,textvariable=etmaxt)
   etmax.place(y=270, x=810)
   if dasel==1:
       value='%.5f'%tALPHA[istart+taille2-1]
   else:
       value='%.5f'%tabtime[taille-1]
   etmaxt.set(value)

   buttonRDA = Tkinter.Button(master=root, text="Replot", command=lambda:_replotDA(-1))
   buttonRDA.place(y=290, x=790)

   if choix==0:
       lymin=Tkinter.Label(root, text="TImin")
   else:
       lymin=Tkinter.Label(root, text="VTmin")
   lymin.place(y=600, x=770)
   eymint=Tkinter.StringVar()
   eymin = Tkinter.Entry(root,width=10,textvariable=eymint)
   eymin.place(y=600, x=810)
   if choix==0:
       value="0."
   else:
       value='%.5f'%np.amin(tCES)
   eymint.set(value)

   if choix==0:
       lymax=Tkinter.Label(root, text="TImax")
   else:
       lymax=Tkinter.Label(root, text="VTmax")
   lymax.place(y=620, x=770)
   eymaxt=Tkinter.StringVar()
   eymax = Tkinter.Entry(root,width=10,textvariable=eymaxt)
   eymax.place(y=620, x=810)

   value='%.5f'%np.amax(tCES)
   eymaxt.set(value)

   lrmin=Tkinter.Label(root, text="Rmin")
   lrmin.place(y=680, x=770)
   ermint=Tkinter.StringVar()
   ermin = Tkinter.Entry(root,width=10,textvariable=ermint)
   ermin.place(y=680, x=810)
   ermint.set(rCES[0][0])

   lrmax=Tkinter.Label(root, text="Rmax")
   lrmax.place(y=700, x=770)
   ermaxt=Tkinter.StringVar()
   ermax = Tkinter.Entry(root,width=10,textvariable=ermaxt)
   ermax.place(y=700, x=810)
   ermaxt.set(rCES[0][taille3-1])

   buttonRR = Tkinter.Button(master=root, text="Replot", command=_replotR)
   buttonRR.place(y=720, x=790)

   buttonG = Tkinter.Button(master=root, text="Generate", command=_generate2)
   buttonG.place(y=770,x=420)

   buttonQ = Tkinter.Button(master=root, text="Quit", command=_quit)
   buttonQ.place(y=770,x=520)

   Tkinter.mainloop()
# If you put root.destroy() here, it will cause an error if the window is
# closed with the window manager.
