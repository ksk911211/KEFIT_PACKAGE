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

CORE_R = [1806, 1826, 1848, 1871, 1894, 1917, 1942, 1966, 1991, 2016, 2041, 2068, 2093, 2120]
EDGE_R = [2124, 2137, 2143, 2149, 2156, 2162, 2177, 2191, 2202, 2216, 2229, 2242, 2257, 2271, 2285, 2297, 2311]

offset_factor = 1.001
shot_c = 21760
shot_d = 24136

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
#                    fig.add_subplot(111,xlim=(rmin,rmax),ylim=(ymin,ymax)).scatter(rTS[j][0:taille3-1], tTS[j][0:taille3-1],c=couleur)
                    a.errorbar(rTS[j][0:taille3-1], tTS[j][0:taille3-1],yerr=tTSerr[j][0:taille3-1],c=couleur,fmt='.')

    if cptj >=0:
        couleur='red'
#        fig.add_subplot(111,xlim=(rmin,rmax),ylim=(ymin,ymax)).scatter(rTS[cptj][0:taille3-1], tTS[cptj][0:taille3-1],c=couleur)
        a.errorbar(rTS[cptj][0:taille3-1], tTS[cptj][0:taille3-1],yerr=tTSerr[cptj][0:taille3-1],c=couleur,fmt='o')
    if choix==0:
        fig.suptitle('TE [keV]')
    else:
        fig.suptitle('NE [$10^{19}$$m^{-3}$]')
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

#    fig2.cla()
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

    if dasel==1:
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
#    ermint.set(str(max(rmin,rTS[0][0])))
#    ermaxt.set(str(min(rmax,rTS[0][taille3-1])))

#    ymin=float(eymint.get())
#    ymax=float(eymaxt.get())
#    eymint.set(str(max(ymin,0.)))
#    eymaxt.set(str(min(ymax,np.amax(tTS))))

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
    global ncore
    print(ncore)
    if choix ==0:
        fname1='te_%d_%d.dat'%(shot,mtime*1000)
        fname2='te_%d_%d_edge.dat'%(shot,mtime*1000)
        f1=open(fname1,'w'); f2=open(fname2,'w');
        f1.write(' R[m]	Z[m]	TE[eV]	Error[eV]\n')
        f2.write(' R[m]  Z[m]    TE[eV]  Error[eV]\n')
        liste=lb2.get(0,Tkinter.END)

        for i in liste:
            for j in range(taille):
                if i == float('%.5f'%tabtime[j]):
                    for k in range(ncore):
                        f1.write(' %.3f  0.000    %.3f   %.3f\n'%(rTS[j][k],1000.*tTS[j][k],1000.*tTSerr[j][k]))
                    for k in range(ncore,taille3):
                        f2.write(' %.3f  0.000    %.3f   %.3f\n'%(rTS[j][k],1000.*tTS[j][k],1000.*tTSerr[j][k]))
        f1.close(); f2.close();
    else:
        fname1='ne_%d_%d.dat'%(shot,mtime*1000)
        fname2='ne_%d_%d_edge.dat'%(shot,mtime*1000)
        f1=open(fname1,'w'); f2=open(fname2,'w');
        f1.write(' R[m]	Z[m]	NE[#E18/m3]	Error[#E18/m3]\n')
        f2.write(' R[m] Z[m]    NE[#E18/m3]     Error[#E18/m3]\n')
        liste=lb2.get(0,Tkinter.END)

        for i in liste:
            for j in range(taille):
                if i == float('%.5f'%tabtime[j]):
                    for k in range(ncore):
                        f1.write(' %.3f  0.000    %.3f   %.3f\n'%(rTS[j][k],10.*tTS[j][k],10.*tTSerr[j][k]))
                    for k in range(ncore,taille3):
                        f2.write(' %.3f  0.000    %.3f   %.3f\n'%(rTS[j][k],10.*tTS[j][k],10.*tTSerr[j][k]))
                    
        f1.close(); f2.close();

    print fname1,fname2,' generated'

def gprofdat(shot,treename,temps,dt):
    import numpy as np
    ncore=1
    tail=-1
    tail2=-1
    DA='\\TOR_HA10' # Dalpha
    TS_CORE='\\TS_CORE'
    TS_EDGE='\\TS_EDGE'
    TS_C=['_TE','_NE']
    TS_CE=['_TERRH','_NERRH']
    TS_conv=[1000.,1.e19]
    

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
# read TS data
                time=eq.get('dim_of(%s)'%'\\TS_CORE2.CORE2_TE').data()
                size=len(time)

                let=-1
                for i in range(size-1):
                    if (let==-1) and (time[i] <= temps) and (temps < time[i+1]):
                        let=i

                if size <2:
                    err=1
                else:
                    taili=0 # radial
# loop on probe (radial)
                    for i in range(nb_ts):
                        dat=TS_CORE+str(i+1)+'.CORE'+str(i+1)+TS_C[choix]
                        try:
                            te=eq.get(dat).data(); ncore=i+1
                            print 'TS core avail  ',i+1
                        except:
                            print 'TS core unavail',i+1
                        else:

                            dat=TS_CORE+str(i+1)+'.CORE'+str(i+1)+TS_CE[choix]
                            try:
                                teerr=eq.get(dat).data()
                            except:
                                teerr=[x / 10. for x in te]
                                print 'No TS err_bar => 10% ',TS_CORE+str(i+1)

                            if te[let] > 0.:
                                dat=TS_CORE+str(i+1)+'.CORE'+str(i+1)+'_POS'
                                if ((float(shot)>=shot_c) and (float(shot)<shot_d)): radius = CORE_R[i]; print 'prescribed R=%3.2f'%radius;
                                else: radius=eq.get(dat).data(); print 'MDS R=%3.2f'%radius
                                tailj=0 # time
# loop on time
                                for j in range(size):
                                    if (temps-dt/2. <= time[j]) and (time[j] < temps+dt/2.):
                                        tabtime[tailj]=time[j]
                                        rTS[tailj][taili]=radius/1000.
                                        tTS[tailj][taili]=te[j]/TS_conv[choix] #keV
                                        tTSerr[tailj][taili]=teerr[j]/TS_conv[choix] #keV
                                        tailj=tailj+1
                                taili=taili+1

                    taili2=0 # radial
# loop on probe (radial)
                    for i in range(nb_ts):
                        dat=TS_EDGE+str(i+1)+'.EDGE'+str(i+1)+TS_C[choix]
                        try:
                            te=eq.get(dat).data()
                            print 'TS edge avail  ',i+1
                        except:
                            print 'TS edge unavail',i+1
                        else:

                            dat=TS_EDGE+str(i+1)+'.EDGE'+str(i+1)+TS_CE[choix]
                            try:
                                teerr=eq.get(dat).data()
                            except:
                                teerr=[x / 10. for x in te]
                                print 'No TS err_bar => 10% ',TS_EDGE+str(i+1)


                            if te[let] > 0.:
                                dat=TS_EDGE+str(i+1)+'.EDGE'+str(i+1)+'_POS'
                                if ((float(shot)>=shot_c) and (float(shot)<shot_d)): radius= EDGE_R[i]; print 'prescribed R=%3.2f'%radius;
                                else:
                                     radius = eq.get(dat).dat(); 
                                     if i==0: radius = radius+2; print 'MDS R=%3.2f'%radius
                                if radius/1000. > rTS[0][taili-1]:
                                    tailj=0 # time
# loop on time
                                    for j in range(size):
                                        if (temps-dt/2. <= time[j]) and (time[j] < temps+dt/2.):
                                            tabtime[tailj]=time[j]
                                            rTS[tailj][taili+taili2]=radius/1000.
                                            tTS[tailj][taili+taili2]=te[j]/TS_conv[choix] #keV
                                            tTSerr[tailj][taili+taili2]=teerr[j]/TS_conv[choix] #keV
                                            tailj=tailj+1
                                    taili2=taili2+1

# remove NaN 
                    for i in range(tailj):  
                        tTS[i]=np.nan_to_num(tTS[i])
                        tTSerr[i]=np.nan_to_num(tTSerr[i])

            except:
                print "Error TS"
            else:
                print "End of TS reading"

            print 'verif', taili+taili2

            if dasel ==1:
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
            else:
                tail2=0

#    print tail
    return (tailj,tail2,taili+taili2,ncore)

    mds.close()


def gprofdatdat(shot,treename,temps,dt):
    import numpy as np

    tail=-1
    tail2=-1

    TS_conv=[1000.,1.e19]

    f=open(savfile+'_size','r')
    stmp=f.readline()
    size=int(stmp.split()[0])
    size2=int(stmp.split()[1])
    nbprobe=int(stmp.split()[2])
    try:ncore=int(stmp.split()[3])
    except: ncore = nbprobe
    f.close()    

    time=np.zeros(size)
    te=np.zeros(size)
    teerr=np.zeros(size)
    for k in range(nbprobe):
        fichier=savfile+str(k+1)+'.npz'
# unzip
        comm='gzip -d '+fichier+'.gz'
        os.system(comm)

        npzfile  = np.load(fichier, mmap_mode='r')
        time=npzfile['arr_0']
        te=npzfile['arr_1']
        teerr=npzfile['arr_2']
        radius=npzfile['arr_3']

# zip again
        comm='gzip '+fichier
        os.system(comm)

        tailj=0 # time
# loop on time
        for j in range(size):
            if (temps-dt/2. <= time[j]) and (time[j] < temps+dt/2.):
                tabtime[tailj]=time[j]
                rTS[tailj][k]=radius/1000.
                tTS[tailj][k]=te[j]/TS_conv[choix] #keV
                tTSerr[tailj][k]=teerr[j]/TS_conv[choix] #keV
                tailj=tailj+1

        for i in range(tailj):  
            tTS[i]=np.nan_to_num(tTS[i])
            tTSerr[i]=np.nan_to_num(tTSerr[i])        

#DA                            
    time=np.zeros(size2)

    Dalpha=np.zeros(size2)
    istart=0
    if dasel==1:
# unzip
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
#    print tabtime
#    print size2,tail2,time[0],time[-1]

                        
    return (tailj,tail2,nbprobe,istart,time,Dalpha,ncore)

def gprofdatmds(shot,treename,temps,dt):
    import numpy as np
    ncore=0
    tail=-1
    tail2=-1
    DA='\\TOR_HA10' # Dalpha
    TS_CORE='\\TS_CORE'
    TS_EDGE='\\TS_EDGE'
    TS_C=['_TE','_NE']
    TS_CE=['_TERRH','_NERRH']
    TS_conv=[1000.,1.e19]

# connection to MDSPlus data
#    with MDS(server="172.17.250.21:8005") as mds:
    with MDS(server="172.17.100.200:8005") as mds:
        try:
            eq=mds.open(shot=shot, tree=treename)
        except: 
            print "Error #2"
        else:
            try:
#            if (1==1):
                print mds.alist

                err=0
# read TS data
                time=eq.get('dim_of(%s)'%'\\TS_CORE2.CORE2_TE').data()
                size=len(time)

                let=-1
                for i in range(size-1):
                    if (let==-1) and (time[i] <= temps) and (temps < time[i+1]):
                        let=i

                if size <2:
                    err=1
                else:
                    taili=0 # radial
# loop on probe (radial)
                    for i in range(nb_ts):
                        dat=TS_CORE+str(i+1)+'.CORE'+str(i+1)+TS_C[choix]
                        try:
                            te=eq.get(dat).data(); ncore = i+1;
                            print 'TS core avail  ',i+1
                        except:
                            print 'TS core unavail',i+1
                        else:

                            dat=TS_CORE+str(i+1)+'.CORE'+str(i+1)+TS_CE[choix]
                            try:
                                teerr=eq.get(dat).data()
                            except:
                                teerr=[x / 10. for x in te]
                                print 'No TS err_bar => 10% ',TS_CORE+str(i+1)

                            if te[let] > 0.:
                                dat=TS_CORE+str(i+1)+'.CORE'+str(i+1)+'_POS'
                                if ((float(shot)>=shot_c) and (float(shot)<shot_d)): radius = CORE_R[i]; print 'Prescribed R=%3.2f'%radius
                                else: radius=eq.get(dat).data(); print 'MDS R=%3.2f'%radius

                                tailj=0 # time
# loop on time
                                for j in range(size):
                                    if (temps-dt/2. <= time[j]) and (time[j] < temps+dt/2.):
                                        tabtime[tailj]=time[j]
                                        rTS[tailj][taili]=radius/1000.
                                        tTS[tailj][taili]=te[j]/TS_conv[choix] #keV
                                        tTSerr[tailj][taili]=teerr[j]/TS_conv[choix] #keV
                                        tailj=tailj+1
                                taili=taili+1

                                fichier=savfile+str(taili)+'.npz'
                                np.savez(fichier,time,te,teerr,radius)

                                comm='gzip '+fichier
                                os.system(comm)

                    taili2=0 # radial
# loop on probe (radial)
                    for i in range(nb_ts):
                        dat=TS_EDGE+str(i+1)+'.EDGE'+str(i+1)+TS_C[choix]
                        try:
                            te=eq.get(dat).data()
                            print 'TS edge avail  ',i+1
                        except:
                            print 'TS edge unavail',i+1
                        else:

                            dat=TS_EDGE+str(i+1)+'.EDGE'+str(i+1)+TS_CE[choix]
                            try:
                                teerr=eq.get(dat).data()
                            except:
                                teerr=[x / 10. for x in te]
                                print 'No TS err_bar => 10% ',TS_EDGE+str(i+1)
                            if te[let] > 0.:
                                dat=TS_EDGE+str(i+1)+'.EDGE'+str(i+1)+'_POS'
                                if ((float(shot)>=shot_c) and (float(shot)<shot_d)): radius = EDGE_R[i]; print 'Prescribed R=%3.2f'%radius
                                else: 
                                     radius=eq.get(dat).data(); 
                                     if i==0: radius = radius + 2; print 'MDS R=%3.2f'%radius
                                if radius/1000. > rTS[0][taili-1]:
                                    tailj=0 # time
# loop on time
                                    for j in range(size):
                                        if (temps-dt/2. <= time[j]) and (time[j] < temps+dt/2.):
                                            tabtime[tailj]=time[j]
                                            rTS[tailj][taili+taili2]=radius/1000.
                                            tTS[tailj][taili+taili2]=te[j]/TS_conv[choix] #keV
                                            tTSerr[tailj][taili+taili2]=teerr[j]/TS_conv[choix] #keV
                                            tailj=tailj+1
                                    taili2=taili2+1

                                    fichier=savfile+str(taili+taili2)+'.npz'
                                    np.savez(fichier,time,te,teerr,radius)
                                    comm='gzip '+fichier
                                    os.system(comm)
# remove NaN 
                    for i in range(tailj):  
                        tTS[i]=np.nan_to_num(tTS[i])
                        tTSerr[i]=np.nan_to_num(tTSerr[i])
                            

            except:
                print "Error TS"
            else:
                print "End of TS reading"

#            print 'verif', taili+taili2

            err=0

            Dalpha=[]
            tstart=0
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

#                print 'tstart=',tstart,size2
#                print 'done'
                if not os.path.isfile(dafile+'.gz'):
#                    f=open(dafile,'w')
#                    for i in range(size2):
#                        f.write('%f %f\n'%(time[i],Dalpha[i]))
#                    f.close()
                    np.savez(dafile,time,Dalpha)
                    comm='gzip '+dafile
                    os.system(comm)


            f=open(savfile+'_size','w')
            f.write('%d %d %d %d'%(size,size2,taili+taili2,ncore))
            f.close()
#    print tail
    mds.close()

    return (tailj,tail2,taili+taili2,tstart,time,Dalpha,ncore)




# python ts5.py shot time dt multiplier percentage DAselect TE/NE
if __name__=='__main__':
   """ """
   import os, sys
   global ncore
   exename=os.path.basename(__file__);
   nargs=len(sys.argv);
   
#   treename = 'HEATING'
   treename='kstar'
   shot=int(sys.argv[1])
   mtime=float(sys.argv[2])
   dt=float(sys.argv[3])
   meanmult=float(sys.argv[4]) # for preselection
   pourcent=float(sys.argv[5]) # for preselection
   dasel=int(sys.argv[6]) # 1 =DA, 0 = no DA
   choix=int(sys.argv[7]) # 0=te, 1=ne
   try: mds_only = sys.argv[8]
   except: mds_only = 'n'

   root = Tkinter.Tk()
   if choix ==0:
       root.wm_title("TE")
   else:
       root.wm_title("NE")
   root.geometry("900x820")
#   root.pack(fill=BOTH, expand=1)

   fig = Figure(figsize=(5, 4), dpi=100)
#   t = np.arange(0, 3, .01)
#fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))

   nb_ts=25
   sizmax=1000
   shape=(sizmax,2*nb_ts)
   tabtime=np.zeros(sizmax)
   rTS=np.zeros(shape)
   tTS=np.zeros(shape)
   tTSerr=np.zeros(shape)
#   tALPHA=np.zeros(sizmax)
#   DALPHA=np.zeros(sizmax)
   if choix ==0:
       savfile='./DATASAVE/'+str(shot)+'/TE'
   else:
       savfile='./DATASAVE/'+str(shot)+'/NE'
   dafile='./DATASAVE/'+str(shot)+'/DA.npz'

# data already saved
   if os.path.isfile(savfile+'1.npz.gz'):
       print 'Data read from DATASAVE'
       (taille,taille2,taille3,istart,tALPHA,DALPHA,ncore)=gprofdatdat(shot,treename,mtime,dt)
   else:
       print 'Data read from KSTAR MDSPlus tree'
#       if mds_only == 'y': exit()
# 1st time...need data from MDSPlus Tree
       (taille,taille2,taille3,istart,tALPHA,DALPHA,ncore)=gprofdatmds(shot,treename,mtime,dt)
       if mds_only == 'y': exit()
# taille = size time TS, taille2= size Dalpha, taille3 = size radial TS
   print taille,taille2,taille3,istart
   if os.path.isfile('./DATASAVE/TS_err'): os.system('rm ./DATASAVE/TS_err')
   if not os.path.isfile(savfile+'1.npz.gz'): 
     f = open('./DATASAVE/TS_err','w')
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
#       print 'debug',len(DALPHA), istart,taille2,moyenne
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

   if choix ==0:
       lymin=Tkinter.Label(root, text="TEmin")
   else:
       lymin=Tkinter.Label(root, text="NEmin")
   lymin.place(y=600, x=770)
   eymint=Tkinter.StringVar()
   eymin = Tkinter.Entry(root,width=10,textvariable=eymint)
   eymin.place(y=600, x=815)
   eymint.set("0.")

   if choix==0:
       lymax=Tkinter.Label(root, text="TEmax")
   else:
       lymax=Tkinter.Label(root, text="NEmax")
   lymax.place(y=620, x=770)
   eymaxt=Tkinter.StringVar()
   eymax = Tkinter.Entry(root,width=10,textvariable=eymaxt)
   eymax.place(y=620, x=815)
   
   value='%.5f'%np.amax(tTS)
   eymaxt.set(value)

   lrmin=Tkinter.Label(root, text="Rmin")
   lrmin.place(y=680, x=770)
   ermint=Tkinter.StringVar()
   ermin = Tkinter.Entry(root,width=10,textvariable=ermint)
   ermin.place(y=680, x=810)
   ermint.set(rTS[0][0])

   lrmax=Tkinter.Label(root, text="Rmax")
   lrmax.place(y=700, x=770)
   ermaxt=Tkinter.StringVar()
   ermax = Tkinter.Entry(root,width=10,textvariable=ermaxt)
   ermax.place(y=700, x=810)
   ermaxt.set(rTS[0][taille3-1])

   buttonRR = Tkinter.Button(master=root, text="Replot", command=_replotR)
   buttonRR.place(y=720, x=790)

   buttonG = Tkinter.Button(master=root, text="Generate", command=_generate2)
   buttonG.place(y=770,x=420)

   buttonQ = Tkinter.Button(master=root, text="Quit", command=_quit)
   buttonQ.place(y=770,x=520)

   Tkinter.mainloop()
# If you put root.destroy() here, it will cause an error if the window is
# closed with the window manager.
