#!/usr/local/anaconda3/bin/python3
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from read_namelist import read_namelist_str as read

class japlot:

    def initialise_variables(self,filename):

        self.currdir = os.getcwd()
        self.input_file = None
        if filename == None:
                try:
                      self.input_file = sys.argv[1]
                except:
                      pass
        else:
                self.input_file = filename

        if (self.input_file == None):
            print('Command should include [plot_file]')
            exit()

        # Stability diagram opt.
        self.use_elite_a=True
        self.use_adj_a = True
        self.use_j2    = False
        self.use_bilinear=True
        self.nq_cut_off=27.7

        self.crit_dia=0.25
        self.crit_alf=0.03

        # Get indiviual plot
        self.plot_spectrum=True
        self.target_i=1
        self.target_j=1

        # Cal. check
        self.qdel_plot=False
        self.fill_up=False

        # cut off vars
        self.cutoff_n=3
        self.grmax=100

        #default variables (do not change)
        self.inf_plot=False
        self.add_n=0
        self.grd_mult=1.0
        self.jphi_mult=1.0

        return

    def read_namelist(self):

        if (os.path.isfile('ja_opt')):
            print('Read ja_plot option...')

            f = open('ja_opt','r')

            while True:
                line = f.readline()
                if not line: break
                self.use_elite_a = read(line,'use_elite_a',self.use_elite_a,4)
                self.use_adj_a = read(line,'use_adj_a',self.use_adj_a,4)
                self.use_j2 = read(line,'use_j2',self.use_j2,4)
                self.use_bilinear = read(line,'use_bilinear',self.use_bilinear,4)
                self.nq_cut_off = read(line,'nq_cut_off',self.nq_cut_off,2)
                self.crit_dia = read(line,'crit_dia',self.crit_dia,2)
                self.crit_alf = read(line,'crit_alf',self.crit_alf,2)
                self.plot_spectrum = read(line,'plot_spectrum',self.plot_spectrum,4)
                self.target_j = read(line,'target_j',self.target_j,1)
                self.target_i = read(line,'target_i',self.target_i,1)
                self.qdel_plot = read(line,'qdel_plot',self.qdel_plot,4)
                self.fill_up = read(line,'fill_up',self.fill_up,4)
                self.cutoff_n = read(line,'cutoff_n',self.cutoff_n,1)
                self.grmax = read(line,'grmax',self.grmax,2)                
                self.inf_plot = read(line,'inf_plot',self.inf_plot,4)
                self.add_n = read(line,'add_n',self.add_n,1)           
                self.grd_mult = read(line,'grd_mult',self.grd_mult,2)
                self.jphi_mult = read(line,'jphi_mult',self.jphi_mult,2) 
            f.close()

        else:
            print('You can add option throught ja_opt')

        return

    def read_result(self,printr=False):

        f1=open(self.input_file,'r')
        line1=f1.readline().split()
        line2=f1.readline().split()
        #data=np.transpose(np.array([line1,line2],float))
        if len(line1) == 3:
            run_id,gridn1,nn=line1
            gridn2 = gridn1
        else:
            run_id,gridn1,gridn2,nn=line1
        self.run_id=run_id; self.gridn1=int(gridn1); self.gridn2=int(gridn2); self.nn=int(nn)
        a1,a2,j1,je = line2
        self.a1 = float(a1); self.a2 = float(a2); self.j1 = float(j1); self.je = float(je); self.j0 = float(j1);
        self.bnd_file='stab_bnd_%s' % run_id
        lines=f1.readlines()
        resultlen=len(lines)
        self.result=np.zeros((resultlen,16))
        for i in range(resultlen):
            self.result[i]=lines[i].split()
        f1.close()
  
        if printr:
            print('=-----------------------------------=')
            print('=---------JA-PLOT Options...--------=')
            print('=-----------------------------------=')
            print('RUN ID         = %i'%run_id)
            print('use_adj_a      = %s'%use_adj_a)
            print('use_j2         = %s'%use_j2)
            print('use_bilinear   = %s'%use_bilinear)
            print('crit_dia       = %4.3f'%crit_dia)
            print('crit_alf       = %4.3f'%crit_alf)
            print('cutoff_n       = %i'%cutoff_n)
            print('nq_cut_off     = %4.2f'%nq_cut_off)
            print('qdel_plot      = %s'%qdel_plot)
            print('plot_spectrum  = %s'%plot_spectrum)
            print('(i,j)          = (%i,%i)'%(target_i,target_j))
            print('=-----------------------------------=')

        return               

    def post_processing(self):

        self.modenn=self.result[0:self.nn,8]
        for i in range(self.nn):
            if self.modenn[i]>=self.cutoff_n:
                break

        if i==self.nn-1:
            self.cutn=0
        else:
            self.cutn=i
    
    #define variables
        self.alphae=np.zeros((self.gridn1,self.gridn2)); self.alphat=np.zeros((self.gridn1,self.gridn2)); self.alphatt=np.zeros((self.gridn1,self.gridn2))
        self.alpha2e=np.zeros((self.gridn1,self.gridn2)); self.alpha2t=np.zeros((self.gridn1,self.gridn2)); self.alpha2tt=np.zeros((self.gridn1,self.gridn2))
        self.jphime=np.zeros((self.gridn1,self.gridn2)); self.jphimt=np.zeros((self.gridn1,self.gridn2)); self.jphimtt=np.zeros((self.gridn1,self.gridn2))
        self.qind=np.zeros((self.gridn1,self.gridn2)); self.qedge=np.zeros((self.gridn1,self.gridn2)); self.jphimed=np.zeros((self.gridn1,self.gridn2));

        self.infwe = np.copy(self.alphae)
        self.infwt = np.copy(self.alphae)
        self.infwtt = np.copy(self.alphae)

        self.gr=np.zeros((self.gridn1,self.gridn2,self.nn)); self.grd=np.copy(self.gr);
        self.gr_n=np.ones((self.gridn1,self.gridn2,self.nn))

        self.grm=np.zeros((self.gridn1,self.gridn2))
        self.grm2=np.copy(self.grm)
        self.grmd=np.zeros((self.gridn1,self.gridn2))
        self.grmd2=np.copy(self.grm)

        self.n_dia=np.ones((self.gridn1,self.gridn2,self.nn)) 
        self.readcount1 = np.zeros((self.gridn1,self.gridn2))
        self.readcount2 = np.zeros((self.gridn1,self.gridn2))

        jset = -1
        for ii in range(len(self.result[:,1])):
            result = self.result[ii,:]
            i = int(result[0])-1
            j = int(result[1])-1
            self.alphae[i,j]=result[2]
            self.alpha2e[i,j]=result[3]
            self.jphime[i,j]=result[4]

            if(jset==j):
                k = k + 1
            else:
                jset = j
                k = 0
                self.readcount1[i,j] = int(ii)

            if not (ii == len(self.result[:,1])-1):
                if not (self.result[ii,1] == self.result[ii+1,1]):
                    self.readcount2[i,j] = k+1
            else:
                self.readcount2[i,j] = k+1
                
            grc=result[9]
            grdc=result[10]
            if grc<self.grmax:
                self.gr[i,j,k]=abs(grc)
                self.gr_n[i,j,k]=result[8]
                self.grd[i,j,k]=abs(grdc)

            ncut_bi=self.nq_cut_off/result[12]
            if self.use_bilinear:
                if (self.gr_n[i,j,k] <=ncut_bi):
                    self.n_dia[i,j,k] = self.gr_n[i,j,k]
                else:
                    self.n_dia[i,j,k] = ncut_bi
            else:
                self.n_dia[i,j,k] = self.gr_n[i,j,k]

        for i in range(self.gridn1):
            for j in range(self.gridn2):

                if (self.alphae[i,j]>0.):
                    X=np.divide(np.multiply(self.grd[i,j,self.cutn:self.nn],self.gr_n[i,j,self.cutn:self.nn]),self.n_dia[i,j,self.cutn:self.nn])
                    pa,b=[max(X),np.argmax(X)]
                    Y=self.gr[i,j,self.cutn:self.nn]
                    pa1,b1=[max(Y),np.argmax(Y)]
                    b=b+self.cutn
                    b1=b1+self.cutn
                    self.grmd[i,j]=self.grd[i,j,b]*self.gr_n[i,j,b]/self.n_dia[i,j,b]*self.grd_mult
                    self.grmd2[i,j]=self.gr_n[i,j,b]
                    nline = int(self.readcount1[i,j])
                    if self.readcount2[i,j] == self.nn:
                        nline1 = nline + b
                        nline2 = nline + b1
                    else:
                        for k in range(int(self.readcount2[i,j])):
                            if self.gr_n[i,j,b] == self.result[nline+k,8]:
                                nline1 = nline + k
                            if self.gr_n[i,j,b1] == self.result[nline+k,8]:
                                nline2 = nline + k                                

                    self.infwe[i,j]=self.result[nline1,13]

                    if self.infwe[i,j]>0:
                        self.infwe[i,j] = 1/self.infwe[i,j]

                    self.qedge[i,j]   = self.result[nline1,12]
                    self.qind[i,j]    = self.result[nline1,11]
                    self.alphat[i,j]  = self.result[nline1,5]
                    self.alpha2t[i,j] = self.result[nline1,6]
                    self.jphimt[i,j]  = self.result[nline1,7]
                    self.alphatt[i,j] = self.result[nline2,5]
                    self.alpha2tt[i,j]= self.result[nline2,6]
                    self.jphimtt[i,j] = self.result[nline2,7]
                    self.infwt[i,j]   = self.result[nline1,14]
                    self.infwtt[i,j]  = self.result[nline2,14]
                    self.jphimed[i,j] = self.result[nline2,15]		

                    if self.infwt[i,j]>0:
                        self.infwt[i,j]=1/self.infwt[i,j]
                    self.grm[i,j]  =self.gr[i,j,b1]
                    self.grm2[i,j] =self.gr_n[i,j,b1]

        if self.use_adj_a:
            self.alpha = np.copy(self.alphat); self.alpha2 = np.copy(self.alpha2t);   self.jphim = np.copy(self.jphimt); self.infw = np.copy(self.infwt); 
            self.alphaa = np.copy(self.alphatt); self.alpha2a = np.copy(self.alpha2tt); self.jphima = np.copy(self.jphimtt); self.infwa = np.copy(self.infwtt);
        else:
            self.alpha = np.copy(self.alphae); self.alpha2 = np.copy(self.alpha2e);   self.jphim = np.copy(self.jphime); self.infw = np.copy(self.infwe); 
            self.alphaa = np.copy(self.alphae); self.alpha2a = np.copy(self.alpha2e); self.jphima = np.copy(self.jphime); self.infwa = np.copy(self.infwe);

        if self.use_j2:
           self.j1    = 0.5*(self.j0+self.je)
           if self.use_adj_a:
            self.jphim = 0.5*(self.jphimt+self.jphimed)
            self.jphima= 0.5*(self.jphimtt+self.jphimed)
           else:
            self.jphim = 0.5*(self.jphime+self.jphimed)
            self.jphima= 0.5*(self.jphime+self.jphimed)
        else: self.j1 = self.j0
        
        for i in range(self.gridn1):
            for j in range(self.gridn2):
                if self.alpha[i,j]==0:
                    if self.fill_up:

                        if (i==0):
                            self.alpha[i,j]=2*self.alpha[i+1,j]-self.alpha[i+2,j]
                            self.alpha2[i,j]=2*self.alpha2[i+1,j]-self.alpha2[i+2,j]
                            self.jphim[i,j]=2*self.jphim[i+1,j]-self.jphim[i+2,j]
                            self.grm[i,j]=2*self.grm[i+1,j]-self.grm[i+2,j]
                            self.grm2[i,j]=2*self.grm2[i+1,j]-self.grm2[i+2,j]
                            self.grmd[i,j]=2*self.grmd[i+1,j]-self.grmd[i+2,j]
                            self.grmd2[i,j]=2*self.grmd2[i+1,j]-self.grmd2[i+2,j]
                        else:
                            self.alpha[i,j]=2*self.alpha[i-1,j]-self.alpha[i-2,j]
                            self.alpha2[i,j]=2*self.alpha2[i-1,j]-self.alpha2[i-2,j] 
                            self.jphim[i,j]=2*self.jphim[i-1,j]-self.jphim[i-2,j]
                            self.grm[i,j]=2*self.grm[i-1,j]-self.grm[i-2,j]
                            self.grm2[i,j]=2*self.grm2[i-1,j]-self.grm2[i-2,j]
                            self.grmd[i,j]=2*self.grmd[i-1,j]-self.grmd[i-2,j]
                            self.grmd2[i,j]=2*self.grmd2[i-1,j]-self.grmd2[i-2,j]

                    else:
                        if i < self.gridn1*0.5:
                            i_ind = self.gridn1
                            i_del = 1
                        else:
                            i_ind = 0
                            i_del = -1                            
                        if j < self.gridn2*0.5:
                            j_ind = self.gridn2
                            j_del = 1
                        else:
                            j_ind = 0
                            j_del = -1  
                        count = False
                        for ii in range(i,i_ind,i_del):
                            for jj in range(j,j_ind,j_del):
                                if self.alpha[ii,jj] > 0.: 
                                    print('>> Error found at [nw,nh]=[%i,%i] --> Ingore...'%(ii,jj))
                                    count = True
                                    break
                            if count:   break

                        self.alpha[i,j]=self.alpha[ii,jj]
                        self.alpha2[i,j]=self.alpha2[ii,jj]
                        self.alphaa[i,j]=self.alphaa[ii,jj]
                        self.alpha2a[i,j]=self.alpha2a[ii,jj]
                        self.jphim[i,j]=self.jphim[ii,jj]
                        self.jphima[i,j]=self.jphima[ii,jj]
                        self.grm[i,j]=self.grm[ii,jj]
                        self.grm2[i,j]=self.grm2[ii,jj]
                        self.grmd[i,j]=self.grmd[ii,jj]
                        self.grmd2[i,j]=self.grmd2[ii,jj]

        self.grmd2=self.grmd2+self.add_n
        self.grm2=self.grm2+self.add_n

        self.j2=self.j1/10**6*self.jphi_mult
        self.jphim=self.jphim*self.jphi_mult
        self.jphima=self.jphima*self.jphi_mult

        if self.use_elite_a:
            self.aa = self.alpha; self.aa2 = self.alphaa; self.aap = self.a1; self.jj=self.jphim; self.jj2 = self.jphima;
        else:
            self.aa = self.alpha2; self.aa2 = self.alpha2a; self.aap = self.a2; self.jj=self.jphim; self.jj2 = self.jphima;

        return

    def get_close_point(self):

        xind = 0;   yind = 0;
        xerr = np.zeros(2)
        yerr = np.zeros(2)
        val0 = (self.aa[0,0]-self.aap)**2 + (self.jphim[0,0]-self.j2)**2
        for i in range(self.gridn1):
            for j in range(self.gridn2):
                val = (self.aa[i,j]-self.aap)**2 + (self.jphim[i,j]-self.j2)**2
                if (val <= val0):
                    xind = i
                    yind = j
                    val0 = val

        xdel = self.aap - self.aa[xind,yind]
        ydel = self.j2 - self.jphim[xind,yind] 


        if xdel >= 0. and ydel >= 0.:
            xerr[0] = xdel
            xerr[1] = self.aa[xind+1,yind] - self.aap 
            yerr[0] = ydel
            yerr[1] = self.jphim[xind,yind+1] - self.j2
        elif xdel >= 0. and ydel <= 0.:
            xerr[0] = xdel
            xerr[1] = self.aa[xind+1,yind] - self.aap 
            yerr[0] = self.j2 -self.jphim[xind,yind-1]
            yerr[1] = abs(ydel)
        elif xdel <= 0. and ydel >= 0.:
            xerr[0] = self.aap - self.aa[xind-1,yind]
            xerr[1] = abs(xdel)
            yerr[0] = ydel
            yerr[1] = self.jphim[xind,yind+1] - self.j2
        elif xdel <= 0. and ydel <= 0.:       
            xerr[0] = self.aap - self.aa[xind-1,yind]
            xerr[1] = abs(xdel)
            yerr[0] = self.j2 -self.jphim[xind,yind-1]
            yerr[1] = abs(ydel)       

        if abs(xdel) <0.1 and abs(ydel) <0.1:
            if xind == len(self.aa[:,0])-1: xind = xind -1
            if yind == len(self.aa[0,:])-1: yind = yind -1
            xerr[0] = (self.aa[xind+1,yind] - self.aa[xind-1,yind])*0.5
            xerr[1] = xerr[0]
            yerr[0] = (self.jphim[xind,yind+1] - self.jphim[xind,yind-1])*0.5
            yerr[1] = yerr[0]

        xerr = 0.5*xerr*1.1
        yerr = 0.5*yerr*1.1


        return (xerr,yerr)


    def draw_plot(self,ftype=1,fig_ex=None):

        xerr,yerr = self.get_close_point()

        if fig_ex==None:
            if (ftype == 1 or ftype ==2):
                fig = plt.figure('J-A diagram')
                ax1 = fig.add_subplot(1,1,1)
            elif (ftype == 3):
                fig = plt.figure('qdel prof')
                ax1 = fig.add_subplot(1,1,1)
            elif (ftype == 4):
                fig = plt.figure('Inf Ballooning')
                ax1 = fig.add_subplot(1,1,1)

            elif (ftype == 5):
                fig = plt.figure('Growth rate spectrum')
                ax1 = fig.add_subplot(1,1,1)

        else:
            fig = fig_ex
            len2 = len(fig.axes)
            if len2 == 1:
                [ax1] = fig.axes
            else:
                [ax1, ax2] = fig.axes
                ax2.cla()

            ax1.cla()

        if (ftype==1 or ftype==2 or ftype==3):
            contour1=ax1.contour(self.aa,self.jphim,self.infw,[1.0])
            try:
                inc2=contour1.allsegs[0][0]
            except:
                inc2=contour1.allsegs[0]
            ill=len(inc2)

            if ftype == 2:   contour2=ax1.contour(self.aa,self.jphim,self.grmd,[self.crit_dia])
            elif ftype==1:   contour2=ax1.contour(self.aa,self.jphim,self.grm,[self.crit_alf])
            ll = len(contour2.allsegs[0])
            c2 = contour2.allsegs[0]
            

        ax1.cla()

        if (ftype ==2):

            A=ax1.pcolormesh(self.aa,self.jphim,self.grmd,shading='gouraud')
#            A=ax1.contourf(self.aa,self.jphim,self.grmd)
            #print(fig.axes)
            plt.set_cmap('viridis')
            #print(fig.axes)
            for i in range(ll):
                 ax1.plot(c2[i][:,0],c2[i][:,1],'lime',linewidth=3)
            ax1.set_title('$\gamma$ / $\omega_{*i}$')
            ax1.set_xlabel('Normalized $\\alpha_{max}$')
            if self.use_j2:
              ax1.set_ylabel('Edge current (j$_{\phi, max}$+j$_{\phi, sep}$)/2 [$MA/m^2$]')
            else:
              ax1.set_ylabel('Edge current j$_{\phi, max}$ [$MA/m^2$]')
            ax1.errorbar(self.aap,self.j2,xerr=np.mean(xerr),yerr=np.mean(yerr),fmt='o',markersize='5',capsize=1.5,capthick=4,ecolor='r',c='r')
           # ax1.scatter(self.aap,self.j1,s=800,c='r',marker='+')
            ax1.scatter(self.aa[self.target_i,self.target_j],self.jphim[self.target_i,self.target_j],c='orange',marker='x',s=50)
            aaa = np.linspace(min(self.aa[:,0]),max(self.aa[:,0]),10)
            bbb = np.copy(aaa * self.j2 / self.aap)
            ax1.plot(aaa,bbb,c='gold',linestyle='dashed')
            A.set_clim(0,1)
            if len2 == 1:
                 fig.colorbar(A,orientation='vertical')
            else:
                 fig.colorbar(A,cax=ax2,orientation='vertical')
            fig.tight_layout()
            f = open('ja_raw.dat','w')
            f.write('%9.6f\t%9.6f\n'%(self.aap,self.j2))
            for i in range(self.gridn1):
                for j in range(self.gridn2):
                    ax1.text(self.aa[i,j],self.jj[i,j],str(int(self.grmd2[i,j])),color='hotpink')
                    f.write('%2i\t%2i\t%13.7e\t%13.7e\t%13.7e\t%i\n'%(i,j,self.aa[i,j],self.jphim[i,j],self.grmd[i,j],self.grmd2[i,j]))
            f.close()
            if self.inf_plot:
                ax1.plot(inc2[0:ill,0],inc2[0:ill,1],'b',linewidth=3)

            
        elif (ftype ==1):

            B=ax1.pcolormesh(self.aa2,self.jj2,self.grm,shading='gouraud')
            #plt.set_cmap('viridis')
            for i in range(ll):
               ax1.plot(c2[i][:,0],c2[i][:,1],'lime',linewidth=3)
            ax1.set_title('$\gamma$ / $\omega_{A}$')
            ax1.set_xlabel('Normalized $\\alpha_{max}$')
            if not self.use_j2: ax1.set_ylabel('Edge current j$_{\phi, max}$ [$MA/m^2$]')
            else: ax1.set_ylabel('Edge current (j$_{\phi, max}$+j$_{\phi, sep}$)/2 [$MA/m^2$]')
            B.set_clim(0,0.1)
            if len2 == 1:
                 fig.colorbar(B,orientation='vertical')
            else:
                 fig.colorbar(B,cax=ax2,orientation='vertical')
            fig.tight_layout()
            aaa = np.linspace(min(self.aa2[:,0]),max(self.aa2[:,0]),10)
            bbb = np.copy(aaa * self.j2 / self.aap)
            f = open('ja_raw.dat','w')
            f.write('%9.6f\t%9.6f\n'%(self.aap,self.j2))
            for i in range(self.gridn1):
                for j in range(self.gridn2):
                    ax1.text(self.aa2[i,j],self.jj2[i,j],str(int(self.grm2[i,j])),color='hotpink')
                    f.write('%2i\t%2i\t%9.6f\t%9.6f\t%9.6f\t%i\n'%(i,j,self.aa2[i,j],self.jj2[i,j],self.grm[i,j],self.grm2[i,j]))
            #ax1.scatter(self.aap,self.j1,s=800,c='r',marker='+')
            f.close()
            ax1.errorbar(self.aap,self.j2,xerr=np.mean(xerr),yerr=np.mean(yerr),fmt='o',markersize='5',capsize=1.5,capthick=4,ecolor='r',c='r')
            ax1.scatter(self.aa[self.target_i,self.target_j],self.jphim[self.target_i,self.target_j],c='orange',marker='x',s=50)
            ax1.plot(aaa,bbb,c='gold',linestyle='dashed')
#            plt.tight_layout()
            if self.inf_plot:
                ax1.plot(inc2[0:ill,0],inc2[0:ill,1],'b',linewidth=3)

        elif (ftype ==3):            
            ax1.pcolor(self.aa,self.jphim,self.qind)

        elif(ftype == 4):
            ax1.set_title('n = $\infinity$ balloning mode')
            ax1.xlabel('Normalized $\\alpha_{max}$')
            if not self.use_j2: ax1.ylabel('Edge current j$_{\phi, max}$ [$MA/m^2$]')
            else: ax1.set_ylabel('Edge current (j$_{\phi, max}$+j$_{\phi, sep}$)/2 [$MA/m^2$]')
            C=ax1.contourf(self.aa,self.jphim,self.infw)
            ax1.plot(inc2[0:ill,0],inc2[0:ill,1],'b',linewidth=3)
            ax1.set_cmap('hot')
            ax1.set_clim(0,2)

        elif (ftype == 5):
            ax1.plot(self.gr_n[self.target_i,self.target_j,:],self.gr[self.target_i,self.target_j,:],'--',c='gold') 
            for k in range(self.nn):
                if self.gr[self.target_i,self.target_j,k]<self.crit_alf:
                    ax1.scatter(self.gr_n[self.target_i,self.target_j,k],self.gr[self.target_i,self.target_j,k],c='b')
                else:
                    ax1.scatter(self.gr_n[self.target_i,self.target_j,k],self.gr[self.target_i,self.target_j,k],c='r')
            

            sa='$\\alpha$'; sjp='j$_{\phi}$'
            titlen='growth rate %s = %3.1f, %s = %3.1f [$MA/m^2$]' %(sa,self.aa[self.target_i,self.target_j],sjp,self.jphim[self.target_i,self.target_j])
            ax1.set_title('%s' %titlen)
            ax1.set_xlabel('mode n')
            ax1.set_ylabel('$\gamma$ / $\omega_{A}$')
            ax1.set_xlim([min(self.gr_n[self.target_i,self.target_j,:])-1,max(self.gr_n[self.target_i,self.target_j,:])+1])
            minx = min(self.gr_n[self.target_i,self.target_j,:])-1
            maxx = max(self.gr_n[self.target_i,self.target_j,:])+1            
            delx = maxx-minx
            delx = int(delx/6)

            ax1.set_xticks(range(int(minx),int(maxx+1),delx))            
        elif (ftype == 6):
            weight = self.gr_n[self.target_i,self.target_j,:]/self.n_dia[self.target_i,self.target_j,:]*self.grd_mult       
            ax1.plot(self.gr_n[self.target_i,self.target_j,:],self.grd[self.target_i,self.target_j,:]*weight,'--',c='gold') 
            for k in range(self.nn):
                weight = self.gr_n[self.target_i,self.target_j,k]/self.n_dia[self.target_i,self.target_j,k]*self.grd_mult
                if self.grd[self.target_i,self.target_j,k]*weight<self.crit_dia:
                    ax1.scatter(self.gr_n[self.target_i,self.target_j,k],self.grd[self.target_i,self.target_j,k]*weight,c='b')
                else:
                    ax1.scatter(self.gr_n[self.target_i,self.target_j,k],self.grd[self.target_i,self.target_j,k]*weight,c='r')

            sa='$\\alpha$'; sjp='j$_{\phi}$'
            titlen='growth rate %s = %3.1f, %s = %3.1f [$MA/m^2$]' %(sa,self.aa[self.target_i,self.target_j],sjp,self.jphim[self.target_i,self.target_j])        
            ax1.set_title('%s' %titlen)
            ax1.set_xlabel('mode n')
            ax1.set_ylabel('$\gamma$ / $\omega_{*i}$')
            ax1.set_xlim([min(self.gr_n[self.target_i,self.target_j,:])-1,max(self.gr_n[self.target_i,self.target_j,:])+1])
            minx = min(self.gr_n[self.target_i,self.target_j,:])-1
            maxx = max(self.gr_n[self.target_i,self.target_j,:])+1            
            delx = maxx-minx
            delx = int(delx/6)

            ax1.set_xticks(range(int(minx),int(maxx+1),delx))

        if (ftype==1 or ftype==2):
    
            f=open(self.bnd_file,'w')
            f.write(' Equilibrium point [alp, jphi[MA/m2]]\n')
            f.write('%f %f \n' %(self.aap,self.j2))

            if (ftype==2):
                s1=' Diamagnetic criterion = %f \n' %self.crit_dia
            else:
                s1=' Alfven criterion = %f \n' %self.crit_alf

            f.write('%s' %s1)
            for j in range(ll):
               for i in range(len(c2[j])):
                   f.write('%f %f \n' %(c2[j][i,0],c2[j][i,1]))

            if (self.inf_plot):
                s1=' n = Inf ball \n'
                f.write('%s' %s1)
                for j in range(ll):
                   for i in range(len(inc2)):
                       f.write('%f %f \n' %(inc2[i,0],inc2[i,1]))
                f.close()
        if(fig==None): 
                plt.show(block=False)
        else:
                plt.draw()

    def __init__(self,filename=None):

        self.initialise_variables(filename)
        self.read_namelist()
        self.read_result()

        return

if __name__ == "__main__":

    import gjaplot

    gjap = gjaplot.japlot()  
    gjap.post_processing()
    fig1,ax0 = plt.subplots(1,1)
    gjap.draw_plot(1,fig1)
    fig2,ax0 = plt.subplots(1,1)
    gjap.draw_plot(2,fig2)
    gjap.plot_spectrum = False
    plt.show(block=False)
    if gjap.plot_spectrum:
        fig3 = plt.figure(3)
        ax3=fig3.add_subplot(1,1,1)
        gjap.draw_plot(5,fig3)
        fig4 = plt.figure(4)
        ax4=fig4.add_subplot(1,1,1)
        gjap.draw_plot(6,fig4)

    input()
