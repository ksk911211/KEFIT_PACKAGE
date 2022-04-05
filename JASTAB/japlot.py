#!/usr/local/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 15:04:32 2018
JA - diagram plot tool (by S.Kim ver.1)
@author: S.Kim, adjusted by C. Lim
"""

import numpy as np
import os, sys
import matplotlib.pyplot as plt
from read_namelist import read_namelist_str as read

currdir = os.getcwd()
input_file = None

#inputs
try:
	input_file = sys.argv[1]
except:
	pass
#input_file='ja_result_6452_mis2'

if (input_file == None):
	print('Command should include [plot_file]')
	exit()

# Stability diagram opt.
use_elite_a=True
use_adj_a = True
use_bilinear=True
nq_cut_off=27.7

crit_dia=0.25
crit_alf=0.03

# Get indiviual plot
plot_spectrum=False
target_i=4
target_j=4

# Cal. check
qdel_plot=False
fill_up=False

# cut off vars
cutoff_n=3
grmax=100

#default variables (do not change)
inf_plot=False
add_n=0
grd_mult=1.0
jphi_mult=1.0

#read ja_opt
if (os.path.isfile('ja_opt')):

    print('Read ja_plot option...')

    f = open('ja_opt','r')

    while True:
        line = f.readline()
        if not line: break
        use_elite_a = read(line,'use_elite_a',use_elite_a,4)
        use_adj_a = read(line,'use_adj_a',use_adj_a,4)
        use_bilinear = read(line,'use_bilinear',use_bilinear,4)
        nq_cut_off = read(line,'nq_cut_off',nq_cut_off,2)
        crit_dia = read(line,'crit_dia',crit_dia,2)
        crit_alf = read(line,'crit_alf',crit_alf,2)
        plot_spectrum = read(line,'plot_spectrum',plot_spectrum,4)
        target_j = read(line,'target_j',target_j,1)
        target_i = read(line,'target_i',target_i,1)
        qdel_plot = read(line,'qdel_plot',qdel_plot,4)
        fill_up = read(line,'fill_up',fill_up,4)
        cutoff_n = read(line,'cutoff_n',cutoff_n,1)
        grmax = read(line,'grmax',grmax,2)                
        inf_plot = read(line,'inf_plot',inf_plot,4)
        add_n = read(line,'add_n',add_n,1)           
        grd_mult = read(line,'grd_mult',grd_mult,2)
        jphi_mult = read(line,'jphi_mult',jphi_mult,2) 
    f.close()

else:
    print('You can add option throught ja_opt')

#read file
#try:
f1=open(input_file,'r')
line1=f1.readline().split()
line2=f1.readline().split()
data=np.transpose(np.array([line1,line2],float))
run_id,gridn,nn=data[:,0]
run_id=int(run_id); gridn=int(gridn); nn=int(nn)
a1,a2,j1=data[:,1]                                  #ref_equilibrium point
bnd_file='stab_bnd_%d' % run_id
lines=f1.readlines()
resultlen=len(lines)
result=np.zeros((resultlen,15))
for i in range(resultlen):
	result[i]=lines[i].split()
f1.close()

#except:
#print('Input file is not avaliable \n')

print('=-----------------------------------=')
print('=---------JA-PLOT Options...--------=')
print('=-----------------------------------=')
print('RUN ID         = %i'%run_id)
print('use_adj_a      = %s'%use_adj_a)
print('use_bilinear   = %s'%use_bilinear)
print('crit_dia       = %4.3f'%crit_dia)
print('crit_alf       = %4.3f'%crit_alf)
print('cutoff_n       = %i'%cutoff_n)
print('nq_cut_off     = %4.2f'%nq_cut_off)
print('qdel_plot      = %s'%qdel_plot)
print('plot_spectrum  = %s'%plot_spectrum)
print('(i,j)          = (%i,%i)'%(target_i,target_j))
print('=-----------------------------------=')

modenn=result[0:nn,8]
for i in range(nn):
    if modenn[i]>=cutoff_n:
        break

if i==nn-1:
    cutn=0
else:
    cutn=i
    
#define variables
alphae=np.zeros((gridn,gridn)); alphat=np.zeros((gridn,gridn)); alphatt=np.zeros((gridn,gridn))
alpha2e=np.zeros((gridn,gridn)); alpha2t=np.zeros((gridn,gridn)); alpha2tt=np.zeros((gridn,gridn))
jphime=np.zeros((gridn,gridn)); jphimt=np.zeros((gridn,gridn)); jphimtt=np.zeros((gridn,gridn))
qind=np.zeros((gridn,gridn)); qedge=np.zeros((gridn,gridn))

infwe = np.copy(alphae)
infwt = np.copy(alphae)
infwtt = np.copy(alphae)

gr=np.zeros((gridn,gridn,nn)); grd=np.copy(gr);
gr_n=np.ones((gridn,gridn,nn))

grm=np.zeros((gridn,gridn))
grm2=np.copy(grm)
grmd=np.zeros((gridn,gridn))
grmd2=np.copy(grm)

readcount=0
jset=0
jset1=jset
for i in range(gridn):
    j=0
    while j<gridn:
        alphae[i,j]=result[readcount,2]
        alpha2e[i,j]=result[readcount,3]
        jphime[i,j]=result[readcount,4]
        
        k=0
        while(jset==jset1):
            jset=result[readcount,1]-1
            nncut_bi=np.ones((nn,1))
            
            grc=result[readcount,9]
            grdc=result[readcount,10]
            if grc<grmax:
                gr[i,j,k]=abs(grc)
                gr_n[i,j,k]=result[readcount,8]
                grd[i,j,k]=abs(grdc)
            readcount=readcount+1
            k=k+1
            if readcount>=len(result):
                break
            jset1=result[readcount,1]-1
        n_dia_temp=np.zeros((1,1,nn))
        if use_bilinear:
            ncut_bi=nq_cut_off/result[readcount-1,12]
            for kk in range(nn):
                if gr_n[i,j,kk]<=ncut_bi:
                    n_dia_temp[0,0,kk]=gr_n[i,j,kk]
                else:
                    n_dia_temp[0,0,kk]=ncut_bi
        else:
            n_dia_temp=gr_n[i,j,:]
        
        X=np.divide(np.multiply(grd[i,j,cutn:nn],gr_n[i,j,cutn:nn]),n_dia_temp[0,0,cutn:nn])
        pa,b=[max(X),np.argmax(X)]
        Y=gr[i,j,cutn:nn]
        pa1,b1=[max(Y),np.argmax(Y)]
        b=b+cutn
        b1=b1+cutn
        grmd[i,j]=grd[i,j,b]*gr_n[i,j,b]/n_dia_temp[0,0,b]*grd_mult
        grmd2[i,j]=gr_n[i,j,b]
        infwe[i,j]=result[readcount-k+b,13]
        if infwe[i,j]>0:
            infwe[i,j]=1/infwe[i,j]
        qedge[i,j]=result[readcount-k+b,12]
        qind[i,j]=result[readcount-k+b,11]
        alphat[i,j]=result[readcount-k+b,5]
        alpha2t[i,j]=result[readcount-k+b,6]
        jphimt[i,j]=result[readcount-k+b,7]
        alphatt[i,j]=result[readcount-k+b1,5]
        alpha2tt[i,j]=result[readcount-k+b1,6]
        jphimtt[i,j]=result[readcount-k+b1,7]
        infwt[i,j]=result[readcount-k+b,14]
        infwtt[i,j]=result[readcount-k+b1,14]
        if infwt[i,j]>0:
            infwt[i,j]=1/infwt[i,j]
        grm[i,j]=gr[i,j,b1]
        grm2[i,j]=gr_n[i,j,b1]
        j=j+1
        if readcount>=len(result):
            break
        jset1=result[readcount,1]-1
        jset=jset1
        if jset!=j:
            j=j+1
    j=0

if use_adj_a:
    alpha = np.copy(alphat); alpha2 = np.copy(alpha2t);   jphim = np.copy(jphimt); infw = np.copy(infwt); 
    alphaa = np.copy(alphatt); alpha2a = np.copy(alpha2tt); jphima = np.copy(jphimtt); infwa = np.copy(infwtt);
else:
    alpha = np.copy(alphae); alpha2 = np.copy(alpha2e);   jphim = np.copy(jphime); infw = np.copy(infwe); 
    alphaa = np.copy(alphae); alpha2a = np.copy(alpha2e); jphima = np.copy(jphime); infwa = np.copy(infwe);

if fill_up:
    for i in range(gridn):
        for j in range(gridn):
            
            if alpha[i,j]==0:
                
                if i==0:
                    alpha[i,j]=2*alpha[i+1,j]-alpha[i+2,j]
                    alpha2[i,j]=2*alpha2[i+1,j]-alpha2[i+2,j]
                    jphim[i,j]=2*jphim[i+1,j]-jphim[i+2,j]
                    grm[i,j]=2*grm[i+1,j]-grm[i+2,j]
                    grm2[i,j]=2*grm2[i+1,j]-grm2[i+2,j]
                    grmd[i,j]=2*grmd[i+1,j]-grmd[i+2,j]
                    grmd2[i,j]=2*grmd2[i+1,j]-grmd2[i+2,j]
                else:
                    alpha[i,j]=2*alpha[i-1,j]-alpha[i-2,j]
                    alpha2[i,j]=2*alpha2[i-1,j]-alpha2[i-2,j]
                    jphim[i,j]=2*jphim[i-1,j]-jphim[i-2,j]
                    grm[i,j]=2*grm[i-1,j]-grm[i-2,j]
                    grm2[i,j]=2*grm2[i-1,j]-grm2[i-2,j]
                    grmd[i,j]=2*grmd[i-1,j]-grmd[i-2,j]
                    grmd2[i,j]=2*grmd2[i-1,j]-grmd2[i-2,j]


grmd2=grmd2+add_n
grm2=grm2+add_n

grmd=grmd

j1=j1/10**6*jphi_mult
jphim=jphim*jphi_mult
jphima=jphima*jphi_mult

if use_elite_a:
    aa = alpha; aa2 = alphaa; aap = a1; jj=jphim; jj2 = jphima;
else:
    aa = alpha2; aa2 = alpha2a; aap = a2; jj=jphim; jj2 = jphima;

plt.figure(1)
contour1=plt.contour(aa,jphim,infw,[1.0])
try:
    inc2=contour1.allsegs[0][0]
except:
    inc2=contour1.allsegs[0]
ill=len(inc2)

contour2=plt.contour(aa,jphim,grmd,[crit_dia])
try:
   c2=contour2.allsegs[0][0]
except:
    c2=contour2.allsegs[0]
ll=len(c2)
plt.close(1)
fig = plt.figure('J-A diagram')
fig.set_size_inches(13,6)
#plt.clf()
plt.subplot(1,2,2)
A=plt.pcolormesh(aa,jphim,grmd,shading='gouraud')
plt.set_cmap('viridis')
#plt.hold('on')
plt.plot(c2[0:ll,0],c2[0:ll,1],'lime',linewidth=3)
plt.title('$\gamma$ / $\omega_{*i}$')
plt.xlabel('Normalized $\\alpha_{max}$')
plt.ylabel('Edge current j$_{\phi, max}$ [MA/$m^2$]')
plt.scatter(aap,j1,s=800,c='r',marker='+')
aaa = np.linspace(min(aa[:,0]),max(aa[:,0]),10)
bbb = np.copy(aaa * j1 / aap)
A.set_clim(0,1)
plt.colorbar(A)

for i in range(gridn):
    for j in range(gridn):
        plt.text(aa[i,j],jj[i,j],str(int(grmd2[i,j])),color='hotpink')
if inf_plot:
    plt.plot(inc2[0:ill,0],inc2[0:ill,1],'b',linewidth=3)
plt.plot(aaa,bbb,c='gold',linestyle='dashed')
plt.figure(10)
contour3=plt.contour(aa,jphim,grm,[crit_alf])
try:
    c3=contour3.allsegs[0][0]
except:
    c3=contour3.allsegs[0]
plt.close(10)

plt.figure('J-A diagram')
#plt.clf()
plt.subplot(1,2,1)
B=plt.pcolormesh(aa2,jj2,grm,shading='gouraud')
plt.set_cmap('viridis')
#plt.hold('on')
plt.plot(c3[0:ll,0],c3[0:ll,1],'lime',linewidth=3)
plt.title('$\gamma$ / $\omega_{A}$')
plt.xlabel('Normalized $\\alpha_{max}$')
plt.ylabel('Edge current j$_{\phi, max}$ [MA/$m^2$]')
B.set_clim(0,0.1)
plt.colorbar(B)
aaa = np.linspace(min(aa2[:,0]),max(aa2[:,0]),10)
bbb = np.copy(aaa * j1 / aap)
for i in range(gridn):
    for j in range(gridn):
        plt.text(aa2[i,j],jj2[i,j],str(int(grm2[i,j])),color='hotpink')
plt.scatter(aap,j1,s=800,c='r',marker='+')
plt.plot(aaa,bbb,c='gold',linestyle='dashed')
plt.tight_layout()
if inf_plot:
    plt.plot(inc2[0:ill,0],inc2[0:ill,1],'b',linewidth=3)
if qdel_plot:
    plt.figure('qdel prof')
    plt.clf()
    plt.pcolor(aa,jphim,qind)

if inf_plot:
    plt.figure('Inf Ballooning')
    plt.clf()
 #   plt.hold('on')
    plt.title('n = $\infinity$ balloning mode')
    plt.xlabel('Normalized $\\alpha_{max}$')
    plt.ylabel('Edge current j$_{\phi, max}$ [MA/$m^2$]')
    C=plt.contourf(aa,jphim,infw)
    plt.plot(inc2[0:ill,0],inc2[0:ill,1],'b',linewidth=3)
    plt.set_cmap('hot')
    C.set_clim(0,2)
    
f=open(bnd_file,'w')
f.write(' Equilibrium point [alp, jphi[MA/m2]]\n')
f.write('%f %f \n' %(aap,j1))
s1=' Diamagnetic criterion = %f \n' %crit_dia
f.write('%s' %s1)
for i in range(len(c2)):
    f.write('%f %f \n' %(c2[i,0],c2[i,1]))
s1=' Alfven criterion = %f \n' %crit_alf
f.write('%s' %s1)
for i in range(len(c3)):
    f.write('%f %f \n' %(c3[i,0],c3[i,1]))
if (inf_plot):
	s1=' n = Inf ball \n'
	f.write('%s' %s1)
	for i in range(len(inc2)):
	    f.write('%f %f \n' %(inc2[i,0],inc2[i,1]))
f.close()
if plot_spectrum:
    plt.figure('Growth rate spectrum')
    plt.clf()
    plt.subplot(1,2,1)
   # plt.hold('on')
    for k in range(nn):
        if gr[target_i,target_j,k]<crit_alf:
            plt.scatter(gr_n[target_i,target_j,k],gr[target_i,target_j,k],c='b')
        else:
            plt.scatter(gr_n[target_i,target_j,k],gr[target_i,target_j,k],c='r')
    sa='$\\alpha$'; sjp='j$_{\phi}$'
    titlen='growth rate %s = %3.1f, %s = %3.1f [MA/$m^2$]' %(sa,aa[target_i,target_j],sjp,jphim[target_i,target_j])
    plt.title('%s' %titlen)
    plt.xlabel('mode n')
    plt.ylabel('$\gamma$ / $\omega_{A}$')
    plt.subplot(1,2,2)
  #  plt.hold('on')
    for k in range(nn):
        if grd[target_i,target_j,k]<crit_dia:
            plt.scatter(gr_n[target_i,target_j,k],grd[target_i,target_j,k],c='b')
        else:
            plt.scatter(gr_n[target_i,target_j,k],grd[target_i,target_j,k],c='r')
    plt.title('%s' %titlen)
    plt.xlabel('mode n')
    plt.ylabel('$\gamma$ / $\omega_{*i}$')

try:
    alpw=np.load('alp_omega_3032')
    alp2=np.zeros((gridn,gridn))
    for i in range(gridn):
        for j in range(gridn):
            alp2[i,j]=1/alpw[gridn*i+j,3]
    plt.figure(10)
    contour4=plt.contour(aa,jphim,alp2,[1.0])
    try:
        c4=contour4.allsegs[0][0]
    except:
        c4=contour4.allsegs[0]
    ll=len(c3)
    plt.close(10)
    plt.figure(6)
    D=plt.pcolormesh(aa,jphim,alp2,shading='gouraud')
    #plt.hold('on')
    plt.set_cmap('hot')
    D.set_clim(0,2)
    plt.plot(c3[0:ll,0],c3[0:ll,1],'b',linewidth=3)
    plt.close(6)
except:
    pass

plt.show(block=False)
input('Press Enter to close the figure')
