# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 13:01:15 2014

@author: root
"""

import os
import csv
import numpy as np
from scipy.integrate import quad
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import leastsq
import scipy.constants as c
import random
import time

#from PlotNManyElectronsT2_24Nov2014 import *

import cProfile, pstats, StringIO
pr = cProfile.Profile()
pr.enable()

#N=1
#
#Q1 = 1e9#np.float('inf')
#Qe = np.ones(N)*Q1
#Q0=1000
#beta = 1.0
#sigma = 1e-5
#betas = beta*np.random.normal(0,sigma,N)#*np.ones(N)
#
#gamma = 1.0
#k = 0
#
#wz = 300e6 # collission / axial bounce frequency
#wo = 30e9   # cyclotron frequency
#wr = wo/gamma # reference frequency
#
#Tcol = gamma*wo/wz
#print Tcol
#
#dt_ave = np.pi/1000.0


#Original def sin1_sin2(t,f1,f2,P1,P2):
#    t1 = np.arange(t-np.pi,t+np.pi,dt_ave)
#    sb = np.sin(f1*t1 + P1)
#    s = np.sin(f2*t1 +P2)
#    
#    mean = np.sum(sb*s)*dt_ave/np.pi/2.0
#    
#    return mean
    
def sin1_sin2(t,f1,f2,P1,P2):
    #t1 = np.arange(t-np.pi,t+np.pi,dt_ave)
    f = lambda x: np.sin(f1*x+P1)*np.sin(f2*x+P2)
    #sb = np.sin(f1*t1 + P1)
    #s = np.sin(f2*t1 +P2)
    
    mean = quad(f,t-np.pi*gamma,t+np.pi*gamma)[0]/2/np.pi/gamma
    #print mean
    
    
    #mean = np.mean(np.sin(f2*t1 +P2)*np.sin(f1*t1 + P1))
    
    return mean
    
def sin1_cos2(t,f1,f2,P1,P2):
    #t1 = np.arange(t-np.pi,t+np.pi,dt_ave)
    #sb = np.sin(f1*t1 + P1)
    #s = np.cos(f2*t1 +P2)
    
    #mean = np.sum(sb*s)*dt_ave/np.pi/2.0
    f = lambda x: np.sin(f1*x+P1)*np.cos(f2*x+P2)
    #sb = np.sin(f1*t1 + P1)
    #s = np.sin(f2*t1 +P2)
    
    mean = quad(f,t-np.pi*gamma,t+np.pi*gamma)[0]/2/np.pi/gamma
    return mean

def cos1_sin2(t,f1,f2,P1,P2):
#    t1 = np.arange(t-np.pi,t+np.pi,dt_ave)
#    sb = np.cos(f1*t1 + P1)
#    s = np.sin(f2*t1 +P2)
#    
#    mean = np.sum(sb*s)*dt_ave/np.pi/2.0
    f = lambda x: np.cos(f1*x+P1)*np.sin(f2*x+P2)
    #sb = np.sin(f1*t1 + P1)
    #s = np.sin(f2*t1 +P2)
    
    mean = quad(f,t-np.pi*gamma,t+np.pi*gamma)[0]/2/np.pi/gamma
    
    return mean
    
def cos1_cos2(t,f1,f2,P1,P2):
#    t1 = np.arange(t-np.pi,t+np.pi,dt_ave)
#    sb = np.cos(f1*t1 + P1)
#    s = np.cos(f2*t1 +P2)
#    
#    mean = np.sum(sb*s)*dt_ave/np.pi/2.0
    f = lambda x: np.cos(f1*x+P1)*np.cos(f2*x+P2)
    #sb = np.sin(f1*t1 + P1)
    #s = np.sin(f2*t1 +P2)
    
    mean = quad(f,t-np.pi*gamma,t+np.pi*gamma)[0]/2/np.pi/gamma
    
    return mean

def sinsin_matrix(t,phi,theta): # this matrix symmetric!?
    M = np.zeros((N+1,N+1))
    for k in range(0,len(theta)):
        M[0][k+1] = sin1_sin2(t,gamma,betas[k],phi,theta[k])
    

    for i in range(0,len(theta)):
        for j in range(0,len(theta)):
            entry = sin1_sin2(t,betas[i],betas[j],theta[i],theta[j])
            M[i+1][j+1] = entry
    return M
        
def coscos_matrix(t,phi,theta): # this matrix symmetric!?
    M = np.zeros((N+1,N+1))
    M[0][0] = 0.5
    for k in range(0,len(theta)):
        M[0][k+1] = cos1_cos2(t,gamma,betas[k],phi,theta[k])

    for i in range(0,len(theta)):
        for j in range(0,len(theta)):
            entry = cos1_cos2(t,betas[i],betas[j],theta[i],theta[j])
            M[i+1][j+1] = entry
    return M
    
def sin1cos2_matrix(t,phi,theta): # this matrix symmetric!?
    M = np.zeros((N+1,N+1))
    M[0][0] = 0.0
    for k in range(0,len(theta)):
        M[0][k+1] = sin1_cos2(t,gamma,betas[k],phi,theta[k])
        
    for k in range(0,len(theta)):
        M[k+1][0] = sin1_cos2(t,betas[k],gamma,theta[k],phi)

    for i in range(0,len(theta)):
        for j in range(0,len(theta)):
            if i==j:
                entry = 0
                M[i+1][j+1] = entry
                continue
            entry = sin1_cos2(t,betas[i],betas[j],theta[i],theta[j])
            M[i+1][j+1] = entry
    return M
    
################


#################    

#def ridot(t,i,rho,phi,r,theta,Mss,Msc,Mcc):
#    
#    cri = -betas[i]*(1.0 + k**2/(1.0-N*k**2) ) /Qe[i]/2.0
#    
#    crho = k*gamma**2/(1-N*k**2)/betas[i]*( Mss[0][i+1]/Q0 - Msc[0][i+1])
#
########## SOMasdfa
#asdfasdf
#asdfasdfasd
#something wrong with index for Mss etc for N=1 anyways
##
#    cij =  betas**2*( - Mss[i+1][1:]/Qe + Msc[i+1][1:] )  # Mss ith row (or colum, its symettric) but delete the gamma entry with [:1]
##                                             #  Msc ith row will have sin(betai + theta i ) cos(betaj +thetaj)
##                                             
#    cij = k**2/ (1-N*k**2) * r *cij
#    cij = np.delete(cij,i)
##    
#    cij = np.sum(cij)/betas[i]
#                                        
##    cij =0
##    
##    
##    for j in range(0,N):
##        if j==i:
##            continue
##       
##        ss = sin1_sin2(t,betas[j],betas[i],theta[j],theta[i])
##        cs = cos1_sin2(t,betas[j],betas[i],theta[j],theta[i])
##        
##        
##        
##        cij += r[j]*betas[j]**2*( -ss/Qe[j]  + cs )
#    
#    rdot = cri*r[i] + crho*rho + cij
#    
#    return rdot
    

def ridot(t,i,rho,phi,r,theta):
    
    ss = sin1_sin2(t,betas[i],gamma,theta[i],phi)
    sc = sin1_cos2(t,betas[i],gamma,theta[i],phi)
    
    
    cri = -betas[i]*(1.0 + k**2/(1.0-N*k**2) ) /Qe[i]/2.0
    
    crho = k*gamma**2/(1-N*k**2)/betas[i]*( ss/Q0 - sc)
    
    cij =0
    for j in range(0,N):
        if j==i:
            continue
       
        ss = sin1_sin2(t,betas[j],betas[i],theta[j],theta[i])
        cs = cos1_sin2(t,betas[j],betas[i],theta[j],theta[i])
        
        cij += r[j]*betas[j]**2*( -ss/Qe[j]  + cs )
    
    rdot = cri*r[i] + crho*rho + k**2/(1-N*k**2)/betas[i] *cij
    
    return rdot


#def rhodot(t,rho,phi,r,theta,Mss,Msc,Mcc):
#    
#    crho = -gamma/2.0/(1-N*k**2)/Q0
#    
#    
#    cr = np.sum( k**2/gamma/(1-N*k**2)* r * betas**2*( Mss[0][1:]/Qe - Msc[0][1:] ) )
#    
##    cr = 0
##    for j in range(0,N):    
##        ss = sin1_sin2(t,gamma,betas[j],phi,theta[j])
##        sc = sin1_cos2(t,gamma,betas[j],phi,theta[j])
##    
##        cr += r[j]*betas[j]**2*( ss/Qe[j]  - sc )
##        
#    
#    rhodot = crho*rho + cr
#    
#    return rhodot
def rhodot(t,rho,phi,r,theta):
    
    crho = -gamma/2.0/(1-N*k**2)/Q0
    
    cr = 0
    for j in range(0,N):    
        ss = sin1_sin2(t,gamma,betas[j],phi,theta[j])
        sc = sin1_cos2(t,gamma,betas[j],phi,theta[j])
    
        cr += r[j]*betas[j]**2*( ss/Qe[j]  - sc )
        
    
    rhodot = crho*rho +k/gamma/(1-N*k**2) *cr
    
    return rhodot
    


    
def phidot(t,rho,phi,r,theta,Mss,Msc,Mcc):
    
    const = gamma*N*k**2/2.0/(1-N*k**2)
    
    c2 = np.sum(  k/(1-N*k**2)/rho/gamma *  r *  betas**2 * ( Msc[:,0][1:]/Qe - Mcc[0][1:] ) )
    
#    c1 = 0
#    for j in range(0,N):    
#        sc = cos1_sin2(t,gamma,betas[j],phi,theta[j])
#        cc = cos1_cos2(t,gamma,betas[j],phi,theta[j])
###    
#        c1 += k*betas[j]**2/(1-N*k**2)* r[j]* ( sc/Qe[j] - cc  )
#        
    #print c1/rho/gamma
    #print c2
#    
    phidot = const +  c2#/rho/gamma
    
    return phidot
    
def thetaidot(t,i,rho,phi,r,theta,Mss,Msc,Mcc):
    
    #cs = cos1_sin2(t,betas[i],gamma,theta[i],phi)
    #cc = cos1_cos2(t,betas[i],gamma,theta[i],phi)
    
    cs = Msc[0][i+1]
    cc = Mcc[0][i+1]
    
    const = betas[i]*k**2/2.0/(1-N*k**2)
    
    c1 = k*gamma**2/(1-N*k**2)/betas[i]*( cs/Q0 - cc  )

    
    cij =  r / r[i]  / betas[i] * k**2/(1-N*k**2) *betas**2*( -Msc[:,i+1][1:]/Qe + Mcc[i+1][1:] ) 
    
    cij = np.delete(cij,i)
    cij = np.sum(cij)
    #print cij
    
#    cij =0
#    for j in range(0,N):
#        if j==i:
#            continue
#        
#        sc = sin1_cos2(t,betas[j],betas[i],theta[j],theta[i])
#        cc = cos1_cos2(t,betas[j],betas[i],theta[j],theta[i])
#        
#        cij += r[j]*betas[j]**2*( -sc/Qe[j]  + cc )
#    
    thetadot = const + rho/r[i] * c1 + cij 
    
    return thetadot
    

    ############

    
    ########################
    
### set ICS


#theta0 = 0
#phi0 = 0
#
#rho0 = 1.0
#r0 = 1000.
#
#
#Tl = 1000.0
#dT =1.0
#Time = np.arange(0,Tl,dT)
#
#
#
#ri0 = np.ones(N)*r0
#thetai0 = np.ones(N)*theta0
#
#Ri = np.array([ri0])
#Thetai = [thetai0]
#Rho = [rho0]
#Phi = [phi0]
#
#
#rinow = ri0
#thetainow = thetai0
#rhonow = rho0
#phinow = phi0

print "Simulating in the rotating frame..."
#rinextlist = np.zeros(N)
#thetainextlist = np.zeros(N)

    ### set ICS

#Theres a bug in ridot and rho dot or child functinos! 
#But this isn't in the ridot rhodot functions from RframeNelectrons.py!!! (currently using them...but slow)


theta0 = 0.0
phi0 = 0
rho0 = 1.0
r0 = 100.0


Tl = 10000.0
dT =1.0
Time = np.arange(0,Tl,dT)

tstart = time.time()
for num in range(1,2,1):

    N=num
    print N," electrons..."
    #N=1

    Q1 = 1e10#np.float('inf')
    Qe = np.ones(N)*Q1
    Q0=1000
    beta = 1.0
    sigma = 1e-3
    betas =np.ones(N) #beta+np.random.normal(0,sigma,N)#*np.ones(N)

    gamma = 1.0
    k = 1e-5

    wz = 150e6 # collission / axial bounce frequency
    wo = 30e9   # cyclotron frequency
    wr = wo/gamma # reference frequency

    Tcol = gamma*wo/wz

    ri0 = np.ones(N)*r0
    thetai0 = np.ones(N)*theta0 #np.random.uniform(0,2*np.pi,size=N) #np.ones(N)*theta0 #2*np.pi*np.random.uniform(size=N)
    
    Ri = np.array([ri0])
    Thetai = [thetai0]
    Rho = [rho0]
    Phi = [phi0]


    rinow = ri0
    thetainow = thetai0
    rhonow = rho0
    phinow = phi0

    rinextlist = np.zeros(N)
    thetainextlist = np.zeros(N)
    for t in Time[1:]:
    
        Mss = sinsin_matrix(t,phinow,thetainow)
        Msc = sin1cos2_matrix(t,phinow,thetainow)
        Mcc = coscos_matrix(t,phinow,thetainow)
        
    
        rhonext = rhodot(t,rhonow,phinow,rinow,thetainow)*dT + rhonow
        phinext = phidot(t,rhonow,phinow,rinow,thetainow,Mss,Msc,Mcc)*dT + phinow    
    
        for i in range(0,N):
            #print ridot(t,i,rhonow,phinow,rinow,thetainow)
            #print rdot(t,i,rhonow,phinow,rinow,thetainow)
        
            rinextlist[i] = ridot(t,i,rhonow,phinow,rinow,thetainow)*dT + rinow[i]
            thetainextlist[i] = thetaidot(t,i,rhonow,phinow,rinow,thetainow,Mss,Msc,Mcc)*dT +thetainow[i]
    #print rinextlist

        Ri= np.vstack((Ri,rinextlist))
        #print Ri
        Thetai = np.vstack((Thetai,thetainextlist))
        Rho.append(rhonext) 
        Phi.append(phinext)
        
        rinow = rinextlist
        thetainow = thetainextlist
        rhonow=rhonext
        phinow= phinext
        
       # if t%Tcol ==0:
        #    betas = beta+np.random.normal(0,sigma,N)
        
    
    tfinish = time.time()

    print tstart - tfinish

#x = L.analytical_qit(1)
    #y = L.get_particlei(1)
    #z = L.get_particlei(0)
    
    Ri = np.array(Ri)
    Thetai = np.array(Thetai)
    Rho = np.array(Rho)
    Phi = np.array(Phi)

#t = 3.0*np.pi/2.0#np.arange(0,2*np.pi,np.pi/1000)
#p =0# np.arange(0,2*np.pi,0.1)
#r = np.arange(2.5,3,0.01)
#rho = np.arange(-0.1,0.1,0.01)
#
#(X,Y) = np.meshgrid(r,rho)
#u = rdot(X,p,Y,t)
#v = rhodot(X,p,Y,t)
#q = plt.quiver(X,Y,u,v)
##p = plt.quiverkey(q,coordinates='data',color='r')
#plt.show()

    Data = np.array([Rho,Phi])
    Data = np.vstack((Data,Ri.T))
    Data = np.vstack((Data,Thetai.T))

    fname = "Rframe_N"+str(N)+"_Qo"+str(Q0)+"_Qi"+str(Q1)+"_k"+str(k)+"_rho"+str(rho0)+"_r"+str(r0)+"_phi"+str(phi0)+"_theta"+str(theta0)+"_sigma"+str(sigma)+"_Tcol"+str(Tcol)+"_"+str(tfinish)+".csv"                      
    np.savetxt("Rframe\\accuartebutslow"+fname,Data.T,delimiter=',')

pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print s.getvalue()

R1 = Ri[:,0]
#R2 = Ri[:,1]


Qc = 2*Tl/( np.log(R1[0]) - np.log(R1[-1])  )
ctime = Qc/beta
freetime = Q1/beta
print ctime/freetime

#t = np.arange(0,500,0.001)
#qexact = Col0.ANLib.analytical_qvec(t)[1]


#plt.ion()
#plt.subplot(211)
#
#
#plt.ylabel('Log (Relative Amplitude)')
##plt.plot(L.time,y[0],'ys-')
##plt.plot(L.time,z[0],'gs-')
#
##plt.plot(t,qexact[:,1],'ro-')
##plt.plot(Time,np.log(R1/R1[0]),'bo-')
##plt.plot(Time,np.log(R2/R2[0]),'yo-')
#
##plt.plot(t,qexact[:,0],'yo-')
##plt.plot(Time,Rho,'ko-')
##plt.plot(Time,Ri,'bo-')
#plt.plot(Time,np.log10(Ri/Ri[0]),'ro-')
#plt.plot(Time,np.log(Rho/Rho[0]),'bo-')
#
#plt.subplot(212)
##for i in range(0,N):
#  #  plt.plot(Time,(Thetai[:,i]-Phi)/np.pi,'bo')
#plt.plot(Time,Phi/np.pi,'bo')
#plt.plot(Time,Thetai/np.pi,'ro')
##plt.plot(Time,(Thetai[:,2]-Phi)/np.pi,'bo')
##plt.plot(Time,(Thetai[:,3]-Phi)/np.pi,'bo')
##plt.plot(Time,Thetai[:,1]/np.pi,'yo-')
##plt.plot(Time,Phi/np.pi,'ro-')
#plt.ylabel('Phase/$\pi$')
#plt.xlabel('Time')

#plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(Time, Ri, 'r.',label='Electrons')

ax2 = ax1.twinx()
ax2.plot(Time, Rho, 'b.',label='Cavity')
ax1.set_ylabel('Amplitude')
plt.legend(loc='best')

ax3 = fig.add_subplot(212)
plt.plot(Time,Phi/np.pi,'b.')
plt.plot(Time,Thetai/np.pi,'r.')
ax3.set_ylabel('Phase/$\pi$')
ax3.set_xlabel('Time')

plt.show()
    




