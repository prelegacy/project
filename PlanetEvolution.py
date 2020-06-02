# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:50:17 2020

@author: Chris
"""


import numpy as np
def skipper(fname, linesToSkip):#skip the starting lines for np load (quicker than standard native library)
        with open(fname) as fin:
            for i, line in enumerate(fin): #for each line
                if i > linesToSkip-1: #at least 6th line pleae
                    yield line
#Consider only reading from file what is neccesary
# data = np.loadtxt(skipper('regalltemps.txt',0))  

r = np.arange(0,1000e3,101)
dr = np.mean(np.diff(r))
trange = dr**2 * 0.01
# print('trange =',trange)
t = np.arange(200e6,500e6,trange)
dt = np.mean(np.diff(t))
# print(dr,dt,dt/dr**2)
######## Accretion code 

#Thermal conductivity in J/(yr*m*K)
#3 values of K are for the regolith, normal chondrite, chondrite after melting
K = np.array([3.2e7,3.2e7,12.6e7])

#Regolith thickness (in number of radial space steps)
Reg = 2

#Average bulk density of chondrite
p = 3440

#Specific heat capacities of [silicates, metal, sulfide, conjoined grains]
c = np.array([892,598,699,616])

#Weight fraction of phases in chondrite: [silicates, metal, sulfide, conjoined grains]
P = [0.76,0.05,0.03,0.16]

#Initial temp throughout asteroid and initial temp of newly accreting material
data = np.loadtxt('regalltemps.txt',usecols=0)
data = data[1:-1]/1000000
a = int(np.min(np.argwhere((data > 200.0)&(data <201))))

data = np.loadtxt(skipper('regalltemps.txt',a))
init = np.array(data[0,1:]) 
# print(len(init))
# print('init',init)
#Import remaining values for asteroid 
#######
######
######
#External temp throughout time
bry = 180 

#Latent heat of fusion for phases in J/kg, in steps 
Hin=Hstart = np.array([[2000,20000,40000,120000,200000],
       [252000,0,0,0,0],
       [358000,0,0,0,0],
       [37978,4250,4647,9308,214817]])
#Total latent heat of fusion (Htot) is divided into five parts.
#Each part is the value of Htot multiplied by the percentage of the phase that melts at that step
#So if first melting step caused 5% melting silicates, 5% of Htot is put as Hin(1,1)

#Latent heat of fusion, same as Hstart, but different format. THis is for use in setting up Hin values. Put in as an array, in the order silicates (5 steps), metal, sulfide, conjoined (5 steps)
Hstart_imp = np.array([20000,20000,40000,120000,20000,252000,358000,136282,23313,27815,65331,52260])

#Tempearture (K) at each melting step. Should be 5 for silicates and conjoined grains, and one for isolate metal and sulfide grains
M = np.array([[1353,1393,1483,1753],
              [1809,0,0,0,0],
              [1463,0,0,0,0],
              [1236,1400,1500,1600,1702]])

#Setup accretion conditions 
#no of accretion steps - set to 50 for now?
Z = 1
#Final radius 
R = 1000e3

#Values of radius (m) at each accretion step
Rvals = np.arange(0,R,101)
#FInal number of space steps once radius is at ifnal sice
r_step_tot = len(Rvals)

#input 
#accretion time imported from file 
t_acc = 200e6
#accretion duration myr needs to be discussed
t_dur = 0.35e6
#final time to compute values 
tfin = 500e6
#Values of time (yr) at each accretion step
tvals = np.arange(t_acc,t_acc+t_dur)

#Number of time steps used between accretion steps
t_step_dur = 201
#no of timesteps used after accretion
t_step_fin = (tfin/1e4+1)

#Setup matrix for values of N (space steps) at each accretion step
N = np.zeros((1,Z))

#Setup matrix for values of J (time steps) at each accretion step
J = np.zeros((1,Z))

#Total number of time steps used 
t_step_tot = t_step_dur*Z + t_step_fin

#Matrix temps_time, column 1 has time at each timestep, remaining are temps at each radial position
temps_time = np.zeros((int(t_step_tot),int(r_step_tot)+1))

#Reynolds functions - general algorithm used for each melting step
#Modified from Eynolds (1996)
#T = T(j+1,n)
#M = melting temp
#H = residual heat of fusion
#Hin = initial heat of fusion
#P = wt fraction of melting phase
#c = specific heat capacity of that phase 

def Reynolds(T,M,H,Hin,P,c):
    if T < M:
        if H >= Hin:
            H = Hin
        else:
            T2 = T + P * ((Hin - H)/c)
            if T2 < M:
                T = T2
                H = Hin
            else:
                T = M
        if H <= 0:
            H = 0
        else: 
            T2 = T - P * ((H/c))
            if T2 > M:
                T = T2
                H = 0
            else:
                T = M
    return [T,H]

def HeatEqn_Grad_A(k,K,Reg,p,c,r,t,P,init,bdry,Hin,Hstart,M,acc_con):
    #Already assigned above 
    #N = len(r)
    #J = len(t)
    #dr = np.mean(np.diff(r))
    #dt = np.mean(np.diff(t))
    
    #Set up matrix of J N
    T = np.array([init])
    # print(T.shape,'before')    
    T = np.array([np.append(T,np.ones(len(r)-len(init))*180)])
    
    
    return T
# initi = len(init(1,:))
# for i in :
#     for j in len(init(:,1)):
#         a


# print(T.shape,'after')
# T = np.append(T,[init/180],axis=0)

    
#Heat equation grad a
#Consider having outputs written to a file to read from as well

