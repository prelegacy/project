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

r = np.arange(0,100e3,101)
dr = np.mean(np.diff(r))
trange = dr**2 * 0.01
print('trange =',trange)
t = np.arange(2.853e6,500e6,trange)
dt = np.mean(np.diff(t))
print(dr,dt,dt/dr**2)
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
init = 180
init1 = 180



bry = 180 

Hin = np.array([[2000,20000,40000,120000,200000],
       [252000,0,0,0,0],
       [358000,0,0,0,0],
       [37978,4250,4647,9308,214817]])
M = np.array([[1353,1393,1483,1753],
              [1809,0,0,0,0],
              [1463,0,0,0,0],
              [1236,1400,1500,1600,1702]])
#Consider having outputs written to a file to read from as well

data = np.loadtxt('regalltemps.txt',usecols=0)
data = data[1:-1]/1000000
a = int(np.min(np.argwhere((data > 200.0)&(data <201))))

data = np.loadtxt(skipper('regalltemps.txt',a))
temps = np.array(data[0,1:]) 
print(len(temps))
print('temp',temps)