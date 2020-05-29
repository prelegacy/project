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
                    
# data = np.loadtxt(skipper('regalltemps.txt',0))  
r = np.arange(0,100e3,101)
dr = np.mean(np.diff(r))
trange = dr**2 * 0.01
print('trange =',trange)
t = np.arange(2.853e6,500e6,trange)
dt = np.mean(np.diff(t))
print(dr,dt,dt/dr**2)

bry = 180 

Hin = np.array([[2000,20000,40000,120000,200000],
       [252000,0,0,0,0],
       [358000,0,0,0,0],
       [37978,4250,4647,9308,214817]])
M = np.array([[1353,1393,1483,1753],
              [1809,0,0,0,0],
              [1463,0,0,0,0],
              [1236,1400,1500,1600,1702]])

def q(x):
    return x