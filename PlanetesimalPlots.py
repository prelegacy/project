#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 24 16:40:30 2020

@author: ccnew1
"""

import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator, FormatStrFormatter

def skipper(fname, linesToSkip):#skip the starting lines for np load (quicker than standard native library)
        with open(fname) as fin:
            for i, line in enumerate(fin): #for each line
                if i > linesToSkip-1: #at least 6th line pleae
                    yield line
                    
data = np.loadtxt(skipper('regalltemps.txt',0))               
# data = np.loadtxt(fname='regalltemps.txt',skiprows=0)
radius = data[0,1:]/1000
time = data[1:,0]/1000000
temps = data[1:,1:]

# Set up grid and test data
x = radius
y = time[0:-2]


fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.set_xlabel=('Radius (m)')
# ax.set_ylabel=('Time (years)')
# ax.set_zlabel=('Temperature (K)')
X, Y = np.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
# surf = ax.contourf(X, Y, temps,cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# ax.set_zlim(180,np.max(temps))
surf = plt.contourf(Y,X,temps[0:-2])
plt.title('2.53e23 Al abundance, canonical ratio, 0.35Myr acc time')
plt.xlabel('Time (Myr)')
plt.ylabel('Radius (km)')
fig.colorbar(surf,label='Temperature (K)')
if np.max(temps)>1700:
    plt.clim(1700)

# plt.savefig('2.53e23, canonical ratio.PNG',bbox_inches='tight',dpi=100)





