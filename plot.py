#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 24 16:40:30 2020

@author: ccnew1
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def skipper(fname, linesToSkip):#skip the starting lines for np load (quicker than standard native library)
        with open(fname) as fin:
            for i, line in enumerate(fin): #for each line
                if i > linesToSkip-1: #at least 6th line pleae
                    yield line
                    
data = np.loadtxt(skipper('regalltemps.txt',0))               

radius = data[0,1:]
time = data[1:,0]
temps = data[1:,1:]

# Set up grid and test data
x = radius
y = time


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.set_xlabel=('Radius (m)')
# ax.set_ylabel=('Time (years)')
# ax.set_zlabel=('Temperature (K)')
# X, Y = np.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
# surf = ax.plot_surface(X, Y, temps,cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)

# fig.colorbar(surf, shrink=0.5, aspect=5)
plt.figure(figsize=(10,10))
plt.plot(temps[0,:],temps[1,:])