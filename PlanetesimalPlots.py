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
                    
# data = np.loadtxt(skipper('regalltemps.txt',0))  

big = 3

if big == 1:
    time = 3     
    data = np.loadtxt('Largetmp.txt',delimiter=',')           
    # data = np.loadtxt(fname='regalltemps.txt',skiprows=0)
    T = data[:,:]
    
    Planetdensity = 3440
    Planetradius = 500e3
    p = 3440
    radius = 100e3
    # fAl = fAl* pa*planetesimals(Planetdensity,Planetradius,p,radius)
    # enhancement = pa *planetesimals(Planetdensity,Planetradius,p,radius)
    EAl = 6.4154e-13
    EFe = 4.87e-13
    LifeAl = 1.07e6
    LifeFe = 3.49e6
    r = np.arange(int(0),int(Planetradius),int(1000))
    
    dr = np.mean(np.diff(r))
    trange = dr**2 * 0.01
    # print('trange =',trange)
    t = np.arange(time*1e6,100e6,trange)
    dt = np.mean(np.diff(t))
    pdist = np.array([120/280 *r[-1],80/280 *r[-1],14/280 *r[-1],(280-214)/280 *r[-1]])
    core = int(np.min(np.argwhere((pdist[0] < r))))
    mantle = int(np.min(np.argwhere((pdist[1]+pdist[0] < r))))
    lower_crust = int(np.min(np.argwhere((pdist[2]+pdist[1]+pdist[0] < r))))
    upper_crust = lower_crust+1
    record = 20
    x = r/1000
    y = t/1e6
    
    
    fig = plt.figure(figsize=(10,10))
    X, Y = np.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
    surf = plt.contourf(Y,X,T,15)
    CS = plt.contour(Y, X,T,levels = [1650,1750],
                     colors=('k',),linestyles=('-',),linewidths=(2,))
    plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=14)
    plt.hlines(r[core]/1000, y[record], np.max(Y),color='red')
    plt.hlines(r[mantle]/1000, y[record], np.max(Y),color='red')
    plt.axvline(y[record],0,500,color='red')
    plt.title('2-zone model,Canonical Al & Fe, accretion at '+ str(t[0]/1e6))
    # if case == 2:
    #     plt.title('4-zone model,Canonical Al & Fe, accretion at '+ str(t[0]/1e6))
    plt.xlabel('Time (Myr)')
    plt.ylabel('Radius (km)')
    cbar=fig.colorbar(surf,label='Temperature (K)')
    
    cbar.ax.invert_yaxis()
    # if wvalue == 1:
    #     plt.savefig('init '+str(time)+' myr acc,'+str(pa)+' too little Plot.PNG',bbox_inches='tight',dpi=100)
    # elif wvalue ==2:
    #     plt.savefig('non init '+str(time)+' myr acc,'+str(pa)+' too little Plot.PNG',bbox_inches='tight',dpi=100)
      
    plt.show()
if big == 2:
    data = np.loadtxt('regalltemps.dat') 
    Time = data[1:,0]
    Radius = data [0,1:]
    T = data[1:,1:]
    # r = np.arange(0,100000,101)
    # dr = np.mean(np.diff(r))
    # trange = dr**2 * 0.01
    # t = np.arange(Time[0]*1e6,500e6,trange)
    x = Radius/1000
    y = Time/1e6
    X, Y = np.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
    fig = plt.figure(figsize=(10,10))
    CS = plt.contour(Y, X,T,levels = [1700],
                    colors=('k',),linestyles=('-',),linewidths=(2,))
    plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=14)
    surf = plt.contourf(Y,X,T,100)

    plt.title('Gradual Accretion 500-km model,growing Al & Fe, accretion at '+ str(Time[0]))
    # # if case == 2:
    # #     plt.title('4-zone model,Canonical Al & Fe, accretion at '+ str(t[0]/1e6))
    plt.xlabel('Time (Myr)')
    plt.ylabel('Radius (km)')
    cbar=fig.colorbar(surf,label='Temperature (K)')
    cbar.ax.invert_yaxis()
if big ==3:
    data = np.loadtxt('regalltemps.dat') 
    Time = data[1:,0]
    Radius = data [0,1:]
    T = data[1:,1:]
    # r = np.arange(0,100000,101)
    # dr = np.mean(np.diff(r))
    # trange = dr**2 * 0.01
    # t = np.arange(Time[0]*1e6,500e6,trange)
    x = Radius/1000
    y = Time/1e6
    X, Y = np.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
    fig = plt.figure(figsize=(10,10))
    CS = plt.contour(Y, X,T,levels = [1700],
                    colors=('k',),linestyles=('-',),linewidths=(2,))
    plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=14)
    surf = plt.contourf(Y,X,T,100)

    plt.title('Gradual Accretion 500-km model,scaled Al & Fe, 1 Myr accretion from '+ str(y[0]) +' Myr')
    # # if case == 2:y
    # #     plt.title('4-zone model,Canonical Al & Fe, accretion at '+ str(t[0]/1e6))
    plt.xlabel('Time (Myr)')
    plt.ylabel('Radius (km)')
    cbar=fig.colorbar(surf,label='Temperature (K)')
    cbar.ax.invert_yaxis()
    fig = plt.figure(figsize=(10,10))
    CS = plt.contour(Y, X,T,levels = [1200],
                        colors=('k',),linestyles=('-',),linewidths=(2,))
    plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=14)
    surf = plt.contourf(Y,X,T,100)
    plt.title('Zoom Gradual Accretion 500-km model,scaled Al & Fe, 1 Myr accretion from '+ str(y[0]) + 'Myr')
    # # if case == 2:
    # #     plt.title('4-zone model,Canonical Al & Fe, accretion at '+ str(t[0]/1e6))
    plt.xlabel('Time (Myr)')
    plt.ylabel('Radius (km)')
    cbar=fig.colorbar(surf,label='Temperature (K)')
    cbar.ax.invert_yaxis()
    plt.xlim(2,10)
    



