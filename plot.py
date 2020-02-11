# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def skipper(fname, linesToSkip):#skip the starting lines for np load (quicker than standard native library)
        with open(fname) as fin:
            for i, line in enumerate(fin): #for each line
                if i > linesToSkip-1: #at least 6th line pleae
                    yield line  #return that row please

data = np.loadtxt(skipper('tmaxoutput.dat',0))
dataa = np.loadtxt(skipper('rmaxoutput.dat',0))
r = dataa[:,0]
temp = dataa[:,1]
t = data[:,0]
tempt = data[:,1]

# print('What plot would you like to make?')
# print('(1) Temp vs Radius')
# print('(2) Temp vs Time')
# val = input()


plt.figure()
plt.plot(r,temp)
plt.title('Temp of planetesimal over time - regular case')
plt.xlabel('Time (Myr)')
plt.ylabel('Temp (K)')
plt.savefig('Temp of planetesimal over time - regular case.png')
plt.show()
plt.figure()

plt.plot(r,temp)
plt.title('Temp of planetesimal over radius - regular case')
plt.xlabel('Radius (M)')
plt.ylabel('Temp (K)')
plt.savefig('Temp of planetesimal over radius - regular case.png')
plt.show()
# if val == str(1):
#     plt.figure()
#     plt.plot(r,temp)
#     plt.title('Temp of planetesimal over radius')
#     plt.xlabel('Radius (M)')
#     plt.ylabel('Temp (K)')
#     plt.show()

# elif val == str(2):
#     plt.figure()
#     plt.plot(r,temp)
#     plt.title('Temp of planetesimal over time')
#     plt.xlabel('Time (Myr)')
#     plt.ylabel('Temp (K)')
#     plt.show()

# else:
#     print('invalid input')
#plt.savefig('TotalKeVsTime- vs')


