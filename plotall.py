# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def skipper(fname, linesToSkip):#skip the starting lines for np load (quicker than standard native library)
        with open(fname) as fin:
            for i, line in enumerate(fin): #for each line
                if i > linesToSkip-1: #at least 6th line pleae
                    yield line  #return that row please

# hdata = np.loadtxt(skipper('Htmaxoutput.dat',0))
# hdataa = np.loadtxt(skipper('Hrmaxoutput.dat',0))
# altcdata = np.loadtxt(skipper('AltCtmaxoutput.dat',0))
# altcdataa = np.loadtxt(skipper('AltCrmaxoutput.dat',0))
# NCMdata = np.loadtxt(skipper('NCMtmaxoutput.dat',0))
# NCMdataa = np.loadtxt(skipper('NCMrmaxoutput.dat',0))
# NMSdata = np.loadtxt(skipper('NMStmaxoutput.dat',0))
# NMScdataa = np.loadtxt(skipper('NMSrmaxoutput.dat',0))

hdata = 
hr = hdataa[:,0]
htemp = hdataa[:,1]
ht = hdata[:,0]
htempt = hdata[:,1]

altcr = altcdataa[:,0]
altctemp = altcdataa[:,1]
altct = altcdata[:,0]
altctempt = altcdata[:,1]

NCMr = NCMdataa[:,0]
NCMtemp = NCMdataa[:,1]
NCMt = NCMdata[:,0]
NCMtempt = NCMdata[:,1]

NMSr = hdataa[:,0]
NMStemp = hdataa[:,1]
NMSt = hdata[:,0]
NMStempt = hdata[:,1]


plt.figure()
plt.plot(hr,htemp, label = 'Regular case',color = 'black')
plt.plot(altcr,altctemp, label = 'Alt C case',color = 'red',linestyle='dotted')
plt.plot(NCMr, NCMtemp, label = 'No Conjoined Metal melting case',color = 'green',linestyle='dashed')
plt.plot(NMSr, NMStemp,label = 'No Metal Sulfides case',color = 'blue',linestyle='dashdot')
plt.legend()
plt.ylim((1225,1230))
plt.title('Temp of planetesimal over time')
plt.xlabel('Time (Myr)')
plt.ylabel('Temp (K)')
#plt.savefig('Temp of planetesimal over time - regular case.png')
plt.show()

plt.figure()
plt.plot(ht,htempt, label = 'Regular case',color = 'black')
plt.plot(altct,altctempt, label = 'Alt C case',color = 'red',linestyle='dotted')
plt.plot(NCMt, NCMtempt, label = 'No Conjoined Metal melting case',color = 'green',linestyle='dashed')
plt.plot(NMSt, NMStempt,label = 'No Metal Sulfides case',color = 'blue',linestyle='dashdot')
plt.title('Temp of planetesimal over radius - regular case')
plt.xlabel('Radius (M)')
plt.ylabel('Temp (K)')
#plt.savefig('Temp of planetesimal over radius - regular case.png')
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


