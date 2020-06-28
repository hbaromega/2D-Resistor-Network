import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import subprocess
from matplotlib.cbook import get_sample_data
import os,subprocess
from itertools import cycle


# SETTING PLOT ATTRIBUTES
# --------------------------------------------------------------------------
# Set borderline width
mpl.rcParams['axes.linewidth'] = 2.0


# Line styles
lines = ["-","--","-.",":","-","--","-.",":"]
colors =['r','g','b','m','orange','deepskyblue','brown','grey']
#colors =['k','r','g','b','m','orange','deepskyblue','brown','grey']
#colors =['k','r','g','b','m','orange','brown','grey']
linecycler = cycle(lines)




# NOW PLOT
# ------------------------------------------------------------- 
List=['set_rectangular.dat','set_hexagonal.dat']
geo=['Rectangular','Hexagonal']
R_th=[200,271]


# PLOT 1
#plt.figure(100)
plt.xlabel=('Power supply voltage (Volts)')
plt.ylabel=('$R_{\mathrm{int}}$') #, fontsize=40)
for i in [0,1]:
     datafile=List[i]
     df=pd.read_csv(datafile, delim_whitespace=True)
     print('df=',df)
     print("df['#V']=",df['#V'])
     exit()
     print('i,datafile=',i,datafile)
     dir='.' 
     data = np.loadtxt(datafile)
     R_ex=data[:,1]
     R_int=1.0/(1.0/R_ex-1.0/R_th[i])
     print('R_th=',R_th[i])
     print('R_int=',R_int)
     print('R_ex=',R_ex)
     plt.plot(data[:,0], R_int, linewidth=2.0*(1+i/4), linestyle=lines[i], color=colors[i], label=geo[i])
plt.legend(loc='upper left', prop=dict(size=20))
fname='Rint_vs_V'
plt.tight_layout()
plt.savefig("{}.png".format(fname))
plt.savefig("{}.eps".format(fname))
plt.show()

exit()
# PLOT 2
plt.figure(101)
plt.xlabel=("Power supply voltage (Volts)")#, fontsize=40)
plt.ylabel=('$R_{\mathrm{expt}}$')#, fontsize=40)
for i in [0,1]:
     datafile=List[i]
     print('i,datafile=',i,datafile)
     dir='.' 
     data = np.loadtxt(datafile)
     plt.plot(data[:,0], data[:,1], linewidth=2.0*(1+i/4), linestyle=lines[i], color=colors[i], label=geo[i])
plt.legend(loc='upper left', prop=dict(size=20))
fname='Rexp_vs_V'
plt.savefig("{}.png".format(fname))
plt.savefig("{}.eps".format(fname))


plt.show()

