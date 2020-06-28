import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import subprocess
from matplotlib.cbook import get_sample_data
import os,subprocess
from itertools import cycle
from sklearn.linear_model import LinearRegression 
from sklearn.metrics import mean_squared_error, r2_score





## Set text
# Optionally set font to Computer Modern to avoid common missing font errors
#mpl.rc('font', family='serif', serif='cm10')
#mpl.rc('text', usetex=True)
#plt.rc('text', usetex=True)
#mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
#mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]


# SETTING PLOT ATTRIBUTES
# --------------------------------------------------------------------------
# Set borderline width
mpl.rcParams['axes.linewidth'] = 2.0


# Line styles
lines = ["-","--","-.",":","-","--","-.",":"]
markers = ['o','*', '.', ',', 'x', '+', 'v', '^', '<', '>', 's', 'd']
colors =['r','g','b','m','orange','deepskyblue','brown','grey']
#colors =['k','r','g','b','m','orange','deepskyblue','brown','grey']
#colors =['k','r','g','b','m','orange','brown','grey']
linecycler = cycle(lines)

# Set tick features 
fig, ax = plt.subplots()
plt.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
plt.tick_params(axis='both',which='minor',width=2,length=5)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.xlim(0,15.5)
plt.ylim(0,0.08)


# NOW PLOT
# ------------------------------------------------------------- 
path = 'set*.dat'
List=['set_rectangular.dat','set_hexagonal.dat']
geo=['Square','Hexagonal']
R_th=[200,271]


# PLOT 1
#fig, ax = plt.subplots(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
#plt.figure(100)
sum0=0; sum1=0
print('type=',type(plt.xlabel))
for i in [0,1]:
     datafile=List[i]
     print('i,datafile=',i,datafile)
     dir='.' 
     data = np.loadtxt(datafile)
     x = data[:,0]
     y = data[:,1]
    

     
     # Scikit needs reshaping as if it's 2D (m x n) structure (n=1 here, column vector): 
     x = x.reshape(-1,1)
     y = y.reshape(-1,1)

     print('Vset=',x)
     print('Iset=',y)
     ax.scatter(x, y, marker=markers[i], color=colors[i], label=geo[i])

     # Linear regression:
     model =  LinearRegression()
     model.fit(x, y)

     y_pred = model.predict(x)

     # model evaluation
     rmse = mean_squared_error(y, y_pred)
     r2 = r2_score(y, y_pred)
     slope = model.coef_
     #slope = slope.reshape(-1)
     slope = slope[[0]]
     slope = slope.item()
     #slope = np.asscalar(slope[[0]]) # Works, but deprecated
     slope = float('{:.4f}'.format(slope))
     plt.plot(x, y_pred, color = colors[i])   
     plt.plot(x, y_pred, color = colors[i], label='slope ={}'.format(slope))   


     # printing values
     print('Network:', geo[i])
     print('Slope:', model.coef_)
     print('Intercept:', model.intercept_)
     print('Root mean squared error: ', rmse)
     print('R2 score: ', r2)



plt.legend(loc='upper left', fancybox=True, framealpha=0.5, prop=dict(size=20))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15,  top=0.95)
ax.set_xlabel("Power supply voltage, $V$ (Volts)", fontsize=20)
ax.set_ylabel('Current, $I$ (Amperes)', fontsize=20)
#plt.ylabel=('$R_{\mathrm{int}}$') #, fontsize=40)

fname='I_vs_V_w_linreg_sklearn'
plt.tight_layout()
ax.set_rasterized(True)
plt.savefig("{}.png".format(fname))
plt.savefig("{}.eps".format(fname))
plt.show()



