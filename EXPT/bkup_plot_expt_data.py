import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import subprocess
from matplotlib.cbook import get_sample_data
import os,subprocess
from itertools import cycle
from sklearn.linear_model import LinearRegression 

def estimate_coef(x, y): 
    # number of observations/points 
    n = np.size(x) 
  
    # mean of x and y vector 
    m_x, m_y = np.mean(x), np.mean(y) 
  
    # calculating cross-deviation and deviation about x 
    SS_xy = np.sum(y*x) - n*m_y*m_x 
    SS_xx = np.sum(x*x) - n*m_x*m_x 
  
    # calculating regression coefficients 
    b_1 = SS_xy / SS_xx 
    b_0 = m_y - b_1*m_x 
  
    return(b_0, b_1)

def plot_regression_line(x, y, b): 
    # plotting the actual points as scatter plot 
    #plt.scatter(x, y, color = "m", 
     #          marker = "o", s = 30) 
  
    # predicted response vector 
    y_pred = b[0] + b[1]*x 
    return y_pred

    # plotting the regression line 
    #plt.plot(x, y_pred, color = "g") 
  
    # putting labels 
    #plt.xlabel('x') 
    #plt.ylabel('y') 
  
    # function to show plot 
    #plt.show()  





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



# NOW PLOT
# ------------------------------------------------------------- 
path = 'set*.dat'
List=['set_rectangular.dat','set_hexagonal.dat']
geo=['Rectangular','Hexagonal']
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
     Vset = data[:,0]
     Iset = data[:,1]
     print('Vset=',Vset)
     print('Iset=',Iset)
     R_ex=[x/y for x,y in zip(Vset,Iset)]
     # R=[x/y for x,y in zip(data[:,0],data[:,1])] # Also works
     print('R_ex=',R_ex,data[:,2])
     R_int=[1.0/(1.0/x-1.0/R_th[i]) for x in R_ex]
     #R_int=1.0#/(1.0/R_ex-1.0/R_th[i])
     print('R_th=',R_th[i])
     print('R_int=',R_int)
     ax.scatter(Vset, R_int, marker=markers[i], color=colors[i], label=geo[i])
     #plt.plt(Vset, R_int, marker=markers[i], color=colors[i], label=geo[i])

     # Linear regression:
     #model =  LinearRegression()
     #results = model.fit(Vset,R_int)

     #R_int = np.array(R_int)
     Vset = np.array(Vset)

     # estimating coefficients 
     b = estimate_coef(Vset,R_int)

     print("Estimated coefficients:\n b_0 = {}; b_1 = {}".format(b[0], b[1])) 
  
     # plotting regression line 
     #plot_regression_line(Vset, R_int, b)
     sum0=sum0+b[0]
     sum1=sum1+b[1]
b0=sum0/2.0; b1=sum1/2.0
y_pred = b0 + b1*Vset 
plt.plot(Vset, y_pred, color = "b")

plt.legend(loc='upper left', prop=dict(size=20))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15,  top=0.95)
ax.set_xlabel("Power supply voltage (Volts)")
ax.set_ylabel('$R_{\mathrm{int}}$ (Ohms)') #, fontsize=40)
#plt.ylabel=('$R_{\mathrm{int}}$') #, fontsize=40)

fname='Rint_vs_V'
plt.tight_layout()
plt.savefig("{}.png".format(fname))
plt.savefig("{}.eps".format(fname))
plt.show()

#exit()
# PLOT 2
plt.figure(101)
plt.xlabel("Power supply voltage (Volts)")#, fontsize=40)
plt.ylabel('$R_{\mathrm{expt}}$ (Ohms)')#, fontsize=40)
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

