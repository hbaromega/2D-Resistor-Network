import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import subprocess
from matplotlib.cbook import get_sample_data
import os,subprocess
from itertools import cycle


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



ylabel="$R_{\mathrm{eff}}$"

flag_default = 0
#flag_default=int(input('Enter flag_default '))
print('flag_default=',flag_default)
if (flag_default==1):
    xmin, xmax = [0,100]
    ymin, ymax = [0,20]
elif (flag_default==2):
    xmin,xmax=[float(item) for item in input('Enter x-axis min and max\n').split()] 
    ymin,ymax=[float(item) for item in input('Enter y-axis min and max\n').split()]


# Set plot features 
#fig, ax = plt.subplots()
fig, ax = plt.subplots(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')



# Set tick features 
plt.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
plt.tick_params(axis='both',which='minor',width=2,length=5)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
#ax.set_axis_bgcolor('grey')
#ax.tick_params(axis="ylabel",direction="in", pad=-22)
#ax.tick_params(axis="ylabel",direction="in")
#ax.tick_params(axis="xlabel",direction="in")


# Set axis-limits
if flag_default != 0:
   plt.xlim(xmin,xmax)
   plt.ylim(ymin,ymax)
plt.tick_params(axis='both', which='major', labelsize=30)

#
flag_inset=0
if flag_inset == 1:
 # Inset plot setup
  #left, bottom, width, height = [0.4, 0.5, 0.6, 0.7]
  #left, bottom, width, height = [0.2, 0.55, 0.35, 0.35]
  left, bottom, width, height = [0.2, 0.6, 0.3, 0.3]
  axin = fig.add_axes([left, bottom, width, height])
  xlabel="$1/L_y$"
  xins1, xins2, yins1, yins2 = [0, 100, 0, 100]
  #axin = fig.add_axes([0.2, 0.55, 0.35, 0.35])
  axin.set_xlim(xins1,xins2)
  axin.set_ylim(yins1,yins2)
  xinlabel="$L_y$"
  #yinlabel=ylabel
  yinlabel=""
  axin.set_xlabel('{}'.format(xinlabel), fontsize=30)
  axin.set_ylabel(yinlabel, fontsize=30)
  axin.tick_params(axis='both', which='major', labelsize=20)
  #axin.major_xticks = np.arange(0, 101, 50) #seems doesn't work
  #axin.major_yticks = np.arange(0, 101, 50)  #seems doesn't work
  #axin.xticks(np.arange(min(x), max(x)+1, 1.0))
  # Syntax : axin.xaxis.set_ticks(np.arange(start, end, stepsize))
  axin.xaxis.set_ticks(np.arange(0, 101, 50))
  axin.yaxis.set_ticks(np.arange(0, 101, 50))
else:
  xlabel="$L_y$"

# Main axes labels
plt.xlabel(xlabel, fontsize=40)
plt.ylabel(ylabel, fontsize=40)
#plt.ylabel(r'$R_\mathrm{eff}$', {'color': 'C0', 'fontsize': 20})

# Line styles
lines = ["-","--","-.",":","-","--","-.",":"]
colors =['r','g','b','m','orange','deepskyblue','brown','grey']
#colors =['k','r','g','b','m','orange','deepskyblue','brown','grey']
#colors =['k','r','g','b','m','orange','brown','grey']
linecycler = cycle(lines)




# NOW PLOT
# ------------------------------------------------------------- 
i = 0
#for Lx in [10,20,50,100,500]:
import glob
import re

path = 'Reff_vs_Ly_fixed_Lx*_zighex.dat'
List=[]
for datafile in glob.glob(path):
     Lx = re.sub('\D', '', datafile) # Syntax: re.sub(<match string>,<new string,<whole string>),
                                 # \D: any character that's not a digit
                                 # \d: any character that's a digit
     List.append(int(Lx)) # convert string to integer

List = sorted(List) 
for Lx in List:
     print('Lx=',Lx)
     dir='.' 
     #dir='../DATAFILES'
     datafile = '{}/Reff_vs_Ly_fixed_Lx{}_zighex.dat'.format(dir,Lx)
     data = np.loadtxt(datafile)
     print('datafile =',datafile)
     #ax.plot(data[:,2], data[:,1], linewidth=2.0*(1+i/4), linestyle=lines[i], color=colors[i], label='$L_x$ = {}'.format(Lx))
     ax.plot(data[:,0], data[:,1], linewidth=2.0*(1+i/4), linestyle=lines[i], color=colors[i], label='$L_x$ = {}'.format(Lx))


     # INSET Plot
     # this is an inset axes over the main axes
     #a = plt.axes([.65, .6, .2, .2], facecolor='y')
     # These are in unitless percentages of the figure size. (0,0 is bottom left)
     #axin.plot(data[:,0], data[:,1], linewidth=1.0*(1+i/4), linestyle=lines[i], color=colors[i])

     i += 1


# Gridlines through origin
#ax.axhline(0, linestyle='--', color='k') # horizontal lines
#ax.axvline(0, linestyle='--', color='k') # vertical lines




# LEGEND AND OTHER DECORATIONS
# ------------------------------
# Legend 
#ax.legend(loc='center', bbox_to_anchor = (0.85, 0.9), prop=dict(size=30))
ax.legend(loc='lower right', bbox_to_anchor = (0.98, 0.15), prop=dict(size=30))
#plt.grid(True)
#ax.legend(loc='upper right', prop={'size':6})
#plt.legend(loc='upper right', bbox_to_anchor = (0.1, 0.99), prop={'size':30})
#plt.legend(loc='upper right', borderpad=2)
#plt.legend(loc='upper left', bbox_to_anchor = (0.5, 0.5))
#plt.legend(loc='upper left')
#plt.legend()

# Annotation 
xannote=4; yannote=220 # absolute coordinates
#bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2)
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.9)
ax.text(xannote, yannote, "Hexagonal \n(zigzag)", ha="center", va="center", size=40, bbox=bbox_props)
#ax.text(0, 0.04, "$g$=1.0031, $\mu$=0.7555955, $T_c$=0.2395, $\Delta=$0.5973, $k_0$=4.0", size=20,bbox=bbox_props)

# Adjust layout
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15,  top=0.95)
# Check: https://stackoverflow.com/questions/18619880/matplotlib-adjust-figure-margin

# SAVE AND SHOW PLOT
# -----------------------------------
fname='Reff_vs_Ly_fixed_Lx_zighex'
plotfile='{}.png'.format(fname)
#plt.tight_layout()
#plt.savefig('test.eps', format='eps', dpi=1000)
#plt.tight_layout(True)
plt.savefig("{}".format(plotfile))
ax.set_rasterized(True)
plt.savefig("{}.eps".format(fname))
plt.show()

