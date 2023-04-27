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

xlabel="$L_x$"; ylabel="$R_{\mathrm{eff}}$ $L_y$/$(R L_x)$"

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

# Axes labels
plt.xlabel(xlabel, fontsize=40)
plt.ylabel(ylabel, fontsize=40)
#plt.ylabel(r'$R_\mathrm{eff}$', {'color': 'C0', 'fontsize': 20})

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



# Line styles
lines = ["-","--","-.",":","-","--","-.",":"]
colors =['r','g','b','m','orange','deepskyblue','brown','grey']
#colors =['k','r','g','b','m','orange','deepskyblue','brown','grey']
#colors =['k','r','g','b','m','orange','brown','grey']
linecycler = cycle(lines)




# NOW PLOT
# ------------------------------------------------------------- 
i = 0
#for Ly in [10,20,50,100,500]:
import glob
import re
print('here')
path = 'Rratio_vs_Lx_fixed_Ly*_zighex.dat'
List=[]
for datafile in glob.glob(path):
     Ly = re.sub('\D', '', datafile) # Syntax: re.sub(<match string>,<new string,<whole string>),
                                 # \D: any character that's not a digit
                                 # \d: any character that's a digit
     #print('datafile, Ly=',datafile, Ly)
     List.append(int(Ly)) # convert string to integer

List = sorted(List)
#print('List=',List)
#exit() 
for Ly in List:
     print('Ly=',Ly)
     dir='.' 
     #dir='../DATAFILES'
     datafile = '{}/Rratio_vs_Lx_fixed_Ly{}_zighex.dat'.format(dir,Ly)
     data = np.loadtxt(datafile)
     print('datafile =',datafile)
     ax.plot(data[:,0], data[:,1], linewidth=2.0*(1+i/4), linestyle=lines[i], color=colors[i], label='$L_y$ = {}'.format(Ly))


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
ax.legend(loc='lower left', prop=dict(size=30))
#ax.legend(loc='upper right', bbox_to_anchor = (0.98, 0.15), prop=dict(size=30))
#plt.grid(True)
#ax.legend(loc='upper right', prop={'size':6})
#plt.legend(loc='upper right', bbox_to_anchor = (0.1, 0.99), prop={'size':30})
#plt.legend(loc='upper right', borderpad=2)
#plt.legend(loc='upper left', bbox_to_anchor = (0.5, 0.5))
#plt.legend(loc='upper left')
#plt.legend()

# Annotation 
xannote=150; yannote=2.95 # absolute coordinates
#bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2)
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.9)
ax.text(xannote, yannote, "Hexagonal \n(zigzag)", ha="center", va="center", size=40, bbox=bbox_props)

# Adjust layout
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15,  top=0.95)
# Check: https://stackoverflow.com/questions/18619880/matplotlib-adjust-figure-margin

# SAVE AND SHOW PLOT
# -----------------------------------
fname='Rratio_vs_Lx_fixed_Ly_zighex'
ax.set_rasterized(True)
plotfile='{}.png'.format(fname)
#plt.tight_layout()
#plt.savefig('test.eps', format='eps', dpi=1000)
#plt.tight_layout(True)
plt.savefig("{}".format(plotfile))
plt.savefig("{}.eps".format(fname))
plt.show()

