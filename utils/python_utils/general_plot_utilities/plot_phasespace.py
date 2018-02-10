#!/usr/bin/python
######################################################################
# Name:         plot_phasespace
# Author:       A. Marocchino
# Date:			2017-11-02
# Purpose:      plot phase space for architect
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil
import time, datetime
import scipy
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import pylab as pyl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
plt.style.use(os.path.join(os.path.expanduser('~'),'Codes','Python_general_controllers','python_plot','plot_style_ppth.mplstyle'))
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes/Code_Architect/Architect/utils/python_utils/architect_graphycal_unit'))
from read_architect_bin import *
### --- ###

# --- inputs --- #
path = os.getcwd()
read_architect_bin(path,'PS')

# if(len(sys.argv)<2):
# 	print('not enought input arguments')
# 	print('select bunch')
# 	sys.exit()
# bunch_select = int(sys.argv[1])
#
#

# --- plot --- #
fig = pyl.figure(1)
fig.set_size_inches(3.0,3.0,forward=True)
ax1  = pyl.subplot(111)

select_bunch = (np.asarray(var.bunch_id)==1)
select_dcut  = (np.asarray(var.dcut)==1.0)
selected     = np.asarray(select_bunch) & np.asarray(select_dcut)
z0           = np.mean(var.z[selected])
ax1.plot(var.z[selected]-z0,-var.pz[selected], '.', markersize=0.2, lw=1, label=r"Driver")

select_bunch = (np.asarray(var.bunch_id)==2)
select_dcut  = (np.asarray(var.dcut)==1.0)
selected     = np.asarray(select_bunch) & np.asarray(select_dcut)
ax1.plot(var.z[selected]-z0,-(var.pz[selected]-200.), '.', markersize=0.2, lw=1, label=r"Tr. bunch")

#--- labels ---#
ax1.set_xlabel(r'Z ($\mu$m)')
ax1.set_ylabel(r'$\beta \gamma$')
#
#
# ax1.legend(loc=9, ncol=2, prop={'size':7.5})
# ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=3)
# # pyl.subplots_adjust(bottom=0.34,left=0.210)
#
# # ax1.xaxis.set_major_locator(MultipleLocator(1.))
# # ax1.xaxis.set_minor_locator(MultipleLocator(.5))
#
pyl.savefig(os.path.join(path,'im_longitudinal_phasespace.pdf'), format='pdf')
plt.show()
