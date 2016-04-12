#!/usr/bin/python
######################################################################
# Name:         plot_integrated_parameters
# Author:       A. Marocchino
# Date:			28-10-2015
# Purpose:      plot main bunch integrated parameters and export them in a figure
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
### --- ###



# --- load witness integrated quantities --- #
bunch = np.loadtxt('bunch_integrated_quantity_2.dat')


# --- plot --- #
fig = pyl.figure(1, figsize=(3.25,6.0))
ax1  = pyl.subplot(211)
ax1.plot(bunch[:,0]*.3/1e4,bunch[:,7],      'b-', lw=2, label = '$\sigma_x$ ($\mu$m)')
ax1.plot(bunch[:,0]*.3/1e4,bunch[:,13],     'm-', lw=2, label = '$\epsilon_x$ ($\mu$m)')
ax1.plot(bunch[:,0]*.3/1e4,bunch[:,16]*100, 'c-', lw=2, label = '$\Delta \gamma / \gamma$ * 100')

ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, prop={'size':7.5})
pyl.ylabel('($\mu$m)', fontsize = 8.0)

ax1.xaxis.set_major_locator(MultipleLocator(1.)) 
ax1.xaxis.set_minor_locator(MultipleLocator(.5)) 
pyl.xticks(fontsize=9)
pyl.yticks(fontsize=9)



ax2 = plt.subplot(212)
ax2.plot(bunch[:,0]*.3/1e4,bunch[:,15]/2.,  'b-', label = '$\sigma_x$ ($\mu$m)')
ax2.xaxis.set_major_locator(MultipleLocator(1.)) 
ax2.xaxis.set_minor_locator(MultipleLocator(.5)) 

pyl.xlabel('witness position (cm)', fontsize = 8.0)
pyl.ylabel('Energy (MeV)', fontsize = 8.0)


pyl.subplots_adjust(bottom=0.15,left=0.190)

pyl.savefig('im_witness_integrated_parameter.eps', format='eps')
pyl.savefig('im_witness_integrated_parameter.pdf', format='pdf')
plt.show()





