#!/usr/bin/python
######################################################################
# Name:         plot_integrated_parameters
# Author:       A. Marocchino
# Date:			20-09-2016
# Purpose:      plot main bunch integrated parameters and export them in a figure
# Source:       python
# --- --- --- #
#--- Architect-PYTHON-paths! ---#
#PATH="~/Codes/Code_Architect/Architect/utils/python_utils:$PATH"
#export PATH
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
### --- ###

# --- inputs --- #
if(len(sys.argv)<2):
	print('not enought input arguments')
	print('select bunch')
	sys.exit()
bunch_select = int(sys.argv[1])


# --- path --- #
path = os.getcwd()

# --- plot --- #
for root_dir, sub_dirs, files in os.walk(path):
    if '==started==' in files or '==completed==' in files:
        if(bunch_select==1): bunch = np.loadtxt(os.path.join(root_dir,'out','integrated_diagnostics','bunch_integrated_quantity_1_dcut.dat'))
        if(bunch_select==2): bunch = np.loadtxt(os.path.join(root_dir,'out','integrated_diagnostics','bunch_integrated_quantity_2_dcut.dat'))

        fig = pyl.figure(1)
        fig.set_size_inches(3.25,3.6,forward=True)
        ax1  = pyl.subplot(111)

        #--- sigma ---#
        ax1.plot(bunch[:,0]/1e4,bunch[:,7], lw=1, label=r"$\sigma_x$ ($\mu$m)")
        ax1.plot(bunch[:,0]/1e4,bunch[:,8], lw=1, label=r"$\sigma_y$ ($\mu$m)")

        #--- emittance ---#
        ax1.plot(bunch[:,0]/1e4,bunch[:,13], lw=1, label=r'$\epsilon_x$ ($\mu$m)' )
        ax1.plot(bunch[:,0]/1e4,bunch[:,14], lw=1, label=r'$\epsilon_y$ ($\mu$m)' )

        #--- energy spread ---#
        ax1.plot(bunch[:,0]/1e4,bunch[:,16]*100, lw=1, label=r'$\sigma_E$%' )

        #--- labels ---#
        ax1.set_xlabel('Z (cm)')


# ax1.legend(loc=9, ncol=2, prop={'size':7.5})
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17),fancybox=True, shadow=False, ncol=3)
pyl.subplots_adjust(bottom=0.30)

# ax1.xaxis.set_major_locator(MultipleLocator(1.))
# ax1.xaxis.set_minor_locator(MultipleLocator(.5))

# pyl.savefig(os.path.join(path,'im_witness_integrated_parameter.pdf'), format='pdf')
plt.show()
