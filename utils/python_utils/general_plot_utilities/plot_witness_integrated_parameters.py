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
### --- ###


# --- path --- #
path = os.getcwd()

# --- plot --- #
for root_dir, sub_dirs, files in os.walk(path):
    if '==started==' in files or '==completed==' in files:
        bunch = np.loadtxt(os.path.join(root_dir,'out','integrated_diagnostics','bunch_integrated_quantity_2_dcut.dat'))

        fig = pyl.figure(1)
        fig.set_size_inches(3.25, 3.0, forward=True)
        ax1  = pyl.subplot(211)

        #--- sigma ---#
        ax1.plot(bunch[:,0]/1e4,bunch[:,8], lw=2, label=r"$\sigma_x$ ($\mu$m)")

        #--- emittance ---#
        ax1.plot(bunch[:,0]/1e4,bunch[:,13], lw=2, label=r'$\epsilon_x$ ($\mu$m)' )

        #--- energy spread ---#
        ax1.plot(bunch[:,0]/1e4,bunch[:,16], lw=2, label=r'$\Delta\gamma/\gamma$' )
        # pyl.ylabel('$\Delta\gamma/\gamma$', fontsize = 8.0)
        # pyl.subplots_adjust(bottom=0.10,left=0.230)
        # plt.legend()
        pyl.xticks(fontsize=9)
        pyl.yticks(fontsize=9)
        # pyl.xlim([0,10])
        # pyl.ylim([0,5])

        #--- gamma plot ---#
        ax2  = pyl.subplot(212)
        ax2.plot(bunch[:,0]/1e4,bunch[:,15]/2., lw=2, label=r'$\Delta\gamma/\gamma$' )
        pyl.xticks(fontsize=9)
        pyl.yticks(fontsize=9)
        # pyl.xlim([0,10])
        pyl.xlabel('$Z$ (cm)', fontsize = 8.0)
        pyl.ylabel('Energy (MeV)', fontsize = 8.0)

pyl.subplots_adjust(bottom=0.14,left=0.210)


ax1.legend(loc=9, ncol=2, prop={'size':7.5})
# ax2.legend(loc=9, ncol=2, prop={'size':7.5})
# ax3.legend(loc=9, ncol=2, prop={'size':7.5})

# ax1.xaxis.set_major_locator(MultipleLocator(1.))
# ax1.xaxis.set_minor_locator(MultipleLocator(.5))

# pyl.savefig('im_witness_integrated_parameter.eps', format='eps')
pyl.savefig(os.path.join(os.path.expanduser('~'),'Downloads','im_witness_integrated_parameter.pdf'), format='pdf')
plt.show()
