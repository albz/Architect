#!/usr/bin/python
######################################################################
# Name:         Architect_Plot_Spectrum2D.py
# Author:       F Massimo      
# Date:         2016-02-08
# Purpose:      plots 2D spectrum of Architect Particles data
# Source:       Python
#####################################################################



### loading shell commands
import os,sys
import struct
from scipy import *
import numpy as np
# from matplotlib import *
# from pylab import *
#import matplotlib.pyplot as plt
from Architect_utilities import *
from Architect_Spectrum_utils import *
import matplotlib.pyplot as plt

# Usage: ./Architect_Plot_Spectrum2D.py R/W horiz_axis vert_axis iter nbins/bins_width binsx_input binsy_input
#
# R: reads pre-existing spectrum file - W: writes new spectrum file, overwrites existing one
#
# Valid options for axes:
# X, Y, Z, Px, Py, Pz, Nbunch, cut, dcut
# 
# nbins: bins_inputx and bins_inputy will be number of bins for the axes
# 
# bins_width: bins_inputx and bins_inputy will be the widths of each bin on the axes, in their corresponding units 


# ------ Read or Compute Histogram

inputs = sys.argv
binsx,binsy,horiz_axis,vert_axis,H = Get_Spectrum2D(inputs)

# ----- Plot histogram
cnt = plt.contourf(binsx,binsy,H.T,70 )
plt.xlabel(horiz_axis)
plt.ylabel(vert_axis)
plt.title('dNparticles /d'+horiz_axis+'d'+vert_axis)
plt.axis('tight')
plt.colorbar(ticks=np.linspace(0.,H[:][:].max(),5), format='%0.2f')
for c in cnt.collections: # this command hides white contour lines
    c.set_edgecolor("face") 

plt.show()
