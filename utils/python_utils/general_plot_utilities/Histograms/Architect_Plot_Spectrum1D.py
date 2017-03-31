#!/usr/bin/python
######################################################################
# Name:         Architect_Plot_Spectrum1D.py
# Author:       F Massimo      
# Date:         2016-02-08
# Purpose:      plots 1D spectrum of Architect Particles data
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


# Usage: ./Architect_Plot_Spectrum1D.py R/W horiz_axis iter nbins/bins_width bins_input
# 
# R: reads pre-existing spectrum file - W: writes new spectrum file, overwrites existing one
# 
# Valid options for axis:"
# X, Y, Z, Px, Py, Pz, Nbunch, cut, dcut
# 
# nbins: bins_input will be number of bins for the horizontal axis"
# 
# bins_width: bins_input will be the widths of each bin on the axis, in the corresponding units 


# ------ Read or Compute Histogram

inputs = sys.argv
bins,horiz_axis,H = Get_Spectrum1D(inputs)

# ----- Plot histogram
plt.plot(bins,H[:,0], linewidth=2.0)
plt.ylim(0. , 1.1*max(H[:,0])  )
plt.xlabel(horiz_axis)
plt.ylabel('dN particles/d'+horiz_axis)
plt.axis('tight') 

plt.show()