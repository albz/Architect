#!/usr/bin/python
######################################################################
# Name:         architect_read_PS_bin.py
# Author:       A. Marocchino
# Date:			2016-07-25
# Purpose:      read Architect PS
# Source:       Python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time
import struct
import numpy as np
sys.path.append(os.path.join(os.path.expanduser('~'),'/Users/albz/Codes/Code_Architect/Architect/utils/python_utils'))
from architect_read_PS_bin import *
### --- ###

dir_path  = '/Users/albz/sims/sims_PWFA/random_test'
file_name = '0000000.arch'

x,y,z,px,py,pz,bunch_id,cut,dcut = architect_read_PS_bin(dir_path,file_name)

np.savetxt('PS_ascii.arch',np.column_stack((x,y,z,px,py,pz,bunch_id)), delimiter='\t')
