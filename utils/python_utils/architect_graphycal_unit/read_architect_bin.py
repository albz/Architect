#!/usr/bin/python
######################################################################
# Name:         read_Architect_bin_section.py
# Author:        A. Marocchino
# Date:			  2016-11-10
# Purpose:     Script Manger for binary file reader of Architect
# Source:       Python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time
import numpy as np
import pylab as pyl
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes/Code_Architect/Architect/utils/python_utils/architect_graphycal_unit'))
import global_variables as var
from architect_read_section_bin import *
from architect_read_PS_bin import *
from general_utilities import *
### --- ###


#--- *** ---# shell inputs
def read_architect_bin(*argv):

    read_argv(argv)
    list_outputs(var.path)
    read_output_version(var.path,var.list_outputs[var.frm_number])
    print('output version ::',var.output_version)
    if(var.output_version==4 and var.what == 'section'): read_Architect_bin_section_output_v_4(var.path,var.list_outputs[var.frm_number])
    if(var.output_version==2 and var.what == 'PS'): architect_read_PS_bin_v_2(var.path,var.list_outputs[var.frm_number])
    if(var.output_version==3 and var.what == 'PS'): architect_read_PS_bin_v_3(var.path,var.list_outputs[var.frm_number])
