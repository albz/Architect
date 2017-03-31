#!/usr/bin/python
######################################################################
# Name:         read_Architect_bin_section.py
# Author:        A. Marocchino
# Date:			  2016-11-10
# Purpose:     General utilities for plots
# Source:       Python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time
import numpy as np
import pylab as pyl
import struct
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes/Code_Architect/Architect/utils/python_utils/architect_graphycal_unit'))
import global_variables as var
### --- ###


def read_output_version(dir_path,file_name):
    path  = os.path.join(os.path.join(dir_path,file_name))
    f        = open(path,'rb')
    var.output_version = struct.unpack('i', f.read(4))[0]

def list_outputs(path):
    for root, dirs, files in os.walk(path):
        var.list_outputs=files
        var.last_output=files[-1]

def read_argv(argv):
    var.frm_number=-1
    var.what = None
    var.path = None
    for arg in argv:
        #---frame number---#
        if var.frm_number == -1:
            try:
                var.frm_number=int(arg)
            except:
                var.frm_number=-1
        #---section or PS---#
        if var.what == None:
            if arg == 'section':
                var.what = 'section'
            if arg == 'PS':
                var.what = 'PS'
        #--- path ---#
        if var.path == None:
            if os.path.exists( str(arg) ) == True:
                var.path = arg

        #--- if left undefined ---#
    if( var.path == None ): var.path = os.getcwd()
    if( var.what == None): var.what = 'section'

    if var.what == 'section':
        var.path = os.path.join(var.path,'out','2D')
    if var.what == 'PS':
        var.path = os.path.join(var.path,'out','PS')
