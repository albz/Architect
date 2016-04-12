#!/usr/bin/python
######################################################################
# Name:         Architect_Load_Phase_Space.py
# Author:       F Massimo
# Date:         2016-02-08
# Purpose:      reads Architect Particles data in .arch file
# Source:       Python
#####################################################################


### loading shell commands
import os,sys
import struct
from scipy import *
import numpy as np
# from matplotlib import *
# from pylab import *
import matplotlib.pyplot as plt
from Architect_utilities import *

### --- ###

def Load_PS(iter):

	#-path
	path = os.getcwd()
	
	X                = [ ]
	Y                = [ ]
	Z                = [ ]
	Px               = [ ]
	Py               = [ ]
	Pz               = [ ]
	Nbunch           = [ ]
	cut              = [ ]
	dcut             = [ ]
	

	# Read phase space
	directory_PS     = os.path.join(path,'out')
	directory_PS     = os.path.join(directory_PS,'PS')
	file_name        = str(iter)+'.arch'
	path_file        = os.path.join(directory_PS,file_name)
	print 'Loading Phase Space from ',path_file
	if (not os.path.isfile(path_file)):
		print 'Error: non-existing Phase Space file.'
		exit(0)

	if (os.stat(path_file).st_size == 0): # checks if file is not empty
		print 'Empty Phase Space file'
		exit(0)

	f  = open(path_file,'r') 			# opens Phase Space file
	
	# read header lines
	Output_version = struct.unpack('i', f.read(4))
	avgz           = struct.unpack('i', f.read(4)) # distance traveled along z by beam center
	Np             = struct.unpack('i', f.read(4))
		
	vars=[]

 	for i in range(0,Np[0]):
 		vars=struct.unpack('fffffffff', f.read(4*9))
 		
 		X.append(vars[0])
 		Y.append(vars[1])
 		Z.append(vars[2])
 		Px.append(vars[3])
 		Py.append(vars[4])
 		Pz.append(vars[5])
 		Nbunch.append(vars[6])
 		cut.append(vars[7])
 		dcut.append(vars[7])
 		
	f.close()


	return X,Y,Z,Px,Py,Pz,Nbunch,cut,dcut
