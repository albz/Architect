#!/usr/bin/python
######################################################################
# Name:         Architect_Spectrum_utils.py
# Author:       F Massimo      
# Date:         2016-02-08
# Purpose:      Functions for spectrum visualization for Architect
# Source:       Python
#####################################################################



### loading shell commands
import os,sys
import struct
from scipy import *
import numpy as np
import math
# from matplotlib import *
# from pylab import *
#import matplotlib.pyplot as plt
from Architect_utilities import *
from Architect_Load_Phase_Space import *


########## 1D spectrum ##########
def Get_Spectrum1D(inputs):
	if ( (len(inputs)>6) or (len(inputs) < 4) ): # wrong number of arguments
		print 'Error: wrong number of arguments'
		print_usage_1D_spectrum()
		exit(0)

	options_RW  = ['R','W']
	RW          = inputs[1]
	if (not (RW in options_RW)):
		print_usage_1D_spectrum()
		exit(0)

	if ( RW == 'R' ): # Read from histogram file
		
		if (len(inputs)==4): # Read correctly from histogram file
			horiz_axis  = inputs[2]
			iter        = inputs[3]

			bins,horiz_axis,H = Read_Spectrum1D(horiz_axis,iter)
		else:
			print "Error in number of arguments"
			print_usage_1D_spectrum()
			exit(0)
			
	elif ( RW == 'W' ):
		if ( len(inputs)<6  ): # Writing a new histogram file
			print "Error: Not enough input arguments to compute new spectrum"
			print_usage_1D_spectrum()
			exit(0)
	
		elif (len(inputs)==6):
			horiz_axis  = inputs[2]
			iter        = inputs[3]
			bin_choice  = inputs[4]
			bins_input  = inputs[5]
			
			options_bins= ['nbins','bins_width'] 		
			if (not (bin_choice in options_bins)):
				print "Error in binning inputs"
				print_usage_2D_spectrum()
				exit(0)
	
			bins,horiz_axis,H = Write_Spectrum1D(horiz_axis,iter,bin_choice,bins_input)
		
	return bins,horiz_axis,H


def Write_Spectrum1D(horiz_axis,iter,bin_choice,bins_input):  # loads phase space, writes 1d spectrum on file, returns spectrum data
	
	options_bins= ['nbins','bins_width'] 		
	valid_axes  = ['X','Y','Z','Px','Py','Pz','Nbunch','cut','dcut']
	
	if (iter<0):
		print "Error: Invalid iteration number"
		exit(0)
	
	if (not (horiz_axis in valid_axes)):
		print "Error in axis choice"
		print_usage_1D_spectrum()
		exit(0)
	
	X,Y,Z,Px,Py,Pz,Nbunch,cut,dcut = Load_PS(iter) # loads phase space
	
	axis   = {'X':X,'Y':Y,'Z':Z,'Px':Px,'Py':Py,'Pz':Pz,'Nbunch':Nbunch,'cut':cut,'dcut':dcut} 
	
	h_axis = axis[horiz_axis] # use of axis dictionary to make the axis choice more synthetic	
	
	if (bin_choice=='nbins'): # number of bins in histogram are specified
		# check if nbins are positive integers
		nbins = int(bins_input)		
		
		if (nbins<=0): 
			print ' '
			print ' '
			print 'Error: invalid number of bins'
			print_usage_1D_spectrum()
			exit(0)	
			
		# ----- creates bins vectors
		bins 	   = np.linspace(min(h_axis),max(h_axis),nbins+1)
		delta_bin = bins[1]-bins[0]
	
		
	elif (bin_choice=='bins_width'): # bin widths are specified
		print 'Horiz axis minimum: ',min(h_axis),' \t- Horiz axis maximum: ',max(h_axis)
		x_res_for_100_bins = ( max(h_axis) - min(h_axis) ) / 100.
		print 'To have a 100 bin resolution, choose bin widths as: ',x_res_for_100_bins # advice is always welcome
		
		delta_bin = float(bins_input)
		
		bins      = arange(min(h_axis),max(h_axis),delta_bin)
		bins      = bins[0:-1]
		nbins     = len(bins)-1
		if (nbins<2):
			print "Error: bin size too coarse for horizontal axis"
			exit(0)
		
	# ----- computes histogram H
	H          = np.zeros(shape=(nbins,1))
	
	nparticles = len(h_axis)	
	i          = np.minimum( (h_axis-bins[0])/delta_bin , ones(nparticles)*(nbins-1) )
	i          = i.astype(int)
	
	for k in range(0,nparticles):
		H[i[k]][0] += 1.

	# ----- writes histogram on file
	path = os.getcwd()	
	spectrum_directory = generate_folder(path,horiz_axis)	
	file_name = horiz_axis+'_'+str(iter)+'.dat'	
	path_file = os.path.join(spectrum_directory,file_name)
	
	print 'Writing Histogram'

	f=open(path_file,'w')
	f.close() # to overwrite file if exists
  
	f=open(path_file,'a+')
	
	f.write('# '+horiz_axis+'\t dNparticles'+horiz_axis+'\n') # header with quantities
	f.write('# '+str(nbins)+'\n') # header with number of bins 
	data = np.array([bins[0:-1],H[:,0] ])
	data = data.T
	np.savetxt(f, data, fmt=['%f','%f'])
	
	f.close()
	
	return bins[0:-1],horiz_axis,H

def Read_Spectrum1D(horiz_axis,iter): # reads from spectrum file
	
	# checking validity of inputs
	valid_axes  = ['X','Y','Z','Px','Py','Pz','Nbunch','cut','dcut']
	
	if (not (horiz_axis in valid_axes)):
		print "Error in axis choice"
		print_usage_1D_spectrum()
		exit(0)
		
	if (iter<0):
		print "Error: Invalid iteration number"
		exit(0)
	
	# creating spectrum path and filename
	path = os.getcwd()
	directory = os.path.join(path,'out')
	directory = os.path.join(directory,'Histograms')	
		
	if not os.path.exists( directory ):
		os.makedirs(directory)
		
	spectrum_directory = os.path.join(directory,horiz_axis)
	
	if not os.path.exists( directory ):
		os.makedirs(spectrum_directory)
		
	file_name = horiz_axis+'_'+iter+'.dat'	
	path_file = os.path.join(spectrum_directory,file_name)
	print 'Reading spectrum file ',path_file
	
	if (not os.path.isfile(path_file)): # spectrum file non existing
		print ' '
		print ' '
		print 'Error: Spectrum data file not existing.'
		print_usage_1D_spectrum()
		exit(0)
	
	# Preparing for reading...
	bins   = [ ]
	
	f      = open(path_file,'r')  # opens Spectrum file 
	f.readline() 				  # skip header
	line   = f.readline()		  # reads number of bins
	line   = line.strip()
	vars   = line.split()
	nbins  = int(vars[1])
	H      = np.zeros(shape=(nbins,1))
	 
	# Reads spectrum data and puts them in the correct form 
	i = 0

	for line in f:
		
		line = line.strip()
		vars = line.split()
		vars = [float(x) for x in vars]
		bins.append (float(vars[0]))
		H[i][0]   = (float(vars[1]))	
		i += 1	
		
	f.close()
	
	return bins,horiz_axis,H

########### 2D spectrum #############

def Get_Spectrum2D(inputs):
	
	if ( (len(inputs)>8) or (len(inputs) < 5) ): # wrong number of arguments
		print 'Error: wrong number of arguments'
		print_usage_2D_spectrum()
		exit(0)

	options_RW  = ['R','W']
	RW          = inputs[1]
	if (not (RW in options_RW)):
		print_usage_2D_spectrum()
		exit(0)

	if ( RW == 'R' ): # Read from histogram file
		
		if (len(inputs)==5): # Read correctly from histogram file
			horiz_axis  = inputs[2]
			vert_axis   = inputs[3]
			iter        = inputs[4]

			binsx,binsy,horiz_axis,vert_axis,H = Read_Spectrum2D(horiz_axis,vert_axis,iter)
		else:
			print "Error in number of arguments"
			print_usage_2D_spectrum()
			exit(0)
			
	elif ( RW == 'W' ):
		if ( len(inputs)<7  ): # Writing a new histogram file
			print "Error: Not enough input arguments to compute new spectrum"
			print_usage_2D_spectrum()
			exit(0)
	
		elif ( (len(inputs)==7) or (len(inputs)==8) ):
			horiz_axis  = inputs[2]
			vert_axis   = inputs[3]
			iter        = inputs[4]
			bin_choice  = inputs[5]
			binsx_input = inputs[6]
			
			options_bins= ['nbins','bins_width'] 		
			if (not (bin_choice in options_bins)):
				print "Error in binning inputs"
				print_usage_2D_spectrum()
				exit(0)
		
			if (len(inputs)==7):
				binsy_input = binsx_input
				print 'bins_input for vertical axis not specified, set equal to bins_input for horizontal axis: '
				print 'bins_input for horizontal axis = ',binsx_input,'; bins_input for vertical axis = ',binsy_input
			elif (len(inputs)==8):
				binsy_input = inputs[7]
	
			binsx,binsy,horiz_axis,vert_axis,H = Write_Spectrum2D(horiz_axis,vert_axis,iter,bin_choice,binsx_input,binsy_input)
		
	return binsx,binsy,horiz_axis,vert_axis,H
	
	
def Write_Spectrum2D(horiz_axis,vert_axis,iter,bin_choice,binsx_input,binsy_input):  # loads phase space, writes 2d spectrum on file, returns spectrum data

	valid_axes  = ['X','Y','Z','Px','Py','Pz','Nbunch','cut','dcut']
	options_bins= ['nbins','bins_width'] 		

	if (not (horiz_axis in valid_axes)) or (not (vert_axis in valid_axes)):
		print "Error in axis choice"
		print_usage_2D_spectrum()
		exit(0)
	
	if (iter<0):
		print "Error: Invalid iteration number"
		exit(0)

	X,Y,Z,Px,Py,Pz,Nbunch,cut,dcut = Load_PS(iter) # loads phase space

	axis   = {'X':X,'Y':Y,'Z':Z,'Px':Px,'Py':Py,'Pz':Pz,'Nbunch':Nbunch,'cut':cut,'dcut':dcut} 

	h_axis = axis[horiz_axis] # use of dictionaries to make these assignments more synthetic
	v_axis = axis[vert_axis]
	
	if (bin_choice=='nbins'): # number of bins in histogram are specified
		# check if nbins are positive integers
		nbinsx = int(binsx_input)
		nbinsy = int(binsy_input)
		
		if (nbinsx*nbinsy<=0): 
			print ' '
			print ' '
			print 'Error: invalid number of bins'
			print_usage_2D_spectrum()
			exit(0)	
			
		# ----- creates bins vectors
		binsx 	   = np.linspace(min(h_axis),max(h_axis),nbinsx+1)
		delta_binx = binsx[1]-binsx[0]
	
		binsy 	   = np.linspace(min(v_axis),max(v_axis),nbinsy+1)
		delta_biny = binsy[1]-binsy[0]
		
	elif (bin_choice=='bins_width'): # bin widths are specified
		print 'Horiz axis minimum: ',min(h_axis),' \t- Horiz axis maximum: ',max(h_axis)
		print 'Vertl axis minimum: ',min(v_axis),' \t- Vertl axis maximum: ',max(v_axis)
		x_res_for_100_bins = ( max(h_axis) - min(h_axis) ) / 100.
		y_res_for_100_bins = ( max(v_axis) - min(v_axis) ) / 100.
		print 'To have a 100 bin resolution in both axes, choose bin widths as: ',x_res_for_100_bins,' ',y_res_for_100_bins # advice is always welcome
		
		delta_binx = float(binsx_input)
		delta_biny = float(binsy_input)
		
		binsx      = arange(min(h_axis),max(h_axis),delta_binx)
		binsx      = binsx[0:-1]
		nbinsx     = len(binsx)-1
		if (nbinsx<2):
			print "Error: bin size too coarse for horizontal axis"
			exit(0)
		
		binsy      = arange(min(v_axis),max(v_axis),delta_biny)
		binsy      = binsy[0:-1]
		nbinsy     = len(binsy)-1
		if (nbinsy<2):
			print "Error: bin size too coarse for vertical axis"
			exit(0)
		
	# ----- computes histogram H
	H          = np.zeros(shape=(nbinsx,nbinsy))
	
	nparticles = len(h_axis)	
	i          = np.minimum( (h_axis-binsx[0])/delta_binx , ones(nparticles)*(nbinsx-1) )
	i          = i.astype(int)
	j          = np.minimum( (v_axis-binsy[0])/delta_biny , ones(nparticles)*(nbinsy-1) )
	j          = j.astype(int)
	
	for k in range(0,nparticles):
		H[i[k]][j[k]] += 1.

	# ----- writes histogram on file
	path = os.getcwd()	
	spectrum_directory = generate_folder(path,horiz_axis+'_'+vert_axis)
		
	file_name = horiz_axis+'_'+vert_axis+'_'+str(iter)+'.dat'	
	path_file = os.path.join(spectrum_directory,file_name)

	print 'Writing Histogram'

	f=open(path_file,'w')
	f.close() # to overwrite file if exists
  
	f=open(path_file,'a+')
	f.write('# '+horiz_axis+'\t'+vert_axis+'\t dNparticles / d'+horiz_axis+'d'+vert_axis+'\n') # header with quantities
	f.write('# '+str(nbinsx)+'\t'+str(nbinsy)+'\n') # header with number of bins 
	
	
	for j in range(0,nbinsy):
		data = np.array([binsx[0:-1],ones(nbinsx)*binsy[j],H[:,j] ])
		data = data.T
		np.savetxt(f, data, fmt=['%f','%f','%f'])
		f.write('\n')
	f.close()
	
	return binsx[0:-1],binsy[0:-1],horiz_axis,vert_axis,H


def Read_Spectrum2D(horiz_axis,vert_axis,iter): # reads from spectrum file
	
	# checking validity of inputs
	valid_axes  = ['X','Y','Z','Px','Py','Pz','Nbunch','cut','dcut']
	
	if (not (horiz_axis in valid_axes)) or (not (vert_axis in valid_axes)):
		print "Error in axis choice"
		print_usage_2D_spectrum()
		exit(0)
		
	if (iter<0):
		print "Error: Invalid iteration number"
		exit(0)
	
	# creating spectrum path and filename
	path = os.getcwd()
	directory = os.path.join(path,'out')
	directory = os.path.join(directory,'Histograms')	
		
	if not os.path.exists( directory ):
		os.makedirs(directory)
		
	spectrum_directory = os.path.join(directory,horiz_axis+'_'+vert_axis)
	
	if not os.path.exists( directory ):
		os.makedirs(spectrum_directory)
		
	file_name = horiz_axis+'_'+vert_axis+'_'+iter+'.dat'	
	path_file = os.path.join(spectrum_directory,file_name)
	print 'Reading spectrum file ',path_file
	
	if (not os.path.isfile(path_file)): # spectrum file non existing
		print ' '
		print ' '
		print 'Error: Spectrum data file not existing.'
		print_usage_2D_spectrum()
		exit(0)
	
	# Preparing for reading...
	binsx  = [ ]
	binsy  = [ ]
	
	f      = open(path_file,'r')  # opens Spectrum file 
	f.readline() 				  # skip header
	line   = f.readline()		  # reads number of bins
	line   = line.strip()
	vars   = line.split()
	nbinsx = int(vars[1])
	nbinsy = int(vars[2])
	H      = np.zeros(shape=(nbinsx,nbinsy))
	 
	# Reads spectrum data and puts them in the correct form 
	i = 0
	j = 0

	for line in f:
			
		line = line.strip()
		vars = line.split()
		vars = [float(x) for x in vars]
		
		H[i][j]   = (float(vars[2]))
		
		if ((j==0)&(i<nbinsx)):
			binsx.append(float(vars[0]))
			
		i += 1
			
		if (i%nbinsx==0):
			binsy.append(float(vars[1]))
			i  = 0
			j += 1
		if (i==0):
			f.next() # skip empty line		
	f.close()

	return binsx,binsy,horiz_axis,vert_axis,H
	
	
#################################################
def print_usage_1D_spectrum():
	print ' '
	print ' '
	print "Usage: ./Architect_Plot_Spectrum1D.py R/W horiz_axis iter nbins/bins_width bins_input "
	print ' '
	print "R: reads pre-existing spectrum file - W: writes new spectrum file, overwrites existing one"
	print ' '
	print "Valid options for axis:"
	print "X, Y, Z, Px, Py, Pz, Nbunch, cut, dcut"
	print ' '
	print "nbins: bins_input will be number of bins for the horizontal axis"
	print ' '
	print "bins_width: bins_input will be the widths of each bin on the axis, in the corresponding units "



def print_usage_2D_spectrum():
	print ' '
	print ' '
	print "Usage: ./Architect_Plot_Spectrum2D.py R/W horiz_axis vert_axis iter nbins/bins_width binsx_input binsy_input"
	print ' '
	print "R: reads pre-existing spectrum file - W: writes new spectrum file, overwrites existing one"
	print ' '
	print "Valid options for axes:"
	print "X, Y, Z, Px, Py, Pz, Nbunch, cut, dcut"
	print ' '
	print "nbins: bins_inputx and bins_inputy will be number of bins for the axes"
	print ' '
	print "bins_width: bins_inputx and bins_inputy will be the widths of each bin on the axes, in their corresponding units "





