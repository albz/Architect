#!/usr/bin/python
######################################################################
# Name:         read_Architect_bin_section.py
# Author:       F. Massimo
# Date:			2016-02-09
# Purpose:      reads Architect section binary from output
# Source:       Python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time
import struct
#from scipy import *
import numpy as np
# from matplotlib import *
# from pylab import *
# import matplotlib as plt
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
### --- ###



# Advice: always call these routines with
# f              = open(file_name,'rb')    # filename is the name of the binary 2D file
# output_version = struct.unpack('i', f.read(4))[0]
# f.close()
# 
# if output_version == 1:
# 	dist,r_mesh,z_mesh,rho_b,n_bck,Er,Ez,Bphi,Jbr,Jbckr,Jbz,Jbckz                                                       = read_Architect_bin_section_output_v_1(dir_path,file_name)
# elif output_version == 2:
# 	dist,r_mesh,z_mesh,rho_b,n_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,Jbr,Jbckr,Jbz,Jbckz               = read_Architect_bin_section_output_v_2(dir_path,file_name)
# elif output_version == 3:
# 	dist,r_mesh,z_mesh,rho_b,n_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,B_ex_poloidal,Jbr,Jbckr,Jbz,Jbckz = read_Architect_bin_section_output_v_3(dir_path,file_name)
		






def read_Architect_bin_section_output_v_1(dir_path,file_name):

	# - #
	path     = os.path.join(os.path.join(dir_path,file_name))
	f        = open(path,'rb')

	output_version = struct.unpack('i', f.read(4))[0]
	dist           = struct.unpack('i', f.read(4))[0] # distance traveled
	#- radial and axial dimensions -#
	Nr             = struct.unpack('i', f.read(4))[0] # half-plane axis dim (radial direction      )
	Nz             = struct.unpack('i', f.read(4))[0] # z dim 				(longitudinal direction)

	
	#### Start reading grid data
	
	# r_mesh axis
	r_mesh = np.zeros(Nr)
	for j in range(0,Nr): # read r axis
		r = struct.unpack('f', f.read(4))[0]
		r_mesh[j] = r
		
	# z_mesh axis
	z_mesh = np.zeros(Nz)
	for i in range(0,Nz): # read z axis
		z = struct.unpack('f', f.read(4))[0]
		z_mesh[i] = z
	
	# rho_b  -  bunch electron density			
	rho_b = np.zeros((Nr,Nz))	
	for j in range(0,Nr):
		for i in range(0,Nz):
			rho_b[j,i] = struct.unpack('f', f.read(4))[0] 
 	
	# n_bck - plasma electron density
	n_bck = np.zeros((Nr,Nz))	
	for j in range(0,Nr):
		for i in range(0,Nz):
			n_bck[j,i] = struct.unpack('f', f.read(4))[0] 
	
	# Er - Electric field, transverse
	Er = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Er[j,i] = struct.unpack('f', f.read(4))[0] 

	# Ez - Electric field, longitudinal
	Ez = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Ez[j,i] = struct.unpack('f', f.read(4))[0] 

	# Bphi - Magnetic field, azimuthal
	Bphi = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Bphi[j,i] = struct.unpack('f', f.read(4))[0] 

	# Jbr - Bunch current density, transverse
	Jbr = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbr[j,i] = struct.unpack('f', f.read(4))[0]
			 
	# Jbckr - Plasma current density, transverse
	Jbckr = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbckr[j,i] = struct.unpack('f', f.read(4))[0]
	
	# Jbr - Bunch current density, transverse
	Jbz = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbz[j,i] = struct.unpack('f', f.read(4))[0] 
	# Jbckr - Plasma current density, transverse
	Jbckz = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbckz[j,i] = struct.unpack('f', f.read(4))[0]  
 		
	f.close()
	
	return dist,r_mesh,z_mesh,rho_b,n_bck,Er,Ez,Bphi,Jbr,Jbckr,Jbz,Jbckz
	
	
def read_Architect_bin_section_output_v_2(dir_path,file_name):

	# - #
	path     = os.path.join(os.path.join(dir_path,file_name))
	f        = open(path,'rb')

	output_version = struct.unpack('i', f.read(4))[0]
	dist           = struct.unpack('i', f.read(4))[0] # distance traveled
	#- vector length -#
	Nr             = struct.unpack('i', f.read(4))[0] # half-plane axis dim (radial direction      )
	Nz             = struct.unpack('i', f.read(4))[0] # z dim 				(longitudinal direction)
		
	#### Start reading grid data
	
	# r_mesh axis
	r_mesh = np.zeros(Nr)
	for j in range(0,Nr): # read r axis
		r = struct.unpack('f', f.read(4))[0]
		r_mesh[j] = r
		
	# z_mesh axis
	z_mesh = np.zeros(Nz)
	for i in range(0,Nz): # read z axis
		z = struct.unpack('f', f.read(4))[0]
		z_mesh[i] = z
	
	# rho_b  -  bunch electron density			
	rho_b = np.zeros((Nr,Nz))	
	for j in range(0,Nr):
		for i in range(0,Nz):
			rho_b[j,i] = struct.unpack('f', f.read(4))[0] 
 	
	# n_bck - plasma electron density
	n_bck = np.zeros((Nr,Nz))	
	for j in range(0,Nr):
		for i in range(0,Nz):
			n_bck[j,i] = struct.unpack('f', f.read(4))[0] 
	
	# Er - Electric field, transverse
	Er = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Er[j,i] = struct.unpack('f', f.read(4))[0] 

	Er_bck = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Er_bck[j,i] = struct.unpack('f', f.read(4))[0] 

	Er_b = np.zeros((Nr,Nz))	# bunch field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Er_b[j,i] = struct.unpack('f', f.read(4))[0] 

	# Ez - Electric field, longitudinal
	Ez = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Ez[j,i] = struct.unpack('f', f.read(4))[0] 

	Ez_bck = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Ez_bck[j,i] = struct.unpack('f', f.read(4))[0] 

	Ez_b = np.zeros((Nr,Nz))	# bunch field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Ez_b[j,i] = struct.unpack('f', f.read(4))[0] 

	# Bphi - Magnetic field, azimuthal
	Bphi = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Bphi[j,i] = struct.unpack('f', f.read(4))[0] 

	Bphi_bck = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Bphi_bck[j,i] = struct.unpack('f', f.read(4))[0] 

	Bphi_b = np.zeros((Nr,Nz))	# bunch field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Bphi_b[j,i] = struct.unpack('f', f.read(4))[0] 

	# Jbr - Bunch current density, transverse
	Jbr = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbr[j,i] = struct.unpack('f', f.read(4))[0] 
			
	# Jbckr - Plasma current density, transverse
	Jbckr = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbckr[j,i] = struct.unpack('f', f.read(4))[0]
	
	# Jbr - Bunch current density, transverse
	Jbz = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbz[j,i] = struct.unpack('f', f.read(4))[0] 
			
	# Jbckr - Plasma current density, transverse
	Jbckz = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbckz[j,i] = struct.unpack('f', f.read(4))[0]  


 	f.close()
 	
 	return dist,r_mesh,z_mesh,rho_b,n_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,Jbr,Jbckr,Jbz,Jbckz 		
	
	
	
	
def read_Architect_bin_section_output_v_3(dir_path,file_name):

	# - #
	path     = os.path.join(os.path.join(dir_path,file_name))
	f        = open(path,'rb')

	output_version = struct.unpack('i', f.read(4))[0]
	dist           = struct.unpack('i', f.read(4))[0] # distance traveled
	#- vector length -#
	Nr             = struct.unpack('i', f.read(4))[0] # half-plane axis dim (radial direction      )
	Nz             = struct.unpack('i', f.read(4))[0] # z dim 				(longitudinal direction)
	
	#### Start reading grid data
	
	# r_mesh axis
	r_mesh = np.zeros(Nr)
	for j in range(0,Nr): # read r axis
		r = struct.unpack('f', f.read(4))[0]
		r_mesh[j] = r
		
	# z_mesh axis
	z_mesh = np.zeros(Nz)
	for i in range(0,Nz): # read z axis
		z = struct.unpack('f', f.read(4))[0]
		z_mesh[i] = z
	
	# rho_b  -  bunch electron density			
	rho_b = np.zeros((Nr,Nz))	
	for j in range(0,Nr):
		for i in range(0,Nz):
			rho_b[j,i] = struct.unpack('f', f.read(4))[0] 
 	
	# n_bck - plasma electron density
	n_bck = np.zeros((Nr,Nz))	
	for j in range(0,Nr):
		for i in range(0,Nz):
			n_bck[j,i] = struct.unpack('f', f.read(4))[0] 
	
	# Er - Electric field, transverse
	Er = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Er[j,i] = struct.unpack('f', f.read(4))[0] 

	Er_bck = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Er_bck[j,i] = struct.unpack('f', f.read(4))[0] 

	Er_b = np.zeros((Nr,Nz))	# bunch field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Er_b[j,i] = struct.unpack('f', f.read(4))[0] 

	# Ez - Electric field, longitudinal
	Ez = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Ez[j,i] = struct.unpack('f', f.read(4))[0] 

	Ez_bck = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Ez_bck[j,i] = struct.unpack('f', f.read(4))[0] 

	Ez_b = np.zeros((Nr,Nz))	# bunch field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Ez_b[j,i] = struct.unpack('f', f.read(4))[0] 

	# Bphi - Magnetic field, azimuthal
	Bphi = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Bphi[j,i] = struct.unpack('f', f.read(4))[0] 

	Bphi_bck = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Bphi_bck[j,i] = struct.unpack('f', f.read(4))[0] 

	Bphi_b = np.zeros((Nr,Nz))	# bunch field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Bphi_b[j,i] = struct.unpack('f', f.read(4))[0] 

	B_ex_poloidal = np.zeros((Nr,Nz))	# bunch field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			B_ex_poloidal[j,i] = struct.unpack('f', f.read(4))[0]

	# Jbr - Bunch current density, transverse
	Jbr = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbr[j,i] = struct.unpack('f', f.read(4))[0] 
			
	# Jbckr - Plasma current density, transverse
	Jbckr = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbckr[j,i] = struct.unpack('f', f.read(4))[0]
	
	# Jbr - Bunch current density, transverse
	Jbz = np.zeros((Nr,Nz))      # total field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbz[j,i] = struct.unpack('f', f.read(4))[0] 
			
	# Jbckr - Plasma current density, transverse
	Jbckz = np.zeros((Nr,Nz))	# background field	
	for j in range(0,Nr):
		for i in range(0,Nz):
			Jbckz[j,i] = struct.unpack('f', f.read(4))[0]  


 	f.close()
 	
 	return dist,r_mesh,z_mesh,rho_b,n_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,B_ex_poloidal,Jbr,Jbckr,Jbz,Jbckz 	
	 
 
 		
	
	
	
