#!/usr/bin/python
######################################################################
# Name:         read_Architect_bin_section.py
# Author:        A. Marocchino (last version: I have made the mess with double precision) and F. Massimo
# Date:			  2016-11-10
# Purpose:     read Binary Sections
# Source:       Python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time
import struct
import numpy as np
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes/Code_Architect/Architect/utils/python_utils/architect_graphycal_unit'))
import global_variables as var
### --- ###


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



def read_Architect_bin_section_output_v_4(dir_path,file_name):
	# - #
	path  = os.path.join(os.path.join(dir_path,file_name))
	f        = open(path,'rb')

	var.output_version = struct.unpack('i', f.read(4))[0]
	var.dist_um = struct.unpack('i', f.read(4))[0] # distance traveled
	#- vector length -#
	var.Nr = struct.unpack('i', f.read(4))[0] # half-plane axis dim (radial direction      )
	var.Nz = struct.unpack('i', f.read(4))[0] # z dim 				(longitudinal direction)

	#--- Grid data ---#

	# r_mesh axis
	var.r_mesh = np.zeros(var.Nr)
	for j in range(0,var.Nr): # read r axis
		r = struct.unpack('d', f.read(8))[0]
		var.r_mesh[j] = r

	# z_mesh axis
	var.z_mesh = np.zeros(var.Nz)
	for i in range(0,var.Nz): # read z axis
		z = struct.unpack('d', f.read(8))[0]
		var.z_mesh[i] = z

	# rho_b  -  bunch electron density
	var.rho_b = np.zeros((var.Nr,var.Nz))
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.rho_b[j,i] = struct.unpack('d', f.read(8))[0]

	# rho_bck - plasma electron density
	var.rho_bck = np.zeros((var.Nr,var.Nz))
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.rho_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.rho = np.zeros((var.Nr,var.Nz))
	var.rho = var.rho_b+var.rho_bck

	# Er fields
	var.Er = np.zeros((var.Nr,var.Nz))
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Er[j,i] = struct.unpack('d', f.read(8))[0]

	var.Er_bck = np.zeros((var.Nr,var.Nz))
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Er_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.Er_b = np.zeros((var.Nr,var.Nz))
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Er_b[j,i] = struct.unpack('d', f.read(8))[0]

	# Ez - Electric field, longitudinal
	var.Ez = np.zeros((var.Nr,var.Nz))      # total field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Ez[j,i] = struct.unpack('d', f.read(8))[0]

	var.Ez_bck = np.zeros((var.Nr,var.Nz))	# background field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Ez_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.Ez_b = np.zeros((var.Nr,var.Nz))	# bunch field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Ez_b[j,i] = struct.unpack('d', f.read(8))[0]

	# Bphi - Magnetic field, azimuthal
	var.Bphi = np.zeros((var.Nr,var.Nz))      # total field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Bphi[j,i] = struct.unpack('d', f.read(8))[0]

	var.Bphi_bck = np.zeros((var.Nr,var.Nz))	# background field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Bphi_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.Bphi_b = np.zeros((var.Nr,var.Nz))	# bunch field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Bphi_b[j,i] = struct.unpack('d', f.read(8))[0]

	var.Bphi_ex = np.zeros((var.Nr,var.Nz))	# bunch field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Bphi_ex[j,i] = struct.unpack('d', f.read(8))[0]

	# Radial Current Jr
	var.Jr_b = np.zeros((var.Nr,var.Nz))      # total field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Jr_b[j,i] = struct.unpack('d', f.read(8))[0]

	#Plasma current density, transverse
	var.Jr_bck = np.zeros((var.Nr,var.Nz))	# background field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Jr_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.Jr = np.zeros((var.Nr,var.Nz))
	var.Jr = var.Jr_b+var.Jr_bck

	# Jbr - Bunch current density, transverse
	var.Jz_b = np.zeros((var.Nr,var.Nz))      # total field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Jz_b[j,i] = struct.unpack('d', f.read(8))[0]

	# Jbckr - Plasma current density, transverse
	var.Jz_bck = np.zeros((var.Nr,var.Nz))	# background field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Jz_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.Jz = np.zeros((var.Nr,var.Nz))
	var.Jz = var.Jz_b+var.Jz_bck

	f.close()


def read_Architect_bin_section_output_v_5(dir_path,file_name):
	# - #
	path  = os.path.join(os.path.join(dir_path,file_name))
	f        = open(path,'rb')

	var.output_version = struct.unpack('i', f.read(4))[0]

	#--- parameters ---#
	var.kp      = struct.unpack('d', f.read(8))[0]
	var.wp      = struct.unpack('d', f.read(8))[0]
	var.dist    = struct.unpack('d', f.read(8))[0]
	var.np      = struct.unpack('d', f.read(8))[0]
	var.dist_um = struct.unpack('d', f.read(8))[0]

	#- vector length -#
	var.Nr = struct.unpack('i', f.read(4))[0] # half-plane axis dim (radial direction      )
	var.Nz = struct.unpack('i', f.read(4))[0] # z dim 				(longitudinal direction)

	#--- Grid data ---#

	# r_mesh axis
	var.r_mesh = np.zeros(var.Nr)
	for j in range(0,var.Nr): # read r axis
		r = struct.unpack('d', f.read(8))[0]
		var.r_mesh[j] = r

	# z_mesh axis
	var.z_mesh = np.zeros(var.Nz)
	for i in range(0,var.Nz): # read z axis
		z = struct.unpack('d', f.read(8))[0]
		var.z_mesh[i] = z

	# rho_b  -  bunch electron density
	var.rho_b = np.zeros((var.Nr,var.Nz))
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.rho_b[j,i] = struct.unpack('d', f.read(8))[0]

	# rho_bck - plasma electron density
	var.rho_bck = np.zeros((var.Nr,var.Nz))
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.rho_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.rho = np.zeros((var.Nr,var.Nz))
	var.rho = var.rho_b+var.rho_bck

	# Er fields
	var.Er = np.zeros((var.Nr,var.Nz))
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Er[j,i] = struct.unpack('d', f.read(8))[0]

	var.Er_bck = np.zeros((var.Nr,var.Nz))
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Er_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.Er_b = np.zeros((var.Nr,var.Nz))
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Er_b[j,i] = struct.unpack('d', f.read(8))[0]

	# Ez - Electric field, longitudinal
	var.Ez = np.zeros((var.Nr,var.Nz))      # total field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Ez[j,i] = struct.unpack('d', f.read(8))[0]

	var.Ez_bck = np.zeros((var.Nr,var.Nz))	# background field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Ez_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.Ez_b = np.zeros((var.Nr,var.Nz))	# bunch field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Ez_b[j,i] = struct.unpack('d', f.read(8))[0]

	# Bphi - Magnetic field, azimuthal
	var.Bphi = np.zeros((var.Nr,var.Nz))      # total field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Bphi[j,i] = struct.unpack('d', f.read(8))[0]

	var.Bphi_bck = np.zeros((var.Nr,var.Nz))	# background field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Bphi_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.Bphi_b = np.zeros((var.Nr,var.Nz))	# bunch field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Bphi_b[j,i] = struct.unpack('d', f.read(8))[0]

	var.Bphi_ex = np.zeros((var.Nr,var.Nz))	# bunch field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Bphi_ex[j,i] = struct.unpack('d', f.read(8))[0]

	# Radial Current Jr
	var.Jr_b = np.zeros((var.Nr,var.Nz))      # total field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Jr_b[j,i] = struct.unpack('d', f.read(8))[0]

	#Plasma current density, transverse
	var.Jr_bck = np.zeros((var.Nr,var.Nz))	# background field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Jr_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.Jr = np.zeros((var.Nr,var.Nz))
	var.Jr = var.Jr_b+var.Jr_bck

	# Jbr - Bunch current density, transverse
	var.Jz_b = np.zeros((var.Nr,var.Nz))      # total field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Jz_b[j,i] = struct.unpack('d', f.read(8))[0]

	# Jbckr - Plasma current density, transverse
	var.Jz_bck = np.zeros((var.Nr,var.Nz))	# background field
	for j in range(0,var.Nr):
		for i in range(0,var.Nz):
			var.Jz_bck[j,i] = struct.unpack('d', f.read(8))[0]

	var.Jz = np.zeros((var.Nr,var.Nz))
	var.Jz = var.Jz_b+var.Jz_bck

	f.close()
