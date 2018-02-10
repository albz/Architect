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
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes/Code_Architect/Architect/utils/python_utils/architect_graphycal_unit'))
import global_variables as var
### --- ###

def architect_read_PS_bin(dir_path,file_name):

	# - #
	path     = os.path.join(os.path.join(dir_path,'out','PS',file_name))
	f        = open(path,'rb')

	#-#
	x=[]; y=[]; z=[]; px=[]; py=[]; pz=[]; bunch_id=[]; cut=[]; dcut=[];

	output_version     = struct.unpack('i', f.read(4))[0]
	distance_travelled = struct.unpack('i', f.read(4))[0]
	number_of_particles= struct.unpack('i', f.read(4))[0]

	print('output version of the file >',output_version)
	print('number of particle >',number_of_particles)

	for p in range(0,number_of_particles):
		x.append(struct.unpack('d', f.read(8))[0])
		y.append(struct.unpack('d', f.read(8))[0])
		z.append(struct.unpack('d', f.read(8))[0])
		px.append(struct.unpack('d', f.read(8))[0])
		py.append(struct.unpack('d', f.read(8))[0])
		pz.append(struct.unpack('d', f.read(8))[0])
		bunch_id.append(struct.unpack('d', f.read(8))[0])
		cut.append(struct.unpack('d', f.read(8))[0])
		dcut.append(struct.unpack('d', f.read(8))[0])

	f.close()

	return x,y,z,px,py,pz,bunch_id,cut,dcut



def architect_read_PS_bin_v_2(dir_path,file_name):
	# - #
	path     = os.path.join(os.path.join(dir_path,file_name))
	f        = open(path,'rb')

	#-#
	var.output_version     = struct.unpack('i', f.read(4))[0]
	var.n_bunches           = struct.unpack('i', f.read(4))[0]
	for i in range(0,var.n_bunches): var.bunch_charges.append(struct.unpack('d', f.read(8))[0])
	var.distance_travelled = struct.unpack('i', f.read(4))[0]
	var.number_of_particles= struct.unpack('i', f.read(4))[0]

	print('output version of the file >',var.output_version)
	print('number of bunches >',var.n_bunches)
	print('number of particle >',var.number_of_particles)

	var.x=[]; var.y=[]; var.z=[];
	var.px=[]; var.py=[]; var.pz=[];
	var.bunch_id=[]; var.cut=[]; var.dcut=[];

	for p in range(0,var.number_of_particles):
		var.x.append(struct.unpack('d', f.read(8))[0])
		var.y.append(struct.unpack('d', f.read(8))[0])
		var.z.append(struct.unpack('d', f.read(8))[0])
		var.px.append(struct.unpack('d', f.read(8))[0])
		var.py.append(struct.unpack('d', f.read(8))[0])
		var.pz.append(struct.unpack('d', f.read(8))[0])
		var.bunch_id.append(struct.unpack('d', f.read(8))[0])
		var.cut.append(struct.unpack('d', f.read(8))[0])
		var.dcut.append(struct.unpack('d', f.read(8))[0])

	var.x=np.array(var.x); var.y=np.array(var.y); var.z=np.array(var.z);
	var.px=np.array(var.px); var.py=np.array(var.py); var.pz=np.array(var.pz);
	var.bunch_id=np.array(var.bunch_id)
	f.close()


def architect_read_PS_bin_v_3(dir_path,file_name):
	# - #
	path     = os.path.join(os.path.join(dir_path,file_name))
	f        = open(path,'rb')

	#-#
	var.output_version     = struct.unpack('i', f.read(4))[0]
	var.n_bunches           = struct.unpack('i', f.read(4))[0]
	for i in range(0,var.n_bunches): var.bunch_charges.append(struct.unpack('d', f.read(8))[0])
	var.distance_travelled = struct.unpack('i', f.read(4))[0]
	var.number_of_particles= struct.unpack('i', f.read(4))[0]

	print('output version of the file >',var.output_version)
	print('number of bunches >',var.n_bunches)
	print('number of particle >',var.number_of_particles)
	print('bunch charge(s) >',var.bunch_charges)

	var.x=[]; var.y=[]; var.z=[];
	var.px=[]; var.py=[]; var.pz=[];
	var.bunch_id=[]; var.cut=[]; var.dcut=[];
	var.macro_particle_charge=[]; var.macro_particle_nume=[];

	for p in range(0,var.number_of_particles):
		var.x.append(struct.unpack('d', f.read(8))[0])
		var.y.append(struct.unpack('d', f.read(8))[0])
		var.z.append(struct.unpack('d', f.read(8))[0])
		var.px.append(struct.unpack('d', f.read(8))[0])
		var.py.append(struct.unpack('d', f.read(8))[0])
		var.pz.append(struct.unpack('d', f.read(8))[0])
		var.bunch_id.append(struct.unpack('d', f.read(8))[0])
		var.cut.append(struct.unpack('d', f.read(8))[0])
		var.dcut.append(struct.unpack('d', f.read(8))[0])
		var.macro_particle_charge.append(struct.unpack('d', f.read(8))[0])
		var.macro_particle_nume.append(struct.unpack('d', f.read(8))[0])

	var.x=np.array(var.x); var.y=np.array(var.y); var.z=np.array(var.z);
	var.px=np.array(var.px); var.py=np.array(var.py); var.pz=np.array(var.pz);
	var.bunch_id=np.array(var.bunch_id)
	var.macro_particle_charge=np.array(var.macro_particle_charge)
	var.macro_particle_nume=np.array(var.macro_particle_nume)
	f.close()
