#!/usr/bin/python
######################################################################
# Name:         read_Architect_bin_section.py
# Author:        A. Marocchino
# Date:			 2016-11-10
# Purpose:    Defines the variable for the plots in a Global style
# Source:       Python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time
import matplotlib.pyplot as plt
import numpy as np
### --- ###

### --- Define Global Variables ---###


#--- **************** ---#
#--- Plot variables ---#
#--- *************** ---#
path=''
what=''
frm_number=0
list_outputs=[]
last_output=[]


#--- ***************************************** ---#
#--- ARCHITECT OUTPUT VARIABLES ---#
#--- ***************************************** ---#
#--- Scalars ---#
output_version=0
kp=0.
wp=0.
dist=0.
np=0.
dist_um=0.

Nr=0
Nz=0

#--- Vectors ---#
r_mesh = []
z_mesh = []

#--- Matrixes ---#
rho=[]; rho_b=[]; rho_bck=[]
Er=[]; Er_b=[]; Er_bck=[]
Ez=[]; Ez_b=[]; Ez_bck=[]
Bphi=[]; Bphi_b=[]; Bphi_bck=[]; Bphi_ex=[]
Jr=[]; Jr_b=[]; Jr_bck=[]
Jz=[]; Jz_b=[]; Jz_bck=[]
Zstar=[]
rho_i=[]

#--- *** Phase Space ***---#
# Scalar #
n_bunches=0
number_of_particles=0
#vectorial
bunch_charges=[]
x=[]; y=[]; z=[];
px=[]; py=[]; pz=[]
bunch_id=[]; cut=[];  dcut=[]
macro_particle_charge=[]; macro_particle_nume=[];
