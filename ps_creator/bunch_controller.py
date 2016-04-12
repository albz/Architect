#!/usr/bin/python
######################################################################
# Name:         bunch_controller
# Date:			12-03-2014
# Purpose:      Initial Phase Space Builder for QFluid bunches
# Source:       python
#####################################################################
import os, os.path, sys
import scipy, pylab, math
import matplotlib.pyplot as plt
import numpy as numpy
import random as random
###>>>
sys.path.append(os.getcwd())
from bunch_generator import *
###>>>




#- Generating a single bunch
s_z = 50.#15.0#50. #10.#25.
s_x = 8. #8.0#8.0 #5.#60.
s_y = s_x
x_0 = 0.0
y_0 = 0.0
z_0 = 0.0
ex  = 1.
ey  = 1.
szg = 0.#120.0
gamma = 200. #235.29#2000.
dg  = 0.1 #D(Energy) in %
np  = 1500000 #30000
Verbose = True



#First argument of bunch_generator_function:
#  1: bivariate gaussian, z-Energy uncorrelated 
#  2: bivariate gaussian, z-Energy correlated 



x,px, y,py, z,pz = bunch_generator_function( 1, s_z, s_x, s_y, x_0, y_0, z_0, ex, ey, gamma, dg,szg, np, Verbose) 


print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'rugby_test_sz50um_sx8um_1500k.inp')
#print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'bunch_01.inp')
#print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'driver_test.inp')
#print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'driver_20umsz_13umsx_ALaDyn_1GeV.inp')
#print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'driver_20umsz_13umsx_1GeV.inp')
#print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'witness_bench.inp')
#print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'driver_test_15umsz_25umsx.inp')
#print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'driver_01.inp')
#print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'driver_30k_test1.inp')
#print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'witness_bench_50particles.inp')
#print_bunch(s_z, s_x, s_y, x_0, y_0, z_0, x,px, y,py, z,pz, 'short_driver_1500k_ALaDyn_comparison.inp')
