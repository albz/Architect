#!/usr/bin/python
######################################################################
# Author:       A. Marocchino
# Date:			14-11-2017
# Purpose:      plot sliced analisys
# Source:       python
# --- --- --- #
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil
import scipy
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import pylab as pyl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes','Code_Architect','Architect','utils','python_utils','architect_graphycal_unit'))
# plt.style.use(os.path.join(os.path.expanduser('~'),'Codes','Python_general_controllers','python_plot','plot_style_ppth.mplstyle'))
from read_architect_bin import *
### --- ###

# --- inputs --- #
if(len(sys.argv)<5):
	print('\nnot enought input arguments')
	print('input 1 :: bunch select')
	print('input 2 :: select slice thinkness in um')
	print('input 3 :: select shifting dimension in um')
	print('input 4 :: frame number to analyse\n')
	sys.exit()
bunch_select          = int(sys.argv[1])
slice_thinkness_um    = float(sys.argv[2])
shifting_thinkness_um = float(sys.argv[3])
frm_number            = int(sys.argv[4])

# --- inputs --- #
path = os.getcwd()
read_architect_bin('PS',path,frm_number)

# --- select bunch --- #
select_bunch = (np.asarray(var.bunch_id)==bunch_select)
select_dcut  = (np.asarray(var.dcut)==1.0)
selected     = np.asarray(select_bunch & select_dcut)

x  = var.x[selected]
y  = var.y[selected]
z  = +(var.z[selected]-np.mean(var.z[selected]))
px = var.px[selected]
py = var.py[selected]
pz = +var.pz[selected]
# W  = var.macro_particle_nume[selected]
#--- *** ---#

# --- particle selector --- #
W          = np.full(len(x), 1.0, dtype=float)

fig = pyl.figure(1)
fig.set_size_inches(3.25,3.6,forward=True)
ax1  = pyl.subplot(111)
#--- sigma ---#
ax1.plot(z,x, 'o', markersize=1.)
plt.show()

#--- whole bunch integrated parameters ---#
# Charge  =  var.bunch_charges[bunch_select-1]*1e-9
Charge=np.array([0.030e-9,0.030e-9])
gamma = np.array(np.sqrt(1. + px**2 + py**2 + pz**2))
mu_x=np.average(x,weights=W)
mu_y=np.average(y,weights=W)
mu_z=np.average(z,weights=W)
mu_Px=np.average(px,weights=W)
mu_Py=np.average(py,weights=W)
mu_Pz=np.average(pz,weights=W)
mu_gamma=np.average(gamma,weights=W)
sigma_x=np.sqrt( np.average((x-mu_x)**2, weights=W) )
sigma_y=np.sqrt( np.average((y-mu_y)**2, weights=W) )
sigma_z=np.sqrt( np.average((z-mu_z)**2, weights=W) )
sigma_Px=np.sqrt( np.average((px-mu_Px)**2, weights=W) )
sigma_Py=np.sqrt( np.average((py-mu_Py)**2, weights=W) )
sigma_Pz=np.sqrt( np.average((pz-mu_Pz)**2, weights=W) )
sigma_gamma=np.sqrt( np.average((gamma-mu_gamma)**2, weights=W) )
cov_x_Px=np.average((x-mu_x)*(px-mu_Px), weights=W)
cov_y_Py=np.average((y-mu_y)*(py-mu_Py), weights=W)
cov_z_Pz=np.average((z-mu_z)*(pz-mu_Pz), weights=W)

en_spread = sigma_gamma/mu_gamma
emittance_x = np.sqrt( sigma_x**2*sigma_Px**2-cov_x_Px**2);
emittance_y = np.sqrt( sigma_y**2*sigma_Py**2-cov_y_Py**2);
emittance_z = np.sqrt( sigma_z**2*sigma_Pz**2-cov_z_Pz**2);
# Current = Charge*1e-15*3e8/np.sqrt(2*np.pi)/sigma_x/1e-6

print(' ### ---------------------------- ### ')
print('Diagnostic for the whole bunch')
print('Energy spread: ', ('%25.8e' % en_spread) ,'%')
print('Normalized Emittance X: ', ('%3.2e' % emittance_x), 'mm-mrad')
print('Normalized Emittance Y: ', ('%3.2e' % emittance_y), 'mm-mrad')
print('Mean Energy: ',  ('%3.2e' % (mu_gamma*0.511)), 'MeV')
print('transverse   - Sigmax:',('%3.2e' % sigma_x),'mum')
print('transverse   - Sigmay:',('%3.2e' % sigma_y),'mum')
print('longitudinal - Sigmaz:',('%3.2e' % sigma_z),'mum')
print(' ### ---------------------------- ### ')


# --- slice analisys --- #
first_particle_z =np.min(z)
last_particle_z  =np.max(z)
window_min       =first_particle_z
window_max       =window_min+slice_thinkness_um
slice_counter    =0

slice_pos=[]; emittance_x_plot=[]; emittance_y_plot=[]; energy_spread_plot=[]; current_plot=[]

while (last_particle_z>window_max):
	Z_selected     = np.asarray(np.asarray(z>=window_min) & np.asarray(z<window_max))

	gamma = np.array(np.sqrt(1. + px[Z_selected]**2 + py[Z_selected]**2 + pz[Z_selected]**2))
	mu_x=np.average(x[Z_selected],weights=W[Z_selected])
	mu_y=np.average(y[Z_selected],weights=W[Z_selected])
	mu_z=np.average(z[Z_selected],weights=W[Z_selected])
	mu_Px=np.average(px[Z_selected],weights=W[Z_selected])
	mu_Py=np.average(py[Z_selected],weights=W[Z_selected])
	mu_Pz=np.average(pz[Z_selected],weights=W[Z_selected])
	mu_gamma=np.average(gamma,weights=W[Z_selected])
	sigma_x=np.sqrt( np.average((x[Z_selected]-mu_x)**2, weights=W[Z_selected]) )
	sigma_y=np.sqrt( np.average((y[Z_selected]-mu_y)**2, weights=W[Z_selected]) )
	sigma_z=np.sqrt( np.average((z[Z_selected]-mu_z)**2, weights=W[Z_selected]) )
	sigma_Px=np.sqrt( np.average((px[Z_selected]-mu_Px)**2, weights=W[Z_selected]) )
	sigma_Py=np.sqrt( np.average((py[Z_selected]-mu_Py)**2, weights=W[Z_selected]) )
	sigma_Pz=np.sqrt( np.average((pz[Z_selected]-mu_Pz)**2, weights=W[Z_selected]) )
	sigma_gamma=np.sqrt( np.average((gamma-mu_gamma)**2, weights=W[Z_selected]) )
	cov_x_Px=np.average((x[Z_selected]-mu_x)*(px[Z_selected]-mu_Px), weights=W[Z_selected])
	cov_y_Py=np.average((y[Z_selected]-mu_y)*(py[Z_selected]-mu_Py), weights=W[Z_selected])
	cov_z_Pz=np.average((z[Z_selected]-mu_z)*(pz[Z_selected]-mu_Pz), weights=W[Z_selected])

	en_spread = sigma_gamma/mu_gamma
	emittance_x = np.sqrt( sigma_x**2*sigma_Px**2-cov_x_Px**2);
	emittance_y = np.sqrt( sigma_y**2*sigma_Py**2-cov_y_Py**2);
	emittance_z = np.sqrt( sigma_z**2*sigma_Pz**2-cov_z_Pz**2);

	# charge      = var.bunch_charges[bunch_select-1]*1e-9 / sum(W) * sum(W[Z_selected])
	charge      = 0.030e-9 / sum(W) * sum(W[Z_selected])
	# charge      = 1. #var.bunch_charges[Z_selected-1]*1e-9/len(var.x)*sum(Z_selected)
	current     = charge / (slice_thinkness_um*1e-6) * 3e8

	#print( "%4d %8d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" % (slice_counter,sum(Z_selected),window_min,window_max,sigma_x,sigma_y,emittance_x,emittance_y,en_spread*1e3,current,mu_gamma) )
	print( "%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" % (0.5*(window_min+window_max),sigma_x,sigma_y,emittance_x,emittance_y,en_spread*1e3,current,mu_gamma) )

	slice_pos.append((window_max+window_min)/2.)
	emittance_x_plot.append(emittance_x)
	emittance_y_plot.append(emittance_y)
	energy_spread_plot.append(en_spread)
	current_plot.append(current)

	window_min+= shifting_thinkness_um
	window_max+= shifting_thinkness_um
	slice_counter+=1

# --- plot --- #
fig = pyl.figure(1)
fig.set_size_inches(4.2,4.0,forward=True)
ax1  = pyl.subplot(211)
ax2  = pyl.subplot(212)

ax1.plot(z,x, '.', markersize=0.2, lw=1, label=r"Driver")
ax2.plot(np.array(slice_pos),np.array(emittance_x_plot), lw=1, label=r"$\varepsilon_x$")
ax2.plot(np.array(slice_pos),np.array(emittance_y_plot), lw=1, label=r"$\varepsilon_x$")
ax2.plot(np.array(slice_pos),np.array(current_plot)/1e3, lw=1, label=r"I")
ax2.plot(np.array(slice_pos),np.array(energy_spread_plot)*1e3, lw=1, label=r"$\sigma_E \times 1000$")

ax1.set_ylabel('X ($\mu$m)', labelpad=-5)
ax2.set_ylabel('arb. units')
ax2.set_xlabel('Z ($\mu$m)', labelpad=-3)

ax2.legend(loc='upper center', bbox_to_anchor=(0.3, -0.27),fancybox=True, shadow=False, ncol=2)
pyl.subplots_adjust(top=0.98,bottom=0.3)
# plt.tight_layout()

# pyl.savefig(os.path.join(path,'witness_rollingwindow_in.pdf'), format='pdf')

plt.show()
