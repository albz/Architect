#!/usr/bin/python
#####################################################################
# Author:		F. Massimo
# Date:			2016-02-09
# Purpose:
# Source:       Python
#####################################################################

# This routine must be used inside an Architect simulation directory, 
# preferably with many output files in the out/2D subdirectory.
# All the section files are read and combined in a density colormap. 
# The colormaps act as frames for an animation, saved in the simulation folder 


### loading shell commands
import os, os.path, glob, sys, shutil
import time, datetime
import scipy
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pylab as pyl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as colors


from read_Architect_bin_section import *



#################################################################
def listdir_nohidden(path): # function to ignore hidden files
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f
#################################################################
# Customize the colormap
color_dictionary = {'red':  ((0.0, 0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.5, 0.8, 1.0),
                   (0.75, 1.0, 1.0),
                   (1.0, 0.4, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.5, 0.9, 0.9),
                   (0.75, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.4),
                   (0.25, 1.0, 1.0),
                   (0.5, 1.0, 0.8),
                   (0.75, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }
mycolormap = LinearSegmentedColormap('MyColormap', color_dictionary)
#################################################################
# Or choose subset of preexisting ones
def truncate_colormap(cmap, minval=0.2, maxval=0.5, n=70):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))

cmap = plt.get_cmap('seismic')
new_cmap = truncate_colormap(cmap, 0.2, 0.5)
#################################################################

matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath'] # to make LateX fonts bold

##### Initial setup
home_path = os.path.expanduser('~')
path      = os.getcwd()
path      = os.path.join(path,'out/2D')

Number_of_frames = len([file_name for file_name in listdir_nohidden(path) if os.path.isfile(os.path.join(path, file_name))])
#Number_of_frames = 10

# find matrix dimensions from first file in directory
f              = open(os.path.join(path, '0000000.arch'),'rb')
output_version = struct.unpack('i', f.read(4))[0]
dist           = struct.unpack('i', f.read(4))[0] # distance traveled
# radial and axial dimensions 
Nr             = struct.unpack('i', f.read(4))[0] # half-plane axis dim (radial direction      )
Nz             = struct.unpack('i', f.read(4))[0] # z dim 				(longitudinal direction)
f.close()

# --- Create 3D array for frames
density_at_frame = np.zeros((Number_of_frames,Nr,Nz))
dist_at_frame    = np.zeros( Number_of_frames       )

##### Read section files
i=0
for file_name in listdir_nohidden(path):
	print 'Reading frame ', i,' of ',Number_of_frames
	
	f              = open(os.path.join(path, file_name),'rb')
	output_version = struct.unpack('i', f.read(4))[0]
	f.close()
	
	if output_version == 1:			
		dist,r_mesh,z_mesh,rho_b,n_bck,Er,Ez,Bphi,Jbr,Jbckr,Jbz,Jbckz                                                       = read_Architect_bin_section_output_v_1(path,file_name)
	elif output_version == 2:
		dist,r_mesh,z_mesh,rho_b,n_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,Jbr,Jbckr,Jbz,Jbckz               = read_Architect_bin_section_output_v_2(path,file_name)
	elif output_version == 3:
		dist,r_mesh,z_mesh,rho_b,n_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,B_ex_poloidal,Jbr,Jbckr,Jbz,Jbckz = read_Architect_bin_section_output_v_3(path,file_name)
	

	density     = n_bck+rho_b
	max_density =  density.max()
	cbarlim     =  min(3.,max_density)
	density_at_frame[i,:,:] = np.minimum(density,cbarlim*np.ones((Nr,Nz)) )
	dist_at_frame[i]        = dist
	i += 1
 	#if(i==Number_of_frames): break


z_mesh,r_mesh = np.meshgrid(z_mesh,r_mesh) 
z_mesh=-z_mesh[::-1] # reverse longitudinal axis

##### Frames setup 
fig = plt.figure()

ax = plt.axes(xlim=(-250., 150.), ylim=(-40., 40.))  
plt.yticks([-40,-20,0,20,40], fontsize = 12.0)
plt.xticks([-200,-100,0,100], fontsize = 12.0)

# plt.xlabel(r'Z-ct ($\mu$m)',usetex=True,fontsize = 14.0)
# plt.ylabel(r'[ -R, R ] ($\mu$m)', usetex=True,fontsize = 14.0)
plt.xlabel(r'Z-ct ($\mu$m)',fontsize = 14.0)
plt.ylabel(r'[ -R, R ] ($\mu$m)',fontsize = 14.0)
plt.gcf().subplots_adjust(bottom=0.15)


# animation function, creates the plot for each iteration
def animate(i):
	max_density =  density_at_frame[i,:,:].max()
	cbarlim     =  min(3.,max_density)
	cont = plt.contourf(z_mesh, r_mesh, density_at_frame[i,:,:], np.linspace(0,cbarlim,50),cmap='hot')
	#plt.title(r'$n_e$/$n_0$,  Z = '+'{0:.2f}'.format(dist_at_frame[i]/1e4)+' cm', usetex=True, fontsize = 14.0,y=1.03 )
	plt.title(r'$n_e$/$n_0$,  Z = '+'{0:.2f}'.format(dist_at_frame[i]/1e4)+' cm',fontsize = 14.0,y=1.03 )
	return cont 

# initialization function, sets colorbar, opens the first plot...
def init():
	max_density =  density_at_frame[0,:,:].max()
	cbarlim     =  min(3.,max_density)
	cont = plt.contourf(z_mesh, r_mesh, density_at_frame[0,:,:], np.linspace(0,cbarlim,50),cmap='hot')
# 	cbar = fig.colorbar(cont,ticks=[0,1,2,3],cmap='hot')   
	return [fig]
	
	
# creates the animation
print 'Creating animation ...' 
anim = ani.FuncAnimation(fig, animate, np.arange(1, Number_of_frames),interval=25, blit=False, init_func=init)  

dpi = 200
h_in_inches = 4.
w_in_inches = 8.
fig.set_size_inches(w_in_inches, h_in_inches, forward=True)  # Your resolution is then dpi * h_in_inches X dpi * w_in_inches



##### Save the animation
print 'Saving animation ...'
writer = ani.writers['ffmpeg'](fps=5)
anim.save('Architect_video.mp4',writer=writer,dpi=dpi) # Obvius advice: Higher resolution needs more time to save
print 'Animation successfully saved'

###### Plot 
#plt.show()



