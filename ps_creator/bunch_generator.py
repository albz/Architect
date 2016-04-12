#!/usr/bin/python
######################################################################
# Name:         bunch_generator
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
###>>>


#- dg "is a %-value!

###- function bunch generator -###
def bunch_generator_function( shape, s_z, s_x, s_y, x_0, y_0, z_0, ex, ey, gamma, dg,szg, np, verbose): 

	#- s_x -> is the standard deviation along X: in the code it uses 2*s_x
	#- s_y -> is the standard deviation along Y: in the code it uses 2*s_y
	#- s_z -> is the standard deviation along Y: in the code it uses 2*s_z
	s_x = 2.*s_x; s_y = 2.*s_y; s_z = 2.*s_z
	
	
	#print shape,nlong,ntran,s_z,s_x,s_y,cut,ex,ey,gamma,dg,np
	shape = int(shape); s_z = float(s_z); s_x = float(s_x)
	s_y = float(s_y); ex = float(ex); ey = float(ey); gamma = float(gamma); dg = float(dg); np = int(np)


	### Initialize arrays ##
	pt		= numpy.zeros([6], float)
	a 		= numpy.zeros([6], float)
	bunch 	= numpy.zeros(shape=(np,6))
	fct1 	= 1./(2.*numpy.pi)

	### Creates the 6D Phase Space ##
	# shape == 1 ---> gaussian distribution
        # shape == 2 ---> gaussian distribution, z-Energy correlated


	if shape==1: ##gaussian distribution
		sigs = numpy.zeros([6], float)

		sigs[0] = 0.5*s_x
		sigs[1] = ex/sigs[0]
		sigs[2] = 0.5*s_y
		sigs[3] = ey/sigs[2]
		sigs[4] = 0.5*s_z
		sigs[5] = (0.01)*dg*gamma
	
		cut = 3.0
		p_cut = ( fct1**3/(sigs[0]*sigs[1]*sigs[2]*sigs[3]*sigs[4]*sigs[5]) )*math.exp(-3.0*cut**2.0)
		#pmax = fct1**3.0/(sigs[0]*sigs[1]*sigs[2]*sigs[3]*sigs[4]*sigs[5])
	
		#--- loop generation
		for i in range(0,np):
			if i%1000==0: 
				print i

			#test=1.; p1=0.
			#while test > p1  or p1 < p_cut:
			p1=0.
			while p1 < p_cut:
				pt=numpy.random.randn(6)*sigs
				#test = pmax * numpy.random.random()
				p1 = numpy.abs( numpy.prod( pt ) ) - p_cut
			
			bunch[i,0:6] = pt
			bunch[i,5] = bunch[i,5] + gamma
	
	

        elif shape==2: #gaussian distribution, z-Energy correlated
		

		if (   szg**2. >= (   s_z**2.  *  (  (0.01)*dg*gamma  )**2. )    ):
			print ' ERROR: Input sz, dg, szg are not valid - remember that covariance matrix of z,E must be positive definite '
			sys.exit()


		sigs = numpy.zeros([6], float)

		sigs[0] = 0.5*s_x
		sigs[1] = ex/sigs[0]
		sigs[2] = 0.5*s_y
		sigs[3] = ey/sigs[2]
		sigs[4] = 0.5*s_z
		sigs[5] = (0.01)*dg*gamma
	
		cut = 3.0
		p_cut = ( fct1**3/(sigs[0]*sigs[1]*sigs[2]*sigs[3]*sigs[4]*sigs[5]) )*math.exp(-3.0*cut**2.0)
		#pmax = fct1**3.0/(sigs[0]*sigs[1]*sigs[2]*sigs[3]*sigs[4]*sigs[5])
	

		mean = [x_0,0.,y_0,0.,z_0,0.]
		cov = [[sigs[0]**2.0,0.,0.,0.,0.,0.],[0.,sigs[1]**2.0,0.,0.,0.,0.],[0.,0.,sigs[2]**2.0,0.,0.,0.],[0.,0.,0.,sigs[3]**2.0,0.,0.],[0.,0.,0.,0.,sigs[4]**2.0,szg],[0.,0.,0.,0.,szg,sigs[5]**2.0]  ]

		
		

		#--- loop generation
		for i in range(0,np):
			if i%1000==0: 
				print i


			
			p1=0.
			
			while p1 < p_cut:
	
				x,px,y,py,z,pz=numpy.random.multivariate_normal(mean,cov,1).T			
							
				pt=[x,px,y,py,z,pz]
				p1 = numpy.abs( numpy.prod( pt ) ) - p_cut	

			

			bunch[i,0] = x
			bunch[i,1] = px
			bunch[i,2] = y
			bunch[i,3] = py
			bunch[i,4] = z
			bunch[i,5] = pz
			bunch[i,5] = bunch[i,5] + gamma
		

			
	else: 
		print 'Error: invalid shape parameter'


	##### Computes the mean value within the bunch of every generated coordinate #####
	
	xm 	= numpy.mean(bunch[:,0])
	ym 	= numpy.mean(bunch[:,2])
	zm 	= numpy.mean(bunch[:,4])
	pxm     = numpy.mean(bunch[:,1])
	pym     = numpy.mean(bunch[:,3]) 
	pzm     = numpy.mean(bunch[:,5])


	##### Computes the Co-variance within the bunch of every generated coordinate #####
	cov_x_pxm = numpy.cov(bunch[:,0],bunch[:,1])[0,1]
	cov_y_pym = numpy.cov(bunch[:,2],bunch[:,3])[0,1]
	cov_z_pzm = numpy.cov(bunch[:,4],bunch[:,5])[0,1]

	#### Transverse emittances #####
	emx = math.sqrt( numpy.var(bunch[:,0]) * numpy.var(bunch[:,1])-cov_x_pxm**2.0)
	emy = math.sqrt( numpy.var(bunch[:,2]) * numpy.var(bunch[:,3])-cov_y_pym**2.0)



	### --- output --- ###
	if verbose == True:
		print '--------------------'
		print 'input emittance_x = ',ex,'----generated bunch emittance_x = ',emx 
		print 'input emittance_y = ',ey,'----generated bunch emittance_y = ',emy 
		print 'input sigma_x = ',sigs[0],'----generated bunch sigma_x = ',numpy.std(bunch[:,0])
		print 'input sigma_y = ',sigs[2],'----generated bunch sigma_y = ',numpy.std(bunch[:,2])
		print 'input sigma_z = ',sigs[4],'----generated bunch sigma_z = ',numpy.std(bunch[:,4])
		print 'input sigma_px = ',sigs[1],'----generated bunch sigma_px = ',numpy.std(bunch[:,1])
		print 'input sigma_py = ',sigs[3],'----generated bunch sigma_py = ',numpy.std(bunch[:,3])
		print 'input dE/E = ',dg,'----generated bunch dE/E = ',100.0*math.sqrt(numpy.var(bunch[:,5]))/pzm
		print 'input E = ',gamma,'----generated bunch E = ',pzm

		if shape==2:
			print 'input sigma_z_pz = ',szg,'----generated bunch  sigma_z_pz = ',cov_z_pzm


	#### Centers the bunch in the Phase Space ######
	for i in range(1,np):
		bunch[i,0] = bunch[i,0] - xm
		bunch[i,2] = bunch[i,2] - ym
		bunch[i,4] = bunch[i,4] - zm				
		bunch[i,1] = bunch[i,1] - pxm		
		bunch[i,3] = bunch[i,3] - pym
		bunch[i,5] = bunch[i,5] - (pzm-gamma)


	#### Writes the Phase Space in output file ##########
	x	=	bunch[:,0]
	px	=	bunch[:,1]
	y	=	bunch[:,2]
	py	=	bunch[:,3]
	z	=	bunch[:,4]
	pz	=	bunch[:,5]


	#-end of the program
	return x,px,y,py,z,pz
	

	
	
###--- PRINTING bunch on a file
def print_bunch(s_z,s_x,s_y, xb_0, yb_0, zb_0, x,px,y,py,z,pz,file_name):

	full_file_name = os.path.join(os.getcwd(),file_name)
	
	if os.path.exists( full_file_name ) == True:
		os.remove( full_file_name )
	f = open(full_file_name,'a')
	f.close()

	parameters_header = numpy.column_stack((xb_0,yb_0,zb_0,s_x,s_y,s_z))
	Matrix = numpy.column_stack( (x,y,z, px,py,pz) )
	
	f = open(full_file_name,'a')
	numpy.savetxt( f ,parameters_header,fmt='%15.14e')
	numpy.savetxt( f ,Matrix,fmt='%15.14e')
	f.close()
	
	
	
	
	
