#!/usr/bin/python
######################################################################
# Name:         Architect_utilities.py
# Author:       F Massimo      
# Date:         2016-02-08
# Purpose:      utilities for CALDER-CIRC postprocessing
# Source:       Python
#####################################################################



import os, os.path, glob, sys, shutil, time, datetime, re


def generate_folder(path,output_data):
	directory = os.path.join(path,'out')
	directory = os.path.join(directory,'Histograms')		

	if not os.path.exists( directory ):
		os.makedirs(directory)
	
	directory = os.path.join(directory,output_data)
	if not os.path.exists( directory ):
		os.makedirs(directory)
	
	return directory	
		

		
		
		