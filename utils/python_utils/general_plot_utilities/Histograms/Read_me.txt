Tools_CALDER-CIRC usage:

python Python_scripts_folder/command


List of commands:

--------------

Architect_Plot_Spectrum1D.py R/W horiz_axis iter nbins/bins_width bins_input
R: reads pre-existing spectrum file - W: writes new spectrum file, overwrites existing one

Valid options for axis:
X, Y, Z, Px, Py, Pz, Nbunch, cut, dcut

If W option is chosen, then a hisogram is created and plotted from the phase space file. If nbins is specified the histogram has the chosen number 
of bins in horizontal axis. If bins_width are specified, the width of each bin in the horizontal axis is specified. 
Bad choices for this value could create errors or inconsistent results. Looking at the resolution suggestion printed at screen is strongly adviced. 

--------------

Architect_Plot_Spectrum2D.py R/W horiz_axis vert_axis iter nbins/bins_width binsx_input binsy_input
R: reads pre-existing spectrum file - W: writes new spectrum file, overwrites existing one

Valid options for axes:
X, Y, Z, Px, Py, Pz, Nbunch, cut, dcut

If W option is chosen, then a hisogram is created and plotted from the phase space file. If nbins is specified the histogram has the chosen number 
of bins in horizontal and vertical axes. If the latter is not specified it is taken as the same of the horizontal axis. If bins_width are specified,
the width of each bin in each axis is specified. Bad choices for this value could create errors or inconsistent results. 
Looking at the resolution suggestion printed at screen is strongly adviced. 
