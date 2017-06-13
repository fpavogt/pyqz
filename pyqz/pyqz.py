# -*- coding: utf-8 -*-
'''
 pyqz: a Python module to derive the ionization parameter and oxygen abundance of HII 
 regions, by comparing their strong emission line ratios with MAPPINGS photoionization 
 models.\n
 Copyright (C) 2014-2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains the master pyqz routines. See the dedicated website for more info:
 
 http://fpavogt.github.io/pyqz

 Updated April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
 
 -----------------------------------------------------------------------------------------
  
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
'''
# ----------------------------------------------------------------------------------------

import sys
import warnings
import numpy as np
import scipy
from scipy import interpolate
import scipy.stats as stats
import itertools
import pickle
import multiprocessing
ncpu_max = multiprocessing.cpu_count()

# For the Kernel Density Estimator (don't force it if it's not there, but issue
# a Warning)
try:
    import statsmodels.api as sm 
except:
    warnings.warn("Statsmodels module not found. KDE_method must be set to "+
                  "'gauss' or else I will crash.")
                    
import matplotlib._cntr as cntr # To construct the KDE contours with plotting anything
from matplotlib.path import Path

# For generic things
import os
from datetime import datetime as dt
import time

# Get some important metadata, the code version and nifty tools ...
import pyqz_metadata as pyqzm
from pyqz_metadata import __version__
import pyqz_tools as pyqzt

# ------------------------------------------------------------------------------

# A function to get the reference grid - useful to make plots !
def get_grid(ratios, coeffs=[[1,0],[0,1]],
                Pk = 5.0, kappa=np.inf, struct='pp', sampling = 1):
    ''' 
    Returns a given line ratio diagnostic grid generated using MAPPINGS for a given set of
    parameters.
        
    :Args:
        ratios: string 
                The line ratios defining the grid, e.g. '[NII]/[SII]+;[OIII]/Hb'
        coeffs: list of list [default: [[1,0],[0,1]] ]
                The different coefficients with which to mix the line ratios. 
                The size of each sub-list must be equal to the number of line ratios 
                involved. Used for projected 3D diagnostic grids.
        Pk: float [default: 5.0] 
            MAPPINGS model pressure. This value must match an existing grid file.
        kappa: float [default: np.inf 
               The kappa value. This value must match an existing grid file.
        struct: string [default: 'pp']
                spherical ('sph') or plane-parallel ('pp') HII regions. This value 
                must match an existing reference grid file.
        sampling: int [default: 1]
                  Use a resampled grid ?
        
    :Returns: 
        out: list
             Returns a list [grid,grid_cols,metadata], where:
             grid [numpy array]: the diagnostic grid as a numpy array.
             grid_cols [list]: labels of each column inside grid.
             metadata [list]: basic info from the associated MAPPINGS simulations.
    '''
     
    # 0) Do some basic tests on the input
    if kappa == 'inf':
        kappa = np.inf
  
    # 1) Get the metadata about the file
    fn = pyqzt.get_MVphotogrid_fn(Pk = Pk, kappa = kappa, struct = struct, 
                                       sampling = sampling)
    metadata = pyqzt.get_MVphotogrid_metadata(fn)

    # Extract the individual line ratios
    ratios = ratios.split(';')   
            
    # Check the line ratio name, just in case ...
    for ratio in ratios:
        if not(ratio in metadata['columns']):
            raise Exception('Line ratio unknown: %s' % ratio) 
    # Also check that the coeffs size matches the number of line ratios
    if len(coeffs) !=2 or len(coeffs[0]) != len(ratios) or \
                                            len(coeffs[1]) != len(ratios) :
        raise Exception('Mixing coefficient error (size mismatch): %s' % coeffs)        

    # 2) Get the grid in a numpy array format - by default, maintain the grid
    # nodes last, and fill the first columns with stuff that may be of interest
       
    # For now, export Q, Z tot and Z gas. Could add more later ...
    data_cols = ['LogQ', 'Tot[O]+12','gas[O]+12']+ratios+['Mix_x','Mix_y']
    data = np.loadtxt(fn, comments='c', delimiter=',',skiprows=4+3*(sampling>1), 
                        usecols = [metadata['columns'].index(m) for m in data_cols[:-2]])
                        
    # And add space for the combined line ratios as well
    data = np.append(data,np.zeros((len(data),2)),axis=1)                    

    # 3) Do the mix of line ratio requested 
    for k in range(len(ratios)): 
        data[:,-2] += coeffs[0][k]*data[:,data_cols.index(ratios[k])]
        data[:,-1] += coeffs[1][k]*data[:,data_cols.index(ratios[k])]
    
    # 4) Sort it for clarity (first LogQ, then [O]+12 ...)
    sort_ind = np.lexsort((data[:,1],data[:,0]))
    
    # 5) Send it back
    return [data[sort_ind],data_cols,metadata]
# ------------------------------------------------------------------------------

# This function is designed to check a given line ratio grid, and return the bad segments
# if any is found.
def check_grid(ratios, coeffs = [[1,0],[0,1]],
                Pk = 5.0, kappa=np.inf, struct='pp', sampling=1):
    '''
    Creates the diagram of a given line ratio grid for rapid inspection, and identify the 
    location of "wraps" (= fold-over ambiguous regions).
        
    :Args:
        ratios: string 
                The line ratios defining the grid, e.g. '[NII]/[SII]+;[OIII]/Hb'
        coeffs: list of list [default: [[1,0],[0,1]] ]
                The different coefficients with which to mix the line ratios. 
                The size of each sub-list must be equal to the number of line ratios 
                involved. Used for projected 3D diagnostic grids.
        Pk: float [default: 5.0] 
            MAPPINGS model pressure. This value must match an existing grid file.
        kappa: float [default: np.inf 
               The kappa value. This value must match an existing grid file.
        struct: string [default: 'pp']
                spherical ('sph') or plane-parallel ('pp') HII regions. This value 
                must match an existing reference grid file.
        sampling: int [default: 1]
                  Use a resampled grid ?

    :Returns: 
        out: list
             A list of bad segments with the associated node coordinates.
   ''' 

    # 0) Let's get the data in question
    [grid, grid_cols, metadata] = get_grid(ratios, 
                                           coeffs=coeffs, Pk=Pk, kappa=kappa,     
                                           struct=struct, sampling = sampling)

    # 1) Look for problematic regions
    # Check as a function of LogQ and Tot[O]+12. Where are these ?
    u = grid_cols.index('LogQ')
    v = grid_cols.index('Tot[O]+12')
    grid_segs = []
    bad_segs = []

    for i in [u,v]:  
        # Here, 'var' plays the role of 'q' or 'z' depending on 'i'.
        for var in np.unique(grid[:,i]):

            # Extract the segments, check if they cross any other
            n_seg = len(grid[grid[:,i]==var])-1
            
            for s in range(n_seg):
                # a) check if segment clashes with existing ones [A,B]
                seg_a = [[grid[grid[:,i]==var]
                                    [:,grid_cols.index('Mix_x')][s],
                              grid[grid[:,i]==var]
                                    [:,grid_cols.index('Mix_y')][s]],
                              [grid[grid[:,i]==var]
                                    [:,grid_cols.index('Mix_x')][s+1],
                              grid[grid[:,i]==var]
                                    [:,grid_cols.index('Mix_y')][s+1]],
                        ]

                # Loop through stored segments
                for seg_b in grid_segs:
                    # Check that they are disconnected
                    if not(seg_a[0] in seg_b) and not(seg_a[1] in seg_b):
                        # Check if they intercept
                        if pyqzt.seg_intersect(seg_a[0],seg_a[1],seg_b[0],seg_b[1]):
                            # Yes, then store both as bad segments
                            bad_segs.append(seg_a)
                            bad_segs.append(seg_b)
                            
                # b) add it to the general list 
                grid_segs.append(seg_a)

    # Let's order a bit my list of bad segments
    if len(bad_segs) >0:
        # Make sure each segment is present only once 
        # (mind you, it's a list of list, which requires some fancy magic)
        bad_segs.sort()
        bad_segs = list(bad_segs for bad_segs,_ in itertools.groupby(bad_segs))
            
    return bad_segs
# ------------------------------------------------------------------------------

# The core function - returns 'q' or'z' for a given ratio (and a given grid !) 
# Interpolate case by case for even better results !
def interp_qz ( qz, 
                data,
                ratios,
                coeffs =[[1,0],[0,1]],
                Pk = 5.0, kappa=np.inf, struct='pp',sampling=1,
                method = 'linear',
                ):
    ''' The core function of pyqz.
    
    Returns the 'q' or 'z' value for a given set of line ratios based on a given 2-D \
    diagnostic grid.
    
    :Args:    
        
        qz: string 
            Which estimate to return; 'LogQ', 'Tot[O]+12' or 'gas[O]+12'  
        data: list of numpy array
              List of Arrays of the line ratio values. One array per line ratio.
        ratios: string 
                The line ratios defining the grid, e.g. '[NII]/[SII]+;[OIII]/Hb'
        coeffs: list of list [default: [[1,0],[0,1]] ]
                The different coefficients with which to mix the line ratios. 
                The size of each sub-list must be equal to the number of line ratios 
                involved. Used for projected 3D diagnostic grids.
        Pk: float [default: 5.0] 
            MAPPINGS model pressure. This value must match an existing grid file.
        kappa: float [default: np.inf 
               The kappa value. This value must match an existing grid file.
        struct: string [default: 'pp']
                spherical ('sph') or plane-parallel ('pp') HII regions. This value 
                must match an existing reference grid file.
        sampling: int [default: 1]
                  Use a resampled grid ?
        method: string [default: 'linear']
                'linear' or 'cubic' interpolation method of scipy.interpolate.griddata. 

    :Returns: 
        out: numpy array
             The computed estimates, with nans when the points land outside the grid.  
        
    :Notes:
        WARNING: using method = 'cubic' can lead to estimates slightly "outside" the grid
        area. Proceed with caution. Usually, it is better to use a resampled grid with 
        method = 'linear'.
    '''

    if not(qz in ['LogQ','Tot[O]+12', 'gas[O]+12']) :
        raise Exception('Unknown qz string: %s' % qz)

    for dat in data:
        if type(dat) != np.ndarray:
            # 05.2017 [F.P.A. Vogt]: fixed bug reported by P. Weilbacher.
            raise Exception('Input line ratios must be Numpy arrays and not: %s' % 
                        type(dat))
        if np.shape(dat) != np.shape(data[0]):
            raise Exception('Input data arrays size mismatch !')    
    
    # Save the input shape for later - I'll reshape everything in 
    # 1-dimension for simplicity.
    input_data_shape = np.shape(data[0])
    input_data_size = np.int(np.array(input_data_shape).prod())
    for i in range(len(data)):
            data[i] = data[i].reshape(input_data_size)
                                                                      
    # Compute the combined "observed" ratios
    data_comb = [np.zeros_like(data[0]),np.zeros_like(data[1])]  

    for k in range(len(ratios.split(';'))): 
        data_comb[0] += coeffs[0][k]*data[k]
        data_comb[1] += coeffs[1][k]*data[k]
 
    # 1-1) Load the corresponding data file
    [the_grid, grid_cols, the_metadata] = get_grid(ratios, 
                                                    coeffs=coeffs, Pk=Pk, 
                                                    kappa=kappa,
                                                    struct=struct, 
                                                    sampling=sampling)

    # 1-2) Check for wraps and other bad segments
    bad_segs = check_grid(ratios, coeffs=coeffs, Pk=Pk, kappa=kappa, 
                            struct=struct, sampling=sampling)                     
            
    # 2-1) Now, get ready to do the interpolation-s ...  
    # Figure out over what range we are stepping ...
    var = grid_cols.index(qz)
    Qvar = grid_cols.index('LogQ')
    Zvar = grid_cols.index('Tot[O]+12')
    Qloops = np.sort(np.unique(the_grid[:,Qvar]))
    Zloops = np.sort(np.unique(the_grid[:,Zvar]))

    # 2-3) To store the final results - fill it with nan-s ...
    interp_data = np.zeros_like(data_comb[0])+np.nan

    # 2-4) Now, loop through every grid panel and do the interpolation ...
    # Go panel by panel and skip those with bad segments !
    for (L1,Qloop) in enumerate(Qloops[:-1]):
        for (L2,Zloop) in enumerate(Zloops[:-1]):
            
            # Get the indices of the columns with the line ratios ... 
            # To make our life easier later ...
            rat1_ind = grid_cols.index('Mix_x')
            rat2_ind = grid_cols.index('Mix_y')

            # 2-4a) Select the slice
            this_panel = the_grid[(the_grid[:,Qvar]>=Qloop) * 
                                (the_grid[:,Qvar]<=Qloops[L1+1])*
                                (the_grid[:,Zvar]>=Zloop) * 
                                (the_grid[:,Zvar]<=Zloops[L2+1]),:] 

            # 2-4a) Check if this box has any bad segments and if so, skip it !
            nedges = len(this_panel[:,grid_cols.index(ratios.split(';')[0])])
            if nedges != 4:
                raise Exception('Something went terribly wrong ... (L736)')
            
            seg_keys = {0:[0,1],
                        1:[1,3],
                        2:[3,2],
                        3:[2,0]}
            is_bad_panel = False 

            for i in range(nedges):
                # loop through each segments. Be careful about the node orders !
                # This is where "seg_keys" comes in ... assumes the grid is
                # sorted ! Should be the case, based on the get_grid fct.
                this_seg = [[this_panel[seg_keys[i][0],rat1_ind],
                             this_panel[seg_keys[i][0],rat2_ind],
                            ],
                            [this_panel[seg_keys[i][1],rat1_ind],
                             this_panel[seg_keys[i][1],rat2_ind],
                            ]]
                # Remember to check segments in BOTH directions ...
                if (this_seg in bad_segs) or \
                        ([this_seg[1],this_seg[0]] in bad_segs):
                    is_bad_panel = True
            # At least one segment is bad, and therefore the panel is bad ...
            # Skip it !
            if is_bad_panel:
                continue
            
            # Here, I check for the presence of data points
            # inside the grid panel. If there's nothing, no need to interpolate 
            # anything !
            path_codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,
                            Path.CLOSEPOLY,]
            # Mind the order of the nodes ... !                                              
            path_verts = this_panel[[0,1,3,2,0],:][:,[rat1_ind,rat2_ind]].tolist()            
            panel_path = Path(path_verts, path_codes)
            
            points_in_panel = \
                panel_path.contains_points([ (data_comb[0][k],data_comb[1][k]) for k in 
                                             range(len(data_comb[0]))],
                                           radius=0.001)
                                               
            # Only do the interpolation if I actually have some point inside 
            # this panel ...
            if np.any(points_in_panel):             
                
                # 2-4b) Stretch slice between 0 and 1 for better results
                # 01.2016: Is this causing to have estimates outside the panel bounds ?
                # YES, if the method is 'cubic' ... -> the strongest reason
                # to date to use only "linear" ?
                # Should not matter much, now that the interpolation is done on a 
                # panel-by-panel basis.
                this_stretched_panel = np.zeros_like(this_panel)
                xmin = np.min(this_panel[:,rat1_ind])
                dxmax = np.max(this_panel[:,rat1_ind]-xmin)
                ymin = np.min(this_panel[:,rat2_ind])
                dymax = np.max(this_panel[:,rat2_ind]-ymin)

                this_stretched_panel[:,rat1_ind] = (this_panel[:,rat1_ind]-xmin)/\
                                                                            dxmax
                this_stretched_panel[:,rat2_ind] = (this_panel[:,rat2_ind]-ymin)/\
                                                                            dymax
    
                # 2-4c) Also do it for the input ratios - only use those that 
                # are INSIDE the panel (or near by, to avoid edge effects)
                stretched_data = np.zeros_like([data_comb[0][points_in_panel],
                                                data_comb[1][points_in_panel]
                                               ],dtype = np.float)
                stretched_data[0] = (data_comb[0][points_in_panel]-xmin)/dxmax
                stretched_data[1] = (data_comb[1][points_in_panel]-ymin)/dymax

                # 2-4d) Interpolate at last!
                this_interp_data = interpolate.griddata(this_stretched_panel[:,
                                                                    [rat1_ind,
                                                                    rat2_ind]
                                                                    ],
                                        this_panel[:,var],
                                        (stretched_data[0],stretched_data[1]),
                                        method=method,
                                        fill_value=np.nan)
                
                # This is problematic when used with get_global_qz ...
                # 2-4e) Make sure that the data falls within the allowed range. Avoid
                # unexpected surprises ...
                # this_interp_data[(this_interp_data<pyqzm.QZs_lim[pyqzm.M_version][qz][0])] = -1
                # this_interp_data[(this_interp_data>pyqzm.QZs_lim[pyqzm.M_version][qz][1])] = -1
                
            
                # 3-4) Store the cleaned interpolated values in the final array
                interp_data[np.where(points_in_panel)[0][np.where(this_interp_data>0)[0]]] = \
                            this_interp_data[this_interp_data>0]

    # Finish it all, without forgetting to return an array of the same initial
    # dimension !
    return interp_data.reshape(input_data_shape)
# ------------------------------------------------------------------------------

# The function used by the multiprocessing module. Does everything for a single spectrum.
def get_global_qz_singlespec((j, final_cols, data, data_cols, which_grids,
                              ids, qzs, Pk, kappa, struct, sampling, error_pdf, 
                              srs, flag_level, KDE_method, KDE_qz_sampling,
                              KDE_do_singles, KDE_pickle_loc, verbose)):
    ''' The single-process function associated with get_global_qz().
        
    This function is used to spawn as many process as input spectra sent by the user
    to get_global_qz(). See http://stackoverflow.com/a/3843313 for more details.
    
    :Notes:    
        The data is fed as one big list. See get_global_qz() for details on each item.
    '''
    # ------------------------------------------------------------------------------------
    # 0) Set things up.
    
    # First, create the final storage structure
    final_data = np.zeros([# the individual estimates
                           len(qzs)*len(which_grids)+
                           # the individual KDE estimates and their errors
                           2*len(qzs)*len(which_grids)+
                           # the "direct" combined estimates + std
                           2*len(qzs)+
                           # the combined KDE estimates + errors
                           2*len(qzs)+
                           # the flag value
                           1+
                           # the number of random point landing offgrid
                           1 
                           ])*np.nan
    
    flag = ''
    
    # Do we want to keep "more", in order to make a plot afterwards ?
    # Let's store it in a dictionary.
    if KDE_pickle_loc:
       KDE_topickle = {}
       KDE_topickle['srs'] = srs
       KDE_topickle['qzs'] = qzs
       KDE_topickle['which_grids'] = which_grids
       KDE_topickle['final_cols'] = final_cols
       KDE_topickle['ids'] = ids 
       KDE_topickle['Pk'] = Pk 
       KDE_topickle['kappa'] = kappa
       KDE_topickle['struct'] = struct
       KDE_topickle['sampling'] = sampling
       KDE_topickle['KDE_method'] = KDE_method
    
    # Is srs = 0 ? Then no KDE !
    if srs == 0:
        noKDE = True
    # Are all the errors 0 ? Then this means no KDE !    
    elif np.all( data[[i for i in range(len(data_cols)) if 'std' in data_cols[i]]]==0):
        noKDE = True
    else:
        noKDE = False
    
    # ------------------------------------------------------------------------------------
    # 1) Generate the random data around this point
    
    nlines = len(data_cols)/2
    line_names = []
    rsf = np.zeros((srs+1, nlines))
        
    # Do it for all the lines separately <= Lines and errors are uncorrelated !
    for (l,line) in enumerate(data_cols):
        if not('std' in line):
            # Keep the line names for latter ...
            line_names.append(line)
            nline = len(line_names)-1
            # Add the line's mean
            rsf[0,nline] = data[l]
            
            # Avoid issues if this line is bad - make it a nan everywhere:
            if rsf[0,nline] <= 0.0 or np.isnan(rsf[0,nline]):
                rsf[1:,nline] = np.nan
            else:
                f0 = data[l] # Flux
                df = data[data_cols.index('std'+line)] # 1-std error
                if df > 0.0: # This is a normal error
                # Generate random fluxes following the error distribution
                # Careful here: if the errors are large, I may generate 
                # negative fluxes !
                # Avoid this by using a truncated normal distribution and
                # set the lower bound to 0
                    cut_a = (0.0 - f0)/df
                    # Set the upper bound to infty
                    cut_b = np.infty
                    cut_func = stats.truncnorm(cut_a,cut_b, 
                                            loc = f0,
                                            scale = df)
                    # generate some random numbers ...
                    rsf[1:,nline] = cut_func.rvs(srs)
                elif df == 0.0 :  
                    # If the error is 0, then this value is exact.
                    rsf[1:,nline] = f0
                elif df == -1.0:
                    # This flux is an upper limit -> draw fluxes with a 
                    # uniform distribution
                    rsf[1:,nline] = np.random.uniform(low=0.0,
                                                    high=f0,
                                                    size = srs)
                else:
                    raise Exception('Error: dFlux <0 and dFlux != -1 '+
                                'is not allowed. Check your errors !')      
    
    # Store these rsf if the user said so:
    if KDE_pickle_loc:
        KDE_topickle['rsf'] = rsf
        KDE_topickle['rsf_keys'] = line_names
    
    # ------------------------------------------------------------------------------------
    # 2) Compute individual q & z values for each of the "srs+1" points
    
    # Create a structure where I can store all the real and pseudo-estimates
    # for 'all' the diagnostics in questions.
    discrete_pdf = {}
                    
    for (i,calc) in enumerate(final_cols):    
        
        ratios = calc.split('|')[0]
        
        # Individual diagnostics
        if ratios in pyqzm.diagnostics and not('KDE' in calc) :
            qz = calc.split('|')[1]
            Rval = []

            for (k,ratio) in enumerate(ratios.split(';')):
                l1 = ratio.split('/')[0]
                l2 = ratio.split('/')[1]
            
                Rval.append(np.log10(rsf[:,line_names.index(l1)]/
                                     rsf[:,line_names.index(l2)]))
                
            coeffs = pyqzm.diagnostics[ratios]['coeffs']
            
            # launch the qz interpolation (10% total time)
            qz_values = interp_qz(qz, Rval, ratios, coeffs = coeffs, Pk = Pk, 
                                  kappa = kappa, struct = struct, sampling = sampling)
                                    
            # Store the 'real' line estimate                
            final_data[i] = qz_values[0]
            
            # Keep the other random estimates as well
            discrete_pdf[calc] = qz_values[:]
     
    # If requested, store these discrete values for separate plotting
    if KDE_pickle_loc:
        KDE_topickle['discrete_pdf'] = discrete_pdf
        
    # ------------------------------------------------------------------------------------    
    # 4) Calculate mean 'qz' values
    
    # First, look at the reconstructed joint probability density function
    
    # Create some storage structures
    all_estimates = {} # All the qz estimates
    all_bws = {} # All the bandwidths
    for qz in qzs:
        all_estimates[qz] = np.array([])
        all_bws[qz] = np.array([])
        
    # Total number of outside points and number of diagnostics
    all_my_rs_offgrid = [0.,0.]
    all_my_single_kernels = {}
        
    for (i,item) in enumerate(discrete_pdf.keys()):
        if 'LogQ' in item and not(np.isnan(discrete_pdf[item][0])):
            
            # the random points
            rs_points = {}
            rs_ref = {}
            
            for qz in qzs:
                rs_points[qz] = discrete_pdf[item.split('|')[0]+'|'+qz][1:]
                rs_ref[qz] = discrete_pdf[item.split('|')[0]+'|'+qz][0]

            # How many points are outside the grid ?
            n_off = len(rs_points['LogQ'][np.isnan(rs_points['LogQ'])])
            all_my_rs_offgrid[0] += n_off
            all_my_rs_offgrid[1] += 1.
                
            # Let's remove all nans first and add them to the global storage
            for qz in qzs:
                rs_points[qz] = rs_points[qz][~np.isnan(rs_points[qz])]          
                all_estimates[qz] = np.append(all_estimates[qz],
                                                    rs_points[qz])
            
            # Do we need to calculate the individual KDE ? 
            if KDE_method == 'gauss' and KDE_do_singles and not(noKDE):
                try:
                    values = np.vstack([rs_points[qzs[0]],rs_points[qzs[1]]])                
                    # Work out the kernel magic ...
                    my_kernel = stats.gaussian_kde(values, bw_method='scott')
                    # Save the kernel for later
                    all_my_single_kernels[item.split('|')[0]] = my_kernel
                except:
                      warnings.warn('Unexpected issue with single KDE: not enough points ?')
                      all_my_single_kernels[item.split('|')[0]] = np.nan
                      
            if KDE_method == 'multiv' and not(noKDE):
                # Calculate the bandwidths (along the 'qz' axis) for the 
                # sample
                # WARNING: calculate the BW manually to match the code
                # sm.nonparametric.bandwidths.bw_scott contains an 
                # 'IQR' thingy that differs from what KDE Multivariate is 
                # actually calling !    
                try:  
                    for qz in qzs:
                        all_bws[qz] = np.append(all_bws[qz],
                                            1.06*np.std(rs_points[qz])*
                                            len(rs_points[qz])**(-1./6.))
                                                
                    if KDE_do_singles:                            
                        my_kernel = sm.nonparametric.KDEMultivariate(
                                            data=[rs_points[qzs[0]],
                                                  rs_points[qzs[1]]],
                                            var_type='cc', 
                                            bw=np.array([all_bws[qzs[0]][-1],
                                                         all_bws[qzs[1]][-1]]) 
                                                                    )
                    all_my_single_kernels[item.split('|')[0]] = my_kernel                                     
                except:
                    warnings.warn('Unexpected issue with single KDE: not enough points ?')
                    all_my_single_kernels[item.split('|')[0]] = np.nan                                                                           
            
    # Store the number of points outside all the grids (cumulative)
    # only do this, if the central point is actually ON the grid !
    if all_my_rs_offgrid[1]>0:
        if noKDE:
            flag = -1
        else:
            final_data[final_cols.index('rs_offgrid')] = \
                    np.round(all_my_rs_offgrid[0]/(all_my_rs_offgrid[1]*srs)*100,1)
                        
    else:
        # No point lands in any of the grids ! Raise a flag for that.
        final_data[final_cols.index('rs_offgrid')] = np.nan
        flag += '8'
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # What is the acceptable range of the diagnostic ?   
    # QZs_lim is now defined in pyqz_metadata for easier access and clarity
    # WARNING: Different diagnostics have different areas ... 
    # I really should account for this eventually ...
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    # Discretize the QZ space ... 201 is a good compromise 
    # (speed vs accuracy)
    QZs_grid = {}
    for qz in qzs:
        QZs_grid[qz] = np.mgrid[pyqzm.QZs_lim[pyqzm.M_version][qz][0]:
                                pyqzm.QZs_lim[pyqzm.M_version][qz][1]:KDE_qz_sampling]

    # WARNING: HERE, I ASSUME THAT qzs has a length of only 2 !
    # Unfortunate, but that's all you get for now ...
    gridX, gridY = np.mgrid[pyqzm.QZs_lim[pyqzm.M_version][qzs[0]][0]:
                            pyqzm.QZs_lim[pyqzm.M_version][qzs[0]][1]:KDE_qz_sampling, 
                            pyqzm.QZs_lim[pyqzm.M_version][qzs[1]][0]:
                            pyqzm.QZs_lim[pyqzm.M_version][qzs[1]][1]:KDE_qz_sampling]
    grid_positions = np.vstack([gridX.ravel(), gridY.ravel()])

    # For speed purposes, we will only reconstruct the KDE over the area where
    # we have some QZ estimates (no need to calculate the PDF over the entire QZ space,
    # it's a waste of time - we know it's almost 0)
    # So, define a subgrid just around all the QZ estimates
    
    # Bug correction (01.2016, F.P.A. Vogt). What if the points is outside 'all'
    # the diagnostics ? Must handle this case as well.
    # Bug correction (05.2017, F.P.A. vogt). Looks like I messed up the boundary conditions
    # which showed up when updating to numpy 1.12 ...
    
    if not(np.all(np.isnan(all_estimates[qzs[0]]))):
        # 09.2016: fixed a crash if I have an estimate outside the allowed range
        if len(np.where(gridX[:,0]<all_estimates[qzs[0]].min())[0])==0:
            gridxmin = 0 #pyqzm.QZs_lim[pyqzm.M_version][qzs[0]][0]
        else:
            gridxmin = np.where(gridX[:,0]<all_estimates[qzs[0]].min())[0].max()
        if len(np.where(gridX[:,0]>all_estimates[qzs[0]].max())[0]) == 0:
            gridxmax = -1 #pyqzm.QZs_lim[pyqzm.M_version][qzs[0]][1]
        else:    
            gridxmax = np.where(gridX[:,0]>all_estimates[qzs[0]].max())[0].min()
        # 
        if len(np.where(gridY[0,:]<all_estimates[qzs[1]].min())[0])==0:
            gridymin = 0 #pyqzm.QZs_lim[pyqzm.M_version][qzs[1]][0]
        else:   
            gridymin = np.where(gridY[0,:]<all_estimates[qzs[1]].min())[0].max()
        if len(np.where(gridY[0,:]>all_estimates[qzs[1]].max())[0])==0:
            gridymax = -1 #pyqzm.QZs_lim[pyqzm.M_version][qzs[1]][1]
        else:
            gridymax = np.where(gridY[0,:]>all_estimates[qzs[1]].max())[0].min()
    
    else: # Hum, seems that I have no estimate at all ...
        noKDE = True
        gridxmin = 0 #gridX.max()
        gridxmax = -1 #gridX.min()
        gridymin = 0 #gridY.max()
        gridymax = -1 #gridY.min()
        
    subgridX = gridX[gridxmin:gridxmax,
                    gridymin:gridymax]
    subgridY = gridY[gridxmin:gridxmax,
                    gridymin:gridymax] 

    subgrid_positions = np.vstack([subgridX.ravel(), subgridY.ravel()])

    # This structure is used to store all the reconstructed (sampled) KDE
    all_my_gridZ = {}
        
    # Compute the global kernel
    if noKDE:
        all_my_single_kernels['all'] = np.nan
                
    elif KDE_method == 'gauss' and not(noKDE):
                        
        values = np.vstack([all_estimates[qzs[0]], all_estimates[qzs[1]]])                
        # Work out the kernel magic ...
        try:
            kernel = stats.gaussian_kde(values, bw_method='scott')
            all_my_single_kernels['all'] = kernel
            
        except:
            if not(noKDE):
                # Could happen if srs is too small, or if the errors are all
                # zero !
                    
                warnings.warn('Not enough valid srs points '+
                                '(spectra #%s) to run the KDE' % j)
                all_my_single_kernels['all'] = np.nan
       
    elif KDE_method == 'multiv' and not(noKDE):

        # Scipy doesn't let me set the bandwidth along both axes ... not 
        # good enough !
        # Use instead the statsmodel KDE routine, see 
        # http://statsmodels.sourceforge.net/stable/index.html
        # First, estimate the bandwith - can't use the default one because 
        # it oversmooths i.e. if the different diagnostics are compact but 
        # away from each other, the default relies on the overall std, 
        # which is too large. Instead, take the mean of the different
        # individual bandwith for each individual dataset). Or min ?
        # bw_min = np.min(np.array(all_my_bw),axis=0)
        bw_mean = np.array([np.nanmean(all_bws[qzs[0]]),
                            np.nanmean(all_bws[qzs[1]])])

        try:
            kernel = sm.nonparametric.KDEMultivariate(
                                            data=[all_estimates[qzs[0]],
                                                  all_estimates[qzs[1]]],
                                                var_type='cc', 
                                                bw=bw_mean)
            all_my_single_kernels['all'] = kernel

        except:
            if not(noKDE):
                # Could happen if srs is too small, or if the errors are all
                # zero !
                warnings.warn('WARNING: not enough valid srs points '+
                                '(spectra #%s)' % j)
                all_my_single_kernels['all'] = np.nan  
            
    # Very well. Now with all the singles and the one global kernel, 
    # let's reconstruct the PDF on the discretized QZ space 
           
    # ----- This is the slowest bit. Could I speed it up ? How ?  ----------                     
    for kern in all_my_single_kernels.keys():
        gridZ = np.zeros_like(gridX)*np.nan
        if type(all_my_single_kernels[kern]) == type(np.nan):
            pass
        elif KDE_method == 'gauss':
            if pyqzm.do_full_KDE_reconstruction:
                gridZ = np.reshape(all_my_single_kernels[kern](grid_positions), 
                                    gridX.shape)
            else:
                gridZ[gridxmin:gridxmax,gridymin:gridymax] = \
                np.reshape(all_my_single_kernels[kern](subgrid_positions),
                            subgridX.shape)
        
        
        elif KDE_method == 'multiv':
            if pyqzm.do_full_KDE_reconstruction:
                gridZ = np.reshape(all_my_single_kernels[kern].pdf(grid_positions), 
                                    gridX.shape)  
            else:                    
                gridZ[gridxmin:gridxmax,gridymin:gridymax] = \
                np.reshape(all_my_single_kernels[kern].pdf(subgrid_positions),
                            subgridX.shape)                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                               
        else:
           raise Exception('Something went terribly wrong ...(L1402)')                     
            
        # And store the grid for later
        all_my_gridZ[kern] = gridZ
    # ----------------------------------------------------------------------        
          
    if KDE_pickle_loc:
    
        KDE_topickle['gridX'] = gridX
        KDE_topickle['gridY'] = gridY
        KDE_topickle['all_gridZ'] = all_my_gridZ                   
                     
    # Now, extract some useful values from the KDE ... 
    for key in all_my_gridZ.keys(): 
        gridZ = all_my_gridZ[key]
        gridZ = np.rot90(gridZ)[::-1,:]
        # Find the location of the peak
        if np.any(~np.isnan(gridZ)):
            peak_loc = np.unravel_index(np.nanargmax(gridZ), gridZ.shape)
            # Normalize the array to the peak value
            peak_val = gridZ[peak_loc]
            gridZ/= peak_val
        else:
            peak_loc = [np.nan,np.nan]         

        if not(noKDE):
            
            # Right now, I want to compute the 0.61 level contour of the KDE, to estimate
            # some "errors" and "mean". THe contour function inside matplotlib is ideal, 
            # but it requirs to create a plot - which is incompatible with multiprocessing
            #. Instead, call directly the core function inside contour. 
            # This is not ideal, as it could change any time, but it is a LOT better than 
            # forcing a non-interactive backend on the user.
            # See e.g.: http://stackoverflow.com/questions/18304722/python-find-contour-lines-from-matplotlib-pyplot-contour
            
            xv,yv = np.meshgrid(QZs_grid[qzs[0]],QZs_grid[qzs[1]])
            # Now the magic function
            c = cntr.Cntr(xv,yv, gridZ) 
            
            paths = c.trace(pyqzm.PDF_cont_level) 
            npaths = len(paths)/2
            paths = [Path(path,codes=paths[npaths+p]) for 
                       (p,path) in enumerate(paths[:npaths])] # Only keep the path nodes.
             
        
        # Get the 'mean location of the 0.61 level contour
        # First, fetch the contour as path; only if contour actually exists
        if np.any(all_my_gridZ[key] == all_my_gridZ[key]):

            # If the path/contour is broken, use the one that contains 
            # the peak of the distribution
            for path in paths:
                if path.contains_point([QZs_grid[qzs[0]][peak_loc[1]],
                                        QZs_grid[qzs[1]][peak_loc[0]]]):
                    break

            vert = path.vertices
            mean_vert = (np.mean(vert[:,0]),np.mean(vert[:,1]))
            # For the error, take the max extent of the ellipse
            err_vert = (np.max(np.abs(vert[:,0]-mean_vert[0])),
                        np.max(np.abs(vert[:,1]-mean_vert[1])))

            # WARNING: what if the contour is not continuous ?
            # I.e. if multiple peaks exist ?
            # For now, simply raise a flag for the combined KDE case
            if key == 'all' and npaths>1:
                flag = '9'
            
        else:
            path = [np.nan,np.nan]
            vert = [np.nan,np.nan]
            mean_vert = [np.nan,np.nan]
            err_vert = [np.nan,np.nan]                
        
        # Save the values as appropriate
        if key =='all':
            for (k,qz) in enumerate(qzs):
                final_data[final_cols.index('<'+qz+'{KDE}>')] = mean_vert[k]
                final_data[final_cols.index('err('+qz+'{KDE})')] = err_vert[k]
        else:
            for (k,qz) in enumerate(qzs):
                final_data[final_cols.index(key+'|'+qz+'{KDE}')] = mean_vert[k]
                final_data[final_cols.index('err('+key+'|'+qz+'{KDE})')] = err_vert[k]

    # Now, look at the real data 'mean' and 'variance'
    which_qz = {}
    for qz in qzs:
        which_qz[qz] = []        
        
    # Start by finding which columns I have to average
    for (i,name) in enumerate(final_cols):
        for qz in qzs:
            try:
                if name.split('|')[1]==qz and not('KDE' in name):
                    which_qz[qz].append(True)
                else:
                    which_qz[qz].append(False)
            except:
                which_qz[qz].append(False)
                                
    # Now get the mean and std for each value
    basic_qz = {}

    for qz in qzs:
        
        mean_qz = np.nanmean(final_data[np.array(which_qz[qz])])
        std_qz = np.nanstd(final_data[np.array(which_qz[qz])])
        basic_qz[qz] = [mean_qz,std_qz]

        # Save these
        final_data[final_cols.index('<'+qz+'>')] = mean_qz
        final_data[final_cols.index('std('+qz+')')] = std_qz   

        # Do a quick check to see if qz and qz_rs are consistent
        # i.e. within each-other's errors
        check1 = np.abs(mean_qz-mean_vert[qzs.index(qz)])/\
                                std_qz<=flag_level
        check2 = np.abs(mean_qz-mean_vert[qzs.index(qz)])/\
                                err_vert[qzs.index(qz)]<=flag_level
            
        for (i,check) in enumerate([check1,check2]):
            if not(check) and not(noKDE):
                flag+=str((len(qzs)*np.int(qzs.index(qz))+i+1))

    if flag == '':
        flag+='0'  

    final_data[final_cols.index('flag')]= np.int(flag)      
    
    if KDE_pickle_loc:
        KDE_topickle['final_data'] = final_data
        
        fn_out = pyqzt.get_KDE_picklefn(ids, Pk, kappa, struct, sampling, KDE_method)
        # Append the path:
        fn_out = os.path.join(KDE_pickle_loc, fn_out)
        
        ftmp = open(fn_out,'w')
        pickle.dump(KDE_topickle,ftmp)
        ftmp.close() 
        
    return [j,final_data]
# ----------------------------------------------------------------------------------------

# Some small function to get the multiprocessing working fine.
# See http://stackoverflow.com/a/3843313
def get_global_qz_singlespec_init(q):
    
    ''' Some function required to get the multiprocessing working fine.
        See http://stackoverflow.com/a/3843313 for details.
    '''
    
    get_global_qz_singlespec.q = q
# ----------------------------------------------------------------------------------------

# Get the "global" Q/Z values from a set of line ratios and diagnostics
def get_global_qz(data, data_cols, which_grids,
                  ids = None,
                  qzs = ['LogQ','Tot[O]+12'],
                  Pk = 5.0, kappa = np.inf, struct = 'pp', sampling = 1,
                  error_pdf = 'normal', 
                  srs = 400,
                  flag_level = 2., # any float >0: flag_level*std = level above 
                                   # which a mismatch between q & z and q_rs & z_rs 
                                   # is flagged 
                  KDE_method='gauss', # Which KDE routine to use:
                                      # 'gauss'=gaussian_kde from Scipy: 
                                      # fast but cannot define 2D bandwidth
                                      # 'multiv' = 
                                      # statsmodel 'KDEMultivariate':
                                      # 100x slower but can fine-tune 2-D 
                                      # bandwidth
                  KDE_qz_sampling = 101j,
                  KDE_do_singles = True,   
                  KDE_pickle_loc = None,
                  verbose = True,
                  nproc = 1,
                  ):

    ''' 
    Derives the LogQ and [O]+12 values for a given set of line fluxes based on a given set 
    of diagnostic grids.
        
    This function only uses "valid" grids defined in pyqz_metadata.py.
          
    For each set of line ratios, this function combines the estimates from invdividual 
    diagnostics into a combined best estimate of LogQ and [O]+12. Observational errors are 
    propagated via the joint probability density function, reconstructed via a Kernel 
    Density Estimation from "srs" realizations of the measured line fluxes, constructed
    randomly according to the observational errors and probability distribution function.
    
    :Args:
        data: numpy array of size Nx2*M  
              N sets of M line fluxes and errors. An error =-1 implies that the associated 
              line flux is an upper limit.
        data_cols: list of Mx2 strings
                   Description of the coumns in "data", e.g. ['[NII]+','stdHa',...]
                   Mind the MAPPINGS + pyqz conventions.
        which_grids: list of strings
                     The list of the diagnostic grids to use for the estimations,
                     e.g. ['[NII]+/Ha;[OII]/Hb', ...]
        ids: list [default: None]
             An optional length-N list of numbers or strings to identify each data point 
             (e.g. spaxel number, source ID, etc ...)
        qzs: list of strings [default: ['LogQ','Tot[O]+12']]
             list of Q/Z values to compute: ['LogQ','Tot[O]+12'] or ['LogQ','gas[O]+12']
        Pk: float [default: 5.0] 
            MAPPINGS model pressure. This value must match an existing grid file.
        kappa: float [default: np.inf] 
               The kappa value. This value must match an existing grid file.
        struct: string [default: 'pp']
                spherical ('sph') or plane-parallel ('pp') HII regions. This value 
                must match an existing reference grid file.
        sampling: int [default: 1]
                  Use a resampled MAPPINGS grid ?
        error_pdf: string [default = 'normal']
                   The shape of the error function for the line fluxes. 
                   Currently, only 'normal' is supported.
        srs: int [default: 400]
             The number of random line fluxes generated to discretize (and reconstruct) 
             the joint probability function. If srs = 0, no KDE is computed.
        flag_level: float [default: 2]
                    A 'flag' is raised (in the output array) when the direct q/z estimates
                    and the KDE q/z estimates (q_rs and z_rs) are offset by more than 
                    flag_level * standard_deviation. Might indicate trouble.
                    A flag is always raised when a point is outside all of the grids (8), 
                    when the KDE PDF is multipeaked (9) or when srs = 0 
                    (-1, no KDE computed).
        KDE_method: string [default: 'gauss']
                    Whether to use scipy.stats.gaussian_kde ('gauss') or 
                    sm.nonparametric.KDEMultivariate ('multiv') to perform the 2-D Kernel 
                    Density Estimation.
        KDE_qz_sampling: complex [default: 101j]
                         The sampling of the QZ plane for the KDE reconstruction. 
        KDE_do_singles: bool [default: True]
                        Weather to compute KDE for single diagnostics (in addition to the 
                        KDE for all the diagnostics together).
        KDE_pickle_loc: string [default: None]
                        If specified, the location to save the pickle files generated for
                        each input points. Useful for subsequent plotting of the output 
                        via pyqz_plots.plot_global_qz()
        verbose: bool [default: True]
                 Do you want to read anything ?   
        nproc: int [default = 1]
               Defines how many process to use. Set to -1 for using as many as possible. 

    :Returns: 
            
        out: [P,r]
              A list where P contains all the estimated Log(Q) and [O]+12
              values, and r contains the columns name   
    '''

    if verbose:    
        print ' '
    
    # Let's keep track of time ...
    starttime = dt.now()
    
    # 0) Do some quick tests to avoid crashes later on ...
    try:
        if flag_level<=0:
            raise Exception('flag_level must be >0: %s' % flag_level)
    except:
        raise Exception('flag_level must be >0: %s' % flag_level)

    try:
        if len(qzs) != 2:
            raise Exception('qzs must be of length 2: %s' % qzs)
        if not('LogQ' in qzs):
            raise Exception("qzs must contain 'LogQ': %s" % qzs)
    except:
        raise Exception('qzs unrecognized: %s' % qzs)

    if not(KDE_method in ['gauss','multiv']):
        raise Exception('KDE_method unrecognized: %s' % KDE_method)

    if not(type(which_grids) == list):
        raise Exception('"which_grids" must be a list, and not: %s' % type(which_grids))
        
    for grid in which_grids:
        if not(grid in pyqzm.diagnostics.keys()): 
            raise Exception('Diagnostic grid unvalid: %s' % grid)
                            
    # Make sure the input array has a vaild 2D size
    if len(np.shape(data)) ==1:
        raise Exception('Bad data shape: "data" should be a 2-D array, but np.shape(data)=%s.' % 
                 np.shape(data))  
        
    # Finally, in nproc = None, set it to the max number of cpus. Helps against
    # weird warning ?
    if nproc == -1:
        nproc = ncpu_max
    elif nproc > ncpu_max:
        raise Exception('nproc set to %i, but I can find only %i cpus.' % (nproc, ncpu_max))
    
    # 3) Create the final storage stucture
    npoints = len(data)
    if verbose:
        if npoints > 1:
            print '--> Received '+np.str(npoints)+' spectra ...'
        else:
            print '--> Received 1 spectrum ...'
    
    # The final_data storage structure ... only contains floats for now ...        
    final_data = np.zeros([npoints,
                           # the individual estimates
                           len(qzs)*len(which_grids)+
                           # the individual KDE estimates and their errors
                           2*len(qzs)*len(which_grids)+
                           # the "direct" combined estimates + std
                           2*len(qzs)+
                           # the combined KDE estimates + errors
                           2*len(qzs)+
                           # the flag value
                           1+
                           # the number of random point landing offgrid
                           1 
                           ])*np.nan

    # The names of the different columns of the final_data matrix
    final_cols = []
    for grid in which_grids:
        for qz in qzs:
            final_cols.append(grid+'|'+qz)
    for qz in qzs:
        final_cols.append('<'+qz+'>')
        final_cols.append('std('+qz+')')        
            
    for grid in which_grids:
        for qz in qzs:
            final_cols.append(grid+'|'+qz+'{KDE}')
            final_cols.append('err('+grid+'|'+qz+'{KDE})')
    for qz in qzs:
        final_cols.append('<'+qz+'{KDE}>')
        final_cols.append('err('+qz+'{KDE})')
            
    # Also add a flag in case of strong mismatch between z & q and z_rs & q_rs
    final_cols.append('flag')
    
    # Add the number of (random) points landing outside the grid
    final_cols.append('rs_offgrid')

    # 3) Start doing the calculation
    # Let's use the multiprocessing package to speed things up
    # Especially when computing multiple points at once.        
    
    jobs = []
    
    # If the user did not provide any ids, let's just use numbers.
    if not(ids):
        ids = np.arange(npoints)+1
    
    # Let's create our list of jobs:
    for j in range(npoints):
        # This need all to be pickl-able ... (?)                   
        jobs.append( (j,final_cols, data[j,:], data_cols, which_grids,
                      ids[j],
                      qzs,
                      Pk, kappa, struct,sampling,
                      error_pdf, 
                      srs,
                      flag_level,
                      KDE_method,
                      KDE_qz_sampling,
                      KDE_do_singles,   
                      KDE_pickle_loc,                      
                      verbose) )                                    
              
    if len(jobs)>0:
        
        if nproc == 1:
            if verbose:
                print '--> Dealing with them one at a time ... be patient now !'
                print '    (no status update until I am done ...)' 
                
            results = map(get_global_qz_singlespec, jobs)    
            
        else:
            
            if verbose:
                print '--> Launching the multiple processes ... ' 
            
            # Build a queue, a pool and spwan processes    
            # (see below and http://stackoverflow.com/a/3843313)
            queue = multiprocessing.Queue()
            mypool = multiprocessing.Pool(nproc, get_global_qz_singlespec_init, [queue])
            rs = mypool.map_async(get_global_qz_singlespec,jobs)

            # We don't add anything else anymore
            mypool.close()
        
            # Let's track the overall progress, see 
            # http://stackoverflow.com/questions/5666576/show-the-progress-of-a-python-multiprocessing-pool-map-call
            chunksize = rs._chunksize
            if verbose:
                while (True):
                    time.sleep(0.5)
                    if (rs.ready()): break
                    remaining = rs._number_left
                    sys.stdout.write("    ... %i job(s) left ...\r" % (remaining*chunksize) )
                    sys.stdout.flush()

            
            # And for safety, make a "join" call ...
            mypool.join()
            
            results = rs._value
            
            if verbose:
                sys.stdout.write("    %i job(s) completed.      \n" % len(results))
                sys.stdout.flush()
            
    # All done. Now let's collect the results ...
    for j in range(len(jobs)):      

        final_data[j,:] = results[j][1]
            
        # Make sure the good spectrum is returned in the good order
        if results[j][0] != j :
            raise Exception('Something went very badly wrong ... (L1203)')

    # and return the final data set.
    if verbose:
        print ' '
        print 'All done in',dt.now()-starttime
    
    return [final_data, final_cols]
# ------------------------------------------------------------------------------

# Get the qz ratios from a file - mind the file structure !
def get_global_qz_ff(fn, 
                     which_grids,
                     qzs = ['LogQ', 'Tot[O]+12'],
                     Pk = 5.0,
                     kappa = np.infty,
                     struct = 'pp',
                     sampling = 1,
                     error_pdf = 'normal', 
                     srs = 400,
                     flag_level = 2., # any float >0: flag_level*std = level above 
                                    # which a mismatch between q & z and q_rs & 
                                    # z_rs is flagged 
                     KDE_method = 'gauss', # Which KDE routine to use:
                                           # 'gauss' from Scipy: fast 
                                           # but cannot define 2-D bandwidth
                                           # statsmodel 'multiv':
                                           # 100x slower but can fine-tune 2-D 
                                           # bandwidth
                     KDE_qz_sampling = 101j,
                     KDE_do_singles = True,     
                     KDE_pickle_loc = False,                 
                     decimals = 5,
                     missing_values='$$$', 
                     suffix_out = '_out',
                     verbose = True,
                     nproc = 1,
                     ):
    
    '''
    The get_global_qz-'from a file' function. Gets the line ratios from a csv 
    file, and not directly as Python arrays.
    
    Save to file the log(q) and 12+log(O/H) values for a given set of line 
    fluxes based on a given set of diagnostic grids. Requires a specific 
    input file, and can batch process several measurements automatically. 
    
    :Args:    
        
        fn: string
            The path+filename of the input file, in CSV format.   
        which_grids: list of strings
                     The list of the diagnostic grids to use for the estimations,
                     e.g. ['[NII]+/Ha;[OII]/Hb', ...]
        qzs: list of strings [default: ['LogQ','Tot[O]+12']]
             list of Q/Z values to compute: ['LogQ','Tot[O]+12'] or ['LogQ','gas[O]+12']
        Pk: float [default: 5.0] 
            MAPPINGS model pressure. This value must match an existing grid file.
        kappa: float [default: np.inf] 
               The kappa value. This value must match an existing grid file.
        struct: string [default: 'pp']
                spherical ('sph') or plane-parallel ('pp') HII regions. This value 
                must match an existing reference grid file.
        sampling: int [default: 1]
                  Use a resampled MAPPINGS grid ?
        error_pdf: string [default = 'normal']
                   The shape of the error function for the line fluxes. 
                   Currently, only 'normal' is supported.
        srs: int [default: 400]
             The number of random line fluxes generated to discretize (and reconstruct) 
             the joint probability function. If srs = 0, no KDE is computed.
        flag_level: float [default: 2]
                    A 'flag' is raised (in the output array) when the direct q/z estimates
                    and the KDE q/z estimates (q_rs and z_rs) are offset by more than 
                    flag_level * standard_deviation. Might indicate trouble.
                    A flag is always raised when a point is outside all of the grids (8), 
                    when the KDE PDF is multipeaked (9) or when srs = 0 
                    (-1, no KDE computed).
        KDE_method: string [default: 'gauss']
                    Whether to use scipy.stats.gaussian_kde ('gauss') or 
                    sm.nonparametric.KDEMultivariate ('multiv') to perform the 2-D Kernel 
                    Density Estimation.
        KDE_qz_sampling: complex [default: 101j]
                         The sampling of the QZ plane for the KDE reconstruction. 
        KDE_do_singles: bool [default: True]
                        Weather to compute KDE for single diagnostics (in addition to the 
                        KDE for all the diagnostics together).
        KDE_pickle_loc: string [default: None]
                        If specified, the location to save the pickle files generated for
                        each input points. Useful for subsequent plotting of the output 
                        via pyqz_plots.plot_global_qz()
        decimals: int [default: 5]
                  The number of decimals to print in the final CSV file.
        missing_values: string [default: '$$$']
                        The symbol used to mark missing values in the CSV file.
        suffix_out: string [default = '_ggQZff-out'}
                    The string to add to the input filename to create the output filename.
        verbose: bool [default: True]
                 Do you want to read anything ?    
        nproc: int [default = 1]
               Defines how many process to use. Set to -1 for using as many as possible.                  
    '''

    # 0) Do some quick tests to avoid crashes later on ...
    if not(os.path.isfile(fn)):
        raise Exception('File unknown: %s' % fn)
    
    if not(fn[-3:] == 'csv'):
        warnings.warn('File extension unknown (%s).' % fn[-3:] +
                    ' Assuming CSV formatted content.')
    
    if suffix_out =='':
        raise Exception("suffix_out='' would overwrite the input file. Specify something else.")
 
    # 1) Get the different ratios and grids names
    data_file = open(fn, 'r')
    rats = data_file.readline()
    data_file.close()
    
    # Clean the end of the line
    if rats[-1:]=='\n' :
        rats=rats[:-1]
    
    rats = rats.split(',') # Assumes CSV !    
    data_range = range(len(rats))
    data_range.pop(rats.index('Id'))
    data = np.genfromtxt(fn,skip_header = 1, missing_values = missing_values, 
                         filling_values = np.nan,
                         usecols = data_range,
                         delimiter = ',',
                         comments='#') 
                         
    ids = np.genfromtxt(fn, skip_header = 1, missing_values = '',
                        filling_values = '',dtype='S15',
                        usecols = (rats.index('Id')),
                        delimiter=',',comments='#')
                        
    # Need to be careful if there's only one spectrum in the file ...
    if len(np.shape(data)) == 1:
        data = data.reshape(1,np.shape(data)[0])
        ids = ids.reshape(1)      
        
    # 2) Alright, ready to launch the machine !
    [P, r] = get_global_qz(data, [rat for rat in rats if rat!='Id'],
                            which_grids=which_grids,qzs=qzs,
                            ids = ids,
                            Pk = Pk,
                            kappa = kappa,
                            struct = struct, 
                            sampling=sampling,
                            error_pdf = error_pdf,
                            srs = srs,
                            flag_level = flag_level,
                            KDE_method = KDE_method,
                            KDE_qz_sampling = KDE_qz_sampling,
                            KDE_do_singles = KDE_do_singles,
                            KDE_pickle_loc = KDE_pickle_loc,
                            nproc = nproc,
                            )    

    
    # 3) And now, save this as another CSV file
    out_file = open(fn[:-4]+suffix_out+'.csv','w')
    # First, the header
    line = ''
    if 'Id' in rats:
        line = 'Id,'
    for item in r:
        line+=item+','
    line = line[:-1]+'\n' 
    out_file.write(line)   
     
    # Then, the data itself
    for i in range(len(P)):
        line = ''
        if 'Id' in rats:
            line = ids[i]+','
        for (j,item) in enumerate(P[i]):
            if r[j] =='flag':
                line +=  np.str(np.int(item))+','
            else:
                line +=  np.str(np.round(item,decimals))+','
        line = line[:-1]+'\n'
        out_file.write(line)      
    out_file.close()    
    
    return P,r
# ----------------------------------------------------------------------------------------