# -*- coding: utf-8 -*-
#
# This program returns the values of log(q) and 12+log(O/H) following the 
# latest MAPPINGS V simulations, and different kappa values.
#
# See the documentation for installation, changelog and usage instructions.
#
# TO DO:
# - account for the different areas of different diagnostics in KDE (hard!)
# - make the allowed diagnostic work for all line setups (e.g. [SII] and [SII]+)
# - speed up interp_qz by doing the whole grid at once. Ok with method='linear' ?
#
# If you find this code useful, please cite the corresponding paper
#
# Dopita et al., ApJ (2013).
#
# Copyright 2014-2015 Frédéric Vogt (frederic.vogt -at- alumni.anu.edu.au)
#
# This file is part of the pyqz Python module.
#
#   The pyqz Python module is free software: you can redistribute it and/or 
#   modify it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, version 3 of the License.
#
#   The pyqz Python module is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License along 
#   with the pyqz Python module.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------------


# ------------------- Import the required modules ------------------------------
print ' '

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
if ncpu_max > 1:
    print ' Multiprocessing possible with: %s cpus max. ' % ncpu_max

# For the Kernel Density Estimator (don't force it if it's not there, but issue
# a Warning)
try:
    import statsmodels.api as sm 
except:
    warnings.warn("Statsmodels module not found. KDE_method must be set to "+
                  "'gauss' or else I will crash.")
                    
# For plotting
import matplotlib
from matplotlib import pyplot as plt

# For now, force a non-interactive backend. Until a better solution is found.
# (Problem: can't open many matplotlib windows with multiple processes at 
#  the same time.)
plt.switch_backend('Agg')
print ' Loaded matplotlib with backend: %s' % matplotlib.get_backend()

from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar

# For generic things
import os
from datetime import datetime as dt
import time

# Get some important metadata, the code version and nifty tools ...
from pyqz_metadata import *
from pyqz_metadata import __version__
from pyqz_tools import *
print ' Loaded pyqz %s' % __version__
print ' '

# ------------------- Make the plots look good ---------------------------------
import matplotlib as mpl
# Use mathtext & Helvetica: avoid usetex for portability's sake ...
mpl.rc('font',**{'family':'sans-serif', 'serif':['Bitstream Vera Serif'], 
                 'sans-serif':['Helvetica'], 'size':20, 
                 'weight':'normal'})
mpl.rc('axes',**{'labelweight':'normal', 'linewidth':1})
mpl.rc('ytick',**{'major.pad':8, 'color':'k'})
mpl.rc('xtick',**{'major.pad':8, 'color':'k'})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm', 
                     'bf':'monospace:bold'})
mpl.rc('text', **{'usetex':False})
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ------------------- And now for the main pyqz functions ----------------------
# ------------------------------------------------------------------------------

# A function to get the reference grid - useful to make plots !
def get_grid(ratios, coeffs=[[1,0],[0,1]],
                Pk = 5.0, kappa=np.inf, struct='sph', sampling = 1):
    ''' 
        Returns a given line ratio diagnostic grid generated using MAPPINGS
        for a given kappa value.
        
        :param ratios: {string} 
                        The line ratios involved in the grid, e.g.
                        '[NII]/[SII]+;[OIII]/Hb'
        :param coeffs: {list; length = 2, default = [[1,0],[0,1]]} 
                        The different coefficients with which to 
                        mix the different line ratios. The size of each list 
                        elements must be equal to the number of line ratios 
                        involved.
        :param Pk: {float;default = 5.0} 
                    MAPPINGS model pressure. 
                    Value must match an existing reference grid file !
        :param kappa: {float; default = np.inf} 
                        The kappa value.
                        Value must match an existing reference grid file !
        :param struct: {string; default = 'sph'}
                        spherical ('sph') or plane-parallel ('pp') HII regions.
                        Value must match an existing reference grid file !
        :param sampling: {int; default = 1}
                        Use a resampled grid ?
        
        :returns: [grid,grid_cols,metadata]
                  grid {numpy array}: the diagnostic grid as a numpy array.
                  grid_cols {list}: labels of each column inside grid.
                  metadata {list}: basic info from the MAPPINGS simulations.
    '''
     
    # 0) Do some basic tests on the input
    if kappa == 'inf':
        kappa = np.inf
    #if sampling > 1 :
    #    resamp = '_samp_'+np.str(sampling)
    #else:
    #    resamp = ''
  
    # 1) Get the metadata about the file
    fn = get_MVphotogrid_fn(Pk = Pk, kappa = kappa, 
                            struct = struct, sampling = sampling)
    metadata = get_MVphotogrid_metadata(fn)

    # Extract the individual line ratios
    ratios = ratios.split(';')   
            
    # Check the line ratio name, just in case ...
    for ratio in ratios:
        if not(ratio in metadata['columns']):
            sys.exit('Line ratio unknown: %s' % ratio) 
    # Also check that the coeffs size matches the number of line ratios
    if len(coeffs) !=2 or len(coeffs[0]) != len(ratios) or \
                                            len(coeffs[1]) != len(ratios) :
        sys.exit('Mixing coefficient error (size mismatch): %s' % coeffs)        

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

# This function is designed to inspect a given line ratio grid
def check_grid(ratios, coeffs = [[1,0],[0,1]],
                Pk = 5.0, kappa=np.inf, struct='sph', sampling=1,
                color_mode = 'Tot[O]+12', show_plot=False,
                nfig=False, save_plot=False):
    '''
        Creates the diagram of a given line ratio grid for rapid inspection.
        
        :param ratios: {string} 
                        The line ratios involved in the grid, e.g.
                        '[NII]/[SII]+;[OIII]/Hb'
        :param coeffs: {list; length = 2, default = [[1,0],[0,1]]} 
                        The different coefficients with which to 
                        mix the different line ratios. The size of each list 
                        elements must be equal to the number of line ratios 
                        involved.
        :param Pk: {float;default = 5.0} 
                    MAPPINGS model pressure. 
                    Value must match an existing reference grid file !
        :param kappa: {float; default = np.inf} 
                        The kappa value.
                        Value must match an existing reference grid file !
        :param struct: {string; default = 'sph'}
                        spherical ('sph') or plane-parallel ('pp') HII regions.
                        Value must match an existing reference grid file !
        :param sampling: {int; default = 1}
                            Use a resampled grid ?
        :param color_mode: {string; default = 'Tot[O]+12'}
                        'Tot[O]+12'(default) or 'gas[O]+12' or 'LogQ' to 
                        color-code the grid  
        :param show_plot: {bool; default = False}
                        whether to show the plot or not      
        :param nfig: {int; default = False}
                        defines the plot window number
        :param save_plot: {string; default = False}
                        'path+name+format' to export the Figure to.
                        
        :returns: a list of bad segments with the node coordinates
   ''' 

    # 0) Let's get the data in question
    [grid, grid_cols, metadata] = get_grid(ratios, 
                                            coeffs=coeffs, Pk=Pk, kappa=kappa,     
                                            struct=struct, sampling = sampling)
                                            
    # 1) Start the plotting
    if show_plot or save_plot:
        if nfig>0:
            plt.close(nfig)
            fig = plt.figure(nfig,figsize=(10,8))
        else:
            fig = plt.figure(figsize=(10,8))
    
        gs = gridspec.GridSpec(1,2, height_ratios=[1], width_ratios=[1,0.05])
        gs.update(left=0.14,right=0.88,bottom=0.14,top=0.92,
                    wspace=0.1,hspace=0.1)    
        ax1 = fig.add_subplot(gs[0,0])  

    # 1-1) Get some useful info for the plot
    if not(color_mode in grid_cols):
        sys.exit('color_mode unknown: %s' % color_mode)

    # 2) Plot the grid points 
    if show_plot or save_plot:
        # Let's make the distinction between the 'TRUE' MAPPINGS point, 
        # and those that were interpolated in a finer grid using Akima splines
        if sampling > 1:
            grid_pts = ax1.scatter(grid[:,grid_cols.index('Mix_x')],
                            grid[:,grid_cols.index('Mix_y')],
                            marker='o',
                            c=grid[:,grid_cols.index(color_mode)],
                            s=30, cmap=pyqz_cmap_0, edgecolor='none', 
                            vmin=np.min(grid[:,grid_cols.index(color_mode)]), 
                            vmax=np.max(grid[:,grid_cols.index(color_mode)]), 
                            zorder=3)
            true_points = ax1.scatter(grid[ [n for n in range(len(grid)) 
                                             if (grid[n,metadata['columns'].index('LogQ')] in metadata['resampled']['LogQ']) and (grid[n,metadata['columns'].index('Tot[O]+12')] in metadata['resampled']['Tot[O]+12'])],grid_cols.index('Mix_x')],
                                     grid[ [n for n in range(len(grid)) 
                                             if (grid[n,metadata['columns'].index('LogQ')] in metadata['resampled']['LogQ']) and (grid[n,metadata['columns'].index('Tot[O]+12')] in metadata['resampled']['Tot[O]+12'])],grid_cols.index('Mix_y')],
                                    marker='o',
                                    c=grid[ [n for n in range(len(grid)) 
                                             if (grid[n,metadata['columns'].index('LogQ')] in metadata['resampled']['LogQ']) and (grid[n,metadata['columns'].index('Tot[O]+12')] in metadata['resampled']['Tot[O]+12'])],grid_cols.index(color_mode)],
                                    s=60, cmap=pyqz_cmap_0, edgecolor='k', facecolor='white',
                                    vmin=np.min(grid[:,grid_cols.index(color_mode)]), 
                                    vmax=np.max(grid[:,grid_cols.index(color_mode)]), 
                                    zorder=5)          
                            
        else:
            grid_pts = ax1.scatter(grid[:,grid_cols.index('Mix_x')],
                            grid[:,grid_cols.index('Mix_y')],
                            marker='o',
                            c=grid[:,grid_cols.index(color_mode)],
                            s=60, cmap=pyqz_cmap_0, edgecolor='k', 
                            vmin=np.min(grid[:,grid_cols.index(color_mode)]), 
                            vmax=np.max(grid[:,grid_cols.index(color_mode)]), 
                            zorder=3)
            

    # 2-1) Draw the grid lines, and check for problematic regions
    # Check as a function of LogQ and Tot[O]+12. Where are these ?
    u = grid_cols.index('LogQ')
    v = grid_cols.index('Tot[O]+12')
    grid_segs = []
    bad_segs = []

    for i in [u,v]:  
        # Here, 'var' plays the role of 'q' or 'z' depending on 'i'.
        for var in np.unique(grid[:,i]):
            # Plot the grid line
            if show_plot or save_plot:
                ax1.plot(grid[grid[:,i]==var][:,grid_cols.index('Mix_x')],
                            grid[grid[:,i]==var][:,grid_cols.index('Mix_y')],
                            'k-', lw = 1, zorder=1)
            
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
                        if seg_intersect(seg_a[0],seg_a[1],seg_b[0],seg_b[1]):
                            # Yes, then store both as bad segments
                            bad_segs.append(seg_a)
                            bad_segs.append(seg_b)
                            
                # b) add it to the general list 
                grid_segs.append(seg_a)

    # I know have a list of bad_segments - plot them all in red to show which
    # they are !
    if len(bad_segs) >0:
        # Make sure each segment is present only once 
        # (mind you, it's a list of list, which requires some fancy magic)
        bad_segs.sort()
        bad_segs = list(bad_segs for bad_segs,_ in itertools.groupby(bad_segs))
    
        if show_plot or save_plot:
            for bad_seg in bad_segs:
                ax1.plot([bad_seg[0][0],bad_seg[1][0]],
                            [bad_seg[0][1],bad_seg[1][1]],
                            'r-',linewidth=4, zorder=0)
      
    if show_plot or save_plot:                          
        # 3) Plot the colorbar
        cb_ax = plt.subplot(gs[0,1])
        cb = Colorbar(ax = cb_ax, mappable = grid_pts, orientation='vertical')
        
        # Colorbar legend
        cb.set_label(color_mode, labelpad = 10)
        
        # 4) Axis names, kappa value, etc ...
        labelx = ''
        labely = ''
        rats = ratios.split(';')
        # Make sure the labels look pretty in ALL cases ...
        for n in range(len(rats)): 
            if coeffs[0][n] !=0:
                if not(coeffs[0][n] in [1,-1]):
                    if n !=0:
                        labelx += '%+.03g ' % coeffs[0][n]
                    else:
                        labelx += '%.03g ' % coeffs[0][n] 
                elif coeffs[0][n] == 1 and labelx != '':
                    labelx += '+ '
                elif coeffs[0][n] == -1:
                    labelx += '- ' 
                labelx += rats[n].replace('+','$^{+}$') + ' '
            if coeffs[1][n] != 0:
                if not(coeffs[1][n] in [1,-1]):
                    if n != 0:
                        labely += '%+.03g ' % coeffs[1][n] 
                    else:
                        labely += '%.03g ' % coeffs[1][n]
                elif coeffs[1][n] == 1 and labely != '':
                    labely += '+ '
                elif coeffs[1][n] == -1:
                    labely += '- '
                labely += rats[n].replace('+','$^{+}$')+' '   
                            
                                                  
        ax1.set_xlabel(labelx,labelpad = 10)
        ax1.set_ylabel(labely,labelpad = 10)
        
        if not(kappa in [np.inf, 'inf']) :
            kappa_str = r'$\kappa$ = '+str(kappa)    
        else :
            kappa_str = r'$\kappa$ = $\infty$'
        
        ax1.text(0.85,0.9,kappa_str, 
                    horizontalalignment='left',verticalalignment='bottom',
                    transform=ax1.transAxes)
        ax1.grid(True)
        
        # All done ...
        plt.show()
        if save_plot :
            fig.savefig(save_plot, bbox_inches='tight')
        if show_plot:
            plt.show()
        else:
            plt.close()
            
    return bad_segs
# ------------------------------------------------------------------------------

# The core function - returns 'q' or'z' for a given ratio (and a given grid !) 
# Interpolate case by case for even better results !
def interp_qz ( qz, 
                data,
                ratios,
                coeffs =[[1,0],[0,1]],
                Pk = 5.0, kappa=np.inf, struct='sph',sampling=1,
                method = 'linear', 
                show_plot = False, n_plot = False, save_plot = False,
                verbose = True, 
                ):
    ''' The core function of pyqz.
    
        Returns the 'q' or 'z' value for a given set of line ratios based on a
        given 2-D diagnostic grid.
        
        :param qz: {string} 
                    Which estimate to return, 'LogQ', 'Tot[O]+12' or 'gas[O]+12'  
        :param data: {list of numpy array} 
                        List of Arrays of the line ratio values. One array per line ratio
                        of size NxM for NxM values).
        :param ratios: {string} 
                        The line ratios involved in the grid, e.g.
                        '[NII]/[SII]+;[OIII]/Hb'
        :param coeffs: {list; length = 2, default = [[1,0],[0,1]]} 
                        The different coefficients with which to 
                        mix the different line ratios. The size of each list 
                        elements must be equal to the number of line ratios 
                        involved.
        :param Pk: {float;default = 5.0} 
                    MAPPINGS model pressure. 
                    Value must match an existing reference grid file !
        :param kappa: {float; default = np.inf} 
                        The kappa value.
                        Value must match an existing reference grid file !
        :param struct: {string; default = 'sph'}
                        spherical ('sph') or plane-parallel ('pp') HII regions.
                        Value must match an existing reference grid file !
        :param sampling: {int; default = 1}
                            Use a resampled grid ?
        :param method: {string; default = 'linear'}
                        'linear' or 'cubic' interpolation
                        method of scipy.interpolate.griddata.
                        WARNING: using'cubic' can lead to estimates slightly 
                        "outside" the allowed range!
        :param show_plot: {bool; default = False}
        :param n_plot: {int; default = False}
                        defines the plot window number
        :param save_plot: {string; default = False}
                        'path+name+format' to export the Figure to.
        :param verbose: {bool; default=True}  

        :returns: {numpy array of NxM}   
    '''

    if not(qz in ['LogQ','Tot[O]+12', 'gas[O]+12']) :
        sys.exit('Unknown qz string: %s' % qz)

    for dat in data:
        if type(dat) != np.ndarray:
            sys.exit('Input line ratios must be Numpy arrays and not: %s, %s' % 
                        type(dat))
        if np.shape(dat) != np.shape(data[0]):
            sys.exit('Input data arrays size mismatch !')    
    
    # Save the input shape for later - I'll reshape everything in 
    # 1-dimension for simplicity.
    input_data_shape = np.shape(data[0])
    input_data_size = np.array(input_data_shape).prod()
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

    # 2-2) Start the plot already ... will plot in the loop ...
    if show_plot or save_plot :
        if n_plot>0:
            plt.close(n_plot)
            fig = plt.figure(n_plot, figsize=(10,8))
        else:
            fig = plt.figure(figsize=(10,8))

        gs = gridspec.GridSpec(1,2, height_ratios=[1], width_ratios=[1,0.05])
        gs.update(left=0.14,right=0.88,bottom=0.14,top=0.92,
                  wspace=0.1,hspace=0.1)
            
        ax1 = fig.add_subplot(gs[0,0])  
        #xlim1 = diagnostics[ratios]['xliml']
        #xlim2 = diagnostics[ratios]['xlimr']
        #ylim1 = diagnostics[ratios]['ylimb']
        #ylim2 = diagnostics[ratios]['ylimt']

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

            # 2-4aa) Check if this box has any bad segments and if so, skip it !
            nedges = len(this_panel[:,grid_cols.index(ratios.split(';')[0])])
            if nedges != 4:
                sys.exit('Something went terribly wrong ... (L736)')
            
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
            
            # Here, I need to add a check for the presence of data points
            # inside the grid panel. If there's nothing, no need to interpolate 
            # anything !
            path_codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,
                            Path.CLOSEPOLY,]
            # Mind the order of the nodes ... !                                              
            path_verts = this_panel[[0,1,3,2,0],:][:,[rat1_ind,rat2_ind]].tolist()            
            panel_path = Path(path_verts, path_codes)
            
            points_in_panel = panel_path.contains_points([(data_comb[0][k],data_comb[1][k]) for k in range(len(data_comb[0]))],
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

                # 2-4d) Interpolate !
                this_interp_data = interpolate.griddata(this_stretched_panel[:,
                                                                    [rat1_ind,
                                                                    rat2_ind]
                                                                    ],
                                        this_panel[:,var],
                                        (stretched_data[0],stretched_data[1]),
                                        method=method,
                                        fill_value=np.nan)
            
                # 3-4) Store the cleaned interpolated values in the final array
                interp_data[np.where(points_in_panel)[0][np.where(this_interp_data>0)[0]]] = \
                            this_interp_data[this_interp_data>0]

            # 4) Plot for test purposes
            if show_plot or save_plot:

                # plot the panel - use the Convex Hull approach for the fun of it

                panel_hull = scipy.spatial.ConvexHull(this_panel[:,
                                                                [rat1_ind,
                                                                 rat2_ind]])
                                                                 
                for simplex in panel_hull.simplices:
                    ax1.plot(this_panel[:,[rat1_ind,rat2_ind]][simplex,0],
                            this_panel[:,[rat1_ind,rat2_ind]][simplex,1],
                            'k-',lw=0.5,zorder=0)                                                      
            
                # 4-2) Grid points
                if sampling > 1:
                    ref_pts = ax1.scatter(this_panel[:,rat1_ind],
                                            this_panel[:,rat2_ind],
                                            marker='o',
                                            c=this_panel[:,var],
                                            s=40, cmap=pyqz_cmap_0, edgecolor='none', 
                                            vmin=np.min(the_grid[:,var]), 
                                            vmax=np.max(the_grid[:,var]),zorder=1)
                    
                    
                    
                    true_pts = ax1.scatter(this_panel[ [n for n in range(len(this_panel)) 
                                             if (this_panel[n,the_metadata['columns'].index('LogQ')] in the_metadata['resampled']['LogQ']) and (this_panel[n,the_metadata['columns'].index('Tot[O]+12')] in the_metadata['resampled']['Tot[O]+12'])],rat1_ind],
                                            this_panel[ [n for n in range(len(this_panel)) 
                                             if (this_panel[n,the_metadata['columns'].index('LogQ')] in the_metadata['resampled']['LogQ']) and (this_panel[n,the_metadata['columns'].index('Tot[O]+12')] in the_metadata['resampled']['Tot[O]+12'])],rat2_ind],
                                    marker='o',
                                    c=this_panel[ [n for n in range(len(this_panel)) 
                                             if (this_panel[n,the_metadata['columns'].index('LogQ')] in the_metadata['resampled']['LogQ']) and (this_panel[n,the_metadata['columns'].index('Tot[O]+12')] in the_metadata['resampled']['Tot[O]+12'])],var],
                                    s=40, cmap=pyqz_cmap_0, edgecolor='k', facecolor='white',
                                    vmin=np.min(the_grid[:,var]), 
                                    vmax=np.max(the_grid[:,var]), 
                                    zorder=1) 
                                    
                else:
                    ref_pts = ax1.scatter(this_panel[:,rat1_ind],this_panel[:,rat2_ind],
                                            marker='o',
                                            c=this_panel[:,var],
                                            s=40, cmap=pyqz_cmap_0, edgecolor='k', 
                                            vmin=np.min(the_grid[:,var]), 
                                            vmax=np.max(the_grid[:,var]))

    if show_plot or save_plot:
        # 4-3) Interpolated points
        ax1.scatter(data_comb[0],data_comb[1],marker='s', c=interp_data,
                        s=60, cmap = pyqz_cmap_0, edgecolor='none',
                        vmin = np.min(the_grid[:,var]),
                        vmax = np.max(the_grid[:,var]),zorder=2)

                                
                                                                        
        # Plot also the points outside the grid ?
        # Which are they ? Need to check if they have a valid input first !
        # mind the "~" to invert the bools ! isn't this cool ?
        my_out_pts = ~np.isnan(data_comb[0]) * ~np.isnan(data_comb[1]) * \
            np.isnan(interp_data)
  
        
        if not(np.all(data[0].reshape(input_data_shape)==fullgrid_x)):
            ax1.scatter(data_comb[0][my_out_pts], 
                        data_comb[1][my_out_pts],
                                  marker = '^', facecolor = 'none', 
                                  edgecolor = 'k', s=60)
        
        # plot the colorbar - only if the point was found in the first place
        #if not(np.all(np.isnan(interp_data))):
        if True:
            cb_ax = plt.subplot(gs[0,1])
            cb = Colorbar(ax = cb_ax, mappable = ref_pts, orientation='vertical')
            # Colorbar legend
            cb.set_label(qz, labelpad=10)
            
        # Axis names
        labelx = ''
        labely = ''
        rats = ratios.split(';')
        for n in range(len(rats)): 
            if coeffs[0][n] !=0:
                if np.abs(coeffs[0][n]) !=1:
                    labelx += '%+.03g ' % coeffs[0][n] 
                elif coeffs[0][n] == 1 and labelx != '':
                    labelx += '+ '
                elif coeffs[0][n] == -1:
                    labelx += '- '
                labelx += rats[n]+ ' '
                
            if coeffs[1][n] !=0:
                if np.abs(coeffs[1][n]) !=1:
                    labely += '%+.03g ' % coeffs[1][n] 
                elif coeffs[1][n] == 1 and labely != '':
                    labely += '+ '
                elif coeffs[1][n] == -1:
                    labely += '- '
                labely += rats[n]+ ' '  
                                                              
        ax1.set_xlabel(labelx,labelpad = 10)
        ax1.set_ylabel(labely,labelpad = 10)
        
        # Axis limits        
        #ax1.set_xlim((xlim1,xlim2))
        #ax1.set_ylim((ylim1,ylim2))
        #ax1.set_aspect((xlim2-xlim1)/(ylim2-ylim1))
        
        if not(kappa in [np.inf, 'inf']) :
            kappa_str = r'$\kappa$ = '+str(kappa)    
        else :
            kappa_str = r'$\kappa$ = $\infty$'

        ax1.text(0.85,0.9,kappa_str, 
                    horizontalalignment='left',verticalalignment='bottom',
                    transform=ax1.transAxes)    
        ax1.grid(True)
        
        
        
        if show_plot:
            plt.show()
        if save_plot :
            fig.savefig(save_plot, bbox_inches='tight')
        
        # if I just want to save them, then close it to save memory
        if not(show_plot) and save_plot : 
            plt.close()

    # Finish it all, without forgetting to return an array of the same initial
    # dimension !
    return interp_data.reshape(input_data_shape)
# ------------------------------------------------------------------------------

# The function used by the multiprocessing module. Does everything (plotting, 
# KDE, QZ estimates, ...) for a single spectrum.
def get_global_qz_singlespec((j,final_cols,data, data_cols, which_grids,
                  ids,
                  qzs,
                  Pk, kappa, struct,sampling,
                  error_pdf, 
                  srs,
                  flag_level,
                  KDE_method,
                  KDE_qz_sampling,
                  KDE_do_singles,   
                  KDE_save_PDFs,                      
                  show_plot,
                  save_plot, 
                  plot_loc,
                  plot_fmt,
                  verbose)):
    ''' The single-process function associated with get_global_qz().
        Used to spawn as many process as input spectrum sent by the user
        to get_global_qz(). See http://stackoverflow.com/a/3843313 for more
        details.
        
        Note that the data is fed as one big list.
        
        :param j: {int}
                   The process number.
        :param final_cols: {list}
                   The labelling of the columns in the final_data array.           
        :param data: {numpy array of size 1x2*M}  
                        M line fluxes and errors. An error =-1 implies
                        the associated line flux is an upper limit.
        :param data_cols: {list of Mx2 strings}
                        Data content, e.g. ['[NII]+','stdHa',...]
                        Syntax MUST match the MAPPINGS + pyqz conventions !
        :param which_grids: {list of strings}
                            list of the model grids to use for the estimations,
                            e.g. ['[NII]+/Ha;[OII]/Hb',...]
        :param ids: {string; default = False}
                    a string to identify 
                    each data point (e.g. spaxel number, source ID, etc ...)
        :param qzs: {list; default = ['LogQ','Tot[O]+12']}
                    list of Q/Z values to compute
        :param Pk: {float;default = 5.0} 
                    MAPPINGS model pressure. 
                    Value must match an existing reference grid file !
        :param kappa: {float; default = np.inf} 
                        The kappa value.
                        Value must match an existing reference grid file !
        :param struct: {string; default = 'sph'}
                        spherical ('sph') or plane-parallel ('pp') HII regions.
                        Value must match an existing reference grid file !
        :param sampling: {int; default = 1}
                        Use a resampled grid ?
        :param error_pdf: {string; default = 'normal'}
                        The shape of the error function for the line fluxes. 
                        Currently, only 'normal' is supported.
        :param srs: {int; default = 400}
                    The number of random line fluxes generated 
                    to discretize (and reconstruct) the joint probability 
                    function. 
        :param flag_level: {float; default = 2}
                           A 'flag' is raised (in the output file) 
                           when the direct q/z estimate and the KDE q/z
                           estimate (q_rs and z_rs) are offset by more than 
                           flag_level * standard_deviation.
                           Might indicate trouble.
                           A flag is always raised when a point is outside all of 
                           the grids (8), when the KDE PDF is multipeaked (9) or 
                           when srs = 0 (-1, no KDE computed)
        :param KDE_method: {string; default = 'gauss'}
                            Whether to use scipy.stats.gaussian_kde ('gauss') or 
                            sm.nonparametric.KDEMultivariate ('multiv') to 
                            perform the 2-D Kernel Density Estimation.
        :param KDE_qz_sampling: {complex; default= 101j}
                            Sampling of the QZ plane for the KDE reconstruction. 
        :param KDE_do_singles: {bool; default=True}
                             Weather to compute KDE for single diagnostics also.
        :param KDE_save_PDFs: {bool/string; default = False}
                           If not False, the location to save all the PDFs 
                           computed by the code.
        :param show_plot: {string/bool; default = False}
                            Show all plots (True), none (False), only the grids 
                            ('grids') or the KDE maps ('KDE').
        :param save_plot: {string/bool; default = False}
                            Save all plots (True), none (False), only the grids 
                            ('grids'), only the KDE maps ('KDE_all') or just the 
                            KDE maps raising a flag ('KDE_flag').
        :param plot_loc: {string; default = './'}
                            Location where figures are saved.
        :param plot_fmt: {string; default = 'png'}
                            'eps', 'pdf', or 'png'.
        :param verbose: {bool; default=True}   

        :returns: [j,P], where P contains all the estimated Log(Q) and [O]+12
                     values  
        
    '''
    
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
        
    # Is srs = 0 ? Then no KDE !
    if srs == 0:
        noKDE = True
    # Are all the errors 0 ? Then this means no KDE !    
    elif np.all( data[[i for i in range(len(data_cols)) if 'std' in data_cols[i]]]==0):
        noKDE = True
    else:
        noKDE = False
        
    # Generate the random data around this point
    nlines = len(data_cols)/2
    line_names = []
    rsf = np.zeros((srs+1, nlines))
        
    # Do it for all the lines separately <= Lines and errors are 
    # uncorrelated !
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
                    sys.exit('Error: dFlux <0 and dFlux != -1 '+
                                'is not allowed. Check your errors !')      

    # Create a structure where I can store all the real and fake estimates
    # for 'all' the diagnostics in questions.
    discrete_pdf = {}
                    
    # Start filling the 'final_data' array with information I already have
    for (i,calc) in enumerate(final_cols):    
    
        # Individual diagnostics
        if calc.split('|')[0] in diagnostics and not('KDE' in calc) :
            qz = calc.split('|')[1]
            ratios = calc.split('|')[0]
            Rval = []
            for k in range(len(ratios.split(';'))):
                this_ratio = calc.split('|')[0].split(';')[k]
                l1 = this_ratio.split('/')[0]
                l2 = this_ratio.split('/')[1]
            
                Rval.append(np.log10(rsf[:,line_names.index(l1)]/
                            rsf[:,line_names.index(l2)]))
                
            coeffs = diagnostics[calc.split('|')[0]]['coeffs']
            
            # Do I want to save the files ?
            if (save_plot in [True,'grids']):
                # TODO?: give the user the ability to feed the custom file
                # name for EACH point ... maybe via optional function param
                if ids:
                    pname = np.str(ids)+'_'
                else:
                    pname = '' 
                plot_name = os.path.join(plot_loc, pname+\
                            ratios.replace(';','_vs_').replace('/','_')+'_'+qz+'_'+\
                            'Pk'+str(np.int(10*Pk))+'_'+\
                            'k'+str(kappa)+'_'+\
                            struct+\
                            '.'+plot_fmt)
            else:
                plot_name=False
            #st = dt.now()
            # launch the qz interpolation (10% total time)
            qz_values = interp_qz(qz,Rval,ratios, 
                                    coeffs = coeffs, 
                                    Pk=Pk,kappa=kappa,struct=struct,sampling=sampling,
                                    show_plot = (show_plot in [True,'grids']),
                                    save_plot=plot_name)
            #print dt.now()-st
            # Store the 'real' line estimate                
            final_data[i] = qz_values[0]
            # Keep the other random estimates as well
            discrete_pdf[calc] = qz_values[:]
        
    # 4) Calculate mean 'qz' values
    # First, look at the reconstructed joint probability density function

    if ( (show_plot in [True,'KDE']) or (save_plot in [True,'KDE_all',
                                                'KDE_flag'])):
        plt.figure(figsize=(13,8))
        gs = gridspec.GridSpec(2,2, height_ratios=[0.05,1.], 
                                width_ratios=[1.,1.])
        gs.update(left=0.1,right=0.95,bottom=0.1,top=0.9,
                    wspace=0.15,hspace=0.07 )
        ax1 = plt.subplot(gs[1,0])
        ax2 = plt.subplot(gs[1,1])

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
                                                
            # Plot these ... 
            if show_plot or save_plot:
              for ax in [ax1,ax2]: 
                if not(noKDE):
                    ax.scatter(rs_points[qzs[0]],rs_points[qzs[1]],
                                marker='o',edgecolor='k',#'darkorange',
                                facecolor='k',linewidth=0.7, s=2, c='k', 
                                zorder=1)  
                                        
                ax.plot(rs_ref[qzs[0]],rs_ref[qzs[1]],c='w', 
                        marker = 's',markeredgewidth=1.5,
                        markerfacecolor='w', 
                        markeredgecolor='firebrick',markersize=11, zorder=10)
                ax.text(rs_ref[qzs[0]],rs_ref[qzs[1]],
                        which_grids.index(item.split('|')[0])+1,size=12,
                        ha='center',va='center', zorder=11,color='firebrick')

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

    # Now, perform the KDE over the entire set of points
        
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
        QZs_grid[qz] = np.mgrid[QZs_lim[qz][0]:QZs_lim[qz][1]:KDE_qz_sampling]

    # WARNING: HERE, I ASSUME THAT qzs has a length of only 2 !
    # Unfortunate, but that's all you get for now ...
    gridX, gridY = np.mgrid[QZs_lim[qzs[0]][0]:QZs_lim[qzs[0]][1]:KDE_qz_sampling, 
                            QZs_lim[qzs[1]][0]:QZs_lim[qzs[1]][1]:KDE_qz_sampling]
    grid_positions = np.vstack([gridX.ravel(), gridY.ravel()])

    # For speed purposes, we will only reconstruct the KDE over the area where
    # we have some QZ estimates (no need to calculate the PDF over the entire QZ space,
    # it's a waste of time - we know it's almost 0)
    # So, define a subgrid just around all the QZ estimates
    
    # Bug correction (01.2016, F.P.A. Vogt). What if the points is outside 'all'
    # the diagnostics ? Must handle this case as well.
    if not(np.all(np.isnan(all_estimates[qzs[0]]))):
        gridxmin = np.where(gridX[:,0]<all_estimates[qzs[0]].min())[0].max()
        gridxmax = np.where(gridX[:,0]>all_estimates[qzs[0]].max())[0].min()
        gridymin = np.where(gridY[0,:]<all_estimates[qzs[1]].min())[0].max()
        gridymax = np.where(gridY[0,:]>all_estimates[qzs[1]].max())[0].min()
    
    else: # Hum, seem that I have no estimate at all ...
        noKDE = True
        gridxmin = gridX.max()
        gridxmax = gridX.min()
        gridymin = gridY.max()
        gridymax = gridY.min()
        
        
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
    
    #st = dt.now()        
    # ----- This is the slowest bit. Could I speed it up ? How ?  ----------                     
    for kern in all_my_single_kernels.keys():
        gridZ = np.zeros_like(gridX)*np.nan
        if type(all_my_single_kernels[kern]) == type(np.nan):
            pass
        elif KDE_method == 'gauss':
            if do_full_KDE_reconstruction:
                gridZ = np.reshape(all_my_single_kernels[kern](grid_positions), 
                                    gridX.shape)
            else:
                gridZ[gridxmin:gridxmax,gridymin:gridymax] = \
                np.reshape(all_my_single_kernels[kern](subgrid_positions),
                            subgridX.shape)
        
        
        elif KDE_method == 'multiv':
            if do_full_KDE_reconstruction:
                gridZ = np.reshape(all_my_single_kernels[kern].pdf(grid_positions), 
                                    gridX.shape)  
            else:                    
                gridZ[gridxmin:gridxmax,gridymin:gridymax] = \
                np.reshape(all_my_single_kernels[kern].pdf(subgrid_positions),
                            subgridX.shape)                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                               
        else:
           sys.exit('Something went terribly wrong ...(L1402)')                     
            
        # And store the grid for later
        all_my_gridZ[kern] = gridZ
    # ----------------------------------------------------------------------        
    #print dt.now()-st
          
    if KDE_save_PDFs:
        # I want to save the different PDF arrays ... how ?
        # pickle would be easy ... I guess I should allow some
        # CSV files as well eventually ... but not before I assess the
        # need for these.
        if ids:
            d_fn = np.str(ids) + '_KDE-sampled-PDF_' + \
                   np.str(KDE_qz_sampling) + '.pkl'
  
        else:    
            d_fn = np.str(j) + '_KDE-sampled-PDF_' + \
                   np.str(KDE_qz_sampling) + '.pkl'
                
        dump_fn = os.path.join(KDE_save_PDFs,d_fn)
        dump_file = open(dump_fn, 'w')
        pickle.dump([gridX,gridY,all_my_gridZ], dump_file)                    
                     
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

        # Plot it
        if not(noKDE):
            if show_plot or save_plot:
                for ax in [ax1,ax2]:
                    if key == 'all':
                        my_c = 'darkorange'
                        my_c2 = 'none'
                        my_lw = 2.5
                        my_ls = '-'
                        my_zo = 7
                        kde = ax.imshow(gridZ,
                                        extent=[QZs_lim[qzs[0]][0],
                                                QZs_lim[qzs[0]][1],
                                                QZs_lim[qzs[1]][0],
                                                QZs_lim[qzs[1]][1],
                                                ], 
                                        cmap=pyqz_cmap_1,#'GnBu_r', 
                                        origin='lower', 
                                        zorder=0, interpolation='nearest')
                        my_zo = 6         
                                        
                    else:
                        my_c = 'lightblue'
                        my_c2 = 'none'
                        my_lw = 1.5
                        my_zo = 3
                        my_ls = '-'
   
                    kde_cont = ax.contour(QZs_grid[qzs[0]],QZs_grid[qzs[1]],
                                            gridZ,[PDF_cont_level], 
                                            colors='k', 
                                            linestyles=my_ls,
                                            linewidths=my_lw+1, 
                                            zorder=my_zo)
   
                    # 0.61 ~ 1-sigma value of a normalized normal distribution
                    kde_cont = ax.contour(QZs_grid[qzs[0]],QZs_grid[qzs[1]],
                                            gridZ,[PDF_cont_level], 
                                            colors=my_c, 
                                            linestyles=my_ls,
                                            linewidths=my_lw, 
                                            zorder=my_zo)
                                               
                                                
                    # Not showing this anymore ... 
                    # Avoid crowding the plots for nothing ...
                    #
                    #ax.plot(QZs_grid[qzs[0]][peak_loc[1]],
                    #        QZs_grid[qzs[1]][peak_loc[0]],
                    #        c=my_c, marker='^',
                    #        markeredgecolor=my_c, 
                    #        markerfacecolor='none', 
                    #        markeredgewidth = 2,markersize=10, zorder=5)
                    #
            else:                   
                # Ok, the user doesn't want a plot, but I still need one to 
                # get the contours !
                plt.figure()
                kde = plt.imshow(gridZ,
                                extent=[QZs_lim[qzs[0]][0],
                                        QZs_lim[qzs[0]][1],
                                        QZs_lim[qzs[1]][0],
                                        QZs_lim[qzs[1]][1],
                                        ],
                                cmap='gist_yarg', origin='lower', zorder=0, 
                                interpolation='nearest')
              
                # 0.61 ~ 1-sigma value of a normalized normal distribution
                kde_cont = plt.contour(QZs_grid[qzs[0]],QZs_grid[qzs[1]],
                                        gridZ,[PDF_cont_level], 
                                        colors='darkorange', linestyles='-',
                                        linewidths=2, zorder=2)
                plt.close()

        # Get the 'mean location of the 0.61 level contour
        # First, fetch the contour as path; only if contour actually exists
        if np.any(all_my_gridZ[key] == all_my_gridZ[key]):
            path = kde_cont.collections[0].get_paths()[0]
            # If the path/contour is broken, use the one that contains 
            # the peak of the distribution
             
            for path in kde_cont.collections[0].get_paths():
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
            if key == 'all' and len(kde_cont.collections[0].get_paths())>1:
                flag = '9'
            
        else:
            path = [np.nan,np.nan]
            vert = [np.nan,np.nan]
            mean_vert = [np.nan,np.nan]
            err_vert = [np.nan,np.nan]                
            
            
        if show_plot or save_plot:
            # Plot it
            if not(np.all(np.isnan(mean_vert))):
                for ax in [ax1,ax2]:
                    ax.errorbar(mean_vert[0],mean_vert[1],
                                xerr=err_vert[0],yerr=err_vert[1],
                                elinewidth=4., ecolor='k', capthick=4., 
                                zorder=my_zo)  
                    ax.plot(np.mean(vert[:,0]),np.mean(vert[:,1]),c='k', 
                                marker='o', markersize=10,
                                markeredgecolor='k',
                                markerfacecolor='none',
                                markeredgewidth=4, zorder=my_zo)
                                    
                    ax.errorbar(mean_vert[0],mean_vert[1],
                                xerr=err_vert[0],yerr=err_vert[1],
                                elinewidth=2., ecolor=my_c, capthick=2., 
                                zorder=my_zo)  
                    ax.plot(np.mean(vert[:,0]),np.mean(vert[:,1]),c=my_c, 
                                marker='o', markersize=10,
                                markeredgecolor=my_c,
                                markerfacecolor=my_c2,
                                markeredgewidth=2, zorder=my_zo)             
 
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
    # Be CAREFUL here: stats.nanstd is by default the 'unbiased' one 
    # (i.e. 1/(n-1)) ! Use instead 1/n for consistency (= np.std)
    basic_qz = {}

    for qz in qzs:
        
        mean_qz = stats.stats.nanmean(final_data[np.array(which_qz[qz])])
        std_qz = stats.stats.nanstd(final_data[np.array(which_qz[qz])], 
                                    bias=True)
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

    # Finalize the plot ...
    if show_plot or save_plot:  
        for ax in [ax1,ax2]:
            ax.errorbar(basic_qz[qzs[0]][0],basic_qz[qzs[1]][0],
                        xerr=basic_qz[qzs[0]][1],yerr=basic_qz[qzs[1]][1], 
                        elinewidth=4, 
                        ecolor='w',capthick=4,zorder=8)
            ax.errorbar(basic_qz[qzs[0]][0],basic_qz[qzs[1]][0],
                        xerr=basic_qz[qzs[0]][1],yerr=basic_qz[qzs[1]][1], 
                        elinewidth=2, 
                        ecolor='firebrick',capthick=2,zorder=8)
            ax.plot(basic_qz[qzs[0]][0],basic_qz[qzs[1]][0],    
                    '*',c='w', markeredgecolor='w',
                    markeredgewidth=4,markerfacecolor='none',
                    markersize=20, zorder=9)
            ax.plot(basic_qz[qzs[0]][0],basic_qz[qzs[1]][0],    
                    '*',c='w', markeredgecolor='firebrick',
                    markeredgewidth=2,markerfacecolor='w',
                    markersize=20, zorder=9)
                        
            ax.set_xlabel(qzs[0])
            # Temprorarily:
            #ax.set_xlabel(r'$\log(Q)$')
                
            ax.grid(True) 
            
        # Set the left plot to cover the full extent of the qz plane
        ax1.set_xlim(QZs_lim[qzs[0]])
        ax1.set_ylim(QZs_lim[qzs[1]])
        ax1.set_ylabel(qzs[1])
        # Temprorarily:
        #ax1.set_ylabel(r'12+$\log(O/H)$')

        # Define the limits for the right plot - a zoomed-in version
        if len(all_estimates[qzs[0]])>1:
            xlims2 = [np.max([np.min(all_estimates[qzs[0]])-0.05,
                            QZs_lim[qzs[0]][0]]),
                        np.min([np.max(all_estimates[qzs[0]])+0.05,
                            QZs_lim[qzs[0]][1]])]
            ylims2 = [np.max([np.min(all_estimates[qzs[1]])-0.05,
                            QZs_lim[qzs[1]][0]]),
                        np.min([np.max(all_estimates[qzs[1]])+0.05,
                            QZs_lim[qzs[1]][1]])]            
        else:
            xlims2 = [np.max([basic_qz[qzs[0]][0]-0.25,
                            QZs_lim[qzs[0]][0]]),
                      np.min([basic_qz[qzs[0]][0]+0.25,
                            QZs_lim[qzs[0]][1]])]
            ylims2 = [np.max([basic_qz[qzs[1]][0]-0.25,
                            QZs_lim[qzs[1]][0]]),
                      np.min([basic_qz[qzs[1]][0]+0.25,  
                            QZs_lim[qzs[1]][1]])]

        ax2.set_xlim(xlims2)
        ax2.set_ylim(ylims2)
        
        # Plot the window of the right plot in the left one
        rectx = [xlims2[0],xlims2[1],xlims2[1],xlims2[0],xlims2[0]]
        recty = [ylims2[0],ylims2[0],ylims2[1], ylims2[1], ylims2[0]]
        ax1.plot(rectx,recty, 'k--', lw = 2, markersize=5,zorder=1)
        #ax2.plot(rectx,recty, 'k--', lw = 2, markersize=5,zorder=1)
            
        ax2.set_aspect('auto')
        ax1.set_aspect('auto')
        ax1.locator_params(axis='x',nbins=5)
        ax2.locator_params(axis='x',nbins=5)

        # Plot the colorbar
        if not(noKDE):
            cb_ax = plt.subplot(gs[0,:])
            cb = Colorbar(ax = cb_ax, mappable = kde, orientation='horizontal')
            # Colorbar legend
            cb.set_label(r'Joint Probability Density (normalized to peak)', 
                        labelpad = -75)
            cb.ax.xaxis.set_ticks_position('top')
            cb.solids.set_edgecolor('face')
            # Draw the 1-sigma level (assuming gaussian = 61% of the peak)
            cb.ax.axvline(x=PDF_cont_level, color = 'darkorange', linewidth = 3, 
                            linestyle = '-')

        if (save_plot in [True, 'KDE_all']) or ((save_plot == 'KDE_flag') and 
                                              (np.int(flag) !=0)):
                
            if ids:
                pname = np.str(ids)+'_'
            else:    
                pname= np.str(j)+'_' 

            plot_name = os.path.join(plot_loc, pname+\
                            qzs[0]+'_'+\
                            qzs[1]+'_'+\
                            KDE_method+'_'+\
                            'srs'+np.str(np.int(srs))+'_'+\
                            'Pk'+str(np.int(10*Pk))+'_'+\
                            'k'+str(kappa)+'_'+\
                            struct+\
                            '.'+plot_fmt)
                
            plt.savefig(plot_name, bbox_inches='tight')            
             
        if '8' in np.str(final_data[final_cols.index('flag')]):
            plt.close()
            if verbose:
                print '  %s: No KDE calculable (bad points/grids ?)' %ids                               
        elif show_plot in [True, 'KDE']:
            plt.show()   
        else:
            plt.close()
     
    return [j,final_data]
# ------------------------------------------------------------------------------
    
# ----------------- not used anymore ... but kept for possible future use ------    
# Some small function to get the multiprocessing working fine.
# See http://stackoverflow.com/a/3843313
def get_global_qz_singlespec_init(q):
    
    ''' Some function required to get the multiprocessing working fine.
        See http://stackoverflow.com/a/3843313 for details.
    '''
    
    get_global_qz_singlespec.q = q
# ------------------------------------------------------------------------------



# Get the "global" Q/Z values from a set of line ratios and diagnostics
def get_global_qz(data, data_cols, which_grids,
                  ids = False,
                  qzs = ['LogQ','Tot[O]+12'],
                  Pk = 5.0, kappa=np.inf, struct='sph',sampling=1,
                  error_pdf = 'normal', 
                  srs = 400,
                  flag_level=2., # any float >0: flag_level*std = level above 
                               # which a mismatch between q & z and q_rs & z_rs 
                               # is flagged 
                  KDE_method='gauss', # Which KDE routine to use:
                                           # 'gauss'=gaussian_kde from Scipy: 
                                           # fast but cannot define 2D bandwidth
                                           # 'multiv' = 
                                           # statsmodel 'KDEMultivariate':
                                           # 100x slower but can fine-tune 2-D 
                                           # bandwidth
                  KDE_qz_sampling=101j,
                  KDE_do_singles = True,   
                  KDE_save_PDFs = False,                      
                  show_plot = False,
                  save_plot = False, 
                  plot_loc = './',
                  plot_fmt = 'png',
                  verbose = True,
                  nproc = 1,
                  ):

    ''' Derives the LogQ and [O]+12 values for a given set of line 
        fluxes based on a given set of diagnostic grids.
        
        This function only uses "valid" grids defined in pyqz_metadata.py.
          
        For each set of line ratios, this function combines the estimates from 
        invdividual diagnostics into a combined best estimate of LogQ and 
        [O]+12. Observational errors are propagated via the joint
        probability density function, reconstructed via a Kernel Density
        Estimation.
        
        :param data: {numpy array of size Nx2*M}  
                        N sets of M line fluxes and errors. An error =-1 implies
                        the associated line flux is an upper limit.
        :param data_cols: {list of Mx2 strings}
                        Data content, e.g. ['[NII]+','stdHa',...]
                        Syntax MUST match the MAPPINGS + pyqz conventions !
        :param which_grids: {list of strings}
                            list of the model grids to use for the estimations,
                            e.g. ['[NII]+/Ha;[OII]/Hb',...]
        :param ids: {list; default = False}
                    an optional length-N list of numbers or string to identify 
                    each data point (e.g. spaxel number, source ID, etc ...)
        :param qzs: {list; default = ['LogQ','Tot[O]+12']}
                    list of Q/Z values to compute
        :param Pk: {float;default = 5.0} 
                    MAPPINGS model pressure. 
                    Value must match an existing reference grid file !
        :param kappa: {float; default = np.inf} 
                        The kappa value.
                        Value must match an existing reference grid file !
        :param struct: {string; default = 'sph'}
                        spherical ('sph') or plane-parallel ('pp') HII regions.
                        Value must match an existing reference grid file !
        :param sampling: {int; default = 1}
                        Use a resampled grid ?
        :param error_pdf: {string; default = 'normal'}
                        The shape of the error function for the line fluxes. 
                        Currently, only 'normal' is supported.
        :param srs: {int; default = 400}
                    The number of random line fluxes generated 
                    to discretize (and reconstruct) the joint probability 
                    function. 
        :param flag_level: {float; default = 2}
                           A 'flag' is raised (in the output file) 
                           when the direct q/z estimate and the KDE q/z
                           estimate (q_rs and z_rs) are offset by more than 
                           flag_level * standard_deviation.
                           Might indicate trouble.
                           A flag is always raised when a point is outside all of 
                           the grids (8), when the KDE PDF is multipeaked (9) or 
                           when srs = 0 (-1, no KDE computed)
        :param KDE_method: {string; default = 'gauss'}
                            Whether to use scipy.stats.gaussian_kde ('gauss') or 
                            sm.nonparametric.KDEMultivariate ('multiv') to 
                            perform the 2-D Kernel Density Estimation.
        :param KDE_qz_sampling: {complex; default= 101j}
                            Sampling of the QZ plane for the KDE reconstruction. 
        :param KDE_do_singles: {bool; default=True}
                             Weather to compute KDE for single diagnostics also.
        :param KDE_save_PDFs: {bool/string; default = False}
                           If not False, the location to save all the PDFs 
                           computed by the code.
        :param show_plot: {string/bool; default = False}
                            Show all plots (True), none (False), only the grids 
                            ('grids') or the KDE maps ('KDE').
        :param save_plot: {string/bool; default = False}
                            Save all plots (True), none (False), only the grids 
                            ('grids'), only the KDE maps ('KDE_all') or just the 
                            KDE maps raising a flag ('KDE_flag').
        :param plot_loc: {string; default = './'}
                            Location where figures are saved.
        :param plot_fmt: {string; default = 'png'}
                            'eps', 'pdf', or 'png'.
        :param verbose: {bool; default=True}   
        :param nproc: {int; default = 1}
                            Defines how many process to use. Set to None for using 
                            as many as possible. Can't be more than the number 
                            available cpus.

        :returns: [P,r], where P contains all the estimated Log(Q) and [O]+12
                     values, and r contains the columns name   
    '''
    print ' '
    # Let's keep track of time ...
    starttime = dt.now()
    
    # 0) Do some quick tests to avoid crashes later on ...
    try:
        if flag_level<=0:
            sys.exit('flag_level must be >0: %s' % flag_level)
    except:
        sys.exit('flag_level must be >0: %s' % flag_level)

    if not(plot_fmt in ['png','eps','pdf']):
            sys.exit('plot_fmt unknown: %s' % plot_fmt) 
    try:
        if len(qzs) != 2:
            sys.exit('qzs must be of length 2: %s' % qzs)
        if not('LogQ' in qzs):
            sys.exit("qzs must contain 'LogQ': %s" % qzs)
    except:
        sys.exit('qzs unrecognized: %s' % qzs)

    if not(KDE_method in ['gauss','multiv']):
        sys.exit('KDE_method unrecognized: %s' % KDE_method)

    if not(type(which_grids) == list):
        sys.exit('"which_grids" must be a list/tuple, and not: %s' % type(which_grids))
        
    for grid in which_grids:
        if not(grid in diagnostics.keys()): 
            sys.exit('Diagnostic grid unvalid: %s' % grid)
    
    if not(show_plot in [True,False,'grids','KDE']):
        sys.exit("show_plot must be [True,False,'grids','KDE'].") 
    
    if not(save_plot in [True,False,'grids','KDE_all','KDE_flag']):
        sys.exit("show_plot must be [True,False,'grids','KDE_all','KDE_flag'].")             
                            
    if srs == 0:
        warnings.warn('srs set to 0. No KDE will be computed.')        

    # Make sure the input array has a vaild 2D size
    if len(np.shape(data)) ==1:
        sys.exit('Data shape mismatch. '+
                 'Should be (x,y), not "%s". '% np.shape(data) +
                 'Use np.reshape() ?')  
                 
    # This is nasty - if nproc>1, make sure we have a non-interactive matplotlib 
    # backend.
    if nproc != 1 and not(matplotlib.get_backend() in ['Agg','agg']):
                            
        warnings.warn('nproc >1 requires a non-interactive backend for matplotlib. Setting it to "agg" now.') 
                                    
        plt.switch_backend('agg')
        if verbose:
            print 'New matplotlib backend: %s ' % matplotlib.get_backend()
            time.sleep(2) # to give time to the info to get through ... !
        
        
    # And if nproc > 1, turn the plotting off so we save time - since no plot will
    # be displayed anyway. Note: this is only for the display - saving figures still
    # works with multiprocessing !
    if nproc != 1 and show_plot:
        warnings.warn('show_plot= %s, but nproc >1 implies show_plot= False.' % show_plot)

        show_plot = False
    # Finally, in nproc = None, set it to the max number of cpus. Helps against
    # weird warning ?
    if nproc == None:
        nproc = ncpu_max
    
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
    # Loop through each value, generate many random ones, etc ...   
     
    # Let's use the multiprocessing package to speed things up
    # Especially when computing multiple points at once.        
    
    jobs = []
    
    if not(ids):
        ids = np.arange(npoints)+1
    
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
                  KDE_save_PDFs,                      
                  show_plot,
                  save_plot, 
                  plot_loc,
                  plot_fmt,
                  verbose) )                                    
              
    if len(jobs)>0:
        
        if nproc == 1:
            if verbose:
                print '--> Dealing with them one at a time ... be patient now !' 
                
            results = map(get_global_qz_singlespec,jobs)    
            
        else:
            
            if verbose:
                print '--> Launching the multiple processes ... be patient now !' 
            
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
            while (True):
                if (rs.ready()): break
                remaining = rs._number_left
                sys.stdout.write("    ... %i job(s) left ...\r" % (remaining*chunksize) )
                sys.stdout.flush()
                time.sleep(0.5)
            
            # And for safety, make a "join" call ...
            mypool.join()
            
            results = rs._value
            
            sys.stdout.write("    %i job(s) completed.      \n" % len(results))
            sys.stdout.flush()
            
    # All done. Now let's collect the results ...
    for j in range(len(jobs)):      

        if nproc == 1:                                                                                                                                                           
            final_data[j,:] = results[j][1]
        else:
            #final_data[j,:] = results.__dict__['_items'][j][1][1]
            final_data[j,:] = results[j][1]
            
            # Make sure the good spectrum is returned in the good order
            if results[j][0] != j :
                sys.exit('Something went very badly wrong ... (L1740)')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
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
                     struct = 'sph',
                     sampling=1,
                     error_pdf = 'normal', 
                     srs = 400,
                     flag_level=2., # any float >0: flag_level*std = level above 
                                    # which a mismatch between q & z and q_rs & 
                                    # z_rs is flagged 
                     KDE_method='gauss', # Which KDE routine to use:
                                           # 'gauss' from Scipy: fast 
                                           # but cannot define 2-D bandwidth
                                           # statsmodel 'multiv':
                                           # 100x slower but can fine-tune 2-D 
                                           # bandwidth
                     KDE_qz_sampling=101j,
                     KDE_do_singles = True,     
                     KDE_save_PDFs = False,                 
                     decimals=5,
                     missing_values='$$$', # How missing values are marked - will be 
                                            # replaced by nan's
                     suffix_out = '_out',
                     show_plot = False,
                     save_plot = False, 
                     plot_loc = './',
                     plot_fmt = 'png',
                     verbose = True,
                     nproc = 1,
                     ):
    
    '''
    The get_global_qz-'from a file' function. Gets the line ratios from a csv 
    file, and not directly as Python arrays.
    
    Save to file the log(q) and 12+log(O/H) values for a given set of line 
    fluxes based on a given set of diagnostic grids. Requires a specific 
    input file, and can batch process several measurements automatically. 
        
    :param fn: {string}
                The path+filename of the input file, in CSV format    
    :param which_grids: {list of strings}
                            list of the model grids to use for the estimations,
                            e.g. ['[NII]/Ha;[OIII]/Hb',...]
    :param qzs: {list; default = ['LogQ','Tot[O]+12']}
                    list of Q/Z values to compute
    :param Pk: {float;default = 5.0} 
                MAPPINGS model pressure. 
                Value must match an existing reference grid file !
    :param kappa: {float; default = np.inf} 
                    The kappa value.
                    Value must match an existing reference grid file !
    :param struct: {string; default = 'sph'}
                    spherical ('sph') or plane-parallel ('pp') HII regions.
                    Value must match an existing reference grid file !
    :param sampling: {int; default = 1}
                    Use a resampled grid ?
    :param error_pdf: {string; default = 'normal'}
                        The shape of the error function for the line fluxes. 
                        Currently, only 'normal' is supported.
    :param srs: {int; default = 400}
                The number of random line fluxes generated 
                to discretize (and reconstruct) the joint probability 
                function. 
    :param flag_level: {float; default = 2}
                       A 'flag' is raised (in the output file) 
                       when the direct q/z estimate and the KDE q/z
                       estimate (q_rs and z_rs) are offset by more than 
                       flag_level * standard_deviation.
                       Might indicate trouble.
                       A flag is always raised when a point is outside all of 
                       the grids (8), when the KDE PDF is multipeaked (9) or 
                       when srs = 0 (-1, no KDE computed)
    :param KDE_method: {string; default = 'gauss'}
                        Whether to use scipy.stats.gaussian_kde ('gauss') or 
                        sm.nonparametric.KDEMultivariate ('multiv') to 
                        perform the 2-D Kernel Density Estimation
    :param KDE_qz_sampling: {complex; default= 101j}
                        Sampling of the QZ plane for the KDE reconstruction. 
    :param KDE_do_singles: {bool; default=True}
                        Weather to compute KDE for single diagnostics also.
    :param KDE_save_PDFs: {bool/string; default = False}
                        If not False, the location to save all the PDFs 
                        computed by the code.
    :param decimals: {int; default = 5}
                        Number of decimals to print in the final CSV file.
    :param missing_values: {string; default = '$$$'}
                            Symbol used to mark missing values.
    :param suffix_out: {string; default = '_out'}
                        String to add to the input filename.
                        Warning: '' will overwrite the input file !!!
    :param show_plot: {string/bool; default = False}
                        Show all plots (True), none (False), only the grids 
                        ('grids') or the KDE maps ('KDE').
    :param save_plot: {string/bool; default = False}
                        Save all plots (True), none (False), only the grids 
                        ('grids'), only the KDE maps ('KDE_all') or just the 
                        KDE maps raising a flag ('KDE_flag').
    :param plot_loc: {string; default = './'}
                        Location where figures are saved.
    :param plot_fmt: {string; default = 'png'}
                        'eps', 'pdf', or 'png'.                    
    :param verbose: {bool; default=True} 
    :param nproc: {int; default = 1}
                            Defines how many process to use. Set to None for using 
                            as many as possible. Can't be more than the number 
                            available cpus.                    
    '''

    # 0) Do some quick tests to avoid crashes later on ...
    if not(os.path.isfile(fn)):
        sys.exit('File unknown: %s' % fn)
    if not(fn[-3:] == 'csv'):
        warnings.warn('File extension unknown (%s).' % fn[-3:] +
                    ' Assuming CSV formatted content.')

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
                            KDE_save_PDFs = KDE_save_PDFs,
                            show_plot = show_plot,
                            save_plot = save_plot,
                            plot_loc = plot_loc,
                            plot_fmt = plot_fmt,
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
# ------------------------------------------------------------------------------
   
# End of the World as we know it. ----------------------------------------------
# ------------------------------------------------------------------------------