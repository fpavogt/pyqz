# -*- coding: utf-8 -*-
'''
 pyqz: a Python module to derive the ionization parameter and oxygen abundance of HII 
 regions, by comparing their strong emission line ratios with MAPPINGS photoionization 
 models.\n
 Copyright (C) 2014-2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains tools for the pyqz routines to create pretty plots seamlessly.

 Updated April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

import os
     
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar

import numpy as np
import pickle

import pyqz
import pyqz_tools as pyqzt
import pyqz_metadata as pyqzm

# Use tex - a pain to setup, but it looks great when it works !
# If you encounter trouble running the code, simply comment out the lines below.

#mpl.rc('MacOSX')
#mpl.rc('font',**{'family':'sans-serif', 'serif':['Computer Modern Serif'], 
#                 'sans-serif':['Helvetica'], 'size':15, 
#                'weight':500, 'variant':'normal'})
#mpl.rc('axes',**{'labelweight':'normal', 'linewidth':1})
#mpl.rc('ytick',**{'major.pad':10, 'color':'k', 'direction':'in'})
#mpl.rc('xtick',**{'major.pad':10,'color':'k', 'direction':'in'})
#mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
#mpl.rc('text', **{'usetex':True})
#mpl.rc('text.latex',preamble=r'\usepackage{cmbright},\usepackage{relsize},'+\
#                             r'\usepackage{upgreek}, \usepackage{amsmath}')
#mpl.rc('contour', **{'negative_linestyle':'solid'}) # dashed | solid

# Avoid using rcParams not to mess up with the user's environment too much. 
# Try via plt.style instead. This should be easier to nullify after the fact by calling
# e.g. plt.style.use('classic')
#
# Note to self: I tried to use plt.style.context to only change the style for the pyqz 
# plots and then revert it, but all the "live displays" were failing (because they revert 
# to the default style eventually).

plt.style.use(os.path.join(pyqzm.pyqz_dir,'pyqz_plots.mplstyle'))

# ----------------------------------------------------------------------------------------      

def get_plot_labels(rats,coeffs):
    ''' 
    Create the x and y label strings for the line ratio diagrams. Make sure they look
    as pretty as can be !
    
    :Args:
        rats: list of string
              The list of line ratios involved.
        coeffs: list of list
                The mixing coeffificents of the line ratios.
    
    :Returns:
        labelx,labely: string, string
                       The pretty labels
    '''             

    labelx = ''
    labely = ''
    
    # Loop through each line ratio, and update the labels as required.
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
                
    return labelx,labely
# ---------------------------------------------------------------------------------------- 

# This function is designed to plot a given grid for inspection.
def plot_grid(ratios, coeffs = [[1,0],[0,1]],
              Pk = 5.0, kappa=np.inf, struct = 'pp', sampling = 1,
              color_mode = 'Tot[O]+12', figname = None, show_plot = True,
              data = None, interp_data = None):
    '''
    Creates the diagram of a given line ratio grid for rapid inspection, including wraps
    (fold-over regions), and of certain line ratios, if "data" is specified.
        
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
        color_mode: string [default: 'Tot[O]+12']
                    Color the grid according to 'Tot[O]+12', 'gas[O]+12' or 'LogQ'.    
        figname: string [default: None]
                  'path+name+format' to save the Figure to.
        show_plot: bool [default: True]
                  Do you want to display the Figure ?
        data: list of numpy array [default: None]
              List of Arrays of the line ratio values. One array per line ratio.
        interp_data: numpy array [default: 'linear']
                     interpolated line ratios (via pyqz.interp_qz)
   ''' 

    # 0) Let's get the data in question
    [grid, grid_cols, metadata] = pyqz.get_grid(ratios, coeffs=coeffs, Pk=Pk, kappa=kappa,     
                                                struct=struct, sampling = sampling)
                                            
    # 1) What are the bad segments in this grid ?
    bad_segs = pyqz.check_grid(ratios, coeffs=coeffs, Pk = Pk, kappa=kappa, struct=struct, 
                               sampling = sampling)                  
                                     
    # 2) Start the plotting
    fig = plt.figure(figsize=(10,8))
    
    gs = gridspec.GridSpec(1,2, height_ratios=[1], width_ratios=[1,0.05])
    gs.update(left=0.14,right=0.88,bottom=0.14,top=0.92,
                wspace=0.1,hspace=0.1)    
    ax1 = fig.add_subplot(gs[0,0])  

    if not(color_mode in grid_cols):
        raise Exception('color_mode unknown: %s' % color_mode)

    # 2) Plot the grid points 
    # Let's make the distinction between the 'TRUE' MAPPINGS point, 
    # and those that were interpolated in a finer grid using Akima splines
    pts = ax1.scatter(grid[:,grid_cols.index('Mix_x')],
                      grid[:,grid_cols.index('Mix_y')],
                      marker='o',
                      c=grid[:,grid_cols.index(color_mode)],
                      s=30, cmap=pyqzm.pyqz_cmap_0, edgecolor='none', 
                      vmin=np.min(grid[:,grid_cols.index(color_mode)]), 
                      vmax=np.max(grid[:,grid_cols.index(color_mode)]), 
                      zorder=3)
           
    # Now mark the "original" points with a black outline. First, which are they ?
    if sampling > 1:
        origin_cond = [n for n in range(len(grid)) if 
                         (grid[n,metadata['columns'].index('LogQ')] in 
                                                      metadata['resampled']['LogQ']) and 
                         (grid[n,metadata['columns'].index('Tot[O]+12')] in 
                                                      metadata['resampled']['Tot[O]+12'])]
    else:
        origin_cond = [n for n in range(len(grid))]
    
    # Now plot the black outlines.
    original_pts = ax1.scatter(grid[origin_cond, grid_cols.index('Mix_x')],
                               grid[origin_cond, grid_cols.index('Mix_y')],
                               marker='o',
                               c=grid[origin_cond, grid_cols.index(color_mode)],
                               s=60, cmap=pyqzm.pyqz_cmap_0, edgecolor='k', facecolor='white',
                               vmin=np.min(grid[:,grid_cols.index(color_mode)]), 
                               vmax=np.max(grid[:,grid_cols.index(color_mode)]), 
                               zorder=5)          


    # 2-1) Draw the grid lines
    # Check as a function of LogQ and Tot[O]+12. Where are these ?
    u = grid_cols.index('LogQ')
    v = grid_cols.index('Tot[O]+12')

    for i in [u,v]:  
        # Here, 'var' plays the role of 'q' or 'z' depending on 'i'.
        for var in np.unique(grid[:,i]):
            # Plot the grid line
            ax1.plot(grid[grid[:,i]==var][:,grid_cols.index('Mix_x')],
                            grid[grid[:,i]==var][:,grid_cols.index('Mix_y')],
                            'k-', lw = 1, zorder=1)
                            
    # Plot the bad segments
    for bad_seg in bad_segs:
        ax1.plot([bad_seg[0][0],bad_seg[1][0]],[bad_seg[0][1],bad_seg[1][1]],  'r-',
                 linewidth=4, zorder=0)
                            
    # Now, also plot the data, if anything was provided !             
    if not(interp_data is None):
        
        # Compute the combined "observed" ratios
        data_comb = [np.zeros_like(data[0]),np.zeros_like(data[1])]    
        for k in range(len(ratios.split(';'))): 
            data_comb[0] += coeffs[0][k]*data[k]
            data_comb[1] += coeffs[1][k]*data[k]
        
        ax1.scatter(data_comb[0][interp_data == interp_data], 
                    data_comb[1][interp_data == interp_data],
                    marker='s', c=interp_data[interp_data == interp_data],
                    s=15, cmap = pyqzm.pyqz_cmap_0, edgecolor='k',
                    vmin = np.min(grid[:,grid_cols.index(color_mode)]),
                    vmax = np.max(grid[:,grid_cols.index(color_mode)]),
                    zorder=2, alpha=0.35)

        # Plot also the points outside the grid ?
        # Which are they ? Need to check if they have a valid input first !
        # mind the "~" to invert the bools ! isn't this cool ?
        my_out_pts = ~np.isnan(data_comb[0]) * ~np.isnan(data_comb[1]) * \
            np.isnan(interp_data)
  
        # Don't do this if this is a test with fullgrid_x ...
        ax1.scatter(data_comb[0][my_out_pts], data_comb[1][my_out_pts],
                    marker = '^', facecolor = 'none', edgecolor = 'k', s=60)           
                               
    # 3) Plot the colorbar
    cb_ax = plt.subplot(gs[0,1])
    cb = Colorbar(ax = cb_ax, mappable = pts, orientation='vertical')
        
    # Colorbar legend
    cb.set_label(color_mode, labelpad = 10)
        
    # 4) Axis names, kappa value, etc ...
    rats = ratios.split(';')
    # Make sure the labels look pretty in ALL cases ...
    labelx,labely = get_plot_labels(rats,coeffs)
                                                  
    ax1.set_xlabel(labelx,labelpad = 10)
    ax1.set_ylabel(labely,labelpad = 10)
        
    if not(kappa in [np.inf, 'inf']) :
        kappa_str = r'$\kappa$ = '+str(kappa)    
    else :
        kappa_str = r'$\kappa$ = $\infty$'
        
    ax1.text(0.85,0.9,kappa_str, horizontalalignment='left',verticalalignment='bottom',
             transform=ax1.transAxes)
    ax1.grid(True)
        
    if figname :
        fig.savefig(figname, bbox_inches='tight')
    if show_plot:
        plt.show()
    else:
        plt.close()
    
# ------------------------------------------------------------------------------


# The function used by the multiprocessing module. Does everything for a single spectrum.
def plot_global_qz(pickle_fn, do_all_diags = True, show_plots = True, save_loc = None):
    ''' A plotting function designed to plot the pickle output from pyqz.get_global_qz().
    
    :Args:
        pickel_fn: string
                   the path+filename of the pickle filed to process.
        do_all_diags: bool [default: True]
                      Whether to construct the diagnostic grid plots for all individual
                      diagnostics or not.
        show_plots: bool [default: True]
                    Whether to display the plots or not.
        save_loc: string [default: None][
                  If set, the location where the plots will be saved. If not set, no plots
                  will be saved.
    :Notes: 
        You must set KDE_pickle_loc to the location of your choice when running 
        pyqz.get_global_qz() and pyqz.get_global_qz_ff(), if you want the associated 
        pickle file to be generated.
    '''
    # 0) Open the pickle file, and check what we are dealing with ...
    ftmp = open(pickle_fn,'r')
    KDE_frompickle = pickle.load(ftmp)
    ftmp.close()

    # What will I call my plots ?

    ids = KDE_frompickle['ids']
    Pk = KDE_frompickle['Pk'] 
    kappa = KDE_frompickle['kappa']
    struct = KDE_frompickle['struct']
    sampling = KDE_frompickle['sampling']
    KDE_method = KDE_frompickle['KDE_method']
    
    fn_in = pyqzt.get_KDE_picklefn(ids, Pk, kappa, struct, sampling, KDE_method)

    # First things first. Generate the grid plots for all the diagnostics. ---------------
    
    discrete_pdf = KDE_frompickle['discrete_pdf']
    rsf = KDE_frompickle['rsf']
    rsf_keys = KDE_frompickle['rsf_keys']
    
    if do_all_diags:
        for key in discrete_pdf.keys():
    
            # Reconstruct the line flux mix from the raw rsf data
            Rval = []
            for (r,ratio) in enumerate(key.split('|')[0].split(';')):
                l1 = ratio.split('/')[0]
                l2 = ratio.split('/')[1]
            
                Rval.append(np.log10(rsf[:,rsf_keys.index(l1)]/
                                     rsf[:,rsf_keys.index(l2)]))
     
            if save_loc:
                plot_name = os.path.join(save_loc, fn_in[:-4] + '_' + 
                                key.replace('/','-').replace(';','_vs_').replace('|','_')
                                + '.pdf')
            else:
                plot_name = None
            
            plot_grid(key.split('|')[0], 
                      coeffs = pyqzm.diagnostics[key.split('|')[0]]['coeffs'],
                      Pk = Pk, kappa=kappa, struct = struct, sampling = sampling,
                      color_mode = key.split('|')[1], 
                      figname = plot_name, show_plot = show_plots,
                      data = Rval, interp_data = discrete_pdf[key])
                  
    # Let's get started with the KDE plot ------------------------------------------------
    plt.figure(figsize=(13,8))
    gs = gridspec.GridSpec(2,2, height_ratios=[0.05,1.], width_ratios=[1.,1.])
    gs.update(left=0.1, right=0.95, bottom=0.1,top=0.9, wspace=0.15, hspace=0.07 )
    ax1 = plt.subplot(gs[1,0])
    ax2 = plt.subplot(gs[1,1])
    
    # First, the direct estimates --------------------------------------------------------
    qzs = KDE_frompickle['qzs']
    which_grids = KDE_frompickle['which_grids']
    
    # Create an array to keep track of how big my zoom window will be for ax2
    ax2_lims = [np.nan, np.nan, np.nan, np.nan]
    
    for (k,key) in enumerate(discrete_pdf.keys()):
        # Make sure we show each point only once
        if not(qzs[0] in key):
            continue

        # Don't do anything if the point lands outside the plot.
        if np.isnan(discrete_pdf[key][0]):
            continue
        
        this_diag = key.split('|')[0]
            
        for ax in [ax1,ax2]: 
            
            # The individual realizations
            ax.scatter(discrete_pdf[this_diag+'|'+qzs[0]][1:],
                       discrete_pdf[this_diag+'|'+qzs[1]][1:],  
                       c='k', marker = 'o', s = 1, zorder = 2)  
            
            # Update the zoom-window limits if required
            ax2_lims = [np.nanmin(( ax2_lims[0], 
                                np.nanmin(discrete_pdf[this_diag+'|'+qzs[0]][1:])
                              )),
                        np.nanmax(( ax2_lims[1], 
                                np.nanmax(discrete_pdf[this_diag+'|'+qzs[0]][1:])
                              )),
                        np.nanmin(( ax2_lims[2], 
                                np.nanmin(discrete_pdf[this_diag+'|'+qzs[1]][1:])
                              )),
                        np.nanmax(( ax2_lims[3], 
                                np.nanmax(discrete_pdf[this_diag+'|'+qzs[1]][1:])
                              )),
                       ]
                        
            # The direct estimate               
            ax.plot(discrete_pdf[this_diag+'|'+qzs[0]][0],
                    discrete_pdf[this_diag+'|'+qzs[1]][0],
                    c='w', marker = 's',markeredgewidth=1.5,
                    markerfacecolor='w', 
                    markeredgecolor='k', markersize=11, zorder=10)
            
            grid_number  = which_grids.index(this_diag)+1
      
            ax.text(discrete_pdf[this_diag+'|'+qzs[0]][0],
                    discrete_pdf[this_diag+'|'+qzs[1]][0],
                    grid_number,
                    size=12, ha='center',va='center', zorder=11, color='k')
            
            if ax == ax1:
                # And also add the name of the diagnostic
                ax.text(pyqzm.QZs_lim[pyqzm.M_version][qzs[0]][0] +0.05, 
                        pyqzm.QZs_lim[pyqzm.M_version][qzs[1]][1] -0.05*grid_number, 
                        '%i: %s' % (grid_number,this_diag),
                        horizontalalignment = 'left', verticalalignment = 'center')
    
    # Now, plot the global KDE and the associated contours -------------------------------
    all_gridZ = KDE_frompickle['all_gridZ']
    gridX = KDE_frompickle['gridX']
    gridY = KDE_frompickle['gridY']
    
    final_data = KDE_frompickle['final_data']
    final_cols = KDE_frompickle['final_cols']
    
    # Loop through all reconstructed KDEs. Plot all the contours, but only the full KDE.
    kde = False
    for gridZ_key in all_gridZ.keys():
        
        # If it is just some nans, go no further
        if np.all(np.isnan(all_gridZ[gridZ_key])):
            continue
        
        # Careful here: gridZ has a weird arrangement ...
        gridZ = np.rot90(all_gridZ[gridZ_key])[::-1,:]

        for ax in [ax1,ax2]:
            if gridZ_key == 'all':
                my_c = 'darkorange'
                my_c2 = 'none'
                my_lw = 2.5
                my_zo = 7
                
                # Plot the KDE
                kde = ax.imshow(gridZ,
                                extent=[pyqzm.QZs_lim[pyqzm.M_version][qzs[0]][0],
                                        pyqzm.QZs_lim[pyqzm.M_version][qzs[0]][1],
                                        pyqzm.QZs_lim[pyqzm.M_version][qzs[1]][0],
                                        pyqzm.QZs_lim[pyqzm.M_version][qzs[1]][1],
                                        ], 
                                        cmap=pyqzm.pyqz_cmap_1,
                                        origin='lower', 
                                        zorder=0, interpolation='nearest')      
            else:
                my_c = 'royalblue'
                my_c2 = 'none'
                my_lw = 1.5
                my_zo = 3
        
            # Plot the contour (two times for the look)
            # Mind the fact that gridX and gridY refer to the "un-derotated" gridZ !
                
            kde_cont = ax.contour(gridX, gridY, all_gridZ[gridZ_key], [pyqzm.PDF_cont_level], 
                                  colors = my_c, linestyles = '-', linewidths = my_lw, 
                                  zorder = my_zo)
                                            
            # Plot the individual and global KDE estimates -------------------------------    
            
            if gridZ_key == 'all':
                final_keys = ['<'+qzs[0]+'{KDE}>','err('+qzs[0]+'{KDE})',
                              '<'+qzs[1]+'{KDE}>','err('+qzs[1]+'{KDE})']
            else:
                final_keys = [gridZ_key+'|'+qzs[0]+'{KDE}',
                              'err('+gridZ_key+'|'+qzs[0]+'{KDE})',
                              gridZ_key+'|'+qzs[1]+'{KDE}',
                              'err('+gridZ_key+'|'+qzs[1]+'{KDE})']
                
            mean_est = [final_data[final_cols.index(final_keys[0])],
                        final_data[final_cols.index(final_keys[2])]]
            err_est = [final_data[final_cols.index(final_keys[1])],
                       final_data[final_cols.index(final_keys[3])]]
                       
            ax.errorbar(mean_est[0], mean_est[1],
                        xerr = err_est[0],yerr = err_est[1],
                        elinewidth = 2., ecolor = my_c, capthick = 0, 
                        zorder = my_zo)  
            ax.plot(mean_est[0], mean_est[1],
                    c = my_c, marker = 'o', markersize = 10, markeredgecolor = my_c,
                    markerfacecolor = my_c2, markeredgewidth = 2, zorder = my_zo)             
    
        
    # Plot the combined direct estimates -------------------------------------------------
    
    final_keys = ['<'+qzs[0]+'>','std('+qzs[0]+')', '<'+qzs[1]+'>','std('+qzs[1]+')']
    
    mean_est = [final_data[final_cols.index(final_keys[0])],
                final_data[final_cols.index(final_keys[2])]]
    err_est = [final_data[final_cols.index(final_keys[1])],
               final_data[final_cols.index(final_keys[3])]]
    
    # Only if that value actually exists !
    if not(np.any(np.isnan(mean_est))):
    
        for ax in [ax1,ax2]:
          
            ax.errorbar(mean_est[0], mean_est[1],
                        xerr = err_est[0], yerr = err_est[1], 
                        elinewidth = 2, 
                        ecolor = 'k', capthick = 0, zorder = 8)

            ax.plot(mean_est[0],mean_est[1],    
                    '*', c = 'w', markeredgecolor = 'k',
                    markeredgewidth = 2, markerfacecolor = 'w',
                    markersize = 20, zorder = 9)
    
    # Very well, finalize the plot now.
    for ax in [ax1, ax2]:               
        ax.set_xlabel(qzs[0])   
        ax.grid(True) 
      
    # Set the left plot to cover the full extent of the qz plane
    ax1.set_xlim(pyqzm.QZs_lim[pyqzm.M_version][qzs[0]])
    ax1.set_ylim(pyqzm.QZs_lim[pyqzm.M_version][qzs[1]])
    ax1.set_ylabel(qzs[1])
    
    
    # Define the limits for the right plot - a zoomed-in version
    if np.any(np.isnan(ax2_lims)):
        ax2.set_xlim(pyqzm.QZs_lim[pyqzm.M_version][qzs[0]])
        ax2.set_ylim(pyqzm.QZs_lim[pyqzm.M_version][qzs[1]])
    else:
        ax2.set_xlim(ax2_lims[:2])
        ax2.set_ylim(ax2_lims[2:])
        
        # Plot the window of the right plot in the left one
        rectx = [ax2_lims[0], ax2_lims[1], ax2_lims[1], ax2_lims[0], ax2_lims[0]]
        recty = [ax2_lims[2], ax2_lims[2], ax2_lims[3], ax2_lims[3], ax2_lims[2]]
        ax1.plot(rectx,recty, 'k--', lw = 1.5, markersize=5,zorder=1)
       
    ax2.set_aspect('auto')
    ax1.set_aspect('auto')
    ax1.locator_params(axis='x',nbins=5)
    ax2.locator_params(axis='x',nbins=5)

    # Plot the colorbar if required
    if kde:
        cb_ax = plt.subplot(gs[0,:])
        cb = Colorbar(ax = cb_ax, mappable = kde, orientation='horizontal')
        # Colorbar legend
        cb.set_label(r'Joint Probability Density (normalized to peak)', labelpad = -60)
        cb.ax.xaxis.set_ticks_position('top')
        cb.solids.set_edgecolor('face')
        # Draw the 1-sigma level (assuming gaussian = 61% of the peak)
        cb.ax.axvline(x = pyqzm.PDF_cont_level, color = 'darkorange', linewidth = 3, 
                      linestyle = '-')
   
    # Get ready to save this:
    if save_loc:
        plot_name = os.path.join(save_loc, fn_in[:-3]+'pdf')
        plt.savefig(plot_name, bbox_inches='tight')            
    if show_plots:     
        plt.show()
    else:
        plt.close()
    
# ----------------------------------------------------------------------------------------