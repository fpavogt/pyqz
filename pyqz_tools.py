# -*- coding: utf-8 -*-
#
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

import numpy as np
import scipy
import sys
import os
import subprocess # To launch shell scripts from Python the proper way ...
from datetime import datetime as dt

from pyqz_metadata import *
from pyqz_metadata import __version__

# ------------------- And now for the tool functions ---------------------------

# A function to calculate if three points are listed in a clockwise direction or
# not. Adapted from Bryce Boe:
# http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
def ccw(A,B,C):
    '''
        Check whether 3 points are listed in a counter-clockwise series or not.
        Adapted from Bryce Boe
        http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    
        :param A: 2-D array for point 1
        :param B: 2-D array for point 2
        :param C: 2-D array for point 3
        
        :returns: True or False                
    '''

    if not(A in [B,C]) and B!=C:  
        return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])
    else:
        sys.exit('A, B & C must be distinct!')
# ------------------------------------------------------------------------------
        
def seg_intersect(A,B,C,D):
    '''
        Check whether 2 segments defined by 4 points intersect.
        Adapted from Bryce Boe
        http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    
        :param A: 2-D array for point 1
        :param B: 2-D array for point 2
        :param C: 2-D array for point 3
        :param D: 2-D array for point 4
        
        :returns: True or False                
    '''  
    if not(A in [B,C,D]) and not(B in [C,D]) and C!=D:   
        return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)
    else:
        sys.exit('A, B, C and D must be distinct!')
# ------------------------------------------------------------------------------

# A function to construct the filename (str) of a given MAPPINGS V grid
def get_MVphotogrid_fn(Pk = 5.0, calibs = 'GCZO', kappa = np.inf, 
                        struct = 'sph', sampling = 1):
    '''
        Returns the filename of a MAPPINGS V grid, given a set of parameters.
        
        :param Pk: {float;default = 5.0} 
                    MAPPINGS model pressure. 
                    Value must match an existing reference grid file !
        :param calibs: {string; default = GCZO}
                    Input models to the MV simulations.
                    Value must match an existing reference grid file !
        :param kappa: {float; default = np.inf} 
                    The kappa value.
                    Value must match an existing reference grid file !
        :param struct: {string; default = 'sph'}
                    spherical ('sph') or plane-parallel ('pp') HII regions.
                    Value must match an existing reference grid file !
        :param sampling: {int; default = 1}
                    Use a resampled grid ?
                    
        :returns: string of the MV filename, incl. the path.
    '''
    # Only allows Pk to have a precision of 1 for now. Issue an error otherwise.
    # If people really want more precision (why?), then they should send
    # frederic.vogt@alumni.anu.edu.au an email when they read this ...
    if np.round(Pk,1) != Pk:
        sys.exit('Error: Pk cnanot have more than 1 decimals. L219')    

    if sampling >1:
        resam = '_samp_'+np.str(sampling)
    else:
        resam = ''
    if kappa == 'inf':
        kappa = np.inf
    elif kappa != np.inf:
        kappa = np.int(kappa)
    fn = os.path.join(pyqz_grid_dir,'grid_QZ_'+struct+'_'+calibs+'_Pk'+
                    np.str(np.int(10*Pk))+'_k'+np.str(kappa)+resam+'.csv')
    return fn
# ------------------------------------------------------------------------------

# A function to fetch the header information from the MAPPINGS files (csv)
def get_MVphotogrid_metadata(fn):
    '''
        Returns the first 4 lines of a MAPPINGS V file (compiled via the 
        Hawk script).
        
        :param fn: {string}
                    filename (inculding path if required) as string.
        
        :returns: dictionnary of strings
    '''
    
    # Create a dictionnary for clarity
    metadata = {}    
    
    # Open the file, and extract the info
    try:
        data_file = open(fn, 'r')
    except:
        sys.exit('Grid file not found: %s' % fn)
         
    # Store the info
    metadata['date'] = data_file.readline().split('\n')[0]
    metadata['MV_id'] = data_file.readline().split('\n')[0]
    metadata['params'] = data_file.readline().split('\n')[0]
    metadata['resampled'] = {}
    if 'samp' in fn:
        metadata['resampled'] = {}
        metadata['resampled']['info'] = data_file.readline().split('\n')[0]
        metadata['resampled']['LogQ'] = [np.float(j) for j in data_file.readline().split(']')[0].split('[')[1].split(' ') if j !='']
        metadata['resampled']['Tot[O]+12'] = [np.float(j) for j in data_file.readline().split(']')[-2].split('[')[1].split(' ') if j!='']
        
    data_content = data_file.readline().split('\n')[0].split(',')
    data_file.close()
    
    # For each column header, remove any '(x)' if they are present ...
    for (i,column_name) in enumerate(data_content):
         data_content[i] = column_name.split('(')[0]

    metadata['columns'] = data_content
    
    # All done. 
    return metadata
# ------------------------------------------------------------------------------

# A function to interpolate the MAPPINGS grid into thiner grid. 
# uses Akyma splines (1-D) in the Q/Z planes (one axis at a time).
def resample_MVphotogrid(fn, sampling = 2):
    ''' 
        A function to interpolate the MAPPINGS grid into thiner grid. 
        uses Akyma splines (1-D) in the Q/Z planes (one axis at a time).
        
        :param fn: {string}
                    filename (inculding path if required) as string.
        :param sampling: {int; default = 2}
                            the refining factor, i.e. number of interpolated
                            points between each native grid node.
    '''
    
    # 1) Get the grid metadata, for later ...
    metadata = get_MVphotogrid_metadata(fn)
    
    # 2) and now, load the raw grid
    orig_grid = np.loadtxt(fn, comments='c', delimiter=',',skiprows=4) 
    
    # 3) Extract the original grid nodes, for later
    LogQ_orig = np.sort(np.unique(orig_grid[:,metadata['columns'].index('LogQ')]))
    TotZ_orig = np.sort(np.unique(orig_grid[:,metadata['columns'].index('Tot[O]+12')]))

    # 4) Create the finer arrays.
    # Round them up to deal with "exact" node points
    TotZ_fine = np.zeros(len(TotZ_orig)+((len(TotZ_orig)-1)*(sampling-1)))
    LogQ_fine = np.zeros(len(LogQ_orig)+((len(LogQ_orig)-1)*(sampling-1)))

    for (k,Q) in enumerate(LogQ_orig[:-1]):
        LogQ_fine[sampling*k:sampling*(k+1)] = \
                    np.round(np.linspace(Q,LogQ_orig[k+1],num=sampling,
                                                    endpoint=False),4)

    for (k,Z) in enumerate(TotZ_orig[:-1]):
        TotZ_fine[sampling*k:sampling*(k+1)] = \
                    np.round(np.linspace(Z,TotZ_orig[k+1],num=sampling,
                                                    endpoint=False),4)
    # Don't forget the last element as well
    LogQ_fine[-1] = LogQ_orig[-1]
    TotZ_fine[-1] = TotZ_orig[-1]
                                                       

    # 4) Create the resampled grid, and start filling it
    fine_grid = np.zeros((len(LogQ_fine)*len(TotZ_fine),len(orig_grid[0])))
    fine_grid[:,metadata['columns'].index('LogQ')] = \
                                np.meshgrid(LogQ_fine,TotZ_fine)[0].flatten()
    fine_grid[:,metadata['columns'].index('Tot[O]+12')] = \
                                np.meshgrid(LogQ_fine,TotZ_fine)[1].flatten()

    # 5) Start the interpolation. 
    # Step 1, along the Z direction, for the original Qs
    for (k,Q) in enumerate(LogQ_orig):
        sub_grid = orig_grid[orig_grid[:,metadata['columns'].index('LogQ')]==Q,:]
        sort_ind = np.lexsort((sub_grid[:,metadata['columns'].index('Tot[O]+12')],sub_grid[:,metadata['columns'].index('LogQ')]))
        sub_grid_sorted = sub_grid[sort_ind]
        for (c,col) in enumerate([cc for cc in metadata['columns'] if
                                    not(cc in ['LogQ','Tot[O]+12'])]):        
            Akyma_func = scipy.interpolate.Akima1DInterpolator(TotZ_orig,
                                                            sub_grid_sorted[:,metadata['columns'].index(col)])
                                              
            fine_grid[fine_grid[:,metadata['columns'].index('LogQ')]==Q,metadata['columns'].index(col)] = Akyma_func(TotZ_fine)
            
    # Step 2, along the Q directions, for all the refined Zs
    for (k,Z) in enumerate(TotZ_fine):
        sub_grid = fine_grid[(fine_grid[:,metadata['columns'].index('Tot[O]+12')]==Z),:]
        sub_grid = sub_grid[[n for n in range(len(sub_grid)) if sub_grid[n,metadata['columns'].index('LogQ')] in LogQ_orig],:]
        sort_ind = np.lexsort((sub_grid[:,metadata['columns'].index('Tot[O]+12')],sub_grid[:,metadata['columns'].index('LogQ')]))
        sub_grid_sorted = sub_grid[sort_ind]
        for (c,col) in enumerate([cc for cc in metadata['columns'] if
                                    not(cc in ['LogQ','Tot[O]+12'])]):   
            Akyma_func = scipy.interpolate.Akima1DInterpolator(LogQ_orig,
                                                              sub_grid_sorted[:,metadata['columns'].index(col)])
                                                              
            fine_grid[fine_grid[:,metadata['columns'].index('Tot[O]+12')]==Z, metadata['columns'].index(col)] = Akyma_func(LogQ_fine)                                                         

    # All done ? I guess I could save it all, then ...
    file_header =  metadata['date'] +'\n'
    file_header += metadata['MV_id'] + '\n'
    file_header += metadata['params'] +'\n'
    file_header += 'Resampled with pyqz '+__version__+' (sampling='+np.str(sampling)+') '+str(dt.now())+'\n'
    file_header += 'Original LogQ: '+np.str(LogQ_orig)+'\n'
    file_header += 'Original Tot[O]+12: '+np.str(TotZ_orig)+'\n'
    for (i,col) in enumerate(metadata['columns']):
        file_header += col+'('+str(i)+'),'
        
    out_fn = fn[:-4]+'_samp_'+np.str(sampling)+'.csv'    
    np.savetxt(out_fn,fine_grid,delimiter=',',header=file_header[:-1],fmt='%0.4f', comments='')
    
    print ' '
    print '  Success: %s resampled by a factor %ix%i and saved as %s' % (fn.split('/')[-1],sampling,sampling,out_fn.split('/')[-1])
# ------------------------------------------------------------------------------    

# A function to run the Awk script provided with MAPPINGS in the loop, to
# easily generate all the grids required for the default pyqz.
def run_awk_loop(Pks = [5.0], kappas = ['inf'], structs = ['sph'], ncpu = 1,
                 awk_loc = '.'):
    '''
         Runs the awk script provided with MAPPINGS in a loop, to easily
         generate all the default MAPPINGS grid required by pyqz. 
         This function is not designed for general use - proceed with caution.
         
         :param Pks: list
                      the values of Pk to compute the grids at.
         :param kappas: list
                      the values of kappa to compute the grids at.
         :param struct: list
                      this values of struct to compute the grids at.
         :param ncpu: int
                      the number of cpus to use.
         :param awk_loc: string
                      path to the awk script provided with MAPPINGS
    '''
    
    # Run some tests first, before crashing after 27 hours ...
    try:
        if np.any(Pks <=0):
            sys.exit('Invalid Pks. Must be >0.')
    except:
        sys.exit('Invalid Pks.')
     
    for kappa in kappas:
        if not(kappa in [2,3,4,6,10,20,50,100,'inf', np.inf]):
            sys.exit('Invalid kappa value(s). Can be [2,3,4,6,10,20,50,100, "inf", np.inf].')
            
    for struct in structs:
        if not struct in ['pp','sph']:
            sys.exit('Invalid structs. Must be in ["pp","sph"].')
            
    # All good, now, let's start the show
    fn_raw = os.path.join(awk_loc,'rungrid.sh')
    fn_tmp = os.path.join(awk_loc,'tmp.sh')
    
    # Start the loop, and the grid constructions
    loop_size = len(Pks)*len(kappas)*len(structs)
    counter = 0
    for Pk in Pks:
        for kappa in kappas:
            for struct in structs:
                
                counter += 1
                
                # Open the original script
                f = open(fn_raw,'r')
                content = f.readlines()
                f.close()
                
                # Replace the specifc lines as required:
                content[57] = 'type="'+struct+'"\n'
                content[70] = 'pres="'+np.str(Pk)+'"\n'
                if not(kappa in ['inf',np.inf]):
                    content[162] = 'kappa="'+np.str(np.round(kappa,1))+'"\n'
                else:
                    content[162] = 'kappa="inf"\n'
                
                # And write this to a temporary .sh file
                f = open(fn_tmp,'w')
                f.writelines(content)
                f.close()
                
                # And make sure it can be run by a Python script ...
                os.system("chmod +x "+fn_tmp)
                
                # perfect, now, I am ready to launch this ...
                # but first, what name do I want to give to this file ?
                fn_out = get_MVphotogrid_fn(Pk=Pk,calibs='GCZO', kappa = kappa,
                                            struct=struct, sampling=1)
                fn_out = fn_out.split('grid_QZ_')[1].split('.csv')[0]
                print ' '
                print np.str(counter)+'/'+np.str(loop_size)+': Running '+fn_out
                print ' '
                sys.stdout.flush()
                
                # And launch it
                subprocess.call([fn_tmp,fn_out,np.str(np.int(ncpu))])
                
    # And don't forget to clean-up the temporary .sh file ...
    os.remove(fn_tmp)
    
    return True           
# ------------------------------------------------------------------------------