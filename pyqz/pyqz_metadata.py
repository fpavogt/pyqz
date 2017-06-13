# -*- coding: utf-8 -*-
'''
 pyqz: a Python module to derive the ionization parameter and oxygen abundance of HII 
 regions, by comparing their strong emission line ratios with MAPPINGS photoionization 
 models.\n
 Copyright (C) 2014-2017,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains some generic metadata used throughout the pyqz module, including the 
 version number, etc ...

 Updated May 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

# import the required modules
import numpy as np
import os
from matplotlib import pylab as plt

# Where are we located ?
pyqz_dir = os.path.dirname(__file__)

# Where are the reference data ?
pyqz_grid_dir = os.path.join(pyqz_dir,'ref_data')

# Define the version of pyqz
__version__ = '0.8.4'

# For test purposes
fullgrid_x, fullgrid_y = np.mgrid[-3:3:0.01,-3:3:0.01]       

# Default colormap
pyqz_cmap_0 = 'Paired'

# Define a custom colorbar for the PDF plot - just because I can.
cbdict = {
'red'  :  ((0.0, 1, 1),(1.00, 0.3, 0.3)),
'green':  ((0.0, 1, 1),(1.00, 0.3, 0.3)),
'blue' :  ((0.0, 1, 1),(1.00, 0.3, 0.3))
}
pyqz_cmap_1 = plt.matplotlib.colors.LinearSegmentedColormap('light_gray', 
                                                            cbdict, 1024)
# and define a color for nan's and other bad points
pyqz_cmap_1.set_bad(color=(1,1,1), alpha=1) 

# Choose which MAPPINGS version is to be used (use only if you know what you are doing!)
M_version = 'MV' # 'MV' or 'MIV'.

# What is the range covered by the MAPPINGS grids in the different spaces ? 
QZs_lim = {'MV':{'LogQ':np.array([6.5,8.5]),
                 'Tot[O]+12':np.array([8.11, 8.985]),
                 'gas[O]+12':np.array([8.00, 8.875]),
                },
           'MIV':{'LogQ':np.array([6.5,8.5]),
                  'Tot[O]+12':np.array([7.39, 9.39]),
                  'gas[O]+12':np.array([7.39, 9.39]),
                 },       
          }
          
# Level of the contours to derive the best qz value from the PDF 
# (normalized to the peak)          
PDF_cont_level = 0.61  

# Reconstructing the Full PDF over the entier QZ plane is very VERY slow. 
# In practice, we already know that outside the QZ estimates (localized),the 
# PDF will be very close to 0 - not worth spending time calculating those 
# estimates then ! 
# But if you really want the full PDF, set the following to True. 
# WARNING: these will result in an exectution time 25 times slower ! 
do_full_KDE_reconstruction = False        

# A list of available diagnostic grids and mixing ratios
diagnostics = {'[NII]/[SII]+;[OIII]/[SII]+':{'coeffs':[[1,0],[0,1]]},
               # ---
               '[NII]/[SII]+;[OIII]/Hb':{'coeffs':[[1,0],[0,1]]},
               # ---
               '[NII]/[SII]+;[OIII]/[OII]+':{'coeffs':[[1,0],[0,1]]},
               # ---
               '[NII]/[OII]+;[OIII]/[OII]+':{'coeffs':[[1,0],[0,1]]},
               # ---
               '[NII]/[OII]+;[OIII]/[SII]+':{'coeffs':[[1,0],[0,1]]},
               # ---
               #'[NII]/[OII]+;[OIII]/Hb':{'coeffs':[[1,0],[0,1]]},
               # ---
               '[NII]/[OII]+;[SII]+/Ha':{'coeffs':[[1,0],[0,1]]},
               # ---
               '[OIII]4363/[OIII];[OIII]/[SII]+':{'coeffs':[[1,0],[0,1]]},
               # ---
               '[OIII]4363/[OIII];[OIII]/[OII]+':{'coeffs':[[1,0],[0,1]]},
               # ---
               '[OIII]4363/[OIII];[SII]+/Ha':{'coeffs':[[1,0],[0,1]]},
               # ---
               #'[NII]/[SII]+;[SII]+/Ha':{'coeffs':[[1,0],[0,1]]},
               # ---
               #'[NII]/Ha;[NII]/[SII]+':{'coeffs':[[1,0],[0,1]]},
               # ---
               ### And now some 3-D line ratios diagnostics
               # From Dopita (2015) Hi-z
               '[NII]/[SII]+;[NII]/Ha;[OIII]/Hb':{'coeffs':[[1.0,0.264,0.0],
                                                            [0.242,-0.910,0.342]
                                                           ]},
               }

