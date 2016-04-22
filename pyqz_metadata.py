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

# import the required modules
import numpy as np
import os
from matplotlib import pylab as plt

# Where are we located ?
pyqz_dir = os.path.dirname(__file__)
# Where are the reference data ?
pyqz_grid_dir = os.path.join(pyqz_dir, 'reference_data')

# Define the version of pyqz
__version__ = '0.7.3'

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

# What is the range covered by the MAPPINGS grids in the different spaces ? 
QZs_lim = {'LogQ':np.array([6.5,8.5]),
            'Tot[O]+12':np.array([8.11, 8.985]),
            'gas[O]+12':np.array([8.00, 8.875]),
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

# A list of available diagnostic grids and mixing ratios (for 3D grids)
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

