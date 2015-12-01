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

# import the required modules
import numpy as np

# Define the version of pyqz
__version__ = '0.7.1'

# For test purposes
fullgrid_x, fullgrid_y = np.mgrid[-3:3:0.01,-3:3:0.01]
testF = np.array([[2.88,0.11,1.0,0.05,2.015,0.045,0.036,0.002,
                    2.72,0.12,0.497,0.017,0.500,0.016]])
testF_names = ['[OII]+','std[OII]+','Hb','stdHb','[OIII]','std[OIII]',
               '[OI]','std[OI]','Ha', 'stdHa','[NII]','std[NII]', 
               '[SII]+','std[SII]+']
test_diag = ['[NII]/[SII]+;[OIII]/[SII]+', 
             '[NII]/[SII]+;[OIII]/Hb',
             '[NII]/[OII]+;[OIII]/[OII]+',
             '[NII]/[SII]+;[NII]/Ha;[OIII]/Hb',
            ]          

# Default colormap
pyqz_cmap_0 = 'Paired'

# What is the range covered by the MAPPINGS grids in the different spaces ? 
QZs_lim = {'LogQ':np.array([6.5,8.5]),
            'Tot[O]+12':np.array([7.459, 9.237]),
            'gas[O]+12':np.array([7.349, 9.127]),
          }
          
# Level of the contours to derive the best qz value from the PDF 
# (normalized to the peak)          
PDF_cont_level = 0.61          


# A list of available diagnostic grids, mixing ratios (for 3D grids) and 
# limiting ranges (the latter carry no scientific meaning, they only help to 
# create decent plots)
diagnostics = {'[NII]/[SII]+;[OIII]/[SII]+':
                    {'coeffs':[[1,0],[0,1]],
                    'xliml':-0.8, 'xlimr':1.0, 'ylimb': -4.0, 'ylimt':3.0},
               '[NII]/[SII]+;[OIII]/Hb':
                    {'coeffs':[[1,0],[0,1]],
                    'xliml':-0.8, 'xlimr':1.0, 'ylimb': -4.0, 'ylimt':2.0},
               '[NII]/[SII]+;[OIII]/[OII]+':
                    {'coeffs':[[1,0],[0,1]],
                    'xliml':-0.8, 'xlimr':1.0, 'ylimb': -2.5, 'ylimt':1.5},
               '[NII]/[OII]+;[OIII]/[OII]+':
                    {'coeffs':[[1,0],[0,1]],
                    'xliml':-2.0, 'xlimr':1.5, 'ylimb': -2.5, 'ylimt':1.0},
               '[NII]/[OII]+;[OIII]/[SII]+':
                    {'coeffs':[[1,0],[0,1]],
                    'xliml':-2.0, 'xlimr':1.5, 'ylimb': -4.0, 'ylimt':3.0},
               '[NII]/[OII]+;[OIII]/Hb':
                    {'coeffs':[[1,0],[0,1]],
                    'xliml':-2.0, 'xlimr':1.5, 'ylimb': -4.0, 'ylimt':2.0},
               '[NII]/Ha;[OIII]/Hb':
                    {'coeffs':[[1,0],[0,1]],
                    'xliml':-3.0, 'xlimr':0.0, 'ylimb': -4.0, 'ylimt':2.0},
               '[NII]/Ha;[OIII]/[OII]+':
                    {'coeffs':[[1,0],[1,0]],
                    'xliml':-3.0, 'xlimr':0.0, 'ylimb': -2.5, 'ylimt':1.5},
                ### And now some 3-D line ratios diagnostics
                '[NII]/[SII]+;[NII]/Ha;[OIII]/Hb': # From Dopita (2015) Hi-z
                    {'coeffs':[[1.0,0.264,0.0],[0.242,-0.910,0.342]],
                    'xliml':-1.5, 'xlimr':0.7, 'ylimb': -0.7, 'ylimt':3.0},
               }

