# -*- coding: utf-8 -*-
#
# This small program is designed to run some self-consistency tests on the pyqz 
# module.
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

# import the required packages
import numpy as np
import pyqz
import unittest
import sys
import os

# ---------------------- Define the basic test functions -----------------------

# Test what happens when I interpolate on a grid node
# Do it only for 1 grid and 2 nodes (on the edge and in the middle)
def interpgridnodes(Pk = 5.0, kappa = np.inf, struct='pp', sampling = 1):
    
    # 1: Extract the "spectra from the grid nodes
    grid_name = '[NII]/[SII]+;[OIII]/[SII]+'
    
    grid = pyqz.get_grid(grid_name, Pk=Pk, struct=struct, kappa=kappa,
                               coeffs=pyqz.diagnostics[grid_name]['coeffs'], 
                               sampling=sampling)

    grid_nodes = grid[0][:,[grid[1].index('LogQ'),grid[1].index('Tot[O]+12'),
                         grid[1].index('Mix_x'),grid[1].index('Mix_y')]]
                         
    # 2: Now, interpolate and checks the output
    interp_qs = pyqz.interp_qz('LogQ', [grid_nodes[:,-2],grid_nodes[:,-1]],
                               grid_name, 
                               coeffs = pyqz.diagnostics[grid_name]['coeffs'],
                               Pk = Pk, kappa = kappa, struct = struct,
                               sampling = sampling, show_plot = False,
                               save_plot = False) 
    interp_zs = pyqz.interp_qz('Tot[O]+12', [grid_nodes[:,-2],grid_nodes[:,-1]],
                               grid_name, 
                               coeffs = pyqz.diagnostics[grid_name]['coeffs'],
                               Pk = Pk, kappa = kappa, struct = struct,
                               sampling = sampling, show_plot = False,
                               save_plot = False) 

    return (np.all(np.round(interp_qs,2) == grid_nodes[:,0]) and
       np.all(np.round(interp_zs,3) == grid_nodes[:,1]))
       
# Test whether points outside the grids are given as nan
def interpoffgrid(Pk = 5.0, kappa = np.inf, struct = 'pp', sampling=1):

    grid_name = pyqz.diagnostics.keys()[0]
    
    interp_q = pyqz.interp_qz('LogQ', [np.array(-1000),np.array(-1000)],
                               grid_name, 
                               coeffs = pyqz.diagnostics[grid_name]['coeffs'],
                               Pk = Pk, kappa = kappa, struct = struct,
                               sampling = sampling, show_plot = False,
                               save_plot = False)           
        
        
    return np.isnan(interp_q)

# When I interpolate a MAPPINGS V simulations done at different Qs, do I get
# the Qs values out correctly ?    
def interp_midMVq(Pk = 5.0, kappa = np.inf, struct = 'pp', sampling = 2,
                  error = 0.1):
        
        # 1: get the intermediate data points from the MV shell script
        # Zs are identical, but at least we have chnaged the Qs
        
        # Where are we ?
        this_dir = os.path.dirname(__file__)
        
        fn = os.path.join(this_dir,'grid_QZ_midQs_pp_GCZO_Pk50_kinf.csv')
        metadata = pyqz.get_MVphotogrid_metadata(fn)
        data = np.loadtxt(fn, comments='c', delimiter=',',skiprows=4)
        
        # Build 'Pseudo' line fluxes - and then assume a 10% error
        Hb = np.ones_like(data[:,0])
        Oiii = 10**data[:, metadata['columns'].index('[OIII]/Hb')]
        Oiip = 10**data[:, metadata['columns'].index('[OII]+/Hb')]
        Nii = Oiip * 10**data[:, metadata['columns'].index('[NII]/[OII]+')]
        Siip = 1./10**data[:, metadata['columns'].index('[NII]/[SII]+')] * Nii
        Ha =  1./10**data[:, metadata['columns'].index('[NII]/Ha')] * Nii              
        
        all_fluxes = np.zeros((len(Hb),12))
        for i in range(len(Hb)):
            all_fluxes[i,0] = 1.0
            all_fluxes[i,1] = all_fluxes[i,0] * error
            all_fluxes[i,2] = Oiii[i]
            all_fluxes[i,3] = all_fluxes[i,2] * error
            all_fluxes[i,4] = Oiip[i]
            all_fluxes[i,5] = all_fluxes[i,4] * error
            all_fluxes[i,6] = Nii[i]
            all_fluxes[i,7] = all_fluxes[i,6] * error
            all_fluxes[i,8] = Siip[i]
            all_fluxes[i,9] = all_fluxes[i,8] * error  
            all_fluxes[i,10] = Ha[i]
            all_fluxes[i,11] = all_fluxes[i,10] * error 
            
        # Launch the interpolation
        nspec = 14 #16
        results = pyqz.get_global_qz(np.reshape(all_fluxes[nspec,:],(1,12)), 
                                        ['Hb','stdHb','[OIII]','std[OIII]',
                                        '[OII]+','std[OII]+','[NII]','std[NII]',
                                        '[SII]+','std[SII]+','Ha','stdHa'],
                                        ['[NII]/[SII]+;[OIII]/[SII]+',
                                        '[NII]/[SII]+;[OIII]/Hb',
                                        '[NII]/[SII]+;[OIII]/[OII]+',
                                        '[NII]/[OII]+;[OIII]/[OII]+',
                                        '[NII]/[OII]+;[OIII]/[SII]+',
                                        '[NII]/[OII]+;[OIII]/Hb'],
                                        Pk = 5.0, kappa = np.inf, struct = 'pp',
                                        show_plot = True,save_plot=True,
                                        plot_loc = './plots',
                                        sampling=1,
                                        KDE_qz_sampling = 101j,
                                        KDE_method='gauss',srs=400,
                                        verbose = False)                                                      
        print ' '
        print results[0][0,results[1].index('<LogQ>')]
        print results[0][0,results[1].index('<LogQ{KDE}>')]
        print data[nspec,metadata['columns'].index('LogQ')]
        print ' '
        print results[0][0,results[1].index('<Tot[O]+12>')]
        print results[0][0,results[1].index('<Tot[O]+12{KDE}>')]
        print data[nspec,metadata['columns'].index('Tot[O]+12')]
        
        return (np.abs(results[0][0,results[1].index('<LogQ>')]/data[nspec,metadata['columns'].index('LogQ')]-1.)<0.01 
                and                                                                             
                np.abs(results[0][0,results[1].index('<LogQ{KDE}>')]/data[nspec,metadata['columns'].index('LogQ')]-1.)<0.01
                and 
                np.abs(results[0][0,results[1].index('<Tot[O]+12>')]/data[nspec,metadata['columns'].index('Tot[O]+12')]-1.)<0.01 
                and
                np.abs(results[0][0,results[1].index('<Tot[O]+12{KDE}>')]/data[nspec,metadata['columns'].index('Tot[O]+12')]-1.)<0.01
               )                                                                                                                                                                                                                            
        
        
# --------------------- Invoke the basic test unit tools -----------------------
class Testpyqz(unittest.TestCase):
  '''
  def test_upper(self):
      self.assertEqual('foo'.upper(), 'FOO')

  def test_isupper(self):
      self.assertTrue('FOO'.isupper())
      self.assertFalse('Foo'.isupper())

  def test_split(self):
      s = 'hello world'
      self.assertEqual(s.split(), ['hello', 'world'])
      # check that s.split fails when the separator is not a string
      with self.assertRaises(TypeError):
          s.split(2)
  '''
  
  def test_interpgridnodes(self):
      
      out = interpgridnodes()      
        
      self.assertTrue(out)
  
  def test_interpoffgrid(self):
      
      out = interpoffgrid()      
        
      self.assertTrue(out)
  
  def test_interp_midMVq(self):
      
      out = interp_midMVq()
      
      self.assertTrue(out)

# launch the testing
print ' '
print ' Starting pyqz tests:'
print' '
sys.stdout.flush()

suite = unittest.TestLoader().loadTestsFromTestCase(Testpyqz)
unittest.TextTestRunner(verbosity=2).run(suite)
