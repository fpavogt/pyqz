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
from .. import pyqz
from .. import pyqz_plots
import unittest
import sys
import os
import glob

# Where do I store figures, etc ... ?
this_dir = os.path.dirname(__file__)
arena = os.path.join(this_dir,'test_arena')

# ---------------------- Define the basic test functions -----------------------

# Test what happens when I interpolate on a grid node
# Do it only for 1 grid and 2 nodes (on the edge and in the middle)
def interpgridnodes(Pk = 5.0, kappa = np.inf, struct='pp', sampling = 1):
    
    # 1: Extract the "spectra from the grid nodes
    grid_name = '[NII]/[SII]+;[OIII]/[SII]+'
    
    grid = pyqz.get_grid(grid_name, Pk=Pk, struct=struct, kappa=kappa,
                               coeffs=pyqz.pyqzm.diagnostics[grid_name]['coeffs'], 
                               sampling=sampling)

    grid_nodes = grid[0][:,[grid[1].index('LogQ'),grid[1].index('Tot[O]+12'),
                         grid[1].index('Mix_x'),grid[1].index('Mix_y')]]
                         
    # 2: Now, interpolate and checks the output
    interp_qs = pyqz.interp_qz('LogQ', [grid_nodes[:,-2],grid_nodes[:,-1]],
                               grid_name, 
                               coeffs = pyqz.pyqzm.diagnostics[grid_name]['coeffs'],
                               Pk = Pk, kappa = kappa, struct = struct,
                               sampling = sampling) 
    interp_zs = pyqz.interp_qz('Tot[O]+12', [grid_nodes[:,-2],grid_nodes[:,-1]],
                               grid_name, 
                               coeffs = pyqz.pyqzm.diagnostics[grid_name]['coeffs'],
                               Pk = Pk, kappa = kappa, struct = struct,
                               sampling = sampling) 

    return (np.all(np.round(interp_qs,2) == grid_nodes[:,0]) and
       np.all(np.round(interp_zs,3) == grid_nodes[:,1]))
       
# Test whether points outside the grids are given as nan
def interpoffgrid(Pk = 5.0, kappa = np.inf, struct = 'pp', sampling=1):

    grid_name = pyqz.pyqzm.diagnostics.keys()[0]
    
    interp_q = pyqz.interp_qz('LogQ', [np.array(-1000),np.array(-1000)],
                               grid_name, 
                               coeffs = pyqz.pyqzm.diagnostics[grid_name]['coeffs'],
                               Pk = Pk, kappa = kappa, struct = struct,
                               sampling = sampling)           
        
    return np.isnan(interp_q)

# When I interpolate a MAPPINGS V simulations done at different Qs, do I get
# the Qs values out correctly ?    
def interp_midMVq(Pk = 5.0, kappa = np.inf, struct = 'pp', sampling = 1,
                  KDE_qz_sampling = 101j, do_single_spec = True, KDE_pickle_loc = None,
                  error = 0.05,nproc=1, diags = ['[NII]/[SII]+;[OIII]/[SII]+']):
        
        # 1: get the intermediate data points from the MV shell script
        # Zs are identical, but at least we have changed the Qs
        
        # Where are we ?
        #this_dir = os.path.dirname(__file__)
        #fn = os.path.join(this_dir,'grid_QZ_midQs_pp_GCZO_Pk50_kinf.csv')
        fn = os.path.join(arena,'grid_QZ_midQs_pp_GCZO_Pk50_kinf.csv')
        
        metadata = pyqz.pyqzt.get_MVphotogrid_metadata(fn)
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
        nspec = 16 #16 on grid,14 = grid edge
        
        if do_single_spec:
            specs = np.reshape(all_fluxes[nspec,:],(1,12))
        else:
            specs = np.ones((24,12))*np.reshape(all_fluxes[nspec,:],(1,12))
            
        results = pyqz.get_global_qz(specs, 
                                        ['Hb','stdHb','[OIII]','std[OIII]',
                                        '[OII]+','std[OII]+','[NII]','std[NII]',
                                        '[SII]+','std[SII]+','Ha','stdHa'],
                                        diags,
                                        Pk = 5.0, kappa = np.inf, struct = 'pp',
                                        KDE_pickle_loc = KDE_pickle_loc,
                                        sampling=sampling,
                                        KDE_qz_sampling = KDE_qz_sampling,
                                        KDE_method='multiv',srs=400,
                                        nproc=nproc,
                                        verbose = False)                                                      
        
        return (np.abs(results[0][0,results[1].index('<LogQ>')]/data[nspec,metadata['columns'].index('LogQ')]-1.)<0.01 
                and                                                                             
                np.abs(results[0][0,results[1].index('<LogQ{KDE}>')]/data[nspec,metadata['columns'].index('LogQ')]-1.)<0.01
                and 
                np.abs(results[0][0,results[1].index('<Tot[O]+12>')]/data[nspec,metadata['columns'].index('Tot[O]+12')]-1.)<0.01 
                and
                np.abs(results[0][0,results[1].index('<Tot[O]+12{KDE}>')]/data[nspec,metadata['columns'].index('Tot[O]+12')]-1.)<0.01
               )                                                                                                                                                                                                                            

# Run get_global_qz with a data points outside all diagnostics
def get_bad_global_qz(Pk = 5.0, kappa = np.inf, struct = 'pp', sampling = 1,
                      nproc=1, error = 0.05, KDE_pickle_loc = None):
                      
        # 1: get the intermediate data points from the MV shell script
        # Zs are identical, but at least we have chnaged the Qs
        
        # Where are we ?
        #this_dir = os.path.dirname(__file__)
        #fn = os.path.join(this_dir,'grid_QZ_midQs_pp_GCZO_Pk50_kinf.csv')
        fn = os.path.join(arena,'grid_QZ_midQs_pp_GCZO_Pk50_kinf.csv')
        
        metadata = pyqz.pyqzt.get_MVphotogrid_metadata(fn)
        data = np.loadtxt(fn, comments='c', delimiter=',',skiprows=4)
        
        # Build 'Pseudo' line fluxes
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
            all_fluxes[i,2] = Oiii[i]*1e3 # Make this an obviously bad line.
            all_fluxes[i,3] = all_fluxes[i,2] * error
            all_fluxes[i,4] = Oiip[i] / 1e3 # Idem
            all_fluxes[i,5] = all_fluxes[i,4] * error
            all_fluxes[i,6] = Nii[i]
            all_fluxes[i,7] = all_fluxes[i,6] * error
            all_fluxes[i,8] = Siip[i]
            all_fluxes[i,9] = all_fluxes[i,8] * error  
            all_fluxes[i,10] = Ha[i]
            all_fluxes[i,11] = all_fluxes[i,10] * error 
            
        # Launch the interpolation
        nspec = 16 #16 on grid,14 = grid edge
        results = pyqz.get_global_qz(np.reshape(all_fluxes[nspec,:],(1,12)), 
                                        ['Hb','stdHb','[OIII]','std[OIII]',
                                        '[OII]+','std[OII]+','[NII]','std[NII]',
                                        '[SII]+','std[SII]+','Ha','stdHa'],
                                        ['[NII]/[SII]+;[OIII]/Hb',
                                         '[NII]/[OII]+;[OIII]/[SII]+'],
                                        Pk = 5.0, kappa = np.inf, struct = 'pp',
                                        KDE_pickle_loc = KDE_pickle_loc,
                                        sampling=1,
                                        KDE_qz_sampling = 201j,
                                        KDE_method='multiv',srs=400,
                                        nproc=nproc,
                                        verbose = False)                                                      
        
        return np.isnan(results[0][0,results[1].index('<LogQ>')])                                                                                                                                                                                                                                               
        
# --------------------- Invoke the basic test unit tools -----------------------
class Testpyqz(unittest.TestCase):
  
  def test01_interpgridnodes(self):
      
      out = interpgridnodes()      
        
      self.assertTrue(out)
  
  def test02_interpoffgrid(self):
      
      out = interpoffgrid()      
        
      self.assertTrue(out)
  
  def test03_interp_midMVq(self):
      
      out = interp_midMVq(nproc = 1, KDE_pickle_loc = arena)
      
      self.assertTrue(out)
  
  def test04_get_bad_global_qz(self):
     
      out = get_bad_global_qz()
          
      self.assertTrue(out)      
              
  def test05_multiprocessing(self):
      
      out = interp_midMVq(nproc=2, KDE_pickle_loc = arena)
          
      self.assertTrue(out) 
         
  def test06_speed_benchmark(self):
      main_diags = ['[NII]/[SII]+;[OIII]/[SII]+',
                    '[NII]/[SII]+;[OIII]/Hb',
                    '[NII]/[SII]+;[OIII]/[OII]+',
                    '[NII]/[OII]+;[OIII]/[OII]+',
                    '[NII]/[OII]+;[OIII]/[SII]+',
                    '[NII]/[OII]+;[SII]+/Ha'
                    ]
                    
      out = interp_midMVq(nproc=1, KDE_pickle_loc = None, diags = [main_diags[0]],
                          sampling=1,KDE_qz_sampling=101j, do_single_spec = True)              
      out = interp_midMVq(nproc=1, KDE_pickle_loc = None, diags = main_diags,
                          sampling=1,KDE_qz_sampling=101j, do_single_spec = True)
      out = interp_midMVq(nproc=1, KDE_pickle_loc = None, diags = main_diags,
                          sampling=2,KDE_qz_sampling=201j, do_single_spec = True)    
      out = interp_midMVq(nproc=1, KDE_pickle_loc = arena, diags = main_diags,
                          sampling=2,KDE_qz_sampling=201j, do_single_spec = True)                                     
      
      out = interp_midMVq(nproc=8, KDE_pickle_loc = None, diags = main_diags,
                          sampling=2,KDE_qz_sampling=201j, do_single_spec = False) 
        
            
      self.assertTrue(out)  
  
  def test07_plots(self):
      
      fns = glob.glob(os.path.join(arena,'*.pkl'))
      
      try:
          for fn in fns:
          
              pyqz_plots.plot_global_qz(fn, do_all_diags = True, show_plots = False, 
                                        save_loc = arena)
      
          out = True
      except:
          out = False
            
      self.assertTrue(out)
  
  
def run_all_tests(test_suite=Testpyqz, cleanup=True):
   # launch the testing
   #print ' '
   #print ' Starting pyqz tests:'
   #print' '
   #sys.stdout.flush()

   suite = unittest.TestLoader().loadTestsFromTestCase(Testpyqz)
   unittest.TextTestRunner(verbosity=2).run(suite)
   
   
   # Cleanup all the mess I just created.
   if cleanup:
      fns = os.listdir(arena)
      for fn in fns:
         if fn[-4:] in ['.pdf','.pkl']:
            os.remove(os.path.join(arena,fn))
  

if __name__ == '__main__':
   # launch the testing
   #print ' '
   #print ' Starting pyqz tests:'
   print ' That was never meant to be printed !!!'
   #sys.stdout.flush()

   #suite = unittest.TestLoader().loadTestsFromTestCase(Testpyqz)
   #unittest.TextTestRunner(verbosity=2).run(suite)
