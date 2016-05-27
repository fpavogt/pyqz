.. pyqz documentation master file, created by
   sphinx-quickstart on Fri Nov 28 18:17:46 2014.

pyqz |release|
==================

The pyqz Python module computes the values of log(Q) [the ionization parameter] and 12+log(O/H) [the oxygen abundance, either total or in the gas phase] for a given set of strong emission lines fluxes from HII regions. 

The log(Q) and 12+log(O/H) values are interpolated from a finite set of **diagnostic line ratio grids** computed with the MAPPINGS code. The grids used by pyqz are chosen to be **flat, without wraps, to decouple the influence of log(Q) and 12+log(O/H)** on the emission line ratios.

pyqz 0.4 was the first publicly released version of the code, which is described in detail in

  Dopita et al., *New Strong Line Abundance Diagnostics for HII Regions: Effects of Updated Atomic Data and kappa-Distributed Electron Energies*, ApJS, 208, 10 (2013). `ADS entry <http://adsabs.harvard.edu/abs/2013ApJS..208...10D>`_ 

pyqz has since been subject to a major overhaul to track the latest heroic developments in 
the MAPPINGS code, support the propagation of observational errors on the emission lines 
fluxes and the construction of the probability density function associated with the estimates 
of log(Q) and 12+log(O/H), auto-detect wraps in the diagnostics grids, and more.

.. note::

   You can also track the latest changes in the code on the dedicated Github repository: 
   https://github.com/fpavogt/pyqz

   See also the :ref:`changelog`.


.. warning:: **Garbage in, garbage out !** 

    Simplicity and ease of use are important design drivers for pyqz. These (hopefully!) imply an easy-to-run code, where all parameters are accessible but most come with a default setting. 

    But be aware that **the default settings will not be optimum for everyone**. Before you use pyqz, read carefully the :ref:`understandingpyqz` section, or chances are, you will be deriving wrong estimates of 12+log(O/H) and log(Q) ! No one wants this to happen.

Contents
---------
.. toctree::
   :maxdepth: 1
   
   Home <self>
   installation
   pyqz_demo_basic/pyqz_demo_basic
   pyqz_demo_advanced/pyqz_demo_advanced
   understanding
   pyqz_demo_param/pyqz_demo_param
   modules/modules
   faq
   changelog
   acknowledge

----

Copyright notice:
 
This file is part of the pyqz Python module.
The pyqz Python module is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
 
The pyqz Python module is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with the pyqz Python module.  If not, see http://www.gnu.org/licenses/ .
 

