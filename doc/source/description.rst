.. _runningpyqz:

Running pyqz
================
From a set of emission line fluxes and errors, ``pyqz`` measures the associated value of the oxygen abundance 12+log(O/H) (a.k.a. "z") and ionization parameters log(Q) (a.k.a. "q"), given a set of MAPPINGS simulations of HII regions. The code uses **flat(-ish)** emission line diagnostic grids to disentangle and interpolate the values of log(Q) and 12+log(O/H).  

This page describes the basic use of the module from a **practical** point of view. It is assumed that the reader is familiar with the `Dopita et al. (2013) <http://adsabs.harvard.edu/abs/2013ApJS..208...10D>`_ article (and references therein) that describes **the theory** behind ``pyqz``. See also Appendix B in `Vogt et al. (2015) <http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1504.03337>`_.


Before we begin: a note on the pyqz syntax
------------------------------------------

The ``pyqz`` module is intimately linked to the MAPPINGS code. While both are stand alone and distinct programs, ``pyqz`` was designed to employ the same notation conventions than that of the MAPPINGS code for clarity, both from a user and programming perspective. 

These conventions, designed to maximise clarity while minimizing the overall character counts, are as follows:

* the ionization parameter is ``logQ``
* the oxygen abundance is ``Tot[O]+12`` (total) or ``gas[O]+12`` (for the gas-phase)
* the Balmer lines from Hydrogen are ``Ha``, ``Hb``, etc ...
* forbidden lines are marked as ``[OIII]``, ``[NII]``, ``[SII]``, ``[OI]``, etc ...
* for the usual strong line doublets, when the doublet line fluxes are considered together (i.e. [OIII]5007 + 4959), a ``+`` is appended to the said emission line, e.g. ``[OIII]+``. By convention, the single line is always the strongest within the doublet. In short, ``[OIII]`` corresponds to [OIII]5007, ``[OIII]+`` corresponds to [OIII]5007+4959, ``[NII]`` corresponds to [NII]6584, ``[NII]+`` corresponds to [NII]6584+6548, etc ...

This syntax must be followed carefully when using ``pyqz``, or errors will arise.
 
The structure of ``pyqz``
-------------------------

The pyqz module is composed of a core function: ``pyqz.get_qz``. 
This function is responsible for interpolating the MAPPINGS IV grid of simulations of HII regions (using ``scipy.interpolate.griddata``) and returns the corresponding value of z or q for a given pair of line ratios. This function is basic, in that it does not propagate errors on its own. You feed it a pair of line ratio, it returns q or z, and that's it. The spirit of the ``pyqz.get_qz`` function really is to be embedded inside a larger Python script to take care of these extra-tasks.

The function ``pyqz.get_qz_ff`` is a wrapper around ``pyqz.get_qz``. It is designed as a top interaction layer for the pyqz module, and can propagate errors or upper-limits on the line flux measurements. It was designed to be self-contained and therefore requires little to no Python coding at all to be run (you still need a Python shell to run it, though!). You feed it an input file (with a specific format) containing all the line fluxes and associated errors, and ``pyqz.get_qz_ff`` returns a text file with all the q and z estimates and the associated errors.

The input text file
"""""""""""""""""""

The input text file that is being fed to ``pyqz.get_qz_ff`` serves two purposes: a) to give it your line fluxes and errors, and b) to tell pyqz which line ratio diagnostics you want to use. The structure of the input file must be as follows:
::
  Name OII dOII Hb dHb OIII dOIII OI dOI Ha dHa NII dNII SII dSII NII/SII;OIII/SII NII/SII;OIII/Hb NII/SII;OIII/OII NII/OII;OIII/SII NII/OII;OIII/Hb NII/OII;OIII/OII NII/Ha;OIII/Hb NII/Ha;OIII/OII
  628-69-208 2.88 0.11 1.0 0.05 2.015 0.045 0.036 0.002 2.72 0.12 0.497 0.017 0.500 0.016
  925-12-066 2.54 0.12 1.0 0.05 3.807 0.114 0.034 0.002 2.84 0.14 0.216 0.009 0.339 0.013
  # This is a commented line

The first line contains the description of the content of each column. It also contains the list of diagnostic diagrams that you want to use. Currently, 8 line ratio diagrams are supported by pqz. These have been selected as they allow to clearly disentangle the values of log(q) and 12+og(O/H) (i.e. the MAPPINGS IV grid does not wrap over itself).
  
1. 'NII/SII;OIII/SII'
2. 'NII/SII;OIII/Hb'
3. 'NII/SII;OIII/OII'
4. 'NII/OII;OIII/SII'
5. 'NII/OII;OIII/Hb'
6. 'NII/OII;OIII/OII'
7. 'NII/Ha;OIII/Hb'
8. 'NII/Ha;OIII/OII'

You do not have to use all of them. For example, if you wanted to use only the diagnostics that do not rely on [OII], the input file would be:
::
  Name OII dOII Hb dHb OIII dOIII OI dOI Ha dHa NII dNII SII dSII NII/SII;OIII/SII NII/SII;OIII/Hb NII/Ha;OIII/Hb
  628-69-208 2.88 0.11 1.0 0.05 2.015 0.045 0.036 0.002 2.72 0.12 0.497 0.017 0.500 0.016
  925-12-066 2.54 0.12 1.0 0.05 3.807 0.114 0.034 0.002 2.84 0.14 0.216 0.009 0.339 0.013

The # symbol is used to write comments. Any line starting with # will be ignored by pyqz. The 'Name' column is not strictly speaking required, but if it exists, figures will be saved with that string inside their filename (avoid spaces at all costs!). Hence, if you are processing multiple spectra at once, providing a 'Name' for each of them will ensure that all the Figures will be saved separately. Missing values should be marked (by default) as '$$$', or alternatively, the 'no-value-string' can be passed as an argument to ``pyqz.get_qz_ff``. Each item must be separated by a single space, or alternatively by a single 'tab', so that **object names should not contain any space !**

Each (non-commented) line in the file corresponds to one set of line fluxes (i.e. 1 spectrum, 1 spaxel, etc …). The file can contain as many lines as you want, and pyqz will process them one at a time (multi-processing is not yet supported). The errors provided are assumed to be the 1-sigma standard deviation from the mean, where the probability density distribution is assumed to be gaussian. All lines errors are also assumed to be uncorrelated. 

.. NOTE:: If your errors are not gaussian, and/or the errors from different lines are correlated, it is in principle possible to update pyqz accordingly. Experienced Python users can look at and update the ``pyqz.get_qz_ff`` function as needed, and/or contact us to discuss about how to implement the update.

In pyqz v0.6.1, support was added for upper limits on the line flux measurements. To mark a given flux measurement as an upper-limit, simply set its error to -1.0 in the input file. In that case, the random sample of additional line fluxes will be derived from a uniform probability distribution between 0 and the flux value provided. 


Launching the routine
"""""""""""""""""""""

Open a IPython shell, cd to the location of your input file, and type:
:: 
    import pyqz
    pyqz.get_qz_ff(kappa,fn) 


pyqz will then process the data provided in the file 'fn', assuming a value of 'kappa'. You can also specify a few keywords to show and/or save some/all/none of the diagrams, to select the output format, etc… The detailed list is listed in :ref:`ref-funcs`.

Some keywords or of particular interest:

* **srs**: the 'size of the random sample' of line fluxes generated by pyqz. This is used to propagate the probability density function associated with each line flux measurements. In other words, this is the number of discrete estimates of the probability density function (in the z-q plane) associated with one diagnostic grid. Hence, the joint probability function density function (combining N diagnostic grids) will be reconstructed via a Kernel Density Estimation routine from N x srs points. srs=400 is the default value, and we suggest a srs=1000 for errors at the 10%-15% level. Larger errors will result in larger probability density peaks, and will require more srs points to be properly discretized - at the cost of additional computation time. 

* **KDE_method**: the Kernel Density Estimation routine used to reconstruct the joint probability density function the z-q plane. Either ``gaussian_kde`` from the ``scipy.stats`` module, or ``KDEMultivariate`` from the ``statsmodels`` package. The former choice is 10-100x faster, but can result is less accurate results if their is disagreement between different diagnostics. The underlying reason is that the kernel bandwidth cannot be explicitly defined individually for the q and z, so that the function tends to over-smooth the distribution. ``KDEMultivariate`` should be preferred as the bandwidth of the kernel is set individually for both the q and z direction using Scott's rule, scaled by the standard deviation of the distribution along the q or z direction.


The output file
"""""""""""""""

By default, the output file is named following the input file, plus ``_out_kxx.txt``, where xx is replaced by the chosen kappa value. It contains all of the information of the input file, plus the q and z estimates derived by pyqz. For example, the output file corresponding to the first example above looks like that:
::
    Name OII dOII Hb dHb OIII dOIII OI dOI Ha dHa NII dNII SII dSII NII/SII;OIII/SII[log_q] NII/SII;OIII/SII[log_z] NII/SII;OIII/Hb[log_q] NII/SII;OIII/Hb[log_z] <q> std[q] <z> std[z] <q_rs> err[q_rs] <z_rs> err[z_rs] flag rs_offgrid 
    628-69-208 2.88 0.11 1.0 0.05 2.015 0.045 0.036 0.002 2.72 0.12 0.497 0.017 0.5 0.016 7.69436 8.59547 8.16709 8.58808 7.61864 0.06147 8.82755 0.01528 7.55806 0.04867 8.84216 0.02063 0 0.0 
    925-12-066 2.54 0.12 1.0 0.05 3.807 0.114 0.034 0.002 2.84 0.14 0.216 0.009 0.339 0.013 7.69436 8.59547 8.16709 8.58808 7.93073 0.23637 8.59177 0.00369 7.69022 0.10611 8.59488 0.0297 2 11.5 

The new columns include:
* for each diagnostic, the estimate of q and z associated with the set of line fluxes (labelled as ``ratio1;ratio2[log_q/z]``

* the mean and standard deviation of these estimates, labelled ``<q>``, ``std[q]``, ``<z>`` and ``std[z]``. Keep in mind that these do not take any error into account. 

* the best estimates of q and z based on the joint probability density function, and associated uncertainty, labelled ``<q_rs>``, ``err[q_rs]``, ``<z_rs>`` and ``err[z_rs]``.

* a flag value (can contain 1,2,3 or 4), raised when the basic and 'KDE' estimates are in disagreement (defined by the ``flag_level`` keyword), and specifically:
   [1] the ``abs(<q>-<q_rs>)/std[q]<=flag_level``

   [2] the ``abs(<q>-<q_rs>)/err[q_rs]<=flag_level``

   [3] the ``abs(<z>-<z_rs>)/std[z]<=flag_level``

   [4] the ``abs(<z>-<z_rs>)/err[z_rs]<=flag_level``

* the number (in %) of the random sample of line fluxes ratios that are located outside of the line ratio diagnostic grids, labelled as ``rs_offgrid``. Line ratio values landing outside of the MAPPINGS IV grid are **NOT** taken into account by pyqz, resulting in a possible **artificial** offset of the ``<q_rs>`` and ``<z_rs>`` values towards the inner part of the q-z plane. ``rs_offgrid>15%`` might therefore indicate trouble. One diagnostic may often be the sole/principal contributor to these off-grid points - a WARNING message is then issued during the processing of the data by ``pyqz.get_qz_ff`` when a specific diagnostic is found to have ``rs_offgrid>15%``. At this time, these individual measurements are **NOT** saved to the final file.



