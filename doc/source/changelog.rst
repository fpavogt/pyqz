.. _changelog:

Changelog
==========
v0.7.2 December 2015, F.P.A. Vogt
 - updated documentation with IPython notebook (static HTML in doc + notebook in ``pyqz/examples/``).
 - added ``run_awk_loop`` function to rapidly create the basic grids for pyqz. This is not intended as a main pyqz feature, but rather as an internal tool to make my life easier. It basically runs the MAPPINGS awk sripts in a loop to avoid many (error prone) manual modifications of the rungrid.sh file. Requires MAPPINGS to be installed properly, as well as the awk scripts - neither provided within ``pyqz`` itself.
 - used ``os.path.join()`` for proper handling across different OS
 
v0.7.1 November 2015, F.P.A. Vogt
 - created "get_MVphotogrid_fn" to construct the filename of the MAPPINGS grids once only (for more portability in future updates).
 - added safety check for the diagnostic grids given by the user in "get_global_qz"
 - updated MAPPINGS models to latest version
 - adjusted the resampling function for more consistency with the "resampling factor"
 - updated FAQ with point about the "paired" colormap, to avoid the wrath of the internet ...
 - improved the looping of the grid inside ``interp_qz`` (by checking the presence of data point inside each panel prior to the interpolation) to speed up the process
 - implemented the automatic reshaping of arrays to a 1xN array inside ``interp_qz`` (and reshaping to the original shape prior to returning the results) 
 - updated the "plot regions" assigned to each diagnostic grid for the latest MAPPINGS grids (no scientific meaning - just to get pretty plots)
 - fixed minor typo in ``get_global_qz_ff``: flag is now written as an integer to the file
 - added the ``KDE_QZ_sampling`` and ``KDE_do_singles`` keywords to ``get_global_qz`` for easier access
 - moved ``QZs_lim`` and ``PDF_cont_level`` to ``pyqz_metadata.py`` for easy access
 - re-instated the ability to calculate the individual KDE alongside the global one
 - created custom colormap for the KDE plot
 - added ability to save all the reconstucted PDFs (individuals AND global) to a pickle file via the ``KDE_save_PDF`` keyword
 - fixed bugs when ``srs = 0`` or all errors are 0: the code now doesn't compute any KDEs in those cases, and issues the flag -1
 - fixed minor "aesthetic" sign issue with axes labelling

v0.7.0 August 2015, F.P.A. Vogt
 - added the function ``refine_MVphotogrid``, which uses Akima splines (in 1-D) to resample a given MAPPINGS grid.
 - updated code to allow the use of resampled MAPPINGS grid (or not)
 - updated plots to differentiate between "genuine" MAPPINGS points and resampled nodes

--------

v0.6.3 July-August 2015, F.P.A. Vogt
 - updated input files to MAPPINGS V files (csv created by Awk scripts, incl. ``Pk``, ``kappa``, ``sph/pp``, ...)
 - added functions to check grids
 - started added proper warnings.warn and sys.exit('error message') 
 - implemented new MAPPINGS V terminology for line ratios
 - fixed several minor/less minor bugs
 - added auto-detection of wraps (not bullet-proof, but sufficient for the task at hand)
 - removed pre-defined "valid" grid regions in favor of auto-detection of the good regions
 - added flag ('9') for when multiple KDE peak (best value) exist 
	Note: 	
		- flag 1 => ('direct' mean -'KDE' mean)/ std 'direct' mean <= flag_level, for qz[0]
		- flag 2 => ('direct' mean - 'KDE' mean)/ std 'KDE' mean <= flag_level, for qz[0]
 		- flag 3 => ('direct' mean -'KDE' mean)/ std 'direct' mean <= flag_level, for qz[1]
		- flag 4 => ('direct' mean - 'KDE' mean)/ std 'KDE' mean <= flag_level, for qz[1]
		- flag 8 => the data lands outside of ALL the grids
		- flag 9 => multiple peak to the KDE map exist. Defaulted to highest peak.
 - updated diagrams design
 - started implementing 3-D diagnostics diagrams
 - added ``get_global_qz_ff``, to feed the data in via CSV files
 - moved project to Github, consolidated the Sphinx Documentation as Github pages.
 
v0.6.2 December 2014, F.P.A. Vogt
  - removed the need for the global PYQZ_DIR variable  
v0.6.1 December 2014, F.P.A. Vogt
  - added number of srs points outside the grid (i.e. how reliable the final estimate is)
  - generate random point via truncated normal function to avoid creating negative fluxes
  - implemented upper limits (marked with flux errors = -1)
  - added decent documentation using the ``Sphinx`` module

v0.6.0 November 2014, F. Vogt
  - added proper support of the observational errors on the line fluxes. Done via propagation of probability density function and KDE reconstruction. KDE can be done via either the ``scipy`` package (fast) or the ``statsmodels`` package (very slow) but more accurate for multi-modal distributions.
  - updated saveplot/savefig commands in ``get_qz_ff`` for more freedom
  - added flag to output file when simple mean and advanced 'KDE' mean disagree by n sigmas
  - created new 'KDE' diagram of the q-z plane.
  - fixed output file header to avoid 'spaces' in column names

--------

v0.5.0 November 2013, F.P.A. Vogt
  - implemented a work around in ``get_qz``, to ensure that ``savefig`` works fine also 
    with eps and pdf (note: the real issue is deep. In some cases (e.g. [NII]/[SII] 
    vs [OIII]/Hb, k=50), the original code would save different images in eps/pdf 
    or png ! It is linked to ``facecolor='none'`` in the path. Work-around: 
    ``facecolor='w'``, and ``zorder=0``.
  - ``get_qzff`` can now handle files with only 1 spectrum properly
  - added ``savefig/plot/plot_loc/save_fmt`` keywords to ``get_qzff`` for better  
    control and ability to save figures
  - fixed a bug in ``get_qzff``: the function can now handle 'extra' column 
    containing numbers and/or text

--------

v0.4.1, F. Vogt
  - fixed a wrong plot label for the [NII] line
v0.4, F. Vogt
  - fixed several bugs related with 2D input arrays
  - limited the number of bad points plotted to 1500 
    (for compatibility with ``grid_x`` and ``grid_y``)
  - clarified the required input structure - must be numpy arrays (1D or 2D)
  - plots are now prettier, and more robust (defined via local rcParams)
  - added possibility to save plot via ``savefig`` keyword
  - improved axis labels in the plots to make them 'publication-ready' (if one wanted to)

--------

v0.3.3, F. Vogt
  - now also displays the points landing outside the grid model with white 
    triangles (only for the 1-D array input type)
v0.3.2b, F. Vogt
  - corrected bug related to integer line ratios (e.g. [0],[0]) 
  - corrected bug related to the step checking if line ratios are on the MAPPINGS IV grid
v0.3.2, June 2013, F. Vogt (following suggestions by D. Nicholls)
  - added 'smart' plot limits (instead of fixed ones)
  - increased grid and data point size
  - added new keyword for choosing the plot window number (n_plot)
  - added 'if' statement to close the plot if all values are NaNs 
    (removed in v0.4)
  - added plot title
v0.3.1b, April 2013, F. Vogt
  - fixed indentation of 4 lines in ``get_qzff``
v0.3.1 April 2013, F. Vogt
  - added .csv output for the ``get_qzff`` (or txt, which ever you like best)
  - changed header column of output files (only 'z' is used for consistency)
v0.3.0 April 2013, F.P.A. Vogt
  - added ``get_qzff`` function to directly work from a txt file
  - corrected ``get_pyqz`` for when points are 'on' the grid.
v0.2.0 April 2013, F.P.A. Vogt
 - modified fitting method to be 'slice-by-slice' for smoother results
 - added the get_grid function
 - added different readable areas for different grids and kappas
v0.1.0 Feb. 2013, F. Vogt
 - created

 

 
  
 
