.. _runningpyqz:

Running ``pyqz``
================

This page is also available as an IPython notebook (``pyqz_demo.ipynb``)
located in ``pyqz/examples/``.

.. note::

    the code syntax in v0.7.0 has changed significantly, and so did the function calls. ``pyqz 0.7.x`` is therefore NOT backward  compatible with older ``pyqz`` versions.

.. warning::
    the examples below will show you how to run the main functions inside ``pyqz``. But these do not exempt you from getting acquainted with the "Understanding ``pyqz``" section of the documentation !

Installing and importing ``pyqz``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installing ``pyqz`` is easy. Download it from Github, unpack it anywhere you
like, and make sure that this location is in your Python path. You
should then be able to import the package and check its version from
within any Python shell:

::

    >>> import pyqz
    >>> pyqz.__version__
    '0.7.2'


Basic use 1: accessing MAPPINGS line ratio diagnostic grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``pyqz`` gives you easy access to the latest MAPPINGS strong line ratio
diagnostic diagrams (and associated info on the MAPPINGS version used to
generate the grid, etc...). This can for example be useful to create
your own line ratio diagnostic plots. You can access the nodes of any
line ratio diagram using ``pyqz.get_grid()``:

::

    >>> a_grid = pyqz.get_grid('[NII]/[SII]+;[OIII]/[SII]+', sampling=1)
    >>> # Uncomment the line below for a look at the structure of a_grid
    >>> # a_grid

The main parameters of the MAPPINGS simulations can be specified via the
following keywords: - Pk let's you define the pressure of the simulated
HII regions, - struct allows you to choose between plane-parallel ('pp')
and spherical ('sph') HII regions, and - kappa lets you define the value
of :math:`\kappa` (from the so-called :math:`\kappa`-distribution).

All these values must match an existing set of MAPPINGS simulations
inside the ``pyqz/reference_data/`` folder, or ``pyqz`` will issue an error. In
other words, ``pyqz`` will not be running new MAPPINGS simulations for you.

So, if one wanted to access the MAPPINGS simulations for plane-parallel
HII regions, with Maxwell-Boltzmann electron density distribution, ``Pk
=5.0``, one should type:

::

    >>> a_grid = pyqz.get_grid('[NII]/[SII]+;[OIII]/[SII]+', struct = 'pp', Pk = 5, kappa = 'inf')
    >>> # Uncomment the line below for a look at the structure of a_grid
    >>> # a_grid

If you want to simply check how a given line ratio diagnostic diagram
looks (and e.g. check whether the MAPPINGS grid is flat, or wrapped) for
line ratios of your choice, you can use ``pyqz.check_grid()``:

::

    >>> bad_segments = pyqz.check_grid('[NII]/[OII]+;[OIII]/[SII]+', show_plot=True)


.. image:: output_11_0.png
   :align: center

An important feature of ``pyqz`` is the auto-detection of wraps in the
diagnostic grids, marked with red segments in the diagram, and returned
as an array by the function ``pyqz.check_grid()``.

By default, the default MAPPINGS grids shipped with ``pyqz`` are corse. For
various reasons better explained elsewhere (see the MAPPINGS
documentation), only a few abundance values have matching stellar tracks
AND stellar atmospheres. Hence, only a few abundance points can be
simulated in a consistent fashion.

Rather than 1) interpolating between stellar tracks and stellar
atmospheres in the abundance space and 2) running extra MAPPINGS models
(which would use inconsistent & interpolated input), ``pyqz`` can directly
resample each diagnostic grid (using the function
``pyqz.refine_MVphotogrid()``, see the docs for more info). The resampling
is performed in the {``LogQ`` and ``Tot[O+12]`` vs line ratio} space for all
line ratios returned by MAPPINGS using Akima splines. Resampled grids
can be accessed via the sampling keyword. Diagnostic grids resampled 2x2
times are shipped in the default ``pyqz`` package and are directly
accessible, e.g.:

::

    >>> bad_segments = pyqz.check_grid('[NII]/[OII]+;[OIII]/[SII]+',show_plot=True, sampling=2)


.. image:: output_13_0.png
   :align: center

In the default ``pyqz`` diagrams, the original MAPPINGS nodes are circled
with a black outline, while the reconstructed nodes are not. For grids
more densely resampled, see the "Advanced use 2" below.

Basic use 2: deriving ``LogQ`` and ``Tot[O+12]`` for a given set of line ratios
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the core of ``pyqz`` lies ``pyqz.interp_qz()``, which is the basic routine
used to interpolate a given line ratio diagnostic grid. The function is
being fed by line ratios stored inside numpy arrays, and will only
return a value for line ratios landing on valid and un-wrapped regions
of the grid:

::

    >>> z = pyqz.interp_qz('Tot[O]+12',[np.array([-0.6]),np.array([-0.1])],'[NII]/[OII]+;[OIII]/[SII]+', sampling=1,struct='pp', show_plot=True)
    >>> z
    array([ 8.697449])

.. image:: output_17_0.png
   :align: center

Of course, one usually wants to compute both ``LogQ`` and ``Tot[O+12]`` or
``gas[O+12]`` for a large set of strong emission line fluxes, combining the
estimates from different line ratio diagnostics diagrams. This is
exactly what the function ``pyqz.get_global_qz()`` allows you to do.

The function is being fed the individual line fluxes and associated
errors in the form of numpy arrays and lists. ID tags for each dataset
can also be given to the function (these are then used if/when saving
the different diagrams to files).

::

    >>> pyqz.get_global_qz(np.array([[2.88,0.05,1.0,0.01,2.015,0.045,0.036,0.002,2.72,0.12,0.497,0.017,0.500,0.016]]),
        ['[OII]+','std[OII]+','Hb','stdHb','[OIII]','std[OIII]','[OI]','std[OI]','Ha','stdHa','[NII]','std[NII]','[SII]+','std[SII]+'], ['[NII]/[SII]+;[OIII]/[SII]+','[NII]/[OII]+;[OIII]/[OII]+'], 
        ids = ['NGC_1234'],
        show_plot='KDE', # set this to True to also see the individual line ratio diagrams
        save_plot=True,  # set this to 'grids', 'KDE_all', False or 'KDE_flags' (i.e. only the problematic points) 
        plot_loc = './example_plots', 
        struct='pp',
        sampling=1)
    --> Processing 1 spectrum ...
    All done in 0:00:04.370502
    [array([[  7.30563675e+00,   8.40752058e+00,   7.58464098e+00,
               8.62834536e+00,   7.44513886e+00,   1.39502113e-01,
               8.51793297e+00,   1.10412390e-01,   7.30007583e+00,
               8.69737797e-03,   8.39739984e+00,   1.69106603e-02,
               7.58765550e+00,   1.58882425e-02,   8.62932458e+00,
               1.19835520e-02,   7.58775759e+00,   3.86829996e-02,
               8.63037803e+00,   3.11968890e-02,   9.24000000e+02,
               0.00000000e+00]]),
     ['[NII]/[SII]+;[OIII]/[SII]+|LogQ',
      '[NII]/[SII]+;[OIII]/[SII]+|Tot[O]+12',
      '[NII]/[OII]+;[OIII]/[OII]+|LogQ',
      '[NII]/[OII]+;[OIII]/[OII]+|Tot[O]+12',
      '<LogQ>',
      'std(LogQ)',
      '<Tot[O]+12>',
      'std(Tot[O]+12)',
      '[NII]/[SII]+;[OIII]/[SII]+|LogQ{KDE}',
      'err([NII]/[SII]+;[OIII]/[SII]+|LogQ{KDE})',
      '[NII]/[SII]+;[OIII]/[SII]+|Tot[O]+12{KDE}',
      'err([NII]/[SII]+;[OIII]/[SII]+|Tot[O]+12{KDE})',
      '[NII]/[OII]+;[OIII]/[OII]+|LogQ{KDE}',
      'err([NII]/[OII]+;[OIII]/[OII]+|LogQ{KDE})',
      '[NII]/[OII]+;[OIII]/[OII]+|Tot[O]+12{KDE}',
      'err([NII]/[OII]+;[OIII]/[OII]+|Tot[O]+12{KDE})',
      '<LogQ{KDE}>',
      'err(LogQ{KDE})',
      '<Tot[O]+12{KDE}>',
      'err(Tot[O]+12{KDE})',
      'flag',
      'rs_offgrid']]


.. image:: output_19_1.png
   :align: center


By default, all line fluxes errors are assumed to be gaussian, where the
input std value is the 1 standard deviation. Alternatively, line fluxes
can be tagged as upper-limits by setting their errors to -1.

Users less keen on using Python extensively can alternatively feed their
data to ``pyqz`` via an appropriately structured .csv file and receive
another .csv file in return:

::

    >>> pyqz.get_global_qz_ff('./example_input.csv', 
                          ['[NII]/[SII]+;[OIII]/[SII]+','[NII]/[OII]+;[OIII]/[OII]+'], 
                          show_plot='KDE', # set this to True to also see the individual line ratio diagrams
                          save_plot=True,  # set this to 'grids', 'KDE_all', False or 'KDE_flags' (i.e. only the problematic points) 
                          plot_loc = './example_plots', 
                          struct='pp',
                          sampling=1)
    --> Processing 1 spectrum ...
    All done in 0:00:04.381982
   (array([[  7.30563675e+00,   8.40752058e+00,   7.58464098e+00,
               8.62834536e+00,   7.44513886e+00,   1.39502113e-01,
               8.51793297e+00,   1.10412390e-01,   7.30009089e+00,
               9.35766576e-03,   8.39731789e+00,   1.43132660e-02,
               7.57782295e+00,   3.05980853e-02,   8.62299975e+00,
               1.85412415e-02,   7.57794077e+00,   5.03249097e-02,
               8.62192008e+00,   3.87429978e-02,   9.24000000e+02,
               0.00000000e+00]]),
     ['[NII]/[SII]+;[OIII]/[SII]+|LogQ',
      '[NII]/[SII]+;[OIII]/[SII]+|Tot[O]+12',
      '[NII]/[OII]+;[OIII]/[OII]+|LogQ',
      '[NII]/[OII]+;[OIII]/[OII]+|Tot[O]+12',
      '<LogQ>',
      'std(LogQ)',
      '<Tot[O]+12>',
      'std(Tot[O]+12)',
      '[NII]/[SII]+;[OIII]/[SII]+|LogQ{KDE}',
      'err([NII]/[SII]+;[OIII]/[SII]+|LogQ{KDE})',
      '[NII]/[SII]+;[OIII]/[SII]+|Tot[O]+12{KDE}',
      'err([NII]/[SII]+;[OIII]/[SII]+|Tot[O]+12{KDE})',
      '[NII]/[OII]+;[OIII]/[OII]+|LogQ{KDE}',
      'err([NII]/[OII]+;[OIII]/[OII]+|LogQ{KDE})',
      '[NII]/[OII]+;[OIII]/[OII]+|Tot[O]+12{KDE}',
      'err([NII]/[OII]+;[OIII]/[OII]+|Tot[O]+12{KDE})',
      '<LogQ{KDE}>',
      'err(LogQ{KDE})',
      '<Tot[O]+12{KDE}>',
      'err(Tot[O]+12{KDE})',
      'flag',
      'rs_offgrid'])


.. image:: output_22_1.png
   :align: center

   
The first line of the input file must contain the name of each column,
following the ``pyqz`` convention. The order itself does not matter, e.g.:

::

    Id,[OII]+,std[OII]+,Hb,stdHb,[OIII],std[OIII],[OI],std[OI],Ha,stdHa,[NII],std[NII],[SII]+,std[SII]+

The Id (optional) can be used to add a tag (i.e. a string) to each set
of line fluxes. This tag will be used in the filenames of the diagrams
(if some are saved) and in the output ``.csv`` file as well.

Commented line begin with #, missing values are marked with
\ :math:`$`\  (set with the missing_values keyword), and the decimal
precison in the output file is set with decimals (default=5).

At this point, it must be stressed that ``pyqz.get_global_qz()`` only
exploit a finite set of diagnostic grids, namely:

::

    >>> pyqz.diagnostics.keys()
    ['[NII]/[OII]+;[OIII]/[SII]+',
     '[NII]/Ha;[OIII]/Hb',
     '[NII]/[OII]+;[OIII]/[OII]+',
     '[NII]/[SII]+;[OIII]/Hb',
     '[NII]/[SII]+;[OIII]/[OII]+',
     '[NII]/[SII]+;[NII]/Ha;[OIII]/Hb',
     '[NII]/[OII]+;[OIII]/Hb',
     '[NII]/Ha;[OIII]/[OII]+',
     '[NII]/[SII]+;[OIII]/[SII]+']


These specific diagnostic diagrams are chosen to be largely flat, i.e.
they are able to cleanly disentangle the influence of ``LogQ`` and
``Tot[O]+12``. One does not need to use all the grids together. For example,
if one knows that an [OII] line flux measurement is corrupted, one ought
to simply use the diagnostic grids that do not rely on this line to
derive the estimates of ``LogQ`` and ``Tot[O]+12``.

Users can easily add new diagnostics to this list (defined inside
``pyqz_metadata.py``), but will do so at their own risk.

Advanced use 1: using custom MAPPINGS grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While ``pyqz`` ships with a default set of HII region simulations from
MAPPINGS, some (all!) users might be interested in using ``pyqz`` with their
own specific sets of MAPPINGS simulations. ``pyqz`` was designed to be
compatible with the grids generated from the awk script provided
alongside MAPPINGS.

If one uses the awk script to create new MAPPINGS grids, the resulting
.csv file must be laced inside ``pyqz/reference_data``. The filename must
match what the function ``pyqz.get_MVphotogrid_fn()`` expects for your
given set of parameters:

::

    >>> pyqz.get_MVphotogrid_fn(Pk=6.7,calibs='GCZO', kappa =10, struct='pp')
    '/Users/fvogt/Tools/Python/fpav_pylib/pyqz/pyqz_dev/reference_data/grid_QZ_pp_GCZO_Pk67_k10.csv'

If one does not use the awk script to generate the custom MAPPINGS grid,
then just make sure your model grid matches the format of existing model
grids located in ``pyqz/reference_data/`` ...

Advanced usage 2: resampling the original MAPPINGS grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, 2 times resampled MAPPINGS grid are shipped with ``pyqz``. These
are generated using the function ``pyqz.resample_MVphotogrid()``, which is
straightforward to use:

::

    >>> grid_fn = pyqz.get_MVphotogrid_fn(Pk=5.0,struct='sph', kappa=np.inf)
    >>> pyqz.resample_MVphotogrid(grid_fn, sampling=2)
    Success: grid_QZ_sph_GCZO_Pk50_kinf.csv resampled by a factor 2x2 and saved as grid_QZ_sph_GCZO_Pk50_kinf_samp_2.csv


More densely resampled grids can then easily be created by varying the
sampling keyword.

Advanced usage 3: "3-D" line ratio diagrams
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``pyqz`` does support 2-D line ratio diagrams constructed from 3 sets of
line ratios (i.e. 3-D line ratio diagrams projected to a given 2-D
plane):

::

    >>> bad_segments = pyqz.check_grid('[NII]/[SII]+;[NII]/Ha;[OIII]/Hb',
                    coeffs = [[1.0,0.264,0.0],[0.242,-0.910,0.342]],
                    show_plot=True,
                    struct='sph',
                    sampling=1)


.. image:: output_37_0.png
   :align: center


