
Installing pyqz
===================

pyqz is available on pypi, which makes its installation easier than ever. 
In a terminal, type:
::
   pip install pyqz


And that should take care of things.

The most recent release of pyqz is also available for download from its `Github repository <https://github.com/fpavogt/pyqz/releases/latest/>`_. 
Interested users can fork the pyqz repository if they want to get access to the 
latest updates not yet released. *Push requests* for bug fixes and new features are 
welcome and will be examined in detail. 
      
Requirements
------------
The following packages are required for pyqz to work properly:

* numpy (1.12.1 or or above)
* scipy (0.19.0 or above)
* matplotlib (2.0.1 or above)

Optional (but strongly recommended): 

* statsmodels (0.6.1 or above)

The statsmodel package is required to perform the Kernel Density Estimations using 
``statsmodel.nonparametric.KDEMultivariate()``, which is often more suitable than the 
alternative ``scipy.stats.gaussian_kde()``.

Testing the installation
------------------------

First, launch a Python shell and check that you can import pyqz, and that it is the intended version:
::
  
  >>> import pyqz
  >>> print pyqz.__version__
 
Next, as a quick test, try to fetch one of the diagnostic grid:
::

  >>> a_grid = pyqz.get_grid('[NII]/[SII]+;[OIII]/[SII]+', sampling=1)
  >>> print a_grid

.. _unittest:

More tests with unittest 
++++++++++++++++++++++++++++++

A more complete set of tests, relying on the Python unittest module, are also available. 
The user willing to run them can do so as follows:
::

  >>> import pyqz.tests
  >>> pyqz.tests.run_all_tests(cleanup=True)
  
  
Note that the last test takes several seconds (~130s or so) to complete. It is
designed to test the functions of pyqz_plots. If it succeeds, the demonstration plots will
be deleted, unless the tests are run with ``cleanup=False``. In that case, the plots will
be stored in: 
::

   >>> print pyqz.tests.arena

.. _troubleshooting:

Troubleshooting
---------------

1. 
If you get the following message when importing pyqz:
::

    WARNING: Statsmodels module not found. KDE_method must be set to 'gauss' or else I will crash.

then pyqz could not import the statsmodels module. This module is required **only if** 
you want to use the ``KDEMultivariate`` function to construct the joint probability 
density function (which we suggest you do). To remove the warning, install statsmodels 
and try reloading pyqz.

2. 
If you encounter other errors when importing the module or running the example above, 
ensure that your numpy, scipy and matplotlib packages are up-to-date and try again. 
  
3. 
If you still encounter errors after doing all that, check the :ref:`faq`.
  
4. 
If the :ref:`faq` doesn't shine some light on your problem, try :ref:`unittest`.
  
5. 
Check if this is a known issue: https://github.com/fpavogt/pyqz/issues
  
6. 
If you still can't figure out what's wrong, please `submit a new issue on the Github 
repository of the project <https://github.com/fpavogt/pyqz/issues>`_. Provide as much detail as possible (error message, minimal example able to reproduce the error, operating system, Python version, etc ...).

.. note::
   Submitting a Github issue is the best way for you to get help rapidly, for us to keep 
   track of the problems that need solving, and for future users to see what changes have 
   been made over time (and look for existing solutions to their problem which may be the 
   same as yours). Submitting a new issue on Github is rapid and easy, but if you are 
   really against doing it (why would you ?), you can always email 
   frederic.vogt@alumni.anu.edu.au for help. 

 