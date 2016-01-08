
Installing ``pyqz``
===================

The most recent release of ``pyqz`` is available for download from its Github repository: https://github.com/fpavogt/pyqz/releases

Installing ``pyqz`` merely requires to let your Python installation know about its existence. Specifically:

1. Unzip the compressed files (downloaded from the link above), 
2. place the ``pyqz/`` folder anywhere you like, and 
3. add its location to your Python path. 
   
In my case (on MAC OSX), my ``.bash_profile`` would look like : 
   ::

      export PYTHONPATH=$PYTHONPATH:/Users/fvogt/Tools/Python/fpav_pylib/pyqz/


Github users are welcome to fork the ``pyqz`` repository if they want to get access to the latest updates not yet released. *Push requests* for bug fixes and new features are welcome and will be examined in detail. 
      
Requirements
------------
The basic packages below are required for ``pyqz`` to work properly:

* ``numpy`` (1.8.1 or or above)
* ``scipy`` (0.14.0 or above)
* ``matplotlib`` (1.4.2 or above)

Optional (but strongly recommended): 

* ``statsmodels`` (0.6.0 or above)

The ``statsmodel`` package is required to perform the Kernel Density Estimations using ``statsmodel.nonparametric.KDEMultivariate()``, which is often more suitable than the alternative ``scipy.stats.gaussian_kde()``.

Testing the installation
------------------------

First, launch a Python shell and check that you can import ``pyqz``:
::
  
  >>> import pyqz
 
  Multiprocessing possible with: 8 cpus max. 
  Loaded matplotlib with backend: agg
  Loaded pyqz 0.7.2

 
If this fails, then the ``pyqz`` folder is not in your Python path. If the command succeeds, ``pyqz`` will tell how many cpus it can find, which ``matplotlib`` backend is currently being used (critical if you plan on using more than 1 cpu), and finally the code version. 

Next, as a quick test, try to fetch one of the diagnostic grid:
::
	
  >>> a_grid = pyqz.get_grid('[NII]/[SII]+;[OIII]/[SII]+', sampling=1)
  >>> print a_grid	

.. _unittest:

More tests with ``unittest`` 
++++++++++++++++++++++++++++++

A more complete set of tests, relying on the Python ``unittest`` module, are located inside ``pyqz/unittest/``. The interested reader willing to check things further can run them as follows:
::

  >>> cd /path-to-pyqz/pyqz/unittest/
  >>> run pyqz_check.py
 
   Starting pyqz tests:
 
  test01_interpgridnodes (__main__.Testpyqz) ... ok
  test02_interpoffgrid (__main__.Testpyqz) ... ok
  test03_interp_midMVq (__main__.Testpyqz) ...  
  --> Received 1 spectrum ...
  --> Dealing with them one at a time ... be patient now !
   
  All done in 0:00:00.340057
  ok
  test04_get_bad_global_qz (__main__.Testpyqz) ...  
  --> Received 1 spectrum ...
  --> Dealing with them one at a time ... be patient now !
    1: No KDE calculable (bad points/grids ?)
   
  All done in 0:00:00.660051
  ok
  test05_multiprocessing (__main__.Testpyqz) ...  
  --> Received 1 spectrum ...
  --> Launching the multiple processes ... be patient now !
      1 job(s) completed.      
   
  All done in 0:00:01.423516
  ok
  test06_speed_benchmark (__main__.Testpyqz) ...  
  --> Received 1 spectrum ...
  --> Dealing with them one at a time ... be patient now !
   
  All done in 0:00:00.409525
   
  --> Received 1 spectrum ...
  --> Dealing with them one at a time ... be patient now !
   
  All done in 0:00:01.623224
   
  --> Received 1 spectrum ...
  --> Dealing with them one at a time ... be patient now !
   
  All done in 0:00:06.755430
   
  --> Received 1 spectrum ...
  --> Dealing with them one at a time ... be patient now !
   
  All done in 0:00:32.537323
   
  --> Received 24 spectra ...
  --> Launching the multiple processes ... be patient now !
      24 job(s) completed.      
   
  All done in 0:00:41.917625
  ok
  
  ----------------------------------------------------------------------
  Ran 6 tests in 85.806s
  
  OK


Be warned - running these tests will change your ``matplotlib`` backend to a non-interactive one for the current Python session. In fact, the mere fact of importing ``pyqz`` will switch your backend to ``agg``.

.. _troubleshooting:

Troubleshooting
---------------

1) If you get the following message when importing ``pyqz``:
   ::
     WARNING: Statsmodels module not found. KDE_method must be set to 'gauss' or else I will crash.

   then ``pyqz`` could not import the ``statsmodels`` module. This module is required **only if** you want to use the ``KDEMultivariate`` function to construct the joint   probability density function (which we suggest you do). To remove the warning, install ``statsmodels`` and try reloading ``pyqz``. See :ref:`runningpyqzbasic` for more details.

2) If you encounter other errors when importing the module or running the example above, ensure that your ``numpy``, ``scipy`` and ``matplotlib`` packages are up-to-date and try again. 

3) If you still encounter errors after doing all that, check the :ref:`faq`.

4) If the :ref:`faq` doesn't shine some light on your problem, try :ref:`unittest`.

5) Check if this is a known issue: https://github.com/fpavogt/pyqz/issues

6) If you still can't figure out what's wrong, please submit a new issue on the Github repository of the project (https://github.com/fpavogt/pyqz/issues) so we can take a look at it. Please provide as much detail as possible (error message, minimal example able to reproduce the error, operating system, Python version, etc É).

.. warning::
   Submitting a Github issue is the best way for you to get help rapidly, for us to keep track of the problems that need solving, and for future users to see what changes have been made over time (and look for existing solutions to their problem which may be the same as yours). Submitting a new issue on Github is rapid and easy, but if you are really against doing it (why would you ?), you can always email frederic.vogt@alumni.anu.edu.au for help. 

 