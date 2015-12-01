
Installing pyqz
===================

The most recent release of ``pyqz`` is available for download from its Github repository: https://github.com/fpavogt/pyqz/releases

Installing ``pyqz`` merely requires to let your Python installation know about its existence. Specifically:

1. Unzip the compressed files (downloaded from the link above), 
2. place the ``./pyqz`` folder anywhere you like, and 
3. add its location to your Python path. 
   
In my case (on MAC OSX), my ``.bash_profile`` would look like : 
   ::

      export PYTHONPATH=$PYTHONPATH:/Users/fvogt/Tools/Python/fpav_pylib/pyqz/
      
Requirements
------------
The basic packages below are required for pyqz to work properly:

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
 
If this fails, then the ``pyqz`` folder is not in your Python path. 

Next, you ought to make sure that you are using the ``pyqz`` you intended:
::
  >>> pyqz.__version__
  '0.7.0'

Then, try to run the example provided with the code:
::

  >>> cd /path-to-pyqz/pyqz/example/	
  >>> pyqz.get_qz_ff(20,'Input.txt', plot = True, savefig='KDE_all')	

This will run the default example, generate (a lot) of figures, and save two of them (the KDE ones). 

.. _troubleshooting:

Troubleshooting
---------------

1) If you get the following message when importing ``pyqz``:
   ::
     WARNING: Statsmodels module not found. KDE_method must be set to 'gauss' or else I will crash.

   then ``pyqz`` could not import the ``statsmodels`` module. This module is required **only if** you want to use the ``KDEMultivariate`` function to construct the joint   probability density function (which we suggest you do). To remove the warning, install ``statsmodels`` and try reloading ``pyqz``. See :ref:`runningpyqz` for more details.

2) If you encounter other errors when importing the module or running the example above, ensure that your ``numpy``, ``scipy`` and ``matplotlib`` packages are up-to-date and try again. 

3) If you still encounter errors after doing all that, check the :ref:`faq`.

4) If you still can't figure out what's wrong, please submit a new issue on the Github repository of the project (https://github.com/fpavogt/pyqz/issues) so we can take a look at it. Please provide as much detail as possible (error message, operating system, Python version, etc É).

.. warning::
   Submitting a Github issue is the best way for you to get help rapidly, for us to keep track of the problems that need solving, and for future users to see what changes have been made over time (and look for existing solutions to their problem which may be the same as yours). Submitting a new issue on Github is rapid and easy, but if you are really (really!) against doing it, you can always email frederic.vogt@alumni.anu.edu.au for help. 

 