
Installing pyqz 
===============

Installing pyqz merely requires to let your Python installation know about its existence. Specifically:

1. Unzip the compressed files, 
2. place the ./pyqz folder anywhere you like, and 
3. add its location to your Python path. 
   
In my case (on MAC OSX), my .bash_profile would   look like : 
   ::

      export PYTHONPATH=$PYTHONPATH:/Users/fvogt/Tools/Python/fpav_pylib/pyqz/
      
Requirements
------------
The basic packages below are required for pyqz to work properly:

* numpy
* scipy
* matplotlib

Optional: 

* statsmodels

pyqz was tested using numpy 1.8.1, scipy 0.14.0, matplotlib 1.4.2 and statsmodel 0.6.0

Testing the installation
------------------------

First, launch a Python shell and check that you can import pyqz:
::
  
  import pyqz
 
If this fails, then the pyqz folder is not in your Python path. 

Next, try to run the example provided with the code:
::

  cd /path-to-pyqz/pyqz/example/	
  pyqz.get_qz_ff(20,'Input.txt', plot = True, savefig='KDE_all')	

This will run the default example included in the pyqz package, generate (a lot) of figures, and save two of them (the KDE ones). 

Troubleshooting
---------------

If you get the following message when importing pyqz:
::
  WARNING: Statsmodels module not found. KDE_method must be set to 'gaussian_kde' or else I will crash.

then pyqz could not find the statsmodels module. This is required ** only if** you want to use the KDEMultivariate function to construct the joint probability density function (which we suggest you do). Install statsmodels and try again. 

If you encounter other errors, check that the required packages listed above are up-to-date, and try again. 

If that doesn't seem to solve the problem, email us at frederic.vogt@anu.edu.au


 