.. _faq:

FAQ 
====

1) **Q**: When importing pyqz, I get the following ``WARNING: Statsmodels module not found. KDE_method must be set to 'gauss' or else I will crash.`` What does this mean ?

  **A**: pyqz could not import the statsmodels module. See :ref:`troubleshooting` for more details.

------

2) **Q**: pyqz is returning a lot of different values for ``LogQ`` and ``Tot[O]+12``. What do they all mean and which one should I use ?

   **A**: See :ref:`estimates`

------

3) **Q**: I don't see any matplotlib/pyplot figures anymore after running pyqz. What is happening ?

   **A**: At this time, pyqz *absolutely requires* a non-interactive backend for maplotlib: 
          pyqz forces this to be ``agg``. This non-interactive backend disables 
          all plot windows, but you will still be able to save them to file. To revert to 
          your default backend, the safest is to restart your Python shell.  You can read 
          more about backends here: http://matplotlib.org/faq/usage_faq.html#what-is-a-backend 
          
          We hope to remove this <<feature>> some time in the future.

------

4) **Q**: I have my own set of simulations of HII regions. Can I use ``pyqz`` with those ?

   **A**: In principle, yes. pyqz is designed to be fed by an awk script shipped with the 
          MAPPINGS code. But if you make your simulations look like what pyqz is expecting 
          (namely, a suitable filename and file structure matching the ones inside 
          ``pyqz/reference_data/``), then you ought to be able to use pyqz just fine with 
          non-MAPPINGS data.

------


5) **Q**: What's with the "Paired" colorbar ? Doesn't that violate all the rules about 
          valid color schemes in scientific diagrams ?

  **A**: Well, yes. But it is also an excellent way of visualizing rapid local changes in 
         the ``LogQ`` and ``Tot[O]+12`` plane, e.g. across the PDF associated with specific 
         ``LogQ`` and ``Tot[O]+12`` estimates, AND throughout the entire extent of the 
         MAPPINGS grids - especially when resampled 2 or 3 times. In any case, you can 
         easily enable a colorbar of your choosing via ``pyqz_cmap_0`` defined in 
         ``pyqz_metadata.py``.