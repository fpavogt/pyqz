.. _faq:

FAQ 
====

1) **Q**: When importing pyqz, I get the following ``WARNING: Statsmodels module not found. KDE_method must be set to 'gauss' or else I will crash.`` What does this mean ?

  **A**: pyqz could not import the statsmodels module. See :ref:`troubleshooting` for more details.

------

2) **Q**: pyqz is returning a lot of different values for ``LogQ`` and ``Tot[O]+12``. What do they all mean and which one should I use ?

   **A**: See :ref:`estimates`

------

3) **Q**: I have my own set of simulations of HII regions. Can I use ``pyqz`` with those ?

   **A**: In principle, yes. pyqz is designed to be fed by an awk script shipped with the MAPPINGS code. But if you make your simulations look like what pyqz is expecting (namely, a suitable filename and file structure matching the ones inside ``pyqz.pyqzm.pyqz_grid_dir``), then you ought to be able to use pyqz just fine with non-MAPPINGS data.

------


4) **Q**: What's with the "Paired" colorbar ? Doesn't that violate all the rules about valid color schemes in scientific diagrams ?

   **A**: Well, yes. But it is also an excellent way of visualizing rapid local changes in the ``LogQ`` and ``Tot[O]+12`` plane, e.g. across the PDF associated with specific ``LogQ`` and ``Tot[O]+12`` estimates, AND throughout the entire extent of the MAPPINGS grids - especially when resampled 2 or 3 times. In any case, you can easily enable a colorbar of your choosing by redefining ``pyqz.pyqzm.pyqz_cmap_0`` to your liking, e.g. ``pyqz.pyqzm.pyqz_cmap_0 = 'Blues'``.
   
------


5) **Q**: My plots look weird after using pyqz ? Why ?

   **A**: To make sure that the label axis are displayed properly, pyqz relies on a custom matplotlib style file (``pyqz_plots.mplstyle``). If you want to revert back to the default style after using pyqz, you can do so via (e.g.): ``plt.style.use('classical')``
   
   
   