.. _faq:

FAQ 
====

1) **Q**: When importing pyqz, I get the following ``WARNING: Statsmodels module not found. KDE_method must be set to 'gauss' or else I will crash.`` What does this mean ?

  **A**: ``pyqz`` could not import the ``statsmodels`` module. See :ref:`troubleshooting` for more details.

------

2) **Q**: ``pyqz`` is returning a lot of different values for log(Q) and 12+log(O/H). What do they all mean and which one should I use ?

   **A**: See :ref:`estimates`

------


3) **Q**: What's with the "Paired" colorbar ? Doesn't that violate all the rules about valid color schemes in scientific diagrams ?

  **A**: Well, yes. But it is also an excellent way of visualizing rapid local changes in the Q and Z plane, e.g. across the PDF associated with specific Q and Z estimates, AND throughout the entire extent of the MAPPINGS V grids - especially when resampled 2 or 3 times. If it hurts your eyes too much, you can easily enable a colorbar of your choosing via ``pyqz_cmap_0`` defined in ``pyqz_metadata.py``.