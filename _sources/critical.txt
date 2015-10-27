.. _critical:

pyqz Do's and Dont's
=========================


.. _estimates:

What estimates of q and z should I use ?
""""""""""""""""""""""""""""""""""""""""""""

[Short answer:] In principle, ``<q_rs>``, ``err[q_rs]``, ``<z_rs>`` and ``err[z_rs]``.

[Long answer:] In principle, ``<q_rs>``, ``err[q_rs]``, ``<z_rs>`` and ``err[z_rs]``. When all goes well (reasonable errors, not close from the edge of the diagnostic grids, consistent diagnostics), these values contain the error associated both with the observational errors and the mismatch between the different MAPPINGS IV diagnostic grids. 

When things get difficult (large observational errors, significant mismatch between the diagnostics, etc…), ``<q_rs>``, ``err[q_rs]``, ``<z_rs>`` and ``err[z_rs]`` **may or may not** be appropriate to use. In these situations, ``<q_rs>``, ``err[q_rs]``, ``<z_rs>`` and ``err[z_rs]`` **could** be a safer bet, but it should then kept in mind that their errors solely account for the mismatch between diagnostics (i.e. no observational errors). 

.. NOTE:: Through the ``flags`` and ``rs_offgrid`` values, pyqz will let you spot these difficult situations, but it is ultimately up to the user to then decide which estimates (if any) are the most appropriate.
