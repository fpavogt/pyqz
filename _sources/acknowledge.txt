.. |DOI_0.8.0| image:: https://zenodo.org/badge/doi/10.5281/zenodo.53201.svg
   :target: http://dx.doi.org/10.5281/zenodo.53201
   
.. |DOI_0.7.1| image:: https://zenodo.org/badge/doi/10.5281/zenodo.34502.svg
   :target: http://dx.doi.org/10.5281/zenodo.34502 
 


Acknowledging pyqz
====================

1) Only use lower case letters when mentioning pyqz, and always include the version number.
   Ideally, you should also include the DOI (Digital Object Identifier) associated with any 
   of the Github release, e.g.:

   - pyqz 0.8.0: |DOI_0.8.0|

   - pyqz 0.7.1: |DOI_0.7.1|
   
-----
   
2) pyqz will be described in detail in 

    **Vogt et al.**, in prep.

   and the associated MAPPINGS V models in

    **Sutherland et al.**, in prep.

   If you find pyqz useful for your research, please cite these references accordingly.

-----

3) Additional references you may wish to cite as well include:

   - the article describing the new projected 3-D diagnostic combining Ha, [NII], [SII]+, [OIII] and Hb:

    **Dopita, Kewley, Sutherland, & Nicholls**, *Chemical abundances in high-redshift galaxies: a powerful new emission line diagnostic*, Ap&SS, 361, 61 (2016). `ADS entry <http://adsabs.harvard.edu/abs/2016Ap%26SS.361...61D>`_ 

   - the article first introducing pyqz 

    **Dopita, Sutherland, Nicholls, Kewley & Vogt**, *New Strong Line Abundance Diagnostics for HII Regions: Effects of Updated Atomic Data and kappa-Distributed Electron Energies*, ApJS, 208, 10 (2013). `ADS entry <http://adsabs.harvard.edu/abs/2013ApJS..208...10D>`_ 

   - the initial description of the updated error propagation mechanism (first implemented in pyqz 0.6.0), described in Appendix B of

    **Vogt, Dopita, Borthakur, Verdes-Montenegro, Heckman, Yun & Chambers**, *Galaxy interactions in compact groups - II. Abundance and kinematic anomalies in HCG 91c*, MNRAS, 450, 2593 (2015). `ADS entry <http://adsabs.harvard.edu/abs/2015MNRAS.450.2593V>`_

-----

4) pyqz also uses several packages that **should also be acknowledged in their own right.** 
   The following Tex-formatted acknowledgment is one way to do so::

    This research has made use of \textsc{pyqz} (Dopita et al., 2013), a Python module to 
    derive the ionization parameter and oxygen abundance of HII regions from their strong 
    emission line ratios hosted at \url{http://http://fpavogt.github.io/pyqz}. \textsc{pyqz} 
    relies on \textsc{statsmodel} (Seabold & Perktold 2010) and \textsc{matplotlib} 
    (Hunter 2007).

   References:
    - `Hunter (2007) <http://cdsads.u-strasbg.fr/abs/2007CSE.....9...90H>`_    
    - Seabold & Perktold (2010)::
 
    @inproceedings{seabold2010,
        title={Statsmodels: Econometric and statistical modeling with python},
        author={Seabold, Skipper and Perktold, Josef},
        booktitle={9th Python in Science Conference},
        year={2010},
    }