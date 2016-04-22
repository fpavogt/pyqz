.. _understandingpyqz:

Understanding pyqz
======================
For a set of nebular emission line fluxes and errors, pyqz measures the associated value 
of the oxygen abundance 12+log(O/H) and ionization parameters log(Q), given a set of 
MAPPINGS simulations of HII regions. The code uses **flat(-ish)** emission line diagnostic 
grids to disentangle and interpolate the values of log(Q) and 12+log(O/H).  

As pyqz wraps around MAPPINGS simulations, it can provide estimates of the total abundance 
(``Tot[O]+12``) or the gas-phase abundance of the HII region (``gas[O]+12``). In the 
reminder of this document, whenever the former is used, it is understood that it is 
replaceable by the latter.

If you have read this doc from the start, you probably have pyqz installed on your 
machine by now, and managed to run the basic examples described in :ref:`runningpyqzbasic`. 
But before you move on to process your own data, there are a few critical elements that 
you cannot ignore any longer. 

.. warning::

   We're serious here - read this page or be doomed !

A note on the pyqz syntax
------------------------------------------

The pyqz module is intimately linked to the MAPPINGS code. While both are stand alone 
and distinct programs, pyqz was designed to employ the same notation conventions than that 
of the MAPPINGS code for clarity, both from a user and programming perspective. 

These conventions, designed to maximise clarity while minimizing the overall character 
counts, are as follows:

* the ionization parameter is ``LogQ``
* the oxygen abundance is ``Tot[O]+12`` (total) or ``gas[O]+12`` (for the gas-phase)
* the Balmer lines from Hydrogen are ``Ha``, ``Hb``, etc ...
* the main forbidden lines are marked as ``[OIII]``, ``[NII]``, ``[SII]``, ``[OI]``, etc ...
* other strong lines are tagged with their wavelength, i.e. ``4363``, ``3726``, ``3729``, etc ...
* for the usual strong line doublets, when the doublet line fluxes are considered together (
  i.e. [OIII]5007 + 4959), a ``+`` is appended to the said emission line, e.g. ``[OIII]+``. 
  By convention, the single line is always the strongest within the doublet. In short, 
  ``[OIII]`` corresponds to [OIII]5007, ``[OIII]+`` corresponds to [OIII]5007+4959,
  ``[NII]`` corresponds to [NII]6584, ``[NII]+`` corresponds to [NII]6584+6548, etc ...

This syntax must be followed carefully when using pyqz, or errors will arise.
 
The spirit of pyqz
-------------------------

The pyqz module is composed of a core function: ``pyqz.interp_qz``. 
This function is responsible for interpolating the MAPPINGS V grid of simulations of HII 
regions (using ``scipy.interpolate.griddata``) and returns the corresponding value of z or 
q for a given pair of line ratios. This function is basic, in that it does not propagate 
errors on its own. You feed it a pair of line ratio, it returns ``LogQ``, ``Tot[O]+12`` 
or ``gas[O]+12``, and that's it.

The function ``pyqz.get_global_qz`` is a wrapper around ``pyqz.interp_qz``. It is designed 
as a top interaction layer for the ``pyqz`` module, and can propagate errors or 
upper-limits on the line flux measurements. You feed it your measured line fluxes and 
associated errors, and it returns all the ``LogQ`` and ``Tot[O]+12`` or ``gas[O]+12`` 
estimates and associated errors.

Yep, that's right: estimateS. What are these ? 

Direct estimates
""""""""""""""""

pyqz uses a well defined set of line ratio diagnostic grids (the list of which can be 
seen using ``pyqz.diagnostics.keys()``) to interpolate ``LogQ`` and ``Tot[O]+12``. Given 
a set of line fluxes, pyqz can therefore compute 1 estimate of ``LogQ`` and ``Tot[O]+12`` 
per diagnostic diagram chosen by the user, e.g. ``[NII]/[SII]+;[OIII]/[SII]+``. These 
**single direct estimates** (labelled with ``|LogQ`` and ``|Tot[O]+12`` for each 
diagnostic diagram, e.g. ``[NII]/[SII]+;[OIII]/[SII]+|LogQ``) are the most straightforward 
ones computed by pyqz.

Of course, because all line ratio diagnostic grids are constructed from the same set of 
MAPPINGS simulations, all these individual direct estimates ought to be consistent, so 
that computing their mean value is a sensible thing to do. These 
**global direct estimates** are labelled ``<LogQ>``, ``<Tot[O]+12>``, etc. and the 
associated standard deviations are labelled ``std(LogQ)``, ``std(Tot[O]+12)``, etc.

KDE estimates
"""""""""""""

As we do not live in a perfect world, some errors are usually associated with the 
measurement of line fluxes (sigh!). The direct estimates do not take any errors into 
account - the **KDE estimates** (KDE = Kernel Density Estimation) do.

The idea is as follows. First, a set of ``srs`` (where ``srs=400`` is the default) 
random flux values (for each emission line) sampling the probability density function of 
each measurement is generated. Each of these ``srs`` pseudo-sets of line fluxes are fed 
through ``pyqz.interp_qz()``, which returns ``srs`` random estimates of ``LogQ`` and 
``Tot[O]+12``. ``pyqz`` then uses a Kernel Density Estimation tool to reconstruct 

a) the probability density function (PDF) in the ``LogQ`` and ``Tot[O]+12`` plane for 
   every single diagnostic grid selected by the user, and 

b) the full probability density function in the ``LogQ`` and ``Tot[O]+12`` plane resulting 
  from the combination of all ``srs`` estimates for all chosen diagnostic grids. 
  
Python users have the ability to pickle these (individual and global) reconstructed PDFs 
for external use (via the ``KDE_save_PDFs`` keyword). 

From the reconstructed probability density functions, pyqz computes the 0.61% 
(i.e. the :math:`1-{\sigma}` contour for a log normal distribution) level contour in 
the ``LogQ`` vs ``Tot[O]+12`` plane, with respect to the peak. pyqz subsequently 
returns as an (individual or global) KDE estimate the mean of the 0.61% contour and its 
associated half spatial extent along the ``LogQ`` and ``Tot[O]+12`` directions.  

These **single KDE estimates** are referred to (accordingly) using ``|LogQ{KDE}`` and 
``|Tot[O]+12{KDE}`` for the individual diagnostic grids (e.g. 
``[NII]/[SII]+;[OIII]/[SII]+|LogQ{KDE}`` with an error 
``err([NII]/[SII]+;[OIII]/[SII]+|LogQ{KDE})``). The **global KDE estimates** are labelled 
as ``<LogQ{KDE}>`` and ``<Tot[O]+12>``, with associated errors ``err(LogQ{KDE})`` and 
``err(Tot[O]+12{KDE})``.

At this point, things are most likely more confused than ever, and one may be wondering ...

.. _estimates:

What estimates of ``LogQ`` and ``Tot[O]+12`` should one use ?
--------------------------------------------------------------

Unfortunately, there is no definite answer to this question. If all goes well (i.e. your 
measurements are reliable and have reasonable errors), the global KDE estimates 
(``<LogQ{KDE}>`` and ``<Tot[O]+12>``) are the values one should use: these combine all 
requested diagnostic grids estimates and observational errors down to one number. 

But many things can go wrong: one (or more) of your line fluxes might be unknowingly off, 
or perhaps the choice of MAPPINGS simulations is not quite appropriate for the HII regions 
one may be working with (in terms of pressure, abundances, structure, depletion, etc.), or 
perhaps real HII regions may simply not behave quite like MAPPINGS is predicting (sigh!). 

**In all those cases, one must use extreme caution with the global KDE estimates.** A lot 
of information lies in the individual estimates of ``LogQ`` and ``Tot[O]+12``, and 
especially in bad cases. 

So, how does one identify the *good* cases from the *bad* cases ?

Comparing the averaged direct estimates (e.g. ``<LogQ>``) with the global KDE estimates 
(e.g. ``<LogQ{KDE}>``) is a good way to spot problem. For each set of line ratios fed to 
``pyqz.get_global_qz()``, the code checks how similar those estimates are, and issues a 
flag if they are not. The possible flag values are as follows:
  - 9: the PDF is multipeaked. This indicates a likely mismatch between some of the 
       diagnostic grids in their estimates of ``LogQ`` and ``Tot[O]+12``.
  - 8: the observed set of line fluxes is located outside the valid region of one or 
       more of the chosen diagnostic grids.
  - -1: no KDE was computed (either ``srs`` was set to 0, or a line flux errors was 
        set to 0).
  - 1 to 4: these flags are raised when the averaged direct estimates are offset by 
            more than ``flag_level`` times their standard deviations, e.g.:

    * 1 :math:`{\leftrightarrow}` :math:`{|}` ``<LogQ>`` - ``<LogQ{KDE}>`` :math:`{|}` :math:`{<}` ``std(LogQ)`` :math:`{\cdot}` ``flag_level``
    * 2 :math:`{\leftrightarrow}` :math:`{|}` ``<LogQ>`` - ``<LogQ{KDE}>`` :math:`{|}`  :math:`{<}` ``err(LogQ{KDE})`` :math:`{\cdot}` ``flag_level``
    * 3 :math:`{\leftrightarrow}` :math:`{|}` ``<Tot[O]+12>`` - ``<Tot[O]+12{KDE}>`` :math:`{|}` :math:`{<}` ``std(Tot[O]+12)`` :math:`{\cdot}` ``flag_level``
    * 4 :math:`{\leftrightarrow}` :math:`{|}` ``<Tot[O]+12>`` - ``<Tot[O]+12{KDE}>`` :math:`{|}` :math:`{<}` ``err(Tot[O]+12{KDE})`` :math:`{\cdot}` ``flag_level``


Looking at the flags can be helpful in identifying potentially problematic sets of line 
fluxes and (maybe?) the cause. Is one diagnostic grid estimates consistently off ? 
Then maybe some errors in one of the associated line ratio measurements is not properly 
accounted for.

In the end, it remains to the user to decide which estimate(s) to use. The final choice 
will significantly depend on the intended usage, the importance given to the ``LogQ`` and 
``Tot[O]+12`` estimates in a subsequent analysis, and the ability to construct a precise 
model of the said HII region in the first place. 

**It cannot be stressed enough that choosing appropriate HII regions parameters (in terms 
of pressure, spatial structure, abundances, etc.) for the MAPPINGS simulations can and 
will influence the final estimates of ``LogQ`` and ``Tot[O]+12``, both single and global 
ones**. 

If you are using pyqz, chances are that you do not possess enough information 
to define these elements with certainty, and simply use the default diagnostic grids 
provided. This is fine. But in case of estimates mismatch, one must then keep this fact 
in mind.











