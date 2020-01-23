.. image:: docs/logo_text.png
    :scale: 50

Conduct Fast Radio Burst Population Synthesis | **David Gardenier** | ASTRON

****************
What's frbpoppy?
****************
Establishing the origin and properties of Fast Radio Bursts (FRBs) is one of the biggest goals in radio astronomy. ``frbpoppy`` was called into life to help astronomers determine which FRB properties to expect. Designed to be simple in use and easy to adapt, ``frbpoppy`` continues the work of ``psrpop`` and ``psrpoppy`` in the realm of FRBs.

*********************
How can I install it?
*********************
1. Ensure ``gfortran`` is installed on your system (used for NE2001).
2. Get the files from the github repository:
   ::

    $ git clone https://github.com/davidgardenier/frbpoppy

3. It's important that frbpoppy is installed locally. Why? It means that you can play around with the code without having to dig into a system wide install. It also allows frbpoppy to create data files more easily. Ubuntu is supported, as should be Macs, however no tests have been done on Windows. Install frbpoppy locally by going to the downloaded directory and running
   ::

    $ python3 setup.py develop

4.  Run frbpoppy for the first time. Frbpoppy will automatically create lookup tables for complex calculations. Note this can take up to 2h on a modern machine (4 cores). Subsequent runs will be orders of magnitude faster.
    ::

     $ python3 examples/_starting_with_frbpoppy_.py


******************
How do I use it?
******************
Check out the `examples` directory, the `tests` directory, or else ``frbpoppy``'s `webpage <https://davidgardenier.github.io/frbpoppy/>`_!

****************************************
Which dependencies does `frbpoppy` have?
****************************************
All requirements can be found in `setup.py <https://github.com/davidgardenier/frbpoppy/blob/master/setup.py>`_ but are also expanded upon in the following list:

 - `bokeh >= 1.3.4` for interactive plotting
 - `numpy >= 1.17.0` for array calculations
 - `pandas >= 0.23.4` for interactive plotting and easy import of csvs
 - `scipy >= 1.1.0` for Bessel functions and integrations
 - `SQLAlchemy >= 1.3.0` for creating and querying  cosmological databases
 - `matplotlib >= 2.2.3,<3.1` for plotting
 - `requests >= 2.20.0.` for downloading new versions of frbcat
 - `dill >= 0.3.1.1` for saving populations in pickled files
 - `tqdm >= 4.35.0` for nice progress bars during long calculations
 - `joblib >= 0.13.2` for parallel processing of long calculations

And if using an old version of Python (<v3.6):

 - `future-fstrings >= 1.2.0` for using f-strings in old Python versions

***************
How can I help?
***************
Spotted a bug, or want to add some functionality? Simply `open an issue <https://github.com/davidgardenier/frbpoppy/issues/new>`_ on github, or `contact me <gardenier@astron.nl>`_ directly. Any help is appreciated!

********************
Who are the authors?
********************
* **David Gardenier**
   As a PhD student at the Netherlands Institute for Radio Astronomy (ASTRON) and the University of Amsterdam (UvA), I'm working with the APERTIF/ALERT team to establish the properties of FRBs. Get in touch with me via `email <gardenier@astron.nl>`_, drop past my `office <http://davidgardenier.com/#slide=4>`_ or say hello at any `FRB conferences <http://davidgardenier.com/activities.html#slide=3>`_!

*****************
Looking for more?
*****************
Check out ``frbpoppy``'s `webpage <https://davidgardenier.github.io/frbpoppy/>`_ for more information!
