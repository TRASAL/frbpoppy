.. image:: docs/logo_text.png
    :scale: 50

Conduct Fast Radio Burst Population Synthesis | **David Gardenier, Joeri van Leeuwen & Yuyang Wang** | ASTRON & UvA

****************
What's frbpoppy?
****************
Establishing the origin and properties of Fast Radio Bursts (FRBs) is one of the biggest goals in radio astronomy. ``frbpoppy`` was called into life to help astronomers determine which FRB properties to expect. Designed to be simple in use and easy to adapt.

*********************
How can I install it?
*********************
1. Determine if you can and want to install frbpoppy system wide (e.g. with ``sudo``). Otherwise, first create and enter a virtual environment in python with e.g.: 
   ::

    $ python3 -m venv <foo>      # where <foo> is the install directory, e.g. ~/frbpoppy/virtualenv
    $ source <foo>/bin/activate  # where <foo> is the same directory as above
   
2. Ensure ``gfortran`` is installed on your system (used for NE2001).
3. Get the files from the github repository:
   ::

    $ git clone https://github.com/TRASAL/frbpoppy

4. It's important that frbpoppy is installed locally. Why? It means that you can play around with the code without having to dig into a system wide install. It also allows frbpoppy to create data files more easily. Ubuntu is supported, as should be Macs, however no tests have been done on Windows. Install frbpoppy locally by going to the downloaded directory and running
   ::

    $ python3 setup.py develop

5.  Run frbpoppy for the first time. Frbpoppy will automatically create lookup tables for complex calculations. Note this can take up to 2h on a modern machine (4 cores). Subsequent runs will be orders of magnitude faster.
    ::

     $ python3 examples/_starting_with_frbpoppy_.py


******************
How do I use it?
******************
Check out the `examples` directory or the `tests` directory.

****************************************
Which dependencies does `frbpoppy` have?
****************************************
All requirements can be found in `setup.py <https://github.com/TRASAL/frbpoppy/blob/master/setup.py>`_ but are also expanded upon in the following list:

 - `bokeh` for interactive plotting
 - `numpy` for array calculations
 - `pandas` for interactive plotting and easy import of csvs
 - `frbcat` for importing frbcat data
 - `scipy` for Bessel functions and integrations
 - `SQLAlchemy` for creating and querying  cosmological databases
 - `matplotlib` for plotting
 - `requests` for downloading new versions of frbcat
 - `dill` for saving populations in pickled files
 - `tqdm` for nice progress bars during long calculations
 - `joblib` for parallel processing of long calculations

And if using an old version of Python (<v3.6):

 - `future-fstrings` for using f-strings in old Python versions

***************
How can I help?
***************
Spotted a bug, or want to add some functionality? Simply `open an issue <https://github.com/TRASAL/frbpoppy/issues/new>`_ on github, or contact `Joeri van Leeuwen <leeuwen@astron.nl>`_ or `Yuyang Wang <y.wang3@uva.nl>`_ directly. Any help is appreciated!

********************
Who are the authors?
********************
* **David Gardenier**
   During his PhD, David worked at the Netherlands Institute for Radio Astronomy (ASTRON) and the University of Amsterdam (UvA), with the APERTIF/ALERT team, to establish the properties of FRBs. He wrote  frbpoppy.  He works in industry now, drop him a line through `his webpage <https://davidgardenier.github.io/>`_ if you enjoy working with frbpoppy.

* **Yuyang Wang**
   Yuyang currently pursues a PhD at the University of Amsterdam (UvA) and the Netherlands Institute for Radio Astronomy (ASTRON). His research focuses on FRB population studies and pulsar timing. Get in touch with him via `email <y.wang3@uva.nl>`_.

* **Joeri van Leeuwen**
   Joeri is a faculty member at the Netherlands Institute for Radio Astronomy (ASTRON) and conceived frbpoppy.  `Email him <leeuwen@astron.nl>`_ or talk to him at a conference. 

   
