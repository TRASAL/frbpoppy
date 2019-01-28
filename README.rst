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


3. Install frbpoppy locally on Ubuntu by going to the downloaded directory and running:
   ::
   
    $ pip3 install -r requirements.txt
    $ python3 setup.py develop

   Macs should also be supported, however no tests have been done on Windows.

   *If receiving an ASCII error, this is an error from upstream plotting library: contact me on how to fix it.*


******************
How do I use it?
******************
Well, that's where things get more interesting.

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
