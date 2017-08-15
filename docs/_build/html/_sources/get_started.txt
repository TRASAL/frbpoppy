===========
Get Started
===========

**********************
How can I install it?
**********************
1. Get the files from the github repository:
   ::

    $ git clone https://github.com/davidgardenier/frbpoppy

2. Install it on your system:
   ::

    $ python setup.py install

****************
How do I use it?
****************
Population synthesis always involves three steps:

1. **Model population**

   Generating a population in ``frbpoppy`` can be done using the ``generate`` function:
   ::

    from frbpoppy.do_populate import generate
    pop = generate()

   The generate function takes a number of parameters, allowing you generate the exact

2. **Survey modeled population**

   In ``frbpoppy``, a population can be observed with a number of different surveys. While many current survey parameters have been included within ``frbpoppy``, it is possible to define your own survey parameters.

3. **Compare obtained results with actual survey results**

   This is left as an exercise for the reader



*****
Oops!
*****
It seems this page is still under construction...

***********
What to do?
***********
This website is being built to guide beginners in the use of frbpoppy. If you've come this far - you must be an expert! If you arrived here wanting to dig deeper into frbpoppy, you might find the `github repository <https://github.com/davidgardenier/frbpoppy>`_ to be of help. All code should be well documented, and if there's something amiss, just `open an issue <https://github.com/davidgardenier/frbpoppy/issues/new>`_ or `drop me a note <gardenier@astron.nl>`_!
