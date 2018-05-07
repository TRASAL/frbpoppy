===========
Get Started
===========

**********************
How can I install it?
**********************
1. Ensure ``gfortran`` is installed on your system (used for NE2001).



2. Get the files from the github repository:
   ::

    $ git clone https://github.com/davidgardenier/frbpoppy



3. Install frbpoppy locally on Ubuntu by going to the downloaded directory
   ::

     $ cd ./frbpoppy


4. Run the following to get a local install - it will allow to make changes to the frbpoppy code base and having them all instantly available across all your other scripts.
   ::

    $ sudo python3 setup.py develop

   Macs should also be supported, however no tests have been done on Windows.

   *If receiving an ASCII error, this is* `a bug <https://github.com/bokeh/bokeh/issues/7272>`_ *in Bokeh.*

5. Test whether frbpoppy is working with:
   ::

    $ python3
    >>> import frbpoppy

   If you don't get any errors - hurray, all should be working!


****************
How do I use it?
****************
Population synthesis always involves three steps:

1. **Model population**

   Generating a population in ``frbpoppy`` can be done using the ``generate`` function:
   ::

    from frbpoppy.do_populate import generate
    pop = generate()

   The generate function takes a number of parameters, allowing you generate the exact population you want, from a steep luminosity function to ultra-long intrinsic pulses. It's all possible. The population you've just generated is called your initial population. The next step will be to observe it with whichever survey takes your fancy.

2. **Survey modeled population**

   In ``frbpoppy``, a population can be observed with a number of different surveys. While many current surveys have been included within ``frbpoppy``, it is possible to define your own survey, or adapt the ones available. To observe a population simply add the following lines

   ::

    from frbpoppy.do_survey import observe
    surv_pop = observe(pop, 'APERTIF')


3. **Compare obtained results with actual survey results**
   One of the easiest ways to compare the observed results to actual observations is to use the built in interactive viewer. Simply pass the populations as arguments in to the plot function, and explore away!
   ::

    from frbpoppy.do_plot import plot
    plot(pop, surv_pop)



*********
What now?
*********
These are the basics, but ``frbpoppy`` offers much more functionality than given in this brief guide. Feel free to pursue the :doc:`documentation`, the `tests <https://github.com/davidgardenier/frbpoppy/tree/master/tests>`_ directory in ``frbpoppy``, or even the `code base <https://github.com/davidgardenier/frbpoppy>`_ itself to find out more.
