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



4. Get a local installation of frbpoppy - it will allow to make changes to the frbpoppy code base and having them all instantly available across all your other scripts.
   ::

    $ sudo python3 setup.py develop

   Macs should also be supported, however no tests have been done on Windows.

5. Test whether frbpoppy is working with:
   ::

    $ python3
    >>> import frbpoppy

   If you don't get any errors - hurray, all should be working!

6. Put frbpoppy through its paces:
   ::

    $ python3 examples/_starting_with_frbpoppy_.py

   The first time you run frbpoppy it will generate lookup tables to speed up future runs. This can take up to 2h on a standard 4 core laptop. Subsequent runs will be mere seconds.

****************
How do I use it?
****************
Population synthesis always involves three steps:

1. **Model a population**

   Generating a population in ``frbpoppy`` can be done using the ``generate`` function:
   ::

    from frbpoppy import CosmicPopulation
    pop = CosmicPopulation()

   The generate function takes a number of parameters, allowing you generate the exact population you want, from a steep luminosity function to ultra-long intrinsic pulses. It's all possible. The population you've just generated is called a cosmic population. The next step will be to observe it with whichever survey takes your fancy.

2. **Model a survey**

   In ``frbpoppy``, a population can be observed with a number of different surveys. While many current surveys have been included within ``frbpoppy``, it is possible to define your own survey, or adapt the ones available. To use a survey included in frbpoppy, simply use the code below. Options for adapting the survey, for instance its beam pattern, can also be given as options to the Survey class.

   ::

    from frbpoppy import Survey
    survey = Survey('wrst-apertif')

2. **Create a survey population**

   Use the both the CosmicPopulation and Survey you've set up to create a surveyed population, incorporating all of the selection effects you didn't even know about. Simply

   ::

    from frbpoppy import SurveyPopulation
    surv_pop = SurveyPopulation(pop, survey)

3. **Compare obtained results with actual survey results**
   One of the easiest ways to compare the observed results to actual observations is to use the built in interactive viewer. Simply pass the populations as arguments in to the plot function, and explore away!
   ::

    from frbpoppy import plot
    plot(pop, surv_pop)



*********
What now?
*********
These are the basics, but ``frbpoppy`` offers much more functionality than given in this brief guide. Feel free to pursue the :doc:`documentation`, the `tests <https://github.com/davidgardenier/frbpoppy/tree/master/tests>`_ directory in ``frbpoppy``, or even the `code base <https://github.com/davidgardenier/frbpoppy>`_ itself to find out more.
