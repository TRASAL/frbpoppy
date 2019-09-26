========
frbpoppy
========

|

Conduct Fast Radio Burst Population Synthesis | **David Gardenier** | ASTRON


****************
What's frbpoppy?
****************
Establishing the origin and properties of Fast Radio Bursts (FRBs) is one of the biggest goals in radio astronomy. ``frbpoppy`` was called into life to help astronomers determine which FRB properties to expect. Designed to be simple in use and easy to adapt, ``frbpoppy`` continues the work of ``psrpop`` and ``psrpoppy`` in the realm of FRBs.

**********************
How can I install it?
**********************
.. code-block:: bash

    git clone https://github.com/davidgardenier/frbpoppy
    cd ./frbpoppy/
    python3 setup.py develop
    python3 examples/_starting_with_frbpoppy_.py

Things not going quite as smoothly? Check the :doc:`get_started` guide for additional information on common installation problems.


******************
How do I use it?
******************
Feeling confident? Start off with this:

.. literalinclude:: ./../examples/_starting_with_frbpoppy_.py
   :linenos:

Curious about additional functionality in ``frbpoppy``? Check the :doc:`get_started` guide, the :doc:`documentation`, or the `tests <https://github.com/davidgardenier/frbpoppy/tree/master/tests>`_ directory in ``frbpoppy``.


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
.. toctree::
   :titlesonly:
   :maxdepth: 1

   get_started
   documentation

*****************
Or try:
*****************
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
