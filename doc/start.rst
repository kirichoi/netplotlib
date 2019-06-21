===========
Quick Start
===========

To start using netplotlib, `install the package using pip <https://netplotlib.readthedocs.io/en/latest/installation.html>`_ or use `Tellurium <http://tellurium.analogmachine.org/>`_. 

This section demonstrates how to load in a model and draw network diagrams. To start with, lets import netplotlib package and define a model. Here, we use a simple feed forward loop as an example.

.. code-block:: python

    import netplotlib as npl

    AntimonyStr = '''
    $Xi -> S1; k0*Xi
    S1 -> S2; k1*S1
    S2 -> S3; k2*S2
    S1 -> S3; k3*S1
    S3 -> $Xo; k4*S3

    Xi = 3; Xo = 2
    k0 = 0.46; k1 = 0.73; k2 = 0.64;
    k3 = 0.51; k4 = 0.22
    '''
    
Next, create an Network object.

.. code-block:: python

    net = npl.Network(AntimonyStr)
    
To generate network diagrams, simply run `draw() <https://netplotlib.readthedocs.io/en/latest/API.html#netplotlib.Network.draw>`_.

.. code-block:: python

    net.draw()

.. image:: ../images/ffl.png
    :width: 55%



