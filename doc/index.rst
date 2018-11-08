.. NetworkEX documentation master file, created by
   sphinx-quickstart on Thu Nov  8 15:04:29 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NetworkEX's documentation!
=====================================

NetworkEX is a simple extension to NetworkX and matplotlib to draw reaction network diagrams from SBML or Antimony strings with ease. NetworkEX also supports drawing weighted reaction network diagrams from an ensemble of models based on reaction frequency or user-supplied weights::

    import networkex as nex

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
    net = nex.NetworkEX(AntimonyStr)
    net.draw()


.. toctree::
   :maxdepth: 2
   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. highlight:: python
