# netplotlib
A simple package for visualizing reaction network models

Copyright 2018 Kiri Choi

## Introduction

Netplotlib is a simple extension to [NetworkX](https://networkx.github.io/) and [matplotlib](https://matplotlib.org/) to draw reaction network diagrams from SBML or Antimony strings with ease. Netplotlib also supports drawing weighted reaction network diagrams from an ensemble of models based on reaction frequency or user-supplied weights.

```{python}
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
net = npl.Network(AntimonyStr)
net.draw()
```

## Installation

To install, run the following command:

```
pip install netplotlib
```

## Documentation

Documentation is available at https://netplotlib.readthedocs.io/en/latest

## Examples

<p align="center">
  <img src="/images/mapk.png" width="400px" />
  <img src="/images/weighted.png" width="400px" /> 
</p>
