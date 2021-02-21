# netplotlib

[![Documentation Status](https://readthedocs.org/projects/netplotlib/badge/?version=latest)](http://netplotlib.readthedocs.org/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/netplotlib.svg)](https://badge.fury.io/py/netplotlib)

A package for visual analysis of biochemical reaction network models

Copyright 2018-2021 Kiri Choi

## Introduction

Netplotlib is an extension to [NetworkX](https://networkx.github.io/) and [matplotlib](https://matplotlib.org/) for visual analysis of biochemical reaction network diagrams written in SBML or Antimony strings. Netplotlib supports visualization of quantities such as flux and species rate of change. Netplotlib provides functions to visualize an ensemble of biochemical reaction networks.

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
