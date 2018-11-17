# NetworkEX
A simple extension to NetworkX for visualizing reaction network models.

Copyright 2018 Kiri Choi

## Introduction

NetworkEX is a simple extension to [NetworkX](https://networkx.github.io/) and [matplotlib](https://matplotlib.org/) to draw reaction network diagrams from SBML or Antimony strings with ease. NetworkEX also supports drawing weighted reaction network diagrams from an ensemble of models based on reaction frequency or user-supplied weights.

```{python}
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
```

## Documentation

Documentation is available at https://networkex.readthedocs.io/en/latest

## Examples

<p align="center">
  <img src="/images/mapk.png" width="400px" />
  <img src="/images/weighted.png" width="400px" /> 
</p>
