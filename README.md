[![PyPI version](https://badge.fury.io/py/bravais.svg)](https://badge.fury.io/py/bravais)

# bravais
A simple package for representing [Bravais lattices](https://en.wikipedia.org/wiki/Bravais_lattice). Primarily useful to check the passed parameters represent a valid Bravais lattice. If a lattice parameter is not specified, it will be assigned randomly (such that all lattice parameters remain compatible with the specified lattice system).

## Installation

`pip install bravais`

## Examples

Import the `BravaisLattice` class:

```python
from bravais import BravaisLattice
```

Quickly generate a monoclinic Bravais lattice without specifying any lattice parameters:

```python
mon_lat = BravaisLattice('monoclinic')
print(mont_lat)
```
```
P-centred monoclinic lattice (a=5.9417, b=4.7245, c=5.7335, alpha=90.00, beta=90.00, gamma=51.01)
```

Generate a body-centred tetragonal Bravais lattice with particular lattice parameters:

```python
tet_lat = BravaisLattice('tetragonal', 'I', a=3)
```

```
I-centred tetragonal lattice (a=3.0000, b=3.0000, c=3.5708, alpha=90.00, beta=90.00, gamma=90.00)
```

Note that the following single-digit codes are used to specify centring-types:

```
P -> primitive
B -> base
I -> body
F -> face
R -> rhombohedral
```
