# Tectonic Utils

A simple library with utility functions for "everyday use" in active tectonics research. There are geodesy and seismology utility functions. There are also a few read-write functions to access commonly used formats.

## Dependencies
* Python 3
* numpy
* scipy 
* netCDF4 (pip install netCDF4)


## Installation
On Unix (other systems not tested yet), install from PyPI using: 

```pip install Tectonic-Utils```

To test that it works, import a utility in Python. I typically use an import structure like: 
```
$ python
>>> from Tectonic_Utils.geodesy import haversine
```

## Documentation
I am beginning the process of documenting this library using sphinx. The preliminary html docs can be found here: https://kmaterna.github.io/Tectonic_Utils/index.html
