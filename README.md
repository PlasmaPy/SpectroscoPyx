# SpectroscoPy
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](./LICENSE.md)
[![Build Status](https://travis-ci.org/PlasmaPy/SpectroscoPy.svg?branch=master)](https://travis-ci.org/PlasmaPy/SpectroscoPy)
[![Coverage Status](https://coveralls.io/repos/github/PlasmaPy/SpectroscoPy/badge.svg?branch=master)](https://coveralls.io/github/PlasmaPy/SpectroscoPy?branch=master)

An open source community developed Python 3.6+ package for spectroscopy in 
the early stages of development. SpectroscoPy intends to be for spectroscopy what
[Astropy](https://github.com/astropy/astropy) is for astronomy - a 
collection of commonly used programs to be shared between spectroscopists 
and researchers globally, running within and leveraging the open source 
scientific Python ecosystem. 

Goals include:
* Creating a generalized syntax for describing atomic energy levels and transitions
* Creating a single interface for multiple spectroscopic databases (e.g., NIST, CHIANTI), opacity databases, and lineshape databases
* Maintaining a system of unique identifiers to precisely share and cite which databases and versions were used for a particular analysis run
* Generating tools convenient for common tasks in simulation and analysis of spectroscopic data
* Creating a spectroscopic simulation suite which can use the aforementioned interfaces to readily compare spectra generated using different databases


SpectroscoPy requires Python 3.6 and is
[not compatible with Python 2](https://pythonclock.org/).

## License

SpectroscoPy is licensed under a 3-clause BSD license with added protections
against software patents - see the [``LICENSE.md``](LICENSE.md) file in
the top-level directory.

