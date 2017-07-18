# README #

libWallModelledLES is a library based on OpenFOAM® technology, extending the
capabilties of OpenFOAM in the area of wall-modelled LES.

The library provides a set of new models, both based on non-linear algebraic
equations (laws-of-the-wall) and ordinary differential equations.
Fine grain control over the models' behaviour is given to the user.
The library also provides develpoers a convenient framework to quickly add new 
models.

## Compatibility ##

The library has been developed using OpenFOAM version 3.0.1. It has also been
compiled using OpenFOAM 2.3.1. Versions 4.x are currently not supported due
to breaking changes in the API (dereferencing, in particular).

## Installing ##

Clone the repository to the directory of your choice and run wmake inside.
This should be it!

## Disclaimer ##

This offering is not approved or endorsed by OpenCFD Limited, producer and
distributor of the OpenFOAM software via www.openfoam.com, and owner of the
OPENFOAM®  and OpenCFD®  trade marks.
