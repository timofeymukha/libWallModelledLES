# README #

libWallModelledLES is a library based on OpenFOAM® technology, extending the
capabilties of OpenFOAM in the area of wall-modelled LES (WMLES). In
particular, so called wall-stress models are considered. These aim at
correctly predicting the wall shear stress at the wall without the need
for the LES mesh to resolve the inner part of the turbulent boundary layers.
Note, that, unlike some other  approaches (e.g. hybrid LES/RANS), the LES
domain here extends all the way to the wall.

To simplify application to general geometries the models in the library
predict the magnitude of the wall shear stress instead of its individual
components.
Similarly to OpenFOAM's native wall functions, the value of the shear stress
is enforced by setting a non-zero value to the turbulent viscosity at the
wall.
Therefore, the models are chosen and configured individually for each wall
boundary in the 0/nut file.

The library provides a set of new models, both based on non-linear algebraic
equations (laws-of-the-wall) and ordinary differential equations.
Fine grain control over the models' behaviour is given to the user.
The library also provides developers a convenient framework to quickly add new 
models.

## Compatibility ##

The library has been developed using OpenFOAM version 3.0.1. It has also been
compiled using OpenFOAM 2.3.1. Versions 4.x are currently not supported due
to breaking changes in the API (dereferencing, in particular).

## Installing ##

Clone the repository to the directory of your choice and run wmake inside.
This should be it!

If you want to build the source code documentation with doxygen, got into the
docs folder and run "doxygen config". This will create an html folder that
can be reat using a browser.

## Documentation ##

Each class is documented in the corresponding .H header file. This includes
usage instructions, and, where applicable, formulas and references to
literature. A tiny toy case can be found under tests/testCases/channel\_flow
where the 0/nut file provides an example of setting up the boundary
conditions.


## Disclaimer ##

This offering is not approved or endorsed by OpenCFD Limited, producer and
distributor of the OpenFOAM software via www.openfoam.com, and owner of the
OPENFOAM®  and OpenCFD®  trade marks.
