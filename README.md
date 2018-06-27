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

The library has been developed using OpenFOAM version 3.0.1.
A branch for version 2.3.1 is available, but does not contain all the latest featrues.
Versions 4.x and up are currently not supported due to breaking changes in the API (dereferencing, in particular).
Setting up a more complicated compilation procedure to support more version is on the todo list, currently no promises can be made regarding when it will be done.
Contributions are highly welcome :).

## Installing ##

Clone the repository to the directory of your choice and run wmake inside.
This should be it!

If you want to build the source code documentation with doxygen, go into the
docs folder and run "doxygen config".
This will create an html folder that can be read using a browser.

In the *tests* folder there is a toy channel flow case that you can try running to make sure that things compiled well.

## Case se-up ##
Assume that you've already set up a case for the classical wall-resolved LES. To convert it to WMLES you need to do the following:

- Add libWallModelledLES.so to the loaded libraries in the controlDict.
- Go into the *nut* file, and set up wall models as the boundary conditions at the walls.
  This is similar to what you would do with OpenFOAM's built in model based on Spalding's law, *nutUSpaldingWallFunction*.
  The set-up and parameters for each model in the library can be found in the header of the associated .H file, see list of files below.
- In your *0* directory, you should add a new volScalarField, *h*.
  This will hold the distance from the faces to the cell-center which will be  used for sampling data into the wall model.
  The value of the internalField is irrelevant, you can put it to some constant scalar, e.g. 0.
  At boundaries where wall modelling is not applied, *zeroGradient* can be used.
  At the walls where wall modelling is applied, the value of h should be provided.
  The value 0 is reserved to indicate sampling from the wall-adjacent cell.
  You can either provide a uniform value for the whole patch or a list of scalars, with separate values for each face.
  To do the latter conveniently based on some criteria, using funkySetFields is recommended, it is a utility, which is part of swak4foam.

## Documentation ##

Each class is documented in the corresponding .H header file. This includes
usage instructions, and, where applicable, formulas and references to
literature. A tiny toy case can be found under tests/testCases/channel\_flow
where the 0/nut file provides an example of setting up the boundary
conditions.

## Source files' contents

The contents of the files in each folder is briefly described below.
Most classes are implemented in a pair of files with the same name ending with .C and .H, as is customary in C++.
Each such pair is treated as one item in the list below, without providing the file extention.

- eddyViscosities
    * Duprat/DupratEddyViscosity Class for eddy viscosity based on (Duprat et al, Physics of Fluids, 2011)
    * EddyViscosity/EddyViscosity Base abstract calss for eddy viscosity models used by ODE wall models.
    * JohnsonKing/JohnsonKingEddyViscosity Class for eddy visocosity based on the mixing length model with vanDriest damping (van Driest, Journal of the Aeronautical Sciences, 1956)

- lawsOfTheWall
    * IntegratedReichardtLawOfTheWall/IntegratedReichardtLawOfTheWall Class for the integrated formulation of the algebraic model based on Reichardts' law. (Reichardt, Zeitschrift für Angewandte Mathematik und Mechanik, 1951)
    * IntegratedWernerWengleLawOfTheWall/IntegratedWernerWengleLawOfTheWall
    * LawOfTheWall/LawOfTheWall
    * ReichardLawOfTheWall/ReichardLawOfTheWall
	* SpaldingLawOfTheWall/SpaldingLawOfTheWall
	* WernerWengleLawOfTheWall/WernerWengleLawOfTheWall
- Make
    * files
	* options
- rootFinding
    * BisectionRootFinder/BisectionRootFinder
	* NewtonRootFinder/NewtonRootFinder
	* RootFinder/RootFinder
- samplers
	* SampledField/SampledField
	* SampledField/SampledPGradField
	* SampledField/SampledVelocityField
	* SampledField/SampledWallGradUField
	* Sampler/Sampler
- sgsModels
	* makeSGSModel.C
	* NoModel
- tests
- wallModels
    * EquilibriumODEWallModelFvPatchScalarField
    * KnownWallShearStressWallModelFvPatchScalarField
    * LOTWWallModelFvPatchScalarField
    * ODEWallModelFvPatchScalarField
    * PGradODEWallModelFvPatchScalarField
    * wallModelFvPatchScalarField



## Disclaimer ##

This offering is not approved or endorsed by OpenCFD Limited, producer and
distributor of the OpenFOAM software via www.openfoam.com, and owner of the
OPENFOAM®  and OpenCFD® trade marks.
