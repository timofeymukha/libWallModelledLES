# README #

libWallModelledLES is a library based on OpenFOAM� technology, extending the capabilities of OpenFOAM in the area of wall-modelled LES (WMLES).
In particular, so-called wall-stress models are considered. These aim at correctly predicting the wall shear stress at the wall without the need for the LES mesh to resolve the inner part of the turbulent boundary layers.
Note that, unlike some other approaches (e.g. hybrid LES/RANS), the LES domain here extends all the way to the wall and **only the inner layer is modelled**, whereas the outer layer of TBLs is fully-resolved.
Chapter 5 of the following thesis may be of interest for getting further acquainted with the methodology of WMLES, see also the publication list in the end of the README.

http://www.diva-portal.org/smash/record.jsf?pid=diva2%3A1236761

To simplify application to general geometries the models in the library predict the magnitude of the wall shear stress instead of its individual components.
Similarly to OpenFOAM's native wall functions, the value of the shear stress is enforced by setting a non-zero value to the turbulent viscosity at the wall.
Therefore, the models are chosen and configured individually for each wall boundary in the 0/nut file, see below.

The library provides a set of new models, both based on non-linear algebraic equations (laws-of-the-wall) and ordinary differential equations.
Fine grain control over the models' behaviour is given to the user.
The library also provides developers a convenient framework to quickly add new models.

If you use the library, please cite the following publication, which fully describes the implemented functionality.

https://doi.org/10.1016/j.cpc.2019.01.016

**This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM� and OpenCFD� trademarks.**

## News ##

- **2019-08-01** Version 0.5 released.
- **2019-02-23** Version 0.4.1 released, containing a small bugfix.
- **2018-11-17** Version 0.4 released, see CHANGELOG.md for list of changes.

## Compatibility ##

The library compiles for versions 3.0.x to 7.x from the OpenFOAM Foundation, and versions 3.0+ to 1906 from OpenCFD. 
A special branch on the repository for version 2.3.1 is available but is now very outdated and will not be further supported.

## Installing ##

Clone the repository to the directory of your choice and run the Allwmake file inside.
This should be it!
If you want to get a specific version of the library, go to Downloads in the menu on the left, then to Tags, and download the associated archive.

If you want to build the source code documentation with doxygen, go into the
docs folder and run "doxygen config".
This will create an html folder that can be read using a browser.

## Key features ##

- Provides a number of wall models, based on both non-linear algebraic and ordinary differential equations, see the class-headers in the wallModels folder.
- Makes it possible to specify the distance to the wall model's sampling point, h, on a per-face basis.
- Allows the user to control all the other parameters of wall modelling, e.g. model constants, iterative solver settings etc.
- Serves a as a convenient framework for implementing new models without a lot of code duplication.

## Validation cases ##
In the *tests* folder there is a toy channel flow case that you can try running to make sure that things compiled well.
No more simulation cases are shipped with the library.
However, OpenFOAM cases for turbulent channel flow and flow over a backward-facing step can be found using the following DOI: 10.6084/m9.figshare.6790013
These thus serve both as tutorials for case set-up and as validation cases.
Further results obtained using the library can be found in the publications listed below.
The cases considered in these publication can thus also be used for validation.

## Running inlcuded  unit and integration tests
To run the tests, Google Test and Google Mock should first be installed, see https://github.com/google/googletest.
The environmental variable `GTEST_DIR` should point to the root directory containing both Google Test and Google Mock (i.e. the directory to which you've clone the `googletest` repo).

The unit tests are located in the `tests` directory, are compiled with `wmake` producing the file `testRunner`, which could be executed to run the tests.
The integration tests are located in `test/integrationTests`, are also compiled with `wmake`, and the produced executable is called `testIntegration`.

## Case set-up ##
Assume that you've already set up a case for the classical wall-resolved LES. To convert it to WMLES you need to do the following:

- Add libWallModelledLES.so to the loaded libraries in the controlDict.
- Go into the *nut* file, and set up wall models as the boundary conditions at the walls.
  This is similar to what you would do with OpenFOAM's built-in model based on Spalding's law, *nutUSpaldingWallFunction*.
  The set-up and parameters for each model in the library can be found in the header of the associated .H file, see list of files below.
- In your *0* directory, you should add a new volScalarField, *h*.
  This will hold the distance from the faces to the cell-centre which will be used for sampling data into the wall model.
  The value of the internalField is irrelevant, you can put it to some constant scalar, e.g. 0.
  At boundaries where wall modelling is not applied, *zeroGradient* can be used.
  At the walls where wall modelling is applied, the value of h should be provided.
  The value 0 is reserved to indicate sampling from the wall-adjacent cell.
  You can either provide a uniform value for the whole patch or a list of scalars, with separate values for each face.
  To do the latter conveniently based on some criteria, using funkySetFields is recommended, it is a utility, which is part of swak4foam.

## Documentation ##

Each class is documented in the corresponding .H header file.
This includes usage instructions, and, where applicable, formulas and references to literature.
A tiny toy case can be found under tests/testCases/channel\_flow where the 0/nut file provides an example of setting up the boundary conditions.
The article linked to in the beginning of the README also serves as documentation.

## Best practice guidelines ##
This is intended to be a summary of tips based on the experience of the developers and users of the library.
The intention is to give a good starting point for new users.
Naturally, results may vary heavily depending on the case in question.

- In the boundary layer, define your grid density as the number of cells per delta^3, where delta is the thickness of the boundary layer.
  A good number is 27000 cells, but you can get good results with less.
  This will need an a-priori knowledge of the distribution of delta across the wall.
  A RANS precursor can do the job.
- Use an isotropic grid in the boundary layer.
  In particular, no need to refine it towards the wall, which may at first seem weird for practitioners of hybrid RANS/LES methods.
  For some inspiration on unstructured meshing strategies see (Mukha, Johansson & Liefvendahl, in ECFD 7, Glasgow, UK, 2018).
- In regions where the TBL is attached, set *h* to be the distance to the second consecutive off-the wall cell centre. In other regions, use *h=0*,
  i.e. sample from the wall-adjacent cell.
- Use a mildly diffusive numerical scheme, e.g. LUST. Tips regarding what other schemes worked well are welcome :).
- The WALE model is a good first choice for SGS modelling.
- If your simulation crashes because of the wall model (you can usually see that in the log), make sure you have a good initial condition.
- If your simulation crashed anyway, use *h =0*, this is pretty much guaranteed to be stable.
- Large values of *h* are known to sometimes lead to a crash, in particular, if the grid below *h* is refined.
- If you use *h=0*, use an algebraic wall model in integral formulation, i.e. the *LOTWWallModel* with e.g. the *IntegratedReichardt* law.
- Use a low tolerance and a decent number of iteration for the Newton solver, this will remove occasional spikes in *nut* that may occur when the solver is not converged but have no impact on the performance in general.
  A tolerance of 0.0001 and 20-30 iterations is usually a good choice.
- Similarly, for ODE models, use a relatively dense 1D grid, e.g. 50 points.
  There is no large impact on performance either.

## Source files' contents

The contents of the files in each folder is briefly described below.
Most classes are implemented in a pair of files with the same name ending with .C and .H, as is customary in C++ and OpenFOAM.
Each such pair is treated as one item in the list below, without providing the file extension.

- eddyViscosities
    * Duprat/DupratEddyViscosity Class for eddy viscosity based on (Duprat et al, Physics of Fluids, 2011).
    * EddyViscosity/EddyViscosity Base abstract class for eddy viscosity models used by ODE wall models.
    * JohnsonKing/JohnsonKingEddyViscosity Class for eddy viscosity based on the mixing length model with vanDriest damping (van Driest, Journal of the Aeronautical Sciences, 1956).

- lawsOfTheWall
    * IntegratedReichardtLawOfTheWall/IntegratedReichardtLawOfTheWall Class for the integrated formulation of Reichardt's law, (Reichardt, Zeit. f�r Ang. Math. und Mech., 1951).
    * IntegratedWernerWengleLawOfTheWall/IntegratedWernerWengleLawOfTheWall Class for the integrated formulation of the law of Werner and Wengle, (Werner & Wengle, Turb. Shear Flows 8, 1991).
    * LawOfTheWall/LawOfTheWall Base abstract class for laws of the wall.
    * ReichardLawOfTheWall/ReichardLawOfTheWall Class for Reichardt's law of the wall, (Reichardt, Zeit. f�r Ang. Math. und Mech., 1951).
    * SpaldingLawOfTheWall/SpaldingLawOfTheWall Class for Spalding's law of the wall, (Spalding, J. of Applied Mechanics, 1961).
    * WernerWengleLawOfTheWall/WernerWengleLawOfTheWall Class for Werner and Wengle's law of the wall, (Werner & Wengle, Turb. Shear Flows 8, 1991).
- Make
    * files File used by wmake to determine what source files to compile.
    * options File used by wmkae to determine what libraries and headers to include at compilation.
- rootFinding
    * BisectionRootFinder/BisectionRootFinder Class for a root finder implementing the bisection method.
    * NewtonRootFinder/NewtonRootFinder Class for a root finder implementing Newton's method.
    * RootFinder/RootFinder Base abstract class for root finders.
- samplers
    * SampledField/SampledField Base abstract class for a field to be sampled by the wall models.
    * SampledField/SampledPGradField Class for sampling the pressure gradient
    * SampledField/SampledVelocityField Class for sampling the velocity
    * SampledField/SampledWallGradUField Class for sampling the wall-normal gradient of velocity.
    * Sampler/Sampler
- sgsModels
    * makeSGSModel.C Helper file to create a new turbulence model
    * NoModel Class for an SGS model with zero SGS viscosity in the internal field.
- tests
- versionRules
    * codeRules.H Defines macros based on the version of OpenFOAM which is used.
    * libraryRules.H Defines locations of libraries included in Make/options depending on the OpenFOAM version used.
    * makeFoamVersionHeader.py A Python script that determines the version of OpenFOAM which is used and writes-out associated data to foamVersion4wmles.H
- wallModels
    * EquilibriumODEWallModelFvPatchScalarField Class for the ODE-based wall model with a zero source term.
    * KnownWallShearStressWallModelFvPatchScalarField Class for wall model that reads a-priori known wall shear stress from disk.
    * LOTWWallModelFvPatchScalarField Class for algebraic (law of the wall based) wall models.
    * ODEWallModelFvPatchScalarField Base class for ODE-based wall models.
    * PGradODEWallModelFvPatchScalarField Class for ODE-based wall model with a source term equal to the pressure gradient.
    * wallModelFvPatchScalarField Base abstract class for wall models.

## Published works using the library

- Mukha, T., Rezaeiravesh, S., & Liefvendahl, M. (2017). An OpenFOAM library for wall-modelled Large-Eddy Simulation. In proceedings of the 12th OpenFOAM Workshop, Exeter, UK.
- Mukha, T., Johansson, M., & Liefvendahl, M. (2018). Effect of wall-stress model and mesh-cell topology on the predictive accuracy of LES of turbulent boundary layer flows.
  In 7th European Conference on Computational Fluid Dynamics, Glasgow, UK.
- Mukha, T., Rezaeiravesh, S., & Liefvendahl, M. (2018). Wall-modelled large-eddy simulation of the flow over a backward-facing step. In proceedings of 13th OpenFOAM Workshop, Shanghai, China. Shanghai, China.
- Liefvendahl, M., & Johansson, M. (2018). Wall-Modeled LES for Ship Hydrodynamics in Model Scale. In proceedings of the 32nd Symposium on Naval Hydrodynamics, Hamburg, Germany.
- Bezinge, G. (2018) Wall-unresolved large eddy simulation of turbulent flow at high Reynolds number: Performance and computational cost investigation. Master thesis.
  Department of Mathematics University of Wyoming (UW) Laramie Institute of Fluid Dynamics Swiss Federal Institute of Technology (ETH) Z�rich.
- Rezaeiravesh, S., Mukha, T., & Liefvendahl, M. (2018). Systematic study of accuracy of wall-modeled large eddy simulation using uncertainty quantification techniques. Available: https://arxiv.org/abs/1810.05213
- Mukha, T., Rezeeiravesh, S., & Liefvendahl, M. (2019). A library for wall-modelled large-eddy simulation based on OpenFOAM technology. Computer Physics Communications. DOI: 10.1016/j.cpc.2019.01.016. Preprint: https://arxiv.org/abs/1807.11786
