# README #

libWallModelledLES is a library based on OpenFOAM® technology, extending the capabilities of OpenFOAM in the area of wall-modelled LES (WMLES).
In particular, so-called wall-stress models are considered. These aim at correctly predicting the wall shear stress at the wall without the need for the LES mesh to resolve the inner part of the turbulent boundary layers.
Note, that, unlike some other  approaches (e.g. hybrid LES/RANS), the LES domain here extends all the way to the wall.

To simplify application to general geometries the models in the library predict the magnitude of the wall shear stress instead of its individual
components.
Similarly to OpenFOAM's native wall functions, the value of the shear stress is enforced by setting a non-zero value to the turbulent viscosity at the wall.
Therefore, the models are chosen and configured individually for each wall boundary in the 0/nut file.

The library provides a set of new models, both based on non-linear algebraic equations (laws-of-the-wall) and ordinary differential equations.
Fine grain control over the models' behaviour is given to the user.
The library also provides developers a convenient framework to quickly add new models.

**This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD® trademarks.**

## Compatibility ##

The library has been developed using OpenFOAM version 3.0.1.
A branch for version 2.3.1 is available but does not contain all the latest features.
Versions 4.x and up are currently not supported due to breaking changes in the API (dereferencing, in particular).
Setting up a more complicated compilation procedure to support more version is on the todo list, currently, no promises can be made regarding when it will be done.
Contributions are highly welcome :).

## Installing ##

Clone the repository to the directory of your choice and run wmake inside.
This should be it!

If you want to build the source code documentation with doxygen, go into the
docs folder and run "doxygen config".
This will create an html folder that can be read using a browser.

In the *tests* folder there is a toy channel flow case that you can try running to make sure that things compiled well.

## Case set-up ##
Assume that you've already set up a case for the classical wall-resolved LES. To convert it to WMLES you need to do the following:

- Add libWallModelledLES.so to the loaded libraries in the controlDict.
- Go into the *nut* file, and set up wall models as the boundary conditions at the walls.
  This is similar to what you would do with OpenFOAM's built-in model based on Spalding's law, *nutUSpaldingWallFunction*.
  The set-up and parameters for each model in the library can be found in the header of the associated .H file, see list of files below.
- In your *0* directory, you should add a new volScalarField, *h*.
  This will hold the distance from the faces to the cell-centre which will be  used for sampling data into the wall model.
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

## Best practice guidelines ##
This is intended to be a summary of tips based on the experience of the developers and users of the library.
The intention is to give a good starting point for new users.
Naturally, results may vary heavily depending on the case in questions.

- In the boundary layer, define your grid density as the number of cells per delta^3, where delta is the thickness of the boundary layer.
  A good number is 27000 cells, but you can get good results with less.
  This will need an apriori knowledge of the distribution of delta across the wall.
  A RANS precursor can do the job.
- Use an isotropic grid in the boundary layer.
  In particular, no need to refine it towards the wall, which may seem weird after wall-resolved LES.
  For some inspiration on unstructured meshing strategies see (Mukha, Johansson & Liefvendahl, in ECFD 7, Glasgow, UK, 2018).
- In regions where the TBL is attached, set *h* to be the distance to the second consecutive off-the wall cell centre. In other regions, use *h=0*,
  i.e. sample from the wall-adjacent cell.
- Use a mildly diffusive numerical scheme, e.g. LUST. Tips regading what other schemes worked well are welcome :).
- The WALE model is a good first choice for SGS modelling.
- If your simulation crashes because of the wall model (you can usually see that in the log), make sure you have a good initial condition.
- If your simulation crashed anyway, use *h =0*, this is the most stable regime.
- Large values of *h* are known to sometimes lead to a crash, in particular, if the grid below *h* is refined.
- If you use *h=0*, use an algebraic wall model in integral formulation, i.e. the *LOTWWallModel* with e.g. the *IntegratedReichardt* law.
- Use a low tolerance and a decent number of iteration for the Newton solver, this will remove spikes in *nut* that may occur when the solver is not converged.
  A tolerance of 0.0001 and 20-30 iterations is usually a good choice.
  All te iterations will be taken only if the the tolerance is not achived earlier.
- Similarly, for ODE models, use a relatively dense 1D grid, e.g. 50 points.
  There is no large impact on performance.

## Source files' contents

The contents of the files in each folder is briefly described below.
Most classes are implemented in a pair of files with the same name ending with .C and .H, as is customary in C++.
Each such pair is treated as one item in the list below, without providing the file extension.

- eddyViscosities
    * Duprat/DupratEddyViscosity Class for eddy viscosity based on (Duprat et al, Physics of Fluids, 2011).
    * EddyViscosity/EddyViscosity Base abstract class for eddy viscosity models used by ODE wall models.
    * JohnsonKing/JohnsonKingEddyViscosity Class for eddy viscosity based on the mixing length model with vanDriest damping (van Driest, Journal of the Aeronautical Sciences, 1956).

- lawsOfTheWall
    * IntegratedReichardtLawOfTheWall/IntegratedReichardtLawOfTheWall Class for the integrated formulation of Reichardt's law, (Reichardt, Zeit. für Ang. Math. und Mech., 1951).
    * IntegratedWernerWengleLawOfTheWall/IntegratedWernerWengleLawOfTheWall Class for the integrated formulation of the law of Werner and Wengle, (Werner & Wengle, Turb. Shear Flows 8, 1991).
    * LawOfTheWall/LawOfTheWall Base abstract class for laws of the wall.
    * ReichardLawOfTheWall/ReichardLawOfTheWall Class for Reichardt's law of the wall, (Reichardt, Zeit. für Ang. Math. und Mech., 1951).
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