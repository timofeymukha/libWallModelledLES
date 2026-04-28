# libWallModelledLES

This is a OpenFOAM-based library for wall-modelled large-eddy simulation.
It mostly implements wall models as boundary conditions, with related
functionality like finding sampling cells.

## Minimal theoretical background

- Wall models fundamentally work by setting the shear stress value at the
boundary. In OpenFOAM this is done by setting an eddy viscosity with a value
such that the desired stress is realized.
- Wall models act in 3 stages for each wall boundary face for each timestep
  1. Sample the LES solution at some distance h away from the wall face, along
     the normal direction.
  2. Compute the wall stress using the equations constituting the wall model.
  3. Assign the corresponding value of eddy viscosity at the wall face.
- To predict the stress the models in the library either use
  - An algebraic equation (explicit or implicit).
  - An ODE equation (usually pre-integrated once).


## Repo structure

- `cellFinders`, algorithms for finding cells from which to sample the solution
  into the wall model.
- `eddyViscosities`, 1D RANS models for ODE-based wall models. NB: *not* LES
  SGS models.
- `lawsOfTheWall`, velocity profiles in inner scaling underpinning algebraic
  wall models.
- `rootFinding`, algebraic equation solvers used in algebraic wall models.
- `samplers`, code for sampling the LES solution into the wall models.
- `wallModels`, the implementation of the wall models as OpenFOAM boundary
   conditions.
- `tests`, unit and integration tests.

## Building

- The `wmake` build system is used.
- The `Allwmake` / `Allwclean` scripts build and clean up the build, respectively. Allways use these scripts and not `wmake` directly.
- Pass `-j` to `Allwmake` to speed up compilation.

## Testing
- `tests` contains unit tests structured in folders matching the code in the
  repo.
- Google Test is used for unit tests. It is found via $GTEST_DIR.
- `fixtures.*` contains useful helper routines.
- `testRunner.C` binds the tests together into a single executable.
- `tests/inegrationTests` contains integration testsa, also using Google Test,
  but leveraging a real OpenFOAM test case, which is found under
  `tests/testCases/channelFlow`.
- `tests/testCases/channelFlow` can also be used for ad hoc manual tests.

## Documentation

- The code itself used Doxygen, with custom OpenFOAM formatting.
- A conversion to readthedocs is performed afterwards, and the source
  documentation is bundled together with resturturedText files, which serve as
  a user guide.
- Any breaking change should be put into the CHANGELOG.md. The same for new
  functionality.

## Code style
- Use the OpenFOAM style guidelines and coding conventions.
- Always prefer using OpenFOAM-internal types, e.g. smart pointers, containers,
  etc.

## Misc rules
- Never stage or commit code to git. This will be done manually.
