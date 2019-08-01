# CHANGELOG

## v0.5.0

This release contains multiple improvements to the code structure, adds a unit test suite and
better handling of sampling cell search for large cases.

### For users

- **BREAKING** The sampled fields have changed type from vectorFields to scalarListLists.
  To run the same case with the new version delete the wallModelSampling directory for the
  time-step you want to start running from.
  This will not affect the simulation unless you've used time-averaging of sampled data.

- The fill value in the samplingCells field is now -1 instead of 0, to allow for a wall-patch
  with index 0.

- Default tolerance for all `RootFinder`s is now 0.001 and the maximum number of iterations is 30.

- When a field is sampled, an associated entry to the log is made.

- The wall model outputs the time it consumed in percent of the total CPU runtime.

- An issue with compiling the library for version 3.0+ for OpenCFD has been fixed.

### For developers

- Library now includes a Google Test unit test suite covering a significant part of the functionality.
  However, currently the tests cannot be run on all the OpenFOAM versions supported by the library.

- Library now includes an integration test suite, that checks running in serial, running in parallel,
  decomposing, and reconstructing the case for each wall model.

- Classes `scalarListListList` and `scalarListListIOList` are introduced to hold sampled data from multiple cells.

- Multiple changes in the Sampler classes
    - Raw pointers are now longer used to store the list of sampled fields.
    - Samplers are now RTS, although users are not offered to select a sampler yet.
    - A heirarchy is established, currently with `Sampler` as the base class and `SingleCellSampler` and `MultiCellSampler` as children, corresponding to sampling from a single or multiple cells, repsectively.
    - The `Sampler` class now compute a distance field in function `distanceField` unless one can be read from disk.
    - The `Sampler` class now searches for candidate sampling cells in `findSearchCellLabels()` function.


- Multiple changes in `SampledField` classes
    - They are now RTS, although this is currently not used.
    - Support sampling from multiple cells into a `scalarListListList` and registering associated
      fields.
    - Can be cloned.
    - `db()` now returns the local wallModel registry, and `mesh()` returns the mesh.

- The base wallModel class no longer holds a `cellIndexList`, `wallGradU` or the `Sampler` as members.
  A sampler of appropriate type is owned only by models that need it, like the LOTW or ODE model,
  and the `wallGradU` and the `cellIndexList` can be accessed through the sampler. 

- Multiple changes in LOTW and EddyViscosity classes.
    - Can be cloned.
    - No longer own a `Sampler` pointer, instead functions like `value()` accept a `SingleCellSampler`
      as a const reference.
    - Have `addFieldsToSampler()` function to add `SampledFields` to a sampler.
    - Functions like `value()` accepting necessary data as parameters are also present.
      For example, `SpaldingLawOfTheWall::value(sampler, index, uTau, nu)` grabs u and y from the sampler and calls `SpaldingLawOfTheWall::value(u, y, uTau, nu)`.
    - In LOTWs, the index parameter to `value()` and `derivative()` is now a label, as it should be.
    - Have direct constructors from relevant model constants, e.g. `SpaldingLawOfTheWall(kappa, B)`.
- Constructors for wallModels use `clone()` methods for members that are `autoPtr<T>`, where
  appropriate.
  This fixes a bug that led to some copy constructors leaving the original in an unusable state.

## v0.4.1
This a hotfix release, squashing a bug in one of the constructors of the LOTW model.

## v0.4.0

### For users
- The code should now be installed by running the Allwmake file.

- Python is now necessary to compile the library. Only standard library packages are used.

- The library now compiles under several versions of OpenFOAM from OpenCFD and the OpenFOAM Foundation. This includes versions 3.x, 4.x, 5.x, and 6.x., 1606+, 1612+ and 1806. If makeSGSModel.C is excluded from the list of files (corresponds to the NoModel SGS model), compilation works for foam-extend 3.2 and 4.0 as well, although the code does not actually seem to work.


### For developers
- The multiple versions are accounted for using pre-processor macros, with code borrowed from swak4foam, see the versionRules folder.

- The Sampler now explicitely constructs the sub-registry corresponding to data sampled from the patch. This is instead of using the `subRegistry()` method, which does not support creating a new sub-registry in older OF versions. 

- For the same purpose, accessing the subregistry via `subRegistry()` is used without a second argument to the method.

- The turbulence.H header is now only used for the `NoModel` SGS model code, and not for the wall models. The viscosity field is grabbed directly from the registry and not via the turbulence model.
