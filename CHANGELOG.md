# CHANGELOG

## v0.7.0

### For users
- The main repository is now on Github. Use Github issues instead of Bitbucket to get help!

- The library works with up to release v2212 from OpenCFD.

- Ambition to support Foundation version is dropped. Pull requests are much welcome though.

- The field for the sampling distance height is now called `hSampler` to avoid name collisions with enthalpy.
  However, if `hSampler` is not found, `h` will be picked up.

- It is now possible to use the library with compressible solvers.

- Intracell interpolation should work properly now.

- New parameter "silent", which suppresses time consumption output.

- Messages like "Sampling field ... for patch ..." are now only printed in debug mode.
  Use `DebugSwitches` in the `controlDict` to enable debug.
  
- A new law of the wall is added, the `RoughLogLaw`, suitable for rough surfaces.

- As a beta release, multi-cell sampling is implemented.
  No simulations have been run using it as of yet.
  
- The online documentation is improved.

### For developers
- The SampledField classes now have `TypeName`.

- 'excludeWallAdjacent' is now in the base Sampler class.

- Slightly different methods to get the viscosity internal field and patch field in the base wall model class (to accommodate compressible solvers).


## v0.6.1
### For users

- The library now works with up to release v2206 from OpenCFD.

## v0.6.0

### For users

- This is the last version to support versions 3.x.

- Documentation portal now available at https://libwmles.readthedocs.io

- It is now possible to choose between two different samplers that differ in
how sampling cells are located. This is done by sampling the `sampler` keyword
to `Tree` and `Crawling`. The default is `Tree`, which corresponds to the
behaviour in previous versions. For mor information, see the documentation.

- It is now possible to tell the software to treat the values of the h field
as the consecutive index of the cell to be sampled from instead of the distance
to be sampled from. This is controlled by the hIsIndex keyword. Default is 0,
corresponding to behaviour in previous versions.

- It is now possible to interpolate the sampled values to a particular point
within the sampling cell. This means that the sampling distance does not
have to strictly correspond to a distance to the cell centre.
To that end, the `interpolationType` is used, with the default value of `cell`, 
which corresponds to using the cell-centre value.
This functionality uses the `interpolation` class available in OpenFOAM,
with corresponding options for the interpolation procedure.
See documentation for more details.

- The NoModel SGS model is removed (to make compilation with the Foundation
  version of OF easier to manage)

### For developers

- Searching for sampling cells is now removed from the Sampler classes and
put into CellFinder classes.

- Useful helper functions now present in helpers.C/H files.
  
- SampledField classes no longer RTS.

- There is now a MultiCellSampler class for building models based on data
sampled from multiple consecutive cells. SampledFields also have support for
this.
  


## v0.5.2

Hot fix for the default value of  `averagingTime` not being 0.

## v0.5.1

Hot fix for `copyToPatchInternalField` not getting written to the `nut` file by `decomposePar`.

## v0.5.0

This release contains multiple improvements to the code structure, adds a unit test suite and
better handling of sampling cell search for large cases.
Stability of wall modelling is also improved.

### For users

- **BREAKING** The sampled fields have changed type from vectorFields to scalarListLists.
  To run the same case with the new version delete the wallModelSampling directory for the
  time-step you want to start running from.
  This will not affect the simulation unless you've used time-averaging of sampled data.

- Optionally, the values of `nut` predicted by the wall model can be copied to the wall-adjacent cells.
  This improves the stability of wall modelling. The copying is controlled by the `copyToPatchInternalField`
  switch in the input dictionary.
  Default behaviour is not to copy, meaning that old cases will behave as before without modification.

- The fill value in the samplingCells field is now -1 instead of 0, to allow for a wall-patch
  with index 0.

- Default tolerance for all `RootFinder`s is now 0.001 and the maximum number of iterations is 30.

- When a field is sampled, an associated entry to the log is made.

- The wall model outputs the time it consumed in percent of the total CPU runtime.

- An issue with compiling the library for version 3.0+ from OpenCFD has been fixed.

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

- The Sampler now explicitly constructs the sub-registry corresponding to data sampled from the patch. This is instead of using the `subRegistry()` method, which does not support creating a new sub-registry in older OF versions. 

- For the same purpose, accessing the subregistry via `subRegistry()` is used without a second argument to the method.

- The turbulence.H header is now only used for the `NoModel` SGS model code, and not for the wall models. The viscosity field is grabbed directly from the registry and not via the turbulence model.
