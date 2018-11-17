# CHANGELOG

## v0.4.0 ##

- The code should now be installed by running the Allwmake file.

- Python is now necessary to compile the library. Only standard library packages are used.

- The library now compiles under several versions of OpenFOAM from OpenCFD and the OpenFOAM Foundation. This includes versions 3.x, 4.x, 5.x, and 6.x., 1606+, 1612+ and 1806. If makeSGSModel.C is excluded from the list of files (correspond sto the NoModel SGS model), compilation works for foam-extend 3.2 and 4.0 as well, although the code does not actually seem to work.

- The multiple versions are accounted for using pre-processor macros, with code borrowed from swak4foam, see the versionRules folder.

- The Sampler now explicitely constructs the sub-registry corresponding to data sampled from the patch. This is instead of using the subRegistry() method, which does not support creating a new sub-registry in older OF versions. 

- For the same purpose, accessing the subregistry via subRegistry() is used without a second argument to the method.

- The turbulence.H header is now only used for the NoModel SGS model code, and not for the wall models. The viscosity field is grabbed directly from the registry and not via the turbulence model.
