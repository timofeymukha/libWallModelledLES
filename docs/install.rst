Installation
============

OpenFOAM version compatibility
------------------------------

The libWMLES library supports two forks of OpenFOAM, developed by `OpenCFD <https://openfoam.com/>`_, and the
`OpenFOAM Foundation <https://openfoam.org/>`_, repsectively.
Curently, all releases from OpenCFD up to v2212 are supported, as well as releases 3.0.x to 7 from the Foundation.
It should be stressed that *supported* means that it is checked that the library compiles.
This should in principle entail that it also works properly, but due to time constraints and the ambitious number of
versions suported, it is impossible to actually run test cases for all of them.
That being said, libWMLES has been extensively used with version 3.0.x, 4.1, v1806 and v1912 of OpenFOAM.
Reports regarding difficulties with running a particular version are welcome: please open an issue on Bitbucket.

The best level of testing is done for latest available release from OpenCFD, for which the test harness is run.
This consists of unit and integration tests that aim to cover all the functionality of the library.
A single system test is also run, which is a simulation of channel flow on a coarse grid, using the Spalding-law wall
model and otherwise default parameters in the :code:`nut` dictionary.

The aim is to support all new OpenCFD, meaning that the amount of supported versions grows with
two per year.
This is likely to become unsustainable, leading to deprecation of support for the oldest versions.

The support for new Foundation versions would be nice, but currently dropped due to lack of time and familiarity with the codebase.
Contributions from users of the Foundation version to provide support are most welcome. 

Compilation
-----------

Clone the repository with git or download it as an archive by navigating Downloads in the menu on Bitbucket.
A prerequisite for installing is having Python installed, but no packages are needed and any modern Python version
should do.
To compile, run ``Allwmake``.
Consider using the ``-j`` flag with a number of processors, to speed up the process.
If you get compilation errors, please make sure your OpenFOAM environment is properly set up before opening a Bitbucket
issue.
In particular, take notice of the first couple of lines in the output for ``Allwmake``, which are of the following
form::

   Current OpenFOAM version is v1806.
   This is a clean install
   OpenFOAM-version: Major 1806 Minor 0 Patch 0 (-1 == x / 0) Fork: com

As you can see, they state which version of OpenFOAM the install script has picked up.
If it doesn't correspond to your expectations, there is almost surely something fishy with you setup.

When recompiling for a different version of OpenFOAM than what was previously used, you should frist run ``Àllwclean``.
If you make several compilation attempts and things don't work, it can be a good idea to delete ``Make/linux*`` and
``lnInclude`` to make sure you start from a clean slate.

Running tests
-------------

To run the tests, Google Test and Google Mock should first be installed, see https://github.com/google/googletest.
The environmental variable :code:`GTEST_DIR` should point to the root directory containing both Google Test and Google
Mock (i.e. the directory to which you've clone the `googletest` repo).
Don't forget to copy the compiled libraries to your :code:`FOAM_USER_LIBBIN`.

The unit tests are located in the :code:`tests` directory, are compiled with :code:`wmake` producing the file `testRunner`,
which could be executed to run the tests.
The integration tests are located in :code:`test/integrationTests`, are also compiled with `wmake`, and the produced
executable is called `testIntegration`.
