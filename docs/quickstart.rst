Quick start
===========

Case setup
----------

Assume that you've already set up a case for the classical wall-resolved LES. To convert it to WMLES you need to do the
following:

* Add :code:`libWallModelledLES.so` to the loaded libraries in the :code:`controlDict`.
* Go into the :code:`nut` file, and set up wall models as the boundary conditions at the walls.
  A minimalistic setup for a Spaling law-based algebraic wall model is given below.

.. code-block::
  :linenos:

  type           LOTWWallModel;
  RootFinder
  {
    type         Newton;
  }
  Law
  {
    type         Spalding;
  }

* In your :code:`0` directory, you should add a new volScalarField, :code:`hSampler`, see the :ref:`sampling` section for details.
  For a quick start, set the value of :code:`hSampler` to :code:`uniform 0` at the wall, and use :code:`zeroGradient` at all
  non-wall patch boundaries.
  This will lead to sampling from the wall adjacent-cell, which is very robust, but inaccurate.

The settings above are not optimal, but should get your case running.
Of course, you should never run your WMLES on a wall-resolving mesh.
Instead, we recommend using a meshing strategy presented in :ref:`grid-construction`.

Miscellaneous tips
------------------

* In regions where the TBL is attached, set :code:`hSampler` to be the distance to the second consecutive off-the-wall cell centre.
  In other regions, set it to 0, i.e. sample from the wall-adjacent cell.
* Use a mildly diffusive numerical scheme, e.g. :code:`LUST`. Tips regarding what other schemes worked well are welcome :).
* The WALE model is a good first choice for SGS modelling. Don't use implicit LES on a WMLES mesh.
* If you use :math:`h = 0`, use an algebraic wall model in integral formulation, i.e. the :code:`LOTWWallModel` with e.g.
  the :code:`IntegratedReichardt` law.

Publically available cases
--------------------------

There is a number of cases that use the library available on the web.
These can serve as good examples on how to setup your simulation!

- WMLES of channel and flat-plate TBL flow on unstructured grids.
  https://doi.org/10.6084/m9.figshare.12482438.v2 
- WMLES of channel flow and flow over a backward-facing step.
  https://doi.org/10.6084/m9.figshare.6790013.v1 
- WMLES of a flat-plate TBL using unstructured grids. 
  https://doi.org/10.6084/m9.figshare.6061298.v2 
