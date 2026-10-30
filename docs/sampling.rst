.. _sampling:

Sampling
========

Sampling refers to collecting the relevant flow quantities for input into the
wall model. The main parameter here is the distance to the point where the
sampling occurs, :math:`h`. The library provides means for setting the value of
the distance individually for each wall face.

Selecting the Sampler
---------------------

The library provides two samplers, :code:`NormalRay` and :code:`Crawling`. By
default, the sampler is selected by how :math:`h` is prescribed, see the
:code:`hIsIndex` keyword below: :code:`NormalRay` is used when :math:`h` is a
distance, and :code:`Crawling` when it is the index of the off-wall cell.
Normally, there is no need to select the sampler explicitly. It is possible,
however, using the :code:`sampler` keyword in the wall model's dictionary, in
which case the keyword is also written out with the boundary condition.

The :code:`NormalRay` sampler follows the line segment from the centre of the
wall face along the inward face normal, up to the distance :math:`h`
prescribed by the user for this face. Starting from the wall-adjacent cell, it
finds the face through which the segment leaves the current cell and continues
in the cell on the other side of that face. The sampling cell is the cell
containing the point at the distance :math:`h`. The faces are decomposed into
triangles formed by their edges and their centres for this purpose. The
:code:`NormalRay` sampler does not rely on mesh structure in any way, i.e. it
works for arbitrary polyhedral cells. Its cost is proportional to the number of
cells crossed, so it is about as fast as the :code:`Crawling` sampler. It only
supports prescribing :math:`h` as a distance.

The :code:`Crawling` sampler, by contrast, relies on the presence of prismatic
or hexahedral mesh layers adjacent to the wall. Starting from the wall face, it
searches for the opposite face within the current cell, and continues going in
the wall-normal direction cell by cell, until the opposite face of the cell it
is currently in is located at a distance :math:`> h` from the wall face, or
until the requested off-wall cell index is reached. For other types of cells,
there is no opposite face, and the crawling is stopped with a warning. On
skewed layers, the column of cells followed by the crawl can deviate from the
wall normal, which is why the :code:`NormalRay` sampler is used for distances
by default.

If the wall-normal line or the crawl reaches a boundary before the distance
:math:`h`, both samplers use the last cell before the boundary. The
:code:`NormalRay` sampler then also reduces the sampling distance to the middle
of the part of the line inside that cell, so that the sampling point remains
inside the domain. Cyclic boundaries are not crossed.

**NB**: The :code:`Tree` sampler of earlier versions has been removed. The
:code:`NormalRay` sampler finds the same cell, i.e. the cell at the distance
:math:`h` along the wall normal, and is much faster. Remove
:code:`sampler Tree;` from existing cases.

Sampling Across Processor Boundaries
------------------------------------

In parallel runs, both samplers continue the search on the neighbouring
processor whenever they reach a processor boundary. The selected sampling cells
are thus independent of the decomposition and the same as in a serial run. The
field values are interpolated on the processor holding the sampling cell and
sent to the processor holding the wall face at each time step. The price is
some additional communication, both during the setup and at each time step.

Note the following.

- Processor boundaries originating from cyclic patches are not crossed, the same
  as cyclic patches are not crossed in serial.
- Utilities that run in serial on the decomposed case, such as
  :code:`reconstructPar`, construct the boundary condition separately for each
  processor, where the processor boundaries cannot be crossed. It is
  recommended to set :code:`LIBWMLES_SKIP_SAMPLING_SETUP=true` for these, see
  below.
- This only applies to wall models that sample from a single cell. The
  :code:`MulticellLOTWWallModel` searches only among the cells of the processor
  holding the wall face, and the selected cells thus depend on the
  decomposition.

When solution data is written to disk, the library will write out a field called
:code:`samplingCells`, which can be examined in order to see, which cells are
used for sampling. The default value of this field is -1, but the cells, which
are used for sampling the value will be set to the index of the corresponding
patch. We encourage the users to examine :code:`samplingCells` to confirm that
the cell selection worked as expected.

Bypassing Sampling Setup
------------------------

Some OpenFOAM utilities, such as :code:`decomposePar`, and :code:`renumberMesh`,
construct the wall model boundary condition without evaluating it.  For very
large cases the sampling-cell search can dominate the utility runtime even
though the result is not used.  This setup can be skipped by setting the
environment variable :code:`LIBWMLES_SKIP_SAMPLING_SETUP=true` before running
the utility. **For renumberMesh setting this is mandatory**.

This switch must not be used for solver runs.  If a wall model tries to sample
while the switch is enabled, the library exits with an error instead of using
uninitialized sampling data.


Prescribing :math:`h`
---------------------

The values of :math:`h` should be set in a field called :code:`hSampler`. The
setting of appropriate values is done in the same way as for any other solution
field. To that end, at the wall boundaries, the boundary condition for
:code:`hSampler` should be set to :code:`fixedValue`. The desired values are
then either set for the whole patch using the :code:`uniform` keyword or
alternatively prescribed on a face-by-face basis using an OpenFOAM list. On
other boundaries of type :code:`patch`, the :code:`zeroGradient` boundary
condition can be used.

By default, the values in the :code:`hSampler` will be interpreted as the
desired distance to the sampling point, and the :code:`NormalRay` sampler is
used. Note, that the value 0 is reserved for sampling from the wall-adjacent
cell. Alternatively, one can let the sampler interpret the set values as the
index of the consecutive off-wall cell, from which to do the sampling. So, for
example, :code:`hSampler` equal to 2 will refer to sampling from the second
off-wall cell. In order to do this, the :code:`hIsIndex` keyword should be set
to :code:`yes` in the dictionary of the wall model in :code:`nut`. The
:code:`Crawling` sampler is then used, since only it supports cell indices.

Interpolation
-------------

By default, the wall model input will be sampled from the cell center of the
cell found by the sampler. This means that if :math:`h` is prescribed as
distance, the actual sampling distance may be different. To avoid that, it is
possible to interpolate the field values within the cell. This functionality is
`built into OpenFOAM
<https://gitlab.com/openfoam/core/openfoam/-/tree/master/src/finiteVolume/interpolation/interpolation>`_,
and several interpolation algorithms are available. See also `This OpenFOAM Wiki
article
<https://openfoamwiki.net/index.php/OpenFOAM_guide/Interpolation_(by_cell)>`_.
Perhaps the most suitable choice is `cellPointFace` and also `cellPoint`. Note
that *interpolation within the wall-adjacent cell is not possible*.

A summary of the parameters pertaining to sampling are given below.

.. code-block::
  :linenos:

  hIsIndex 0; // 1
  sampler NormalRay; // Crawling, optional, selected by hIsIndex by default
  interpolation cell; // cellPoint, cellPointFace ...
