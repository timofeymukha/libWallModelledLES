.. _sampling:

Sampling
========

Sampling refers to collecting the relevant flow quantities for input into the wall model.
The main parameter here is the distance to the point where the sampling occurs, :math:`h`.
The library provides means for setting the value of the distance individually for each wall face.

Selecting the Sampler
---------------------

The library provides two samplers, which are selected using the :code:`sampler` keyword in the wall model's dictionary.
The alternatives are :code:`Tree` and :code:`Crawling`, the former being the default setting.
The :code:`Tree` sampler uses the :code:`indexedOctree` class available in OpenFOAM in order to search for the cell
containing any given point.
For each face, the point to search for is computed by following the face normal for a distance :math:`h` prescribed by
the user for this face.
The advantage of the :code:`Tree` sampler is that it does not rely on mesh structure in any way.
**NB**: the :code:`Tree` sampler is known to crash :code:`reconstructPar`, when :code:`cyclicAMI` boundaries are present
in the case.

The :code:`Crawling` sampler, by contrast, relies on the presence of prismatic mesh layers adjacent to the wall.
Starting from the wall face, it searches for the opposite face within the current cell, and continues going
in the wall-normal direction cell by cell, until the cell centre of the cell it is currently in is located at a distance
:math:`\geq h`. 
The :code:`Crawling` sampler is the recommended choice if prismatic layers are present in the mesh, because it is
significantly faster for large cases, and allows more flexibility in how to prescribe :math:`h`, see next section for
details.

There is one situation in which the :code:`Tree` and :code:`Crawling` samplers behave inconsistently: when the distance
:math:`h`is set too large and the point lies outside the domain.
The :code:`Tree` sampler will revert to using the wall-adjacent cell in this case.
The :code:`Crawling` sampler will simply continue to crawl up cell by cell until it it hits a patch (typically a boundary
between processors) and then take the last valid cell's centre to use for sampling.

When solution data is written to disk, the library will write out a field called :code:`SamplingCells`, which can be
examined in order to see, which cells are used for sampling.
The default value of this field is -1, but the cells, which are used for sampling the value will be set to the index of
the corresponding patch.
We encourage the users to examine :code:`SamplingCells` to confirm that the cell selection worked as expected.


Prescribing :math:`h`
---------------------

The values of :math:`h` should be set in a field called :code:`hSampler`.
The setting of appropriate values is done in the same way as for any other solution field.
To that end, at the wall boundaries, the boundary condition for :code:`hSampler` should be set to :code:`fixedValue`.
The desired values are then either set for the whole patch using the :code:`uniform` keyword or alternatively prescribed
on a face-by-face basis using an OpenFOAM list.
On other boundaries of type :code:`patch`, the :code:`zeroGradient` boundary condition can be used.

When the :code:`Tree` sampler is used, the values in the :code:`hSampler` will be interpreted as the desired distance to the
sampling point.
Note, that the value 0 is reserved for sampling from the wall-adjacent cell. 
By default, the same will be done by the :code:`Crawling` sampler, however alternatively one can let the sampler interpret
the set values as the index of the consecutive off-wall cell, from which to do the sampling.
So, for example, :code:`hSampler` equal to 2 will refer to sampling from the second off-wall cell. 
In order to do this, the :code:`hIsIndex` keyword should be set to :code:`yes` in the dictionary of the wall model in
:code:`nut`.

Interpolation
-------------

By default, the wall model input will be sampled from the cell center of the cell found by the sampler.
This means that if :math:`h` is prescribed as distance, the actual sampling distance may be different.
To avoid that, it is possible to interpolate the field values within the cell.
This functionality is `built into OpenFOAM <https://develop.openfoam.com/Development/openfoam/-/tree/master/src/finiteVolume/interpolation/interpolation>`_,
and several interpolation algorithms are available.
See also `This OpenFOAM Wiki article <https://openfoamwiki.net/index.php/OpenFOAM_guide/Interpolation_(by_cell)>`_.
Perhaps the most suitable choice is `cellPointFace` and also `cellPoint`.
Note that *interpolation within the wall-adjacent cell is not possible*.

A summary of the parameters pertaining to sampling are given below.

.. code-block::
  :linenos:

  sampler Tree; //Crawling
  hIsIndex 0; // 1
  interpolation cell; // cellPoint, cellPointFace ...

