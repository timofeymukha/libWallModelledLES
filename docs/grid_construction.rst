.. _grid-construction:

Grid construction
=================

The fact that the inner layer is modelled has a very important implication for the grid construction strategy.
**The grid size is no longer controlled by the viscous length scale**.
In other words, we are no longer interested in the size of the grid in "plus"-units.
This includes the distance to the first off-wall cell centres and the associated :math:`y^+` values. 
This is quite unique to WMLES and people are so used to measuring :math:`y^+` that this sometimes becomes difficult
to accept!

Since in WMLES the aim is to resolve the outer layer only, the relevant length scale for grid construction is the
boundary layer thickness :math:`\delta`.
The cells size is thus conveniently defined as the number of cells per :math:`\delta`.
Alternatively, sometimes the number of cells per a :math:`\delta^3`-cube is used.
Typical cell sizes that can be found in the literature are :math:`\delta/15`-:math:`\delta/30`.
Results from channel flow simulations in this span of cell sizes can be found in :cite:`Mukha2019` and :cite:`Rezaeiravesh2019a`.
In order to construct the mesh in relation to :math:`\delta`, it's value needs to be known a priori.
The recommendation is to use a RANS simulation to that end. 

Some works suggest different resolutions for the streamwise and spanwise wall-parallel directions.
However, here we assume a complex geometry flow is simulated, meaning that these directions are not readily distinguished
at the pre-processing stage. 
There is also no strict rules regarding wether to refine the mesh towards the wall.
For example, in :cite:`Kawai2012`, the authors argue that the different refinement levels should be tested for the
wall-normal spacing below :math:`h`, i.e the distance to the sampling point.
The predictions of the wall shear stress should eventually converge upon successive refinement.
However, uniform spacing in the wall normal direction has been successfully used, see e.g. :cite:`Mukha2019`
:cite:`Rezaeiravesh2019a`.
Our subjective recommendation is therefore to use an isotropic grid, aiming for the same characteristic length of the
grid in all three coordinate directions.

The following grid construction strategy for meshing boundary layers has been demonstrated in :cite:`Mukha2021`.

* Estimate the distribution of :math:`\delta` by running a RANS precursor or using approximative analytical relations, when such are available.
* Construct a surface mesh with the chosen distribution of :math:`d/\delta` on the wall boundary.
  Here, :math:`d` is the characteristic size of the mesh cell.
* Extrude prismatic layers between the wall and :math:`\delta`, with the number of layers equal to :math:`\delta/d`.
  The height of each prism is thus equal to :math:`d`.
* Apply suitable discretization to the rest of the computational domain.

Flat plate turbulent boundary simulation cases with grids constructed according to this strategy can be found on
`Figshare <https://doi.org/10.6084/m9.figshare.12482438.v2>`_.

There are no concrete guidelines for how to mesh the near wall region when the boundary layer separates.
This is currently a topic of research.
In the literature, typically the grid used for the boundary layer prior to separation is extended into the region
occupied by the separation bubble.
