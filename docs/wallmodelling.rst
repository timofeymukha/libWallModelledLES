
Wall modelling theory
=====================

General considerations
----------------------

Generally, WMLES implies modelling of the dynamics of the inner layer of turbulent boundary layers instead of resolving
them with the grid.
This entails significant computational savings at high Reynolds numbers, compared to performing a wall-resolved LES.
A word of warning: **WMLES is not cheap**.
Quite the opposite, WMLES is still expensive and is, for that reason, only starting to be used in applications,
primarily in aerodynamical.

Note that the definition is rather general in the sense that it doesn't define how the inner layer is to be modelled.
On the other hand, it is quite specific since it applies only to attached boundary layers with a clear separation
between outer and inner scales.
This library only provides functionality for so-called *wall stress modelling*.
This is a particular type of WMLES, in which the lack of inner layer resolution is compensated by providing the correct
value of the wall shear stress.
Providing this value is the task of the wall model.
Wall stress modelling fits naturally into the numerical framework of OppenFOAM, i.e. finite volume discretization.

A schematic of how a wall-stress model works is shown in :numref:`fig-model-scheme`.
For each wall face we consider some point away from the wall, located at a distance :math:`h`.
At this point, we can sample the solution to the LES equations, for example, the velocity and pressure gradient vectors.
These sampled values are given to the wall model, which uses them to predict the wall shear stress :math:`\bar \tau_w`.

.. _fig-model-scheme:

.. figure:: /figures/wm_scheme.png
   :align: center
   :width: 75%

   Modus operandi of a wall-stress model.

This library provides associated functionality:

* Possibility to define an arbitrary :math:`h` for each face.
* Sampling relevant fields and temporal averaging of the sampled values.
* A set of wall models implemented as boundary conditions for :code:`nut`.
  

Further reading
---------------
* Paper presenting this library :cite:`Mukha2019`. *Please cite this if you use the library*.
* Review papers on WMLES :cite:`Larsson2016`, :cite:`Bose2018a`.
* Youtube playlist on WMLES theory: https://www.youtube.com/playlist?list=PLrwFJPCcTaPWl00ZE04a2o2A99NmF7-lt
