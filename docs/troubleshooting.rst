.. _troubleshooting

Troubleshooting
===============

The number 1 cause of the case crashing is the divergence in the Newton solver
used in algebraic models.
Look at the error trace, does crash occur in a function called :code:`value` or
:code:`derivative`? 
Then it is the Newton solver.

The divergence occurs due to a bad velocity value being sampled to the wall model.
This typically occurs due to bad (e.g. uniform) initial conditions.
The remedy is therefore to run your case for a bit without wall modelling to get
a reasonable initial field and then turn on the wall model.

However, often one can also run with :code:`hSampler` set to 0, i.e. sampling
from the near-wall cell. 
This tends to be very stable.

If your case crashes anyway, try setting :code:`copyToPatchInternalField` to 1
in the boundary condition configuration. 
This copies the wall value of the viscosity to the near-wall cell.
This a trick, but it tends to stabilize the velocity value in the near-wall cell.
This should certainly be used id you see that your wall :code:`nut` values are
starting to climb up unreasonably. 
The latter is, in fact, how ODE models "diverge", even if the case will be
able to run for a long time.

Finally, don't refine the mesh towards the wall, at least not a lot.
This tends to lead to less stable simulations and won't improve the accuracy.