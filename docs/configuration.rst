.. _configuration:

Configuration
=============

There are some shared options that apply to all all models.
The most important ones relate to the sampling of the data to the wall model,
which are discussed in a separate section, see :ref:`sampling`.
Here, we list the others.


- :code:`silent`. This controls the verbosity of the wall model in the log file.
  Currently, if one sets this to 1 the wall model will no longer print information
  on time consumption.
- :code:`copyToPatchInternalField`.
  If true, the wall values of :code:`nut` will be copied to the wall-adjacent 
  cell.
  This promotes the stability of the model.

  

