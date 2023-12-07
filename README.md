# README #

libWallModelledLES is a library based on OpenFOAM® technology, extending the capabilities of OpenFOAM in the area of wall-modelled LES (WMLES).
This is a turbulence modelling methodology, which allows to make LES cheaper by not resolving the inner region of turbulent boundary layers.

If you use the library, please cite the following publication. This is also a good source for understanding the theory behind the models.

https://doi.org/10.1016/j.cpc.2019.01.016

**This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.**

## News ##

- **2023-04-28** Version 0.7.0 released.
- **2023-04-25** Development moves to [Github](https://github.com/timofeymukha/libWallModelledLES/), Bitbucket remains as a mirror.
- **2023-01-05** Version 0.6.1 released.
- **2021-08-30** Version 0.6.0 released.
- **2019-10-28** Version 0.5.1 released.
- **2019-08-01** Version 0.5.0 released.
- **2019-02-23** Version 0.4.1 released, containing a small bugfix.
- **2018-11-17** Version 0.4.0 released, see CHANGELOG.md for list of changes.

## Documentation
[https://libwmles.readthedocs.io](https://libwmles.readthedocs.io)

## Compatibility ##

See "Installation" section on the documentation portal. In short: the latest ESI versions should work, Foundation version 7  and below should work.

## Getting help

Please first read the troubleshooting section in the documentation.
If that does not help, please open [an issue on Github](https://github.com/timofeymukha/libWallModelledLES/issues)!

## Where this code lives
This code is available on several public repositories:
- [Github](https://github.com/timofeymukha/libWallModelledLES) --- the main repository, where all the development happens, and where you should open issues to get help.
- [Bitbucket](https://bitbucket.org/lesituu/libwallmodelledles/) --- mirror, which only gets update upon new releases.
 
- [Gitlab](https://gitlab.com/chalmers-marine-technology/libwallmodelledles) --- mirror, which only gets updated upon new releases.

## Published works using the library

If your works is missing from this glorious list and you want it here, [open an issue](https://github.com/timofeymukha/libWallModelledLES/issues)!


- Mukha, T., Rezaeiravesh, S., & Liefvendahl, M. (2017). An OpenFOAM library for wall-modelled Large-Eddy Simulation. In proceedings of the 12th OpenFOAM Workshop, Exeter, UK.
- Mukha, T., Johansson, M., & Liefvendahl, M. (2018). Effect of wall-stress model and mesh-cell topology on the predictive accuracy of LES of turbulent boundary layer flows.
  In 7th European Conference on Computational Fluid Dynamics, Glasgow, UK.
- Mukha, T., Rezaeiravesh, S., & Liefvendahl, M. (2018). Wall-modelled large-eddy simulation of the flow over a backward-facing step. In proceedings of 13th OpenFOAM Workshop, Shanghai, China. Shanghai, China.
- Liefvendahl, M., & Johansson, M. (2018). Wall-Modeled LES for Ship Hydrodynamics in Model Scale. In proceedings of the 32nd Symposium on Naval Hydrodynamics, Hamburg, Germany.
- Bezinge, G. (2018) Wall-unresolved large eddy simulation of turbulent flow at high Reynolds number: Performance and computational cost investigation. Master thesis.
  Department of Mathematics University of Wyoming (UW) Laramie Institute of Fluid Dynamics Swiss Federal Institute of Technology (ETH) Z�rich.
- Mukha, T., Rezeeiravesh, S., & Liefvendahl, M. (2019). A library for wall-modelled large-eddy simulation based on OpenFOAM technology. Computer Physics Communications. DOI: 10.1016/j.cpc.2019.01.016. Preprint: https://arxiv.org/abs/1807.11786
- Rezaeiravesh, S., Mukha, T., & Liefvendahl, M. (2019). Systematic study of accuracy of wall-modeled large eddy simulation using uncertainty quantification techniques. Computers & Fluids. DOI: 10.1016/j.compfluid.2019.03.025. Preprint: https://arxiv.org/abs/1810.05213
- Mukha, T. (2019) The effect of numerical dissipation on the predictive accuracy of wall-modelled large-eddy simulation. Trudy ISP RAN/Proc. ISP RAS. DOI: 10.15514/ISPRAS-2019-31(6)-11
- Malkus, T., &  Belloni, C. (2020) Wall-modeled large-eddy simulations of airfoil trailing edge noise. In proceedings of the 8th ESI OpenFOAM Conference. URL: https://www.esi-group.com/sites/default/files/resource/other/1682/8th_OpenFOAM_Conference_Ohio_State_University_Malkus.pdf
- Mukha, T., Bensow, R.E., Liefvendahl, M., 2021. Predictive accuracy of wall-modelled large-eddy simulation on unstructured grids. Comput. Fluids 221, 104885. https://doi.org/10.1016/j.compfluid.2021.104885 (Open access)
