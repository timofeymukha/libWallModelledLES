# README #

libWallModelledLES is a library based on OpenFOAM® technology, extending the
capabilities of OpenFOAM in the area of wall-modelled LES (WMLES). This is a
turbulence modelling methodology, which allows to make LES cheaper by not
resolving the inner region of turbulent boundary layers.

If you use the library, please cite the following publication. This is also a
good source for understanding the theory behind the models.

https://doi.org/10.1016/j.cpc.2019.01.016

**This offering is not approved or endorsed by OpenCFD Limited, producer and
distributor of the OpenFOAM software via www.openfoam.com, and owner of the
OPENFOAM(R) and OpenCFD(R) trademarks.**

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

See "Installation" section on the documentation portal. In short: the latest ESI
versions should work, Foundation version 7  and below should work.

## Getting help

Please first read the troubleshooting section in the documentation.
If that does not help, please open [an issue on Github](https://github.com/timofeymukha/libWallModelledLES/issues)!

## Where this code lives
This code is available on several public repositories:
- [Github](https://github.com/timofeymukha/libWallModelledLES) --- the main repository, where all the development happens, and where you should open issues to get help.
- [Bitbucket](https://bitbucket.org/lesituu/libwallmodelledles/) --- mirror, which only gets update upon new releases.

- [Gitlab](https://gitlab.com/chalmers-marine-technology/libwallmodelledles) --- mirror, which only gets updated upon new releases.

## Published works using the library

If your work is missing from this glorious list and you want it here, [open an
issue](https://github.com/timofeymukha/libWallModelledLES/issues)!

### Journal papers
1. Mukha, T., Rezeeiravesh, S., & Liefvendahl, M. (2019). A library for
   wall-modelled large-eddy simulation based on OpenFOAM technology. Computer
   Physics Communications. https://doi.org/10.1016/j.cpc.2019.01.016
2. Rezaeiravesh, S., Mukha, T., & Liefvendahl, M. (2019). Systematic study of
   accuracy of wall-modeled large eddy simulation using uncertainty
   quantification techniques. Computers & Fluids.
   https://doi.org/10.1016/j.compfluid.2019.03.025
3. Mukha, T. (2019). The effect of numerical dissipation on the predictive
   accuracy of wall-modelled large-eddy simulation. Trudy ISP RAN / Proceedings
   of ISP RAS. https://doi.org/10.15514/ISPRAS-2019-31(6)-11
4. Mukha, T., Bensow, R. E., & Liefvendahl, M. (2021). Predictive accuracy of
   wall-modelled large-eddy simulation on unstructured grids. Computers &
   Fluids, 221, 104885. https://doi.org/10.1016/j.compfluid.2021.104885
5. Ren, X., Su, H., Yu, H.-H., & Yan, Z. (2022). Wall-modeled large eddy
   simulation and detached eddy simulation of wall-mounted separated flow via
   OpenFOAM. Aerospace, 9(12), 759. https://doi.org/10.3390/aerospace9120759
6. He, K., Zhou, F., Zhao, W., et al. (2023). Numerical analysis of turbulent
   fluctuations around an axisymmetric body of revolution based on wall-modeled
   large eddy simulations. Journal of Hydrodynamics, 35, 1041–1051.
   https://doi.org/10.1007/s42241-024-0077-8
7. Chen, S., Yang, L., Zhao, W., et al. (2023). Wall-modeled large eddy
   simulation for the flows around an axisymmetric body of revolution. Journal
   of Hydrodynamics, 35, 199–209. https://doi.org/10.1007/s42241-023-0026-y
8. Taghvaei, M., & Amani, E. (2023). Wall-modeled large-eddy simulation of
   turbulent non-Newtonian power-law fluid flows. Journal of Non-Newtonian Fluid
   Mechanics, 322, 105136. https://doi.org/10.1016/j.jnnfm.2023.105136
9. Jiang, P., Liao, S., & Xie, B. (2024). Large-eddy simulation of flow noise
   from turbulent flows past an axisymmetric hull using high-order schemes.
   Ocean Engineering, 312(Part 2), 119150.
   https://doi.org/10.1016/j.oceaneng.2024.119150
10. Fazeli, M., Emdad, H., Alishahi, M. M., & Rezaeiravesh, S. (2024).
    Wall-modeled large eddy simulation of 90° bent pipe flows with/without
    particles: A comparative study. International Journal of Heat and Fluid
    Flow. https://doi.org/10.1016/j.ijheatfluidflow.2023.109268
11. Nuca, R., Mukha, T., & Parsani, M. (2025). Explicit formulations of widely
    used wall models for large-eddy simulation. Physics of Fluids, 37, 035215.
    https://doi.org/10.1063/5.0253882
12. Mayoral, S., & Massis, A. (2025). A numerical investigation of the
    longitudinal vortex pair structure in underbody diffuser flows. Journal of
    Fluids Engineering, 147(7), 071105. https://doi.org/10.1115/1.4068036
13. Cato, A. S., Kozul, M., & Sandberg, R. (2026). Development of explicit
    algebraic LES wall models using consistent CFD-driven machine learning.
    International Journal of Heat and Fluid Flow, 117(Part B), 110157.
    https://doi.org/10.1016/j.ijheatfluidflow.2025.110157
14. Hansen, C., Yang, X. I. A., & Abkar, M. (2026). Wall-modeled large eddy
    simulation of turbulent smooth body separation using the OpenFOAM flow
    solver. Journal of Fluids Engineering, 148(1), 011502.
    https://doi.org/10.1115/1.4069033
15. He, K., Zhou, F., Zhang, J., & Wan, D. (2025). Wall-Modeled Large Eddy
    Simulation of Turbulent Boundary Layer Flows over an Axisymmetric Body of
    Revolution. International Journal of Offshore and Polar Engineering, 35(04),
    407–414. https://onepetro.org/IJOPE/article-abstract/35/04/407/794625

### Conference proceedings
1. Mukha, T., Rezaeiravesh, S., & Liefvendahl, M. (2017). An OpenFOAM library
   for wall-modelled Large-Eddy Simulation. In Proceedings of the 12th OpenFOAM
   Workshop, Exeter, UK.
2. Mukha, T., Johansson, M., & Liefvendahl, M. (2018). Effect of wall-stress
   model and mesh-cell topology on the predictive accuracy of LES of turbulent
   boundary layer flows. In 7th European Conference on Computational Fluid
   Dynamics, Glasgow, UK.
3. Mukha, T., Rezaeiravesh, S., & Liefvendahl, M. (2018). Wall-modelled
   large-eddy simulation of the flow over a backward-facing step. In Proceedings
   of the 13th OpenFOAM Workshop, Shanghai, China.
4. Liefvendahl, M., & Johansson, M. (2018). Wall-modeled LES for ship
   hydrodynamics in model scale. In Proceedings of the 32nd Symposium on Naval
   Hydrodynamics, Hamburg, Germany.
5. Malkus, T., & Belloni, C. (2020). Wall-modeled large-eddy simulations of
airfoil trailing edge noise. In Proceedings of the 8th ESI OpenFOAM Conference.
https://www.esi-group.com/sites/default/files/resource/other/1682/8th_OpenFOAM_Conference_Ohio_State_University_Malkus.pdf
6. Mayoral, S., & Massis, A. (2023). Wall-modeled large eddy simulation of flow
   past an Ahmed body with a 25° slant angle. In Proceedings of the ASME 2023
   International Mechanical Engineering Congress and Exposition, Volume 9:
   Fluids Engineering. ASME. https://doi.org/10.1115/IMECE2023-113847
7. He, K., Zhou, F., Zhao, W., & Wan, D. (2024). Wall-modeled large eddy
   simulation for a highly decelerated axisymmetric turbulent boundary layer. In
   ISOPE International Ocean and Polar Engineering Conference. ISOPE-I-24-273.
   https://onepetro.org/ISOPEIOPEC/proceedings-abstract/ISOPE24/ISOPE24/ISOPE-I-24-273/546645

### Theses
1. Bezinge, G. (2018). Wall-unresolved large eddy simulation of turbulent flow
   at high Reynolds number: Performance and computational cost investigation
   [Master’s thesis]. Department of Mathematics, University of Wyoming / Laramie
   Institute of Fluid Dynamics, ETH Zurich.
2. Mukha, T. (2018). Modelling techniques for large-eddy simulation of
   wall-bounded turbulent flows [Doctoral thesis, Uppsala University]. Acta
   Universitatis Upsaliensis.
3. Rezaeiravesh, S. (2018). Effect of grid resolution on large eddy simulation
   of wall-bounded turbulence [Doctoral thesis, Uppsala University]. Acta
   Universitatis Upsaliensis.
4. Fernandez, I. J. (2023). Explicit wall model for LES of turbulent flows over
   a smooth separated body [Master’s thesis, California State University,
   Fullerton]. ScholarWorks.
5. Charisoudis, N. (2023). Wall-resolved and wall-modelled LES for a NACA4412
   wing profile in OpenFOAM [Master’s thesis, KTH Royal Institute of
   Technology]. DiVA.
6. Sikirica, A. (2025). Adaptive mesh refinement for computationally efficient
   large eddy simulations [Doctoral thesis, University of Rijeka, Faculty of
   Engineering]. Repository of the Faculty of Engineering, University of Rijeka.
   https://repository.riteh.uniri.hr/en/object/riteh:5124