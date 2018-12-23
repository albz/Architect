# ![Architect Logo](./logo/logo.png) rchitect

[![Build Status Master](https://travis-ci.org/albz/Architect.svg?branch=master)](https://travis-ci.org/albz/Architect "master")

-------

`Architect` _/ˈɑːkɪtɛkt/_ a hybrid code for Beam Driven Plasma Acceleration.

-------

Architect is a time explicit serial hybrid code design to perform quick simulations for beam driven plasma wakefield acceleration. The code has been mainly used and tested (versus experimental results) with electron beams but it could evolve a bunch with arbitrary charge. Architect uses a _hybrid approach_: relativistic bunches are treated kinetically as in a Particle-in-cell (PIC) code while the background plasma is modelled as a fluid.

The _fluid_ and _Maxwell's_ equations are solved in cylindrical symmetry to further reduce the computational costs. Bunch particles can, instead, be solved either in a 6D phase-space or with weighted particles in cylindrical symmetry.

The fluid approach has shown to be a powerful physical-mathematical approach that by drastically reducing the computational costs allows for a larger number of runs, or eventually for longer simulations.

Some of these results, schemes, intuitions are being transported in the full PIC code ALaDyn.

-------

Architect’s development depends on its visibility from publications or presentations featuring its usage or results. Please, we gently invite you, when
  + publishing an article
  + publishing a conference proceedings
  + publishing any manuscripts that had benefited from Architet's results at any levels, from preliminary runs to more important results
  + conference presentations
  + office presentations
  + private communications

to acknowledge the code usage by citing the two following articles:
_[1]_ **Efficient modeling of plasma wakefield acceleration in quasi-non-linear-regimes with the hybrid code Architect**, A. Marocchino, F. Massimo, A. R. Rossi, E. Chiadroni, M. Ferrario - _Nuclear Instruments and Methods in Physics Research A (2016)_ http://dx.doi.org/10.1016/j.nima.2016.03.005

_[2]_ **Comparisons of time explicit hybrid kinetic-fluid code Architect for Plasma Wakefield Acceleration with a full PIC code** F. Massimo, S. Atzeni and A. Marocchino - _Journal of Computational Physics 327 (2016) 841–850_
http://dx.doi.org/10.1016/j.jcp.2016.09.067

_[3]_ and eventually the _Zenodo DOI_ of this latest release, v1.0.0-alpha: [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.49572.svg)](http://dx.doi.org/10.5281/zenodo.49572)


If help or changes in the code were obtained from Architect's main developers (A. Marocchino and F. Massimo), please acknowledge their participation in any subsequent publication or presentation.

If your publication makes significant use of Architect, we will gladly list it in the Publications section.

---

This newer version is released as an alpha version as is, it is distributed in the hope that it will be useful
but without any warranty. It will be maintained here on GitHub.

Architect is protected by the open-source Gnu GPL 3 license. Copyright on the code is by A. Marocchino, F. Massimo
