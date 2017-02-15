# SSGW: Periodic Surface Gravity Waves

This Matlab code computes irrotational 2D periodic steady surface pure gravity waves of arbitrary length in arbitrary depth. The formulation is based on the so-called Babenko equation and pseudo-spectral discretization in the conformal domain. The resulting equation is solved using Petviashvili iteration method.

MANDATORY INPUT PARAMETERS:
* kd = k*d : relative depth (wavenumber "k" times mean water depth "d").
* kH2 = k*H/2 : steepness (half the total wave height "H" times the wavenumber "k").

OPTIONAL INPUT PARAMETERS:
* N : number of positive Fourier modes (default, N=2048).
* tol : tolerance (default, tol=1e-14).

OUTPUT PARAMETERS:
* zs = complex abscissas at the free surface (at the computational nodes).
* ws = complex velocity at the free surface (at the computational nodes).
* PP = Physical Parameters: PP(1)=depth, PP(2)=wavenumber, PP(3)=wavelenght,
* PP(4)=celerity c_e, PP(5)=celerity c_s, PP(6)=Bernoulli constant,
* PP(7)=crest height, PP(8)=trough height, PP(9)=impulse,
* PP(10)=potential energy, pp(11)=kinetic energy, PP(12)=radiation stress,
* PP(13)=momentum flux, PP(14)=energy flux, PP(16)=group velocity.

For more information on the method, please, read our preprint:

* D. Clamond & D. Dutykh. [Accurate fast computation of steady two-dimensional surface gravity waves in arbitrary depth.](https://hal.archives-ouvertes.fr/hal-01465813/) Submitted, 2017