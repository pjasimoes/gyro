# gyro
gyrosynchrotron radiation code
Full equations from Ramaty ApJ 158 1969
http://adsabs.harvard.edu/abs/1969ApJ...158..753R
This version is based on IDL version of Ramaty's code original FORTRAN.
Improvements:
added terms for anisotropic pitch-angle distribution;
Gauss-Legendre integration;
arbitrary distributions in energy and pitch-angle (array input);
OpenMP implementation;
Wild-Hill approximations for Bessel functions; http://adsabs.harvard.edu/abs/1971AuJPh..24...43W
Includes free-free emission/absorption;
