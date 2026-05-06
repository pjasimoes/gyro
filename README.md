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

## Build the Binary

From the project folder, run:

```bash
make all
```

This generates the executable binary `gyro` in the same folder.

On macOS, the default `g++` is often Apple Clang and may not support OpenMP.
The current `makefile` handles this automatically and falls back to a build
without `-fopenmp` when needed.

Other available targets:

```bash
make single
make debug
make clean
```

## Run the Python Code

After building the binary, run:

```bash
python3 gyro.py
```

The Python wrapper expects the compiled `gyro` binary to be in the same folder
as `gyro.py`.
