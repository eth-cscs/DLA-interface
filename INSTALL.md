# Dependencies

- BLAS/LAPACK
- MPI
- ScaLAPACK
- Parsec/DPlasma (Optional)

For building DLA interface the following packages are required, if building on
a local machine with GNU toolchain.

- Debian/Ubuntu:
  - sudo apt-get install build-essential git cmake doxygen openmpi-bin
- CentOS7:
  - sudo yum groupinstall "Development Tools"
  - sudo yum install git cmake doxygen openmpi-devel

Along these packages a scalable LAPACK and BLAS/LAPACK implementations are
required. A set of scripts are prepared for compiling OpenBLAS and ScaLAPACK
and finally DLA interface using the previous packages to build the software
into directory BUILD.

  $ cd DLA-interface
  $ ./tools/scripts/build-OpenBLAS.sh
  $ ./tools/scripts/build-ScaLAPACK.sh
  $ ./tools/scripts/build-DLAI.sh

On a HPC the required packages must be loaded or installed. For CMake options
check the tools/scripts/build-DLAI.sh to see how to provide BLAS/LAPACK and
ScaLAPACK libraries.

# Cmake options:

## BLAS, LAPACK and SCALAPACK

Using MKL as provider for all of them is possible by specifying the CMake option `DLAI_WITH_MKL=on`.
See FindMKL module for more options, moreover you can select the specific variant of MKL via the variables:
- `MKL_LAPACK_TARGET`
- `MKL_SCALAPACK_TARGET`

Otherwise, if you want to use implementations provided implicitly by the compiler, or you want to use a specific implementation, you can
specify how to link to the specific implementation with the variables:
- `LAPACK_LIBRARY` (for both BLAS and LAPACK)
- `SCALAPACK_LIBRARY`

## Parsec/DPlasma

`PARSEC_DIR` has to point to the directory where `ParsecConfig.cmake` is,
e.g. `<parsec_install_prefix>/lib/cmake`.
Parsec has to be configured to build DPlasma (Need COREBLAS from Plasma-2.8.0), and the option `PARSEC_WITH_DEVEL_HEADERS=ON` should be set.

## Other options

- `DLAI_WITH_FORTRAN` default `ON` (`OFF` if ScaLAPACK is not available): If `ON` builds the fortran iso C interface. (Only available if ScaLAPACK is enabled)
- `DLAI_PRINT_DEBUG_INFO` default `OFF`: If `ON` prints some extra debug informations. (E.g. thread binding.)
- `DLAI_PRINT_DEBUG_CALL_PARAM` default `OFF`: If `ON` DLA interface routines print the call arguments.

# Testing

`make test` execute unit tests to check the correctnes of the installation.
