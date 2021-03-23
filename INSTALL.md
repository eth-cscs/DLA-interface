# Dependencies

- BLAS/LAPACK
- MPI
- ScaLAPACK
- Parsec/DPlasma (Optional)

# Cmake options:

## BLAS / LAPACK

BLAS and LAPACK can be configured in 3 ways setting `DLA_LAPACK_TYPE`
- `"Compiler"`: The compiler already link with the correct libraries.
- `"MKL"`: MKL located in `MKL_ROOT` (or env variable `MKLROOT` if not defined) is used.
- `"Custom"`: Link the libraries specified by `DLA_BLAS_LAPACK_LIB`.

For the MKL case the following options apply:
- `MKL_THREADING`:
  - `"Sequential"`
  - `"GNU OpenMP"` (Default)
  - `"Intel OpenMP"` (Default with Apple)

## Scalapack

ScaLAPACK can be configured in 3 ways setting `DLA_SCALAPACK_TYPE`
- `"Compiler"`: The compiler already link with the correct libraries
- `"MKL"`: MKL ScaLAPACK is used (Only if `DLA_LAPACK_TYPE="MKL"` as well).
- `"Custom"`: Link the libraries specified by `DLA_SCALAPACK_LIB`

For the MKL case the MPI interface has to be defined iwith `MKL_MPI_TYPE` as well:
- `"IntelMPI"`: Supports IntelMPI, MPICH, MVAPICH
- `"OpenMPI"`: Supports OpenMPI

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
