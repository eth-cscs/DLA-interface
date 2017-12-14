# Dependencies

- BLAS/LAPACK
- MPI
- Boost: program options library
- ScaLAPACK (Optional)
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

## MPI

Options:
- `TEST_RUNNER` The command used to run the tests. (e.g. mpirun, srun, ...)
- `TEST_RUNNER_NP_OPT` The name of the option used to specify the number of ranks.

If compiler wrappers which include mpi are used (e.g. `mpicc` and `mpicxx`, Cray's `cc` and `CC`)
no libraries and include paths are added and the following default are used:
- `TEST_RUNNER = mpirun`
- `TEST_RUNNER_NP_OPT = -n`

Otherwise `FindMPI` is used to finde the MPI libraries, and the following defaults are used:
- `TEST_RUNNER = ${MPIEXEC}`
- `TEST_RUNNER_NP_OPT = ${MPIEXEC_NUMPROC_FLAG}`

## Boost program options

`find_package(Boost)` is used.
`BOOST_ROOT` can be used to specify the Boost directory.

## Scalapack

BLAS and LAPACK can be configured in 3 ways setting `DLA_LAPACK_TYPE`
- `"None"`: ScaLAPACK is disabled.
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

# Testing

`make test` execute unit tests to check the correctnes of the installation.
