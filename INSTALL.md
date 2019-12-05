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
- `DLA_ALL_TESTS_USE_RUNNER` default `OFF`: If `ON` single rank tests are invoked with `TEST_RUNNER` too.

If compiler wrappers which include mpi are used (e.g. `mpicc` and `mpicxx`, Cray's `cc` and `CC`)
no libraries and include paths are added and the following default are used:
- `TEST_RUNNER = mpirun`
- `TEST_RUNNER_NP_OPT = -n`

Otherwise `FindMPI` is used to find the MPI libraries, and the following defaults are used:
- `TEST_RUNNER = ${MPIEXEC}`
- `TEST_RUNNER_NP_OPT = ${MPIEXEC_NUMPROC_FLAG}`

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

- `DLA_WITH_FORTRAN` default `ON` (`OFF` if ScaLAPACK is not available): If `ON` builds the fortran iso C interface. (Only available if ScaLAPACK is enabled)
- `DLA_PRINT_DEBUG_INFO` default `OFF`: If `ON` prints some extra debug informations. (E.g. thread binding.)
- `DLA_PRINT_DEBUG_CALL_PARAM` default `OFF`: If `ON` DLA interface routines print the call arguments.
- `DLA_COVERAGE_TEST` default `OFF`: If `ON` enables coverage test (Requires GCC). The library is built with the `--coverage` flags and enables `make coverage` which make the coverage test summary with lcov.


# Testing

`make test` execute unit tests to check the correctnes of the installation.

## Coverage testing

To run the coverage test:
- Configure the library with `DLA_COVERAGE_TEST=ON`
- Build the library and the tests: `make`
- Execute the unit tests: `make test`
- Generate the summary: `make coverage`
