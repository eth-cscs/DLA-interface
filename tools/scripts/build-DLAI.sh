#!/usr/bin/bash

# Variables
BUILDROOT=${BUILDROOT:-$(cd ${0%/*} && echo ${PWD%/tools/scripts})}
MAKE_JOBS=${MAKE_JOBS:-$(nproc)}

if [ ! -d ${BUILDROOT}/BUILD ]; then
    mkdir ${BUILDROOT}/BUILD
fi

cd ${BUILDROOT}/BUILD

cmake .. \
   -DCMAKE_BUILD_TYPE=Debug \
   -DDLA_LAPACK_TYPE:STRING="Custom" \
   -DDLA_BLAS_LAPACK_LIB="${BUILDROOT}/OpenBLAS/BUILD/lib/libopenblas.so" \
   -DDLA_BLAS_LAPACK_INC="${BUILDROOT}/OpenBLAS/BUILD/include" \
   -DDLA_SCALAPACK_TYPE="Custom" \
   -DDLA_SCALAPACK_LIB="${BUILDROOT}/ScaLAPACK/BUILD/lib/libscalapack.so" \
   -DDLA_WITH_TEST="ON" \
   -DDLA_COVERAGE_TEST="ON" \
   -DDLA_WITH_FORTRAN="ON"

make -j${MAKE_JOBS}
make test