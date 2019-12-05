#!/bin/sh -x
set -e


# Variables
BUILDROOT=${BUILDROOT:-$(cd ${0%/*} && echo ${PWD%/tools/scripts})}
MAKE_JOBS=${MAKE_JOBS:-$(nproc)}

DOWNLOAD_DIR=${BUILDROOT}/download

# Package variables
VERSION=v2.1.0
GIT=git@github.com:Reference-ScaLAPACK/scalapack.git
SRC_DIR="${BUILDROOT}/ScaLAPACK"


# Clone sources
if [ ! -d ${SRC_DIR} ]; then
    git clone --branch ${VERSION} --single-branch ${GIT} ${SRC_DIR}
fi

if [ ! -e ${SRC_DIR}/.configured ]; then
    mkdir -p ${SRC_DIR}/BUILD
    cd ${SRC_DIR}/BUILD

    # Add the OpenBLAS library so that ScaLAPACK can find it
    CPATH=${BUILDROOT}/OpenBLAS/BUILD/include:${CPATH} \
    LD_LIBRARY_PATH=${BUILDROOT}/OpenBLAS/BUILD/lib:${LD_LIBRARY_PATH} \
    LIBRARY_DIR=${BUILDROOT}/OpenBLAS/BUILD/lib:${LIBRARY_DIR} \
    PKG_CONFIG_PATH=${BUILDROOT}/OpenBLAS/BUILD/lib/pkgconfig:${PKG_CONFIG_PATH} \
    cmake .. -DUSE_OPTIMIZED_BLAS=1 -DBUILD_SHARED_LIBS=On
    touch ${SRC_DIR}/.configured
fi


cd ${SRC_DIR}/BUILD

# Build
if [ ! -e ${SRC_DIR}/.built ]; then
    make -j${MAKE_JOBS}
    touch ${SRC_DIR}/.built
fi
