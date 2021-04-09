#!/bin/sh -x
set -e


# Variables
BUILDROOT=${BUILDROOT:-$(cd ${0%/*} && echo ${PWD%/tools/scripts})}
MAKE_JOBS=${MAKE_JOBS:-$(nproc)}

DOWNLOAD_DIR=${BUILDROOT}/download

# Package variables
VERSION=v0.3.7
GIT=git@github.com:xianyi/OpenBLAS.git
SRC_DIR="${BUILDROOT}/OpenBLAS"

# Clone sources
if [ ! -d ${SRC_DIR} ]; then
    git clone --branch ${VERSION} --single-branch ${GIT} ${SRC_DIR}
fi


cd ${SRC_DIR}

# Build
if [ ! -e ${SRC_DIR}/.built ]; then
    # Support multiple architectures with dynamic arch.
    DYNAMIC_ARCH=${DYNAMIC_ARCH:-1}
    DYNAMIC_ARCH=${DYNAMIC_ARCH} make -j${MAKE_JOBS}
    touch ${SRC_DIR}/.built
fi

# Install to OpenBLAS/BUILD
make install PREFIX=${BUILDROOT}/OpenBLAS/BUILD