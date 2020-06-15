ARG BASE_IMAGE=ubuntu:20.04

FROM $BASE_IMAGE

WORKDIR /root

ENV DEBIAN_FRONTEND noninteractive
ENV FORCE_UNSAFE_CONFIGURE 1

# Install basic tools
RUN apt-get update -qq && apt-get install -qq -y --no-install-recommends \
    software-properties-common \
    build-essential gfortran binutils lcov \
    git tar wget curl gpg-agent jq tzdata bison flex python python-dev pkg-config && \
    rm -rf /var/lib/apt/lists/*

# Install cmake
RUN wget -qO- "https://cmake.org/files/v3.17/cmake-3.17.0-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C /usr/local

# Install MPICH ABI compatible with Cray's lib on Piz Daint
ARG MPICH_VERSION=3.3.2
ARG MPICH_PATH=/usr/local/mpich
ENV MPICH_PATH=${MPICH_PATH}
RUN wget -q https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz && \
    tar -xzf mpich-${MPICH_VERSION}.tar.gz && \
    cd mpich-${MPICH_VERSION} && \
    ./configure \
      --prefix=$MPICH_PATH && \
    make install -j$(nproc) && \
    rm -rf /root/mpich-${MPICH_VERSION}.tar.gz /root/mpich-${MPICH_VERSION}

# Install OpenBLAS
ARG OPENBLAS_VERSION=0.3.9
ARG OPENBLAS_PATH=/usr/local/openblas
ENV OPENBLAS_PATH=${OPENBLAS_PATH}
RUN wget -qO - https://github.com/xianyi/OpenBLAS/archive/v${OPENBLAS_VERSION}.tar.gz -O openblas.tar.gz && \
    tar -xzf openblas.tar.gz && \
    cd OpenBLAS-${OPENBLAS_VERSION}/ && \
    make USE_OPENMP=1 USE_THREAD=1 USE_LOCKING=1 DEBUG=1 -j$(nproc) && \
    make install NO_STATIC=1 PREFIX=${OPENBLAS_PATH} && \
    rm -rf /root/openblas.tar.gz /root/OpenBLAS-${OPENBLAS_VERSION}/ && \
    echo ${OPENBLAS_PATH}/lib >> /etc/ld.so.conf.d/openblas.conf && \
    ldconfig

# Install ScaLAPACK
ARG SCALAPACK_VERSION=2.1.0
ARG SCALAPACK_PATH=/usr/local/scalapack
ENV SCALAPACK_PATH=${SCALAPACK_PATH}
RUN wget -qO - http://www.netlib.org/scalapack/scalapack-${SCALAPACK_VERSION}.tgz -O scalapack.tgz && \
    tar -xzf scalapack.tgz && \
    cd scalapack-${SCALAPACK_VERSION}/ && \
    mkdir build && \
    cd build && \
    CC=/usr/local/mpich/bin/mpicc CXX=/usr/local/mpich/bin/mpicxx FC=/usr/local/mpich/bin/mpif90 cmake .. \
      -DCMAKE_INSTALL_PREFIX=${SCALAPACK_PATH} \
      -DCMAKE_BUILD_TYPE=Debug \
      -DBLAS_LIBRARIES="-L${OPENBLAS_PATH}/lib -lopenblas" \
      -DLAPACK_LIBRARIES="-L${OPENBLAS_PATH}/lib -lopenblas" \
      -DBUILD_SHARED_LIBS=ON && \
    make -j$(nproc) && \
    make install && \
    rm -rf /root/scalapack.tgz /root/scalapack-${SCALAPACK_VERSION}/ && \
    echo ${SCALAPACK_PATH}/lib >> /etc/ld.so.conf.d/scalapack.conf && \
    ldconfig

# Install hwloc
ARG HWLOC_MAJOR=2
ARG HWLOC_MINOR=2
ARG HWLOC_PATCH=0
ENV HWLOC_VERSION=${HWLOC_MAJOR}.${HWLOC_MINOR}.${HWLOC_PATCH}
ARG HWLOC_PATH=/usr/local/hwloc
ENV HWLOC_PATH=${HWLOC_PATH}
RUN wget -q https://download.open-mpi.org/release/hwloc/v${HWLOC_MAJOR}.${HWLOC_MINOR}/hwloc-${HWLOC_VERSION}.tar.gz -O hwloc.tar.gz && \
    tar -xzf hwloc.tar.gz && \
    cd hwloc-${HWLOC_VERSION} && \
    ./configure --prefix=$HWLOC_PATH && \
    make -j$(nproc) install && \
    rm -rf /root/hwloc.tar.gz /root/hwloc-${HWLOC_VERSION}

# TODO: Install ELPA

# TODO: Install PaRSEC + DPLASMA

# Add deployment tooling
RUN wget -q https://github.com/haampie/libtree/releases/download/v1.2.0/libtree_x86_64.tar.gz && \
    tar -xzf libtree_x86_64.tar.gz && \
    rm libtree_x86_64.tar.gz && \
    ln -s /root/libtree/libtree /usr/local/bin/libtree
