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

RUN cd /usr/local/bin && \
    curl -Ls https://codecov.io/bash > codecov.sh && \
    echo "d6aa3207c4908d123bd8af62ec0538e3f2b9f257c3de62fad4e29cd3b59b41d9 codecov.sh" | sha256sum --check --quiet && \
    chmod +x codecov.sh

# Install cmake
RUN wget -qO- "https://cmake.org/files/v3.18/cmake-3.18.6-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C /usr/local

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

# Install ELPA
ARG ELPA_VERSION=2020.05.001
ENV ELPA_VERSION=${ELPA_VERSION}
ARG ELPA_PATH=/usr/local/elpa
ENV ELPA_PATH=${ELPA_PATH}
RUN wget -q https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/${ELPA_VERSION}/elpa-${ELPA_VERSION}.tar.gz && \
    tar -xzf elpa-${ELPA_VERSION}.tar.gz && \
    rm elpa-${ELPA_VERSION}.tar.gz && \
    cd elpa-${ELPA_VERSION} && \
    CC=/usr/local/mpich/bin/mpicc FC=/usr/local/mpich/bin/mpif90 ./configure \
      FCFLAGS="-g3 -mavx2 -mavx -mfma" \
      CFLAGS="-g3 -mavx2 -mavx -mfma -funsafe-loop-optimizations -funsafe-math-optimizations -ftree-vect-loop-version -ftree-vectorize" \
      --enable-option-checking=fatal \
      SCALAPACK_LDFLAGS="-L${SCALAPACK_PATH}/lib -lscalapack -L${OPENBLAS_PATH}/lib -lopenblas" \
      MPI_BINARY=/usr/local/mpich/bin/mpirun \
      --enable-single-precision \
      --disable-avx512 \
      --prefix=$ELPA_PATH && \
    make install -j$(nproc) && \
    rm -rf /root/elpa-${ELPA_VERSION} && \
    echo ${ELPA_PATH}/lib >> /etc/ld.so.conf.d/elpa.conf && \
    ldconfig

# Install Boost
ARG BOOST_MAJOR=1
ARG BOOST_MINOR=72
ARG BOOST_PATCH=0
ARG BOOST_PATH=/usr/local/boost
RUN wget -q https://dl.bintray.com/boostorg/release/${BOOST_MAJOR}.${BOOST_MINOR}.${BOOST_PATCH}/source/boost_${BOOST_MAJOR}_${BOOST_MINOR}_${BOOST_PATCH}.tar.gz -O boost.tar.gz && \
    tar -xzf boost.tar.gz && \
    cd boost_${BOOST_MAJOR}_${BOOST_MINOR}_${BOOST_PATCH} && \
    ./bootstrap.sh --prefix=$BOOST_PATH && \
    ./b2 toolset=gcc variant=debug -j$(nproc) install && \
    rm -rf /root/boost.tar.gz /root/boost_${BOOST_MAJOR}_${BOOST_MINOR}_${BOOST_PATCH}

# Install tcmalloc (their version tagging is a bit inconsistent; patch version is not always included)
ARG GPERFTOOLS_VERSION=2.7
ARG GPERFTOOLS_PATH=/usr/local/gperftools
RUN wget -q https://github.com/gperftools/gperftools/releases/download/gperftools-${GPERFTOOLS_VERSION}/gperftools-${GPERFTOOLS_VERSION}.tar.gz -O gperftools.tar.gz && \
    tar -xzf gperftools.tar.gz && \
    cd gperftools-${GPERFTOOLS_VERSION} && \
    ./configure \
      --prefix=${GPERFTOOLS_PATH} && \
    make -j$(nproc) && \
    make install && \
    rm -rf /root/gperftools.tar.gz /root/gperftools-${GPERFTOOLS_VERSION}

# Install HPX
ARG HPX_FORK=STEllAR-GROUP
ARG HPX_VERSION=1.6.0
ARG HPX_WITH_CUDA=OFF
ARG HPX_PATH=/usr/local/hpx
RUN wget -q https://github.com/${HPX_FORK}/hpx/archive/${HPX_VERSION}.tar.gz -O hpx.tar.gz && \
    tar -xzf hpx.tar.gz && \
    cd hpx-${HPX_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. \
      -DBOOST_ROOT=$BOOST_PATH \
      -DHWLOC_ROOT=$HWLOC_PATH \
      -DTCMALLOC_ROOT=$GPERFTOOLS_PATH \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=$HPX_PATH \
      -DCMAKE_CXX_FLAGS_DEBUG="-g -Og -fno-omit-frame-pointer" \
      -DHPX_WITH_SANITIZERS=ON \
      -DHPX_WITH_STACK_OVERFLOW_DETECTION=OFF \
      -DHPX_WITH_MAX_CPU_COUNT=128 \
      -DHPX_WITH_NETWORKING=OFF \
      -DHPX_WITH_CUDA=$HPX_WITH_CUDA \
      -DHPX_WITH_TESTS=OFF \
      -DHPX_WITH_EXAMPLES=OFF && \
    make -j$(nproc) && \
    make install && \
    rm -rf /root/hpx.tar.gz /root/hpx-${HPX_VERSION}

# Install BLASPP
ARG BLASPP_VERSION=2020.10.02
ARG BLASPP_PATH=/usr/local/blaspp
RUN wget -q https://bitbucket.org/icl/blaspp/downloads/blaspp-${BLASPP_VERSION}.tar.gz -O blaspp.tar.gz && \
    tar -xzf blaspp.tar.gz && \
    cd blaspp-${BLASPP_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. \
      -Dbuild_tests=OFF \
      -DBLAS_LIBRARIES="-L${OPENBLAS_PATH}/lib;-lopenblas" \
      -Dblas=OpenBLAS \
      -DCMAKE_BUILD_TYPE=Debug \
      -Duse_openmp=OFF \
      -Duse_cuda=OFF \
      -DCMAKE_CXX_FLAGS_DEBUG="-g -Og -fno-omit-frame-pointer" \
      -DCMAKE_INSTALL_PREFIX=$BLASPP_PATH && \
    make -j$(nproc) && \
    make install && \
    rm -rf /root/blaspp.tar.gz /root/blaspp-${BLASPP_VERSION}

# Install LAPACKPP
ARG LAPACKPP_VERSION=2020.10.02
ARG LAPACKPP_PATH=/usr/local/lapackpp
RUN wget -q https://bitbucket.org/icl/lapackpp/downloads/lapackpp-$LAPACKPP_VERSION.tar.gz -O lapackpp.tar.gz && \
    tar -xzf lapackpp.tar.gz && \
    cd lapackpp-${LAPACKPP_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. \
      -DCMAKE_BUILD_TYPE=Debug \
      -Dbuild_tests=OFF \
      -DCMAKE_CXX_FLAGS_DEBUG="-g -Og -fno-omit-frame-pointer" \
      -DCMAKE_INSTALL_PREFIX=$LAPACKPP_PATH && \
    make -j$(nproc) install && \
    rm -rf /root/lapackpp.tar.gz /root/lapackpp-${LAPACKPP_VERSION}

# Install DLA-Future
ARG DLAF_VERSION=1312dfc
ARG DLAF_PATH=/usr/local/dlaf
ENV DLAF_PATH=${DLAF_PATH}
ARG DLAF_WITH_CUDA=OFF
RUN wget -q https://github.com/eth-cscs/DLA-Future/tarball/${DLAF_VERSION} -O dlaf.tar.gz && \
    tar -xzf dlaf.tar.gz && \
    cd eth-cscs-DLA-Future-${DLAF_VERSION} && \
    mkdir build && cd build && \
    CC=/usr/local/mpich/bin/mpicc CXX=/usr/local/mpich/bin/mpicxx cmake .. \
      -DMKL_ROOT=/opt/intel/compilers_and_libraries/linux/mkl \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=${DLAF_PATH} \
      -DLAPACK_LIBRARY="-L/${OPENBLAS_PATH}/lib;openblas" \
      -DDLAF_WITH_CUDA=${DLAF_WITH_CUDA} \
      -DCMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES=/usr/local/cuda/targets/x86_64-linux/include \
      -DCMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES="/usr/local/cuda/targets/x86_64-linux/lib/stubs;/usr/local/cuda/targets/x86_64-linux/lib;/usr/lib/gcc/x86_64-linux-gnu/7;/usr/lib/x86_64-linux-gnu;/usr/lib;/lib/x86_64-linux-gnu;/lib;/usr/local/cuda/lib64/stubs" \
      -DDLAF_WITH_MKL=OFF \
      -DDLAF_BUILD_TESTING=OFF \
      -DDLAF_BUILD_MINIAPPS=OFF && \
    make -j$(nproc) install && \
    rm -rf /root/dlaf.tar.gz /root/eth-cscs-DLA-Future-${DLAF_VERSION}

# TODO: Install PaRSEC + DPLASMA

# Add deployment tooling
RUN wget -q https://github.com/haampie/libtree/releases/download/v1.2.0/libtree_x86_64.tar.gz && \
    tar -xzf libtree_x86_64.tar.gz && \
    rm libtree_x86_64.tar.gz && \
    ln -s /root/libtree/libtree /usr/local/bin/libtree
