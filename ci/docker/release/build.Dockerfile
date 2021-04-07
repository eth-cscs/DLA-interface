ARG BASE_IMAGE=ubuntu:20.04

FROM $BASE_IMAGE

WORKDIR /root

SHELL ["/bin/bash", "-l", "-c"]

ENV DEBIAN_FRONTEND noninteractive
ENV FORCE_UNSAFE_CONFIGURE 1

# Install basic tools
RUN apt-get update -qq && apt-get install -qq -y --no-install-recommends \
    software-properties-common \
    build-essential gfortran \
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

# Install MKL
ARG MKL_VERSION=2020.0-088
RUN wget -qO - https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB 2>/dev/null | apt-key add - && \
    apt-add-repository 'deb https://apt.repos.intel.com/mkl all main' && \
    apt-get install -y -qq --no-install-recommends intel-mkl-64bit-${MKL_VERSION} && \
    rm -rf /var/lib/apt/lists/* && \
    echo "/opt/intel/lib/intel64\n/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64" >> /etc/ld.so.conf.d/intel.conf && \
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
#ENV ELPA_VERSION=2020.05.001
#ENV ELPA_PATH=/usr/local/elpa
#RUN wget -q https://elpa.mpcdf.mpg.de/html/Releases/${ELPA_VERSION}/elpa-${ELPA_VERSION}.tar.gz && \
#    tar -xzf elpa-${ELPA_VERSION}.tar.gz && \
#    rm elpa-${ELPA_VERSION}.tar.gz && \
#    cd elpa-${ELPA_VERSION} && \
#    source /opt/intel/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64 && \
#    CC=/usr/local/mpich/bin/mpicc FC=/usr/local/mpich/bin/mpif90 ./configure \
#      FCFLAGS="-O3 -march=native -mavx512f -mavx2 -mavx -mfma" \
#      CFLAGS="-O3 -march=native -mavx512f -mavx2 -mavx -mfma -funsafe-loop-optimizations -funsafe-math-optimizations -ftree-vect-loop-version -ftree-vectorize" \
#      --enable-option-checking=fatal \
#      SCALAPACK_LDFLAGS="-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl" \
#      MPI_BINARY=/usr/local/mpich/bin/mpirun \
#      --enable-single-precision \
#      --prefix=$ELPA_PATH && \
#    make install -j$(nproc) && \
#    rm -rf /root/elpa-${ELPA_VERSION} && \
#    echo ${ELPA_PATH}/lib >> /etc/ld.so.conf.d/elpa.conf && \
#    ldconfig

# Install PaRSEC
#ARG PaRSEC_VERSION=bccb6e190553
#ENV PaRSEC_PATH=/usr/local/parsec
#RUN wget -q https://bitbucket.org/icldistcomp/parsec/get/${PaRSEC_VERSION}.tar.gz -O parsec.tar.gz && \
#    tar -xzf parsec.tar.gz && \
#    rm parsec.tar.gz && \
#    cd icldistcomp-parsec-${PaRSEC_VERSION} && \
#    mkdir build && \
#    cd build && \
#    cmake .. \
#      -DMPI_ROOT=${MPICH_PATH} \
#      -DCMAKE_INSTALL_PREFIX=${PaRSEC_PATH} && \
#    make install -j$(nproc) && \
#    rm -rf icldistcomp-parsec-${PaRSEC_VERSION} && \
#    echo ${PaRSEC_PATH}/lib >> /etc/ld.so.conf.d/parsec.conf && \
#    ldconfig

# Install DPLASMA
#ARG DPLASMA_VERSION=51858e3183ea
#ENV DPLASMA_PATH=/usr/local/dplasma
#RUN wget -q https://bitbucket.org/icldistcomp/dplasma/get/${DPLASMA_VERSION}.tar.gz -O dplasma.tar.gz && \
#    tar -xzf dplasma.tar.gz && \
#    rm dplasma.tar.gz && \
#    cd icldistcomp-dplasma-${DPLASMA_VERSION} && \
#    source /opt/intel/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64 && \
#    mkdir build && \
#    cd build && \
#    cmake .. \
#      -DPython_EXECUTABLE=`which python` \
#      -DPaRSEC_ROOT=${PaRSEC_PATH} \
#      -DMPI_ROOT=${MPICH_PATH} \
#      -DCMAKE_INSTALL_PREFIX=${DPLASMA_PATH} && \
#    make install -j$(nproc) && \
#    rm -rf icldistcomp-dplasma-${DPLASMA_VERSION} && \
#    echo ${DPLASMA_PATH}/lib >> /etc/ld.so.conf.d/dplasma.conf && \
#    ldconfig

# Add deployment tooling
RUN wget -q https://github.com/haampie/libtree/releases/download/v1.2.0/libtree_x86_64.tar.gz && \
    tar -xzf libtree_x86_64.tar.gz && \
    rm libtree_x86_64.tar.gz && \
    ln -s /root/libtree/libtree /usr/local/bin/libtree
