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
ARG ELPA_VERSION=2020.05.001
ENV ELPA_VERSION=${ELPA_VERSION}
ARG ELPA_PATH=/usr/local/elpa
ENV ELPA_PATH=${ELPA_PATH}
RUN wget -q https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/${ELPA_VERSION}/elpa-${ELPA_VERSION}.tar.gz && \
    tar -xzf elpa-${ELPA_VERSION}.tar.gz && \
    rm elpa-${ELPA_VERSION}.tar.gz && \
    cd elpa-${ELPA_VERSION} && \
    source /opt/intel/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64 && \
    CC=/usr/local/mpich/bin/mpicc FC=/usr/local/mpich/bin/mpif90 ./configure \
      FCFLAGS="-O3 -march=native -mavx512f -mavx2 -mavx -mfma" \
      CFLAGS="-O3 -march=native -mavx512f -mavx2 -mavx -mfma -funsafe-loop-optimizations -funsafe-math-optimizations -ftree-vect-loop-version -ftree-vectorize" \
      --enable-option-checking=fatal \
      SCALAPACK_LDFLAGS="-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl" \
      MPI_BINARY=/usr/local/mpich/bin/mpirun \
      --enable-single-precision \
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
    ./b2 -j$(nproc) debug-symbols=on install && \
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
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_INSTALL_PREFIX=$HPX_PATH \
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
RUN source /opt/intel/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64 && \
    wget -q https://bitbucket.org/icl/blaspp/downloads/blaspp-${BLASPP_VERSION}.tar.gz -O blaspp.tar.gz && \
    tar -xzf blaspp.tar.gz && \
    cd blaspp-${BLASPP_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. \
      -Dbuild_tests=OFF \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -Duse_openmp=OFF \
      -Duse_cuda=OFF \
      -Dblas='Intel MKL' \
      -Dblas_threaded=OFF \
      -DCMAKE_INSTALL_PREFIX=$BLASPP_PATH && \
    make -j$(nproc) && \
    make install && \
    rm -rf /root/blaspp.tar.gz /root/blaspp-${BLASPP_VERSION}

# Install LAPACKPP
ARG LAPACKPP_VERSION=2020.10.02
ARG LAPACKPP_PATH=/usr/local/lapackpp
RUN source /opt/intel/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64 && \
    wget -q https://bitbucket.org/icl/lapackpp/downloads/lapackpp-$LAPACKPP_VERSION.tar.gz -O lapackpp.tar.gz && \
    tar -xzf lapackpp.tar.gz && \
    cd lapackpp-${LAPACKPP_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -Dbuild_tests=OFF \
      -DCMAKE_INSTALL_PREFIX=$LAPACKPP_PATH && \
    make -j$(nproc) install && \
    rm -rf /root/lapackpp.tar.gz /root/lapackpp-${LAPACKPP_VERSION}

# Install DLA-Future
ARG DLAF_VERSION=19b83c0
ARG DLAF_PATH=/usr/local/dlaf
ENV DLAF_PATH=${DLAF_PATH}
ARG DLAF_WITH_CUDA=OFF
RUN source /opt/intel/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64 && \
    wget -q https://github.com/eth-cscs/DLA-Future/tarball/${DLAF_VERSION} -O dlaf.tar.gz && \
    tar -xzf dlaf.tar.gz && \
    cd eth-cscs-DLA-Future-${DLAF_VERSION} && \
    mkdir build && cd build && \
    CC=/usr/local/mpich/bin/mpicc CXX=/usr/local/mpich/bin/mpicxx cmake .. \
      -DMKL_ROOT=/opt/intel/compilers_and_libraries/linux/mkl \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_INSTALL_PREFIX=${DLAF_PATH} \
      -DDLAF_WITH_CUDA=${DLAF_WITH_CUDA} \
      -DCMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES=/usr/local/cuda/targets/x86_64-linux/include \
      -DCMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES="/usr/local/cuda/targets/x86_64-linux/lib/stubs;/usr/local/cuda/targets/x86_64-linux/lib;/usr/lib/gcc/x86_64-linux-gnu/7;/usr/lib/x86_64-linux-gnu;/usr/lib;/lib/x86_64-linux-gnu;/lib;/usr/local/cuda/lib64/stubs" \
      -DDLAF_WITH_MKL=ON \
      -DDLAF_BUILD_TESTING=OFF \
      -DDLAF_BUILD_MINIAPPS=OFF && \
    make -j$(nproc) install && \
    rm -rf /root/dlaf.tar.gz /root/eth-cscs-DLA-Future-${DLAF_VERSION}

# Install PaRSEC
#ARG PARSEC_VERSION=aa7947ff3cd4
#ARG PARSEC_PATH=/usr/local/parsec
#ENV PARSEC_PATH=${PARSEC_PATH}
#ARG PARSEC_WITH_CUDA=OFF # TODO: ADD TO BUILD_ARGS
#RUN wget -q https://bitbucket.org/icldistcomp/parsec/get/${PARSEC_VERSION}.tar.gz -O parsec.tar.gz && \
#    tar -xzf parsec.tar.gz && \
#    rm parsec.tar.gz && \
#    cd icldistcomp-parsec-${PARSEC_VERSION} && \
#    mkdir build && \
#    cd build && \
#    cmake .. \
#      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
#      -DMPI_ROOT=${MPICH_PATH} \
#      -DPARSEC_GPU_WITH_CUDA=${PARSEC_WITH_CUDA} \
#      -DCMAKE_INSTALL_PREFIX=${PARSEC_PATH} && \
#    make install -j$(nproc) && \
#    rm -rf icldistcomp-parsec-${PARSEC_VERSION} && \
#    echo ${PARSEC_PATH}/lib >> /etc/ld.so.conf.d/parsec.conf && \
#    ldconfig

# Install DPLASMA
#ARG DPLASMA_VERSION=1ddaa8eb3cc5
#ARG DPLASMA_PATH=/usr/local/dplasma
#ENV DPLASMA_PATH=${DPLASMA_PATH}
#RUN wget -q https://bitbucket.org/icldistcomp/dplasma/get/${DPLASMA_VERSION}.tar.gz -O dplasma.tar.gz && \
#    tar -xzf dplasma.tar.gz && \
#    rm dplasma.tar.gz && \
#    cd icldistcomp-dplasma-${DPLASMA_VERSION} && \
#    source /opt/intel/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64 && \
#    mkdir build && \
#    cd build && \
#    cmake .. \
#      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
#      -DPython_EXECUTABLE=`which python` \
#      -DPaRSEC_ROOT=${PARSEC_PATH} \
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
