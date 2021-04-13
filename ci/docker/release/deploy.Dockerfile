ARG BUILD_IMAGE

# This is the folder where the project is built
ARG BUILD=/DLA-Interface-build

# This is where we copy the sources to
ARG SOURCE=/DLA-Interface

# Where "make install" should go / we copy the bare minimum
# of binaries to here
ARG DEPLOY=/root/DLA-Interface.bundle

FROM $BUILD_IMAGE as builder

ARG BUILD
ARG SOURCE
ARG DEPLOY

# Build DLA-Interface
COPY . ${SOURCE}

SHELL ["/bin/bash", "-c"]

RUN mkdir ${BUILD} && cd ${BUILD} && \
    source /opt/intel/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64 && \
    CC=${MPICH_PATH}/bin/mpicc CXX=${MPICH_PATH}/bin/mpicxx FC=${MPICH_PATH}/bin/mpif90 cmake ${SOURCE} \
      -DCMAKE_BUILD_TYPE=$build_type \
      -DDLAI_CI_RUNNER_USES_MPIRUN=ON \
      -DDLAI_WITH_MKL=1 \
      -DHWLOC_ROOT=${HWLOC_PATH} \
      -DDLAI_WITH_ELPA=ON \
      -DELPA_MODULE_SPEC=$ELPA_PATH/lib/pkgconfig/elpa-${ELPA_VERSION}.pc \
      -DDLAI_WITH_DLAF=ON \
      -DDLAF_ROOT=$DLAF_PATH \
      -DDLAI_BUILD_TESTING=ON \
      -DDLAI_BUILD_MINIAPPS=ON && \
      make -j$(nproc)

# Prune and bundle binaries
RUN mkdir ${BUILD}-tmp && cd ${BUILD} && \
    source /opt/intel/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64 && \
    export TEST_BINARIES=`ctest --show-only=json-v1 | jq '.tests | map(.command[0]) | .[]' | tr -d \"` && \
    libtree -d ${DEPLOY} ${TEST_BINARIES} && \
    rm -rf ${DEPLOY}/usr/bin && \
    libtree -d ${DEPLOY} $(which ctest addr2line) && \
    cp -L ${SOURCE}/ci/mpi-ctest ${DEPLOY}/usr/bin && \
    echo "$TEST_BINARIES" | xargs -I{file} find -samefile {file} -exec cp --parents '{}' ${BUILD}-tmp ';' && \
    find -name CTestTestfile.cmake -exec cp --parent '{}' ${BUILD}-tmp ';' && \
    rm -rf ${BUILD} && \
    mv ${BUILD}-tmp ${BUILD}

# Deploy MKL separately, since it dlopen's some libs
RUN source /opt/intel/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64 && \
    export MKL_LIB=$MKLROOT/lib/intel64 && \
    libtree -d ${DEPLOY} \
    --chrpath \
    ${MKL_LIB}/libmkl_avx.so \
    ${MKL_LIB}/libmkl_avx2.so \
    ${MKL_LIB}/libmkl_core.so \
    ${MKL_LIB}/libmkl_def.so \
    ${MKL_LIB}/libmkl_intel_thread.so \
    ${MKL_LIB}/libmkl_mc.so \
    ${MKL_LIB}/libmkl_mc3.so \
    ${MKL_LIB}/libmkl_sequential.so \
    ${MKL_LIB}/libmkl_tbb_thread.so \
    ${MKL_LIB}/libmkl_vml_avx.so \
    ${MKL_LIB}/libmkl_vml_avx2.so \
    ${MKL_LIB}/libmkl_vml_cmpt.so \
    ${MKL_LIB}/libmkl_vml_def.so \
    ${MKL_LIB}/libmkl_vml_mc.so \
    ${MKL_LIB}/libmkl_vml_mc3.so

FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

ARG BUILD
ARG DEPLOY

# tzdata is needed to print correct time
RUN apt-get update -qq && \
    apt-get install -qq -y --no-install-recommends \
      tzdata && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder ${BUILD} ${BUILD}
COPY --from=builder ${DEPLOY} ${DEPLOY}

# Make it easy to call our binaries.
ENV PATH="${DEPLOY}/usr/bin:$PATH"

# Automatically print stacktraces on segfault
ENV LD_PRELOAD=/lib/x86_64-linux-gnu/libSegFault.so

RUN echo "${DEPLOY}/usr/lib/" > /etc/ld.so.conf.d/dlai.conf && ldconfig

WORKDIR ${BUILD}
