module switch PrgEnv-cray PrgEnv-gnu
module load daint-mc        
module load CMake
module load cray-python
module unload cray-libsci   # should have no effect but it makes sure that libsci is not used if MKL is linked incorrectly.
module load intel           # defines $MKLROOT
module load hwloc
module load Boost

export CRAYPE_LINK_TYPE=dynamic

# Use MKL 2018u2 (Fix LU bug)
source /apps/daint/UES/sandbox/rasolca/intel/mkl/bin/mklvars.sh intel64
