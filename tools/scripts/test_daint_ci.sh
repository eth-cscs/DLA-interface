OPTIONS=dfo:v
LONGOPTIONS=debug,force,output:,verbose

OPTIONS=`getopt -o h::t:b:p:l:s:e:d:x: -l help::,tag:,build-type,partition:,lapack:,scalapack:,elpa:,dplasma:,hpx_linalg: --name "$0" -- "$@"`
ret=$?
if [ $ret -ne 0 ]; then
    exit 1
fi

# default values
tag="notag"
build_type="Release"
partition="mc"
lapack="MKLst"
scalapack="MKL"
elpa="Yes"
dplasma="Yes"
hpx_linalg="No"

# Print help statement and exit.
print_help()
{
  printf "Usage: $0 [options]\n\n"
  printf "Options:\n"

  # --tag
  printf "  %-35s %s\n" \
    "-t, --tag" \
    "Build tag [default: ${tag}]."

  # --build-type
  printf "  %-35s %s\n" \
    "-b, --build-type [Release|Debug]" \
    "Build type [default: ${build_type}]."

  # --partition
  printf "  %-35s %s\n" \
    "-p, --partition [mc|gpu]" \
    "Build type [default: ${partition}]."

  # --lapack
  printf "  %-35s %s\n" \
    "-l, --lapack [MKLst|MKLmt]" \
    "Lapack type [default: ${lapack}]."

  # --scalapack
  printf "  %-35s %s\n" \
    "-s, --scalapack [MKL]" \
    "Scalapack type [default: ${scalapack}]."

  # --elpa
  printf "  %-35s %s\n" \
    "-e, --elpa [No|Yes]" \
    "Build with ELPA [default: ${elpa}]."

  # --dplasma
  printf "  %-35s %s\n" \
    "-d, --dplasma [No|Yes]" \
    "Build with DPlasma [default: ${dplasma}]."

  # --hpx_linalg
  printf "  %-35s %s\n" \
         "-x, --hpx_linalg [No|Yes]" \
         "Build with HPX_LINALG [default: ${hpx_linalg}]."

  # --help
  printf "  %-35s %s\n" \
    "-h, --help" \
    "Print this help statement."
}

eval set -- "$OPTIONS"

# now enjoy the options in order and nicely split until we see --
while true; do
  case "$1" in
    -h|--h*)         print_help ; exit 0 ;;
    -b|--build-type) build_type="$2" ; shift 2 ;;
    -t|--tag)        tag="$2"        ; shift 2 ;;
    -p|--partition)  partition="$2"  ; shift 2 ;;
    -l|--lapack)     lapack="$2"     ; shift 2 ;;
    -s|--scalapack)  scalapack="$2"  ; shift 2 ;;
    -e|--elpa)       elpa="$2"       ; shift 2 ;;
    -d|--dplasma)    dplasma="$2"    ; shift 2 ;;
    -x|--hpx_linalg) hpx_linalg="$2" ; shift 2 ;;
    --) shift ; break ;;
    *) echo "Options internal error. $1" ; exit 1 ;;
  esac
done

case $partition in
  mc|gpu) ;;
  *) echo "Wrong --partition option: $partition" ; print_help ; exit 1 ;;
esac

case $lapack in
  MKLst) OPT_LAPACK=(-DDLA_LAPACK_TYPE=MKL -DMKL_THREADING=Sequential) ;;
  MKLmt) OPT_LAPACK=(-DDLA_LAPACK_TYPE=MKL -DMKL_THREADING="Intel OpenMP") ;;
  *) echo "Wrong --lapack option: $lapack" ; print_help ; exit 1 ;;
esac

case $scalapack in
  MKL)   OPT_SCALAPACK=(-DDLA_SCALAPACK_TYPE=MKL) ;;
  *) echo "Wrong --scalapack option: $scalapack" ; print_help ; exit 1 ;;
esac

ELPA_VERS=2018.05.001
ELPA_LIB_DIR=/apps/daint/UES/sandbox/rasolca/elpa

case $elpa in
  No|no)   OPT_ELPA=() ;;
  Yes|yes) OPT_ELPA=(-DELPA_ROOT="${ELPA_LIB_DIR}" -DELPA_VERSION="${ELPA_VERS}") ;;
  *)  echo "Wrong --elpa option: $elpa" ; print_help ; exit 1 ;;
esac

PLASMA_LIB_DIR=/apps/daint/UES/sandbox/rasolca/plasma-${partition}
PSEC_LIB_DIR=/apps/daint/UES/sandbox/rasolca/parsec-${partition}

case $dplasma in
  No|no)   OPT_DPLASMA=() ;;
  Yes|yes) OPT_DPLASMA=(-DParsec_DIR="${PSEC_LIB_DIR}/lib/cmake/" -DPLASMA_DIR="${PLASMA_LIB_DIR}") ;;
  *)  echo "Wrong --dplasma option: $dplasma" ; print_help ; exit 1 ;;
esac

HPX_LINALG_DIR=/apps/daint/UES/simbergm/jenkins/DLA-interface/hpx_linalg/lib/cmake/HPX_LINALG/

case $hpx_linalg in
    No|no)   OPT_HPX_LINALG=() ;;
    Yes|yes) OPT_HPX_LINALG=(-DHPX_LINALG_DIR="${HPX_LINALG_DIR}") ;;
    *)  echo "Wrong --hpx_linalg option: $hpx_linalg" ; print_help ; exit 1 ;;
esac

echo -n "Running with options:"
echo -n " tag: ${tag},"
echo -n " build_type: ${build_type},"
echo -n " partition: ${partition},"
echo -n " lapack: ${lapack},"
echo -n " scalapack: ${scalapack},"
echo -n " elpa: ${elpa},"
echo -n " dplasma: ${dplasma}"
echo    " hpx_linalg: ${hpx_linalg}"

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
SRC_DIR=$SCRIPT_DIR/../..

source $SCRIPT_DIR/daint-${partition}_env.sh

cd $SRC_DIR

BUILD_DIR=$SRC_DIR/build_$tag

rm -rf $BUILD_DIR
mkdir $BUILD_DIR
cd $BUILD_DIR

OPT_CMAKE=(\
  -DCMAKE_BUILD_TYPE=$build_type \
  -DTEST_RUNNER="srun" \
  -DDLA_ALL_TESTS_USE_RUNNER=ON \
  "${OPT_LAPACK[@]}" \
  "${OPT_SCALAPACK[@]}" \
  "${OPT_ELPA[@]}" \
  "${OPT_DPLASMA[@]}" \
  "${OPT_HPX_LINALG[@]}" \
  -DHWLOC_ROOT=$EBROOTHWLOC \
  ..)
echo -n "executing cmake with options: ${OPT_CMAKE[@]}"
CC=cc CXX=CC FC=ftn cmake "${OPT_CMAKE[@]}"
ret=$?

if [ $ret -ne 0 ]; then
  echo "Configuration failed."
  exit $ret
fi

make -j 8
ret=$?

if [ $ret -ne 0 ]; then
  echo "Building failed."
  exit $ret
fi

make test
