#include "ref_scalapack_tools.h"

namespace reference_scalapack_tools {

  int numroc(int n, int nb, int iproc, int isrcproc, int nprocs) {
    // -- ScaLAPACK tools routine (version 1.7) --
    //    University of Tennessee, Knoxville, Oak Ridge National Laboratory,
    //    and University of California, Berkeley.
    //    May 1, 1997
    //
    // Purpose
    // =======
    //
    // NUMROC computes the NUMber of Rows Or Columns of a distributed
    // matrix owned by the process indicated by IPROC.
    //
    // Arguments
    // =========
    //
    // N         (global input) INTEGER
    //           The number of rows/columns in distributed matrix.
    //
    // NB        (global input) INTEGER
    //           Block size, size of the blocks the distributed matrix is
    //           split into.
    //
    // IPROC     (local input) INTEGER
    //           The coordinate of the process whose local array row or
    //           column is to be determined.
    //
    // ISRCPROC  (global input) INTEGER
    //           The coordinate of the process that possesses the first
    //           row or column of the distributed matrix.
    //
    // NPROCS    (global input) INTEGER
    //           The total number processes over which the matrix is
    //           distributed.

    // Figure PROC's distance from source process
    int mydist = (nprocs + iproc - isrcproc) % nprocs;

    // Figure the total number of whole NB blocks N is split up into
    int nblocks = n / nb;

    // Figure the minimum number of rows/cols a process can have
    int return_value = (nblocks / nprocs) * nb;

    // See if there are any extra blocks
    int extrablks = nblocks % nprocs;

    // If I have an extra block
    if (mydist < extrablks)
      return_value += nb;

    // If I have last block, it may be a partial block
    else if (mydist == extrablks)
      return_value += n % nb;

    return return_value;
  }

  void infog1l(int gindx, int nb, int nprocs, int myroc, int isrcproc, int& lindx, int& rocsrc) {
    // -- ScaLAPACK tools routine (version 1.7) --
    //    University of Tennessee, Knoxville, Oak Ridge National Laboratory,
    //    and University of California, Berkeley.
    //    May 1, 1997
    //
    // Purpose
    // =======
    //
    // INFOG1L computes the starting local indexes LINDX corresponding to
    // the distributed submatrix starting globally at the entry pointed by
    // GINDX.  This routine returns the coordinates of the process in the
    // grid owning the submatrix entry of global index GINDX: ROCSRC.
    // INFOG1L is a 1-dimensional version of INFOG2L.
    //
    // Arguments
    // =========
    //
    // GINDX     (global input) INTEGER
    //           The global starting index of the submatrix.
    //
    // NB        (global input) INTEGER
    //           The block size.
    //
    // NPROCS    (global input) INTEGER
    //           The total number of processes over which the distributed
    //           submatrix is distributed.
    //
    // MYROC     (local input) INTEGER
    //           The coordinate of the process calling this routine.
    //
    // ISRCPROC  (global input) INTEGER
    //           The coordinate of the process having the first entry of
    //           the distributed submatrix.
    //
    // LINDX     (local output) INTEGER
    //           The local starting indexes of the distributed submatrix.
    //
    // ROCSRC    (global output) INTEGER
    //           The coordinate of the process that possesses the first
    //           row and column of the submatrix.

    int gcpy = gindx - 1;
    int iblk = gcpy / nb;
    rocsrc = (iblk + isrcproc) % nprocs;

    lindx = (iblk / nprocs + 1) * nb + 1;

    if ((myroc + nprocs - isrcproc) % nprocs >= iblk % nprocs) {
      if (myroc == rocsrc) {
        lindx += gcpy % nb;
      }
      lindx = lindx - nb;
    }
    return;
  }

}  // namespace reference_scalapack
