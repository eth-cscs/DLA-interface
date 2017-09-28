#ifndef DLA_INTERFACE_TEST_INCLUDE_REF_SCALAPACK_TOOLS_H
#define DLA_INTERFACE_TEST_INCLUDE_REF_SCALAPACK_TOOLS_H

namespace reference_scalapack_tools {

  int numroc(int n, int nb, int iproc, int isrcproc, int nprocs);

  void infog1l(int gindx, int nb, int nprocs, int myroc, int isrcproc, int& lindx, int& rocsrc);

  int indxl2g(int indxloc, int nb, int iproc, int isrcproc, int nprocs);

}  // namespace reference_scalapack

#endif  // DLA_INTERFACE_TEST_INCLUDE_REF_SCALAPACK_TOOLS_H
