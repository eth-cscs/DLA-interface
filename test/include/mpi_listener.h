#ifndef DLA_INTERFACE_TEST_INCLUDE_MPI_LISTENER_H
#define DLA_INTERFACE_TEST_INCLUDE_MPI_LISTENER_H

#include <fstream>

namespace testing {
  std::ofstream& setMPIListener(std::string f_base);
  // Implemented in pretty_mpi_unit_test_result_printer.ipp
}

#endif  // DLA_INTERFACE_TEST_INCLUDE_MPI_LISTENER_H
