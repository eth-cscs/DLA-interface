#include <array>
#include <iostream>
#include <random>
#include <tuple>
#include <vector>
#include "cxxopts/cxxopts.hpp"
#include "c_dla_interface.h"
#include "distributed_matrix.h"
#include "dla_interface.h"

inline double square(double x) {
  return x * x;
}

double overlapElement(double alpha, const std::array<double, 3>& p1, const std::array<double, 3>& p2) {
  double tmp = 0;
  for (auto i : {0, 1, 2}) {
    tmp += square(p1[i] - p2[i]);
  }
  return std::exp(-tmp * alpha / 2.);
}

using namespace dla_interface;

int main(int argc, char** argv) {
  // clang-format off
  cxxopts::Options desc(argv[0], "Allowed Options");
  desc.add_options()
      ("help", "produce help message")
      ("n,size", "The size of the square matrix.", cxxopts::value<int>()->default_value("1024"))
      ("nb", "The block size for the block cyclic distribution. ", cxxopts::value<int>()->default_value("128"))
      ("alpha", "", cxxopts::value<double>()->default_value(".5"))
      ("a", "", cxxopts::value<double>()->default_value("20"))
      ("seed", "Seed for random number generator.", cxxopts::value<unsigned int>()->default_value("2017"))
      ("U,upper", "Use upper triangular part of the matrix. (If not specified the lower triangular part is used.)")
      ("scalapack", "Run test with ScaLAPACK.")
      ("dplasma", "Run test with DPlasma.")
      ("p,row_procs", "The number of rows in the 2D communicator.", cxxopts::value<int>()->default_value("1"))
      ("q,col_procs", "The number of cols in the 2D communicator.", cxxopts::value<int>()->default_value("1"))
      ("nr_threads", "The number of threads per rank.", cxxopts::value<int>()->default_value("1"))
      ("no_check", "Disable correctness checking.")
      ;
  // clang-format on

  auto vm = desc.parse(argc, argv);

  if (vm.count("help")) {
    std::cout << desc.help() << std::endl;
    return 1;
  }

  const int n = vm["size"].as<int>();
  const int nb = vm["nb"].as<int>();
  const double alpha = vm["alpha"].as<double>();
  const double a = vm["a"].as<double>();
  const unsigned int seed = vm["seed"].as<unsigned int>();
  const int p = vm["row_procs"].as<int>();
  const int q = vm["col_procs"].as<int>();
  const int nr_threads = vm["nr_threads"].as<int>();

  const bool check = !vm.count("no_check");
  const UpLo uplo = vm.count("upper") ? Upper : Lower;

  std::vector<SolverType> solvers;
  if (vm.count("scalapack"))
    solvers.push_back(ScaLAPACK);
  if (vm.count("dplasma"))
    solvers.push_back(DPlasma);

  int init_mpi = 1;
  MPI_Fint mpi_comm_f = MPI_Comm_c2f(MPI_COMM_WORLD);

  //////////////////////////////////////////////////////////////
  // C interface
  //////////////////////////////////////////////////////////////
  dlai_initialize_arg(&nr_threads, &argc, &argv, &init_mpi);
  int ictxt = dlai_create_2d_grid(&mpi_comm_f, &p, &q, "R");

  int timer_opt = 2;
  dlai_set_print_timer_option(&timer_opt);
  //////////////////////////////////////////////////////////////

  auto& comm_grid = comm::CommunicatorManager::getCommunicator2DGridFromBlacsContext(ictxt);

  std::vector<std::array<double, 3>> point;
  point.reserve(n);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> dist(-a, a);
  for (int i = 0; i < n; ++i) {
    point.emplace_back(std::array<double, 3>({dist(rng), dist(rng), dist(rng)}));
  }

  DistributedMatrix<double> mat(n, n, nb, nb, comm_grid, scalapack_dist);
  std::unique_ptr<DistributedMatrix<double>> mat_copy;

  for (auto solver : solvers) {
    // Set the elements.
    for (int j = 0; j < mat.localSize().second; ++j) {
      for (int i = 0; i < mat.localSize().first; ++i) {
        Local2DIndex local_index(i, j);
        Global2DIndex global_index = mat.getGlobal2DIndex(local_index);
        mat(local_index) = overlapElement(alpha, point[global_index.row], point[global_index.col]);
      }
    }
    if (check) {
      // create a copy for checking.
      mat_copy = std::make_unique<DistributedMatrix<double>>(mat);
    }

    auto mat_sca = mat.getScalapackDescription();
    auto desca_array = std::get<3>(mat_sca);

    std::string solver_str = util::getSolverString(solver);
    const char uplo_c = static_cast<char>(uplo);
    double* a_ptr = std::get<0>(mat_sca);
    int ia = std::get<1>(mat_sca);
    int ja = std::get<2>(mat_sca);
    int* desca = &desca_array[0];
    const char* solver_c = solver_str.c_str();

    //////////////////////////////////////////////////////////////
    // C interface
    //////////////////////////////////////////////////////////////
    dlai_d_cholesky_factorization(&uplo_c, &n, a_ptr, &ia, &ja, desca, solver_c);
    //////////////////////////////////////////////////////////////

    if (check) {
      auto mat_copy_sca = mat_copy->getScalapackDescription();
      auto descc_array = std::get<3>(mat_copy_sca);
      double* c_ptr = std::get<0>(mat_copy_sca);
      int ic = std::get<1>(mat_copy_sca);
      int jc = std::get<2>(mat_copy_sca);
      int* descc = &descc_array[0];

      if (uplo == Lower) {
        for (int j = 0; j < mat.localSize().second; ++j) {
          for (int i = 0; i < mat.localSize().first; ++i) {
            Local2DIndex local_index(i, j);
            Global2DIndex global_index = mat.getGlobal2DIndex(local_index);
            if (global_index.row < global_index.col)
              mat(local_index) = 0;
          }
        }

        //////////////////////////////////////////////////////////////
        // C interface
        //////////////////////////////////////////////////////////////
        dlai_d_matrix_multiplication("N", "T", &n, &n, &n, DLA_D_M_ONE, a_ptr, &ia, &ja, desca, a_ptr,
                                     &ia, &ja, desca, DLA_D_ONE, c_ptr, &ic, &jc, descc, solver_c);
        //////////////////////////////////////////////////////////////
      }
      else {
        for (int j = 0; j < mat.localSize().second; ++j) {
          for (int i = 0; i < mat.localSize().first; ++i) {
            Local2DIndex local_index(i, j);
            Global2DIndex global_index = mat.getGlobal2DIndex(local_index);
            if (global_index.row > global_index.col)
              mat(local_index) = 0;
          }
        }

        //////////////////////////////////////////////////////////////
        // C interface
        //////////////////////////////////////////////////////////////
        dlai_d_matrix_multiplication("T", "N", &n, &n, &n, DLA_D_M_ONE, a_ptr, &ia, &ja, desca, a_ptr,
                                     &ia, &ja, desca, DLA_D_ONE, c_ptr, &ic, &jc, descc, solver_c);
        //////////////////////////////////////////////////////////////
      }

      for (int j = 0; j < mat_copy->localSize().second; ++j) {
        for (int i = 0; i < mat_copy->localSize().first; ++i) {
          Local2DIndex local_index(i, j);
          if (std::abs((*mat_copy)(local_index)) > 1000 * n * std::numeric_limits<double>::epsilon()) {
            std::cout << comm_grid.id2D() << " Correctness test failed! "
                      << mat.getGlobal2DIndex(local_index) << std::endl;
            return 1;
          }
          if (std::abs((*mat_copy)(local_index)) > 100 * n * std::numeric_limits<double>::epsilon()) {
            std::cout << comm_grid.id2D() << " Warning: Correctnes test suspicious! "
                      << mat.getGlobal2DIndex(local_index) << std::endl;
            return 1;
          }
        }
      }
    }
  }

  dlai_finalize();

  return 0;
}
