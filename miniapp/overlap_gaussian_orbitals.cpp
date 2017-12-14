#include <array>
#include <iostream>
#include <random>
#include <vector>
#include <boost/program_options.hpp>
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
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("size,n", boost::program_options::value<int>()->default_value(1024), "The size of the square matrix. (Default = 1024)")
      ("nb", boost::program_options::value<int>()->default_value(128), "The block size for the block cyclic distribution. (Default = 128)")
      ("alpha", boost::program_options::value<double>()->default_value(.5), "(Default = .5)")
      ("a", boost::program_options::value<double>()->default_value(20), "(Default = 20.)")
      ("seed", boost::program_options::value<unsigned int>()->default_value(2017), "Seed for random number generator. (Default = 2017)")
      ("upper,U", "Use upper triangular part of the matrix. (If not specified the lower triangular part is used.)")
      ("scalapack", "Run test with ScaLAPACK.")
      ("dplasma", "Run test with DPlasma.")
      ("row_procs,p", boost::program_options::value<int>()->default_value(1), "The number of rows in the 2D communicator. (Default = 1)")
      ("col_procs,q", boost::program_options::value<int>()->default_value(1), "The number of cols in the 2D communicator. (Default = 1)")
      ("nr_threads", boost::program_options::value<int>()->default_value(1), "The number of threads per rank. (Default 1)")
      ("no_check", "Disable correctness checking.")
      ;
  // clang-format on

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
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

  comm::CommunicatorManager::initialize(nr_threads, &argc, &argv, true);

  auto& comm_grid =
      comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, p, q, RowMajor);

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

    choleskyFactorization(uplo, mat, solver, 2);

    if (check) {
      if (uplo == Lower) {
        for (int j = 0; j < mat.localSize().second; ++j) {
          for (int i = 0; i < mat.localSize().first; ++i) {
            Local2DIndex local_index(i, j);
            Global2DIndex global_index = mat.getGlobal2DIndex(local_index);
            if (global_index.row < global_index.col)
              mat(local_index) = 0;
          }
        }

        matrixMultiply(NoTrans, Trans, -1., mat, mat, 1., *mat_copy, solver, 2);
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

        matrixMultiply(Trans, NoTrans, -1., mat, mat, 1., *mat_copy, solver, 2);
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
  return 0;
}
