#include <array>
#include <iostream>
#include <random>
#include <vector>
#include "cxxopts/cxxopts.hpp"
#include "distributed_matrix.h"
#include "dla_interface.h"
#include "matrix_index.h"
#include "types.h"
#include "util_types.h"

extern "C" void dgemm_(const char* transa, const char* transb, const int* m, const int* n,
                       const int* k, const double* alpha, const double* a, const int* lda,
                       const double* b, const int* ldb, const double* beta, double* c,
                       const int* ldc);

using namespace dla_interface;

template <class T, class F_T>
void fill_random(DistributedMatrix<T>& mat, F_T random) {
  int m = mat.localSize().first;
  int n = mat.localSize().second;

  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < m; ++i) {
      mat(Local2DIndex(i, j)) = random();
    }
  }
}

int main(int argc, char** argv) {
  // clang-format off
  cxxopts::Options desc(argv[0], "Allowed Options");
  desc.add_options()
      ("help", "produce help message")
      ("m,size_m", "Size m", cxxopts::value<int>()->default_value("1024"))
      ("n,size_n", "Size n", cxxopts::value<int>()->default_value("1024"))
      ("k,size_k", "Size k", cxxopts::value<int>()->default_value("1024"))
      ("def_nb", "The default block size for the block cyclic distribution.", cxxopts::value<int>()->default_value("128"))
      ("a_mb", "The row block size of a for the block cyclic distribution. If <= 0 def_nb is used.", cxxopts::value<int>()->default_value("-1"))
      ("a_nb", "The column block size of a for the block cyclic distribution. If <= 0 def_nb is used.", cxxopts::value<int>()->default_value("-1"))
      ("b_mb", "The row block size of b for the block cyclic distribution. If <= 0 def_nb is used.", cxxopts::value<int>()->default_value("-1"))
      ("b_nb", "The column block size of b for the block cyclic distribution. If <= 0 def_nb is used.", cxxopts::value<int>()->default_value("-1"))
      ("c_mb", "The row block size of c for the block cyclic distribution. If <= 0 def_nb is used.", cxxopts::value<int>()->default_value("-1"))
      ("c_nb", "The column block size of c for the block cyclic distribution. If <= 0 def_nb is used.", cxxopts::value<int>()->default_value("-1"))
      ("seed", "Seed for random number generator.", cxxopts::value<unsigned int>()->default_value("2017"))
      ("transa", "a is stored transposed.")
      ("transb", "b is stored transposed.")
      ("scalapack", "Run test with ScaLAPACK.")
      ("dplasma", "Run test with DPlasma.")
      ("p,row_procs", "The number of rows in the 2D communicator.", cxxopts::value<int>()->default_value("1"))
      ("q,col_procs", "The number of cols in the 2D communicator.", cxxopts::value<int>()->default_value("1"))
      ("nr_threads", "The number of threads per rank.", cxxopts::value<int>()->default_value("1"))
      ("r,repetitions", "The number of time to repeat the multiplication.", cxxopts::value<int>()->default_value("1"))
      ;
  // clang-format on

  auto vm = desc.parse(argc, argv);

  if (vm.count("help")) {
    std::cout << desc.help() << std::endl;
    return 1;
  }

  const int m = vm["m"].as<int>();
  const int n = vm["n"].as<int>();
  const int k = vm["k"].as<int>();
  const int def_nb = vm["def_nb"].as<int>();
  auto bsize = [def_nb](int val) { return val > 0 ? val : def_nb; };
  const int a_mb = bsize(vm["a_mb"].as<int>());
  const int a_nb = bsize(vm["a_nb"].as<int>());
  const int b_mb = bsize(vm["b_mb"].as<int>());
  const int b_nb = bsize(vm["b_nb"].as<int>());
  const int c_mb = bsize(vm["c_mb"].as<int>());
  const int c_nb = bsize(vm["c_nb"].as<int>());
  const unsigned int seed = vm["seed"].as<unsigned int>();
  const int p = vm["row_procs"].as<int>();
  const int q = vm["col_procs"].as<int>();
  const int nr_threads = vm["nr_threads"].as<int>();
  const int rep = std::max(1, vm["repetitions"].as<int>());

  const OpTrans transa = vm.count("transa") ? Trans : NoTrans;
  const OpTrans transb = vm.count("transb") ? Trans : NoTrans;

  const int a_m = transa == NoTrans ? m : k;
  const int a_n = transa == NoTrans ? k : m;
  const int b_m = transb == NoTrans ? k : n;
  const int b_n = transb == NoTrans ? n : k;

  std::vector<SolverType> solvers;
  for (size_t i = 0; i < vm.count("scalapack"); ++i)
    solvers.push_back(ScaLAPACK);
  for (size_t i = 0; i < vm.count("dplasma"); ++i)
    solvers.push_back(DPlasma);

  comm::CommunicatorManager::initialize(nr_threads, &argc, &argv, true);

  auto& comm_grid =
      comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, p, q, RowMajor);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> dist(-1, 1);
  auto random = [&dist, &rng]() { return dist(rng); };

  double alpha = dist(rng);
  double beta = dist(rng);

  DistributedMatrix<double> mat_a(a_m, a_n, a_mb, a_nb, comm_grid, scalapack_dist);
  DistributedMatrix<double> mat_b(b_m, b_n, b_mb, b_nb, comm_grid, scalapack_dist);

  fill_random(mat_a, random);
  fill_random(mat_b, random);

  for (auto solver : solvers) {
    double min_elapsed = 3e9;  // ~100 years

    for (int id_rep = 0; id_rep < rep; ++id_rep) {
      DistributedMatrix<double> mat_c(m, n, c_mb, c_nb, comm_grid, scalapack_dist);

      fill_random(mat_c, random);

      util::Timer<> timer(comm_grid.rowOrderedMPICommunicator());

      matrixMultiplication(transa, transb, alpha, mat_a, mat_b, beta, mat_c, solver, 2);
      auto elapsed = timer.elapsed(0, timer.save_time());
      min_elapsed = std::min(min_elapsed, elapsed);
    }
    if (comm_grid.id2D() == std::make_pair(0, 0)) {
      double mnk = static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(k);
      std::cout << util::getSolverString(solver) << ": Best " << min_elapsed << " s ,"
                << util::nrOps<double>(mnk, mnk) / min_elapsed / 1e9 << " GFlop/s" << std::endl;
    }
  }

  return 0;
}
