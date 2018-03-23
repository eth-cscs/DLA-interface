#include <array>
#include <iostream>
#include <random>
#include <vector>
#include "cxxopts/cxxopts.hpp"
#include "distributed_matrix.h"
#include "dla_interface.h"
#include "types.h"
#include "util_local_matrix.h"
#include "util_distributed_matrix.h"

extern "C" void dgemm_(const char* transa, const char* transb, const int* m, const int* n,
                       const int* k, const double* alpha, const double* a, const int* lda,
                       const double* b, const int* ldb, const double* beta, double* c,
                       const int* ldc);

using namespace dla_interface;
using namespace testing;

int main(int argc, char** argv) {
  // clang-format off
  cxxopts::Options desc(argv[0], "Allowed Options");
  desc.add_options()
      ("help", "produce help message")
      ("m,m", "Size m", cxxopts::value<int>()->default_value("1024"))
      ("n,n", "Size n", cxxopts::value<int>()->default_value("1024"))
      ("k,k", "Size k", cxxopts::value<int>()->default_value("1024"))
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

  const OpTrans transa = vm.count("transa") ? Trans : NoTrans;
  const OpTrans transb = vm.count("transb") ? Trans : NoTrans;

  const int a_m = transa == NoTrans ? m : k;
  const int a_n = transa == NoTrans ? k : m;
  const int b_m = transb == NoTrans ? k : n;
  const int b_n = transb == NoTrans ? n : k;

  std::vector<SolverType> solvers;
  if (vm.count("scalapack"))
    solvers.push_back(ScaLAPACK);
  if (vm.count("dplasma"))
    solvers.push_back(DPlasma);

  comm::CommunicatorManager::initialize(nr_threads, &argc, &argv, true);

  auto& comm_grid =
      comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, p, q, RowMajor);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> dist(-1, 1);
  auto random = [&dist, &rng](int, int) { return dist(rng); };

  double alpha = dist(rng);
  double beta = dist(rng);
  LocalMatrix<double> full_mat_a(a_m, a_n);
  fillLocalMatrix(full_mat_a, random);
  LocalMatrix<double> full_mat_b(b_m, b_n);
  fillLocalMatrix(full_mat_b, random);
  LocalMatrix<double> full_mat_c(m, n);
  fillLocalMatrix(full_mat_c, random);
  LocalMatrix<double> full_mat_c_res(full_mat_c);

  int lda = full_mat_a.leadingDimension();
  int ldb = full_mat_b.leadingDimension();
  int ldc = full_mat_c_res.leadingDimension();

  dgemm_((char*)&transa, (char*)&transb, &m, &n, &k, &alpha, full_mat_a.ptr(), &lda,
         full_mat_b.ptr(), &ldb, &beta, full_mat_c_res.ptr(), &ldc);

  DistributedMatrix<double> mat_a(a_m, a_n, a_mb, a_nb, comm_grid, scalapack_dist);
  DistributedMatrix<double> mat_b(b_m, b_n, b_mb, b_nb, comm_grid, scalapack_dist);
  DistributedMatrix<double> mat_c(m, n, c_mb, c_nb, comm_grid, scalapack_dist);

  auto val_a = [&full_mat_a](int i, int j) { return full_mat_a(i, j); };
  auto val_b = [&full_mat_b](int i, int j) { return full_mat_b(i, j); };
  auto val_c = [&full_mat_c](int i, int j) { return full_mat_c(i, j); };
  auto val_c_res = [&full_mat_c_res](int i, int j) { return full_mat_c_res(i, j); };

  fillDistributedMatrix(mat_a, val_a);
  fillDistributedMatrix(mat_b, val_b);

  bool status_all = true;

  for (auto solver : solvers) {
    fillDistributedMatrix(mat_c, val_c);

    matrixMultiplication(transa, transb, alpha, mat_a, mat_b, beta, mat_c, solver, 2);

    auto status = checkNearDistributedMatrix(mat_c, val_c_res, 1e-10, 1e-2);
    if (not status) {
      std::cout << "***" << util::getSolverString(solver) << " Failed! ***" << std::endl;
    }
    status_all = status_all && status;
  }
  return status_all ? 0 : 1;
}
