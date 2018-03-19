#include <array>
#include <iostream>
#include <random>
#include <vector>
#include <boost/program_options.hpp>
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
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("m,m", boost::program_options::value<int>()->default_value(1024), "Size m")
      ("n,n", boost::program_options::value<int>()->default_value(1024), "Size n")
      ("k,k", boost::program_options::value<int>()->default_value(1024), "Size k")
      ("def_nb", boost::program_options::value<int>()->default_value(128), "The default block size for the block cyclic distribution. (Default = 128)")
      ("a_mb", boost::program_options::value<int>()->default_value(-1), "The row block size of a for the block cyclic distribution. If <= 0 def_nb is used. (Default = -1 i.e. def_nb)")
      ("a_nb", boost::program_options::value<int>()->default_value(-1), "The column block size of a for the block cyclic distribution. If <= 0 def_nb is used. (Default = -1 i.e. def_nb)")
      ("b_mb", boost::program_options::value<int>()->default_value(-1), "The row block size of b for the block cyclic distribution. If <= 0 def_nb is used. (Default = -1 i.e. def_nb)")
      ("b_nb", boost::program_options::value<int>()->default_value(-1), "The column block size of b for the block cyclic distribution. If <= 0 def_nb is used. (Default = -1 i.e. def_nb)")
      ("c_mb", boost::program_options::value<int>()->default_value(-1), "The row block size of c for the block cyclic distribution. If <= 0 def_nb is used. (Default = -1 i.e. def_nb)")
      ("c_nb", boost::program_options::value<int>()->default_value(-1), "The column block size of c for the block cyclic distribution. If <= 0 def_nb is used. (Default = -1 i.e. def_nb)")
      ("seed", boost::program_options::value<unsigned int>()->default_value(2017), "Seed for random number generator. (Default = 2017)")
      ("transa", "a is stored transposed.")
      ("transb", "b is stored transposed.")
      ("scalapack", "Run test with ScaLAPACK.")
      ("dplasma", "Run test with DPlasma.")
      ("row_procs,p", boost::program_options::value<int>()->default_value(1), "The number of rows in the 2D communicator. (Default = 1)")
      ("col_procs,q", boost::program_options::value<int>()->default_value(1), "The number of cols in the 2D communicator. (Default = 1)")
      ("nr_threads", boost::program_options::value<int>()->default_value(1), "The number of threads per rank. (Default 1)")
      ;
  // clang-format on

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
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
