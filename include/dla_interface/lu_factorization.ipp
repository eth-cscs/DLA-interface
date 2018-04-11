template <class ElType>
void LUFactorization(DistributedMatrix<ElType>& mat, std::vector<int>& ipiv, SolverType solver,
                     int print_timers) {
  ipiv.resize(mat.localBaseIndex().row + mat.localSize().first + mat.blockSize().first);
  LUFactorization(mat, &ipiv[0], solver, print_timers);
}

template <class ElType>
void LUFactorization(DistributedMatrix<ElType>& mat, int* ipiv, SolverType solver, int print_timers) {
  dlai__util__debug_print_call_param(mat, solver);
  auto& comm_grid = mat.commGrid();
  util::Timer<> timer_full(comm_grid.rowOrderedMPICommunicator(), print_timers > 0);

  dlai__util__checkBlocksAreSquare(mat);
  dlai__util__checkBaseIndexAtBlock(mat);

  double x = std::min(mat.size().first, mat.size().second);
  double d = std::max(mat.size().first, mat.size().second) - x;
  double ops = x * x * (x / 3 + d / 2);
  double flop = util::nrOps<ElType>(ops, ops);

  switch (solver) {
#ifdef DLA_HAVE_SCALAPACK
    case ScaLAPACK: {
      std::array<int, 4> timer_index;
      util::Timer<> timer_part(comm_grid.rowOrderedMPICommunicator(), print_timers > 1);
      timer_index[0] = 0;
      int info = 0;
      {
        DistributedMatrix<ElType> mat_scalapack(scalapack_dist, mat);
        timer_index[1] = timer_part.save_time();
        auto matrix_info = mat_scalapack.getScalapackDescription();

        scalapack_wrappers::pgetrf(mat_scalapack.size().first, mat_scalapack.size().second,
                                   std::get<0>(matrix_info), std::get<1>(matrix_info),
                                   std::get<2>(matrix_info), &std::get<3>(matrix_info)[0], ipiv,
                                   info);

        timer_index[2] = timer_part.save_time();
      }
      timer_index[3] = timer_part.save_time();
      if (comm_grid.id2D() == std::make_pair(0, 0)) {
        timer_part.print_elapsed(timer_index[0], timer_index[1], "Conversion a: ");

        timer_part.print_elapsed(timer_index[1], timer_index[2],
                                 "LU Factorization (ScaLAPACK) time: ", flop);
        timer_part.print_elapsed(timer_index[2], timer_index[3], "Back conversion a: ");
      }
      if (info != 0) {
        if (info < 0)
          throw std::invalid_argument(errorMessage("Argument ", -info, " is wrong."));
        else
          throw std::invalid_argument(errorMessage("Matrix is singular (", info, ")"));
      }
      break;
    }
#endif

#ifdef DLA_HAVE_DPLASMA
    case DPlasma: {
      std::array<int, 5> timer_index;
      util::Timer<> timer_part(comm_grid.rowOrderedMPICommunicator(), print_timers > 1);
      timer_index[0] = 0;
      int info = 0;
      {
        DistributedMatrix<ElType> mat_tile(tile_dist, mat);
        timer_index[1] = timer_part.save_time();
        auto matrix_info = mat_tile.getDPlasmaDescription();
        parsec_tiled_matrix_dc_t* dp_mat =
            reinterpret_cast<parsec_tiled_matrix_dc_t*>(&std::get<0>(matrix_info));

        int ipiv_size = std::min(mat_tile.size().first, mat_tile.size().second);
        int nb = mat_tile.blockSize().second;
        auto ipiv_tile = DistributedMatrix<int>(1, mat_tile.baseIndex().col + ipiv_size, 1, nb,
                                                comm_grid, tile_dist)
                             .subMatrix(1, ipiv_size, 0, mat_tile.baseIndex().col);
        auto ipiv_info = ipiv_tile->getDPlasmaDescription();
        parsec_tiled_matrix_dc_t* dp_ipiv =
            reinterpret_cast<parsec_tiled_matrix_dc_t*>(&std::get<0>(ipiv_info));

        dplasma_wrappers::dplasma_run<dplasma_wrappers::pgetrf<ElType>>(std::get<1>(matrix_info),
                                                                        dp_mat, dp_ipiv, &info);

        timer_index[2] = timer_part.save_time();

        int rank;
        int size;
        MPI_Comm_rank(comm_grid.rowOrderedMPICommunicator(), &rank);
        MPI_Comm_size(comm_grid.rowOrderedMPICommunicator(), &size);
        int g_base_index = mat_tile.baseIndex().row;
        int l_base_index = mat_tile.localBaseIndex().row;

        std::memset(ipiv, 0, (g_base_index + mat_tile.size().first) * sizeof(int));
        for (int j = 0; j < ipiv_tile->localSize().second; ++j) {
          int g_index = ipiv_tile->getGlobal2DIndex(Local2DIndex(0, j)).col;
          ipiv[mat_tile.baseIndex().row + g_index] =
              ipiv_tile->operator()(Local2DIndex(0, j)) + (l_base_index + g_index) / nb * nb;
        }
        MPI_Allreduce(MPI_IN_PLACE, ipiv + g_base_index, ipiv_size, MPI_INT, MPI_SUM,
                      comm_grid.rowOrderedMPICommunicator());

        timer_index[3] = timer_part.save_time();
      }
      timer_index[4] = timer_part.save_time();
      if (comm_grid.id2D() == std::make_pair(0, 0)) {
        timer_part.print_elapsed(timer_index[0], timer_index[1], "Conversion a: ");

        timer_part.print_elapsed(timer_index[1], timer_index[2],
                                 "LU Factorization (DPlasma) time: ", flop);
        timer_part.print_elapsed(timer_index[2], timer_index[3], "Conversion ipiv: ");
        timer_part.print_elapsed(timer_index[3], timer_index[4], "Back conversion a: ");
      }
      if (info != 0) {
        throw std::invalid_argument(errorMessage("Matrix is not positive definite (", info, ")"));
      }
      break;
    }
#endif

    default:
      throw std::invalid_argument(
          errorMessage("LU factorization is not available for solver ", solver));
  }
  int index_end = timer_full.save_time();
  if (comm_grid.id2D() == std::make_pair(0, 0)) {
    timer_full.print_elapsed(0, index_end, "DLA LU Factorization time: ", flop);
  }

  return;
}
