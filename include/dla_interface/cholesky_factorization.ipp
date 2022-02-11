template <class ElType>
void choleskyFactorization(UpLo uplo, DistributedMatrix<ElType>& mat, SolverType solver,
                           int print_timers) {
  dlai__util__debug_print_call_param(uplo, mat, solver);
  auto& comm_grid = mat.commGrid();
  util::Timer<> timer_full(comm_grid.rowOrderedMPICommunicator(), print_timers > 0);

  if (uplo == All)
    uplo = Lower;

  dlai__util__checkIsSquare(mat);
  dlai__util__checkBlocksAreSquare(mat);
  dlai__util__checkBaseIndexAtBlock(mat);

  solver = dlai__util__fallbackCommunicator(comm_grid, solver);

  double n = mat.size().first;
  double n3 = n * n * n;
  double flop = util::nrOps<ElType>(n3 / 6, n3 / 6);
  
  std::cout << " using solver type: " << solver << std::endl;
  
  switch (solver) {

// Implementation for DLA-Future library call
#ifdef DLAI_WITH_DLAF
  case DLAF: {

	// Test compile/linking
	// std::cout << " -> DLAF: Inside Cholesky_Factorization()" << std::endl;


	std::array<int, 4> timer_index;
	util::Timer<> timer_part(comm_grid.rowOrderedMPICommunicator(), print_timers > 1);
	timer_index[0] = 0;
	int info = 0;

	DistributedMatrix<ElType> mat_scalapack(scalapack_dist, mat);
	Global2DIndex base_index = mat_scalapack.baseIndex();

	if (Global2DIndex(0, 0) != base_index) {
	  throw std::invalid_argument(
			  errorMessage(
			  "DLAF Cholesky requires baseIndex == (0, 0) for the input matrix, ",
			  base_index, "given")
			  );
	}
	timer_index[1] = timer_part.save_time();

	// *** Wrap communicator create part ***
	auto& comm_grid = mat_scalapack.commGrid();
	dlaf_wrappers::CommunicatorGrid dlaf_comm_grid = dlaf_wrappers::comm_grid(comm_grid);

	// *** Wrap matrix create part ***
	dlaf::Matrix<ElType, dlaf::Device::CPU> dlaf_matrix = dlaf_wrappers::matrix<ElType, dlaf::Device::CPU>(mat_scalapack, dlaf_comm_grid);

	// std::cout << "dlaf_wrappers::cholesky() - START" << std::endl;
	// *** Wrap Cholesky solver part ***
	dlaf_wrappers::cholesky<ElType, dlaf::Device::CPU>(dlaf_comm_grid, dlaf_matrix);
	// std::cout << "dlaf_wrappers::cholesky() - END" << std::endl;

	timer_index[2] = timer_part.save_time();
	timer_index[3] = timer_part.save_time();

	if (comm_grid.id2D() == std::make_pair(0, 0)) {
		timer_part.print_elapsed(timer_index[0], timer_index[1], "Conversion a: ");

		timer_part.print_elapsed(timer_index[1], timer_index[2],
								 "Cholesky Factorization (DLA-Future) time: ", flop);
		timer_part.print_elapsed(timer_index[2], timer_index[3], "Back conversion a: ");
	}

	if (info != 0) {
		throw std::invalid_argument(errorMessage("Matrix is not positive definite (", info, ")"));
	}

	break;
  }
#endif


#ifdef DLAI_WITH_SCALAPACK
    case ScaLAPACK: {
      std::array<int, 4> timer_index;
      util::Timer<> timer_part(comm_grid.rowOrderedMPICommunicator(), print_timers > 1);
      timer_index[0] = 0;
      int info = 0;
      {
        DistributedMatrix<ElType> mat_scalapack(scalapack_dist, mat);
        timer_index[1] = timer_part.save_time();
        auto matrix_info = mat_scalapack.getScalapackDescription();

        scalapack_wrappers::ppotrf(uplo, mat_scalapack.size().first, std::get<0>(matrix_info),
                                   std::get<1>(matrix_info), std::get<2>(matrix_info),
                                   &std::get<3>(matrix_info)[0], info);

        timer_index[2] = timer_part.save_time();
      }
      timer_index[3] = timer_part.save_time();
      if (comm_grid.id2D() == std::make_pair(0, 0)) {
        timer_part.print_elapsed(timer_index[0], timer_index[1], "Conversion a: ");

        timer_part.print_elapsed(timer_index[1], timer_index[2],
                                 "Cholesky Factorization (ScaLAPACK) time: ", flop);
        timer_part.print_elapsed(timer_index[2], timer_index[3], "Back conversion a: ");
      }
      if (info != 0) {
        if (info < 0)
          throw std::invalid_argument(errorMessage("Argument ", -info, " is wrong."));
        else
          throw std::invalid_argument(errorMessage("Matrix is not positive definite (", info, ")"));
      }
      break;
    }
#endif

#ifdef DLAI_WITH_DPLASMA
    case DPlasma: {
      std::array<int, 4> timer_index;
      util::Timer<> timer_part(comm_grid.rowOrderedMPICommunicator(), print_timers > 1);
      timer_index[0] = 0;
      int info = 0;
      {
        DistributedMatrix<ElType> mat_tile(tile_dist, mat);
        timer_index[1] = timer_part.save_time();
        auto matrix_info = mat_tile.getDPlasmaDescription();
        parsec_tiled_matrix_dc_t* dp_mat =
            reinterpret_cast<parsec_tiled_matrix_dc_t*>(&std::get<0>(matrix_info));

        info = dplasma_wrappers::dplasma_run<dplasma_wrappers::ppotrf<ElType>>(
            std::get<1>(matrix_info), dplasma_wrappers::plasmaUpLo(uplo), dp_mat);

        timer_index[2] = timer_part.save_time();
      }
      timer_index[3] = timer_part.save_time();
      if (comm_grid.id2D() == std::make_pair(0, 0)) {
        timer_part.print_elapsed(timer_index[0], timer_index[1], "Conversion a: ");

        timer_part.print_elapsed(timer_index[1], timer_index[2],
                                 "Cholesky Factorization (DPlasma) time: ", flop);
        timer_part.print_elapsed(timer_index[2], timer_index[3], "Back conversion a: ");
      }
      if (info != 0) {
        if (info < 0)
          throw std::invalid_argument(errorMessage("DPLASMA Error ", info));
        else
          throw std::invalid_argument(errorMessage("Matrix is not positive definite (", info, ")"));
      }
      break;
    }
#endif

    default:
      throw std::invalid_argument(
          errorMessage("Cholesky factorization is not available for solver ", solver));
  }

  int index_end = timer_full.save_time();
  if (comm_grid.id2D() == std::make_pair(0, 0)) {
    timer_full.print_elapsed(0, index_end, "DLA Cholesky Factorization time: ", flop);
  }

  return;
}
