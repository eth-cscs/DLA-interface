template <class ElType>
void hermitianEigenvectors(UpLo uplo, DistributedMatrix<ElType>& mat,
                           std::vector<BaseType<ElType>>& evalues,
                           DistributedMatrix<ElType>& evectors, SolverType solver, int print_timers) {
  evalues.resize(mat.size().first);
  hermitianEigenvectors(uplo, mat, &evalues[0], evectors, solver, print_timers);
}

template <class ElType>
void hermitianEigenvectors(UpLo uplo, DistributedMatrix<ElType>& mat, BaseType<ElType>* evalues,
                           DistributedMatrix<ElType>& evectors, SolverType solver, int print_timers) {
  dlai__util__debug_print_call_param(uplo, mat, evectors, solver, print_timers);
  auto& comm_grid = mat.commGrid();
  util::Timer<> timer_full(comm_grid.rowOrderedMPICommunicator(), print_timers > 0);

  dlai__util__checkIsSquare(mat);
  dlai__util__checkBlocksAreSquare(mat);
  dlai__util__checkBaseIndexAtBlock(mat);

#ifdef DLA_HAVE_ELPA
  if (solver == ELPA) {
    solver =
        dlai__util__fallbackScaLAPACKCondition(mat.baseIndex() != Global2DIndex(0, 0), comm_grid,
                                               solver, "ELPA supports only baseIndex = (0, 0).");
  }
#endif

  switch (solver) {
#ifdef DLA_HAVE_SCALAPACK
    case ScaLAPACK: {
      // use the lower part if A is fully set.
      if (uplo == All)
        uplo = Lower;
      std::array<int, 6> timer_index;
      util::Timer<> timer_part(comm_grid.rowOrderedMPICommunicator(), print_timers > 1);
      timer_index[0] = 0;
      int info = 0;
      {
        DistributedMatrix<ElType> mat_scalapack(scalapack_dist, mat);
        timer_index[1] = timer_part.save_time();
        DistributedMatrix<ElType> mat_ev_scalapack(scalapack_dist, evectors);
        timer_index[2] = timer_part.save_time();

        auto matrix_info = mat_scalapack.getScalapackDescription();
        auto matrix_ev_info = mat_ev_scalapack.getScalapackDescription();

        scalapack_wrappers::pheevd(  //
            uplo, mat_scalapack.size().first, std::get<0>(matrix_info), std::get<1>(matrix_info),
            std::get<2>(matrix_info), &std::get<3>(matrix_info)[0], evalues,
            std::get<0>(matrix_ev_info), std::get<1>(matrix_ev_info), std::get<2>(matrix_ev_info),
            &std::get<3>(matrix_ev_info)[0], info);

        timer_index[3] = timer_part.save_time();
      }
      timer_index[4] = timer_part.save_time();
      if (comm_grid.id2D() == std::make_pair(0, 0)) {
        timer_part.print_elapsed(timer_index[0], timer_index[1], "Conversion a: ");
        timer_part.print_elapsed(timer_index[1], timer_index[2], "Conversion ev: ");
        timer_part.print_elapsed(timer_index[2], timer_index[3],
                                 "D&C Eigensolver (ScaLAPACK) time: ");
        timer_part.print_elapsed(timer_index[3], timer_index[4], "Back conversion a, ev: ");
      }
      if (info != 0) {
        if (info < 0)
          throw std::invalid_argument(errorMessage("Argument ", -info, " is wrong."));
        else
          throw std::invalid_argument(errorMessage("D&C Eigensolver returned (", info, ")"));
      }

      break;
    }
#endif
#ifdef DLA_HAVE_ELPA
    case ELPA: {
      std::array<int, 6> timer_index;
      util::Timer<> timer_part(comm_grid.rowOrderedMPICommunicator(), print_timers > 1);
      timer_index[0] = 0;
      int info = 0;
      {
        DistributedMatrix<ElType> mat_scalapack(scalapack_dist, mat);
        timer_index[1] = timer_part.save_time();
        DistributedMatrix<ElType> mat_ev_scalapack(scalapack_dist, evectors);
        timer_index[2] = timer_part.save_time();

        auto matrix_info = mat_scalapack.getScalapackDescription();
        auto matrix_ev_info = mat_ev_scalapack.getScalapackDescription();

        if (uplo != All) {
          UpLo inv_uplo = uplo == Lower ? Upper : Lower;
          int ij = uplo == Lower ? 2 : 1;
          int ji = uplo == Lower ? 1 : 2;
          // matrix size - 1
          int n1 = mat_scalapack.size().first - 1;
          scalapack_wrappers::ptradd(inv_uplo, ConjTrans, n1, n1, 1., std::get<0>(matrix_info), ij,
                                     ji, &std::get<3>(matrix_info)[0], 0., std::get<0>(matrix_info),
                                     ji, ij, &std::get<3>(matrix_info)[0]);
        }

        elpa::Handle handle(mat_scalapack, mat_scalapack.size().first);
        elpa::eigenvectors(handle, std::get<0>(matrix_info), evalues, std::get<0>(matrix_ev_info),
                           &info);

        timer_index[3] = timer_part.save_time();
      }
      timer_index[4] = timer_part.save_time();
      if (comm_grid.id2D() == std::make_pair(0, 0)) {
        timer_part.print_elapsed(timer_index[0], timer_index[1], "Conversion a: ");
        timer_part.print_elapsed(timer_index[1], timer_index[2], "Conversion ev: ");
        timer_part.print_elapsed(timer_index[2], timer_index[3], "ELPA Eigensolver time: ");
        timer_part.print_elapsed(timer_index[3], timer_index[4], "Back conversion a, ev: ");
      }
      if (info != 0) {
        throw std::invalid_argument(
            errorMessage("ELPA Eigensolver returned (", info, " (", elpa_strerr(info), "))"));
      }

      break;
    }
#endif

    default:
      throw std::invalid_argument(errorMessage("Eigensolver is not available for solver ", solver));
  }

  int timer_index_end = timer_full.save_time();
  if (comm_grid.id2D() == std::make_pair(0, 0)) {
    timer_full.print_elapsed(0, timer_index_end, "DLA Eigensolver time: ");
  }
}
