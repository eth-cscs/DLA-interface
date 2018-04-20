template <class ElType>
void matrixMultiplication(OpTrans trans_a, OpTrans trans_b, ElType alpha,
                          const DistributedMatrix<ElType>& mat_a,
                          const DistributedMatrix<ElType>& mat_b, ElType beta,
                          DistributedMatrix<ElType>& mat_c, SolverType solver, int print_timers) {
  dlai__util__debug_print_call_param(trans_a, trans_b, alpha, mat_a, mat_b, beta, mat_c, solver);
  auto& comm_grid = mat_c.commGrid();
  util::Timer<> timer_full(comm_grid.rowOrderedMPICommunicator(), print_timers > 0);

  solver = dlai__util__fallbackCommunicator(comm_grid, solver);

  // For real Types OpTrans 'C' and 'T' are equivalent.
  // Since DPlasma doesn't support 'C' for real types replace 'C' with 'T' for reals.
  trans_a = util::RemoveConjOpTypeIfReal<ElType>(trans_a);
  trans_b = util::RemoveConjOpTypeIfReal<ElType>(trans_b);

  int m = mat_c.size().first;
  int n = mat_c.size().second;
  int m2 = trans_a == NoTrans ? mat_a.size().first : mat_a.size().second;
  int k = trans_a == NoTrans ? mat_a.size().second : mat_a.size().first;
  int k2 = trans_b == NoTrans ? mat_b.size().first : mat_b.size().second;
  int n2 = trans_b == NoTrans ? mat_b.size().second : mat_b.size().first;
  if (m != m2 || n != n2 || k != k2)
    throw std::invalid_argument(errorMessage("Incompatible sizes in matrix multiplication: ", m,
                                             " != ", m2, " || ", n, " != ", n2, " || ", k, " != ",
                                             k2));

  double mnk = static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(k);
  double flop = util::nrOps<ElType>(mnk, mnk);
  switch (solver) {
#ifdef DLA_HAVE_SCALAPACK
    case ScaLAPACK: {
      std::array<int, 6> timer_index;
      util::Timer<> timer_part(comm_grid.rowOrderedMPICommunicator(), print_timers > 1);
      timer_index[0] = 0;
      {
        auto mat_a_scalapack_ptr = mat_a.convertConst(scalapack_dist);
        timer_index[1] = timer_part.save_time();
        decltype(mat_a_scalapack_ptr) mat_b_scalapack_ptr;
        // if mat_a and mat_b are the same matrix (same memory) only mat_a is converted,
        // otherwise convert both matrices.
        if (mat_a.isSameMatrix(mat_b))
          mat_b_scalapack_ptr = mat_a_scalapack_ptr;
        else
          mat_b_scalapack_ptr = mat_b.convertConst(scalapack_dist);

        timer_index[2] = timer_part.save_time();
        DistributedMatrix<ElType> mat_c_scalapack(scalapack_dist, mat_c);
        timer_index[3] = timer_part.save_time();
        auto matrix_a_info = mat_a_scalapack_ptr->getScalapackDescription();
        auto matrix_b_info = mat_b_scalapack_ptr->getScalapackDescription();
        auto matrix_c_info = mat_c_scalapack.getScalapackDescription();

        scalapack_wrappers::pgemm(
            trans_a, trans_b, m, n, k, alpha, std::get<0>(matrix_a_info),
            std::get<1>(matrix_a_info), std::get<2>(matrix_a_info), &std::get<3>(matrix_a_info)[0],
            std::get<0>(matrix_b_info), std::get<1>(matrix_b_info), std::get<2>(matrix_b_info),
            &std::get<3>(matrix_b_info)[0], beta, std::get<0>(matrix_c_info),
            std::get<1>(matrix_c_info), std::get<2>(matrix_c_info), &std::get<3>(matrix_c_info)[0]);

        timer_index[4] = timer_part.save_time();
      }
      timer_index[5] = timer_part.save_time();
      if (comm_grid.id2D() == std::make_pair(0, 0)) {
        timer_part.print_elapsed(timer_index[0], timer_index[1], "Conversion a: ");
        timer_part.print_elapsed(timer_index[1], timer_index[2], "Conversion b: ");
        timer_part.print_elapsed(timer_index[2], timer_index[3], "Conversion c: ");
        timer_part.print_elapsed(timer_index[3], timer_index[4],
                                 "Matrix Matrix Multiplication (ScaLAPACK) time: ", flop);
        timer_part.print_elapsed(timer_index[4], timer_index[5], "Back conversion c: ");
      }

      break;
    }
#endif

#ifdef DLA_HAVE_DPLASMA
    case DPlasma: {
      std::array<int, 6> timer_index;
      util::Timer<> timer_part(comm_grid.rowOrderedMPICommunicator(), print_timers > 1);
      timer_index[0] = 0;
      {
        auto mat_a_tile_ptr = mat_a.convertConst(tile_dist);
        timer_index[1] = timer_part.save_time();
        decltype(mat_a_tile_ptr) mat_b_tile_ptr;
        // if mat_a and mat_b are the same matrix (same memory) only mat_a is converted,
        // otherwise convert both matrices.
        if (mat_a.isSameMatrix(mat_b))
          mat_b_tile_ptr = mat_a_tile_ptr;
        else
          mat_b_tile_ptr = mat_b.convertConst(tile_dist);

        timer_index[2] = timer_part.save_time();
        DistributedMatrix<ElType> mat_c_tile(tile_dist, mat_c);
        timer_index[3] = timer_part.save_time();
        auto matrix_a_info = mat_a_tile_ptr->getDPlasmaDescription();
        auto matrix_b_info = mat_b_tile_ptr->getDPlasmaDescription();
        auto matrix_c_info = mat_c_tile.getDPlasmaDescription();
        const parsec_tiled_matrix_dc_t* dp_mat_a =
            reinterpret_cast<const parsec_tiled_matrix_dc_t*>(&std::get<0>(matrix_a_info));
        const parsec_tiled_matrix_dc_t* dp_mat_b =
            reinterpret_cast<const parsec_tiled_matrix_dc_t*>(&std::get<0>(matrix_b_info));
        parsec_tiled_matrix_dc_t* dp_mat_c =
            reinterpret_cast<parsec_tiled_matrix_dc_t*>(&std::get<0>(matrix_c_info));

        dplasma_wrappers::dplasma_run<dplasma_wrappers::pgemm<ElType>>(
            std::get<1>(matrix_c_info), dplasma_wrappers::plasmaTrans(trans_a),
            dplasma_wrappers::plasmaTrans(trans_b), util::castToC(alpha), dp_mat_a, dp_mat_b,
            util::castToC(beta), dp_mat_c);

        timer_index[4] = timer_part.save_time();
      }
      timer_index[5] = timer_part.save_time();
      if (comm_grid.id2D() == std::make_pair(0, 0)) {
        timer_part.print_elapsed(timer_index[0], timer_index[1], "Conversion a: ");
        timer_part.print_elapsed(timer_index[1], timer_index[2], "Conversion b: ");
        timer_part.print_elapsed(timer_index[2], timer_index[3], "Conversion c: ");
        timer_part.print_elapsed(timer_index[3], timer_index[4],
                                 "Matrix Matrix Multiplication (DPlasma) time: ", flop);
        timer_part.print_elapsed(timer_index[4], timer_index[5], "Back conversion c: ");
      }

      break;
    }
#endif

    default:
      throw std::invalid_argument(
          errorMessage("Matrix multiplication is not available for solver ", solver));
  }

  int timer_index_end = timer_full.save_time();
  if (comm_grid.id2D() == std::make_pair(0, 0)) {
    timer_full.print_elapsed(0, timer_index_end, "DLA Matrix Matrix Multiplication time: ", flop);
  }
}
