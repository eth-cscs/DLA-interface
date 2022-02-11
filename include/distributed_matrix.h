#ifndef DLA_INTERFACE_DISTRIBUTED_MATRIX_H
#define DLA_INTERFACE_DISTRIBUTED_MATRIX_H

#include <array>
#include <mpi.h>
#include <tuple>
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "matrix_index.h"
#include "local_matrix.h"
#include "scalapack.h"
#include "types.h"

namespace dla_interface {

/// A class DistributedMatrix creates a (m x n) matrix distributed in (mb x nb) blocks.

///<br><br>
/// \note
/// <b>(1)</b> \anchor note_01<br>
/// If the matrix distribution matches the distribution of the original matrix:
/// <ul>
///     <li>the matrix reference directly the elements given by the ptr,</li>
///     <li>the matrix element and the original matrix elements are the same,</li>
///     <li>the leading dimension and leading number of blocks remain unchanged.</li>
/// </ul>
/// If the matrix distribution does not match the distribution of the original matrix:
/// <ul>
///     <li>new memory is allocated for the constructed matrix,</li>
///     <li>the leading dimension and leading number of blocks are not the same,</li>
///     <li>the elements of the original matrix are copied to the new matrix (the elements are redistributed),</li>
///     <li>On destruction of the DistributedMatrix object, the elements are copied back to the original matrix (redistributing them).</li>
/// </ul>
/// <b>Requirement</b>: The lifetime of the original matrix has to be longer than the lifetime of the DistributedMatrix object.
///
/// <br><br>
/// \note
/// <b>(2)</b> \anchor note_02<br>
/// The returned Local/DistributedMatrix will use part of the memory allocated for the original matrix.
/// <ul>
/// 		<li>If the original matrix is destroyed, the memory of the original matrix is not deallocated
///          until all non zero size submatrices are destroyed.</li>
/// 		<li>If elements are changed in the original matrix, the submatrix will change as well and
///          vice-versa.</li>
/// 		<li>If the original matrix is re-assigned using the assignement or move-assignement operator
///          the submatrix object is not modified.
///  	<li><b>Exception:</b> if the memory has been provided by the user (see \ref note_01 "Note (1)"), the user has to
///          guarantee that the memory is valid for the lifetime of the submatrix as well.</li>
/// </ul>

  template <class ElType>
  class DistributedMatrix {
    public:
    using ElementType = ElType;

    /// Creates a (0 x 0) matrix with blocksize (1 x 1).
    DistributedMatrix();

	/// Creates a (m x n) matrix (or a (0 x 0) matrix if m == 0 or n == 0),
	/// distributed on the rank grid described by comm according the given distribution
	/// using the blocksize (mb x nb).
	///
	/// @throws std::invalid_argument if m < 0, n < 0, mb < 1, nb < 1.
	///
	/// @param m Number of rows.
	/// @param n Number of columns
	/// @param mb Distributed blocksize number of rows
	/// @param nb Distributed blocksize number of columns
	/// @param comm Described by Communicator2DGrid
	/// @param distribution enum DistributionType
	/// <ul>
	/// <li>scalapack_dist - Scalapack is 2D block cyclic distribution</li>
	/// <li>tile_dist - Tile is 2D block cyclic tile distribution</li>
	/// </ul>
	///
    DistributedMatrix(SizeType m, SizeType n, SizeType mb, SizeType nb,
                      const comm::Communicator2DGrid& comm, DistributionType distribution);

	/// Creates a (m x n) matrix (or a (0 x 0) matrix if m == 0 or n == 0),
	/// distributed on the rank grid described by comm according the given distribution
	/// using the blocksize (mb x nb) and leading dimension ld.
	///
	/// @throws std::invalid_argument
	/// <ul>
	/// 	<li>if m < 0, n < 0, mb < 1, nb < 1, or ld < 1</li>
	/// 	<li>if ld < local_m and n > 0 and distribution == scalapack_dist</li>
	/// 	<li>if ld < mb and distribution == tile_dist</li>
	/// </ul>
	///
	/// @param m Number of rows.
	/// @param n Number of columns
	/// @param mb Distributed blocksize number of rows
	/// @param nb Distributed blocksize number of columns
	/// @param ld SizeType
	/// @param comm Described by Communicator2DGrid
	/// @param distribution enum DistributionType
	/// <ul>
	/// <li>scalapack_dist - Scalapack is 2D block cyclic distribution</li>
	/// <li>tile_dist - Tile is 2D block cyclic tile distribution</li>
	/// </ul>
    DistributedMatrix(SizeType m, SizeType n, SizeType mb, SizeType nb, SizeType ld,
                      const comm::Communicator2DGrid& comm, DistributionType distribution);

	/// Creates a (m x n) matrix (or a (0 x 0) matrix if m == 0 or n == 0),
	/// distributed on the rank grid described by comm according the given distribution
	/// using the blocksize (mb x nb) and leading dimension ld.
	/// This matrix uses the elements stored in the memory provided by ptr (of length len elements)
	/// with leading dimension ld, and leading number of blocks leading_nr_blocks and
	/// distributed according to original_distribution (see \ref note_01 "Note (1)").
	///
	/// <b>Note:</b> if original_distribution == scalapack_dist leading_nr_blocks is ignored.
	///
	/// @throws std::invalid_argument
	/// <ul>
	/// 	<li>if m < 0, n < 0, mb < 1, nb < 1, or ld < 1</li>
	/// 	<li>if ld < local_m and n > 0 and distribution == scalapack_dist</li>
	/// 	<li>if ld < mb and distribution == tile_dist</li>
	/// 	<li>if leading_nr_blocks < ceil(local_m / mb) and distribution == tile_dist</li>
	/// </ul>
	///
	///
	/// @param m Number of rows.
	/// @param n Number of columns
	/// @param mb Distributed blocksize number of rows
	/// @param nb Distributed blocksize number of columns
	/// @param comm Described by Communicator2DGrid
	/// @param new_distribution enum DistributionType
	/// <ul>
	/// <li>scalapack_dist - Scalapack is 2D block cyclic distribution</li>
	/// <li>tile_dist - Tileis 2D block cyclic tile distribution</li>
	/// </ul>
	/// @param ptr Pointer to matrix in memory
	/// @param len Length of matrix elements in memory
	/// @param ld Leading dimension of matrix in memory
	/// @param leading_nr_blocks Leading number of blocks of matrix in memory
	/// @param original_distribution DistributionType of matrix in memory
	/// <ul>
	/// <li>scalapack_dist - Scalapack is 2D block cyclic distribution</li>
	/// <li>tile_dist - Tile is 2D block cyclic tile distribution</li>
	/// </ul>
	/// <b>Note:</b> if original_distribution == scalapack_dist leading_nr_blocks is ignored.
    DistributedMatrix(SizeType m, SizeType n, SizeType mb, SizeType nb,
                      const comm::Communicator2DGrid& comm, DistributionType new_distribution,
                      ElementType* ptr, std::size_t len, SizeType ld, SizeType leading_nr_blocks,
                      DistributionType original_distribution);

	/// Creates a distributed matrix which has the same size, block size and rank grid
	/// of the given DistributedMatrix distributed according distribution
	/// and which reference the same memory (see \ref note_01 "Note (1)").
	///
	/// <b>Note:</b> if original_distribution == scalapack_dist leading_nr_blocks is ignored.
	///
	/// @param distribution DistributionType of matrix in memory
	/// <ul>
	/// <li>scalapack_dist - Scalapack is 2D block cyclic distribution</li>
	/// <li>tile_dist - Tile is 2D block cyclic tile distribution</li>
	/// </ul>
	/// <b>Note:</b> if original_distribution == scalapack_dist leading_nr_blocks is ignored.
	/// @param mat Reference to matrix in memory
    DistributedMatrix(DistributionType distribution, DistributedMatrix& mat);

    #ifdef DLAI_WITH_SCALAPACK
	/// <b>Enabled</b> only if <b>DLAI_WITH_SCALAPACK</b> is defined!
	///
	/// Creates a distributed matrix which has the same size, block size and rank grid
	/// of the given Scalapack matrix and which reference the same memory (see \ref note_01 "Note (1)").
	///
	/// @throws std::invalid_argument
	/// <ul>
	/// 	<li>if the scalapack context was not created via the CommunicatorManager interface</li>
	/// 	<li>if the scalapack descriptor does not describe a valid Scalapack matrix</li>
	/// 	<li>if the indices i < 1, j < 1</li>
	/// 	<li>if the sizes m < 0, n < 0</li>
	/// 	<li>if the indices i + m - 1, j + n - 1 are larger than the corresponding matrix global sizes</li>
	/// 	<li>if desc[6] (rsrc) != 0 or desc[7] (csrc) != 0</li>	 ///
	/// </ul>
	///
	/// @param distribution DistributionType of matrix in memory
	/// <ul>
	/// <li>scalapack_dist - Scalapack is 2D block cyclic distribution</li>
	/// <li>tile_dist - Tile is 2D block cyclic tile distribution</li>
	/// </ul>
	/// @param m Matrix number of rows.
	/// @param n Matrix number of columns
	/// @param ptr Pointer to matrix in memory
	/// @param i Starting row of matrix
	/// @param j Starting column of matrix
	/// @param desc ScaLapack descriptor
    DistributedMatrix(DistributionType distribution, int m, int n, ElementType* ptr,
                      ScalapackIndex i, ScalapackIndex j, constScalapackDescriptor desc);
#endif
#ifdef DLAI_WITH_DPLASMA
	/// <b>Enabled</b> only if <b>DLAI_WITH_DPLASMA</b> is defined!
	///
	/// @todo check this and test it.
	///
	/// Creates a distributed matrix which has the same size, block size and rank grid
	/// of the given D-PLASMA matrix and which reference the same memory (see \ref note_01 "Note (1)").
	///
	/// DistributedMatrix(DistributionType distribution,   desc);
	///

#endif

#ifdef DLAI_WITH_DLAF
	/// <b>Enabled</b> only if <b>DLAI_WITH_DLAF</b> is defined!
	///
	/// @todo check this and test it.
	///
	/// Creates a distributed matrix which has the same size, block size and rank grid
	/// of the given DLAF matrix and which reference the same memory (see \ref note_01 "Note (1)").
	///
	/// DistributedMatrix(DistributionType distribution, DPlasmaDescriptor desc);
	///
#endif

	/// Creates a copy of rhs, i.e. a distributed matrix with the same size, block size,
	/// istribution, leading dimension, rank grid and a copy of the elements.
	///
	/// @param rhs Const reference to DistributedMatrix.
    DistributedMatrix(const DistributedMatrix& rhs);

    ~DistributedMatrix();

    /// Creates a distributed matrix of the same size, block size, distribution,
    /// leading dimension and rank grid moving the elements of rhs.
    /// <b>Postcondition:</b> rhs is in an unspecified state.
    ///
    /// @param rhs
    DistributedMatrix(DistributedMatrix&& rhs) noexcept;

    /// <i>Copy assignement operator</i><br>
    /// Changes the size, block size, distribution, leading dimension and rank grid
    /// of *this to match rhs values and copies the elements of rhs.

    /// \note
    /// <ul>
    /// 	<li>new memory is allocated,</li>
    ///     <li>old memory is deallocated if:</li>
    /// 	<ul>
    ///     	<li>the memory was allocated via the DistributedMatrixClass constructors and</li>
    ///         <li>the memory is not used by any other Distributed/LocalMatrixClass objects.</li>
    /// 	</ul>
    /// </ul>
    /// All the pointers and references to elements of *this are invalidated.
    ///
    /// @param rhs constant reference to DistributedMatrix
    DistributedMatrix& operator=(const DistributedMatrix& rhs);

    /// <i>Move-assignement operator</i><br>
    /// Changes the size, block size, distribution, leading dimension and rank grid
    /// of *this to match rhs values and moves the elements of rhs.
    /// All the pointers and references to elements of *this and rhs are invalidated.
    /// <b>Postcondition:</b> rhs is in a unspecified state.
    ///
    /// @param rhs reference to DistributedMatrix
    DistributedMatrix& operator=(DistributedMatrix&& rhs) noexcept;

    /// Returns a ptr to a new const distributed matrix which has the same size, block size and rank
    /// grid of the given DistributedMatrix distributed according distribution
    /// and which reference the same memory (see \ref note_01 "Note (1)").
    ///
    /// @param distribution DistributionType
    /// @return ptr to DistributedMatrix
    std::shared_ptr<const DistributedMatrix> convertConst(DistributionType distribution) const {
      DistributedMatrix<ElementType>& tmp_mat = const_cast<DistributedMatrix<ElementType>&>(*this);
      auto res = std::make_shared<DistributedMatrix<ElementType>>(distribution, tmp_mat);
      res->remove_referenced();
      return res;
    }

    /// Returns a ptr to a new const (m x n) matrix (or a (0 x 0) matrix if m == 0 or n == 0),
    /// distributed on the rank grid described by comm according to the distribution
    /// current_distribution using the blocksize (mb x nb).
    /// This matrix uses the elements stored in the memory provided by ptr (of length len elements)
    /// with leading dimension ld, and leading number of blocks leading_nr_blocks and
    /// distributed according to original_distribution.
    /// (see \ref note_01 "Note (1)").
    /// \throws std::invalid_argument
    /// <ul>
    /// 	<li>if m < 0, n < 0, mb < 1, nb < 1, or ld < 1</li>
    ///     </li>if ld < local_m and n > 0 and distribution == scalapack_dist</li>
    ///     </li>if ld < mb and distribution == tile_dist</li>
    ///     </li>if leading_nr_blocks < ceil(local_m / mb) and distribution == tile_dist</li>
    /// </ul>
    /// <br>
    ///
    /// \note if original_distribution == scalapack_dist leading_nr_blocks is ignored.
    ///
    /// @param m
    /// @param n
    /// @param mb
    /// @param nb
    /// @param comm
    /// @param new_distribution
    /// @param ptr
    /// @param len
    /// @param ld
    /// @param leading_nr_blocks
    /// @param original_distribution
    /// @return
    static std::shared_ptr<const DistributedMatrix> convertConst(
        SizeType m, SizeType n, SizeType mb, SizeType nb, const comm::Communicator2DGrid& comm,
        DistributionType new_distribution, const ElementType* ptr, std::size_t len, SizeType ld,
        SizeType leading_nr_blocks, DistributionType original_distribution) {
      ElementType* tmp_ptr = const_cast<ElementType*>(ptr);
      auto res = std::make_shared<DistributedMatrix<ElementType>>(
          m, n, mb, nb, comm, new_distribution, tmp_ptr, len, ld, leading_nr_blocks,
          original_distribution);
      res->remove_referenced();
      return res;
    }

#ifdef DLAI_WITH_SCALAPACK
    /// <b>Enabled</b> only if <b>DLAI_WITH_SCALAPACK</b> is defined!<br>
    ///
    /// Returns a ptr to a new const distributed matrix which has the same size, block size and rank
    /// grid of the given Scalapack matrix and which reference the same memory (see \ref note_01 "Note (1)")..
    /// \throws std::invalid_argument
    /// <ul>
    /// 	<li>if the scalapack context was not created via the CommunicatorManager interface.</li>
    ///     <li>if the scalapack descriptor does not describe a valid Scalapack matrix</li>
    ///     <li>if the indices i < 1, j < 1,</li>
    ///     <li>if the sizes m < 0, n < 0,</li>
    ///     <li>if the indices i + m - 1, j + n - 1 are larger than the corresponding matrix global sizes.</li>
    ///     <li>if desc[6] (rsrc) != 0 or desc[7] (csrc) != 0.</li>
    /// </ul>
    static std::shared_ptr<const DistributedMatrix> convertConst(DistributionType distribution,
                                                                 int m, int n, const ElementType* ptr,
                                                                 ScalapackIndex i, ScalapackIndex j,
                                                                 constScalapackDescriptor desc) {
      ElementType* tmp_ptr = const_cast<ElementType*>(ptr);
      auto res =
          std::make_shared<DistributedMatrix<ElementType>>(distribution, m, n, tmp_ptr, i, j, desc);
      res->remove_referenced();
      return res;
    }
#endif

    /// Returns true if this and rhs are the same matrix, i.e.
    /// <ul>
    /// 	<li>*this and rhs have the same shape,</li>
    ///     <li>*this and rhs are distributed in the same way,</li>
    ///     <li>*this and rhs reference the same memory.</li>
    /// </ul>
    ///
    /// @param rhs Reference to DistributedMatrix
    /// @return bool true if rhs == *this
    bool isSameMatrix(const DistributedMatrix& rhs) const;

    /// Copies the value of the elements of rhs.
    /// \throws std::invalid_argument if *this and rhs do not have the same size.<br>
    ///
    /// If size() and rhs.size() != (0, 0):
    /// \throws std::invalid_argument if *this and rhs do not have the same block size, rank grid,
    ///     element node distribution (i.e. if the (i, j)-th element of this and (i, j)-th element of
    ///     rhs are on different nodes for any (i, j)).
    ///
    /// @param rhs const ref to DistributedMatrix
    /// @return ref to DistributedMatrix
    DistributedMatrix& copy(const DistributedMatrix& rhs);

    /// Creates a DistributedMatrix which represent the (m, n) submatrix starting at global index (i,j),
    /// and returns a shared_ptr to it.
    ///
    /// \throws std::invalid_argument if m, n, i, j < 0, i + m > size().first, j + n > size().second
    /// (see \ref note_02 "Note (2)").
    ///
    /// @param m
    /// @param n
    /// @param i
    /// @param j
    /// @return
    std::shared_ptr<const DistributedMatrix<ElType>> subMatrix(SizeType m, SizeType n, IndexType i,
                                                               IndexType j) const {
      return std::shared_ptr<const DistributedMatrix<ElementType>>(
          subMatrixInternal(__func__, m, n, Global2DIndex(i, j)));
    }
    std::shared_ptr<DistributedMatrix<ElType>> subMatrix(SizeType m, SizeType n, IndexType i,
                                                         IndexType j) {
      return std::shared_ptr<DistributedMatrix<ElementType>>(
          subMatrixInternal(__func__, m, n, Global2DIndex(i, j)));
    }

    /// Creates a LocalMatrix which represent the (m, n) submatrix starting at global index (i, j),
    /// and returns a shared_ptr to it.
    ///
    /// \throws std::invalid_argument
    /// <ul>
    /// 	<li>if m, n, i, j < 0, i + m > size().first, j + n > size().second</li>
    ///     <li>if (baseIndex().first + i) % mb + m > mb or (baseIndex().second + j) % nb + n > nb.</li>
    ///	</ul>
    /// (see \ref note_02 "Note (2)").
    ///
    /// @param m
    /// @param n
    /// @param i
    /// @param j
    /// @return
    std::shared_ptr<const LocalMatrix<ElType>> localSubMatrix(SizeType m, SizeType n, IndexType i,
                                                              IndexType j) const {
      return std::shared_ptr<const LocalMatrix<ElementType>>(
          localSubMatrixInternal(__func__, m, n, Global2DIndex(i, j)));
    }
    std::shared_ptr<LocalMatrix<ElType>> localSubMatrix(SizeType m, SizeType n, IndexType i,
                                                        IndexType j) {
      return std::shared_ptr<LocalMatrix<ElementType>>(
          localSubMatrixInternal(__func__, m, n, Global2DIndex(i, j)));
    }

    /// Returns matrix base index
    ///
    /// @return Global2DIndex
    Global2DIndex baseIndex() const {
      return base_index_;
    }

    /// Returns the matrix size.
    ///
    /// @return matrix size as std::pair
    std::pair<SizeType, SizeType> size() const {
      return size_;
    }

    /// Returns matrix local base index.
    ///
    /// @return Local2DIndex
    Local2DIndex localBaseIndex() const {
      return local_base_index_;
    }

    /// Returns local size matrix.
    ///
    /// @return local size as std::pair
    std::pair<SizeType, SizeType> localSize() const {
      return local_size_;
    }

    /// Returns block size.
    ///
    /// @return block size as std::pair
    std::pair<SizeType, SizeType> blockSize() const {
      return block_size_;
    }

    /// Returns the number of elements allocated in the memory.
    /// Note: Different matrices can share the same memory.
    ///
    /// @return number
    std::size_t localStorageSize() const {
      return ptr_->size();
    }

    /// Returns leading matrix dimensions.
    ///
    /// @returns
    /// <td>
    /// 	<li>the leading dimension of the local matrix if the matrix is a scalapack matrix,</li>
    /// 	<li>the leading dimension of each tile if the matrix is a tile matrix otherwise.</li>
    /// </td>
    SizeType leadingDimension() const {
      return ld_;
    }

    /// Returns leading number of blocks.
    ///
    /// @returns
    /// <td>
    /// 	<li>1 if the matrix is a scalapack matrix,</li>
    /// 	<li>the number of blocks stored in each column of tiles otherwise.</li>
    /// </td>
    SizeType leadingNumberOfBlocks() const {
      return leading_nr_blocks_;
    }

    /// Returns matrix distribution.
    ///
    /// @return DistributionType
    DistributionType distribution() const {
      return distribution_;
    }

    /// Returns the matrix communicator.
    ///
    /// @return comm::Communicator2DGrid
    /// <b>Precondition:</b> the matrix is a non default constructed matrix.
    const comm::Communicator2DGrid& commGrid() const {
      assert(comm_grid_ != nullptr);
      return *comm_grid_;
    }

    /// Returns the index of the position in memory of the global (i, j) element.
    /// It holds ptr() + getLocalStorageIndex(index) == ptr(index).
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if i < 0 or i > size().first,</li>
    ///     <li>if j < 0 or j > size().second,</li>
    ///     <li>if the element is not owned by this rank.</li>
    /// <ul>
    ///
    /// @return local index
    std::size_t getLocalStorageIndex(IndexType i, IndexType j) const {
      return getLocalStorageIndex(Global2DIndex(i, j));
    }

    /// Same as getLocalStorageIndex(IndexType i, IndexType j) with i = index.row and j = index.col.
    std::size_t getLocalStorageIndex(Global2DIndex index) const;

    /// Returns the index of the position in memory of the local (index.row, index.col) element.
    /// It holds ptr() + getLocalStorageIndex(index) == ptr(index).
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if index.row < 0 or index.row > localSize().first,</li>
    ///     <li>if index.col < 0 or index.col > localSize().second.</li>
    /// <ul>
    ///
    /// @return local index
    std::size_t getLocalStorageIndex(Local2DIndex index) const;

    /// Returns the local 2D index of the global (i, j) element if owned by this rank.
    /// Otherwise returns the local 2D index of the global (i1, j1) element,
    /// where i1 is the smallest index > i and j1 is the smallest index > j
    /// such that the global element (i1, j1) is owned by this rank.
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if i < 0 or i > size().first,</li>
    ///     <li>if j < 0 or j > size().second.</li>
    /// </ul>
    ///
    /// @return local 2D index
    Local2DIndex getLocal2DIndex(IndexType i, IndexType j) const {
      return getLocal2DIndex(Global2DIndex(i, j));
    }

    /// Same as getLocal2DIndex(IndexType i, IndexType j) with i = index.row and j = index.col.
    Local2DIndex getLocal2DIndex(Global2DIndex index) const;

    /// Returns the 2D id of the rank which own the global (i, j) element.
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if i < 0 or i > size().first,</li>
    ///     <li>if j < 0 or j > size().second.</li>
    /// </ul>
    std::pair<int, int> getRankId2D(IndexType i, IndexType j) const {
      return getRankId2D(Global2DIndex(i, j));
    }
    /// Same as getRankId2D(IndexType i, IndexType j) with i = index.row and j = index.col.
    std::pair<int, int> getRankId2D(Global2DIndex index) const;

    /// Returns the global 2D index of the local (index.row, index.col) element.
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if index.row < 0 or index.row > size().first,</li>
    ///     <li>if index.col < 0 or index.col > size().second.</li>
    /// </ul>
    ///
    /// @return global 2D index
    Global2DIndex getGlobal2DIndex(Local2DIndex index) const;

    /// Returns a reference to the global (i, j) element.
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if i < 0 or i > size().first,</li>
    ///     <li>if j < 0 or j > size().second,</li>
    ///     <li>if the element is not owned by this rank.</li>
    /// </ul>
    ///
    /// @return const ref ElementType
    const ElementType& operator()(IndexType i, IndexType j) const {
      return *ptr(i, j);
    }

    /// Same as operator()(IndexType i, IndexType j) with i = index.row and j = index.col.
    ElementType& operator()(IndexType i, IndexType j) {
      return *ptr(i, j);
    }
    /// Same as operator()(IndexType i, IndexType j) with i = index.row and j = index.col.
    const ElementType& operator()(Global2DIndex index) const {
      return *ptr(index);
    }
    /// Same as operator()(IndexType i, IndexType j) with i = index.row and j = index.col.
    ElementType& operator()(Global2DIndex index) {
      return *ptr(index);
    }

    /// Returns a reference to the local (index.row, index.col) element.
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if index.row < 0 or index.row > localSize().first,</li>
    ///     <li>if index.col < 0 or index.col > localSize().second.</li>
    /// </ul>
    const ElementType& operator()(Local2DIndex index) const {
      return *ptr(index);
    }
    /// Same as operator()(Local2DIndex index) with index = local2Dindex.
    ElementType& operator()(Local2DIndex index) {
      return *ptr(index);
    }

    /// Returns nullptr if localSize().first == 0 and localSize.second == 0.
    /// Returns the pointer to the first element allocated in the memory otherwise.
    /// If ptr() != nullptr it holds ptr() + getLocalStorageIndex(index) == ptr(index),
    /// where index is either a GlobalIndex or a LocaIndex.
    const ElementType* ptr() const;
    ElementType* ptr();

    /// Returns the pointer to the global (i, j) element.
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if i < 0 or i > size().first,</li>
    ///     <li>if j < 0 or j > size().second,</li>
    ///     <li>if the element is not owned by this rank.</li>
    /// </ul>
    const ElementType* ptr(IndexType i, IndexType j) const {
      return ptr(Global2DIndex(i, j));
    }
    ElementType* ptr(IndexType i, IndexType j) {
      return ptr(Global2DIndex(i, j));
    }
    /// Same as ptr(IndexType i, IndexType j) with i = index.row and j = index.col.
    const ElementType* ptr(Global2DIndex index) const;
    ElementType* ptr(Global2DIndex index);

    /// Returns the pointer to the local (index.row, index.col) element.
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if index,row < 0 or index.row > localSize().first,</li>
    ///     <li>if index.col < 0 or index.col > localSize().second.</li>
    /// </ul>
    const ElementType* ptr(Local2DIndex index) const;
    ElementType* ptr(Local2DIndex index);

#ifdef DLAI_WITH_SCALAPACK
    /// <b>Enabled</b> only if <b>DLAI_WITH_SCALAPACK</b> is defined!<br>
    ///
    /// Returns a tuple containing ptr, i, j and a std::array containing the descriptor, needed for
    /// ScaLAPACK calls.
    ///
    /// @throws std::invalid_argument
    /// if the matrix is not distributed in the ScaLAPACK way.
    ///
    std::tuple<ElementType*, IndexType, IndexType, std::array<int, 9>> getScalapackDescription();
    std::tuple<const ElementType*, IndexType, IndexType, std::array<int, 9>> getScalapackDescription() const;
#endif
#ifdef DLAI_WITH_DPLASMA
    /// <b>Enabled</b> only if <b>DLAI_WITH_DPLASMA</b> is defined!<br>
    private:
    DPlasmaDescriptor getDPlasmaDescriptionInternal() const;

    public:
    /// Returns a tuple containing the DPLASMA descriptor, the MPI comunicator
    /// and a flag to indicate if the matrix is const.
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if the matrix is not tile distributed,</li>
    ///     <li>if leadingDimension() != blockSize().first.</li>
    /// </ul>
    std::tuple<DPlasmaDescriptor, MPI_Comm> getDPlasmaDescription();
    std::tuple<const DPlasmaDescriptor, MPI_Comm> getDPlasmaDescription() const;
#endif

#ifdef DLAI_WITH_DLAF
    /*
    /// <b>Enabled</b> only if <b>DLAI_WITH_DLAF</b> is defined!<br>
    private:
    DLAFDescriptor getDLAFDescriptionInternal() const;

    public:
    /// Returns a tuple containing the DLAF descriptor, the MPI comunicator
    /// and a flag to indicate if the matrix is const.
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if the matrix is not tile distributed,</li>
    ///     <li>if leadingDimension() != blockSize().first.</li>
    /// </ul>
    std::tuple<DLAFDescriptor, MPI_Comm> getDLAFDescription();
    std::tuple<const DLAFDescriptor, MPI_Comm> getDLAFDescription() const;
    */
#endif

    template <class Out>
    void debugDump(Out& out) const;

    template <class Out>
    void debugInfo(Out& out) const;

    /*************************/
    /*** Private interface ***/
    /*************************/

    private:
    // new storage
    DistributedMatrix(const char* func, std::pair<SizeType, SizeType> size,
                      std::pair<SizeType, SizeType> block_size, bool ld_set, SizeType ld,
                      const comm::Communicator2DGrid* comm, DistributionType distribution);

    // reference to allocated storage
    DistributedMatrix(const char* func, std::pair<SizeType, SizeType> size,
                      Global2DIndex base_index, std::pair<SizeType, SizeType> block_size,
                      const comm::Communicator2DGrid* comm, DistributionType distribution,
                      std::shared_ptr<memory::MemoryAllocator<ElementType>> original_ptr,
                      SizeType original_ld, SizeType original_leading_nr_blocks,
                      DistributionType original_distribution);

    void checkLocalIndex(const char* func, const Local2DIndex& index) const;
    void checkGlobalIndex(const char* func, bool check_rank, const Global2DIndex& index) const;

    Global2DIndex globalBaseIndexFromLocalBaseIndex(const char* func, const Local2DIndex& index) const;
    std::size_t storageBaseIndexFromLocalBaseIndex(const char* func, const Local2DIndex& index) const;
    static std::size_t storageBaseIndexFromLocalBaseIndex(
        const char* func, const Local2DIndex& index, const std::pair<SizeType, SizeType>& block_size,
        SizeType ld, SizeType leading_nr_blocks, DistributionType dist);
    Local2DIndex localBaseIndexFromGlobalBaseIndex(const char* func, const Global2DIndex& index) const;
    std::size_t storageBaseIndexFromGlobalBaseIndex(const char* func,
                                                    const Global2DIndex& index) const;
    std::pair<int, int> rankFromBaseGlobalIndex(const char* func, const Global2DIndex& index) const;

    static IndexType index1DBaseGlobalFromBaseLocal(IndexType index, SizeType block_size,
                                                    int rank_src, int rank, int comm_size);
    static int rank1DFromBaseGlobalIndex(IndexType index, SizeType block_size, int rank_src,
                                         int comm_size);
    static IndexType index1DBaseLocalFromBaseGlobal(IndexType index, SizeType block_size,
                                                    int rank_src, int rank, int comm_size);

    void checkAndComputeLocalParam(const char* func, bool ld_set, bool ld_nr_bl_set);
    void checkOrSetLeadingDims(const char* func, SizeType& ld, bool ld_set, SizeType& ld_nr_blks,
                               bool ld_nr_bl_set, DistributionType distribution);
    std::size_t allocationSize(const char* func) const;
    std::size_t allocationSize(const char* func, SizeType ld, SizeType leading_nr_blocks,
                               DistributionType distribution) const;

    void copyInternal(bool copy_back, std::shared_ptr<memory::MemoryAllocator<ElementType>> rhs_ptr,
                      SizeType rhs_ld, SizeType rhs_leading_nr_blocks,
                      DistributionType rhs_dist) noexcept;

    DistributedMatrix<ElType>* subMatrixInternal(const char* func, SizeType m, SizeType n,
                                                 Global2DIndex index) const;

    LocalMatrix<ElType>* localSubMatrixInternal(const char* func, SizeType m, SizeType n,
                                                Global2DIndex index) const;

    void reset();

    // Remove the copy-back step for const matrices.
    void remove_referenced();

    std::shared_ptr<memory::MemoryAllocator<ElementType>> ptr_;
    std::pair<SizeType, SizeType> size_;
    std::pair<SizeType, SizeType> block_size_;
    std::pair<SizeType, SizeType> local_size_;
    Global2DIndex base_index_;
    Local2DIndex local_base_index_;
    SizeType ld_;
    SizeType leading_nr_blocks_;
    DistributionType distribution_;
    const comm::Communicator2DGrid* comm_grid_;
    std::shared_ptr<memory::MemoryAllocator<ElementType>> referenced_ptr_ = nullptr;
    SizeType referenced_ld_;
    SizeType referenced_leading_nr_blocks_;
    DistributionType referenced_distribution_;
  };

  // ***********************************
  // *** Constructors Implementation ***
  // ***********************************
  //
  // Implementation of DistributedMatrix constructors is in the separate file
  // !!! Must include file IPP !!!
  //
  #include "distributed_matrix.ipp"

  namespace util {
    template <class ElType>
    std::ostream& operator<<(std::ostream& out, const DistributedMatrix<ElType>& mat) {
      mat.debugInfo(out);
      return out;
    }
  }
}

using dla_interface::util::operator<<;


#endif  // DLA_INTERFACE_DISTRIBUTED_MATRIX_H
