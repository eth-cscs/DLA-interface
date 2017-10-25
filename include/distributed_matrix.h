#ifndef DLA_INTERFACE_DISTRIBUTED_MATRIX_H
#define DLA_INTERFACE_DISTRIBUTED_MATRIX_H

#include <array>
#include <tuple>
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "matrix_index.h"
#include "local_matrix.h"
#include "scalapack.h"
#include "types.h"

namespace dla_interface {

  template <class ElType>
  class DistributedMatrix {
    public:
    using ElementType = ElType;

    // Creates a (0 x 0) matrix with blocksize (1 x 1).
    DistributedMatrix();

    // Creates a (m x n) matrix (or a (0 x 0) matrix if m == 0 or n == 0),
    // distributed on the rank grid described by comm according the given distribution
    // using the blocksize (mb x nb).
    // Throws a std::invalid_argument if m < 0, n < 0, mb < 1, nb < 1.
    DistributedMatrix(SizeType m, SizeType n, SizeType mb, SizeType nb,
                      const comm::Communicator2DGrid& comm, DistributionType distribution);

    // Creates a (m x n) matrix (or a (0 x 0) matrix if m == 0 or n == 0),
    // distributed on the rank grid described by comm according the given distribution
    // using the blocksize (mb x nb) and leading dimension ld.
    // Throws a std::invalid_argument
    //     - if m < 0, n < 0, mb < 1, nb < 1, or ld < 1
    //     - if ld < local_m and n > 0 and distribution == scalapack_dist
    //     - if ld < mb and distribution == tile_dist
    DistributedMatrix(SizeType m, SizeType n, SizeType mb, SizeType nb, SizeType ld,
                      const comm::Communicator2DGrid& comm, DistributionType distribution);

    // Creates a (m x n) matrix (or a (0 x 0) matrix if m == 0 or n == 0),
    // distributed on the rank grid described by comm according to the distribution
    // current_distribution using the blocksize (mb x nb).
    // This matrix uses the elements stored in the memory provided by ptr (of length len elements)
    // with leading dimension ld, and leading number of blocks leading_nr_blocks and
    // distributed according to original_distribution.
    // See note (1).
    // Throws a std::invalid_argument
    //     - if m < 0, n < 0, mb < 1, nb < 1, or ld < 1
    //     - if ld < local_m and n > 0 and distribution == scalapack_dist
    //     - if ld < mb and distribution == tile_dist
    //     - if leading_nr_blocks < ceil(local_m / mb) and distribution == tile_dist
    // Note: if original_distribution == scalapack_dist leading_nr_blocks is ignored.
    DistributedMatrix(SizeType m, SizeType n, SizeType mb, SizeType nb,
                      const comm::Communicator2DGrid& comm, DistributionType new_distribution,
                      ElementType* ptr, std::size_t len, SizeType ld, SizeType leading_nr_blocks,
                      DistributionType original_distribution);

    // Creates a distributed matrix which has the same size, block size and rank grid
    // of the given DistributedMatrix distributed according distribution
    // and which reference the same memory (See note (1)).
    DistributedMatrix(DistributionType distribution, DistributedMatrix& mat);
#ifdef DLA_HAVE_SCALAPACK
    // Creates a distributed matrix which has the same size, block size and rank grid
    // of the given Scalapack matrix and which reference the same memory (See note (1)).
    // Throws std::invalid_argument
    //     - if the scalapack context was not created via the CommunicatorManager interface.
    //     - if the scalapack descriptor does not describe a valid Scalapack matrix
    //     - if the indices i < 1, j < 1,
    //     - if the sizes m < 0, n < 0,
    //     - if the indices i + m - 1, j + n - 1 are larger than the corresponding matrix global
    //     sizes.
    //     - if desc[6] (rsrc) != 0 or desc[7] (csrc) != 0.
    DistributedMatrix(DistributionType distribution, int m, int n, ElementType* ptr,
                      ScalapackIndex i, ScalapackIndex j, ScalapackDescriptor desc);
#endif
#ifdef DLA_HAVE_DPLASMA
// TODO: check this and test it.
// Creates a distributed matrix which has the same size, block size and rank grid
// of the given D-PLASMA matrix and which reference the same memory (See note (1)).
// DistributedMatrix(DistributionType distribution, DPlasmaDescriptor desc);
#endif

    // (1) If the matrix distribution matches the distribution of the original matrix:
    //         - the matrix reference directly the elements given by the ptr,
    //         - the matrix element and the original matrix elements are the same,
    //         - the leading dimension and leading number of blocks remain unchanged.
    //     If the matrix distribution does not match the distribution of the original matrix:
    //         - new memory is allocated for the constructed matrix,
    //         - the leading dimension and leading number of blocks are not the same,
    //         - the elements of the original matrix are copied to the new matrix
    //           (the elements are redistributed),
    //         - On destruction of the DistributedMatrix object, the elements are copied back
    //           to the original matrix (redistributing them).
    //     Requirement: The lifetime of the original matrix has to be longer than the lifetime
    //     of the DistributedMatrix object.

    // Creates a copy of rhs, i.e. a distributed matrix with the same size, block size,
    // distribution, leading dimension, rank grid and a copy of the elements.
    DistributedMatrix(const DistributedMatrix& rhs);

    ~DistributedMatrix();

    // Creates a distributed matrix of the same size, block size, distribution,
    // leading dimension and rank grid moving the elements of rhs.
    // Postcondition: rhs is in an unspecified state.
    DistributedMatrix(DistributedMatrix&& rhs) noexcept;

    // Copy assignement operator:
    // Changes the size, block size, distribution, leading dimension and rank grid
    // of *this to match rhs values and copies the elements of rhs.
    // Note: - new memory is allocated,
    //       - old memory is deallocated if:
    //             - the memory was allocated via the DistributedMatrixClass constructors and
    //             - the memory is not used by any other Distributed/LocalMatrixClass objects.
    // All the pointers and references to elements of *this are invalidated.
    DistributedMatrix& operator=(const DistributedMatrix& rhs);

    // Move-assignement operator:
    // Changes the size, block size, distribution, leading dimension and rank grid
    // of *this to match rhs values and moves the elements of rhs.
    // All the pointers and references to elements of *this and rhs are invalidated.
    // Postcondition: rhs is in a unspecified state.
    DistributedMatrix& operator=(DistributedMatrix&& rhs) noexcept;

    // Copies the value of the elements of rhs.
    // Throws std::invalid_argument if *this and rhs do not have the same size.
    // If size() and rhs.size() != (0, 0):
    //   Throws std::invalid_argument if *this and rhs do not have the same block size, rank grid,
    //     element node distribution (i.e. if the (i, j)-th element of this and (i, j)-th element of
    //     rhs are on different nodes for any (i, j)).
    DistributedMatrix& copy(const DistributedMatrix& rhs);

    // Creates a DistributedMatrix which represent the (m, n) submatrix starting at global index (i,
    // j),
    // and returns a shared_ptr to it.
    // Throws std::invalid_argument if m, n, i, j < 0, i + m > size().first, j + n > size().second
    // See note (2).
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
    // Creates a LocalMatrix which represent the (m, n) submatrix starting at global index (i, j),
    // and returns a shared_ptr to it.
    // Throws std::invalid_argument
    //     - if m, n, i, j < 0, i + m > size().first, j + n > size().second
    //     - if (baseIndex().first + i) % mb + m > mb or (baseIndex().second + j) % nb + n > nb.
    // See note (2).
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

    // (2) The returned Local/DistributedMatrix will use part of the memory allocated for the
    // original matrix.
    //     If the original matrix is destroyed, the memory of the original matrix is not deallocated
    //     until all non zero size submatrices are destroyed.
    //     If elements are changed in the original matrix, the submatrix will change as well and
    //     vice-versa.
    //     If the original matrix is re-assigned using the assignement or move-assignement operator
    //     the submatrix object is not modified.
    //     Exception: if the memory has been provided by the user (see note (1)), the user has to
    //                guarantee that the memory is valid for the lifetime of the submatrix as well.

    Global2DIndex baseIndex() const {
      return base_index_;
    }
    std::pair<SizeType, SizeType> size() const {
      return size_;
    }
    Local2DIndex localBaseIndex() const {
      return local_base_index_;
    }
    std::pair<SizeType, SizeType> localSize() const {
      return local_size_;
    }
    std::pair<SizeType, SizeType> blockSize() const {
      return block_size_;
    }
    // Returns the number of elements allocated in the memory.
    // Note: Different matrices can share the same memory.
    std::size_t localStorageSize() const {
      return ptr_->size();
    }

    // Returns
    // - the leading dimension of the local matrix if the matrix is a scalapack matrix,
    // - the leading dimension of each tile if the matrix is a tile matrix otherwise.
    SizeType leadingDimension() const {
      return ld_;
    }
    // Returns
    // - 1 if the matrix is a scalapack matrix,
    // - the number of blocks stored in each column of tiles otherwise.
    SizeType leadingNumberOfBlocks() const {
      return leading_nr_blocks_;
    }

    DistributionType distribution() const {
      return distribution_;
    }

    // Precondition: the matrix is a non default constructed matrix.
    const comm::Communicator2DGrid& commGrid() const {
      assert(comm_grid_ != nullptr);
      return *comm_grid_;
    }

    // Returns the index of the position in memory of the global (i, j) element.
    // It holds ptr() + getLocalStorageIndex(index) == ptr(index).
    // Throws std::invalid_argument
    //            - if i < 0 or i > size().first,
    //            - if j < 0 or j > size().second,
    //            - if the element is not owned by this rank.
    std::size_t getLocalStorageIndex(IndexType i, IndexType j) const {
      return getLocalStorageIndex(Global2DIndex(i, j));
    }
    // Same as getLocalStorageIndex(IndexType i, IndexType j) with i = index.row and j = index.col.
    std::size_t getLocalStorageIndex(Global2DIndex index) const;

    // Returns the index of the position in memory of the local (index.row, index.col) element.
    // It holds ptr() + getLocalStorageIndex(index) == ptr(index).
    // Throws std::invalid_argument
    //            - if index.row < 0 or index.row > localSize().first,
    //            - if index.col < 0 or index.col > localSize().second.
    std::size_t getLocalStorageIndex(Local2DIndex index) const;

    // Returns the local 2D index of the global (i, j) element if owned by this rank.
    // Otherwise returns the local 2D index of the global (i1, j1) element,
    // where i1 is the smallest index > i and j1 is the smallest index > j
    // such that the global element (i1, j1) is owned by this rank.
    // Throws std::invalid_argument
    //            - if i < 0 or i > size().first,
    //            - if j < 0 or j > size().second.
    Local2DIndex getLocal2DIndex(IndexType i, IndexType j) const {
      return getLocal2DIndex(Global2DIndex(i, j));
    }
    // Same as getLocal2DIndex(IndexType i, IndexType j) with i = index.row and j = index.col.
    Local2DIndex getLocal2DIndex(Global2DIndex index) const;

    // Returns the 2D id of the rank which own the global (i, j) element.
    // Throws std::invalid_argument
    //            - if i < 0 or i > size().first,
    //            - if j < 0 or j > size().second.
    std::pair<int, int> getRankId2D(IndexType i, IndexType j) const {
      return getRankId2D(Global2DIndex(i, j));
    }
    // Same as getRankId2D(IndexType i, IndexType j) with i = index.row and j = index.col.
    std::pair<int, int> getRankId2D(Global2DIndex index) const;

    // Returns the global 2D index of the local (index.row, index.col) element.
    // Throws std::invalid_argument
    //            - if index.row < 0 or index.row > size().first,
    //            - if index.col < 0 or index.col > size().second.
    Global2DIndex getGlobal2DIndex(Local2DIndex index) const;

    // Returns a reference to the global (i, j) element.
    // Throws std::invalid_argument
    //            - if i < 0 or i > size().first,
    //            - if j < 0 or j > size().second,
    //            - if the element is not owned by this rank.
    const ElementType& operator()(IndexType i, IndexType j) const {
      return *ptr(i, j);
    }
    ElementType& operator()(IndexType i, IndexType j) {
      return *ptr(i, j);
    }
    // Same as operator()(IndexType i, IndexType j) with i = index.row and j = index.col.
    const ElementType& operator()(Global2DIndex index) const {
      return *ptr(index);
    }
    ElementType& operator()(Global2DIndex index) {
      return *ptr(index);
    }

    // Returns a reference to the local (index.row, index.col) element.
    // Throws std::invalid_argument
    //            - if index.row < 0 or index.row > localSize().first,
    //            - if index.col < 0 or index.col > localSize().second.
    const ElementType& operator()(Local2DIndex index) const {
      return *ptr(index);
    }
    ElementType& operator()(Local2DIndex index) {
      return *ptr(index);
    }

    // Returns nullptr if localSize().first == 0 and localSize.second == 0.
    // Returns the pointer to the first element allocated in the memory otherwise.
    // If ptr() != nullptr it holds ptr() + getLocalStorageIndex(index) == ptr(index),
    // where index is either a GlobalIndex or a LocaIndex.
    const ElementType* ptr() const;
    ElementType* ptr();

    // Returns the pointer to the global (i, j) element.
    // Throws std::invalid_argument
    //            - if i < 0 or i > size().first,
    //            - if j < 0 or j > size().second.
    //            - if the element is not owned by this rank.
    const ElementType* ptr(IndexType i, IndexType j) const {
      return ptr(Global2DIndex(i, j));
    }
    ElementType* ptr(IndexType i, IndexType j) {
      return ptr(Global2DIndex(i, j));
    }
    // Same as ptr(IndexType i, IndexType j) with i = index.row and j = index.col.
    const ElementType* ptr(Global2DIndex index) const;
    ElementType* ptr(Global2DIndex index);

    // Returns the pointer to the local (index.row, index.col) element.
    // Throws std::invalid_argument
    //            - if index,row < 0 or index.row > localSize().first,
    //            - if index.col < 0 or index.col > localSize().second.
    const ElementType* ptr(Local2DIndex index) const;
    ElementType* ptr(Local2DIndex index);

#ifdef DLA_HAVE_SCALAPACK
    // Returns a tuple containing ptr, i, j and a std::array containing the descriptor, needed for
    // ScaLAPACK calls.
    // Throws std::invalid_argument if the matrix is not distributed in the ScaLAPACK way.
    std::tuple<ElementType*, IndexType, IndexType, std::array<int, 9>> getScalapackDescription();
#endif
#ifdef DLA_HAVE_DPLASMA
    // Returns the DPLASMA descriptor.
    DPlasmaDescriptor getDPlasmaDescription();
#endif

    template <class Out>
    void debugDump(Out& out);

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

#include "distributed_matrix.ipp"
}

#endif  // DLA_INTERFACE_DISTRIBUTED_MATRIX_H
