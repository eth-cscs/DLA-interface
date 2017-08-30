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
    //     - if ld < m and n > 0 and distribution == scalapack_dist
    //     - if ld < mb and distribution == tile_dist
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
    std::size_t allocationSize(const char* func);

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
  };

  template <class ElType>
  DistributedMatrix<ElType>::DistributedMatrix()
   : DistributedMatrix(__func__, std::make_pair(0, 0), std::make_pair(1, 1), false, 0, nullptr,
                       scalapack_dist) {}

  template <class ElType>
  DistributedMatrix<ElType>::DistributedMatrix(SizeType m, SizeType n, SizeType mb, SizeType nb,
                                               const comm::Communicator2DGrid& comm,
                                               DistributionType distribution)
   : DistributedMatrix(__func__, std::make_pair(m, n), std::make_pair(mb, nb), false, 0, &comm,
                       distribution) {}

  template <class ElType>
  DistributedMatrix<ElType>::DistributedMatrix(SizeType m, SizeType n, SizeType mb, SizeType nb,
                                               SizeType ld, const comm::Communicator2DGrid& comm,
                                               DistributionType distribution)
   : DistributedMatrix(__func__, std::make_pair(m, n), std::make_pair(mb, nb), true, ld, &comm,
                       distribution) {}

  // Creates a (m x n) matrix (or a (0 x 0) matrix if m==0 or n == 0),
  // distributed on the rank grid described by comm according to the distribution
  // current_distribution
  // using the blocksize (mb x nb) and leading dimension ld.
  // This matrix uses the elements stored in the memory provided using ptr
  // (at least len elements allocated) and distributed according to original_distribution.
  // See note (1).
  // Throws a std::invalid_argument
  //     - if m < 0, n < 0, mb < 1, nb < 1, or ld < 1
  //     - if ld < m and n > 0 and distribution == scalapack_dist
  //     - if ld < mb and distribution == tile_dist
  //     - if len is smaller than the local elements needed to store the matrix.
  template <class ElType>
  DistributedMatrix<ElType>::DistributedMatrix(SizeType m, SizeType n, SizeType mb, SizeType nb,
                                               const comm::Communicator2DGrid& comm,
                                               DistributionType new_distribution, ElementType* ptr,
                                               std::size_t len, SizeType ld,
                                               SizeType leading_nr_blocks,
                                               DistributionType original_distribution)

   : DistributedMatrix(__func__, std::make_pair(m, n), Global2DIndex(0, 0), std::make_pair(mb, nb),
                       &comm, new_distribution,
                       std::make_shared<memory::MemoryAllocator<ElementType>>(ptr, len), ld,
                       leading_nr_blocks, original_distribution) {}

  // Creates a distributed matrix which has the same size, block size and rank grid
  // of the given DistributedMatrix and which reference the same memory (See note (1)).
  template <class ElType>
  DistributedMatrix<ElType>::DistributedMatrix(DistributionType distribution, DistributedMatrix& mat)
   : DistributedMatrix(__func__, mat.size_, mat.base_index_, mat.block_size_, mat.comm_grid_,
                       distribution, mat.ptr_, mat.ld_, mat.leading_nr_blocks_, mat.distribution_) {}
#ifdef DLA_HAVE_SCALAPACK
  template <class ElType>
  DistributedMatrix<ElType>::DistributedMatrix(DistributionType distribution, int m, int n,
                                               ElementType* ptr, ScalapackIndex i, ScalapackIndex j,
                                               ScalapackDescriptor desc)
      // use comma operator to check {index, size and total size}, rsrc and csrc.
      : DistributedMatrix(
            (({
               if (m + i - 1 > desc[2] && n > 0)
                 throw std::invalid_argument(errorMessage(" row size + row index (", m + i - 1,
                                                          ") > total row size (", desc[2], ")."));
               if (n + j - 1 > desc[3] && m > 0)
                 throw std::invalid_argument(errorMessage(" col size + col index (", n + j - 1,
                                                          ") > total col size (", desc[3], ")."));
               if (desc[6] != 0 || desc[7] != 0)
                 throw std::invalid_argument(errorMessage("Not supported: rsrc (", desc[6],
                                                          ") != 0 or csrc (", desc[7], ") != 0."));
             }),
             __func__),
            std::make_pair(m, n), Global2DIndex(i - 1, j - 1), std::make_pair(desc[4], desc[5]),
            &comm::CommunicatorManager::getCommunicator2DGridFromBlacsContext(desc[1]), distribution,
            std::make_shared<memory::MemoryAllocator<ElementType>>(
                ptr, util::multSize(
                         desc[8],
                         index1DBaseLocalFromBaseGlobal(
                             desc[3], desc[5], desc[7],
                             comm::CommunicatorManager::getCommunicator2DGridFromBlacsContext(desc[1])
                                 .id2D()
                                 .second,
                             comm::CommunicatorManager::getCommunicator2DGridFromBlacsContext(desc[1])
                                 .size2D()
                                 .second))),
            desc[8], 1, scalapack_dist) {}
#endif

  template <class ElType>
  DistributedMatrix<ElType>::DistributedMatrix(const DistributedMatrix& rhs)
   : size_(rhs.size_), block_size_(rhs.block_size_), local_size_(rhs.local_size_),
     base_index_(rhs.base_index_), local_base_index_(rhs.local_base_index_), ld_(rhs.ld_),
     leading_nr_blocks_(rhs.leading_nr_blocks_), distribution_(rhs.distribution_),
     comm_grid_(rhs.comm_grid_) {
    // TODO: can local_base_index_ be optimized?
    ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(allocationSize(__func__));
    copy(rhs);
  }

  template <class ElType>
  DistributedMatrix<ElType>::DistributedMatrix(DistributedMatrix&& rhs) noexcept
   : ptr_(rhs.ptr_), size_(rhs.size_), block_size_(rhs.block_size_), local_size_(rhs.local_size_),
     base_index_(rhs.base_index_), local_base_index_(rhs.local_base_index_), ld_(rhs.ld_),
     leading_nr_blocks_(rhs.leading_nr_blocks_), distribution_(rhs.distribution_),
     comm_grid_(rhs.comm_grid_) {
    rhs.reset();
  }

  template <class ElType>
  DistributedMatrix<ElType>& DistributedMatrix<ElType>::operator=(const DistributedMatrix& rhs) {
    if (this != &rhs) {
      size_ = rhs.size_;
      block_size_ = rhs.block_size_;
      local_size_ = rhs.local_size_;
      base_index_ = rhs.base_index_;
      // TODO: can local_base_index_ be optimized?
      local_base_index_ = rhs.local_base_index_;
      ld_ = rhs.ld_;
      leading_nr_blocks_ = rhs.leading_nr_blocks_;
      distribution_ = rhs.distribution_;
      comm_grid_ = rhs.comm_grid_;

      ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(allocationSize(__func__));
      copy(rhs);
    }
    return *this;
  }

  template <class ElType>
  DistributedMatrix<ElType>& DistributedMatrix<ElType>::operator=(DistributedMatrix&& rhs) noexcept {
    ptr_ = rhs.ptr_;
    size_ = rhs.size_;
    block_size_ = rhs.block_size_;
    local_size_ = rhs.local_size_;
    base_index_ = rhs.base_index_;
    local_base_index_ = rhs.local_base_index_;
    ld_ = rhs.ld_;
    leading_nr_blocks_ = rhs.leading_nr_blocks_;
    distribution_ = rhs.distribution_;
    comm_grid_ = rhs.comm_grid_;

    rhs.reset();
    return *this;
  }

  template <class ElType>
  DistributedMatrix<ElType>& DistributedMatrix<ElType>::copy(const DistributedMatrix& rhs) {
    if (size_ != rhs.size_)
      throw std::invalid_argument(errorMessage("Sizes do not match: ", size_, " != ", rhs.size_));
    if (size_ == std::make_pair(0, 0))
      return *this;
    if (comm_grid_ != rhs.comm_grid_)
      throw std::invalid_argument(errorMessage("Matrix has different communicator grids."));
    if (block_size_ != rhs.block_size_)
      throw std::invalid_argument(
          errorMessage("Block sizes do not match: ", block_size_, " != ", rhs.block_size_));
    if (local_size_ != rhs.local_size_)
      throw std::invalid_argument(
          errorMessage("Local sizes do not match: ", block_size_, " != ", rhs.block_size_));
    if (getRankId2D(0, 0) != rhs.getRankId2D(0, 0))
      throw std::invalid_argument(errorMessage("Element (0, 0) is on different ranks: ",
                                               getRankId2D(0, 0), " != ", rhs.getRankId2D(0, 0)));
    if (getGlobal2DIndex(Local2DIndex(0, 0)) != rhs.getGlobal2DIndex(Local2DIndex(0, 0)))
      throw std::invalid_argument(errorMessage("Elements are on different nodes: Global index of ",
                                               Local2DIndex(0, 0), " are different: ",
                                               getGlobal2DIndex(Local2DIndex(0, 0)), " != ",
                                               rhs.getGlobal2DIndex(Local2DIndex(0, 0))));
    if (local_base_index_.row % block_size_.first !=
            rhs.local_base_index_.row % rhs.block_size_.first ||
        local_base_index_.col % block_size_.second !=
            rhs.local_base_index_.col % rhs.block_size_.second)
      throw std::invalid_argument(
          errorMessage("Elements are on different nodes: Local base index are: ", local_base_index_,
                       " and ", rhs.local_base_index_));

    if (distribution_ != rhs.distribution_)
      throw std::invalid_argument(errorMessage("Not yet implemented"));

    switch (distribution_) {
      case scalapack_dist:
        util::memory::copy2D(local_size_, rhs.ptr(), rhs.ld_, ptr(), ld_);
        break;
      case tile_dist: {
        // Copy tile per tile
        int j_tile = 0;
        int n_tile = block_size_.second - local_base_index_.col % block_size_.second;
        while (j_tile < local_size_.second) {
          int i_tile = 0;
          int m_tile = block_size_.first - local_base_index_.row % block_size_.first;
          while (i_tile < local_size_.first) {
            Local2DIndex index(i_tile, j_tile);
            util::memory::copy2D(std::make_pair(m_tile, n_tile), rhs.ptr(index), rhs.ld_,
                                 ptr(index), ld_);

            i_tile = i_tile + m_tile;
            m_tile = std::min(local_size_.first - i_tile, block_size_.first);
          }
          j_tile = j_tile + n_tile;
          n_tile = std::min(local_size_.second - j_tile, block_size_.second);
        }
        break;
      }
      default:
        throw std::invalid_argument(
            errorMessageFunc(__func__, "Invalid distribution: ", distribution_));
    }

    return *this;
  }

  template <class ElType>
  std::size_t DistributedMatrix<ElType>::getLocalStorageIndex(Global2DIndex index) const {
    checkGlobalIndex(__func__, true, index);
    return storageBaseIndexFromLocalBaseIndex(
               __func__, localBaseIndexFromGlobalBaseIndex(__func__, index + base_index_)) -
           storageBaseIndexFromLocalBaseIndex(__func__, local_base_index_);
  }

  template <class ElType>
  std::size_t DistributedMatrix<ElType>::getLocalStorageIndex(Local2DIndex index) const {
    checkLocalIndex(__func__, index);
    return storageBaseIndexFromLocalBaseIndex(__func__, index + local_base_index_) -
           storageBaseIndexFromLocalBaseIndex(__func__, local_base_index_);
  }

  template <class ElType>
  Local2DIndex DistributedMatrix<ElType>::getLocal2DIndex(Global2DIndex index) const {
    checkGlobalIndex(__func__, false, index);
    return localBaseIndexFromGlobalBaseIndex(__func__, index + base_index_) - local_base_index_;
  }

  template <class ElType>
  std::pair<int, int> DistributedMatrix<ElType>::getRankId2D(Global2DIndex index) const {
    checkGlobalIndex(__func__, false, index);
    return rankFromBaseGlobalIndex(__func__, index + base_index_);
  }

  template <class ElType>
  Global2DIndex DistributedMatrix<ElType>::getGlobal2DIndex(Local2DIndex index) const {
    checkLocalIndex(__func__, index);
    return globalBaseIndexFromLocalBaseIndex(__func__, index + local_base_index_) - base_index_;
  }

  template <class ElType>
  const ElType* DistributedMatrix<ElType>::ptr() const {
    if (local_size_.first == 0)
      return nullptr;
    else
      return ptr_->ptr(storageBaseIndexFromLocalBaseIndex(__func__, local_base_index_));
  }
  template <class ElType>
  ElType* DistributedMatrix<ElType>::ptr() {
    if (local_size_.first == 0)
      return nullptr;
    else
      return ptr_->ptr(storageBaseIndexFromLocalBaseIndex(__func__, local_base_index_));
  }

  template <class ElType>
  const ElType* DistributedMatrix<ElType>::ptr(Global2DIndex index) const {
    checkGlobalIndex(__func__, true, index);
    return ptr_->ptr(storageBaseIndexFromGlobalBaseIndex(__func__, index + base_index_));
  }
  template <class ElType>
  ElType* DistributedMatrix<ElType>::ptr(Global2DIndex index) {
    checkGlobalIndex(__func__, true, index);
    return ptr_->ptr(storageBaseIndexFromGlobalBaseIndex(__func__, index + base_index_));
  }

  template <class ElType>
  const ElType* DistributedMatrix<ElType>::ptr(Local2DIndex index) const {
    checkLocalIndex(__func__, index);
    return ptr_->ptr(storageBaseIndexFromLocalBaseIndex(__func__, index + local_base_index_));
  }
  template <class ElType>
  ElType* DistributedMatrix<ElType>::ptr(Local2DIndex index) {
    checkLocalIndex(__func__, index);
    return ptr_->ptr(storageBaseIndexFromLocalBaseIndex(__func__, index + local_base_index_));
  }

#ifdef DLA_HAVE_SCALAPACK
  template <class ElType>
  std::tuple<ElType*, IndexType, IndexType, std::array<int, 9>> DistributedMatrix<
      ElType>::getScalapackDescription() {
    if (distribution_ != scalapack_dist) {
      throw std::invalid_argument(errorMessage(
          "Cannot return ScaLAPACK descriptor of this matrix. distribution() != scalapack_dist"));
    }
    int total_m = size_.first + base_index_.row;
    int total_n = size_.second + base_index_.col;
    int i = base_index_.row + 1;
    int j = base_index_.col + 1;
    auto ictxt = comm_grid_->blacsContext();
    std::array<int, 9> desc;
    int info;

    scalapack::descinit_(&desc[0], &total_m, &total_n, &block_size_.first, &block_size_.second,
                         DLA_I_ZERO, DLA_I_ZERO, &ictxt, &ld_, &info);
    if (info != 0) {
      throw std::runtime_error(
          errorMessage("Internal error: Wrong argument in descriptor initialization: ", info));
    }

    return std::make_tuple(ptr_->ptr(), i, j, desc);
  }
#endif

  template <class ElType>
  DistributedMatrix<ElType>::DistributedMatrix(const char* func, std::pair<SizeType, SizeType> size,
                                               std::pair<SizeType, SizeType> block_size, bool ld_set,
                                               SizeType ld, const comm::Communicator2DGrid* comm,
                                               DistributionType distribution)
   : size_(size), block_size_(block_size), base_index_(0, 0), ld_(ld), distribution_(distribution),
     comm_grid_(comm) {
    checkAndComputeLocalParam(func, ld_set, false);
    ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(allocationSize(func));
  }

  // if distribution == scalapack_distribution, leading_nr_blocks is not used.
  template <class ElType>
  DistributedMatrix<ElType>::DistributedMatrix(
      const char* func, std::pair<SizeType, SizeType> size, Global2DIndex base_index,
      std::pair<SizeType, SizeType> block_size, const comm::Communicator2DGrid* comm,
      DistributionType distribution,
      std::shared_ptr<memory::MemoryAllocator<ElementType>> original_ptr, SizeType original_ld,
      SizeType original_leading_nr_blocks, DistributionType original_distribution)
   : size_(size), block_size_(block_size), base_index_(base_index), ld_(original_ld),
     leading_nr_blocks_(distribution == scalapack_dist ? 1 : original_leading_nr_blocks),
     distribution_(distribution), comm_grid_(comm) {
    if (distribution_ == original_distribution) {
      checkAndComputeLocalParam(func, true, true);
      std::size_t len = allocationSize(func);
      if (len > original_ptr->size())
        throw std::invalid_argument(errorMessageFunc(func,
                                                     "Local matrix needs more elements: Required ",
                                                     len, " > given ", original_ptr->size()));
      ptr_ = original_ptr;
    }
    else {
      throw std::invalid_argument(
          errorMessageFunc(func, "Not implemented yet, change of distribution."));
    }
  }

  template <class ElType>
  void DistributedMatrix<ElType>::checkLocalIndex(const char* func, const Local2DIndex& index) const {
    if (index.row < 0 || index.row >= local_size_.first)
      throw std::invalid_argument(errorMessageFunc(func, "Row local index ", index.row,
                                                   " out of range [0, ", local_size_.first, ")"));
    if (index.col < 0 || index.col >= local_size_.second)
      throw std::invalid_argument(errorMessageFunc(func, "Column local index ", index.col,
                                                   " out of range [0, ", local_size_.second, ")"));
  }

  template <class ElType>
  void DistributedMatrix<ElType>::checkGlobalIndex(const char* func, bool check_rank,
                                                   const Global2DIndex& index) const {
    if (index.row < 0 || index.row >= size_.first)
      throw std::invalid_argument(errorMessageFunc(func, "Row global index ", index.row,
                                                   " out of range [0, ", size_.first, ")"));
    if (index.col < 0 || index.col >= size_.second)
      throw std::invalid_argument(errorMessageFunc(func, "Column global index ", index.col,
                                                   " out of range [0, ", size_.second, ")"));
    if (check_rank && comm_grid_->id2D() != getRankId2D(index))
      throw std::invalid_argument(
          errorMessageFunc(func, "Global index ", index, " is not on rank ", comm_grid_->id2D()));
  }

  template <class ElType>
  Global2DIndex DistributedMatrix<ElType>::globalBaseIndexFromLocalBaseIndex(
      const char* func, const Local2DIndex& index) const {
    IndexType index_i = index1DBaseGlobalFromBaseLocal(
        index.row, block_size_.first, 0, comm_grid_->id2D().first, comm_grid_->size2D().first);
    IndexType index_j = index1DBaseGlobalFromBaseLocal(
        index.col, block_size_.second, 0, comm_grid_->id2D().second, comm_grid_->size2D().second);
    return Global2DIndex(index_i, index_j);
  }

  template <class ElType>
  std::size_t DistributedMatrix<ElType>::storageBaseIndexFromLocalBaseIndex(
      const char* func, const Local2DIndex& index) const {
    switch (distribution_) {
      case scalapack_dist:
        return index.row + ld_ * index.col;
      case tile_dist: {
        IndexType i_blk = index.row / block_size_.first;
        IndexType j_blk = index.col / block_size_.second;
        std::size_t block_id = util::multSize(leading_nr_blocks_, j_blk) + i_blk;
        std::size_t block_index = util::multSize(ld_, block_size_.second) * block_id;
        return block_index + util::multSize(ld_, index.col % block_size_.second) +
               index.row % block_size_.first;
      }
      default:
        throw std::invalid_argument(errorMessageFunc(func, "Invalid distribution: ", distribution_));
    }
    return 0;
  }

  template <class ElType>
  Local2DIndex DistributedMatrix<ElType>::localBaseIndexFromGlobalBaseIndex(
      const char* func, const Global2DIndex& index) const {
    IndexType local_index_i = index1DBaseLocalFromBaseGlobal(
        index.row, block_size_.first, 0, comm_grid_->id2D().first, comm_grid_->size2D().first);
    IndexType local_index_j = index1DBaseLocalFromBaseGlobal(
        index.col, block_size_.second, 0, comm_grid_->id2D().second, comm_grid_->size2D().second);
    return Local2DIndex(local_index_i, local_index_j);
  }

  template <class ElType>
  std::size_t DistributedMatrix<ElType>::storageBaseIndexFromGlobalBaseIndex(
      const char* func, const Global2DIndex& index) const {
    return storageBaseIndexFromLocalBaseIndex(func, localBaseIndexFromGlobalBaseIndex(func, index));
  }

  template <class ElType>
  std::pair<int, int> DistributedMatrix<ElType>::rankFromBaseGlobalIndex(
      const char* func, const Global2DIndex& index) const {
    int rank_i =
        rank1DFromBaseGlobalIndex(index.row, block_size_.first, 0, comm_grid_->size2D().first);
    int rank_j =
        rank1DFromBaseGlobalIndex(index.col, block_size_.second, 0, comm_grid_->size2D().second);
    return std::make_pair(rank_i, rank_j);
  }

  // Checks size_, block_size_, leading_nr_blocks_, ld_, distribution_, base_index_
  // and sets local_size_, local_base_index_.
  // The following operations are done before checking:
  // If ld_set == false, ld_ is set.
  // If ld_nr_bl_set == false, leading_nr_blocks_ is set.
  template <class ElType>
  void DistributedMatrix<ElType>::checkAndComputeLocalParam(const char* func, bool ld_set,
                                                            bool ld_nr_bl_set) {
    constexpr int chunk = 64;
    if (size_.first < 0 || size_.second < 0) {
      throw std::invalid_argument(errorMessageFunc(func, "Invalid matrix size ", size_));
    }

    if (block_size_.first <= 0 || block_size_.second <= 0) {
      throw std::invalid_argument(errorMessageFunc(func, "Invalid matrix block_size ", block_size_));
    }

    if (base_index_.row < 0 || base_index_.col < 0) {
      throw std::invalid_argument(errorMessageFunc(func, "Invalid base_index ", base_index_));
    }

    // Needed to compute the minimum of ld or of ld_nr_blocks.
    int loc_m = 0;
    if (size_.first == 0 || size_.second == 0) {
      size_ = std::make_pair(0, 0);
      local_size_ = std::make_pair(0, 0);
      local_base_index_ = Local2DIndex(0, 0);
    }
    else {
      IndexType local_index_i =
          index1DBaseLocalFromBaseGlobal(base_index_.row, block_size_.first, 0,
                                         comm_grid_->id2D().first, comm_grid_->size2D().first);
      IndexType local_index_j =
          index1DBaseLocalFromBaseGlobal(base_index_.col, block_size_.second, 0,
                                         comm_grid_->id2D().second, comm_grid_->size2D().second);
      SizeType local_m =
          index1DBaseLocalFromBaseGlobal(size_.first + base_index_.row, block_size_.first, 0,
                                         comm_grid_->id2D().first, comm_grid_->size2D().first) -
          local_index_i;
      loc_m = static_cast<int>(local_m);
      SizeType local_n =
          index1DBaseLocalFromBaseGlobal(size_.second + base_index_.col, block_size_.second, 0,
                                         comm_grid_->id2D().second, comm_grid_->size2D().second) -
          local_index_j;
      if (local_m == 0 || local_n == 0) {
        local_base_index_ = Local2DIndex(0, 0);
        local_size_ = std::make_pair(0, 0);
      }
      else {
        local_base_index_ = Local2DIndex(local_index_i, local_index_j);
        local_size_ = std::make_pair(local_m, local_n);
      }
    }

    int ld_min = 1;
    int leading_nr_blocks_min = 1;
    switch (distribution_) {
      case scalapack_dist:
        ld_min = std::max(1, loc_m);
        break;
      case tile_dist:
        ld_min = block_size_.first;
        leading_nr_blocks_min = std::max(1, util::ceilDiv(loc_m, block_size_.first));
        break;
      default:
        throw std::invalid_argument(errorMessageFunc(func, "Invalid distribution: ", distribution_));
    }
    if (!ld_set) {
      ld_ = ld_min == 1 ? 1 : util::ceilDiv(ld_min, chunk) * chunk;
    }
    if (!ld_nr_bl_set) {
      leading_nr_blocks_ = leading_nr_blocks_min;
    }
    if (ld_ < ld_min) {
      throw std::invalid_argument(
          errorMessageFunc(func, "ld (", ld_, " < ", ld_min, ") is too small."));
    }
    if (distribution_ == scalapack_dist && leading_nr_blocks_ != 1) {
      throw std::invalid_argument(errorMessageFunc(func, "leading_nr_block is not 1."));
    }
    if (leading_nr_blocks_ < leading_nr_blocks_min) {
      throw std::invalid_argument(errorMessageFunc(func, "leading_nr_block (", leading_nr_blocks_,
                                                   " < ", leading_nr_blocks_min,
                                                   ") is too small."));
    }
  }

  // Compute the minimum size of the allocated storage.
  // Precondition: The global parameters and the local parameter has to be set, in particular:
  // - ld_, leading_nr_blocks_, local_size_, local_base_index_, block_size_ and distribution has to
  // be set correctly.
  template <class ElType>
  std::size_t DistributedMatrix<ElType>::allocationSize(const char* func) {
    switch (distribution_) {
      case scalapack_dist:
        return util::multSize(ld_, local_base_index_.col + local_size_.second);
      case tile_dist:
        return util::multSize(ld_, block_size_.second) *
               util::multSize(
                   leading_nr_blocks_,
                   util::ceilDiv(local_base_index_.col + local_size_.second, block_size_.second));
      default:
        throw std::invalid_argument(errorMessageFunc(func, "Invalid distribution: ", distribution_));
    }
    return 0;
  }

  template <class ElType>
  IndexType DistributedMatrix<ElType>::index1DBaseGlobalFromBaseLocal(IndexType index,
                                                                      SizeType block_size,
                                                                      int rank_src, int rank,
                                                                      int comm_size) {
    int my_shifted_rank = (rank - rank_src + comm_size) % comm_size;
    return block_size * (index / block_size * comm_size + my_shifted_rank) + index % block_size;
  }

  template <class ElType>
  int DistributedMatrix<ElType>::rank1DFromBaseGlobalIndex(IndexType index, SizeType block_size,
                                                           int rank_src, int comm_size) {
    return (index / block_size + rank_src) % comm_size;
  }

  template <class ElType>
  IndexType DistributedMatrix<ElType>::index1DBaseLocalFromBaseGlobal(IndexType index,
                                                                      SizeType block_size,
                                                                      int rank_src, int rank,
                                                                      int comm_size) {
    int my_shifted_rank = (rank - rank_src + comm_size) % comm_size;
    int block_id = index / block_size;
    int shifted_rank = block_id % comm_size;
    int local_block_id = block_id / comm_size;
    int local_index = local_block_id * block_size;
    if (my_shifted_rank == shifted_rank)
      return local_index + index % block_size;
    else if (my_shifted_rank < shifted_rank)
      return local_index + block_size;
    else
      return local_index;
  }

  template <class ElType>
  DistributedMatrix<ElType>* DistributedMatrix<ElType>::subMatrixInternal(  //
      const char* func, SizeType m, SizeType n, Global2DIndex index) const {
    checkGlobalIndex(func, false, index);

    if (m < 0 || n < 0)
      throw std::invalid_argument(
          errorMessageFunc(func, "Invalid submatrix size (", m, ", ", n, ")"));
    if (index.row + m > size_.first || index.col + n > size_.second)
      throw std::invalid_argument(errorMessageFunc(func, "Submatrix is exceeding matrix borders (",
                                                   index.row + m, ", ", index.col + n,
                                                   ") larger than ", size_));
    if (m == 0 || n == 0)
      return new DistributedMatrix<ElementType>(func, std::make_pair(0, 0), block_size_, true, ld_,
                                                comm_grid_, distribution_);

    return new DistributedMatrix<ElementType>(func, std::make_pair(m, n), index, block_size_,
                                              comm_grid_, distribution_, ptr_, ld_,
                                              leading_nr_blocks_, distribution_);
  }

  // Creates a LocalMatrix which represent the (m, n) submatrix starting at position (i ,j),
  // and returns a shared_ptr to it.
  // Throws std::invalid_argument
  //     - if m, n, i, j < 0, i + m > size().first, j + n > size().second
  //     - if (baseIndex.first() + i) % mb + m > mb or (baseIndex.second() + j) % nb + n > nb.
  // See note (2).
  template <class ElType>
  LocalMatrix<ElType>* DistributedMatrix<ElType>::localSubMatrixInternal(  //
      const char* func, SizeType m, SizeType n, Global2DIndex index) const {
    checkGlobalIndex(func, false, index);

    if (m < 0 || n < 0)
      throw std::invalid_argument(
          errorMessageFunc(func, "Invalid submatrix size (", m, ", ", n, ")"));
    if (index.row + m > size_.first || index.col + n > size_.second)
      throw std::invalid_argument(errorMessageFunc(func, "Submatrix is exceeding matrix borders (",
                                                   index.row + m, ", ", index.col + n,
                                                   ") larger than ", size_));
    if ((base_index_.row + index.row) % block_size_.first + m > block_size_.first ||
        (base_index_.col + index.col) % block_size_.second + n > block_size_.second)
      throw std::invalid_argument(
          errorMessageFunc(func, "Submatrix is exceeding matrix block borders."));
    if (comm_grid_->id2D() != getRankId2D(index) || m == 0 || n == 0)
      return new LocalMatrix<ElementType>(0, 0, ld_);

    return new LocalMatrix<ElementType>(
        func, m, n, ld_, ptr_, storageBaseIndexFromGlobalBaseIndex(__func__, index + base_index_));
  }

  template <class ElType>
  void DistributedMatrix<ElType>::reset() {
    ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(0);
    size_ = std::make_pair(0, 0);
    // block_size_ not changed
    local_size_ = std::make_pair(0, 0);
    base_index_ = Global2DIndex(0, 0);
    local_base_index_ = Local2DIndex(0, 0);
    ld_ = 1;
    leading_nr_blocks_ = 1;
    // distribution_ not changed
    // comm_grid not changed
  }
}

#endif  // DLA_INTERFACE_DISTRIBUTED_MATRIX_H
