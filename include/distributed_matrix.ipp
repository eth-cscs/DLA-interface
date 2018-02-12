// public interface
#include "util_output.h"

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

template <class ElType>
DistributedMatrix<ElType>::DistributedMatrix(SizeType m, SizeType n, SizeType mb, SizeType nb,
                                             const comm::Communicator2DGrid& comm,
                                             DistributionType new_distribution, ElementType* ptr,
                                             std::size_t len, SizeType ld, SizeType leading_nr_blocks,
                                             DistributionType original_distribution)

 : DistributedMatrix(
       __func__, std::make_pair(m, n), Global2DIndex(0, 0), std::make_pair(mb, nb), &comm,
       new_distribution, std::make_shared<memory::MemoryAllocator<ElementType>>(ptr, len), ld,
       original_distribution == scalapack_dist ? 1 : leading_nr_blocks, original_distribution) {}

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
DistributedMatrix<ElType>::~DistributedMatrix() {
  if (referenced_ptr_) {
    copyInternal(true, referenced_ptr_, referenced_ld_, referenced_leading_nr_blocks_,
                 referenced_distribution_);
  }
}

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
bool DistributedMatrix<ElType>::isSameMatrix(const DistributedMatrix& rhs) const {
  if (this == &rhs)
    return true;

  if (ptr_->ptr() == rhs.ptr_->ptr() && size_ == rhs.size_ && block_size_ == rhs.block_size_ &&
      base_index_ == rhs.base_index_ && ld_ == rhs.ld_ &&
      leading_nr_blocks_ == rhs.leading_nr_blocks_ && distribution_ == rhs.distribution_ &&
      comm_grid_ == rhs.comm_grid_)
    return true;

  return false;
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
    throw std::invalid_argument(errorMessage(
        "Elements are on different nodes: Global index of ", Local2DIndex(0, 0), " are different: ",
        getGlobal2DIndex(Local2DIndex(0, 0)), " != ", rhs.getGlobal2DIndex(Local2DIndex(0, 0))));
  if (local_base_index_.row % block_size_.first != rhs.local_base_index_.row % rhs.block_size_.first ||
      local_base_index_.col % block_size_.second != rhs.local_base_index_.col % rhs.block_size_.second)
    throw std::invalid_argument(
        errorMessage("Elements are on different nodes: Local base index are: ", local_base_index_,
                     " and ", rhs.local_base_index_));

  if (distribution_ == scalapack_dist && rhs.distribution_ == scalapack_dist)
    util::memory::copy2D(local_size_, rhs.ptr(), rhs.ld_, ptr(), ld_);
  else {
    // Copy tile per tile
    int j_tile = 0;
    int n_tile = block_size_.second - local_base_index_.col % block_size_.second;
    while (j_tile < local_size_.second) {
      int i_tile = 0;
      int m_tile = block_size_.first - local_base_index_.row % block_size_.first;
      while (i_tile < local_size_.first) {
        Local2DIndex index(i_tile, j_tile);
        util::memory::copy2D(std::make_pair(m_tile, n_tile), rhs.ptr(index), rhs.ld_, ptr(index),
                             ld_);

        i_tile = i_tile + m_tile;
        m_tile = std::min(local_size_.first - i_tile, block_size_.first);
      }
      j_tile = j_tile + n_tile;
      n_tile = std::min(local_size_.second - j_tile, block_size_.second);
    }
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

  return std::make_tuple(ptr_->ptr(), i, j, std::move(desc));
}

template <class ElType>
std::tuple<const ElType*, IndexType, IndexType, std::array<int, 9>> DistributedMatrix<
    ElType>::getScalapackDescription() const {
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

  return std::make_tuple(ptr_->ptr(), i, j, std::move(desc));
}
#endif
#ifdef DLA_HAVE_DPLASMA
template <class ElType>
DPlasmaDescriptor DistributedMatrix<ElType>::getDPlasmaDescriptionInternal() const {
  if (distribution_ != tile_dist) {
    throw std::invalid_argument(errorMessage(
        "Cannot return DPlasma descriptor of this matrix. distribution() != tile_dist"));
  }
  if (ld_ != block_size_.first) {
    throw std::invalid_argument(errorMessage(
        "leadingDimension() ", ld_, " != mb (i.e. blockSize().first())", block_size_.first));
  }

  DPlasmaDescriptor desc;
  MPI_Comm row_ordered_comm = comm_grid_->rowOrderedMPICommunicator();
  int rank;
  int size;
  MPI_Comm_rank(row_ordered_comm, &rank);
  MPI_Comm_size(row_ordered_comm, &size);

  // The value of the following expression is different among the ranks (!!!).
  // However when used in DPlasma to compute the leading number of blocks (nb_elem_r) gives the
  // expected result on each rank.
  int m_full_estimate = std::max(
      base_index_.row + size_.first,
      ((leading_nr_blocks_ - 1) * comm_grid_->size2D().first + comm_grid_->id2D().first + 1) *
          block_size_.first);

  two_dim_block_cyclic_init(&desc, TypeInfo<ElType>::dplasma_type, matrix_Tile, size, rank,
                            block_size_.first, block_size_.second, m_full_estimate,
                            base_index_.col + size_.second, base_index_.row, base_index_.col,
                            size_.first, size_.second, 1, 1, comm_grid_->size2D().first);
  return desc;
}

template <class ElType>
std::tuple<DPlasmaDescriptor, MPI_Comm> DistributedMatrix<ElType>::getDPlasmaDescription() {
  DPlasmaDescriptor desc = getDPlasmaDescriptionInternal();
  MPI_Comm row_ordered_comm = comm_grid_->rowOrderedMPICommunicator();

  desc.mat = ptr_->ptr();

  return std::make_tuple(desc, row_ordered_comm);
}

template <class ElType>
std::tuple<const DPlasmaDescriptor, MPI_Comm> DistributedMatrix<ElType>::getDPlasmaDescription() const {
  DPlasmaDescriptor desc = getDPlasmaDescriptionInternal();
  MPI_Comm row_ordered_comm = comm_grid_->rowOrderedMPICommunicator();

  desc.mat = ptr_->ptr();

  return std::make_tuple(desc, row_ordered_comm);
}
#endif

template <class ElType>
template <class Out>
void DistributedMatrix<ElType>::debugDump(Out& out) {
  for (int j = 0; j < local_size_.second; ++j) {
    for (int i = 0; i < local_size_.first; ++i) {
      Local2DIndex ind(i, j);
      util::dumpDistributedMatrixElement(
          out, getGlobal2DIndex(ind), ind, getLocalStorageIndex(ind),
          storageBaseIndexFromLocalBaseIndex(__func__, ind + local_base_index_), operator()(ind));
    }
  }
  out << "\n";
}
// private interface

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

// if distribution == scalapack_distribution, leading_nr_blocks has to be 1.
template <class ElType>
DistributedMatrix<ElType>::DistributedMatrix(
    const char* func, std::pair<SizeType, SizeType> size, Global2DIndex base_index,
    std::pair<SizeType, SizeType> block_size, const comm::Communicator2DGrid* comm,
    DistributionType distribution,
    std::shared_ptr<memory::MemoryAllocator<ElementType>> original_ptr, SizeType original_ld,
    SizeType original_leading_nr_blocks, DistributionType original_distribution)
 : size_(size), block_size_(block_size), base_index_(base_index), ld_(original_ld),
   leading_nr_blocks_(original_leading_nr_blocks), distribution_(distribution), comm_grid_(comm) {
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
    checkAndComputeLocalParam(func, false, false);  // overwrites ld_ and ld_nr_blocks_.
    std::size_t original_len =
        allocationSize(func, original_ld, original_leading_nr_blocks, original_distribution);
    checkOrSetLeadingDims(func, original_ld, true, original_leading_nr_blocks, true,
                          original_distribution);
    if (original_len > original_ptr->size())
      throw std::invalid_argument(errorMessageFunc(func,
                                                   "Local matrix needs more elements: Required ",
                                                   original_len, " > given ", original_ptr->size()));
    referenced_ptr_ = original_ptr;
    referenced_ld_ = original_ld;
    referenced_leading_nr_blocks_ = original_leading_nr_blocks;
    referenced_distribution_ = original_distribution;

    ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(allocationSize(func));

    copyInternal(false, referenced_ptr_, referenced_ld_, referenced_leading_nr_blocks_,
                 referenced_distribution_);
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
  return storageBaseIndexFromLocalBaseIndex(func, index, block_size_, ld_, leading_nr_blocks_,
                                            distribution_);
}

template <class ElType>
std::size_t DistributedMatrix<ElType>::storageBaseIndexFromLocalBaseIndex(
    const char* func, const Local2DIndex& index, const std::pair<SizeType, SizeType>& block_size,
    SizeType ld, SizeType leading_nr_blocks, DistributionType dist) {
  switch (dist) {
    case scalapack_dist:
      return index.row + ld * index.col;
    case tile_dist: {
      IndexType i_blk = index.row / block_size.first;
      IndexType j_blk = index.col / block_size.second;
      std::size_t block_id = util::multSize(leading_nr_blocks, j_blk) + i_blk;
      std::size_t block_index = util::multSize(ld, block_size.second) * block_id;
      return block_index + util::multSize(ld, index.col % block_size.second) +
             index.row % block_size.first;
    }
    default:
      throw std::invalid_argument(errorMessageFunc(func, "Invalid distribution: ", dist));
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
  int rank_i = rank1DFromBaseGlobalIndex(index.row, block_size_.first, 0, comm_grid_->size2D().first);
  int rank_j =
      rank1DFromBaseGlobalIndex(index.col, block_size_.second, 0, comm_grid_->size2D().second);
  return std::make_pair(rank_i, rank_j);
}

template <class ElType>
IndexType DistributedMatrix<ElType>::index1DBaseGlobalFromBaseLocal(  //
    IndexType index, SizeType block_size, int rank_src, int rank, int comm_size) {
  int my_shifted_rank = (rank - rank_src + comm_size) % comm_size;
  return block_size * (index / block_size * comm_size + my_shifted_rank) + index % block_size;
}

template <class ElType>
int DistributedMatrix<ElType>::rank1DFromBaseGlobalIndex(IndexType index, SizeType block_size,
                                                         int rank_src, int comm_size) {
  return (index / block_size + rank_src) % comm_size;
}

template <class ElType>
IndexType DistributedMatrix<ElType>::index1DBaseLocalFromBaseGlobal(  //
    IndexType index, SizeType block_size, int rank_src, int rank, int comm_size) {
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

// Checks size_, block_size_, leading_nr_blocks_, ld_, distribution_, base_index_
// and sets local_size_, local_base_index_.
// The following operations are done before checking:
// If ld_set == false, ld_ is set.
// If ld_nr_bl_set == false, leading_nr_blocks_ is set.
template <class ElType>
void DistributedMatrix<ElType>::checkAndComputeLocalParam(const char* func, bool ld_set,
                                                          bool ld_nr_bl_set) {
  if (size_.first < 0 || size_.second < 0) {
    throw std::invalid_argument(errorMessageFunc(func, "Invalid matrix size ", size_));
  }

  if (!comm_grid_ && size_ != std::make_pair(0, 0)) {
    throw std::invalid_argument(
        errorMessageFunc(func, "Only (0, 0) matrices can have no communicator."));
  }

  if (block_size_.first <= 0 || block_size_.second <= 0) {
    throw std::invalid_argument(errorMessageFunc(func, "Invalid matrix block_size ", block_size_));
  }

  if (base_index_.row < 0 || base_index_.col < 0) {
    throw std::invalid_argument(errorMessageFunc(func, "Invalid base_index ", base_index_));
  }

  if (size_.first == 0 || size_.second == 0) {
    size_ = std::make_pair(0, 0);
    local_size_ = std::make_pair(0, 0);
    local_base_index_ = Local2DIndex(0, 0);
  }
  else {
    IndexType local_index_i = index1DBaseLocalFromBaseGlobal(
        base_index_.row, block_size_.first, 0, comm_grid_->id2D().first, comm_grid_->size2D().first);
    IndexType local_index_j =
        index1DBaseLocalFromBaseGlobal(base_index_.col, block_size_.second, 0,
                                       comm_grid_->id2D().second, comm_grid_->size2D().second);
    SizeType local_m =
        index1DBaseLocalFromBaseGlobal(size_.first + base_index_.row, block_size_.first, 0,
                                       comm_grid_->id2D().first, comm_grid_->size2D().first) -
        local_index_i;

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
  checkOrSetLeadingDims(func, ld_, ld_set, leading_nr_blocks_, ld_nr_bl_set, distribution_);
}

// Checks if ld and ld_nr_blks are correct for a matrix of local_size_, block_size_, local_index_
// and distribution.
// Throws if one the checks fails.
// The following operations are done before checking:
// If ld_set == false, ld is set.
// If ld_nr_bl_set == false, ld_nr_blks is set.
// Precondition: local_size_, block_size_ and local_base_index_ has to be correctly set.
template <class ElType>
void DistributedMatrix<ElType>::checkOrSetLeadingDims(const char* func, SizeType& ld, bool ld_set,
                                                      SizeType& ld_nr_blks, bool ld_nr_bl_set,
                                                      DistributionType distribution) {
  constexpr int chunk = 64;
  SizeType ld_min = 1;
  SizeType ld_val = 1;
  SizeType leading_nr_blocks_min = 1;

  SizeType local_total_m = 0;
  if (comm_grid_) {
    local_total_m =
        index1DBaseLocalFromBaseGlobal(size_.first + base_index_.row, block_size_.first, 0,
                                       comm_grid_->id2D().first, comm_grid_->size2D().first);
  }

  switch (distribution) {
    case scalapack_dist:
      ld_min = std::max(1, local_total_m);
      ld_val = ld_min == 1 ? 1 : util::ceilDiv(ld_min, chunk) * chunk;
      break;
    case tile_dist:
      ld_min = block_size_.first;
      ld_val = block_size_.first;
      leading_nr_blocks_min = std::max(1, util::ceilDiv(local_total_m, block_size_.first));
      break;
    default:
      throw std::invalid_argument(errorMessageFunc(func, "Invalid distribution: ", distribution));
  }
  if (!ld_set) {
    ld = ld_val;
  }
  if (!ld_nr_bl_set) {
    ld_nr_blks = leading_nr_blocks_min;
  }
  if (ld < ld_min) {
    throw std::invalid_argument(
        errorMessageFunc(func, "ld (", ld, " < ", ld_min, ") is too small."));
  }
  if (distribution == scalapack_dist && ld_nr_blks != 1) {
    throw std::invalid_argument(errorMessageFunc(func, "leading_nr_block is not 1."));
  }
  else if (ld_nr_blks < leading_nr_blocks_min) {
    throw std::invalid_argument(errorMessageFunc(func, "leading_nr_block (", ld_nr_blks, " < ",
                                                 leading_nr_blocks_min, ") is too small."));
  }
}

// Compute the minimum size of the allocated storage.
// Precondition: The global parameters and the local parameter has to be set, in particular:
// - ld_, leading_nr_blocks_, local_size_, local_base_index_, block_size_ and distribution has to
// be set correctly.
template <class ElType>
std::size_t DistributedMatrix<ElType>::allocationSize(const char* func) const {
  return allocationSize(func, ld_, leading_nr_blocks_, distribution_);
}

// Compute the minimum size of the allocated storage for the given leadind dimensions and
// distribution.
// Precondition: The global parameters and the local parameter has to be set, in particular:
// - local_size_, local_base_index_, block_size_ and distribution has to be set correctly.
// - ld, leading_nr_blocks has to be correct.
template <class ElType>
std::size_t DistributedMatrix<ElType>::allocationSize(  //
    const char* func, SizeType ld, SizeType leading_nr_blocks, DistributionType distribution) const {
  switch (distribution) {
    case scalapack_dist:
      return util::multSize(ld, local_base_index_.col + local_size_.second);
    case tile_dist:
      return util::multSize(ld, block_size_.second) *
             util::multSize(
                 leading_nr_blocks,
                 util::ceilDiv(local_base_index_.col + local_size_.second, block_size_.second));
    default:
      throw std::invalid_argument(errorMessageFunc(func, "Invalid distribution: ", distribution));
  }
  return 0;
}

// Copies element of rhs to this if copy_back == false.
// Copies element of this to rhs if copy_back == true.
template <class ElType>
void DistributedMatrix<ElType>::copyInternal(
    bool copy_back, std::shared_ptr<memory::MemoryAllocator<ElementType>> rhs_ptr, SizeType rhs_ld,
    SizeType rhs_leading_nr_blocks, DistributionType rhs_dist) noexcept {
  // Copy tile per tile
  int j_tile = 0;
  int n_tile = block_size_.second - local_base_index_.col % block_size_.second;
  while (j_tile < local_size_.second) {
    int i_tile = 0;
    int m_tile = block_size_.first - local_base_index_.row % block_size_.first;
    while (i_tile < local_size_.first) {
      Local2DIndex index(i_tile, j_tile);
      std::size_t rhs_index = storageBaseIndexFromLocalBaseIndex(
          __func__, index + local_base_index_, block_size_, rhs_ld, rhs_leading_nr_blocks, rhs_dist);
      if (copy_back)
        util::memory::copy2D(std::make_pair(m_tile, n_tile), ptr(index), ld_,
                             rhs_ptr->ptr(rhs_index), rhs_ld);
      else
        util::memory::copy2D(std::make_pair(m_tile, n_tile), rhs_ptr->ptr(rhs_index), rhs_ld,
                             ptr(index), ld_);

      i_tile = i_tile + m_tile;
      m_tile = std::min(local_size_.first - i_tile, block_size_.first);
    }
    j_tile = j_tile + n_tile;
    n_tile = std::min(local_size_.second - j_tile, block_size_.second);
  }
}

template <class ElType>
DistributedMatrix<ElType>* DistributedMatrix<ElType>::subMatrixInternal(  //
    const char* func, SizeType m, SizeType n, Global2DIndex index) const {
  checkGlobalIndex(func, false, index);

  if (m < 0 || n < 0)
    throw std::invalid_argument(errorMessageFunc(func, "Invalid submatrix size (", m, ", ", n, ")"));
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

template <class ElType>
LocalMatrix<ElType>* DistributedMatrix<ElType>::localSubMatrixInternal(  //
    const char* func, SizeType m, SizeType n, Global2DIndex index) const {
  checkGlobalIndex(func, false, index);

  if (m < 0 || n < 0)
    throw std::invalid_argument(errorMessageFunc(func, "Invalid submatrix size (", m, ", ", n, ")"));
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

template <class ElType>
void DistributedMatrix<ElType>::remove_referenced() {
  referenced_ptr_ = nullptr;
  referenced_ld_ = {};
  referenced_leading_nr_blocks_ = {};
  referenced_distribution_ = {};
}
