#ifndef DLA_INTERFACE_LOCAL_MATRIX_H
#define DLA_INTERFACE_LOCAL_MATRIX_H

#include <memory>
#include <utility>
#include "memory_allocator.h"
#include "types.h"
#include "error_message.h"
#include "util_math.h"
#include "util_memory.h"

namespace dla_interface {

  template <class ElType>
  class LocalMatrix {
    public:
    using ElementType = ElType;

    // Creates a (0 x 0) matrix with leading dimension 1.
    LocalMatrix();

    // Creates a new (m x n) local matrix with leading dimension >= m, if m > 0 and n > 0.
    // If m = 0 or n = 0 a (0 x 0) matrix with leading dimension 1 is created.
    // Throws std::invalid_argument if m < 0 or n < 0.
    LocalMatrix(SizeType m, SizeType n);

    // Creates a new (m x n) local matrix with leading dimension ld, if m > 0 and n > 0.
    // If m == 0 or n == 0 a (0 x 0) matrix with leading dimension ld is created.
    // Throws std::invalid_argument if m < 0 or n < 0 or ld < 1.
    // Throws std::invalid_argument if ld < m and n > 0.
    LocalMatrix(SizeType m, SizeType n, SizeType ld);

    // If m > 0 and n > 0 creates a local matrix (m x n) matrix using the memory provided using ptr.
    // See note (1).
    // If m == 0 or n == 0 a (0 x 0) matrix with leading dimension ld is created.
    // Throws std::invalid_argument if m < 0 or n < 0 or ld < 1.
    // Throws std::invalid_argument if ld < m and n > 0.
    // Throws std::invalid_argument if ptr == nullptr and (m == 0 or n == 0).
    LocalMatrix(SizeType m, SizeType n, ElementType* ptr, SizeType ld);

    // Creates a copy of rhs, i.e. a matrix with the same size,
    // the same leading dimension and a copy of the elements.
    LocalMatrix(const LocalMatrix& rhs);

    // Creates a matrix of the same size and leading dimension moving the elements of rhs.
    // Postcondition: rhs is in an unspecified state.
    LocalMatrix(LocalMatrix&& rhs);

    // Copy assignement operator:
    // Changes the size and leading dimension of *this to match rhs values
    // and copies the elements of rhs.
    // Note: - new memory is allocated,
    //       - old memory is deallocated if:
    //             - the memory was allocated via the Local/DistributedMatrixClass constructors and
    //             - the memory is not used by any other Local/DistributedMatrixClass objects.
    // All the pointers and references to elements of *this are invalidated.
    LocalMatrix& operator=(const LocalMatrix& rhs);

    // Move-assignement operator:
    // Changes the size and leading dimension of *this to match rhs values
    // and moves the elements of rhs.
    // All the pointers and references to elements of *this and rhs are invalidated.
    // Postcondition: rhs is in a unspecified state.
    LocalMatrix& operator=(LocalMatrix&& rhs);

    // Copies the value of the elements of rhs.
    // Throws std::invalid_argument if *this and rhs do not have the same sizes.
    LocalMatrix& copy(const LocalMatrix& rhs);

    // Creates a LocalMatrix which represent the (m, n) submatrix starting at position (i ,j),
    // and returns a shared_ptr to it.
    // (If m == 0 or n == 0 a new (0 x 0) matrix is returned.)
    // Throws std::invalid_argument if m, n, i, j < 0, i + m > size().first, j + n > size().second
    // See note (2).
    std::shared_ptr<LocalMatrix> subMatrix(SizeType m, SizeType n, IndexType i, IndexType j) {
      return std::shared_ptr<LocalMatrix<ElementType>>(subMatrixInternal(__func__, m, n, i, j));
    }
    std::shared_ptr<const LocalMatrix> subMatrix(SizeType m, SizeType n, IndexType i,
                                                 IndexType j) const {
      return std::shared_ptr<const LocalMatrix<ElementType>>(subMatrixInternal(__func__, m, n, i, j));
    }

    // Returns a reference to the element at position (i, j).
    // Throws std::invalid_argument if i < 0, i > size().first, j < 0 or j > size().second.
    const ElementType& operator()(IndexType i, IndexType j) const {
      return *ptr(i, j);
    }
    ElementType& operator()(IndexType i, IndexType j) {
      return *ptr(i, j);
    }

    // Returns nullptr if size().first == 0 or size().second == 0.
    // Returns ptr(0, 0) otherwise.
    const ElementType* ptr() const {
      return ptr_->ptr(offset_);
    }
    ElementType* ptr() {
      return ptr_->ptr(offset_);
    }

    // Returns the pointer to the element at position (i, j).
    // Throws std::invalid_argument if i < 0, i > size().first, j < 0 or j > size().second.
    const ElementType* ptr(IndexType i, IndexType j) const {
      checkIndex(__func__, i, j);
      return ptr_->ptr(offset_ + i + util::multSize(ld_, j));
    }
    ElementType* ptr(IndexType i, IndexType j) {
      checkIndex(__func__, i, j);
      return ptr_->ptr(offset_ + i + util::multSize(ld_, j));
    }

    // Returns a pair containing the size of the matrix.
    std::pair<SizeType, SizeType> size() const {
      return size_;
    }

    // Returns the leading dimension.
    SizeType leadingDimension() const {
      return ld_;
    }

    // (1) The returned LocalMatrix will use the memory provided by the user through the pointer.
    //     Requirement: The memory region [ptr, ptr + m + ld * (n-1)) has to be valid for the
    //     lifetime of
    //     the LocalMatrix object.
    // (2) The returned LocalMatrix will use part of the memory allocated for the original matrix.
    //     If the original matrix is destroyed, the memory of the original matrix is not deallocated
    //     until all submatrices are destroyed.
    //     If elements are changed in the original matrix, the submatrix will change as well and
    //     vice-versa.
    //     If the original matrix is re-assigned using the assignement or move-assignement operator
    //     the submatrix object is not modified.
    //     Exception: if the memory has been provided by the user (see note (1)), the user has to
    //                guarantee that the memory is valid for the lifetime of the submatrix as well.

    private:
    LocalMatrix(const char* func, SizeType m, SizeType n, SizeType ld,
                std::shared_ptr<memory::MemoryAllocator<ElementType>> ptr, size_t offset);

    void checkSizeAndLd(const char* func);

    void checkIndex(const char* func, int i, int j) const;

    LocalMatrix* subMatrixInternal(const char* func, SizeType m, SizeType n, IndexType i,
                                   IndexType j) const;

    static int getLd(int m, int n);

    std::shared_ptr<memory::MemoryAllocator<ElementType>> ptr_;
    std::size_t offset_;
    std::pair<SizeType, SizeType> size_;
    SizeType ld_;
  };

  // Implementation of public member functions
  template <class ElType>
  LocalMatrix<ElType>::LocalMatrix() : LocalMatrix(0, 0, 1) {}

  template <class ElType>
  LocalMatrix<ElType>::LocalMatrix(SizeType m, SizeType n) : LocalMatrix(m, n, getLd(m, n)){};

  template <class ElType>
  LocalMatrix<ElType>::LocalMatrix(SizeType m, SizeType n, SizeType ld)
   : offset_(0), size_(m, n), ld_(ld) {
    checkSizeAndLd(__func__);
    ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(util::multSize(ld_, size_.second));
  }

  template <class ElType>
  LocalMatrix<ElType>::LocalMatrix(SizeType m, SizeType n, ElementType* ptr, SizeType ld)
   : offset_(0), size_(m, n), ld_(ld) {
    checkSizeAndLd(__func__);
    std::size_t len = util::multSize(ld_, size_.second);
    if (len == 0) {
      ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>();
    }
    else {
      ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(ptr, len);
    }
  }

  template <class ElType>
  LocalMatrix<ElType>::LocalMatrix(const LocalMatrix& rhs)
   : offset_(0), size_(rhs.size_), ld_(rhs.ld_) {
    ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(util::multSize(ld_, size_.second));
    util::memory::copy2D(size_, rhs.ptr(), rhs.ld_, ptr(), ld_);
  }

  template <class ElType>
  LocalMatrix<ElType>::LocalMatrix(LocalMatrix&& rhs)
   : ptr_(rhs.ptr_), offset_(rhs.offset_), size_(rhs.size_), ld_(rhs.ld_) {
    rhs.ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(0);
    rhs.offset_ = 0;
    rhs.size_ = std::make_pair(0, 0);
    rhs.ld_ = 0;
  }

  template <class ElType>
  LocalMatrix<ElType>& LocalMatrix<ElType>::operator=(const LocalMatrix& rhs) {
    if (this != &rhs) {
      offset_ = 0;
      size_ = rhs.size_;
      ld_ = rhs.ld_;
      ptr_ =
          std::make_shared<memory::MemoryAllocator<ElementType>>(util::multSize(ld_, size_.second));
      util::memory::copy2D(size_, rhs.ptr(), rhs.ld_, ptr(), ld_);
    }

    return *this;
  }

  template <class ElType>
  LocalMatrix<ElType>& LocalMatrix<ElType>::operator=(LocalMatrix&& rhs) {
    if (this != &rhs) {
      offset_ = rhs.offset_;
      size_ = rhs.size_;
      ld_ = rhs.ld_;
      ptr_ = rhs.ptr_;
      rhs.ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(0);
      rhs.offset_ = 0;
      rhs.size_ = std::make_pair(0, 0);
      rhs.ld_ = 0;
    }

    return *this;
  }

  template <class ElType>
  LocalMatrix<ElType>& LocalMatrix<ElType>::copy(const LocalMatrix& rhs) {
    if (size_ != rhs.size_)
      throw std::invalid_argument(errorMessage("Sizes do not match: ", size_, " != ", rhs.size_));
    util::memory::copy2D(size_, rhs.ptr(), rhs.ld_, ptr(), ld_);
    return *this;
  }

  // Implementation of private member functions
  template <class ElType>
  LocalMatrix<ElType>::LocalMatrix(const char* func, SizeType m, SizeType n, SizeType ld,
                                   std::shared_ptr<memory::MemoryAllocator<ElementType>> ptr,
                                   size_t offset)
   : ptr_(ptr), offset_(offset), size_(m, n), ld_(ld) {
    checkSizeAndLd(func);
    if (offset + util::multSize(ld_, size_.second) > ptr_->size())
      throw std::invalid_argument(errorMessageFunc(
          func, "Allocated memory (", ptr_->size(), ") is smaller than needed memory (",
          offset + util::multSize(ld_, size_.second), ")."));
  }

  template <class ElType>
  void LocalMatrix<ElType>::checkSizeAndLd(const char* func) {
    if (size_.first < 0 || size_.second < 0)
      throw std::invalid_argument(errorMessageFunc(func, "Invalid matrix size ", size_));
    int ld_min = size_.first;
    if (size_.first == 0 || size_.second == 0) {
      size_ = std::make_pair(0, 0);
      ld_min = 1;
      offset_ = 0;
    }
    if (ld_ < ld_min)
      throw std::invalid_argument(
          errorMessageFunc(func, "ld (", ld_, " < ", ld_min, ") is too small."));
  }

  template <class ElType>
  void LocalMatrix<ElType>::checkIndex(const char* func, int i, int j) const {
    if (i < 0 || i >= size_.first)
      throw std::invalid_argument(
          errorMessageFunc(func, "Index i ", i, " out of range [0, ", size_.first, ")"));
    if (j < 0 || j >= size_.second)
      throw std::invalid_argument(
          errorMessageFunc(func, "Index j ", j, " out of range [0, ", size_.second, ")"));
  }

  template <class ElType>
  LocalMatrix<ElType>* LocalMatrix<ElType>::subMatrixInternal(const char* func, SizeType m,
                                                              SizeType n, IndexType i,
                                                              IndexType j) const {
    checkIndex(func, i, j);
    if (m < 0 || n < 0)
      throw std::invalid_argument(
          errorMessageFunc(func, "Invalid submatrix size (", m, ", ", n, ")"));
    if (i + m > size_.first || j + n > size_.second)
      throw std::invalid_argument(errorMessageFunc(func, "Submatrix is exceeding matrix borders (",
                                                   i + m, ", ", j + n, ") larger than ", size_));
    if (m == 0 || n == 0)
      return new LocalMatrix<ElementType>(0, 0, ld_);

    return new LocalMatrix<ElementType>(func, m, n, ld_, ptr_, i + util::multSize(ld_, j));
  }

  template <class ElType>
  int LocalMatrix<ElType>::getLd(int m, int n) {
    if (m == 0 || n == 0)
      return 1;
    constexpr int chunk = 64;
    return util::ceilDiv(m, chunk) * chunk;
  }
}

#endif  // DLA_INTERFACE_LOCAL_MATRIX_H
