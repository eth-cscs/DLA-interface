//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2019, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

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
  /// Forward declaration for friend.
  template <class ElType>
  class DistributedMatrix;

  template <class ElType>
  class LocalMatrix {

	///<br><br>
	/// \note
	/// <b>(1)</b> \anchor note_01<br>
	/// The returned LocalMatrix will use the memory provided by the user through the pointer.
	/// <b>Requirement:</b> The memory region [ptr, ptr + m + ld * (n-1)) has to be valid for the
	/// lifetime of the LocalMatrix object.
	///
	///<br>
	/// <b>(2)</b> \anchor note_02<br>
	/// The returned LocalMatrix will use part of the memory allocated for the original matrix.
	/// If the original matrix is destroyed, the memory of the original matrix is not deallocated
	/// until all submatrices are destroyed.
	/// If elements are changed in the original matrix, the submatrix will change as well and
	/// vice-versa.
	/// If the original matrix is re-assigned using the assignement or move-assignement operator
	/// the submatrix object is not modified.
	/// <b>Exception:</b> if the memory has been provided by the user (see \ref note_01 "Note (1)"), the user has to
	/// guarantee that the memory is valid for the lifetime of the submatrix as well.

    public:
    using ElementType = ElType;

    /// Creates a (0 x 0) matrix with leading dimension 1.
    LocalMatrix();

    /// Creates a new (m x n) local matrix with leading dimension >= m, if m > 0 and n > 0.
    /// If m = 0 or n = 0 a (0 x 0) matrix with leading dimension 1 is created.
    ///
    /// @throws std::invalid_argument if m < 0 or n < 0.
    LocalMatrix(SizeType m, SizeType n);

    /// Creates a new (m x n) local matrix with leading dimension ld, if m > 0 and n > 0.
    /// If m == 0 or n == 0 a (0 x 0) matrix with leading dimension ld is created.
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if m < 0 or n < 0 or ld < 1,</li>
    /// 	<li>if ld < m and n > 0.</li>
    /// <ul>
    LocalMatrix(SizeType m, SizeType n, SizeType ld);

    /// If m > 0 and n > 0 creates a local matrix (m x n) matrix using the memory provided using ptr.
    /// (see \ref note_01 "Note (1)").
    /// If m == 0 or n == 0 a (0 x 0) matrix with leading dimension ld is created.
    ///
    /// @throws std::invalid_argument
    /// <ul>
    /// 	<li>if m < 0 or n < 0 or ld < 1,</li>
    /// 	<li>if ld < m and n > 0,</li>
    ///     <li>if ptr == nullptr and (m == 0 or n == 0).</li>
    /// </ul>
    LocalMatrix(SizeType m, SizeType n, ElementType* ptr, SizeType ld);

    /// Creates a copy of rhs, i.e. a matrix with the same size,
    /// the same leading dimension and a copy of the elements.
    LocalMatrix(const LocalMatrix& rhs);

    /// Creates a matrix of the same size and leading dimension moving the elements of rhs.
    /// <b>Postcondition:</b> rhs is in an unspecified state.
    LocalMatrix(LocalMatrix&& rhs);

    /// Copy assignement operator:
    /// Changes the size and leading dimension of *this to match rhs values
    /// and copies the elements of rhs.
    ///
    /// @note
    /// <ul>
    /// 	<li>new memory is allocated,</li>
    /// 	<li>old memory is deallocated if:</li>
    /// 	<ul>
    /// 		<li>the memory was allocated via the Local/DistributedMatrixClass constructors and</li>
    /// 		<li>the memory is not used by any other Local/DistributedMatrixClass objects.</li>
    /// 	</ul>
    /// <ul>
    ///
    /// All the pointers and references to elements of *this are invalidated.
    LocalMatrix& operator=(const LocalMatrix& rhs);

    /// Move-assignement operator:<br>
    /// Changes the size and leading dimension of *this to match rhs values
    /// and moves the elements of rhs.
    /// All the pointers and references to elements of *this and rhs are invalidated.
    ///
    /// <b>Postcondition:</b> rhs is in a unspecified state.
    LocalMatrix& operator=(LocalMatrix&& rhs);

    /// Copies the value of the elements of rhs.
    ///
    /// @throws std::invalid_argument if *this and rhs do not have the same sizes.
    LocalMatrix& copy(const LocalMatrix& rhs);

    /// Creates a LocalMatrix which represent the (m, n) submatrix starting at position (i ,j),
    /// and returns a shared_ptr to it.
    /// (If m == 0 or n == 0 a new (0 x 0) matrix is returned.)
    ///
    /// @throws std::invalid_argument if m, n, i, j < 0, i + m > size().first, j + n > size().second
    /// (see \ref note_02 "Note (2)").
    std::shared_ptr<LocalMatrix> subMatrix(SizeType m, SizeType n, IndexType i, IndexType j) {
      return std::shared_ptr<LocalMatrix<ElementType>>(subMatrixInternal(__func__, m, n, i, j));
    }
    std::shared_ptr<const LocalMatrix> subMatrix(SizeType m, SizeType n, IndexType i,
                                                 IndexType j) const {
      return std::shared_ptr<const LocalMatrix<ElementType>>(subMatrixInternal(__func__, m, n, i, j));
    }

    /// Returns a reference to the element at position (i, j).
    ///
    /// @throws std::invalid_argument if i < 0, i > size().first, j < 0 or j > size().second.
    const ElementType& operator()(IndexType i, IndexType j) const {
      return *ptr(i, j);
    }
    ElementType& operator()(IndexType i, IndexType j) {
      return *ptr(i, j);
    }

    /// Returns nullptr if size().first == 0 or size().second == 0, or returns ptr(0, 0) otherwise.
    ///
    /// @returns
    /// <ul>
    /// 	<li>nullptr if size().first == 0 or size().second == 0,</li>
    /// 	<li>ptr(0, 0) otherwise.</li>
    /// </ul>
    const ElementType* ptr() const {
      return ptr_->ptr(offset_);
    }
    ElementType* ptr() {
      return ptr_->ptr(offset_);
    }

    /// Returns the pointer to the element at position (i, j).
    ///
    /// @throws std::invalid_argument if i < 0, i > size().first, j < 0 or j > size().second.
    const ElementType* ptr(IndexType i, IndexType j) const {
      checkIndex(__func__, i, j);
      return ptr_->ptr(offset_ + i + util::multSize(ld_, j));
    }
    ElementType* ptr(IndexType i, IndexType j) {
      checkIndex(__func__, i, j);
      return ptr_->ptr(offset_ + i + util::multSize(ld_, j));
    }

    /// Returns a pair containing the size of the matrix.
    std::pair<SizeType, SizeType> size() const {
      return size_;
    }

    /// Returns the leading dimension.
    SizeType leadingDimension() const {
      return ld_;
    }



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

    friend class DistributedMatrix<ElementType>;
  };


#include "local_matrix.ipp"
}

#endif  // DLA_INTERFACE_LOCAL_MATRIX_H
