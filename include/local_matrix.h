#ifndef DLA_INTERFACE_LOCAL_MATRIX_H
#define DLA_INTERFACE_LOCAL_MATRIX_H

#include <memory>
#include <utility>
#include "types.h"

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
    std::shared_ptr<LocalMatrix> subMatrix(SizeType m, SizeType n, IndexType i, IndexType j);
    std::shared_ptr<const LocalMatrix> subMatrix(SizeType m, SizeType n, IndexType i,
                                                 IndexType j) const;

    // Returns a reference to the element at position (i, j).
    // Throws std::invalid_argument if i < 0, i > size().first, j < 0 or j > size().second.
    const ElementType& operator()(IndexType i, IndexType j) const;
    ElementType& operator()(IndexType i, IndexType j);

    // Returns nullptr if size().first == 0 or size().second == 0.
    // Returns ptr(0, 0) otherwise.
    const ElementType* ptr() const;
    ElementType* ptr();

    // Returns the pointer to the element at position (i, j).
    // Throws std::invalid_argument if i < 0, i > size().first, j < 0 or j > size().second.
    const ElementType* ptr(IndexType i, IndexType j) const;
    ElementType* ptr(IndexType i, IndexType j);

    // Returns a pair containing the size of the matrix.
    std::pair<SizeType, SizeType> size() const;

    // Returns the leading dimension.
    SizeType leadingDimension() const;

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
  };
}

#endif  // DLA_INTERFACE_LOCAL_MATRIX_H
