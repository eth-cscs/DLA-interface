// Implementation of public member functions
template <class ElType>
LocalMatrix<ElType>::LocalMatrix() : LocalMatrix(0, 0, 1) {}

template <class ElType>
LocalMatrix<ElType>::LocalMatrix(SizeType m, SizeType n) : LocalMatrix(m, n, getLd(m, n)) {}

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
    ptr_ = std::make_shared<memory::MemoryAllocator<ElementType>>(util::multSize(ld_, size_.second));
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
    throw std::invalid_argument(errorMessageFunc(func, "Allocated memory (", ptr_->size(),
                                                 ") is smaller than needed memory (",
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
LocalMatrix<ElType>* LocalMatrix<ElType>::subMatrixInternal(const char* func, SizeType m, SizeType n,
                                                            IndexType i, IndexType j) const {
  checkIndex(func, i, j);
  if (m < 0 || n < 0)
    throw std::invalid_argument(errorMessageFunc(func, "Invalid submatrix size (", m, ", ", n, ")"));
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
