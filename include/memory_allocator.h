//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2021, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_MEMORY_ALLOCATOR_H
#define DLA_INTERFACE_MEMORY_ALLOCATOR_H

#include <cassert>
#include <cstdlib>
#include "error_message.h"

namespace dla_interface {
  namespace memory {

    template <class ElType>
    class MemoryAllocator {
      public:
      using ElementType = ElType;

      /// Creates nullptr object.
      MemoryAllocator() : ptr_(nullptr), size_(0), allocated_(false) {}
      // Allocates size ElementType and cretes an object.
      MemoryAllocator(size_t size) : ptr_(nullptr), size_(size), allocated_(false) {
        allocate();
      }

      MemoryAllocator(const MemoryAllocator&) = delete;
      MemoryAllocator(MemoryAllocator&& rhs)
       : ptr_(rhs.ptr_), size_(rhs.size_), allocated_(rhs.allocated_) {
        rhs.ptr_ = nullptr;
        rhs.size_ = 0;
        rhs.allocated_ = false;
      }

      MemoryAllocator& operator=(const MemoryAllocator&) = delete;
      MemoryAllocator& operator=(MemoryAllocator&& rhs) {
        deallocate();
        ptr_ = rhs.ptr_;
        size_ = rhs.size_;
        allocated_ = rhs.allocated_;
        rhs.ptr_ = nullptr;
        rhs.size_ = 0;
        rhs.allocated_ = false;
        return *this;
      }

      /// Use the memory pointed by ptr and creates an object.
      ///
      /// <b>Precondition:</b> ptr[i] is valid for i < size.
      MemoryAllocator(ElementType* ptr, size_t size) : ptr_(ptr), size_(size), allocated_(false) {
        if (size > 0 && ptr == nullptr)
          throw std::invalid_argument(errorMessage("ptr == nullptr for size (", size, ") > 0."));
      }

      /// Deallocates the elements if they were allocated by this class.
      ~MemoryAllocator() {
        deallocate();
      }

      /// Returns the i-th element.
      ///
      /// <b>Precondition:</b> i == 0 or i < size_;
      const ElementType* ptr(std::size_t i = 0) const {
        assert(i < size_ || i == 0);
        return ptr_ + i;
      }
      ElementType* ptr(std::size_t i = 0) {
        assert(i < size_ || i == 0);
        return ptr_ + i;
      }

      std::size_t size() const {
        return size_;
      }

      private:
      void allocate() {
        if (size_ > 0) {
          ptr_ = reinterpret_cast<ElementType*>(std::malloc(sizeof(ElementType) * size_));
          allocated_ = true;
        }
      }

      void deallocate() {
        if (allocated_) {
          std::free(ptr_);
          allocated_ = false;
        }
        ptr_ = nullptr;
      }

      ElementType* ptr_;
      std::size_t size_;
      bool allocated_;
    };
  }
}

#endif  // DLA_INTERFACE_MEMORY_ALLOCATOR_H
