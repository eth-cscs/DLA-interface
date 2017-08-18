#include "memory_allocator.h"

#include "gtest/gtest.h"

using namespace dla_interface;
using namespace testing;

template <typename Type>
class MemoryAllocatorTest : public ::testing::Test {};

typedef ::testing::Types<int, float, double, unsigned long long> TestTypes;
TYPED_TEST_CASE(MemoryAllocatorTest, TestTypes);

TYPED_TEST(MemoryAllocatorTest, DefaultConstructor) {
  using Type = TypeParam;
  memory::MemoryAllocator<Type> mem;

  EXPECT_EQ(nullptr, mem.ptr());
  EXPECT_EQ(0u, mem.size());
}

TYPED_TEST(MemoryAllocatorTest, Constructor1) {
  using Type = TypeParam;
  memory::MemoryAllocator<Type> mem(7);

  EXPECT_NE(nullptr, mem.ptr());
  EXPECT_EQ(7u, mem.size());
}

TYPED_TEST(MemoryAllocatorTest, Constructor2) {
  using Type = TypeParam;
  std::vector<Type> src(5);

  memory::MemoryAllocator<Type> mem(&src[0], 5);

  EXPECT_EQ(&src[0], mem.ptr());
  EXPECT_EQ(5u, mem.size());
}

TYPED_TEST(MemoryAllocatorTest, MoveConstructor) {
  using Type = TypeParam;
  memory::MemoryAllocator<Type> mem(6);

  Type* ptr = mem.ptr();
  std::size_t size = mem.size();

  memory::MemoryAllocator<Type> mem2(std::move(mem));
  EXPECT_EQ(ptr, mem2.ptr());
  EXPECT_EQ(size, mem2.size());
}

TYPED_TEST(MemoryAllocatorTest, MoveAssignementOperator) {
  using Type = TypeParam;
  memory::MemoryAllocator<Type> mem2(10);
  {
    memory::MemoryAllocator<Type> mem1(5);
    *mem1.ptr() = static_cast<Type>(10);
    Type* ptr = mem1.ptr();
    std::size_t size = mem1.size();
    memory::MemoryAllocator<Type>* obj_ptr = &mem2;

    EXPECT_EQ(obj_ptr, &(mem2 = std::move(mem1)));
    EXPECT_EQ(ptr, mem2.ptr());
    EXPECT_EQ(size, mem2.size());
  }
  EXPECT_EQ(static_cast<Type>(10), *mem2.ptr());
}

TYPED_TEST(MemoryAllocatorTest, Ptr) {
  using Type = TypeParam;
  {
    memory::MemoryAllocator<Type> mem(5);
    Type* ptr = mem.ptr();
    std::size_t size = mem.size();
    const memory::MemoryAllocator<Type>& const_mem(mem);

    for (std::size_t i = 0; i < size; ++i) {
      EXPECT_EQ(ptr + i, mem.ptr(i));
      EXPECT_EQ(ptr + i, const_mem.ptr(i));
    }
  }

  {
    std::vector<Type> src(5);
    memory::MemoryAllocator<Type> mem(&src[0], 3);
    std::size_t size = mem.size();
    const memory::MemoryAllocator<Type>& const_mem(mem);

    for (std::size_t i = 0; i < size; ++i) {
      EXPECT_EQ(&src[i], mem.ptr(i));
      EXPECT_EQ(&src[i], const_mem.ptr(i));
    }
  }
}
