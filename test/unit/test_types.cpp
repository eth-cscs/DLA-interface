#include "types.h"

#include "gtest/gtest.h"

using namespace dla_interface;
using namespace testing;

TEST(TypesTest, BaseComplexTypes) {
  // ElementType check
  EXPECT_TRUE((std::is_same<float, TypeInfo<float>::ElementType>::value));
  EXPECT_TRUE((std::is_same<double, TypeInfo<double>::ElementType>::value));
  EXPECT_TRUE((std::is_same<std::complex<float>, TypeInfo<std::complex<float>>::ElementType>::value));
  EXPECT_TRUE(
      (std::is_same<std::complex<double>, TypeInfo<std::complex<double>>::ElementType>::value));

  // BaseType check
  EXPECT_TRUE((std::is_same<float, TypeInfo<float>::BaseType>::value));
  EXPECT_TRUE((std::is_same<double, TypeInfo<double>::BaseType>::value));
  EXPECT_TRUE((std::is_same<float, TypeInfo<std::complex<float>>::BaseType>::value));
  EXPECT_TRUE((std::is_same<double, TypeInfo<std::complex<double>>::BaseType>::value));
  EXPECT_TRUE((std::is_same<float, BaseType<float>>::value));
  EXPECT_TRUE((std::is_same<double, BaseType<double>>::value));
  EXPECT_TRUE((std::is_same<float, BaseType<std::complex<float>>>::value));
  EXPECT_TRUE((std::is_same<double, BaseType<std::complex<double>>>::value));

  // ComplexType
  EXPECT_TRUE((std::is_same<std::complex<float>, TypeInfo<float>::ComplexType>::value));
  EXPECT_TRUE((std::is_same<std::complex<float>, TypeInfo<std::complex<float>>::ComplexType>::value));
  EXPECT_TRUE((std::is_same<std::complex<double>, TypeInfo<double>::ComplexType>::value));
  EXPECT_TRUE(
      (std::is_same<std::complex<double>, TypeInfo<std::complex<double>>::ComplexType>::value));
  EXPECT_TRUE((std::is_same<std::complex<float>, ComplexType<float>>::value));
  EXPECT_TRUE((std::is_same<std::complex<float>, ComplexType<std::complex<float>>>::value));
  EXPECT_TRUE((std::is_same<std::complex<double>, ComplexType<double>>::value));
  EXPECT_TRUE((std::is_same<std::complex<double>, ComplexType<std::complex<double>>>::value));
}

TEST(TypesTest, Zero) {
  EXPECT_EQ(0, I_ZERO);
  EXPECT_EQ(0.f, S_ZERO);
  EXPECT_EQ(0., D_ZERO);
  EXPECT_EQ(std::complex<float>(0.f), C_ZERO);
  EXPECT_EQ(std::complex<double>(0.f), Z_ZERO);

  EXPECT_EQ(0, *DLA_I_ZERO);
  EXPECT_EQ(0.f, *DLA_S_ZERO);
  EXPECT_EQ(0., *DLA_D_ZERO);
  EXPECT_EQ(std::complex<float>(0.f), *DLA_C_ZERO);
  EXPECT_EQ(std::complex<double>(0.f), *DLA_Z_ZERO);
}

TEST(TypesTest, One) {
  EXPECT_EQ(1, I_ONE);
  EXPECT_EQ(1.f, S_ONE);
  EXPECT_EQ(1., D_ONE);
  EXPECT_EQ(std::complex<float>(1.f), C_ONE);
  EXPECT_EQ(std::complex<double>(1.f), Z_ONE);

  EXPECT_EQ(1, *DLA_I_ONE);
  EXPECT_EQ(1.f, *DLA_S_ONE);
  EXPECT_EQ(1., *DLA_D_ONE);
  EXPECT_EQ(std::complex<float>(1.f), *DLA_C_ONE);
  EXPECT_EQ(std::complex<double>(1.f), *DLA_Z_ONE);
}
