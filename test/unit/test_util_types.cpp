#include "util_types.h"

#include <exception>
#include <string>
#include "gtest/gtest.h"
#include "types.h"

using namespace dla_interface;
using namespace testing;

TEST(UtilTypesTest, NrOps) {
  EXPECT_DOUBLE_EQ(3., util::nrOps<float>(3, 0));
  EXPECT_DOUBLE_EQ(4., util::nrOps<double>(4, 0));
  EXPECT_DOUBLE_EQ(6., util::nrOps<std::complex<float>>(3, 0));
  EXPECT_DOUBLE_EQ(8., util::nrOps<std::complex<double>>(4, 0));

  EXPECT_DOUBLE_EQ(6., util::nrOps<float>(0, 6));
  EXPECT_DOUBLE_EQ(5., util::nrOps<double>(0, 5));
  EXPECT_DOUBLE_EQ(36., util::nrOps<std::complex<float>>(0, 6));
  EXPECT_DOUBLE_EQ(30., util::nrOps<std::complex<double>>(0, 5));

  EXPECT_DOUBLE_EQ(13., util::nrOps<float>(6, 7));
  EXPECT_DOUBLE_EQ(12., util::nrOps<double>(9, 3));
  EXPECT_DOUBLE_EQ(54., util::nrOps<std::complex<float>>(6, 7));
  EXPECT_DOUBLE_EQ(36., util::nrOps<std::complex<double>>(9, 3));
}

TEST(UtilTypesTest, getOpTrans) {
  EXPECT_EQ(NoTrans, util::getOpTrans('N'));
  EXPECT_EQ(NoTrans, util::getOpTrans('n'));
  EXPECT_EQ(Trans, util::getOpTrans('T'));
  EXPECT_EQ(Trans, util::getOpTrans('t'));
  EXPECT_EQ(ConjTrans, util::getOpTrans('C'));
  EXPECT_EQ(ConjTrans, util::getOpTrans('c'));

  EXPECT_THROW(util::getOpTrans('U'), std::invalid_argument);
  EXPECT_THROW(util::getOpTrans('u'), std::invalid_argument);
  EXPECT_THROW(util::getOpTrans('L'), std::invalid_argument);
  EXPECT_THROW(util::getOpTrans('l'), std::invalid_argument);
  EXPECT_THROW(util::getOpTrans('R'), std::invalid_argument);
  EXPECT_THROW(util::getOpTrans('r'), std::invalid_argument);

  EXPECT_THROW(util::getOpTrans('3'), std::invalid_argument);
  EXPECT_THROW(util::getOpTrans('Q'), std::invalid_argument);
  EXPECT_THROW(util::getOpTrans('{'), std::invalid_argument);
}

TEST(UtilTypesTest, getUpLo) {
  EXPECT_EQ(Upper, util::getUpLo('U'));
  EXPECT_EQ(Upper, util::getUpLo('u'));
  EXPECT_EQ(Lower, util::getUpLo('L'));
  EXPECT_EQ(Lower, util::getUpLo('l'));

  EXPECT_THROW(util::getUpLo('N'), std::invalid_argument);
  EXPECT_THROW(util::getUpLo('n'), std::invalid_argument);
  EXPECT_THROW(util::getUpLo('T'), std::invalid_argument);
  EXPECT_THROW(util::getUpLo('t'), std::invalid_argument);
  EXPECT_THROW(util::getUpLo('C'), std::invalid_argument);
  EXPECT_THROW(util::getUpLo('c'), std::invalid_argument);
  EXPECT_THROW(util::getUpLo('R'), std::invalid_argument);
  EXPECT_THROW(util::getUpLo('r'), std::invalid_argument);

  EXPECT_THROW(util::getUpLo('3'), std::invalid_argument);
  EXPECT_THROW(util::getUpLo('Q'), std::invalid_argument);
  EXPECT_THROW(util::getUpLo('{'), std::invalid_argument);
}

TEST(UtilTypesTest, getOrdering) {
  EXPECT_EQ(RowMajor, util::getOrdering('R'));
  EXPECT_EQ(RowMajor, util::getOrdering('r'));
  EXPECT_EQ(ColMajor, util::getOrdering('C'));
  EXPECT_EQ(ColMajor, util::getOrdering('c'));

  EXPECT_THROW(util::getOrdering('N'), std::invalid_argument);
  EXPECT_THROW(util::getOrdering('n'), std::invalid_argument);
  EXPECT_THROW(util::getOrdering('T'), std::invalid_argument);
  EXPECT_THROW(util::getOrdering('t'), std::invalid_argument);
  EXPECT_THROW(util::getOrdering('U'), std::invalid_argument);
  EXPECT_THROW(util::getOrdering('u'), std::invalid_argument);
  EXPECT_THROW(util::getOrdering('L'), std::invalid_argument);
  EXPECT_THROW(util::getOrdering('l'), std::invalid_argument);

  EXPECT_THROW(util::getOrdering('3'), std::invalid_argument);
  EXPECT_THROW(util::getOrdering('Q'), std::invalid_argument);
  EXPECT_THROW(util::getOrdering('{'), std::invalid_argument);
}

TEST(UtilTypesTest, SolverString) {
  EXPECT_EQ(std::string("ScaLAPACK"), util::getSolverString(ScaLAPACK));
  EXPECT_EQ(std::string("ELPA"), util::getSolverString(ELPA));
  EXPECT_EQ(std::string("DPLASMA"), util::getSolverString(DPlasma));
  EXPECT_EQ(std::string("Chameleon"), util::getSolverString(Chameleon));

  EXPECT_EQ(ScaLAPACK, util::getSolverType("ScaLAPACK"));
  EXPECT_EQ(ELPA, util::getSolverType("ELPA"));
  EXPECT_EQ(DPlasma, util::getSolverType("DPLASMA"));
  EXPECT_EQ(Chameleon, util::getSolverType("Chameleon"));
}
