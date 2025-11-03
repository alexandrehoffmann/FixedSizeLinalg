#include <gtest/gtest.h>

#include <FSLinalg/Vector.hpp>
#include <FSLinalg/Matrix.hpp>

TEST(lazy, matrixVectorProduct)
{
	FSLinalg::UnitVector<4> e0(0);
	FSLinalg::RealVector<3> a({2, 3, 0});
	
	const auto expr1 = 0.5*FSLinalg::outer(e0, a)*a;
	const FSLinalg::RealVector<4> expected1({6.5, 0, 0, 0});
	
	EXPECT_EQ(FSLinalg::RealVector<4>(expr1), expected1);
	EXPECT_TRUE  (expr1.createTemporaryMatrix);
	EXPECT_FALSE (expr1.createTemporaryVector);
	
	const auto expr2 = 0.5*FSLinalg::transpose(FSLinalg::outer(e0, a))*e0;
	const FSLinalg::RealVector<3> expected2({1, 1.5, 0});
	
	EXPECT_EQ(FSLinalg::RealVector<3>(expr2), expected2);
	EXPECT_TRUE  (expr2.createTemporaryMatrix);
	EXPECT_FALSE (expr2.createTemporaryVector);
	
	FSLinalg::RealMatrix<4, 3> A(FSLinalg::outer(e0, a));
	
	const auto expr3 = 0.5*FSLinalg::transpose(A)*e0;
	
	EXPECT_EQ(FSLinalg::RealVector<3>(expr3), expected2);
	EXPECT_FALSE(expr3.createTemporaryMatrix);
	EXPECT_FALSE(expr3.createTemporaryVector);
}
