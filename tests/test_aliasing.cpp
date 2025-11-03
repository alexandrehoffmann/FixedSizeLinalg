#include <gtest/gtest.h>

#include <FSLinalg/Vector.hpp>
#include <FSLinalg/Matrix.hpp>

TEST(aliasing, cross)
{
	FSLinalg::RealVector<3> a({2, 3, 1});
	FSLinalg::RealVector<3> b({4, 6, 5});
	
	const auto expr2 = -2.*FSLinalg::cross(a, b);
	
	const FSLinalg::RealVector<3> c(expr2);
	const FSLinalg::RealVector<3> expected({-18, 12, 0});
	
	EXPECT_EQ(c, expected);
	
	a = expr2;
	
	EXPECT_EQ(a, expected);
}

TEST(aliasing, outer)
{
	FSLinalg::UnitVector<3> e0(0);
	FSLinalg::RealVector<3> a({2, 3, 0});
	
	const auto expr = 0.5*FSLinalg::outer(e0, a)*a;
	const FSLinalg::RealVector<3> expected({6.5, 0, 0});
	
	EXPECT_EQ(FSLinalg::RealVector<3>(expr), expected);
	
	a = expr;
	
	EXPECT_EQ(a, expected);
}

TEST(aliasing, transposed)
{
	FSLinalg::RealMatrix<6,6> A({
		{2, 1, 3, 4, 6, 5},
		{7, 8, 4, 8, 5, 9},
		{2, 4, 2, 2, 0, 5},
		{1, 8, 4, 5, 2, 0},
		{1, 9, 5, 4, 6, 9},
		{3, 5, 2, 6, 8, 9}});
		
	const auto expr = FSLinalg::transpose(A);
	
	const FSLinalg::RealMatrix<6,6> expected({
		{2, 7, 2, 1, 1, 3},
		{1, 8, 4, 8, 9, 5},
		{3, 4, 2, 4, 5, 2},
		{4, 8, 2, 5, 4, 6},
		{6, 5, 0, 2, 6, 8},
		{5, 9, 5, 0, 9, 9}});
	
	EXPECT_EQ(expr, expected);
	
	A = expr;
	
	EXPECT_EQ(A, expected);
}
