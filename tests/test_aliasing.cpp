#include <gtest/gtest.h>

#include <FSLinalg/Matrix.hpp>

TEST(aliasing, cross)
{
	FSLinalg::RealRowVector<3> a({2, 3, 1});
	FSLinalg::RealRowVector<3> b({4, 6, 5});
	
	const auto expr2 = -2.*FSLinalg::cross(a, b);
	
	const FSLinalg::RealRowVector<3> c(expr2);
	const FSLinalg::RealRowVector<3> expected({-18, 12, 0});
	
	EXPECT_EQ(c, expected);
	
	a = expr2;
	
	EXPECT_EQ(a, expected);
}

TEST(aliasing, outer)
{
	FSLinalg::UnitRowVector<3> e0(0);
	FSLinalg::RealRowVector<3> a({2, 3, 0});
	
	const auto expr = 0.5*FSLinalg::outer(e0, a)*a;
	const FSLinalg::RealRowVector<3> expected({6.5, 0, 0});
	
	EXPECT_EQ(FSLinalg::RealRowVector<3>(expr), expected);
	
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

TEST(aliasing, product)
{
	FSLinalg::RealMatrix<6,6> A({
		{2, 1, 3, 4, 6, 5},
		{7, 8, 4, 8, 5, 9},
		{2, 4, 2, 2, 0, 5},
		{1, 8, 4, 5, 2, 0},
		{1, 9, 5, 4, 6, 9},
		{3, 5, 2, 6, 8, 9}});
		
	const FSLinalg::RealMatrix<6,6> B({
		{6, 2, 3, 4, 7, 1},
		{7, 6, 7, 9, 6, 3},
		{7, 2, 6, 1, 4, 5},
		{2, 4, 0, 7, 1, 5},
		{3, 2, 5, 8, 2, 0},
		{7, 0, 9, 3, 3, 7}});
		
	const auto expr = A*B;
	
	const FSLinalg::RealMatrix<6,6> expected({
		{101, 44, 106, 111, 63, 75},
		{220, 112, 207, 227, 158, 154},
		{93, 40, 91, 75, 63, 69},
		{106, 82, 93, 131, 80, 70},
		{193, 94, 207, 193, 124, 136},
		{166, 80, 177, 192, 108, 121}});
	
	A = expr;
	
	EXPECT_EQ(A, expected);
}
