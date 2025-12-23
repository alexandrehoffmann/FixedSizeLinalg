#include <gtest/gtest.h>

#include <FSLinalg/Tensor.hpp>

template<size_t rank> using Shape = std::array<unsigned int, rank>;

TEST(tensor, shape_to_strides)
{	
	{
		Shape<1> shape{12};
		Shape<1> expected{1};
		
		EXPECT_EQ(FSLinalg::TensorUtils::getStrides(shape), expected);
	}
	{
		Shape<2> shape{3, 4};
		Shape<2> expected{4, 1};
		
		EXPECT_EQ(FSLinalg::TensorUtils::getStrides(shape), expected);
	}
	{
		Shape<3> shape{3, 2, 4};
		Shape<3> expected{8, 4, 1};
		
		EXPECT_EQ(FSLinalg::TensorUtils::getStrides(shape), expected);
	}
	{
		Shape<4> shape{3, 2, 4, 8};
		Shape<4> expected{64, 32, 8, 1};
		
		EXPECT_EQ(FSLinalg::TensorUtils::getStrides(shape), expected);
	}
}

TEST(tensor, binary_ops)
{
	FSLinalg::RealTensor<3> a({2, 4, 1});
	FSLinalg::RealTensor<3> b({1, 5, -1});
	
	const auto expr1 = a + b;
	const auto expr2 = a - b;
	
	FSLinalg::RealTensor<3> expected1({3, 9, 0});
	FSLinalg::RealTensor<3> expected2({1, -1, 2});
	
	EXPECT_EQ(expr1, expected1);
	EXPECT_EQ(expr2, expected2);
	
	const auto expr3 = expr1 > expr2;
	
	FSLinalg::BoolTensor<3> expected3({true, true, false});
	
	EXPECT_EQ(expr3, expected3);
}
