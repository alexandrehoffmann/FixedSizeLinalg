#include <gtest/gtest.h>

#include <FSLinalg/Matrix.hpp>

TEST(chain, trivial)
{
	const FSLinalg::RealMatrix<6,5> A({
		{2, 1, 3, 4, 6},
		{7, 8, 4, 8, 5},
		{2, 4, 2, 2, 0},
		{1, 8, 4, 5, 2},
		{1, 9, 5, 4, 6},
		{3, 5, 2, 6, 8}});
		
	const FSLinalg::RealRowVector<5> x({2, 6, 3, 7, 5});
	
	{
		const auto expr = A*x;
		
		using ProdAnalyzer = FSLinalg::MatrixProductAnalyzer<std::decay_t<decltype(expr)>>;
		using DimArray     = typename ProdAnalyzer::DimArray;
		
		EXPECT_EQ(ProdAnalyzer::getLength(), 2);
		
		constexpr DimArray dims = ProdAnalyzer::getDims();
		
		EXPECT_EQ(dims, DimArray({6, 5, 1}));
		
		EXPECT_EQ(ProdAnalyzer::getOptimalCost(), 30);
		EXPECT_EQ(ProdAnalyzer::getOptimalSplit(), 1);
		
		constexpr bool b1 = std::is_same<std::decay_t<decltype(expr)>, typename ProdAnalyzer::OptimalBracketing>::value;
		
		EXPECT_TRUE(b1);
	}
	{
		const FSLinalg::RealRowVector<6> y({2, 6, 3, 7, 5, 1});
		
		const auto expr = 0.5*FSLinalg::transpose(A)*y;
		
		using ProdAnalyzer = FSLinalg::MatrixProductAnalyzer<std::decay_t<decltype(expr)>>;
		using DimArray     = typename ProdAnalyzer::DimArray;
		
		EXPECT_EQ(ProdAnalyzer::getLength(), 2);
		
		constexpr DimArray dims = ProdAnalyzer::getDims();
		
		EXPECT_EQ(dims, DimArray({5, 6, 1}));
		
		EXPECT_EQ(ProdAnalyzer::getOptimalCost(), 30);
		EXPECT_EQ(ProdAnalyzer::getOptimalSplit(), 1);
		
		constexpr bool b1 = std::is_same<std::decay_t<decltype(expr)>, typename ProdAnalyzer::OptimalBracketing>::value;
		
		EXPECT_TRUE(b1);
	}
}

TEST(chain, noBracketing)
{
	FSLinalg::UnitRowVector<4> e0(0);
	FSLinalg::RealRowVector<3> a({2, 3, 0});
	
	const auto expr1 = FSLinalg::transpose(a)*FSLinalg::outer(a, e0)*FSLinalg::outer(e0, a)*a;
	
	using ProdAnalyzer = FSLinalg::MatrixProductAnalyzer<std::decay_t<decltype(expr1)>>;
	using DimArray     = typename ProdAnalyzer::DimArray;
	
	EXPECT_EQ(ProdAnalyzer::getLength(), 6);
	
	constexpr DimArray dims = ProdAnalyzer::getDims();
	
	EXPECT_EQ(dims, DimArray({1, 3, 1, 4, 1, 3, 1}));
	EXPECT_EQ(ProdAnalyzer::getOptimalCost(), 12);
	EXPECT_EQ(ProdAnalyzer::getOptimalSplit(), 2);
	
	using ExpectedExpr = decltype( (FSLinalg::transpose(a)*a)*((FSLinalg::transpose(e0)*e0)*(FSLinalg::transpose(a)*a)) );
	
	constexpr bool b1 = std::is_same<ExpectedExpr, typename ProdAnalyzer::OptimalBracketing>::value;
	
	EXPECT_TRUE(b1);
	
	using ReBracketType = typename ProdAnalyzer::ReBracketType;
	
	ReBracketType rebrackedExpr = ProdAnalyzer::reBracket(expr1);
	
	const FSLinalg::RealMatrix<1,1> result(rebrackedExpr);
	const FSLinalg::RealMatrix<1,1> expected({169.});
	
	EXPECT_EQ(result, expected);
}

TEST(chain, wrongBracketing)
{
	FSLinalg::UnitRowVector<4> e0(0);
	FSLinalg::RealRowVector<3> a({2, 3, 0});
	
	const auto expr1 = FSLinalg::transpose(a)*(FSLinalg::outer(a, e0)*FSLinalg::outer(e0, a))*a;
	
	using ProdAnalyzer = FSLinalg::MatrixProductAnalyzer<std::decay_t<decltype(expr1)>>;
	using DimArray     = typename ProdAnalyzer::DimArray;
	
	EXPECT_EQ(ProdAnalyzer::getLength(), 6);
	
	constexpr DimArray dims = ProdAnalyzer::getDims();
	
	EXPECT_EQ(dims, DimArray({1, 3, 1, 4, 1, 3, 1}));
	EXPECT_EQ(ProdAnalyzer::getOptimalCost(), 12);
	EXPECT_EQ(ProdAnalyzer::getOptimalSplit(), 2);
	
	using ExpectedExpr = decltype( (FSLinalg::transpose(a)*a)*((FSLinalg::transpose(e0)*e0)*(FSLinalg::transpose(a)*a)) );
	
	constexpr bool b1 = std::is_same<ExpectedExpr, typename ProdAnalyzer::OptimalBracketing>::value;
	
	EXPECT_TRUE(b1);
	
	using ReBracketType = typename ProdAnalyzer::ReBracketType;
	
	ReBracketType rebrackedExpr = ProdAnalyzer::reBracket(expr1);
	
	const FSLinalg::RealMatrix<1,1> result(rebrackedExpr);
	const FSLinalg::RealMatrix<1,1> expected({169.});
	
	EXPECT_EQ(result, expected);
}

TEST(chain, general)
{
	const FSLinalg::RealMatrix<12, 3> A({
		{2, 7, 0},
		{3, 3, 4},
		{9, 3, 7},
		{9, 6, 4},
		{7, 1, 1},
		{9, 0, 4},
		{8, 8, 7},
		{1, 1, 9},
		{4, 1, 3},
		{7, 4, 2},
		{4, 2, 2},
		{2, 3, 8}});
	
	const FSLinalg::RealMatrix<3, 8> B({
		{9, 0, 4, 4, 2, 0, 6, 5},
		{7, 0, 6, 4, 1, 4, 1, 3},
		{3, 9, 2, 2, 5, 4, 2, 9}});
		
	const FSLinalg::RealMatrix<8, 5> C({
		{8, 1, 4, 3, 5},
		{0, 3, 0, 7, 2},
		{5, 2, 3, 1, 9},
		{1, 7, 5, 0, 5},
		{7, 7, 3, 8, 2},
		{1, 6, 6, 1, 2},
		{6, 4, 1, 7, 3},
		{2, 4, 5, 4, 8}});
		
	const FSLinalg::RealMatrix<5, 2> D({
		{2, 4},
		{1, 8},
		{6, 1},
		{2, 9},
		{1, 9}});
		
	const FSLinalg::RealMatrix<12,2> expected({
		{11504, 30045},
		{14138, 40118},
		{27308, 78125},
		{26402, 73571},
		{12756, 35823},
		{18986, 54677},
		{32062, 89869},
		{16504, 49153},
		{11554, 33154},
		{18002, 49937},
		{11252, 31636},
		{18864, 54785}});
	
	const auto expr = A*B*C*D;
	const FSLinalg::RealMatrix<12,2> result = expr;
	
	EXPECT_EQ(expected, result);
	
	using ProdAnalyzer = FSLinalg::MatrixProductAnalyzer<std::decay_t<decltype(expr)>>;
	using DimArray     = typename ProdAnalyzer::DimArray;
	
	EXPECT_FALSE(expr.isOptimallyBracked());
	EXPECT_TRUE(ProdAnalyzer::reBracket(expr).isOptimallyBracked());
	
	// A*B*C*D // automatic re-bracketing
	// (A*B*C*D).reBracket()
	// keepBracketing(A*B*C*D)
	EXPECT_EQ(ProdAnalyzer::getLength(), 4);
	
	constexpr DimArray dims = ProdAnalyzer::getDims();
	
	EXPECT_EQ(dims, DimArray({12, 3, 8, 5, 2}));
	EXPECT_EQ(ProdAnalyzer::getOptimalCost(), 200);
	EXPECT_EQ(ProdAnalyzer::getOptimalSplit(), 1);
	
	using ReBracketType = typename ProdAnalyzer::ReBracketType;
	
	ReBracketType rebrackedExpr = ProdAnalyzer::reBracket(expr);
	
	const FSLinalg::RealMatrix<12,2> result2 = FSLinalg::keepBrackets(expr);
	const FSLinalg::RealMatrix<12,2> rebrackedResult = rebrackedExpr;
	
	EXPECT_EQ(expected, result2);
	EXPECT_EQ(expected, rebrackedResult);
}
