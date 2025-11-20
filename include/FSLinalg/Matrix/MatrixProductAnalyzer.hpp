#ifndef FSLINALG_MATRIX_PRODUCT_TRAITS_HPP
#define FSLINALG_MATRIX_PRODUCT_TRAITS_HPP

#include <cstddef>
#include <type_traits>

#include <FSLinalg/Matrix/MatrixProductChain.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> class MatrixProduct;
	
namespace detail
{

template<class Expr>
struct MatrixProductAnalyzerImpl
{
	static constexpr size_t length = 1;
	
	template<size_t n> requires(n == 0)
	using NthMatrix = Expr;
	
	template<size_t n>
	static constexpr const NthMatrix<n>& getMatrix(const Expr& expr, std::integral_constant<size_t, n>) { return expr; }
	
	template<size_t idx, size_t chainLenP1> requires (idx+1 <chainLenP1) 
	static constexpr void fillDims(std::integral_constant<size_t, idx>, std::array<size_t, chainLenP1>& dims) { dims[idx] = Expr::nRows; dims[idx+1] = Expr::nCols; }
};

template<class Lhs, class Rhs>
struct MatrixProductAnalyzerImpl< MatrixProduct<Lhs, Rhs> >
{
	static constexpr size_t lhsLength = MatrixProductAnalyzerImpl<Lhs>::length;
	static constexpr size_t rhsLength = MatrixProductAnalyzerImpl<Rhs>::length;
	static constexpr size_t length    = lhsLength + rhsLength;
	
	template<size_t n, class Enable = void>
	struct NthMatrixHelper;
	
	template<size_t n> 
	struct NthMatrixHelper<n, std::enable_if_t<n < lhsLength>>
	{
		using Type = typename MatrixProductAnalyzerImpl<Lhs>::template NthMatrix<n>;
	};
	
	template<size_t n> 
	struct NthMatrixHelper<n, std::enable_if_t<lhsLength <= n and n < length>>
	{
		using Type = typename MatrixProductAnalyzerImpl<Rhs>::template NthMatrix<n - lhsLength>;
	};

	template<size_t n> 
	using NthMatrix = typename NthMatrixHelper<n>::Type;
	
	template<size_t n>
	static constexpr const NthMatrix<n>& getMatrix(const MatrixProduct<Lhs, Rhs>& expr, std::integral_constant<size_t, n>);

	template<size_t idx, size_t chainLenP1> requires (idx+1 <chainLenP1) 
	static constexpr void fillDims(std::integral_constant<size_t, idx>, std::array<size_t, chainLenP1>& dims);
};

} // namespace detail

template<class Expr>
struct MatrixProductAnalyzer
{
	static_assert(IsMatrix<Expr>::value, "Expr must be a matrix");
	
	using Impl = detail::MatrixProductAnalyzerImpl<Expr>;
	using DimArray = std::array<size_t, Impl::length+1>;
	
	template<size_t n> using NthMatrix = typename Impl::template NthMatrix<n>;
	
	template<size_t n>
	static const NthMatrix<n>& getMatrix(const Expr& expr, std::integral_constant<size_t, n> ic) { return Impl::getMatrix(expr, ic); }
	
	static constexpr size_t getLength() { return Impl::length; }
	
	static constexpr DimArray getDims();
	
	// we use an external template class as to not recompute everything for every product
	// if the optimal splitting for an chain with the same dims has already been computed 
	// we can re-use it.
	static constexpr size_t getOptimalCost  () { return MatrixProductChain<getDims()>::template minCostAndSplit<0, getLength()>.first;  }
	static constexpr size_t getOptimalSplit () { return MatrixProductChain<getDims()>::template minCostAndSplit<0, getLength()>.second; }
private:
	template<size_t start, size_t end> requires(start <= end and end <= getLength())
	struct OptimalBracketingHelper
	{
		static constexpr size_t split = MatrixProductChain<getDims()>::template minCostAndSplit<start, end>.second;
		
		static_assert(start <= split and split+1 < end+1, "invalid split");
		
		using LhsBracketing = OptimalBracketingHelper<start, split>;
		using RhsBracketing = OptimalBracketingHelper<split, end>;
		
		using Lhs = typename LhsBracketing::Type;
		using Rhs = typename RhsBracketing::Type;
		
		using Type          = MatrixProduct<Lhs, Rhs>;
		using ReBracketType = Type; 
		
		static ReBracketType reBracket(const Expr& expr) { return ReBracketType(LhsBracketing::reBracket(expr), RhsBracketing::reBracket(expr)); }
	};
	
	template<size_t idx> requires(idx < getLength())
	struct OptimalBracketingHelper<idx,idx+1>
	{
		using Type          = NthMatrix<idx>;
		using ReBracketType = const Type&;
		
		static ReBracketType reBracket(const Expr& expr) { return MatrixProductAnalyzer<Expr>::getMatrix(expr, std::integral_constant<size_t, idx>{}); }
	};
public:
	using OptimalBracketing = typename OptimalBracketingHelper<0, getLength()>::Type;
	using ReBracketType     = typename OptimalBracketingHelper<0, getLength()>::ReBracketType;
	
	static ReBracketType reBracket(const Expr& expr) { return OptimalBracketingHelper<0, getLength()>::reBracket(expr); }
};

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_PRODUCT_TRAITS_HPP
