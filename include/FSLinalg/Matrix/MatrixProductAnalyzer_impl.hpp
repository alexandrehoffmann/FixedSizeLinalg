#ifndef FSLINALG_MATRIX_PRODUCT_TRAITS_IMPL_HPP
#define FSLINALG_MATRIX_PRODUCT_TRAITS_IMPL_HPP

#include <FSLinalg/Matrix/MatrixProductAnalyzer.hpp>

namespace FSLinalg
{
namespace detail
{

template<class Lhs, class Rhs> template<size_t idx, size_t chainLenP1> requires (idx+1 <chainLenP1)
constexpr void MatrixProductAnalyzerImpl< MatrixProduct<Lhs, Rhs> >::fillDims(std::integral_constant<size_t, idx> ic, std::array<size_t, chainLenP1>& dims)
{
	if constexpr (IsMatrixProduct<Lhs>::value and IsMatrixProduct<Rhs>::value)
	{		
		MatrixProductAnalyzerImpl<Lhs>::fillDims(ic, dims);
		MatrixProductAnalyzerImpl<Rhs>::fillDims(std::integral_constant<size_t, idx + lhsLength>{}, dims);
	}
	else if constexpr (IsMatrixProduct<Lhs>::value)
	{	
		MatrixProductAnalyzerImpl<Lhs>::fillDims(ic, dims);	
		dims[idx + lhsLength + 1] = Rhs::nCols;
	}
	else if constexpr (IsMatrixProduct<Rhs>::value)
	{		
		dims[idx] = Lhs::nRows;
		MatrixProductAnalyzerImpl<Rhs>::fillDims(std::integral_constant<size_t, idx + 1>{}, dims);
	}
	else
	{		
		dims[idx]     = Lhs::nRows;
		dims[idx + 1] = Lhs::nCols;
		dims[idx + 2] = Rhs::nCols;	
	}
}

template<class Lhs, class Rhs> template<size_t n>
constexpr auto MatrixProductAnalyzerImpl< MatrixProduct<Lhs, Rhs> >::getMatrix(const MatrixProduct<Lhs, Rhs>& expr, std::integral_constant<size_t, n> ic) -> const NthMatrix<n>&
{
	if constexpr (n < lhsLength)
	{
		return MatrixProductAnalyzerImpl<Lhs>::getMatrix(expr.m_lhs, ic);
	}
	else
	{
		return MatrixProductAnalyzerImpl<Rhs>::getMatrix(expr.m_rhs, std::integral_constant<size_t, n - lhsLength>{});
	}
}

} // namespace detail

template<class Expr>
constexpr auto MatrixProductAnalyzer<Expr>::getDims() -> DimArray
{ 
	DimArray dims; 
	Impl::fillDims(std::integral_constant<size_t, 0>{}, dims); 
	return dims; 
}

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_PRODUCT_TRAITS_IMPL_HPP
