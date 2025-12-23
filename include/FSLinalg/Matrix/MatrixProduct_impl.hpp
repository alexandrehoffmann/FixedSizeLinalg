#ifndef FSLINALG_MATRIX_PRODUCT_IMPL_HPP
#define FSLINALG_MATRIX_PRODUCT_IMPL_HPP

#include <FSLinalg/Matrix/MatrixProduct.hpp>
#include <FSLinalg/BasicLinalg/GeneralMatrixMatrixProduct.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> 
constexpr bool MatrixProduct<Lhs,Rhs>::isOptimallyBracked()
{
	using OptimallyBracketedSelf = typename MatrixProductAnalyzer<Self>::OptimalBracketing;
	return std::is_same<Self, OptimallyBracketedSelf>::value;
}

template<class Lhs, class Rhs> 
template<typename Bool, typename Alpha, class Dst, bool keepBracketing>
void MatrixProduct<Lhs,Rhs>::assignToHelper(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst, BIC::Fixed<bool, keepBracketing>) const 
	requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	if constexpr (isOptimallyBracked() or keepBracketing)
	{
		using GemmAssign = BasicLinalg::GeneralMatrixMatrixProduct<isLhsTransposed, isLhsConjugated, LhsNRows, LhsNCols, isRhsTransposed, isRhsConjugated, RhsNRows, RhsNCols, false>;
		
		StripSymbolsAndEvalMatrix<Lhs> strippedLhs(m_lhs);
		StripSymbolsAndEvalMatrix<Rhs> strippedRhs(m_rhs);
		
		const auto beta = alpha*strippedLhs.getAlpha()*strippedRhs.getAlpha();
		
		const StrippedLhs& A = strippedLhs.getMatrix();
		const StrippedRhs& B = strippedRhs.getMatrix();
		
		if (checkAliasing and (A.isAliasedTo(dst) or B.isAliasedTo(dst)))
		{
			Matrix<typename Dst::Scalar, Dst::nRows, Dst::nCols> tmp(dst);
			GemmAssign::run(beta, A, B, tmp);
			for (Size i=0; i!=size; ++i) { dst[i] = tmp[i]; }
		}
		else
		{
			GemmAssign::run(beta, A, B, dst.derived());
		}
	}
	else
	{
		MatrixProductAnalyzer<Self>::reBracket(*this).assignToImpl(checkAliasing, alpha, dst);
	}
}

template<class Lhs, class Rhs> 
template<typename Bool, typename Alpha, class Dst, bool keepBracketing>
void MatrixProduct<Lhs,Rhs>::incrementHelper(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst, BIC::Fixed<bool, keepBracketing>) const 
	requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{	
	if constexpr (isOptimallyBracked() or keepBracketing)
	{
		using GemmAssign    = BasicLinalg::GeneralMatrixMatrixProduct<isLhsTransposed, isLhsConjugated, LhsNRows, LhsNCols, isRhsTransposed, isRhsConjugated, RhsNRows, RhsNCols, false>;
		using GemmIncrement = BasicLinalg::GeneralMatrixMatrixProduct<isLhsTransposed, isLhsConjugated, LhsNRows, LhsNCols, isRhsTransposed, isRhsConjugated, RhsNRows, RhsNCols, true>;
		
		StripSymbolsAndEvalMatrix<Lhs> strippedLhs(m_lhs);
		StripSymbolsAndEvalMatrix<Rhs> strippedRhs(m_rhs);
		
		const auto beta      = alpha*strippedLhs.getAlpha()*strippedRhs.getAlpha();
		const StrippedLhs& A = strippedLhs.getMatrix();
		const StrippedRhs& B = strippedRhs.getMatrix();
		
		if (checkAliasing and (A.isAliasedTo(dst) or B.isAliasedTo(dst)))
		{
			Matrix<typename Dst::Scalar, Dst::nRows, Dst::nCols> tmp(dst);
			GemmAssign::run(beta, A, B, tmp);
			for (Size i=0; i!=size; ++i) { dst[i] += tmp[i]; }
		}
		else
		{
			GemmIncrement::run(beta, A, B, dst.derived());
		}
	}
	else
	{
		MatrixProductAnalyzer<Self>::reBracket(*this).incrementImpl(checkAliasing, alpha, dst);
	}
}

template<class Lhs, class Rhs> 
template<typename Bool, typename Alpha, class Dst, bool keepBracketing>
void MatrixProduct<Lhs,Rhs>::decrementHelper(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst, BIC::Fixed<bool, keepBracketing>) const 
	requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	if constexpr (isOptimallyBracked() or keepBracketing)
	{
		using GemmAssign    = BasicLinalg::GeneralMatrixMatrixProduct<isLhsTransposed, isLhsConjugated, LhsNRows, LhsNCols, isRhsTransposed, isRhsConjugated, RhsNRows, RhsNCols, false>;
		using GemmIncrement = BasicLinalg::GeneralMatrixMatrixProduct<isLhsTransposed, isLhsConjugated, LhsNRows, LhsNCols, isRhsTransposed, isRhsConjugated, RhsNRows, RhsNCols, true>;
		
		StripSymbolsAndEvalMatrix<Lhs> strippedLhs(m_lhs);
		StripSymbolsAndEvalMatrix<Rhs> strippedRhs(m_rhs);
		
		const auto beta     = alpha*strippedLhs.getAlpha()*strippedRhs.getAlpha();
		const StrippedLhs& A = strippedLhs.getMatrix();
		const StrippedRhs& B = strippedRhs.getMatrix();
		
		if (checkAliasing and (A.isAliasedTo(dst) or B.isAliasedTo(dst)))
		{
			Matrix<typename Dst::Scalar, Dst::nRows, Dst::nCols> tmp(dst);
			GemmAssign::run(beta, A, B, tmp);
			for (Size i=0; i!=size; ++i) { dst[i] -= tmp[i]; }
		}
		else
		{
			GemmIncrement::run(-beta, A, B, dst.derived());
		}
	}
	else
	{
		MatrixProductAnalyzer<Self>::reBracket(*this).decrementImpl(checkAliasing, alpha, dst);
	}
}

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_PRODUCT_IMPL_HPP
