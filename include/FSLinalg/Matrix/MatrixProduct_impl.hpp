#ifndef FSLINALG_MATRIX_PRODUCT_IMPL_HPP
#define FSLINALG_MATRIX_PRODUCT_IMPL_HPP

#include <FSLinalg/Matrix/MatrixProduct.hpp>
#include <FSLinalg/BasicLinalg/GeneralMatrixMatrixProduct.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> 
template<typename Alpha, class Dst, bool checkAliasing>
void MatrixProduct<Lhs,Rhs>::assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	using GemmAssign = BasicLinalg::GeneralMatrixMatrixProduct<isLhsTransposed, isLhsConjugated, LhsNRows, LhsNCols, isRhsTransposed, isRhsConjugated, RhsNRows, RhsNCols, false>;
	
	StripSymbolsAndEvalMatrix<Lhs> strippedLhs(m_lhs);
	StripSymbolsAndEvalMatrix<Rhs> strippedRhs(m_rhs);
	
	const auto beta      = alpha*strippedLhs.getAlpha()*strippedRhs.getAlpha();
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

template<class Lhs, class Rhs> 
template<typename Alpha, class Dst, bool checkAliasing>
void MatrixProduct<Lhs,Rhs>::incrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
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

template<class Lhs, class Rhs> 
template<typename Alpha, class Dst, bool checkAliasing>
void MatrixProduct<Lhs,Rhs>::decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
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

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_PRODUCT_IMPL_HPP
