#ifndef FSLINALG_MATRIX_VECTOR_PRODUCT_IMPL_HPP
#define FSLINALG_MATRIX_VECTOR_PRODUCT_IMPL_HPP

#include <FSLinalg/Vector/MatrixVectorProduct.hpp>

#include <FSLinalg/BasicLinalg/GeneralMatrixVectorProduct.hpp>

namespace FSLinalg
{
	
template<class Lhs, class Rhs> 
template<typename Beta, class Dst, bool checkAliasing>
void MatrixVectorProduct<Lhs,Rhs>::assignToImpl(const Beta& beta, VectorBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Beta>::value)
{
	using GemvAssign = BasicLinalg::GeneralMatrixVectorProduct<isMatTransposed, isMatConjugated, MatNRows, MatNCols, isVecConjugated, false>;
	
	StripSymbolsAndEvalMatrix<Lhs> strippedLhs(m_matrix);
	StripSymbolsAndEvalVector<Rhs> strippedRhs(m_vector);
	
	const auto alpha     = beta*strippedLhs.getAlpha()*strippedRhs.getAlpha();
	const StrippedMat& A = strippedLhs.getMatrix();
	const StrippedVec& x = strippedRhs.getVector();
	
	if (checkAliasing and (A.isAliasedTo(dst) or x.isAliasedTo(dst)))
	{
		Vector<typename Dst::Scalar, Dst::size> tmp(dst);
		GemvAssign::run(alpha, A, x, tmp);
		for (Size i=0; i!=size; ++i) { dst[i] = tmp[i]; }
	}
	else
	{
		GemvAssign::run(alpha, A, x, dst.derived());
	}
}

template<class Lhs, class Rhs> 
template<typename Beta, class Dst, bool checkAliasing>
void MatrixVectorProduct<Lhs,Rhs>::incrementImpl(const Beta& beta, VectorBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Beta>::value)
{
	using GemvAssign    = BasicLinalg::GeneralMatrixVectorProduct<isMatTransposed, isMatConjugated, MatNRows, MatNCols, isVecConjugated, false>;
	using GemvIncrement = BasicLinalg::GeneralMatrixVectorProduct<isMatTransposed, isMatConjugated, MatNRows, MatNCols, isVecConjugated, true>;
	
	StripSymbolsAndEvalMatrix<Lhs> strippedLhs(m_matrix);
	StripSymbolsAndEvalVector<Rhs> strippedRhs(m_vector);
	
	const auto alpha     = beta*strippedLhs.getAlpha()*strippedRhs.getAlpha();
	const StrippedMat& A = strippedLhs.getMatrix();
	const StrippedVec& x = strippedRhs.getVector();
	
	if (checkAliasing and (A.isAliasedTo(dst) or x.isAliasedTo(dst)))
	{
		Vector<typename Dst::Scalar, Dst::size> tmp(dst);
		GemvAssign::run(alpha, A, x, tmp);
		for (Size i=0; i!=size; ++i) { dst[i] += tmp[i]; }
	}
	else
	{
		GemvIncrement::run(alpha, A, x, dst.derived());
	}
}

template<class Lhs, class Rhs> 
template<typename Beta, class Dst, bool checkAliasing>
void MatrixVectorProduct<Lhs,Rhs>::decrementImpl(const Beta& beta, VectorBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Beta>::value)
{
	using GemvAssign    = BasicLinalg::GeneralMatrixVectorProduct<isMatTransposed, isMatConjugated, MatNRows, MatNCols, isVecConjugated, false>;
	using GemvIncrement = BasicLinalg::GeneralMatrixVectorProduct<isMatTransposed, isMatConjugated, MatNRows, MatNCols, isVecConjugated, true>;
	
	StripSymbolsAndEvalMatrix<Lhs> strippedLhs(m_matrix);
	StripSymbolsAndEvalVector<Rhs> strippedRhs(m_vector);
	
	const auto alpha     = beta*strippedLhs.getAlpha()*strippedRhs.getAlpha();
	const StrippedMat& A = strippedLhs.getMatrix();
	const StrippedVec& x = strippedRhs.getVector();
	
	if (checkAliasing and (A.isAliasedTo(dst) or x.isAliasedTo(dst)))
	{
		Vector<typename Dst::Scalar, Dst::size> tmp(dst);
		GemvAssign::run(alpha, A, x, tmp);
		for (Size i=0; i!=size; ++i) { dst[i] -= tmp[i]; }
	}
	else
	{
		GemvIncrement::run(-alpha, A, x, dst.derived());
	}
}
	
} // namespace FSLinalg

#endif // FSLINALG_MATRIX_VECTOR_PRODUCT_IMPL_HPP
