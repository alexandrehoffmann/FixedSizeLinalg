#ifndef FSLINALG_VECTOR_OUTER_PRODUCT_IMPL_HPP
#define FSLINALG_VECTOR_OUTER_PRODUCT_IMPL_HPP

#include <FSLinalg/Matrix/VectorOuterProduct.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> template<typename Alpha, class Dst, bool checkAliasing>
void VectorOuterProduct<Lhs,Rhs>::assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	TmpLhs lhs(m_lhs);
	TmpRhs rhs(m_rhs);
	
	if (checkAliasing and (lhs.isAliasedTo(dst) or rhs.isAliasedTo(dst)))
	{
		Matrix<typename Dst::Scalar, Dst::nRows, Dst::nCols> tmp_dst;
		
		for (Size i=0; i!=Dst::nRows; ++i) { for (Size j=0; j!=Dst::nCols; ++j) { tmp_dst(i,j) = lhs[i]*rhs[i]; } }
		for (Size i=0; i!=size; ++i) { dst[i] = alpha*tmp_dst[i]; }
	}
	else
	{
		for (Size i=0; i!=Dst::nRows; ++i) { for (Size j=0; j!=Dst::nCols; ++j) { dst(i,j) = alpha*lhs[i]*rhs[i]; } }
	}
}

template<class Lhs, class Rhs> template<typename Alpha, class Dst, bool checkAliasing>
void VectorOuterProduct<Lhs,Rhs>::incrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	TmpLhs lhs(m_lhs);
	TmpRhs rhs(m_rhs);
	
	if (checkAliasing and (lhs.isAliasedTo(dst) or rhs.isAliasedTo(dst)))
	{
		Matrix<typename Dst::Scalar, Dst::nRows, Dst::nCols> tmp_dst;
		
		for (Size i=0; i!=Dst::nRows; ++i) { for (Size j=0; j!=Dst::nCols; ++j) { tmp_dst(i,j) = lhs[i]*rhs[i]; } }
		for (Size i=0; i!=size; ++i) { dst[i] += alpha*tmp_dst[i]; }
	}
	else
	{
		for (Size i=0; i!=Dst::nRows; ++i) { for (Size j=0; j!=Dst::nCols; ++j) { dst(i,j) += alpha*lhs[i]*rhs[i]; } }
	}
}

template<class Lhs, class Rhs> template<typename Alpha, class Dst, bool checkAliasing>
void VectorOuterProduct<Lhs,Rhs>::decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	TmpLhs lhs(m_lhs);
	TmpRhs rhs(m_rhs);
	
	if (checkAliasing and (lhs.isAliasedTo(dst) or rhs.isAliasedTo(dst)))
	{
		Matrix<typename Dst::Scalar, Dst::nRows, Dst::nCols> tmp_dst;
		
		for (Size i=0; i!=Dst::nRows; ++i) { for (Size j=0; j!=Dst::nCols; ++j) { tmp_dst(i,j) = lhs[i]*rhs[i]; } }
		for (Size i=0; i!=size; ++i) { dst[i] -= alpha*tmp_dst[i]; }
	}
	else
	{
		for (Size i=0; i!=Dst::nRows; ++i) { for (Size j=0; j!=Dst::nCols; ++j) { dst(i,j) -= alpha*lhs[i]*rhs[i]; } }
	}
}

} // namespace FSLinalg

#endif // FSLINALG_VECTOR_OUTER_PRODUCT_IMPL_HPP
