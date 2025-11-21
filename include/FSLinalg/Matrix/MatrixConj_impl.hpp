#ifndef FSLINALG_MATRIX_CONJ_IMPL_HPP
#define FSLINALG_MATRIX_CONJ_IMPL_HPP

#include <FSLinalg/Matrix/MatrixConj.hpp>

namespace FSLinalg
{

template<class Expr> template<typename Bool, typename Alpha, class Dst>
void MatrixConj<Expr>::assignToImpl(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	m_expr.assignTo(checkAliasing, conj(alpha), dst);
	if constexpr (IsComplexScalar<Scalar>::value)
	{
		if constexpr (Dst::hasFlatRandomAccess)
		{
			for (Size i=0; i!=size; ++i)
			{
				conjInPlace(dst[i]);
			}
		}
		else
		{
			for (Size i=0; i!=nRows; ++i)
			{
				for (Size j=0; j!=nRows; ++j)
				{
					conjInPlace(dst(i, j));
				}
			}
		}
	}
}
	
template<class Expr> template<typename Bool, typename Alpha, class Dst>
void MatrixConj<Expr>::incrementImpl(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	if constexpr (IsRealScalar<Scalar>::value)
	{
		m_expr.increment(checkAliasing, alpha, dst);
	}
	else
	{
		Matrix<Scalar, nRows, nCols> tmp(m_expr);
		if constexpr (Dst::hasFlatRandomAccess)
		{
			for (Size i=0; i!=size; ++i)
			{
				dst[i] += alpha*conj(tmp[i]);
			}
		}
		else
		{
			for (Size i=0; i!=nRows; ++i)
			{
				for (Size j=0; j!=nRows; ++j)
				{
					dst(i,j) += alpha*conj(tmp(i,j));
				}
			}
		}
	}
}
	
template<class Expr> template<typename Bool, typename Alpha, class Dst>
void MatrixConj<Expr>::decrementImpl(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	if constexpr (IsRealScalar<Scalar>::value)
	{
		m_expr.decrement(checkAliasing, alpha, dst);
	}
	else
	{
		Matrix<Scalar, nRows, nCols> tmp(m_expr);
		if constexpr (Dst::hasFlatRandomAccess)
		{
			for (Size i=0; i!=size; ++i)
			{
				dst[i] -= alpha*conj(tmp[i]);
			}
		}
		else
		{
			for (Size i=0; i!=nRows; ++i)
			{
				for (Size j=0; j!=nRows; ++j)
				{
					dst(i,j) -= alpha*conj(tmp(i,j));
				}
			}
		}
	}
}

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_CONJ_IMPL_HPP
