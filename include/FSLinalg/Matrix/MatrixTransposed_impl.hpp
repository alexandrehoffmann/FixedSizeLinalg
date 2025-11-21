#ifndef FSLINALG_MATRIX_TRANSPOSE_IMPL_HPP
#define FSLINALG_MATRIX_TRANSPOSE_IMPL_HPP

#include <FSLinalg/Matrix/MatrixTransposed.hpp>

namespace FSLinalg
{
	
template<typename Expr> 
auto MatrixTransposed<Expr>::operator*=(const RealScalar& alpha) -> MatrixTransposed& requires(isScalarComplex)
{
	if constexpr (hasFlatRandomAccess)
	{
		for (Size i=0; i!=size; ++i) { m_expr[i] *= alpha; }
	}
	else
	{
		for (Size i=0; i!=nRows; ++i)
		{
			for (Size j=0; j!=nRows; ++j)
			{
				m_expr(i,j) *= alpha;
			}
		}
	}
}

template<typename Expr> 
auto MatrixTransposed<Expr>::operator/=(const RealScalar& alpha) -> MatrixTransposed& requires(isScalarComplex)
{
	if constexpr (hasFlatRandomAccess)
	{
		for (Size i=0; i!=size; ++i) { m_expr[i] /= alpha; }
	}
	else
	{
		for (Size i=0; i!=nRows; ++i)
		{
			for (Size j=0; j!=nRows; ++j)
			{
				m_expr(i,j) /= alpha;
			}
		}
	}
}

template<typename Expr> 
auto MatrixTransposed<Expr>::operator*=(const Scalar& alpha) -> MatrixTransposed&
{
	if constexpr (hasFlatRandomAccess)
	{
		for (Size i=0; i!=size; ++i) { m_expr[i] *= alpha; }
	}
	else
	{
		for (Size i=0; i!=nRows; ++i)
		{
			for (Size j=0; j!=nRows; ++j)
			{
				m_expr(i,j) *= alpha;
			}
		}
	}
}

template<typename Expr> 
auto MatrixTransposed<Expr>::operator/=(const Scalar& alpha) -> MatrixTransposed&
{
	if constexpr (hasFlatRandomAccess)
	{
		for (Size i=0; i!=size; ++i) { m_expr[i] /= alpha; }
	}
	else
	{
		for (Size i=0; i!=nRows; ++i)
		{
			for (Size j=0; j!=nRows; ++j)
			{
				m_expr(i,j) /= alpha;
			}
		}
	}
}

template<typename Expr>  template<typename Bool, typename Alpha, class Dst>
void MatrixTransposed<Expr>::assignToImpl(const Bool /* checkAliasing */, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	Matrix<Scalar, nCols, nRows> tmp(m_expr);
	
	for (Size i=0; i!=nRows; ++i)
	{
		for (Size j=0; j!=nRows; ++j)
		{
			dst(i,j) = alpha*tmp(j,i);
		}
	}
}
	
template<typename Expr>  template<typename Bool, typename Alpha, class Dst>
void MatrixTransposed<Expr>::incrementImpl(const Bool /* checkAliasing */, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	Matrix<Scalar, nCols, nRows> tmp(m_expr);
	
	for (Size i=0; i!=nRows; ++i)
	{
		for (Size j=0; j!=nRows; ++j)
		{
			dst(i,j) += alpha*tmp(j,i);
		}
	}
}
	
template<typename Expr> template<typename Bool, typename Alpha, class Dst>
void MatrixTransposed<Expr>::decrementImpl(const Bool /* checkAliasing */, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	Matrix<Scalar, nCols, nRows> tmp(m_expr);
	
	for (Size i=0; i!=nRows; ++i)
	{
		for (Size j=0; j!=nRows; ++j)
		{
			dst(i,j) -= alpha*tmp(j,i);
		}
	}
}

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_TRANSPOSE_IMPL_HPP
