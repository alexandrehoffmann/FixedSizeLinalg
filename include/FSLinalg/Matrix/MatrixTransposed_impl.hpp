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

template<typename Expr>  template<typename Alpha, class Dst, bool checkAliasing>
void MatrixTransposed<Expr>::assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
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
	
template<typename Expr>  template<typename Alpha, class Dst, bool checkAliasing>
void MatrixTransposed<Expr>::incrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
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
	
template<typename Expr> template<typename Alpha, class Dst, bool checkAliasing>
void MatrixTransposed<Expr>::decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
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
