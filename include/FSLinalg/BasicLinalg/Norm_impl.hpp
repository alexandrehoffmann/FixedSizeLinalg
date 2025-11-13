#ifndef FSLINALG_NORM_IMPL_HPP
#define FSLINALG_NORM_IMPL_HPP

#include <FSLinalg/BasicLinalg/Norm.hpp>
#include <FSLinalg/Scalar.hpp>

namespace FSLinalg
{

template<typename Expr>
typename Expr::RealScalar squaredNorm(const VectorBase<Expr>& base_expr)
{
	using TmpExpr    = std::conditional_t<Expr::hasReadRandomAccess, const Expr&, Vector<typename Expr::Scalar, Expr::size> >;
	using Size       = typename Expr::Size;
	using RealScalar = typename Expr::RealScalar;
	
	TmpExpr expr(base_expr.derived());
	
	RealScalar res(0);
	
	for (Size i=0; i!=expr.size; ++i)
	{
		res += abs2(expr[i]);
	}
	
	return res;
}

template<typename Expr>
typename Expr::RealScalar squaredNorm(const MatrixBase<Expr>& base_expr)
{
	using TmpExpr    = std::conditional_t<Expr::hasReadRandomAccess, const Expr&, MatrixBase<typename Expr::Scalar, Expr::nRows, Expr::nCols> >;
	using Size       = typename Expr::Size;
	using RealScalar = typename Expr::RealScalar;
	
	TmpExpr expr(base_expr.derived());
	
	RealScalar res(0);
	
	if constexpr (expr.hasFlatRandomAccess)
	{
		for (Size i=0; i!=expr.size; ++i)
		{
			res += abs2(expr[i]);
		}
	}
	else
	{
		for (Size i=0; i!=expr.nRows; ++i)
		{
			for (Size j=0; j!=expr.nCols; ++j)
			{
				res += abs2(expr(i,j));
			}
		}
	}
	
	return res;
}

} // namespace FSLinalg

#endif // FSLINALG_NORM_IMPL_HPP
