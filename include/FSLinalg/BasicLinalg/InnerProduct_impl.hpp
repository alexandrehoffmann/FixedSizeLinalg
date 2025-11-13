#ifndef FSLINALG_INNER_PRODUCT_IMPL_HPP
#define FSLINALG_INNER_PRODUCT_IMPL_HPP

#include <FSLinalg/Scalar.hpp>
#include <FSLinalg/BasicLinalg/Product.hpp>
#include <FSLinalg/BasicLinalg/InnerProduct.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs>
InnerProductScalar<Lhs,Rhs> inner(const MatrixBase<Lhs>& base_lhs, const MatrixBase<Rhs>& base_rhs)
{
	static_assert(Lhs::nRows == Rhs::nCols, "Matrices sizes must match");
	static_assert(Lhs::nCols == Rhs::nCols, "Matrices sizes must match");
	
	using TmpLhs = std::conditional_t<Lhs::hasReadRandomAccess, const Lhs&, Matrix<typename Lhs::Scalar, Lhs::nRows, Lhs::nCols> >;
	using TmpRhs = std::conditional_t<Rhs::hasReadRandomAccess, const Rhs&, Matrix<typename Rhs::Scalar, Rhs::nRows, Rhs::nCols> >;
	
	using Size   = std::common_type_t<typename Lhs::Size, typename Rhs::Size>;
	using Scalar = InnerProductScalar<Lhs,Rhs>;
	
	constexpr BasicLinalg::Product<true,false> prod;
	
	TmpLhs lhs(base_lhs.derived());
	TmpRhs rhs(base_rhs.derived());
	
	Scalar res(0);
	
	if constexpr (lhs.hasFlatRandomAccess and rhs.hasFlatRandomAccess)
	{
		for (Size i=0; i!=Lhs::size; ++i)
		{
			res += prod(lhs[i], rhs[i]);
		}
	}
	else
	{
		for (Size i=0; i!=lhs.nRows; ++i)
		{
			for (Size j=0; j!=rhs.nCols; ++j)
			{
				res += prod(lhs(i,j), rhs(i,j));
			}
		}
	}
	
	return res;
}
	
} // namespace FSLinalg

#endif // FSLINALG_INNER_PRODUCT_IMPL_HPP
