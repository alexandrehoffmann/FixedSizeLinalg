#ifndef FSLINALG_VECTOR_CONJ_IMPL_HPP
#define FSLINALG_VECTOR_CONJ_IMPL_HPP

#include <FSLinalg/Vector/VectorConj.hpp>
#include <FSLinalg/Vector/Vector.hpp>

namespace FSLinalg
{

template<class Expr> template<typename Alpha, class Dst, bool checkAliasing>
void VectorConj<Expr>::assignToImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	m_expr.assignTo(conj(alpha), dst, bc);
	if constexpr (IsComplexScalar<Scalar>::value)
	{
		for (Size i=0; i!=size; ++i)
		{
			conjInPlace(dst[i]);
		}
	}
}

template<class Expr> template<typename Alpha, class Dst, bool checkAliasing>
void VectorConj<Expr>::incrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	if constexpr (IsRealScalar<Scalar>::value)
	{
		m_expr.increment(alpha, dst, bc);
	}
	else
	{
		Vector<Scalar, size> tmp(m_expr);
		for (Size i=0; i!=size; ++i)
		{
			dst[i] += alpha*conj(tmp[i]);
		}
	}
}

template<class Expr> template<typename Alpha, class Dst, bool checkAliasing>
void VectorConj<Expr>::decrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	if constexpr (IsRealScalar<Scalar>::value)
	{
		m_expr.decrement(alpha, dst, bc);
	}
	else
	{
		Vector<Scalar, size> tmp(m_expr);
		for (Size i=0; i!=size; ++i)
		{
			dst[i] -= alpha*conj(tmp[i]);
		}
	}
}

} // namespace FSLinalg

#endif // FSLINALG_VECTOR_CONJ_IMPL_HPP
