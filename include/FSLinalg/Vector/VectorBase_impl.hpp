#ifndef FSLINALG_VECTOR_BASE_IMPL_HPP
#define FSLINALG_VECTOR_BASE_IMPL_HPP

#include <FSLinalg/Vector/VectorBase.hpp>
#include <FSLinalg/Vector/Vector.hpp>

namespace FSLinalg
{

template<class Derived> template<typename Alpha, class Dst, bool checkAliasing>
void VectorBase<Derived>::assignTo(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::size == size, "Vector sizes must match");
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			Vector<Scalar, size> tmp(*this);
			
			for (Size i=0; i!=getSize(); ++i) { dst[i] = alpha*tmp[i]; }
		}
		else
		{
			for (Size i=0; i!=getSize(); ++i) {	dst[i] = alpha*Base::derived().getImpl(i); }
		}
	}
	else
	{
		Base::derived().assignToImpl(alpha, dst, bc);
	}
}

template<class Derived> template<typename Alpha, class Dst, bool checkAliasing>
void VectorBase<Derived>::increment(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::size == size, "Vector sizes must match");
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			Vector<Scalar, size> tmp(*this);

			for (Size i=0; i!=getSize(); ++i) { dst[i] += alpha*tmp[i]; }
		}
		else
		{
			for (Size i=0; i!=getSize(); ++i) { dst[i] += alpha*Base::derived().getImpl(i); }
		}
	}
	else
	{
		Base::derived().incrementImpl(alpha, dst, bc);
	}
}

template<class Derived> template<typename Alpha, class Dst, bool checkAliasing>
void VectorBase<Derived>::decrement(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::size == size, "Vector sizes must match");
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			Vector<Scalar, size> tmp(*this);

			for (Size i=0; i!=getSize(); ++i) {	dst[i] -= alpha*tmp[i]; }
		}
		else
		{
			for (Size i=0; i!=getSize(); ++i) {	dst[i] -= alpha*Base::derived().getImpl(i);	}
		}
	}
	else
	{
		Base::derived().decrementImpl(alpha, dst, bc);
	}
}

template<typename Lhs, typename Rhs> requires(Lhs::hasReadRandomAccess and Rhs::hasReadRandomAccess) 
bool operator==(const FSLinalg::VectorBase<Lhs>& lhs, const FSLinalg::VectorBase<Rhs>& rhs)
{
	static_assert(Lhs::size == Rhs::size);
	
	using Size = std::common_type_t<typename Lhs::Size, typename Rhs::Size>;
	
	bool areEquals = true;
	for (Size i=0; i!=Lhs::size; ++i)
	{
		areEquals = areEquals and (lhs[i] == rhs[i]);
	}
	return areEquals;
}

} // namespace FSLinalg

#endif // FSLINALG_VECTOR_BASE_IMPL_HPP
