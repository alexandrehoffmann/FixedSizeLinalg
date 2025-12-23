#ifndef FSLINALG_MISC_NESTED_LOOP_IMPL_HPP
#define FSLINALG_MISC_NESTED_LOOP_IMPL_HPP

#include <FSLinalg/misc/NestedLoop.hpp>

namespace FSLinalg
{
namespace misc
{
	
template<size_t d, size_t dim> template<typename Size, size_t rank, class UnaryFunc> requires(rank >= dim)
constexpr UnaryFunc&& NestedLoop<d,dim>::run(const std::array<Size, rank>& shape, std::array<Size, rank>& index, UnaryFunc&& func)
{
	if constexpr (d == dim)
	{
		func(std::as_const(index));
	}
	else
	{
		for (index[d]=0; index[d] != shape[d]; ++index[d])
		{
			NestedLoop<d+1,dim>::run(shape, index, std::forward<UnaryFunc>(func));
		}
	}
	return std::forward<UnaryFunc>(func);
}  

template<typename Size, size_t dim, class UnaryFunc>
constexpr UnaryFunc&& nestedLoop(const std::array<Size, dim>& shape, UnaryFunc&& func)
{
	std::array<Size, dim> index;
	return NestedLoop<0, dim>::run(shape, index, std::forward<UnaryFunc>(func));
}

} // namespace misc
} // namespace FSLinalg

#endif // FSLINALG_MISC_NESTED_LOOP_IMPL_HPP
