#ifndef FSLINALG_MISC_NESTED_LOOP_HPP
#define FSLINALG_MISC_NESTED_LOOP_HPP

#include <cstddef>
#include <array>

#include <BIC/Core.hpp>

namespace FSLinalg
{
namespace misc
{
	
template<size_t d, size_t dim>
struct NestedLoop
{
	static_assert(d <= dim);
	
	template<typename Size, size_t rank, class UnaryFunc> requires(rank >= dim)
	static constexpr UnaryFunc&& run(const std::array<Size, rank>& shape, std::array<Size, rank>& index, UnaryFunc&& func);
};  

template<typename Size, Size dim, class UnaryFunc>
constexpr UnaryFunc&& nestedLoop(const std::array<Size, dim>& shape, UnaryFunc&& func);

} // namespace misc
} // namespace FSLinalg

#include <FSLinalg/misc/NestedLoop_impl.hpp>

#endif // FSLINALG_MISC_NESTED_LOOP_HPP
