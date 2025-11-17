#ifndef FSLINALG_MISC_STATIC_FOR_HPP
#define FSLINALG_MISC_STATIC_FOR_HPP

namespace FSLinalg
{
namespace misc
{
	
	template<size_t i, size_t n>
	struct StaticFor
	{
		static_assert(i < n);
		
		template<class UnaryFunc>
		static constexpr UnaryFunc&& run(UnaryFunc&& func) { func(std::integral_constant<size_t, i>{}); return StaticFor<i+1,n>::run(std::forward<UnaryFunc>(func)); }
	};
	
	template<size_t n>
	struct StaticFor<n,n>
	{
		template<class UnaryFunc>
		static constexpr UnaryFunc&& run(UnaryFunc&& func) { return std::forward<UnaryFunc>(func); }
	};
	
} // namespace misc
} // namespace FSLinalg

#endif // FSLINALG_MISC_STATIC_FOR_HPP
