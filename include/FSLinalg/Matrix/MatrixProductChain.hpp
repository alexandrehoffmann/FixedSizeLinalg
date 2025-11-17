#ifndef FSLINALG_MATRIX_PRODUCT_CHAIN_HPP
#define FSLINALG_MATRIX_PRODUCT_CHAIN_HPP

#include <array>

namespace FSLinalg
{

template<std::array dims> 
struct MatrixProductChain
{
	template<size_t i, size_t j>
	static constexpr std::pair<size_t, size_t> minMulCostAndSplitRec(std::integral_constant<size_t, i>, std::integral_constant<size_t, j>);
	
	template<size_t i, size_t j> 
	static constexpr std::pair<size_t, size_t> minCostAndSplit = minMulCostAndSplitRec(std::integral_constant<size_t, i>{}, std::integral_constant<size_t, j>{});
};

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_PRODUCT_CHAIN_HPP
