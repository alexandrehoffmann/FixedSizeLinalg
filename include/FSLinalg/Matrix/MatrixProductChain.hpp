#ifndef FSLINALG_MATRIX_PRODUCT_CHAIN_HPP
#define FSLINALG_MATRIX_PRODUCT_CHAIN_HPP

#include <array>

#include <BIC/Core.hpp>

namespace FSLinalg
{

template<std::array dims> 
struct MatrixProductChain
{
	template<size_t i, size_t j>
	static constexpr std::pair<size_t, size_t> minMulCostAndSplitRec(BIC::Fixed<size_t, i> fixed_i, BIC::Fixed<size_t, j> fixed_j);
	
	template<size_t i, size_t j> 
	static constexpr std::pair<size_t, size_t> minCostAndSplit = minMulCostAndSplitRec(BIC::fixed<size_t, i>, BIC::fixed<size_t, j>);
};

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_PRODUCT_CHAIN_HPP
