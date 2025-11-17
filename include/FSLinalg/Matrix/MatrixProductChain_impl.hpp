#ifndef FSLINALG_MATRIX_PRODUCT_CHAIN_IMPL_HPP
#define FSLINALG_MATRIX_PRODUCT_CHAIN_IMPL_HPP

#include <FSLinalg/Matrix/MatrixProductChain.hpp>
#include <FSLinalg/misc/StaticFor.hpp>

namespace FSLinalg
{

template<std::array dims> template<size_t i, size_t j>
constexpr std::pair<size_t, size_t> MatrixProductChain<dims>::minMulCostAndSplitRec(std::integral_constant<size_t, i>, std::integral_constant<size_t, j>)
{
	if constexpr (i + 1 == j) 
	{ 
		return std::make_pair(0u, i + 1);
	}
	else
	{
		size_t minCost     = std::numeric_limits<size_t>::max();
		size_t optSpliting = i + 1; 
		
		misc::StaticFor<i + 1, j>::run([&minCost, &optSpliting](const auto k) -> void
		{
			constexpr size_t curr = minCostAndSplit<i, k>.first + minCostAndSplit<k, j>.first + dims[i]*dims[k]*dims[j];
			if (curr < minCost)
			{
				minCost = curr;
				optSpliting = k;
			}
		});
		
		return std::make_pair(minCost, optSpliting);
	}
}

} //namespace FSLinalg

#endif // FSLINALG_MATRIX_PRODUCT_CHAIN_IMPL_HPP
