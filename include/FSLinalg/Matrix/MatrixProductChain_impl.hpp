#ifndef FSLINALG_MATRIX_PRODUCT_CHAIN_IMPL_HPP
#define FSLINALG_MATRIX_PRODUCT_CHAIN_IMPL_HPP

#include <FSLinalg/Matrix/MatrixProductChain.hpp>

namespace FSLinalg
{

template<std::array dims> template<size_t I, size_t J>
constexpr std::pair<size_t, size_t> MatrixProductChain<dims>::minMulCostAndSplitRec(BIC::Fixed<size_t, I> i, BIC::Fixed<size_t, J> j)
{
	if constexpr (i + 1 == j) 
	{ 
		return std::make_pair(0u, i + 1);
	}
	else
	{
		size_t minCost     = std::numeric_limits<size_t>::max();
		size_t optSpliting = i + 1; 
		
		BIC::foreach(BIC::next(i), j, [i, j, &minCost, &optSpliting](const auto k) -> void
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
