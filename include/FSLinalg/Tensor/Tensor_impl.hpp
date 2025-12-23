#ifndef FSLINALG_TENSOR_IMPL_HPP
#define FSLINALG_TENSOR_IMPL_HPP

#include <FSLinalg/Tensor/Tensor.hpp>

#include <cassert>
#include <random>

namespace FSLinalg
{

template<typename T, unsigned int... dims> template<typename U, unsigned int d>
void Tensor<T,dims...>::initFromNestedInitializerList(misc::NestedInitializerList<U, d> values, Scalar* data)
{
	assert(values.size() == shape[d-1]);
	if constexpr (d == 1)
	{
		std::copy(std::cbegin(values), std::cend(values), data); 
	}
	else
	{
		for (const misc::NestedInitializerList<U, d-1>& inner_values : values)
		{
			initFromNestedInitializerList<U, d-1>(inner_values, data);
			data += strides[d-1];
		}
	}
}

template<typename T, unsigned int... dims>
auto Tensor<T,dims...>::random(const RealScalar& lb, const RealScalar& ub) -> Tensor
{
	std::random_device rd; 
	std::mt19937 gen(rd()); 
	std::uniform_real_distribution<RealScalar> dist(lb, ub);
	
	Tensor ret;
	
	for (Size i=0; i!=size; ++i) { ret[i] = dist(gen); }
	
	return ret;
}

} // namespace FSLinalg

#endif // FSLINALG_TENSOR_IMPL_HPP
