#ifndef FSLINALG_TENSOR_UTILS_IMPL_HPP
#define FSLINALG_TENSOR_UTILS_IMPL_HPP

namespace FSLinalg
{
namespace TensorUtils
{

template<typename Size, size_t rank> 
constexpr std::array<Size, rank> getStrides(const std::array<Size, rank>& shape)
{
	std::array<Size, rank> strides;
	
	for (Size i=0; i!=rank; ++i)
	{
		strides[i] = Size(1);
		for (Size j=i+1; j!=rank; ++j)
		{
			strides[i] *= shape[j];
		}
	}
	
	return strides;
}
	
} // namespace TensorUtils
} // namespace FSLinalg

#endif // FSLINALG_TENSOR_UTILS_IMPL_HPP
