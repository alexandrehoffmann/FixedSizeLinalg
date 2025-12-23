#ifndef FSLINALG_TENSOR_UTILS_HPP
#define FSLINALG_TENSOR_UTILS_HPP

#include <array>
#include <cstddef>

namespace FSLinalg
{
namespace TensorUtils
{

template<typename Size, size_t rank> constexpr std::array<Size, rank> getStrides(const std::array<Size, rank>& idx);
	
} // namespace TensorUtils
} // namespace FSLinalg

#include <FSLinalg/Tensor/TensorUtils_impl.hpp>

#endif // FSLINALG_TENSOR_UTILS_HPP
