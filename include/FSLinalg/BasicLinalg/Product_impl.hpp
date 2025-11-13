#ifndef FSLINALG_BASIC_LINALG_PRODUCT_IMPL_HPP
#define FSLINALG_BASIC_LINALG_PRODUCT_IMPL_HPP

#include <FSLinalg/BasicLinalg/Product.hpp>

namespace FSLinalg
{
namespace BasicLinalg
{

template<bool conjugateA, bool conjugateB> template<Scalar_concept A, Scalar_concept B>
constexpr auto Product<conjugateA,conjugateB>::operator()(const A& a, const B& b) const -> CpxReturnType<A,B> 
{
	if      constexpr (not conjugateA and not conjugateB) { return a*b;             }
	else if constexpr (not conjugateA and     conjugateB) { return a*conj(b);       }
	else if constexpr (    conjugateA and not conjugateB) { return conj(a)*b;       }
	else if constexpr (    conjugateA and     conjugateB) { return conj(a)*conj(b); }
}

} //namespace BasicLinalg
} // namespace FSLinalg


#endif // FSLINALG_BASIC_LINALG_PRODUCT_IMPL_HPP
