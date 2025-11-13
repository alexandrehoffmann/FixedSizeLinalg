#ifndef FSLINALG_BASIC_LINALG_TRIPLE_PRODUCT_IMPL_HPP
#define FSLINALG_BASIC_LINALG_TRIPLE_PRODUCT_IMPL_HPP

#include <FSLinalg/BasicLinalg/TripleProduct.hpp>

namespace FSLinalg
{
namespace BasicLinalg
{

template<bool conjugateA, bool conjugateB, bool conjugateC> template<Scalar_concept A, Scalar_concept B, Scalar_concept C>
constexpr auto TripleProduct<conjugateA,conjugateB,conjugateC>::operator()(const A& a, const B& b, const C& c) const -> CpxReturnType<A,B,C> 
{
	if      constexpr (not conjugateA and not conjugateB and not conjugateC) { return a*b*c;                   }
	else if constexpr (not conjugateA and not conjugateB and     conjugateC) { return a*b*conj(c);             }
	else if constexpr (not conjugateA and     conjugateB and not conjugateC) { return a*conj(b)*c;             }
	else if constexpr (not conjugateA and     conjugateB and     conjugateC) { return a*conj(b)*conj(c);       }
	else if constexpr (    conjugateA and not conjugateB and not conjugateC) { return conj(a)*b*c;             }
	else if constexpr (    conjugateA and not conjugateB and     conjugateC) { return conj(a)*b*conj(c);       }
	else if constexpr (    conjugateA and     conjugateB and not conjugateC) { return conj(a)*conj(b)*c;       }
	else if constexpr (    conjugateA and     conjugateB and     conjugateC) { return conj(a)*conj(b)*conj(c); }
}

} //namespace BasicLinalg
} // namespace FSLinalg


#endif // FSLINALG_BASIC_LINALG_TRIPLE_PRODUCT_IMPL_HPP
