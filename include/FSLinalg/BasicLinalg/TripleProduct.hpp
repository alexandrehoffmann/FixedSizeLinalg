#ifndef FSLINALG_BASIC_LINALG_TRIPLE_PRODUCT_HPP
#define FSLINALG_BASIC_LINALG_TRIPLE_PRODUCT_HPP

#include <FSLinalg/Scalar.hpp>

namespace FSLinalg
{
namespace BasicLinalg
{

template<bool conjugateA, bool conjugateB, bool conjugateC>
struct TripleProduct
{
	template<RealScalar_concept A, RealScalar_concept B, RealScalar_concept C>
	using RealReturnType = std::common_type_t<A, B, C>;
	
	template<Scalar_concept A, Scalar_concept B, Scalar_concept C>
	using CpxReturnType = std::complex<RealReturnType<
		typename NumTraits<A>::Real,
		typename NumTraits<B>::Real,
		typename NumTraits<C>::Real>>;
	
	template<Scalar_concept A, Scalar_concept B, Scalar_concept C>
	constexpr CpxReturnType<A,B,C> operator()(const A& a, const B& b, const C& c) const;
	
	template<RealScalar_concept A, RealScalar_concept B, RealScalar_concept C>
	constexpr RealReturnType<A,B,C> operator()(const A& a, const B& b, const C& c) const { return a*b*c; }
};

} //namespace BasicLinalg
} // namespace FSLinalg

#include <FSLinalg/BasicLinalg/TripleProduct_impl.hpp>

#endif // FSLINALG_BASIC_LINALG_TRIPLE_PRODUCT_HPP
