#ifndef FSLINALG_BASIC_LINALG_PRODUCT_HPP
#define FSLINALG_BASIC_LINALG_PRODUCT_HPP

#include <FSLinalg/Scalar.hpp>

namespace FSLinalg
{
namespace BasicLinalg
{

template<bool conjugateA, bool conjugateB>
struct Product
{
	template<RealScalar_concept A, RealScalar_concept B>
	using RealReturnType = decltype(std::declval<A>() * std::declval<B>());
	
	template<Scalar_concept A, Scalar_concept B>
	using CpxReturnType = std::complex<RealReturnType<
		typename NumTraits<A>::Real,
		typename NumTraits<B>::Real>>;
	
	template<Scalar_concept A, Scalar_concept B>
	constexpr CpxReturnType<A,B> operator()(const A& a, const B& b) const;
	
	template<RealScalar_concept A, RealScalar_concept B>
	constexpr RealReturnType<A,B> operator()(const A& a, const B& b) const { return a*b; }
};

} //namespace BasicLinalg
} // namespace FSLinalg

#include <FSLinalg/BasicLinalg/Product_impl.hpp>

#endif // FSLINALG_BASIC_LINALG_PRODUCT_HPP
