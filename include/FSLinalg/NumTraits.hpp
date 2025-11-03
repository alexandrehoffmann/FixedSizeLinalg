#ifndef FSLINALG_NUM_TRAITS_HPP
#define FSLINALG_NUM_TRAITS_HPP

#include <complex>
#include <concepts>

namespace FSLinalg
{

template<typename T> struct NumTraits;

template<std::floating_point T>
struct NumTraits<T>
{	
	using Real = T;
	
	static constexpr bool IsComplex = false;
	static constexpr Real epsilon   = std::numeric_limits<Real>::epsilon();
	static constexpr Real max       = std::numeric_limits<Real>::max();
	static constexpr Real min       = std::numeric_limits<Real>::min();
	static constexpr Real infinity  = std::numeric_limits<Real>::infinity();
};
	
template<std::floating_point T>
struct NumTraits< std::complex<T> >
{
	using Real = T;
	
	static constexpr bool IsComplex = true;
	static constexpr Real epsilon   = std::numeric_limits<Real>::epsilon();
	static constexpr Real max       = std::numeric_limits<Real>::max();
	static constexpr Real min       = std::numeric_limits<Real>::min();
	static constexpr Real infinity  = std::numeric_limits<Real>::infinity();
};

//// for unit vectors

template<>
struct NumTraits<bool>
{
	using Real = bool;
	
	static constexpr bool IsComplex = false;
	static constexpr Real epsilon   = std::numeric_limits<bool>::epsilon();
	static constexpr Real max       = std::numeric_limits<bool>::max();
	static constexpr Real min       = std::numeric_limits<bool>::min();
	static constexpr Real infinity  = std::numeric_limits<bool>::infinity();
};

} // namespace FSLinalg

#endif // FSLINALG_NUM_TRAITS_HPP
