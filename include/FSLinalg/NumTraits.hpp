#ifndef FSLINALG_NUM_TRAITS_HPP
#define FSLINALG_NUM_TRAITS_HPP

#include <complex>
#include <concepts>

#include <cstdint>

namespace FSLinalg
{

template<typename T> struct NumTraits;

template<std::floating_point T>
struct NumTraits<T>
{	
	using Real = T;
	
	static constexpr bool isComplex = false;
	static constexpr Real epsilon   = std::numeric_limits<Real>::epsilon();
	static constexpr Real max       = std::numeric_limits<Real>::max();
	static constexpr Real min       = std::numeric_limits<Real>::min();
	static constexpr Real infinity  = std::numeric_limits<Real>::infinity();
};
	
template<std::floating_point T>
struct NumTraits< std::complex<T> >
{
	using Real = T;
	
	static constexpr bool isComplex = true;
	static constexpr Real epsilon   = std::numeric_limits<Real>::epsilon();
	static constexpr Real max       = std::numeric_limits<Real>::max();
	static constexpr Real min       = std::numeric_limits<Real>::min();
	static constexpr Real infinity  = std::numeric_limits<Real>::infinity();
};

//// for unit vectors

template<>
struct NumTraits<int>
{
	using Real = int;
	
	static constexpr bool isComplex = false;
	static constexpr Real epsilon   = std::numeric_limits<uint8_t>::epsilon();
	static constexpr Real max       = std::numeric_limits<uint8_t>::max();
	static constexpr Real min       = std::numeric_limits<uint8_t>::min();
	static constexpr Real infinity  = std::numeric_limits<uint8_t>::infinity();
};

} // namespace FSLinalg

#endif // FSLINALG_NUM_TRAITS_HPP
