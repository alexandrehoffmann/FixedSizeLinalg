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

template<typename T, T value>
struct NumTraits< BIC::Fixed<T, value> >
{
	using Real = typename NumTraits<T>::Real;
	
	static constexpr bool isComplex = NumTraits<T>::isComplex;
	static constexpr Real epsilon   = NumTraits<T>::epsilon;
	static constexpr Real max       = NumTraits<T>::max;
	static constexpr Real min       = NumTraits<T>::min;
	static constexpr Real infinity  = NumTraits<T>::infinity;
};

//// for unit vectors

template<>
struct NumTraits<int>
{
	using Real = int;
	
	static constexpr bool isComplex = false;
	static constexpr Real epsilon   = std::numeric_limits<int>::epsilon();
	static constexpr Real max       = std::numeric_limits<int>::max();
	static constexpr Real min       = std::numeric_limits<int>::min();
	static constexpr Real infinity  = std::numeric_limits<int>::infinity();
};

//// for tensor comparison operators

template<>
struct NumTraits<bool>
{
	using Real = bool;
	
	static constexpr bool isComplex = false;
	static constexpr Real epsilon   = std::numeric_limits<bool>::epsilon();
	static constexpr Real max       = std::numeric_limits<bool>::max();
	static constexpr Real min       = std::numeric_limits<bool>::min();
	static constexpr Real infinity  = std::numeric_limits<bool>::infinity();
};

} // namespace FSLinalg

#endif // FSLINALG_NUM_TRAITS_HPP
