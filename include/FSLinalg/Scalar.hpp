#ifndef FSLINALG_SCALAR_HPP
#define FSLINALG_SCALAR_HPP

#include <complex>
#include <type_traits>

#include <FSLinalg/NumTraits.hpp>

namespace FSLinalg
{

template<typename T> concept        Scalar_concept = requires { typename NumTraits<T>::Real; };
template<typename T> concept    RealScalar_concept = Scalar_concept<T> and (not NumTraits<T>::IsComplex);
template<typename T> concept ComplexScalar_concept = Scalar_concept<T> and NumTraits<T>::IsComplex;

template<typename T> struct IsScalar        : std::bool_constant< Scalar_concept<T> > {};
template<typename T> struct IsRealScalar    : std::bool_constant< RealScalar_concept<T> > {};
template<typename T> struct IsComplexScalar : std::bool_constant< ComplexScalar_concept<T> > {};

template<RealScalar_concept T> const T&    real(const T& v) { return v; }
template<RealScalar_concept T> constexpr T imag(const T&  ) { return 0; }
template<RealScalar_concept T> const T&    conj(const T& v) { return v; }

template<RealScalar_concept T> T& conjInPlace(T& v) { return v; }

template<RealScalar_concept T> const T&        real(const std::complex<T>& v) { return std::get<0>(v); }
template<RealScalar_concept T> const T&        imag(const std::complex<T>& v) { return std::get<1>(v); }
template<RealScalar_concept T> std::complex<T> conj(const std::complex<T>& v) { return v.conj();       }

template<RealScalar_concept T> T& conjInPlace(std::complex<T>& v) { std::get<1>(v) = -std::get<1>(v); return v; }

} // namespace FSLinalg

#endif // FSLINALG_SCALAR_HPP
