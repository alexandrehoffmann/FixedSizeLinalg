#ifndef FSLINALG_SCALAR_HPP
#define FSLINALG_SCALAR_HPP

#include <complex>
#include <type_traits>

#include <BIC/Core.hpp>

#include <FSLinalg/NumTraits.hpp>

namespace FSLinalg
{

template<typename T> concept        Scalar_concept = requires { typename NumTraits<T>::Real; };
template<typename T> concept    RealScalar_concept = Scalar_concept<T> and (not NumTraits<T>::isComplex);
template<typename T> concept ComplexScalar_concept = Scalar_concept<T> and      NumTraits<T>::isComplex;

template<typename T> struct IsScalar        : BIC::Fixed<bool, Scalar_concept<T> > {};
template<typename T> struct IsRealScalar    : BIC::Fixed<bool, RealScalar_concept<T> > {};
template<typename T> struct IsComplexScalar : BIC::Fixed<bool, ComplexScalar_concept<T> > {};

template<RealScalar_concept T> constexpr const T& real (const T& v) { return v;           }
template<RealScalar_concept T> constexpr       T  imag (const T&  ) { return 0;           }
template<RealScalar_concept T> constexpr const T& conj (const T& v) { return v;           }
template<RealScalar_concept T>           const T& abs  (const T& v) { return std::abs(v); }
template<RealScalar_concept T>           const T& abs2 (const T& v) { return v*v;         }

template<RealScalar_concept T> T& conjInPlace(T& v) { return v; }

template<RealScalar_concept T> constexpr const T&        real (const std::complex<T>& z) { return reinterpret_cast<const T(&)[2]>(z)[0]; }
template<RealScalar_concept T> constexpr const T&        imag (const std::complex<T>& z) { return reinterpret_cast<const T(&)[2]>(z)[1]; }
template<RealScalar_concept T>           std::complex<T> conj (const std::complex<T>& z) { return z.conj();                              }
template<RealScalar_concept T>           T               abs  (const std::complex<T>& z) { return std::abs(z);                           }
template<RealScalar_concept T>           T               abs2 (const std::complex<T>& z) { return real(z)*real(z) + imag(z)*imag(z);     }

template<RealScalar_concept T> T& conjInPlace(std::complex<T>& z) { reinterpret_cast<const T(&)[2]>(z)[1] = -reinterpret_cast<const T(&)[2]>(z)[1]; return z; }

} // namespace FSLinalg

#endif // FSLINALG_SCALAR_HPP
