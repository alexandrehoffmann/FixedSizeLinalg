#ifndef FSLINALG_VECTOR_HPP
#define FSLINALG_VECTOR_HPP

#include <FSLinalg/Vector/VectorBase.hpp>

#include <cassert>
#include <array>

namespace FSLinalg
{

template<typename T, unsigned int N> class Vector;

template<typename T, unsigned int N>
struct VectorTraits< Vector<T, N> >
{	
	using Scalar     = T;
	using Size       = unsigned int;
	
	static constexpr bool hasReadRandomAccess  = true;
	static constexpr bool hasWriteRandomAccess = true;
	static constexpr bool causesAliasingIssues = false;
	static constexpr bool isLeaf               = true;
	
	static constexpr Size size = N;   
};

template<typename T, unsigned int N>
class Vector : public VectorBase< Vector<T, N> >
{
public:
	using Self = Vector<T, N>;
	FSLINALG_DEFINE_VECTOR
	
	static constexpr bool IsScalarComplex = IsComplexScalar<Scalar>::value;
	
	Vector(const RealScalar& value = RealScalar(0))  requires(IsScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] = value; } }
	Vector(std::initializer_list<RealScalar> values) requires(IsScalarComplex) { assert(values.size() == size); std::copy(std::cbegin(values), std::cend(values), std::begin(m_data)); }
	
	Vector(const Scalar& value = Scalar(0))      { m_data.fill(value); }
	Vector(std::initializer_list<Scalar> values) { assert(values.size() == size); std::copy(std::cbegin(values), std::cend(values), std::begin(m_data)); }
	
	Vector(const Vector& other) : m_data(other.m_data) {}
	
	template<class Expr> Vector(const VectorBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.assignTo(1., *this, std::false_type{}); }
	
	template<class Expr> Vector& operator= (const VectorBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.assignTo  (1., *this, std::true_type{}); return *this; }
	template<class Expr> Vector& operator+=(const VectorBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.increment (1., *this, std::true_type{}); return *this; }
	template<class Expr> Vector& operator-=(const VectorBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.decrement (1., *this, std::true_type{}); return *this; }
	
	Vector& operator*=(const RealScalar& alpha) requires(IsScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] *= alpha; } return *this; }
	Vector& operator/=(const RealScalar& alpha) requires(IsScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] /= alpha; } return *this; }
	
	Vector& operator*=(const Scalar& alpha) { for (Size i=0; i!=size; ++i) { m_data[i] *= alpha; } return *this; }
	Vector& operator/=(const Scalar& alpha) { for (Size i=0; i!=size; ++i) { m_data[i] /= alpha; } return *this; }
	
	const_ReturnType getImpl(const Size i) const { return m_data[i]; }
	      ReturnType getImpl(const Size i)       { return m_data[i]; }
	      
	template<class Dst> 
	bool isAliasedToImpl(const VectorBase<Dst>& dst) const 
	{ 
		if constexpr (size == Dst::size) { return std::addressof(dst.derived()) == this; }
		else                             { return false;                                 }
	}
	
	template<class Dst> 
	constexpr bool isAliasedToImpl(const MatrixBase<Dst>&) const { return false; }
private:
	std::array<Scalar, size> m_data;
};

template<unsigned int N> using RealVector = Vector<double, N>;
template<unsigned int N> using CpxVector  = Vector<std::complex<double>, N>;

} // namespace FSLinalg

#endif // FSLINALG_VECTOR_HPP
