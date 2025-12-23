#ifndef FSLINALG_TENSOR_HPP
#define FSLINALG_TENSOR_HPP

#include <FSLinalg/Tensor/TensorBase.hpp>
#include <FSLinalg/Tensor/TensorUtils.hpp>
#include <FSLinalg/misc/NestedInitializerList.hpp>

#include <array>
#include <numeric>

namespace FSLinalg
{

template<typename T, unsigned int... dims> class Tensor;

template<typename T, unsigned int... dims>
struct TensorTraits< Tensor<T, dims...> >
{		
	static_assert(sizeof...(dims) > 0);
	
	using Scalar = T;
	using Size   = unsigned int;
	using Shape  = std::array<Size, sizeof...(dims)>;
	
	static constexpr bool hasReadRandomAccess  = true;
	static constexpr bool hasWriteRandomAccess = true;
	static constexpr bool hasFlatRandomAccess  = true;
	static constexpr bool causesAliasingIssues = false;
	static constexpr bool isLeaf               = true;
	
	static constexpr Shape shape = Shape({dims...});
};

template<typename T, unsigned int... dims> 
class Tensor : public TensorBase< Tensor<T, dims...> >
{
public:
	using Self = Tensor<T, dims...>;
	FSLINALG_DEFINE_TENSOR
	
	template<class Dst>
	struct CanBeAlisaedTo : BIC::Fixed<bool,  
		    IsTensor<Dst>::value 
		and Base::shape == Dst::shape
		and std::is_same<Scalar, typename Dst::Scalar>::value > {};
	
	static constexpr bool isScalarComplex = IsComplexScalar<Scalar>::value;
	
	static constexpr Shape strides = TensorUtils::getStrides(shape);
	
	Tensor(const RealScalar& value = RealScalar(0))              requires(isScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] = value; } }
	Tensor(misc::NestedInitializerList<RealScalar, rank> values) requires(isScalarComplex) { initFromNestedInitializerList<RealScalar, rank>(values, m_data.data()); }
	
	Tensor(const Scalar& value = Scalar(0)) { m_data.fill(value); }
	Tensor(misc::NestedInitializerList<Scalar, rank> values) { initFromNestedInitializerList<Scalar, rank>(values, m_data.data()); }
	
	Tensor(const Tensor& other) : m_data(other.m_data) {}
	
	template<class Expr> Tensor(const TensorBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.assignTo(BIC::fixed<bool, false>, BIC::fixed<RealScalar, RealScalar(1)>, *this); }
	
	template<class Expr> Tensor& operator= (const TensorBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.assignTo  (BIC::fixed<bool, true>, BIC::fixed<RealScalar, RealScalar(1)>, *this); return *this; }
	template<class Expr> Tensor& operator+=(const TensorBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.increment (BIC::fixed<bool, true>, BIC::fixed<RealScalar, RealScalar(1)>, *this); return *this; }
	template<class Expr> Tensor& operator-=(const TensorBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.decrement (BIC::fixed<bool, true>, BIC::fixed<RealScalar, RealScalar(1)>, *this); return *this; }
	template<class Expr> Tensor& operator*=(const TensorBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.multiply  (BIC::fixed<bool, true>, BIC::fixed<RealScalar, RealScalar(1)>, *this); return *this; }
	template<class Expr> Tensor& operator/=(const TensorBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.divide    (BIC::fixed<bool, true>, BIC::fixed<RealScalar, RealScalar(1)>, *this); return *this; }
	
	Tensor& operator*=(const RealScalar& alpha) requires(isScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] *= alpha; } return *this; }
	Tensor& operator/=(const RealScalar& alpha) requires(isScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] /= alpha; } return *this; }
	
	Tensor& operator*=(const Scalar& alpha) { for (Size i=0; i!=size; ++i) { m_data[i] *= alpha; } return *this; }
	Tensor& operator/=(const Scalar& alpha) { for (Size i=0; i!=size; ++i) { m_data[i] /= alpha; } return *this; }
	
	void setZero() { m_data.fill(Scalar(0)); }
	
	const_ReturnType getImpl(const Size i) const { return m_data[i]; }
	      ReturnType getImpl(const Size i)       { return m_data[i]; }
	      
	template<std::integral... Idx> const_ReturnType getImpl(const Idx... idx) const requires(sizeof...(Idx) == rank) { return m_data[toFlatIndex(idx...)]; }
	template<std::integral... Idx>       ReturnType getImpl(const Idx... idx)       requires(sizeof...(Idx) == rank) { return m_data[toFlatIndex(idx...)]; }
	      
	template<class Dst>           bool isAliasedToImpl(const TensorBase<Dst>& dst) const requires(    CanBeAlisaedTo<Dst>::value) { return std::addressof(dst.derived()) == this; }
	template<class Dst> constexpr bool isAliasedToImpl(const TensorBase<Dst>&    ) const requires(not CanBeAlisaedTo<Dst>::value) { return false; }
	
	static Tensor zero() { return Tensor(RealScalar(0)); }
	static Tensor ones() { return Tensor(RealScalar(1)); }
	
	static Tensor random(const RealScalar& lb = RealScalar(-1), const RealScalar& ub = RealScalar(1));
private:
	template<std::integral... Idx, size_t... Is> Size toFlatIndexHelper(BIC::FixedIndices<Is...>, const Idx... idx) const requires(sizeof...(Idx) == rank and sizeof...(Is) == rank) { return ((idx*strides[Is]) + ...); } 
	
	template<std::integral... Idx> Size toFlatIndex(const Idx... idx) const requires(sizeof...(Idx) == rank) { return toFlatIndexHelper(BIC::indexSeq<0, rank>, idx...); } 
	
	template<typename U, unsigned int d> static void initFromNestedInitializerList(misc::NestedInitializerList<U, d> values, Scalar* data);

	std::array<Scalar, size> m_data;
};

template<unsigned int... dims> using BoolTensor = Tensor<bool, dims...>;
template<unsigned int... dims> using RealTensor = Tensor<double, dims...>;
template<unsigned int... dims> using CpxTensor  = Tensor<std::complex<double>, dims...>;

namespace detail
{

template <typename T, std::array shape, size_t... Is>  constexpr Tensor<T, shape[Is]...> tensorLike(BIC::FixedIndices<Is...>) { return {}; }

} // namespace detail

template<typename T, std::array shape> using TensorFromShape = decltype(detail::tensorLike<T,shape>(BIC::indexSeq<0, shape.size()>));

} // namespace FSLinalg

#endif // FSLINALG_TENSOR_HPP
