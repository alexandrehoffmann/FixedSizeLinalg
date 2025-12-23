#ifndef FSLINALG_TENSOR_BASE_HPP
#define FSLINALG_TENSOR_BASE_HPP

#include <FSLinalg/CRTPBase.hpp>
#include <FSLinalg/Scalar.hpp>
#include <FSLinalg/misc/Logical.hpp>

#include <numeric>
#include <concepts>

namespace FSLinalg
{

template<class Derived> struct TensorTraits;

template<class Derived>
class TensorBase : public CRTPBase<Derived>
{
public:
	using CRTP             = CRTPBase<Derived>;	
	using DerivedTraits    = TensorTraits<Derived>;
	using Scalar           = typename DerivedTraits::Scalar;
	using RealScalar       = typename NumTraits<Scalar>::Real;
	using ReturnType       = std::conditional_t<DerivedTraits::hasWriteRandomAccess,       Scalar&, Scalar>;
	using const_ReturnType = std::conditional_t<DerivedTraits::hasWriteRandomAccess, const Scalar&, Scalar>;
	using Size             = typename DerivedTraits::Size;
	using Shape            = typename DerivedTraits::Shape;
	
	template<class Dst> 
	struct IsConvertibleTo : BIC::Fixed<bool, 
		    std::is_base_of<TensorBase<Dst>, Dst>::value 
		and Dst::hasWriteRandomAccess
		and std::is_convertible<Scalar, typename Dst::Scalar>::value> {};
		
	template<class Src> 
	struct IsConstructibleFrom : BIC::Fixed<bool, 
		    std::is_base_of<TensorBase<Src>, Src>::value 
		and DerivedTraits::hasWriteRandomAccess
		and std::is_convertible<typename Src::Scalar, Scalar>::value> {};
		
	static constexpr bool  hasReadRandomAccess  = DerivedTraits::hasReadRandomAccess;
	static constexpr bool  hasWriteRandomAccess = DerivedTraits::hasWriteRandomAccess;
	static constexpr bool  hasFlatRandomAccess  = DerivedTraits::hasFlatRandomAccess;
	static constexpr bool  causesAliasingIssues = DerivedTraits::causesAliasingIssues;
	static constexpr bool  isLeaf               = DerivedTraits::isLeaf;
	static constexpr Shape shape                = DerivedTraits::shape;
	static constexpr Size  rank                 = shape.size();
	static constexpr Size  size                 = std::reduce(std::begin(shape), std::end(shape), Size(1), std::multiplies{});
	
	constexpr Size getRank  ()             const { return rank;     }
	constexpr Size getShape (const Size d) const { return shape[d]; }
	constexpr Size getSize  ()             const { return size;     }
	
	const_ReturnType operator()(const Shape& idx) const requires(hasReadRandomAccess)  { return getHelper(idx, BIC::indexSeq<0, rank>); }
	      ReturnType operator()(const Shape& idx)       requires(hasWriteRandomAccess) { return getHelper(idx, BIC::indexSeq<0, rank>); }
	
	template<std::integral... Idx> const_ReturnType operator()(const Idx... idx) const requires(sizeof...(Idx) == rank and hasReadRandomAccess)  { return CRTP::derived().getImpl(idx...); }
	template<std::integral... Idx>       ReturnType operator()(const Idx... idx)       requires(sizeof...(Idx) == rank and hasWriteRandomAccess) { return CRTP::derived().getImpl(idx...); }
	      
	const_ReturnType operator[](const Size i) const requires(hasReadRandomAccess  and hasFlatRandomAccess) { return CRTP::derived().getImpl(i); }
	      ReturnType operator[](const Size i)       requires(hasWriteRandomAccess and hasFlatRandomAccess) { return CRTP::derived().getImpl(i); }

	template<class Dst> bool isAliasedTo(const TensorBase<Dst>& other) const { return CRTP::derived().isAliasedToImpl(other); }

	template<typename Bool, typename Alpha, class Dst>
	void assignTo(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Bool, typename Alpha, class Dst>
	void increment(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Bool, typename Alpha, class Dst>
	void decrement(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Bool, typename Alpha, class Dst>
	void multiply(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Bool, typename Alpha, class Dst>
	void divide(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
protected:
	template<size_t... Is> const_ReturnType getHelper(const Shape& idx, BIC::FixedIndices<Is...>) const requires(sizeof...(Is) == rank and hasReadRandomAccess)  { return CRTP::derived().getImpl(idx[Is]...); }
	template<size_t... Is>       ReturnType getHelper(const Shape& idx, BIC::FixedIndices<Is...>)       requires(sizeof...(Is) == rank and hasWriteRandomAccess) { return CRTP::derived().getImpl(idx[Is]...); }
};

template<typename Expr> struct IsTensor : BIC::Fixed<bool,  std::is_base_of<TensorBase<Expr>, Expr>::value > {};

template<typename Expr> concept Tensor_concept         = IsTensor<Expr>::value;
template<typename Expr> concept ReadableTensor_concept = IsTensor<Expr>::value and Expr::hasReadRandomAccess;
template<typename Expr> concept WritableTensor_concept = IsTensor<Expr>::value and Expr::hasWriteRandomAccess;

#define FSLINALG_DEFINE_TENSOR \
	using Base             = TensorBase<Self>; \
	using Scalar           = typename Base::Scalar; \
	using RealScalar       = typename Base::RealScalar; \
	using ReturnType       = typename Base::ReturnType; \
	using const_ReturnType = typename Base::const_ReturnType; \
	using Size             = typename Base::Size; \
	using Shape            = typename Base::Shape; \
	\
	template<class Dst> using IsConvertibleTo     = typename Base::template IsConvertibleTo<Dst>; \
	template<class Src> using IsConstructibleFrom = typename Base::template IsConstructibleFrom<Src>; \
	\
	using Base::hasReadRandomAccess;\
	using Base::hasWriteRandomAccess;\
	using Base::hasFlatRandomAccess;\
	using Base::causesAliasingIssues;\
	using Base::isLeaf;\
	using Base::rank;\
	using Base::shape;\
	using Base::size;\
    \

template<typename Lhs, typename Rhs> requires(Lhs::hasReadRandomAccess and Rhs::hasReadRandomAccess)  bool operator==(const FSLinalg::TensorBase<Lhs>& lhs, const FSLinalg::TensorBase<Rhs>& rhs);

template<class Expr> requires(std::is_same<typename Expr::Scalar, bool>::value and Expr::hasReadRandomAccess) bool all(const TensorBase<Expr>& expr);
template<class Expr> requires(std::is_same<typename Expr::Scalar, bool>::value and Expr::hasReadRandomAccess) bool any(const TensorBase<Expr>& expr);
   
} // namespace FSLinalg

#include <FSLinalg/Tensor/Formater.hpp>

#endif // FSLINALG_TENSOR_BASE_HPP
