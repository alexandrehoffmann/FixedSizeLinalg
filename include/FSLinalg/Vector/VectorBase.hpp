#ifndef FSLINALG_VECTOR_BASE_HPP
#define FSLINALG_VECTOR_BASE_HPP

#include <FSLinalg/CRTPBase.hpp>
#include <FSLinalg/Scalar.hpp>
#include <type_traits>

namespace FSLinalg
{
	
template<class Derived> class MatrixBase;

template<class Derived> struct VectorTraits;

template<class Derived>
class VectorBase : public CRTPBase<Derived>
{
public:
	using Base             = CRTPBase<Derived>;	
	using DerivedTraits    = VectorTraits<Derived>;
	using Scalar           = typename DerivedTraits::Scalar;
	using RealScalar       = typename NumTraits<Scalar>::Real;
	using ReturnType       = std::conditional_t<DerivedTraits::hasWriteRandomAccess,       Scalar&, Scalar>;
	using const_ReturnType = std::conditional_t<DerivedTraits::hasWriteRandomAccess, const Scalar&, Scalar>;
	using Size             = typename DerivedTraits::Size;
	
	template<class Dst> 
	struct IsConvertibleTo : std::bool_constant<
		std::is_base_of<VectorBase<Dst>, Dst>::value 
		and Dst::hasWriteRandomAccess
		and std::is_convertible<Scalar, typename Dst::Scalar>::value> {};
		
	template<class Src> 
	struct IsConstructibleFrom : std::bool_constant<
		std::is_base_of<VectorBase<Src>, Src>::value 
		and DerivedTraits::hasWriteRandomAccess
		and std::is_convertible<typename Src::Scalar, Scalar>::value> {};

	static constexpr bool hasReadRandomAccess  = DerivedTraits::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = DerivedTraits::hasWriteRandomAccess;
	static constexpr bool causesAliasingIssues = DerivedTraits::causesAliasingIssues;
	static constexpr bool isLeaf               = DerivedTraits::isLeaf;
	static constexpr Size size                 = DerivedTraits::size;
	
	constexpr Size getSize() const { return size; }
	
	const_ReturnType operator[](const Size& i) const requires(hasReadRandomAccess)  { return Base::derived().getImpl(i); }
	      ReturnType operator[](const Size& i)       requires(hasWriteRandomAccess) { return Base::derived().getImpl(i); }

	template<class Dst>           bool isAliasedTo(const VectorBase<Dst>& dst) const { return Base::derived().isAliasedToImpl(dst); }
	template<class Dst> constexpr bool isAliasedTo(const MatrixBase<Dst>&    ) const { return false; }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignTo(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void increment(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrement(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
};

template<typename Expr> struct IsVector : std::bool_constant< std::is_base_of<VectorBase<Expr>, Expr>::value > {};

template<typename Expr> concept Vector_concept         = IsVector<Expr>::value;
template<typename Expr> concept ReadableVector_concept = IsVector<Expr>::value and Expr::hasReadRandomAccess;
template<typename Expr> concept WritableVector_concept = IsVector<Expr>::value and Expr::hasWriteRandomAccess;

template<typename Lhs, typename Rhs> requires(Lhs::hasReadRandomAccess and Rhs::hasReadRandomAccess) bool operator==(const FSLinalg::VectorBase<Lhs>& lhs, const FSLinalg::VectorBase<Rhs>& rhs);
template<typename Lhs, typename Rhs> requires(Lhs::hasReadRandomAccess and Rhs::hasReadRandomAccess) bool operator!=(const FSLinalg::VectorBase<Lhs>& lhs, const FSLinalg::VectorBase<Rhs>& rhs) { return not(lhs == rhs); }

} // namespace FSLinalg

#define FSLINALG_DEFINE_VECTOR \
	using Base             = VectorBase<Self>; \
	using Scalar           = typename Base::Scalar; \
	using RealScalar       = typename Base::RealScalar; \
	using ReturnType       = typename Base::ReturnType; \
	using const_ReturnType = typename Base::const_ReturnType; \
	using Size             = typename Base::Size; \
	\
	template<class Dst> using IsConvertibleTo     = typename Base::template IsConvertibleTo<Dst>; \
	template<class Src> using IsConstructibleFrom = typename Base::template IsConstructibleFrom<Src>; \
	\
	using Base::hasReadRandomAccess;\
	using Base::hasWriteRandomAccess;\
	using Base::causesAliasingIssues;\
	using Base::isLeaf;\
	using Base::size;\
	\

#endif // FSLINALG_VECTOR_BASE_HPP
