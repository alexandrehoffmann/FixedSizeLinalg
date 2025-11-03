#ifndef FSLINALG_VECTOR_OUTER_PRODUCT_HPP
#define FSLINALG_VECTOR_OUTER_PRODUCT_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>
#include <FSLinalg/Vector/VectorBase.hpp>
#include <FSLinalg/Vector/Vector.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> class VectorOuterProduct;

template<class Lhs, class Rhs>
struct MatrixTraits< VectorOuterProduct<Lhs,Rhs> >
{
	static_assert(IsVector<Lhs>::value and IsVector<Rhs>::value, "Both LHS and RHS must be vectors");
	
	using LhsTraits = VectorTraits<Lhs>;
	using RhsTraits = VectorTraits<Rhs>;
	
	using Scalar = decltype(std::declval<typename LhsTraits::Scalar>() * std::declval<typename RhsTraits::Scalar>());
	using Size   = std::common_type_t<typename LhsTraits::Size, typename RhsTraits::Size>;
	
	static constexpr bool hasReadRandomAccess  = LhsTraits::hasReadRandomAccess and RhsTraits::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool hasFlatRandomAccess  = false;
	static constexpr bool causesAliasingIssues = LhsTraits::causesAliasingIssues or RhsTraits::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size nRows = LhsTraits::size;   
	static constexpr Size nCols = RhsTraits::size;   
};

template<class Lhs, class Rhs>
class VectorOuterProduct : public MatrixBase< VectorOuterProduct<Lhs,Rhs> >
{
public:
	using Self = VectorOuterProduct<Lhs,Rhs>;
	FSLINALG_DEFINE_MATRIX
	
	VectorOuterProduct(const VectorBase<Lhs>& lhs, const VectorBase<Rhs>& rhs) : m_lhs(lhs.derived()), m_rhs(rhs.derived()) {}
	
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess) { return m_lhs.getImpl(i)*m_rhs.getImpl(j); }
	
	template<class Dst> bool isAliasedToImpl(const VectorBase<Dst>& other) const { return m_lhs.isAliasedToImpl(other) or m_rhs.isAliasedToImpl(other); }
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& other) const { return m_lhs.isAliasedToImpl(other) or m_rhs.isAliasedToImpl(other); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
private:
	using TmpLhs = std::conditional_t<Lhs::hasReadRandomAccess, const Lhs&, Vector<typename Lhs::Scalar, Lhs::size> >;
	using TmpRhs = std::conditional_t<Rhs::hasReadRandomAccess, const Rhs&, Vector<typename Rhs::Scalar, Rhs::size> >;

	std::conditional_t<Lhs::isLeaf, const Lhs&, Lhs> m_lhs;
	std::conditional_t<Rhs::isLeaf, const Rhs&, Rhs> m_rhs;
};

template<class Lhs, class Rhs>
VectorOuterProduct<Lhs,Rhs> outer(const VectorBase<Lhs>& lhs, const VectorBase<Rhs>& rhs) { return VectorOuterProduct<Lhs,Rhs>(lhs, rhs); }

} // FSLinalg

#endif // FSLINALG_VECTOR_OUTER_PRODUCT_HPP
