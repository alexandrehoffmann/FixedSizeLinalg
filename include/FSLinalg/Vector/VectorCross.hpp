#ifndef FSLINALG_VECTOR_CROSS_HPP
#define FSLINALG_VECTOR_CROSS_HPP

#include <FSLinalg/Vector/VectorBase.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> class VectorCross;

template<class Lhs, class Rhs>
struct VectorTraits< VectorCross<Lhs, Rhs> >
{	
	static_assert(IsVector<Lhs>::value and IsVector<Rhs>::value, "Both LHS and RHS must be vectors");
	
	using LhsTraits = VectorTraits<Lhs>;
	using RhsTraits = VectorTraits<Rhs>;
	
	static_assert(LhsTraits::size == 3, "Only defined for 3D vectors");
	static_assert(RhsTraits::size == 3, "Only defined for 3D vectors");
	
	using Scalar = decltype(std::declval<typename LhsTraits::Scalar>() * std::declval<typename RhsTraits::Scalar>());
	using Size   = std::common_type_t<typename LhsTraits::Size, typename RhsTraits::Size>;
	
	static constexpr bool hasReadRandomAccess  = false;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool causesAliasingIssues = true;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size size = 3;   
};

template<class Lhs, class Rhs>
class VectorCross : public VectorBase< VectorCross<Lhs, Rhs> >
{
public:
	using Self = VectorCross<Lhs, Rhs>;
	FSLINALG_DEFINE_VECTOR
	
	VectorCross(const VectorBase<Lhs>& lhs, const VectorBase<Rhs>& rhs) : m_lhs(lhs.derived()), m_rhs(rhs.derived()) {}
	
	template<class Dst> bool isAliasedToImpl(const VectorBase<Dst>& dst) const { return m_lhs.isAliasedToImpl(dst) or m_rhs.isAliasedToImpl(dst); }
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& dst) const { return m_lhs.isAliasedToImpl(dst) or m_rhs.isAliasedToImpl(dst); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
private:
	std::conditional_t<Lhs::isLeaf, const Lhs&, Lhs> m_lhs;
	std::conditional_t<Rhs::isLeaf, const Rhs&, Rhs> m_rhs;
	
	using TmpLhs = std::conditional_t<Lhs::hasReadRandomAccess, const Lhs&, Vector<typename Lhs::Scalar, Lhs::size> >;
	using TmpRhs = std::conditional_t<Rhs::hasReadRandomAccess, const Rhs&, Vector<typename Rhs::Scalar, Rhs::size> >;
	
	template<typename Alpha, class LhsVec, class RhsVec, class DstVec> static void assignCross    (const Alpha& alpha, const LhsVec& lhs, const RhsVec& rhs, DstVec& dst);
	template<typename Alpha, class LhsVec, class RhsVec, class DstVec> static void incrementCross (const Alpha& alpha, const LhsVec& lhs, const RhsVec& rhs, DstVec& dst);
	template<typename Alpha, class LhsVec, class RhsVec, class DstVec> static void decrementCross (const Alpha& alpha, const LhsVec& lhs, const RhsVec& rhs, DstVec& dst);
};

template<class Lhs, class Rhs>
VectorCross<Lhs, Rhs> cross(const VectorBase<Lhs>& lhs, const VectorBase<Rhs>& rhs) { return VectorCross<Lhs, Rhs>(lhs, rhs); }

} // namespace FSLinalg

#endif // FSLINALG_VECTOR_CROSS_HPP
