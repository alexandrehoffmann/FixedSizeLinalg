#ifndef FSLINALG_VECTOR_SUB_HPP
#define FSLINALG_VECTOR_SUB_HPP

#include <FSLinalg/Vector/VectorBase.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> class VectorSub;

template<class Lhs, class Rhs>
struct VectorTraits< VectorSub<Lhs,Rhs> >
{
	static_assert(IsVector<Lhs>::value and IsVector<Rhs>::value, "Both LHS and RHS must be vectors");
	
	using LhsTraits = VectorTraits<Lhs>;
	using RhsTraits = VectorTraits<Rhs>;
	
	static_assert(LhsTraits::size == RhsTraits::size, "Vectors size must match");
	
	using Scalar = decltype(std::declval<typename LhsTraits::Scalar>() - std::declval<typename RhsTraits::Scalar>());
	using Size   = std::common_type_t<typename LhsTraits::Size, typename RhsTraits::Size>;
	
	static constexpr bool hasReadRandomAccess  = LhsTraits::hasReadRandomAccess and RhsTraits::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool causesAliasingIssues = LhsTraits::causesAliasingIssues or RhsTraits::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size size = LhsTraits::size;   
};

template<class Lhs, class Rhs> 
class VectorSub : public VectorBase< VectorSub<Lhs, Rhs> >
{
public:
	using Self = VectorSub<Lhs, Rhs>;
	FSLINALG_DEFINE_VECTOR
	
	VectorSub(const VectorBase<Lhs>& lhs, const VectorBase<Rhs>& rhs) : m_lhs(lhs.derived()), m_rhs(rhs.derived()) {}
	
	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess) { return m_lhs.getImpl(i) - m_rhs.getImpl(i); }
	
	template<class Dst> bool isAliasedToImpl(const VectorBase<Dst>& dst) const { return m_lhs.isAliasedToImpl(dst) or m_rhs.isAliasedToImpl(dst); }
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& dst) const { return m_lhs.isAliasedToImpl(dst) or m_rhs.isAliasedToImpl(dst); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)  { m_lhs.assignTo(alpha, dst, bc); m_rhs.decrement(alpha, dst, bc); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_lhs.increment(alpha, dst, bc); m_rhs.decrement(alpha, dst, bc); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_lhs.decrement(alpha, dst, bc); m_rhs.increment(alpha, dst, bc); }
private:
	std::conditional_t<Lhs::isLeaf, const Lhs&, Lhs> m_lhs;
	std::conditional_t<Rhs::isLeaf, const Rhs&, Rhs> m_rhs;
};

template<class Lhs, class Rhs> 
FSLinalg::VectorSub<Lhs,Rhs> operator-(const FSLinalg::VectorBase<Lhs>& lhs, const FSLinalg::VectorBase<Rhs>& rhs) { return FSLinalg::VectorSub<Lhs,Rhs>(lhs, rhs); }

} // namespace FSLinalg

#endif // FSLINALG_VECTOR_SUB_HPP
