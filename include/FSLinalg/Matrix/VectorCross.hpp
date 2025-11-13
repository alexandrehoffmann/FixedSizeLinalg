#ifndef FSLINALG_VECTOR_CROSS_HPP
#define FSLINALG_VECTOR_CROSS_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> class VectorCross;

template<class Lhs, class Rhs>
struct MatrixTraits< VectorCross<Lhs, Rhs> >
{	
	static_assert(IsMatrix<Lhs>::value and IsMatrix<Rhs>::value, "Both LHS and RHS must be matrices");
	static_assert(Lhs::isRowVector and Lhs::size == 3, "Lhs must be a 3D RowVector");
	static_assert(Rhs::isRowVector and Rhs::size == 3, "Lhs must be a 3D RowVector");
	
	using Scalar = decltype(std::declval<typename Lhs::Scalar>() * std::declval<typename Rhs::Scalar>());
	using Size   = std::common_type_t<typename Lhs::Size, typename Rhs::Size>;

	static constexpr bool hasReadRandomAccess  = false;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool hasFlatRandomAccess  = false;
	static constexpr bool causesAliasingIssues = true;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size nRows = 3;   
	static constexpr Size nCols = 1;  
};

template<class Lhs, class Rhs>
class VectorCross : public MatrixBase< VectorCross<Lhs,Rhs> >
{
public:
	using Self = VectorCross<Lhs,Rhs>;
	FSLINALG_DEFINE_MATRIX
	
	VectorCross(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) : m_lhs(lhs.derived()), m_rhs(rhs.derived()) {}
	
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& dst) const { return m_lhs.isAliasedToImpl(dst) or m_rhs.isAliasedToImpl(dst); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
private:
	std::conditional_t<Lhs::isLeaf, const Lhs&, Lhs> m_lhs;
	std::conditional_t<Rhs::isLeaf, const Rhs&, Rhs> m_rhs;
	
	using TmpLhs = std::conditional_t<Lhs::hasReadRandomAccess, const Lhs&, RowVector<typename Lhs::Scalar, Lhs::size> >;
	using TmpRhs = std::conditional_t<Rhs::hasReadRandomAccess, const Rhs&, RowVector<typename Rhs::Scalar, Rhs::size> >;
	
	template<typename Alpha, class LhsVec, class RhsVec, class DstVec> static void assignCross    (const Alpha& alpha, const LhsVec& lhs, const RhsVec& rhs, DstVec& dst);
	template<typename Alpha, class LhsVec, class RhsVec, class DstVec> static void incrementCross (const Alpha& alpha, const LhsVec& lhs, const RhsVec& rhs, DstVec& dst);
	template<typename Alpha, class LhsVec, class RhsVec, class DstVec> static void decrementCross (const Alpha& alpha, const LhsVec& lhs, const RhsVec& rhs, DstVec& dst);
};

template<class Lhs, class Rhs> requires(Lhs::isRowVector and Lhs::size == 3 and Rhs::isRowVector and Rhs::size == 3)
VectorCross<Lhs, Rhs> cross(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) { return VectorCross<Lhs, Rhs>(lhs, rhs); }

} // namespace FSLinalg

#endif // FSLINALG_VECTOR_CROSS_HPP
