#ifndef FSLINALG_MATRIX_SUM_HPP
#define FSLINALG_MATRIX_SUM_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> class MatrixSum;

template<class Lhs, class Rhs>
struct MatrixTraits< MatrixSum<Lhs,Rhs> >
{
	static_assert(IsMatrix<Lhs>::value and IsMatrix<Rhs>::value, "Both LHS and RHS must be matrices");
	
	using LhsTraits = MatrixTraits<Lhs>;
	using RhsTraits = MatrixTraits<Rhs>;
	
	static_assert(LhsTraits::nRows == RhsTraits::nCols, "Matrices size must match");
	static_assert(LhsTraits::nCols == RhsTraits::nCols, "Matrices size must match");
	
	using Scalar = decltype(std::declval<typename LhsTraits::Scalar>() + std::declval<typename RhsTraits::Scalar>());
	using Size   = std::common_type_t<typename LhsTraits::Size, typename RhsTraits::Size>;
	
	static constexpr bool hasReadRandomAccess  = LhsTraits::hasReadRandomAccess and RhsTraits::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool hasFlatRandomAccess  = LhsTraits::hasFlatRandomAccess and RhsTraits::hasFlatRandomAccess;
	static constexpr bool causesAliasingIssues = LhsTraits::causesAliasingIssues or RhsTraits::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size nRows = LhsTraits::nRows;   
	static constexpr Size nCols = LhsTraits::nCols;   
};

template<class Lhs, class Rhs>
class MatrixSum : public MatrixBase< MatrixSum<Lhs,Rhs> >
{
public:
	using Self = MatrixSum<Lhs,Rhs>;
	FSLINALG_DEFINE_MATRIX
	
	MatrixSum(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) : m_lhs(lhs.derived()), m_rhs(rhs.derived()) {}
	
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess) { return m_lhs.getImpl(i,j) + m_rhs.getImpl(i,j); }
	
	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess and hasFlatRandomAccess) { return m_lhs.getImpl(i) + m_rhs.getImpl(i); }
	
	template<class Dst> bool isAliasedToImpl(const VectorBase<Dst>& other) const { return m_lhs.isAliasedToImpl(other) or m_rhs.isAliasedToImpl(other); }
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& other) const { return m_lhs.isAliasedToImpl(other) or m_rhs.isAliasedToImpl(other); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)  { m_lhs.assignTo(alpha, dst, bc); m_rhs.increment(alpha, dst, bc); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_lhs.increment(alpha, dst, bc); m_rhs.increment(alpha, dst, bc); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_lhs.decrement(alpha, dst, bc); m_rhs.decrement(alpha, dst, bc); }
private:
	std::conditional_t<Lhs::isLeaf, const Lhs&, Lhs> m_lhs;
	std::conditional_t<Rhs::isLeaf, const Rhs&, Rhs> m_rhs;
};

template<class Lhs, class Rhs> 
FSLinalg::MatrixSum<Lhs,Rhs> operator+(const FSLinalg::MatrixBase<Lhs>& lhs, const FSLinalg::MatrixBase<Rhs>& rhs) { return FSLinalg::MatrixSum<Lhs,Rhs>(lhs, rhs); }

} // FSLinalg

#endif // FSLINALG_MATRIX_SUM_HPP
