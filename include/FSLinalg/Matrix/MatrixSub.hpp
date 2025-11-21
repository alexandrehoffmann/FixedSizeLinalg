#ifndef FSLINALG_MATRIX_SUB_HPP
#define FSLINALG_MATRIX_SUB_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> class MatrixSub;

template<class Lhs, class Rhs>
struct MatrixTraits< MatrixSub<Lhs,Rhs> >
{
	static_assert(IsMatrix<Lhs>::value and IsMatrix<Rhs>::value, "Both LHS and RHS must be matrices");	
	static_assert(Lhs::nRows == Rhs::nRows, "Matrices size must match");
	static_assert(Lhs::nCols == Rhs::nCols, "Matrices size must match");
	
	using Scalar = decltype(std::declval<typename Lhs::Scalar>() - std::declval<typename Rhs::Scalar>());
	using Size   = std::common_type_t<typename Lhs::Size, typename Rhs::Size>;
	
	static constexpr bool hasReadRandomAccess  = Lhs::hasReadRandomAccess and Rhs::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool hasFlatRandomAccess  = Lhs::hasFlatRandomAccess and Rhs::hasFlatRandomAccess;
	static constexpr bool causesAliasingIssues = Lhs::causesAliasingIssues or Rhs::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size nRows = Lhs::nRows;   
	static constexpr Size nCols = Lhs::nCols;   
};

template<class Lhs, class Rhs>
class MatrixSub : public MatrixBase< MatrixSub<Lhs,Rhs> >
{
public:
	using Self = MatrixSub<Lhs,Rhs>;
	FSLINALG_DEFINE_MATRIX
	
	MatrixSub(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) : m_lhs(lhs.derived()), m_rhs(rhs.derived()) {}
	
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess) { return m_lhs.getImpl(i,j) - m_rhs.getImpl(i,j); }
	
	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess and hasFlatRandomAccess) { return m_lhs.getImpl(i) - m_rhs.getImpl(i); }
	
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& other) const { return m_lhs.isAliasedToImpl(other) or m_rhs.isAliasedToImpl(other); }

	template<typename Bool, typename Alpha, class Dst>
	void assignToImpl(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)  { m_lhs.assignTo(checkAliasing, alpha, dst); m_rhs.decrement(checkAliasing, alpha, dst); }
	
	template<typename Bool, typename Alpha, class Dst>
	void incrementImpl(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_lhs.increment(checkAliasing, alpha, dst); m_rhs.decrement(checkAliasing, alpha, dst); }
	
	template<typename Bool, typename Alpha, class Dst>
	void decrementImpl(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_lhs.decrement(checkAliasing, alpha, dst); m_rhs.increment(checkAliasing, alpha, dst); }
private:
	std::conditional_t<Lhs::isLeaf, const Lhs&, Lhs> m_lhs;
	std::conditional_t<Rhs::isLeaf, const Rhs&, Rhs> m_rhs;
};

template<class Lhs, class Rhs> 
MatrixSub<Lhs,Rhs> operator-(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) { return MatrixSub<Lhs,Rhs>(lhs, rhs); }

} // FSLinalg

#endif // FSLINALG_MATRIX_SUB_HPP
