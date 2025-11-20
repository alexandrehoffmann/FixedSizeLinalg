#ifndef FSLINALG_KEEP_BRACKETING_HPP
#define FSLINALG_KEEP_BRACKETING_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>
#include <FSLinalg/Matrix/MatrixProduct.hpp>

namespace FSLinalg
{
	
template<class Expr> class KeepBrackets;

template<class Expr>
struct MatrixTraits< KeepBrackets<Expr> >
{
	static_assert(IsMatrixProduct<Expr>::value);
	
	using Scalar = typename Expr::Scalar;
	using Size   = typename Expr::Size;
	
	static constexpr bool hasReadRandomAccess  = Expr::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = Expr::hasWriteRandomAccess;
	static constexpr bool hasFlatRandomAccess  = Expr::hasFlatRandomAccess;
	static constexpr bool causesAliasingIssues = Expr::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size nRows = Expr::nRows;   
	static constexpr Size nCols = Expr::nCols;  
};
	
template<class Expr>
class KeepBrackets : public MatrixBase< KeepBrackets<Expr> >
{
public:
	using Self = KeepBrackets<Expr>;
	FSLINALG_DEFINE_MATRIX
	
	KeepBrackets(const MatrixBase<Expr>& expr) : m_expr(expr.derived()) {}
	
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess)  { return m_expr.getImpl(i, j); }
	      ReturnType getImpl(const Size i, const Size j)       requires(hasWriteRandomAccess) { return m_expr.getImpl(i, j); }
	      
	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess  and hasFlatRandomAccess) { return m_expr.getImpl(i); }
	      ReturnType getImpl(const Size i)       requires(hasWriteRandomAccess and hasFlatRandomAccess) { return m_expr.getImpl(i); }
	      
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& other) const { return m_expr.isAliasedToImpl(other); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.assignToHelper(alpha, dst, bc, std::true_type{}); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.incrementHelper(alpha, dst, bc, std::true_type{}); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.decrementHelper(alpha, dst, bc, std::true_type{}); }
private:
	std::conditional_t<Expr::isLeaf, const Expr&, Expr> m_expr;
};

template<class Expr>
KeepBrackets<Expr> keepBrackets(const MatrixBase<Expr>& expr) { return KeepBrackets<Expr>(expr); }
	
} // namespace FSLinalg

#endif // FSLINALG_KEEP_BRACKETING_HPP
