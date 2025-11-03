#ifndef FSLINALG_MATRIX_CONJ_HPP
#define FSLINALG_MATRIX_CONJ_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>

namespace FSLinalg
{

template<class Expr> class MatrixConj;
	
template<class Expr> 
struct MatrixTraits< MatrixConj<Expr> >
{
	static_assert(IsMatrix<Expr>::value, "Expr must be a Matrix");
	
	using ExprTraits = MatrixTraits<Expr>;
		
	using Scalar = typename ExprTraits::Scalar;
	using Size   = typename ExprTraits::Size;
	
	static constexpr bool hasReadRandomAccess  = ExprTraits::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool hasFlatRandomAccess  = ExprTraits::hasFlatRandomAccess;
	static constexpr bool causesAliasingIssues = ExprTraits::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size nRows = ExprTraits::nRows;   
	static constexpr Size nCols = ExprTraits::nCols;   
};

template<class Expr> 
class MatrixConj : public MatrixBase< MatrixConj<Expr> >
{
public:
	using Self = MatrixConj<Expr>;
	FSLINALG_DEFINE_MATRIX
	
	MatrixConj(const MatrixBase<Expr>& expr) : m_expr(expr.derived()) {}
	
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess)  { return conj(m_expr.getImpl(i,j)); }
	
	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess and hasFlatRandomAccess) { return conj(m_expr.getImpl(i)); }
	
	template<class Dst> bool isAliasedToImpl(const VectorBase<Dst>& other) const { return m_expr.isAliasedToImpl(other); }
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& other) const { return m_expr.isAliasedToImpl(other); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
private:
	std::conditional_t<Expr::isLeaf, const Expr&, Expr> m_expr;
};

template<class Expr> 
FSLinalg::MatrixConj<Expr> conj(const FSLinalg::MatrixBase<Expr>& expr) { return FSLinalg::MatrixConj<Expr>(expr); }

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_CONJ_HPP
