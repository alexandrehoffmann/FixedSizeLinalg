#ifndef FSLINALG_MATRIX_MINUS_HPP
#define FSLINALG_MATRIX_MINUS_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>

namespace FSLinalg
{
	
template<class GeneralMatrix> class StripSymbolsFromVectorOuterProduct;
template<class GeneralMatrix> class StripSymbolsAndEvalMatrix;
template<class Expr>          class MatrixMinus;

template<class Expr> 
struct MatrixTraits< MatrixMinus<Expr> >
{
	static_assert(IsMatrix<Expr>::value, "Expr must be a Matrix");
	
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
class MatrixMinus : public MatrixBase< MatrixMinus<Expr> >
{
public:
	using Self = MatrixMinus<Expr>;
	FSLINALG_DEFINE_MATRIX
	
	friend class StripSymbolsAndEvalMatrix<Self>; 
	friend class StripSymbolsFromVectorOuterProduct<Self>; 
	
	MatrixMinus(const MatrixBase<Expr>&  expr) : m_expr(expr.derived()) {}
	
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess)  { return -m_expr.getImpl(i, j); }

	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess and hasFlatRandomAccess) { return -m_expr.getImpl(i); }
	
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& other) const { return m_expr.isAliasedToImpl(other); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.assignTo(-alpha, dst, bc); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.decrement(alpha, dst, bc); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.increment(alpha, dst, bc); }
private:
	std::conditional_t<Expr::isLeaf, const Expr&, Expr> m_expr;
}; 

template<class Expr> 
MatrixMinus<Expr> operator-(const MatrixBase<Expr>& expr) { return MatrixMinus<Expr>(expr); }

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_MINUS_HPP
