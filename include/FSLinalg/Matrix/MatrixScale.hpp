#ifndef FSLINALG_MATRIX_SCALE_HPP
#define FSLINALG_MATRIX_SCALE_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>

namespace FSLinalg
{
	
template<class GeneralMatrix>        class StripSymbolsFromVectorOuterProduct;
template<class GeneralMatrix>        class StripSymbolsAndEvalMatrix;
template<typename Alpha, class Expr> class MatrixScale;

template<typename Alpha, class Expr> 
struct MatrixTraits< MatrixScale<Alpha,Expr> >
{
	static_assert(IsScalar<Alpha>::value and IsMatrix<Expr>::value, "Alpha must be a Scalar and Expr must be a Vector");
		
	using Scalar = decltype(std::declval<Alpha>() * std::declval<typename Expr::Scalar>());
	using Size   = typename Expr::Size;
	
	static constexpr bool hasReadRandomAccess  = Expr::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool hasFlatRandomAccess  = Expr::hasFlatRandomAccess;
	static constexpr bool causesAliasingIssues = Expr::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size nRows = Expr::nRows;   
	static constexpr Size nCols = Expr::nCols;   
};

template<typename Alpha, class Expr> 
class MatrixScale : public MatrixBase< MatrixScale<Alpha,Expr> >
{
public:
	using Self = MatrixScale<Alpha,Expr>;
	FSLINALG_DEFINE_MATRIX
	
	friend class StripSymbolsFromVectorOuterProduct<Self>; 
	friend class StripSymbolsAndEvalMatrix<Self>; 
	
	MatrixScale(const Alpha& alpha, const MatrixBase<Expr>&  expr) : m_alpha(alpha), m_expr(expr.derived()) {}
	
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess)  { return m_alpha*m_expr.getImpl(i, j); }

	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess and hasFlatRandomAccess) { return m_alpha*m_expr.getImpl(i); }
	
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& other) const { return m_expr.isAliasedToImpl(other); }

	template<typename Beta, class Dst, bool checkAliasing>
	void assignToImpl(const Beta& beta, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.assignTo(beta*m_alpha, dst, bc); }
	
	template<typename Beta, class Dst, bool checkAliasing>
	void incrementImpl(const Beta& beta, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.increment(beta*m_alpha, dst, bc); }
	
	template<typename Beta, class Dst, bool checkAliasing>
	void decrementImpl(const Beta& beta, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.decrement(beta*m_alpha, dst, bc); }
private:
	Alpha m_alpha;
	std::conditional_t<Expr::isLeaf, const Expr&, Expr> m_expr;
}; 

template<typename Alpha, class Expr> 
MatrixScale<Alpha, Expr> operator*(const Alpha& alpha, const MatrixBase<Expr>& expr) requires(IsScalar<Alpha>::value) { return MatrixScale<Alpha,Expr>(alpha, expr); }

template<typename Alpha, class Expr> 
MatrixScale<Alpha, Expr> operator*(const MatrixBase<Expr>& expr, const Alpha& alpha) requires(IsScalar<Alpha>::value) { return MatrixScale<Alpha,Expr>(alpha, expr); }

template<typename Alpha, class Expr> 
MatrixScale<Alpha, Expr> operator/(const MatrixBase<Expr>& expr, const Alpha& alpha) requires(IsScalar<Alpha>::value) { return MatrixScale<Alpha,Expr>(1. / alpha, expr); }

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_SCALE_HPP
