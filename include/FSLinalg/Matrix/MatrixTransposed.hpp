#ifndef FSLINALG_MATRIX_TRANSPOSE_HPP
#define FSLINALG_MATRIX_TRANSPOSE_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>
#include <FSLinalg/Matrix/MatrixConj.hpp>

namespace FSLinalg
{

template<class GeneralMatrix> class StripSymbolsFromVectorOuterProduct;
template<class GeneralMatrix> class StripSymbolsAndEvalMatrix;
template<class Expr>          class MatrixTransposed;

template<class Expr> 
struct MatrixTraits< MatrixTransposed<Expr> >
{
	static_assert(IsMatrix<Expr>::value, "Expr must be a Matrix");
		
	using Scalar = typename Expr::Scalar;
	using Size   = typename Expr::Size;
	
	static constexpr bool hasReadRandomAccess  = Expr::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = Expr::hasWriteRandomAccess;
	static constexpr bool hasFlatRandomAccess  = Expr::hasFlatRandomAccess and (Expr::isRowVector or Expr::isColVector);
	static constexpr bool causesAliasingIssues = not (Expr::isRowVector or Expr::isColVector);
	static constexpr bool isLeaf               = false;
	
	static constexpr Size nRows = Expr::nCols;   
	static constexpr Size nCols = Expr::nRows;   
};

template<class Expr> 
class MatrixTransposed : public MatrixBase< MatrixTransposed<Expr> >
{
public:
	using Self = MatrixTransposed<Expr>;
	FSLINALG_DEFINE_MATRIX
	
	friend class StripSymbolsFromVectorOuterProduct<Self>; 
	friend class StripSymbolsAndEvalMatrix<Self>; 
	
	static constexpr bool isScalarComplex = IsComplexScalar<Scalar>::value;
	
	MatrixTransposed(const MatrixBase<Expr>& expr) : m_expr(expr.derived()) {}
	
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess)  { return m_expr.getImpl(j, i); }
	      ReturnType getImpl(const Size i, const Size j)       requires(hasWriteRandomAccess) { return m_expr.getImpl(j, i); }
	      
	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess  and hasFlatRandomAccess) { return m_expr.getImpl(i); }
	      ReturnType getImpl(const Size i)       requires(hasWriteRandomAccess and hasFlatRandomAccess) { return m_expr.getImpl(i); }

	template<class SrcExpr> MatrixTransposed& operator= (const MatrixBase<SrcExpr>& srcExpr) requires(IsConstructibleFrom<SrcExpr>::value) { srcExpr.assignTo  (1., *this, std::true_type{}); return *this; }
	template<class SrcExpr> MatrixTransposed& operator+=(const MatrixBase<SrcExpr>& srcExpr) requires(IsConstructibleFrom<SrcExpr>::value) { srcExpr.increment (1., *this, std::true_type{}); return *this; }
	template<class SrcExpr> MatrixTransposed& operator-=(const MatrixBase<SrcExpr>& srcExpr) requires(IsConstructibleFrom<SrcExpr>::value) { srcExpr.decrement (1., *this, std::true_type{}); return *this; }
	
	MatrixTransposed& operator*=(const RealScalar& alpha) requires(isScalarComplex);
	MatrixTransposed& operator/=(const RealScalar& alpha) requires(isScalarComplex);
	
	MatrixTransposed& operator*=(const Scalar& alpha);
	MatrixTransposed& operator/=(const Scalar& alpha);

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
MatrixTransposed<Expr> transpose(const MatrixBase<Expr>& expr) { return MatrixTransposed<Expr>(expr); }

template<class Expr> 
MatrixConj< MatrixTransposed<Expr> > adjoint(const MatrixBase<Expr>& expr) { return MatrixConj< MatrixTransposed<Expr> >(MatrixTransposed<Expr>(expr)); }

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_TRANSPOSE_HPP
