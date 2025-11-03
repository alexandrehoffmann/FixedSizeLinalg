#ifndef FSLINALG_MATRIX_TRANSPOSE_HPP
#define FSLINALG_MATRIX_TRANSPOSE_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>
#include <FSLinalg/Matrix/MatrixConj.hpp>

namespace FSLinalg
{

template<class Expr> class MatrixTransposed;

template<class Expr> 
struct MatrixTraits< MatrixTransposed<Expr> >
{
	static_assert(IsMatrix<Expr>::value, "Expr must be a Matrix");
	
	using ExprTraits = MatrixTraits<Expr>;
		
	using Scalar = typename ExprTraits::Scalar;
	using Size   = typename ExprTraits::Size;
	
	static constexpr bool hasReadRandomAccess  = ExprTraits::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = ExprTraits::hasWriteRandomAccess;
	static constexpr bool hasFlatRandomAccess  = false;
	static constexpr bool causesAliasingIssues = true;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size nRows = ExprTraits::nCols;   
	static constexpr Size nCols = ExprTraits::nRows;   
};

template<class Expr> 
class MatrixTransposed : public MatrixBase< MatrixTransposed<Expr> >
{
public:
	using Self = MatrixTransposed<Expr>;
	FSLINALG_DEFINE_MATRIX
	
	friend class StripSymbolsAndEvalMatrix<Self>; 
	
	static constexpr bool IsScalarComplex = IsComplexScalar<Scalar>::value;
	
	MatrixTransposed(const MatrixBase<Expr>& expr) : m_expr(expr.derived()) {}
	
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess)  { return m_expr.getImpl(j, i); }
	      ReturnType getImpl(const Size i, const Size j)       requires(hasWriteRandomAccess) { return m_expr.getImpl(j, i); }

	template<class SrcExpr> MatrixTransposed& operator= (const MatrixBase<SrcExpr>& srcExpr) requires(IsConstructibleFrom<SrcExpr>::value) { srcExpr.assignTo  (1., *this, std::true_type{}); return *this; }
	template<class SrcExpr> MatrixTransposed& operator+=(const MatrixBase<SrcExpr>& srcExpr) requires(IsConstructibleFrom<SrcExpr>::value) { srcExpr.increment (1., *this, std::true_type{}); return *this; }
	template<class SrcExpr> MatrixTransposed& operator-=(const MatrixBase<SrcExpr>& srcExpr) requires(IsConstructibleFrom<SrcExpr>::value) { srcExpr.decrement (1., *this, std::true_type{}); return *this; }
	
	MatrixTransposed& operator*=(const RealScalar& alpha) requires(IsScalarComplex);
	MatrixTransposed& operator/=(const RealScalar& alpha) requires(IsScalarComplex);
	
	MatrixTransposed& operator*=(const Scalar& alpha);
	MatrixTransposed& operator/=(const Scalar& alpha);

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
MatrixTransposed<Expr> transpose(const MatrixBase<Expr>& expr) { return MatrixTransposed<Expr>(expr); }

template<class Expr> 
MatrixConj< MatrixTransposed<Expr> > adjoint(const MatrixBase<Expr>& expr) { return MatrixConj< MatrixTransposed<Expr> >(MatrixTransposed<Expr>(expr)); }

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_TRANSPOSE_HPP
