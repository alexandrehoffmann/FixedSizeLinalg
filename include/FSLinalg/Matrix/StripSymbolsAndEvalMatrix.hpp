#ifndef FSLINALG_STRIP_SYMBOLS_AND_EVAL_MATRIX_HPP
#define FSLINALG_STRIP_SYMBOLS_AND_EVAL_MATRIX_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>
#include <FSLinalg/Matrix/Matrix.hpp>
#include <FSLinalg/Matrix/MatrixScale.hpp>
#include <FSLinalg/Matrix/MatrixMinus.hpp>
#include <FSLinalg/Matrix/MatrixConj.hpp>
#include <FSLinalg/Matrix/MatrixTransposed.hpp>

namespace FSLinalg
{

template<class Expr>
class StripSymbolsAndEvalMatrix
{
public:
	static_assert(IsMatrix<Expr>::value, "Expr must be a matrix");

	using TmpMatrix = FSLinalg::Matrix< typename Expr::Scalar, Expr::nRows, Expr::nCols >;
	using Matrix = std::conditional_t<Expr::isLeaf, Expr, TmpMatrix>;
	using Scalar = BIC::Fixed<typename Expr::RealScalar, typename Expr::RealScalar(1)>;
	
	static constexpr bool isConjugated = false;
	static constexpr bool isTransposed = false;
	
	static constexpr unsigned int nRows = Expr::nRows;
	static constexpr unsigned int nCols = Expr::nCols;
	
	static constexpr bool createsTemporary = not Expr::isLeaf;
	
	StripSymbolsAndEvalMatrix(const MatrixBase<Expr>& expr) : m_matrix(expr.derived()) {}
	
	const     Matrix& getMatrix() const { return m_matrix; }
	constexpr Scalar  getAlpha()  const { return {}; }
private:
	std::conditional_t<Expr::isLeaf, const Expr&, TmpMatrix> m_matrix;
};

template<typename Alpha, class Expr> 
class StripSymbolsAndEvalMatrix< MatrixScale<Alpha,Expr> >
{
public:
	using Matrix = typename StripSymbolsAndEvalMatrix<Expr>::Matrix;
	using Scalar = decltype(std::declval<Alpha>() * std::declval<typename StripSymbolsAndEvalMatrix<Expr>::Scalar>());
	
	static constexpr bool isConjugated = StripSymbolsAndEvalMatrix<Expr>::isConjugated;
	static constexpr bool isTransposed = StripSymbolsAndEvalMatrix<Expr>::isTransposed;
	
	static constexpr unsigned int nRows = StripSymbolsAndEvalMatrix<Expr>::nRows;
	static constexpr unsigned int nCols = StripSymbolsAndEvalMatrix<Expr>::nCols;
	
	static constexpr bool createsTemporary = StripSymbolsAndEvalMatrix<Expr>::createsTemporary;
	
	StripSymbolsAndEvalMatrix(const MatrixScale<Alpha,Expr>& scaled_expr) : m_expr(scaled_expr.m_expr), m_alpha(scaled_expr.m_alpha) {}
	
	const Matrix& getMatrix() const { return m_expr.getMatrix(); }
	
	Scalar getAlpha() const { return m_alpha*m_expr.getAlpha(); }
private:
	StripSymbolsAndEvalMatrix<Expr> m_expr;
	Alpha                           m_alpha;
};

template<class Expr> 
class StripSymbolsAndEvalMatrix< MatrixMinus<Expr> >
{
public:
	using Matrix = typename StripSymbolsAndEvalMatrix<Expr>::Matrix;
	using Scalar = typename StripSymbolsAndEvalMatrix<Expr>::Scalar;
	
	static constexpr bool isConjugated = StripSymbolsAndEvalMatrix<Expr>::isConjugated;
	static constexpr bool isTransposed = StripSymbolsAndEvalMatrix<Expr>::isTransposed;
	
	static constexpr unsigned int nRows = StripSymbolsAndEvalMatrix<Expr>::nRows;
	static constexpr unsigned int nCols = StripSymbolsAndEvalMatrix<Expr>::nCols;
	
	static constexpr bool createsTemporary = StripSymbolsAndEvalMatrix<Expr>::createsTemporary;
	
	StripSymbolsAndEvalMatrix(const MatrixMinus<Expr>& minus_expr) : m_expr(minus_expr.m_expr) {}
	
	const Matrix& getMatrix() const { return m_expr.getMatrix(); }
	
	Scalar getAlpha() const { return -m_expr.getAlpha(); }
private:
	StripSymbolsAndEvalMatrix<Expr> m_expr;
};

template<class Expr>
class StripSymbolsAndEvalMatrix< MatrixConj<Expr> >
{
public:
	using Matrix = typename StripSymbolsAndEvalMatrix<Expr>::Matrix;
	using Scalar = typename StripSymbolsAndEvalMatrix<Expr>::Scalar;
	
	static constexpr bool isConjugated = not StripSymbolsAndEvalMatrix<Expr>::isConjugated;
	static constexpr bool isTransposed = StripSymbolsAndEvalMatrix<Expr>::isTransposed;
	
	static constexpr unsigned int nRows = StripSymbolsAndEvalMatrix<Expr>::nRows;
	static constexpr unsigned int nCols = StripSymbolsAndEvalMatrix<Expr>::nCols;
	
	static constexpr bool createsTemporary = StripSymbolsAndEvalMatrix<Expr>::createsTemporary;
	
	StripSymbolsAndEvalMatrix(const MatrixConj<Expr>& conj_expr) : m_expr(conj_expr.m_expr) {}
	
	const Matrix& getMatrix() const { return m_expr.getMatrix(); }
	
	Scalar getAlpha()  const { return conj(m_expr.getAlpha());  }
private:
	StripSymbolsAndEvalMatrix<Expr> m_expr;
};

template<class Expr>
class StripSymbolsAndEvalMatrix< MatrixTransposed<Expr> >
{
public:
	using Matrix = typename StripSymbolsAndEvalMatrix<Expr>::Matrix;
	using Scalar = typename StripSymbolsAndEvalMatrix<Expr>::Scalar;
	
	static constexpr bool isConjugated = StripSymbolsAndEvalMatrix<Expr>::isConjugated;
	static constexpr bool isTransposed = not StripSymbolsAndEvalMatrix<Expr>::isTransposed;
	
	static constexpr unsigned int nRows = StripSymbolsAndEvalMatrix<Expr>::nRows;
	static constexpr unsigned int nCols = StripSymbolsAndEvalMatrix<Expr>::nCols;
	
	static constexpr bool createsTemporary = StripSymbolsAndEvalMatrix<Expr>::createsTemporary;
	
	StripSymbolsAndEvalMatrix(const MatrixTransposed<Expr>& transposed_expr) : m_expr(transposed_expr.m_expr) {}
	
	const Matrix& getMatrix() const { return m_expr.getMatrix(); }
	
	Scalar getAlpha() const { return m_expr.getAlpha();  }
private:
	StripSymbolsAndEvalMatrix<Expr> m_expr;
};

} // namespace FSLinalg

#endif // FSLINALG_STRIP_SYMBOLS_AND_EVAL_MATRIX_HPP
