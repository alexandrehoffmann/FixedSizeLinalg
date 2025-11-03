#ifndef FSLINALG_STRIP_SYMBOLS_AND_EVAL_VECTOR_HPP
#define FSLINALG_STRIP_SYMBOLS_AND_EVAL_VECTOR_HPP

#include <FSLinalg/Vector/VectorBase.hpp>
#include <FSLinalg/Vector/Vector.hpp>
#include <FSLinalg/Vector/VectorScale.hpp>
#include <FSLinalg/Vector/VectorConj.hpp>

namespace FSLinalg
{

template<class Expr>
class StripSymbolsAndEvalVector
{
public:
	static_assert(IsVector<Expr>::value, "Expr must be a vector");

	using TmpVector = FSLinalg::Vector< typename Expr::Scalar, Expr::size >;
	using Vector = std::conditional_t<Expr::isLeaf, Expr, TmpVector>;
	using Scalar = typename Expr::RealScalar;
	
	static constexpr bool isConjugated = false;
	
	static constexpr unsigned int size = Expr::size;
	
	static constexpr bool createsTemporary = not Expr::isLeaf;
	
	StripSymbolsAndEvalVector(const VectorBase<Expr>& expr) : m_vector(expr.derived()) {}
	
	const     Vector& getVector() const { return m_vector; }
	constexpr Scalar  getAlpha()  const { return 1.; }
private:
	std::conditional_t<Expr::isLeaf, const Expr&, TmpVector> m_vector;
};

template<typename Alpha, class Expr> 
class StripSymbolsAndEvalVector< VectorScale<Alpha,Expr> >
{
public:
	using Vector = typename StripSymbolsAndEvalVector<Expr>::Vector;
	using Scalar = decltype( std::declval<Alpha>()*std::declval<typename StripSymbolsAndEvalVector<Expr>::Scalar>() );
	
	static constexpr bool isConjugated = StripSymbolsAndEvalVector<Expr>::isConjugated;
	
	static constexpr unsigned int size = StripSymbolsAndEvalVector<Expr>::size;
	
	static constexpr bool createsTemporary = StripSymbolsAndEvalVector<Expr>::createsTemporary;
	
	StripSymbolsAndEvalVector(const VectorScale<Alpha,Expr>& scaled_expr) : m_expr(scaled_expr.m_expr), m_alpha(scaled_expr.m_alpha) {}
	
	const Vector& getVector() const { return m_expr.getVector(); }
	
	Scalar getAlpha() const { return m_alpha*m_expr.getAlpha(); }
private:
	StripSymbolsAndEvalVector<Expr> m_expr;
	Alpha                           m_alpha;
};

template<class Expr>
class StripSymbolsAndEvalVector< VectorConj<Expr> >
{
public:
	using Vector = typename StripSymbolsAndEvalVector<Expr>::Vector;
	using Scalar = typename StripSymbolsAndEvalVector<Expr>::Scalar;
	
	static constexpr bool isConjugated = not StripSymbolsAndEvalVector<Expr>::isConjugated;
	
	static constexpr unsigned int size = StripSymbolsAndEvalVector<Expr>::size;
	
	static constexpr bool createsTemporary = StripSymbolsAndEvalVector<Expr>::createsTemporary;
	
	StripSymbolsAndEvalVector(const VectorConj<Expr>& conj_expr) : m_expr(conj_expr.m_expr) {}
	
	const Vector& getVector() const { return m_expr.getVector(); }
	
	Scalar  getAlpha()  const { return m_expr.getAlpha();  }
private:
	StripSymbolsAndEvalVector<Expr> m_expr;
};

} // namespace FSLinalg

#endif // FSLINALG_STRIP_SYMBOLS_AND_EVAL_VECTOR_HPP
