#ifndef FSLINALG_STRIP_SYMBOLS_FROM_VECTOR_OUTER_PRODUCT_HPP
#define FSLINALG_STRIP_SYMBOLS_FROM_VECTOR_OUTER_PRODUCT_HPP

#include <FSLinalg/Matrix/VectorOuterProduct.hpp>
#include <FSLinalg/Matrix/MatrixScale.hpp>
#include <FSLinalg/Matrix/MatrixMinus.hpp>
#include <FSLinalg/Matrix/MatrixConj.hpp>
#include <FSLinalg/Matrix/MatrixTransposed.hpp>

namespace FSLinalg
{

template<typename GeneralMatrix>     struct IsVectorOuterProduct                                : std::false_type {}; 
template<class Lhs, class Rhs>       struct IsVectorOuterProduct< VectorOuterProduct<Lhs,Rhs> > : std::true_type  {}; 
template<typename Alpha, class Expr> struct IsVectorOuterProduct< MatrixScale<Alpha,Expr> >     : std::bool_constant<IsVectorOuterProduct<Expr>::value> {};
template<class Expr>                 struct IsVectorOuterProduct< MatrixMinus<Expr> >           : std::bool_constant<IsVectorOuterProduct<Expr>::value> {};
template<class Expr>                 struct IsVectorOuterProduct< MatrixConj<Expr> >            : std::bool_constant<IsVectorOuterProduct<Expr>::value> {};
template<class Expr>                 struct IsVectorOuterProduct< MatrixTransposed<Expr> >      : std::bool_constant<IsVectorOuterProduct<Expr>::value> {};

template<class GeneralMatrix> class StripSymbolsFromVectorOuterProduct;

template<class Lhs, class Rhs>
class StripSymbolsFromVectorOuterProduct< VectorOuterProduct<Lhs,Rhs> >
{
public:
	using LhsVector = Lhs;
	using RhsVector = Rhs;
	using Scalar    = std::common_type_t<typename LhsVector::RealScalar, typename RhsVector::RealScalar>;
	
	static constexpr bool isConjugated = false;
	static constexpr bool isTransposed = false;
	
	StripSymbolsFromVectorOuterProduct(const VectorOuterProduct<Lhs,Rhs>& outerProduct) : m_outerProduct(outerProduct) {}
	
	const LhsVector& getLhsVector() const { return m_outerProduct.m_lhs; }
	const RhsVector& getRhsVector() const { return m_outerProduct.m_rhs; }
	
	constexpr Scalar getAlpha() const { return 1;  }
private:
	VectorOuterProduct<Lhs,Rhs> m_outerProduct;
};

template<typename Alpha, class Expr>
class StripSymbolsFromVectorOuterProduct< MatrixScale<Alpha,Expr> >
{
	static_assert(IsVectorOuterProduct<Expr>::value, "Expr must be an outer product between two vectors");
public:
	using LhsVector = typename StripSymbolsFromVectorOuterProduct<Expr>::LhsVector;
	using RhsVector = typename StripSymbolsFromVectorOuterProduct<Expr>::RhsVector;
	using Scalar    = decltype(std::declval<typename StripSymbolsFromVectorOuterProduct<Expr>::Scalar>() * std::declval<Alpha>());
	
	static constexpr bool isConjugated = StripSymbolsFromVectorOuterProduct<Expr>::isConjugated;
	static constexpr bool isTransposed = StripSymbolsFromVectorOuterProduct<Expr>::isTransposed;
	
	StripSymbolsFromVectorOuterProduct(const MatrixScale<Alpha,Expr>& scaled_expr) : m_expr(scaled_expr.m_expr), m_alpha(scaled_expr.m_alpha) {}
	
	const LhsVector& getLhsVector() const { return m_expr.getLhsVector(); }
	const RhsVector& getRhsVector() const { return m_expr.getRhsVector(); }
	
	Scalar getAlpha() const { return m_alpha*m_expr.getAlpha(); }
private:
	StripSymbolsFromVectorOuterProduct<Expr> m_expr;
	Alpha                                    m_alpha;
};

template<class Expr>
class StripSymbolsFromVectorOuterProduct< MatrixMinus<Expr> >
{
	static_assert(IsVectorOuterProduct<Expr>::value, "Expr must be an outer product between two vectors");
public:
	using LhsVector = typename StripSymbolsFromVectorOuterProduct<Expr>::LhsVector;
	using RhsVector = typename StripSymbolsFromVectorOuterProduct<Expr>::RhsVector;
	using Scalar    = typename StripSymbolsFromVectorOuterProduct<Expr>::Scalar;
	
	static constexpr bool isConjugated = StripSymbolsFromVectorOuterProduct<Expr>::isConjugated;
	static constexpr bool isTransposed = StripSymbolsFromVectorOuterProduct<Expr>::isTransposed;
	
	StripSymbolsFromVectorOuterProduct(const MatrixMinus<Expr>& minus_expr) : m_expr(minus_expr.m_expr) {}
	
	const LhsVector& getLhsVector() const { return m_expr.getLhsVector(); }
	const RhsVector& getRhsVector() const { return m_expr.getRhsVector(); }
	
	Scalar getAlpha() const { return -m_expr.getAlpha(); }
private:
	StripSymbolsFromVectorOuterProduct<Expr> m_expr;
};

template<class Expr>
class StripSymbolsFromVectorOuterProduct< MatrixConj<Expr> >
{
	static_assert(IsVectorOuterProduct<Expr>::value, "Expr must be an outer product between two vectors");
public:
	using LhsVector = typename StripSymbolsFromVectorOuterProduct<Expr>::LhsVector;
	using RhsVector = typename StripSymbolsFromVectorOuterProduct<Expr>::RhsVector;
	using Scalar    = typename StripSymbolsFromVectorOuterProduct<Expr>::Scalar;
	
	static constexpr bool isConjugated = not StripSymbolsFromVectorOuterProduct<Expr>::isConjugated;
	static constexpr bool isTransposed = StripSymbolsFromVectorOuterProduct<Expr>::isTransposed;
	
	StripSymbolsFromVectorOuterProduct(const MatrixConj<Expr>& conj_expr) : m_expr(conj_expr.m_expr) {}
	
	const LhsVector& getLhsVector() const { return m_expr.getLhsVector(); }
	const RhsVector& getRhsVector() const { return m_expr.getRhsVector(); }
	
	Scalar getAlpha() const { return conj(m_expr.getAlpha()); }
private:
	StripSymbolsFromVectorOuterProduct<Expr> m_expr;
};

template<class Expr>
class StripSymbolsFromVectorOuterProduct< MatrixTransposed<Expr> >
{
	static_assert(IsVectorOuterProduct<Expr>::value, "Expr must be an outer product between two vectors");
public:
	using LhsVector = typename StripSymbolsFromVectorOuterProduct<Expr>::LhsVector;
	using RhsVector = typename StripSymbolsFromVectorOuterProduct<Expr>::RhsVector;
	using Scalar    = typename StripSymbolsFromVectorOuterProduct<Expr>::Scalar;
	
	static constexpr bool isConjugated = StripSymbolsFromVectorOuterProduct<Expr>::isConjugated;
	static constexpr bool isTransposed = not StripSymbolsFromVectorOuterProduct<Expr>::isTransposed;
	
	StripSymbolsFromVectorOuterProduct(const MatrixTransposed<Expr>& transposed_expr) : m_expr(transposed_expr.m_expr) {}
	
	const LhsVector& getLhsVector() const { return m_expr.getLhsVector(); }
	const RhsVector& getRhsVector() const { return m_expr.getRhsVector(); }
	
	Scalar getAlpha() const { return m_expr.getAlpha(); }
private:
	StripSymbolsFromVectorOuterProduct<Expr> m_expr;
};
	
} // namespace FSLinalg

#endif // FSLINALG_STRIP_SYMBOLS_VECTOR_OUTER_PRODUCT_HPP
