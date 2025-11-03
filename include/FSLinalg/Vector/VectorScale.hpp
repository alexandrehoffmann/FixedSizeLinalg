#ifndef FSLINALG_VECTOR_SCALE_HPP
#define FSLINALG_VECTOR_SCALE_HPP

#include <FSLinalg/Vector/VectorBase.hpp>

namespace FSLinalg
{

template<typename UndefinedVector>   class StripSymbolsAndEvalVector;
template<typename Alpha, class Expr> class VectorScale;

template<typename Alpha, class Expr> 
struct VectorTraits< VectorScale<Alpha,Expr> >
{
	static_assert(IsScalar<Alpha>::value and IsVector<Expr>::value, "Alpha must be a Scalar and Expr must be a Vector");
	
	using ExprTraits = VectorTraits<Expr>;
	
	using Scalar     = decltype(std::declval<Alpha>()*std::declval<typename ExprTraits::Scalar>());
	using Size       = typename ExprTraits::Size;
	
	static constexpr bool hasReadRandomAccess  = ExprTraits::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool causesAliasingIssues = ExprTraits::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size size = ExprTraits::size;   
};

template<typename Alpha, class Expr> 
class VectorScale : public VectorBase< VectorScale<Alpha, Expr> >
{
public:
	using Self = VectorScale<Alpha, Expr>;
	FSLINALG_DEFINE_VECTOR
	
	friend class StripSymbolsAndEvalVector<Self>; 
	
	VectorScale(const Alpha& alpha, const VectorBase<Expr>& expr) : m_alpha(alpha), m_expr(expr.derived()) {}

	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess) { return m_alpha*m_expr.getImpl(i); }
	
	template<class Dst> bool isAliasedToImpl(const VectorBase<Dst>& dst) const { return m_expr.isAliasedToImpl(dst); }
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& dst) const { return m_expr.isAliasedToImpl(dst); }

	template<typename Beta, class Dst, bool checkAliasing>
	void assignToImpl(const Beta& beta, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Beta>::value)  { m_expr.assignTo(beta*m_alpha, dst, bc); }
	
	template<typename Beta, class Dst, bool checkAliasing>
	void incrementImpl(const Beta& beta, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Beta>::value) { m_expr.increment(beta*m_alpha, dst, bc); }
	
	template<typename Beta, class Dst, bool checkAliasing>
	void decrementImpl(const Beta& beta, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Beta>::value) { m_expr.decrement(beta*m_alpha, dst, bc); }
private:
	Alpha m_alpha;
	std::conditional_t<Expr::isLeaf, const Expr&, Expr> m_expr;
};

template<typename Alpha, class Expr> 
FSLinalg::VectorScale<Alpha, Expr> operator*(const Alpha& alpha, const FSLinalg::VectorBase<Expr>& expr) requires(FSLinalg::IsScalar<Alpha>::value) { return FSLinalg::VectorScale<Alpha,Expr>(alpha, expr); }

template<typename Alpha, class Expr> 
FSLinalg::VectorScale<Alpha, Expr> operator*(const FSLinalg::VectorBase<Expr>& expr, const Alpha& alpha) requires(FSLinalg::IsScalar<Alpha>::value) { return FSLinalg::VectorScale<Alpha,Expr>(alpha, expr); }

template<typename Alpha, class Expr> 
FSLinalg::VectorScale<Alpha, Expr> operator/(const FSLinalg::VectorBase<Expr>& expr, const Alpha& alpha) requires(FSLinalg::IsScalar<Alpha>::value) { return FSLinalg::VectorScale<Alpha,Expr>(1. / alpha, expr); }

} // namespace FSLinalg

#endif // FSLINALG_VECTOR_SCALE_HPP
