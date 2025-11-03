#ifndef FSLINALG_VECTOR_CONJ_HPP
#define FSLINALG_VECTOR_CONJ_HPP

#include <FSLinalg/Vector/VectorBase.hpp>

namespace FSLinalg
{

template<class Expr> class VectorConj;

template<class Expr> 
struct VectorTraits< VectorConj<Expr> >
{
	static_assert(IsVector<Expr>::value, "Expr must be a Vector");
	
	using ExprTraits = VectorTraits<Expr>;
	
	using Scalar     = typename ExprTraits::Scalar;
	using Size       = typename ExprTraits::Size;
	
	static constexpr bool hasReadRandomAccess  = ExprTraits::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = IsRealScalar<Scalar>::value;
	static constexpr bool causesAliasingIssues = ExprTraits::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size size = ExprTraits::size;   
};

template<class Expr> 
class VectorConj : public VectorBase< VectorConj<Expr> >
{
public:
	using Self = VectorConj<Expr>;
	FSLINALG_DEFINE_VECTOR
	
	friend class StripSymbolsAndEvalVector<Self>; 
	
	VectorConj(const VectorBase<Expr>& expr) : m_expr(expr.derived()) {}
	
	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess) { return conj(m_expr.getImpl(i)); }
	
	template<class Dst> bool isAliasedToImpl(const VectorBase<Dst>& dst) const { return m_expr.isAliasedToImpl(dst); }
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& dst) const { return m_expr.isAliasedToImpl(dst); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
private:
	std::conditional_t<Expr::isLeaf, const Expr&, Expr> m_expr;
};

template<class Expr> 
VectorConj<Expr> conj(const VectorBase<Expr>& expr) { return VectorConj<Expr>(expr); }

} // namespace FSLinalg

#endif // FSLINALG_VECTOR_CONJ_HPP
