#ifndef FSLINALG_VECTOR_MINUS_HPP
#define FSLINALG_VECTOR_MINUS_HPP

#include <FSLinalg/Vector/VectorBase.hpp>

namespace FSLinalg
{

template<class Expr> class VectorMinus;

template<class Expr> 
struct VectorTraits< VectorMinus<Expr> >
{
	static_assert(IsScalar<Alpha>::value and IsVector<Expr>::value, "Alpha must be a Scalar and Expr must be a Vector");
	
	using ExprTraits = VectorTraits<Expr>;
	
	using Scalar     = typename ExprTraits::Scalar;
	using Size       = typename ExprTraits::Size;
	
	static constexpr bool hasReadRandomAccess  = ExprTraits::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool causesAliasingIssues = ExprTraits::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size size = ExprTraits::size;   
};

template<class Expr> 
class VectorMinus : public VectorBase< VectorMinus<Expr> >
{
public:
	using Self = VectorMinus<Expr>;
	FSLINALG_DEFINE_VECTOR
	
	VectorMinus(const VectorBase<Expr>& expr) : m_expr(expr.derived()) {}
	
	const_ReturnType getImpl(const Size i) const requires(hasReadRandomAccess) { return -m_expr.getImpl(i); }
	
	template<class Dst> bool isAliasedToImpl(const VectorBase<Dst>& dst) const { return m_expr.isAliasedToImpl(dst); }
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& dst) const { return m_expr.isAliasedToImpl(dst); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)  { m_expr.assignTo(-alpha, dst, bc); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.decrement(alpha, dst, bc); }
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { m_expr.increment(alpha, dst, bc); }
private:
	std::conditional_t<Expr::isLeaf, const Expr&, Expr> m_expr;
};

template<class Expr> 
FSLinalg::VectorMinus<Expr> operator-(const FSLinalg::VectorBase<Expr>& expr) { return FSLinalg::VectorMinus<Expr>(expr); }

} // namespace FSLinalg


#endif // FSLINALG_VECTOR_MINUS_HPP
