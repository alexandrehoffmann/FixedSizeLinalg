#ifndef FSLINALG_TENSOR_BINARY_OP_HPP
#define FSLINALG_TENSOR_BINARY_OP_HPP

#include <FSLinalg/Tensor/TensorBase.hpp>
#include <FSLinalg/Tensor/TensorBase.hpp>
#include <FSLinalg/misc/BinaryOp.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs, class Op> class TensorBinaryOp;

template<class Lhs, class Rhs, class Op> 
struct TensorTraits< TensorBinaryOp<Lhs, Rhs, Op> >
{		
	static_assert(IsTensor<Lhs>::value and IsTensor<Rhs>::value, "Both LHS and RHS must be tensors");	
	static_assert(Lhs::rank  == Rhs::rank, "Tensors rank must match");
	static_assert(Lhs::shape == Rhs::shape, "Tensors shape must match");
	static_assert(std::is_invocable<Op, typename Lhs::Scalar, typename Rhs::Scalar>::value, "Op must be a binary op");
	
	using Scalar = decltype(std::declval<Op>()(std::declval<typename Lhs::Scalar>(), std::declval<typename Rhs::Scalar>()));
	using Size   = std::common_type_t<typename Lhs::Size, typename Rhs::Size>;
	using Shape  = std::array<Size, Lhs::rank>;
	
	static constexpr bool hasReadRandomAccess  = Lhs::hasReadRandomAccess and Rhs::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool hasFlatRandomAccess  = Lhs::hasFlatRandomAccess and Rhs::hasFlatRandomAccess;
	static constexpr bool causesAliasingIssues = Lhs::causesAliasingIssues or Rhs::causesAliasingIssues;
	static constexpr bool isLeaf               = false;
	
	static constexpr Shape shape = Lhs::shape;
};

template<class Lhs, class Rhs, class Op> 
class TensorBinaryOp : public TensorBase< TensorBinaryOp<Lhs, Rhs, Op> >
{
public:
	using Self = TensorBinaryOp<Lhs, Rhs, Op>;
	FSLINALG_DEFINE_TENSOR	
	
	static constexpr bool isSum = std::is_same<Op, BinaryOp::Add>::value;
	static constexpr bool isSub = std::is_same<Op, BinaryOp::Sub>::value;
	static constexpr bool isMul = std::is_same<Op, BinaryOp::Mul>::value;
	static constexpr bool isDiv = std::is_same<Op, BinaryOp::Div>::value;
	
	TensorBinaryOp(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) : m_lhs(lhs.derived()), m_rhs(rhs.derived()) {}
	
	template<std::integral... Idx> 
	const_ReturnType getImpl(const Idx... idx) const requires(sizeof...(Idx) == rank and hasReadRandomAccess)  { return m_op(m_lhs(idx...), m_rhs(idx...)); }
	const_ReturnType getImpl(const Size i)     const requires(hasReadRandomAccess  and hasFlatRandomAccess)    { return m_op(m_lhs[i], m_rhs[i]);           }
	
	template<class Dst> bool isAliasedToImpl(const TensorBase<Dst>& other) const { return m_lhs.isAliasedToImpl(other) or m_rhs.isAliasedToImpl(other); }
	
	template<typename Bool, typename Alpha, class Dst>
	void assignToImpl(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Bool, typename Alpha, class Dst>
	void incrementImpl(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Bool, typename Alpha, class Dst>
	void decrementImpl(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Bool, typename Alpha, class Dst>
	void multiplyImpl(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Bool, typename Alpha, class Dst>
	void divideImpl(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
private:
	std::conditional_t<Lhs::isLeaf, const Lhs&, Lhs> m_lhs;
	std::conditional_t<Rhs::isLeaf, const Rhs&, Rhs> m_rhs;
	Op                                               m_op;
};

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::Add> operator+(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::Add>(lhs, rhs); }

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::Sub> operator-(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::Sub>(lhs, rhs); }

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::Mul> operator*(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::Mul>(lhs, rhs); }

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::Div> operator/(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::Div>(lhs, rhs); }

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::Mod> operator%(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::Mod>(lhs, rhs); }

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::Equal> operator==(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::Equal>(lhs, rhs); }

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::NotEqual> operator!=(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::NotEqual>(lhs, rhs); }

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::Greater> operator>(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::Greater>(lhs, rhs); }

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::GreaterOrEqual> operator>=(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::GreaterOrEqual>(lhs, rhs); }

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::Lower> operator<(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::Lower>(lhs, rhs); }

template<class Lhs, class Rhs> 
TensorBinaryOp<Lhs,Rhs,BinaryOp::LowerOrEqual> operator<=(const TensorBase<Lhs>& lhs, const TensorBase<Rhs>& rhs) { return TensorBinaryOp<Lhs,Rhs,BinaryOp::LowerOrEqual>(lhs, rhs); }

} // namespace FSLinalg

#endif // FSLINALG_TENSOR_BINARY_OP_HPP
