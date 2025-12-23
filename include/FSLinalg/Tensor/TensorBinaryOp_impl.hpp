#ifndef FSLINALG_TENSOR_BINARY_OP_IMPL_HPP
#define FSLINALG_TENSOR_BINARY_OP_IMPL_HPP

#include <FSLinalg/Tensor/TensorBinaryOp.hpp>
#include <FSLinalg/Tensor/Tensor.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs, class BinaryOp> template<typename Bool, typename Alpha, class Dst>
void TensorBinaryOp<Lhs, Rhs, BinaryOp>::assignToImpl(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{	
	if constexpr (isSum)
	{
		m_lhs.assignTo(checkAliasing, alpha, dst); 
		m_rhs.increment(checkAliasing, alpha, dst);
	}
	else if constexpr (isSub)
	{
		m_lhs.assignTo(checkAliasing, alpha, dst); 
		m_rhs.decrement(checkAliasing, alpha, dst);
	}
	else if constexpr (isMul)
	{
		m_lhs.assignTo(checkAliasing, alpha, dst); 
		m_rhs.multiply(checkAliasing, BIC::fixed<RealScalar, RealScalar(1)>, dst);
	}
	else if constexpr (isDiv)
	{
		m_lhs.assignTo(checkAliasing, alpha, dst); 
		m_rhs.divide(checkAliasing, BIC::fixed<RealScalar, RealScalar(1)>, dst);
	}
	else
	{
		using TmpLhs = std::conditional_t<Lhs::hasReadRandomAccess, const Lhs&, TensorFromShape<typename Lhs::Scalar, Lhs::shape> >;
		using TmpRhs = std::conditional_t<Rhs::hasReadRandomAccess, const Rhs&, TensorFromShape<typename Rhs::Scalar, Rhs::shape> >;
		
		TmpLhs lhs(m_lhs);
		TmpLhs rhs(m_rhs);
		
		if (checkAliasing and (lhs.isAliasedTo(dst) or rhs.isAliasedTo(dst)))
		{
			TensorFromShape<typename Dst::Scalar, Dst::shape> tmp_dst(dst);
			
			if constexpr (TmpLhs::hasFlatRandomAccess and TmpRhs::hasFlatRandomAccess)
			{
				for (Size i=0; i!=size; ++i) { tmp_dst[i] = assignOp(tmp_dst[i], alpha*m_op(lhs[i], rhs[i])); }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					tmp_dst(index) = assignOp(tmp_dst(index), alpha*m_op(lhs(index), rhs(index)));
				});
			}
		}
		else
		{
			if constexpr (TmpLhs::hasFlatRandomAccess and TmpRhs::hasFlatRandomAccess)
			{
				for (Size i=0; i!=size; ++i) { dst[i] = assignOp(dst[i], alpha*m_op(lhs[i], rhs[i])); }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) = assignOp(dst(index), alpha*m_op(lhs(index), rhs(index)));
				});
			}
		}
	}
}
	
template<class Lhs, class Rhs, class BinaryOp> template<typename Bool, typename Alpha, class Dst>
void TensorBinaryOp<Lhs, Rhs, BinaryOp>::incrementImpl(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	if constexpr (isSum)
	{
		m_lhs.increment(checkAliasing, alpha, dst); 
		m_rhs.increment(checkAliasing, alpha, dst);
	}
	else if constexpr (isSub)
	{
		m_lhs.increment(checkAliasing, alpha, dst); 
		m_rhs.decrement(checkAliasing, alpha, dst);
	}
	else
	{
		TensorFromShape<Scalar, shape> tmp;
		assignTo(BIC::fixed<bool,false>, BIC::fixed<RealScalar,RealScalar(1)>, tmp);
		
		tmp.increment(BIC::fixed<bool,false>, alpha, dst);
	}
}
	
template<class Lhs, class Rhs, class BinaryOp> template<typename Bool, typename Alpha, class Dst>
void TensorBinaryOp<Lhs, Rhs, BinaryOp>::decrementImpl(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	if constexpr (isSum)
	{
		m_lhs.decrement(checkAliasing, alpha, dst); 
		m_rhs.decrement(checkAliasing, alpha, dst);
	}
	else if constexpr (isSub)
	{
		m_lhs.decrement(checkAliasing, alpha, dst); 
		m_rhs.increment(checkAliasing, alpha, dst);
	}
	else
	{
		TensorFromShape<Scalar, shape> tmp;
		assignTo(BIC::fixed<bool,false>, BIC::fixed<RealScalar,RealScalar(1)>, tmp);
		
		tmp.decrement(BIC::fixed<bool,false>, alpha, dst);
	}
}
	
template<class Lhs, class Rhs, class BinaryOp> template<typename Bool, typename Alpha, class Dst>
void TensorBinaryOp<Lhs, Rhs, BinaryOp>::multiplyImpl(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	// dst *= alpha*(lhs op rhs)
	if constexpr (isMul)
	{
		m_lhs.multiply(checkAliasing, alpha, dst); 
		m_rhs.multiply(checkAliasing, BIC::fixed<RealScalar,RealScalar(1)>, dst);
	}
	else if constexpr (isDiv)
	{
		m_lhs.multiply(checkAliasing, alpha, dst); 
		m_rhs.divide(checkAliasing, BIC::fixed<RealScalar,RealScalar(1)>, dst);
	}
	else
	{
		TensorFromShape<Scalar, shape> tmp;
		assignTo(BIC::fixed<bool,false>, BIC::fixed<RealScalar,RealScalar(1)>, tmp);
		
		tmp.multiply(BIC::fixed<bool,false>, alpha, dst);
	}
}
	
template<class Lhs, class Rhs, class BinaryOp> template<typename Bool, typename Alpha, class Dst>
void TensorBinaryOp<Lhs, Rhs, BinaryOp>::divideImpl(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	constexpr BIC::Fixed<RealScalar,RealScalar(1)> one;
	
	// dst /= alpha*(lhs op rhs)
	if constexpr (isMul)
	{
		m_lhs.divide(checkAliasing, one / alpha, dst); 
		m_rhs.divide(checkAliasing, BIC::fixed<RealScalar,RealScalar(1)>, dst);
	}
	else if constexpr (isDiv)
	{
		// dst = rhs*dst / (alpha*lhs)
		m_lhs.divide(checkAliasing, one / alpha, dst); 
		m_rhs.multiply(checkAliasing, BIC::fixed<RealScalar,RealScalar(1)>, dst);
	}
	else
	{
		TensorFromShape<Scalar, shape> tmp;
		assignTo(BIC::fixed<bool,false>, BIC::fixed<RealScalar,RealScalar(1)>, tmp);
		
		tmp.divide(BIC::fixed<bool,false>, alpha, dst);
	}
}

} // namespace FSLinalg

#endif // FSLINALG_TENSOR_BINARY_OP_IMPL_HPP
