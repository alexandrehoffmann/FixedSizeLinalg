#ifndef FSLINALG_VECTOR_CROSS_IMPL_HPP
#define FSLINALG_VECTOR_CROSS_IMPL_HPP

#include <FSLinalg/Matrix/VectorCross.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> template<class Alpha, class LhsVec, class RhsVec, class DstVec> 
void VectorCross<Lhs, Rhs>::assignCross(const Alpha& alpha, const LhsVec& lhs, const RhsVec& rhs, DstVec& dst)
{
	dst[0] = alpha*(lhs[1]*rhs[2] - lhs[2]*rhs[1]);
	dst[1] = alpha*(lhs[2]*rhs[0] - lhs[0]*rhs[2]);
	dst[2] = alpha*(lhs[0]*rhs[1] - lhs[1]*rhs[0]);
}

template<class Lhs, class Rhs> template<class Alpha, class LhsVec, class RhsVec, class DstVec> 
void VectorCross<Lhs, Rhs>::incrementCross(const Alpha& alpha, const LhsVec& lhs, const RhsVec& rhs, DstVec& dst)
{
	dst[0] += alpha*(lhs[1]*rhs[2] - lhs[2]*rhs[1]);
	dst[1] += alpha*(lhs[2]*rhs[0] - lhs[0]*rhs[2]);
	dst[2] += alpha*(lhs[0]*rhs[1] - lhs[1]*rhs[0]);
}

template<class Lhs, class Rhs> template<class Alpha, class LhsVec, class RhsVec, class DstVec> 
void VectorCross<Lhs, Rhs>::decrementCross(const Alpha& alpha, const LhsVec& lhs, const RhsVec& rhs, DstVec& dst)
{
	dst[0] -= alpha*(lhs[1]*rhs[2] - lhs[2]*rhs[1]);
	dst[1] -= alpha*(lhs[2]*rhs[0] - lhs[0]*rhs[2]);
	dst[2] -= alpha*(lhs[0]*rhs[1] - lhs[1]*rhs[0]);
}

template<class Lhs, class Rhs> template<typename Alpha, class Dst, bool checkAliasing>
void VectorCross<Lhs, Rhs>::assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{	
	TmpLhs lhs(m_lhs);
	TmpRhs rhs(m_rhs);
	
	if (checkAliasing and (lhs.isAliasedTo(dst) or rhs.isAliasedTo(dst)))
	{
		RowVector<typename Dst::Scalar, Dst::size> tmp_dst;
		assignCross(alpha, lhs, rhs, tmp_dst);
		
		for (Size i=0; i!=size; ++i) { dst[i] = tmp_dst[i]; }
	}
	else
	{
		assignCross(alpha, lhs, rhs, dst);
	}
}

template<class Lhs, class Rhs> template<typename Alpha, class Dst, bool checkAliasing>
void VectorCross<Lhs, Rhs>::incrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	TmpLhs lhs(m_lhs);
	TmpRhs rhs(m_rhs);
	
	if (checkAliasing and (lhs.isAliasedTo(dst) or rhs.isAliasedTo(dst)))
	{
		RowVector<typename Dst::Scalar, Dst::size> tmp_dst;
		assignCross(alpha, lhs, rhs, tmp_dst);
		
		for (Size i=0; i!=size; ++i) { dst[i] += tmp_dst[i]; }
	}
	else
	{
		incrementCross(alpha, lhs, rhs, dst);
	}
}

template<class Lhs, class Rhs> template<typename Alpha, class Dst, bool checkAliasing>
void VectorCross<Lhs, Rhs>::decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	TmpLhs lhs(m_lhs);
	TmpRhs rhs(m_rhs);
	
	if (checkAliasing and (lhs.isAliasedTo(dst) or rhs.isAliasedTo(dst)))
	{
		RowVector<typename Dst::Scalar, Dst::size> tmp_dst;
		assignCross(alpha, lhs, rhs, tmp_dst);
		
		for (Size i=0; i!=size; ++i) { dst[i] -= tmp_dst[i]; }
	}
	else
	{
		decrementCross(alpha, lhs, rhs, dst);
	}
}

} // namespace FSLinalg

#endif // FSLINALG_VECTOR_CROSS_IMPL_HPP
