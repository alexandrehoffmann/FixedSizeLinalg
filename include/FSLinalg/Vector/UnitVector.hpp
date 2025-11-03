#ifndef FSLINALG_UNIT_VECTOR_HPP
#define FSLINALG_UNIT_VECTOR_HPP

#include <FSLinalg/Vector/VectorBase.hpp>

namespace FSLinalg
{

template<unsigned int N> class UnitVector;

template<unsigned int N>
struct VectorTraits< UnitVector<N> >
{	
	using Scalar = bool;
	using Size   = unsigned int;
	
	static constexpr bool hasReadRandomAccess  = true;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool causesAliasingIssues = false;
	static constexpr bool isLeaf               = true;
	
	static constexpr Size size = N;   
};

template<unsigned int N>
class UnitVector : public VectorBase< UnitVector<N> >
{
public:
	using Self = UnitVector<N>;
	FSLINALG_DEFINE_VECTOR
	
	UnitVector(const Size i) : m_i(i) {}
	
	const_ReturnType getImpl(const Size i) const { return const_ReturnType(i == m_i); }
	
	Size getId() const { return m_i; }
	
	template<class Dst> constexpr bool isAliasedToImpl(const VectorBase<Dst>&) const { return false; }
	template<class Dst> constexpr bool isAliasedToImpl(const MatrixBase<Dst>&) const { return false; }
private:
	Size m_i;
};	

} // namespace FSLinalg

#endif // FSLINALG_UNIT_VECTOR_HPP
