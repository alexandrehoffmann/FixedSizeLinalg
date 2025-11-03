#ifndef FSLINALG_UNIT_MATRIX_HPP
#define FSLINALG_UNIT_MATRIX_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>

namespace FSLinalg
{
	
template<unsigned int Nrows, unsigned Ncols> class UnitMatrix;

template<unsigned int Nrows, unsigned Ncols>
struct MatrixTraits< UnitMatrix<Nrows, Ncols> >
{	
	using Scalar     = bool;
	using Size       = unsigned int;
	
	static constexpr bool hasReadRandomAccess  = true;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool hasFlatRandomAccess  = false;
	static constexpr bool causesAliasingIssues = false;
	static constexpr bool isLeaf               = true;
	
	static constexpr Size nRows = Nrows;
	static constexpr Size nCols = Ncols;
};

template<unsigned int Nrows, unsigned Ncols> 
class UnitMatrix : public MatrixBase< UnitMatrix<Nrows, Ncols> >
{
public:
	using Self = UnitMatrix<Nrows, Ncols>;
	FSLINALG_DEFINE_MATRIX
	
	struct Id { Size i; Size j; };
	
	UnitMatrix(const Size a_i, const Size a_j) : m_id({.i = a_i, .j = a_j}) {}
	
	const_ReturnType getImpl(const Size i, const Size j) const { return const_ReturnType(i == m_id.i and j == m_id.j); }
	
	Id getId() const { return m_id; }
	
	template<class Dst> constexpr bool isAliasedToImpl(const VectorBase<Dst>&) const { return false; }
	template<class Dst> constexpr bool isAliasedToImpl(const MatrixBase<Dst>&) const { return false; }
private:
	Id m_id;
};

} // namespace FSLinalg

#endif // FSLINALG_UNIT_MATRIX_HPP
