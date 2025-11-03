#ifndef FSLINALG_MATRIX_VECTOR_PRODUCT_HPP
#define FSLINALG_MATRIX_VECTOR_PRODUCT_HPP

#include <FSLinalg/Vector/VectorBase.hpp>
#include <FSLinalg/Vector/StripSymbolsAndEvalVector.hpp>

#include <FSLinalg/Matrix/MatrixBase.hpp>
#include <FSLinalg/Matrix/StripSymbolsAndEvalMatrix.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> class MatrixVectorProduct;

template<class Lhs, class Rhs>
struct VectorTraits< MatrixVectorProduct<Lhs,Rhs> >
{
	static_assert(IsMatrix<Lhs>::value, "LHS must be a matrix");
	static_assert(IsVector<Rhs>::value, "RHS must be vectors");
	
	using LhsTraits = MatrixTraits<Lhs>;
	using RhsTraits = VectorTraits<Rhs>;
	
	static_assert(LhsTraits::nCols == RhsTraits::size, "Matrix and Vector sizes must match");
	
	using Scalar = decltype(std::declval<typename LhsTraits::Scalar>() * std::declval<typename RhsTraits::Scalar>());
	using Size   = std::common_type_t<typename LhsTraits::Size, typename RhsTraits::Size>;
	
	static constexpr bool hasReadRandomAccess  = false;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool causesAliasingIssues = true;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size size = LhsTraits::nRows;
};

template<class Lhs, class Rhs> 
class MatrixVectorProduct : public VectorBase< MatrixVectorProduct<Lhs, Rhs> >
{
public:
	using Self = MatrixVectorProduct<Lhs, Rhs>;
	FSLINALG_DEFINE_VECTOR
	
	static constexpr bool createTemporaryMatrix = StripSymbolsAndEvalMatrix<Lhs>::createsTemporary;
	static constexpr bool createTemporaryVector = StripSymbolsAndEvalVector<Rhs>::createsTemporary;
	
	MatrixVectorProduct(const MatrixBase<Lhs>& matrix, const VectorBase<Rhs>& vector) : m_matrix(matrix.derived()), m_vector(vector.derived()) {}
	
	template<class Dst> bool isAliasedToImpl(const VectorBase<Dst>& dst) const { return m_matrix.isAliasedToImpl(dst) or m_vector.isAliasedToImpl(dst); }
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& dst) const { return m_matrix.isAliasedToImpl(dst) or m_vector.isAliasedToImpl(dst); }

	template<typename Beta, class Dst, bool checkAliasing>
	void assignToImpl(const Beta& beta, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Beta>::value);
	
	template<typename Beta, class Dst, bool checkAliasing>
	void incrementImpl(const Beta& beta, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Beta>::value);
	
	template<typename Beta, class Dst, bool checkAliasing>
	void decrementImpl(const Beta& beta, VectorBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Beta>::value);
private:
	using StrippedMat = typename StripSymbolsAndEvalMatrix<Lhs>::Matrix;
	using StrippedVec = typename StripSymbolsAndEvalVector<Rhs>::Vector;
	
	static constexpr bool isMatTransposed = StripSymbolsAndEvalMatrix<Lhs>::isTransposed;
	static constexpr bool isMatConjugated = StripSymbolsAndEvalMatrix<Lhs>::isConjugated;
	static constexpr bool isVecConjugated = StripSymbolsAndEvalVector<Rhs>::isConjugated;
	
	static constexpr unsigned int MatNRows = StripSymbolsAndEvalMatrix<Lhs>::nRows;
	static constexpr unsigned int MatNCols = StripSymbolsAndEvalMatrix<Lhs>::nCols;

	std::conditional_t<Lhs::isLeaf, const Lhs&, Lhs> m_matrix;
	std::conditional_t<Rhs::isLeaf, const Rhs&, Rhs> m_vector;
};

template<class Lhs, class Rhs> 
FSLinalg::MatrixVectorProduct<Lhs,Rhs> operator*(const FSLinalg::MatrixBase<Lhs>& A, const FSLinalg::VectorBase<Rhs>& x) { return FSLinalg::MatrixVectorProduct<Lhs,Rhs>(A, x); }

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_VECTOR_PRODUCT_HPP
