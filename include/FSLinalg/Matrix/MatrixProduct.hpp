#ifndef FSLINALG_MATRIX_MATRIX_PRODUCT_HPP
#define FSLINALG_MATRIX_MATRIX_PRODUCT_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>
#include <FSLinalg/Matrix/MatrixTransposed.hpp>
#include <FSLinalg/Matrix/StripSymbolsAndEvalMatrix.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs> class MatrixProduct;

template<class Lhs, class Rhs>
struct MatrixTraits< MatrixProduct<Lhs, Rhs> >
{
	static_assert(IsMatrix<Lhs>::value and IsMatrix<Rhs>::value, "Both LHS and RHS must be matrices");
	static_assert(Lhs::nCols == Rhs::nRows, "Matrices size must match");
		
	using Scalar = decltype(std::declval<typename Lhs::Scalar>() * std::declval<typename Rhs::Scalar>());
	using Size   = std::common_type_t<typename Lhs::Size, typename Rhs::Size>;
	
	static constexpr bool hasReadRandomAccess  = Lhs::isRowVector and Rhs::isColVector and Lhs::hasFlatRandomAccess and Rhs::hasFlatRandomAccess;
	static constexpr bool hasWriteRandomAccess = false;
	static constexpr bool hasFlatRandomAccess  = false;
	static constexpr bool causesAliasingIssues = true;
	static constexpr bool isLeaf               = false;
	
	static constexpr Size nRows = Lhs::nRows;   
	static constexpr Size nCols = Rhs::nCols;   
};

template<class Lhs, class Rhs>
class MatrixProduct : public MatrixBase< MatrixProduct<Lhs,Rhs> >
{
public:
	using Self = MatrixProduct<Lhs,Rhs>;
	FSLINALG_DEFINE_MATRIX
	
	MatrixProduct(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) : m_lhs(lhs.derived()), m_rhs(rhs.derived()) {}
	
	static constexpr bool createTemporaryLhs = StripSymbolsAndEvalMatrix<Lhs>::createsTemporary;
	static constexpr bool createTemporaryRhs = StripSymbolsAndEvalMatrix<Rhs>::createsTemporary;
	
	/**
	 * @brief Only awailable when multiplying row-vector and col-vector
	 * When multiplying a row-vector and a col-vector, we can compute A*B as A(i,0)*B(0,j)
	 */
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess)  { return m_lhs.getImpl(i)*m_rhs.getImpl(j); }
	
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& other) const { return m_lhs.isAliasedToImpl(other) or m_rhs.isAliasedToImpl(other); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignToImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void incrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrementImpl(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
private:
	using StrippedLhs = typename StripSymbolsAndEvalMatrix<Lhs>::Matrix;
	using StrippedRhs = typename StripSymbolsAndEvalMatrix<Rhs>::Matrix;
	
	static constexpr bool isLhsTransposed = StripSymbolsAndEvalMatrix<Lhs>::isTransposed;
	static constexpr bool isLhsConjugated = StripSymbolsAndEvalMatrix<Lhs>::isConjugated;
	static constexpr bool isRhsTransposed = StripSymbolsAndEvalMatrix<Rhs>::isTransposed;
	static constexpr bool isRhsConjugated = StripSymbolsAndEvalMatrix<Rhs>::isConjugated;
	
	static constexpr unsigned int LhsNRows = StripSymbolsAndEvalMatrix<Lhs>::nRows;
	static constexpr unsigned int LhsNCols = StripSymbolsAndEvalMatrix<Lhs>::nCols;
	
	static constexpr unsigned int RhsNRows = StripSymbolsAndEvalMatrix<Rhs>::nRows;
	static constexpr unsigned int RhsNCols = StripSymbolsAndEvalMatrix<Rhs>::nCols;

	std::conditional_t<Lhs::isLeaf, const Lhs&, Lhs> m_lhs;
	std::conditional_t<Rhs::isLeaf, const Rhs&, Rhs> m_rhs;
};

template<typename Expr>        struct IsMatrixProduct                           : std::false_type {};
template<class Lhs, class Rhs> struct IsMatrixProduct< MatrixProduct<Lhs,Rhs> > : std::true_type  {};

template<class Lhs, class Rhs> MatrixProduct<Lhs,Rhs> operator*(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) { return MatrixProduct<Lhs,Rhs>(lhs, rhs); }

template<class Lhs, class Rhs> requires(Lhs::isRowVector and Rhs::isRowVector)
MatrixProduct< Lhs,MatrixTransposed<Rhs> > outer(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) { return MatrixProduct< Lhs,MatrixTransposed<Rhs> >(lhs, MatrixTransposed<Rhs>(rhs)); }

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_MATRIX_PRODUCT_HPP
