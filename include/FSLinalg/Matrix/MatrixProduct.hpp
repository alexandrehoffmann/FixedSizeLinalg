#ifndef FSLINALG_MATRIX_PRODUCT_HPP
#define FSLINALG_MATRIX_PRODUCT_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>
#include <FSLinalg/Matrix/MatrixTransposed.hpp>
#include <FSLinalg/Matrix/StripSymbolsAndEvalMatrix.hpp>
#include <FSLinalg/Matrix/MatrixProductAnalyzer.hpp>

namespace FSLinalg
{

namespace detail
{

template<class Expr> struct MatrixProductAnalyzerImpl;

} // namespace detail

template<class Lhs, class Rhs> class MatrixProduct;
template<class Expr>           class KeepBrackets;

template<class Lhs, class Rhs>
struct MatrixTraits< MatrixProduct<Lhs, Rhs> >
{
	static_assert(IsMatrix<Lhs>::value and IsMatrix<Rhs>::value, "Both LHS and RHS must be matrices");
	static_assert(Lhs::nCols == Rhs::nRows, "Matrices size must match");
		
	using Scalar = decltype(std::declval<typename Lhs::Scalar>() * std::declval<typename Rhs::Scalar>());
	using Size   = std::common_type_t<typename Lhs::Size, typename Rhs::Size>;
	
	static constexpr bool hasReadRandomAccess  = (Lhs::isRowVector and Lhs::hasFlatRandomAccess) or (Rhs::isColVector and Rhs::hasFlatRandomAccess);
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
	
	friend struct detail::MatrixProductAnalyzerImpl< Self >;
	friend class KeepBrackets< Self >;
	
	MatrixProduct(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) : m_lhs(lhs.derived()), m_rhs(rhs.derived()) {}
	
	static constexpr bool createTemporaryLhs = StripSymbolsAndEvalMatrix<Lhs>::createsTemporary;
	static constexpr bool createTemporaryRhs = StripSymbolsAndEvalMatrix<Rhs>::createsTemporary;
	
	static constexpr bool isOptimallyBracked();
	
	/**
	 * @brief Only awailable when multiplying row-vector and col-vector
	 * When multiplying a row-vector and a col-vector, we can compute A*B as A(i,0)*B(0,j)
	 */
	const_ReturnType getImpl(const Size i, const Size j) const requires(hasReadRandomAccess)  { return m_lhs.getImpl(i)*m_rhs.getImpl(j); }
	
	template<class Dst> bool isAliasedToImpl(const MatrixBase<Dst>& other) const { return m_lhs.isAliasedToImpl(other) or m_rhs.isAliasedToImpl(other); }

	template<typename Bool, typename Alpha, class Dst>
	void assignToImpl(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { assignToHelper(checkAliasing, alpha, dst, BIC::fixed<bool, false>); }
	
	template<typename Bool, typename Alpha, class Dst>
	void incrementImpl(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { incrementHelper(checkAliasing, alpha, dst, BIC::fixed<bool, false>); }
	
	template<typename Bool, typename Alpha, class Dst>
	void decrementImpl(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value) { decrementHelper(checkAliasing, alpha, dst, BIC::fixed<bool, false>); }
private:
	template<typename Bool, typename Alpha, class Dst, bool keepBracketing>
	void assignToHelper(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst, BIC::Fixed<bool, keepBracketing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Bool, typename Alpha, class Dst, bool keepBracketing>
	void incrementHelper(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst, BIC::Fixed<bool, keepBracketing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Bool, typename Alpha, class Dst, bool keepBracketing>
	void decrementHelper(const Bool checkAliasing, const Alpha& alpha, MatrixBase<Dst>& dst, BIC::Fixed<bool, keepBracketing>) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);

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

template<typename Expr>        struct IsMatrixProduct                           : BIC::Fixed<bool, false> {};
template<class Lhs, class Rhs> struct IsMatrixProduct< MatrixProduct<Lhs,Rhs> > : BIC::Fixed<bool, true>  {};

template<class Lhs, class Rhs> MatrixProduct<Lhs,Rhs> operator*(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) { return MatrixProduct<Lhs,Rhs>(lhs, rhs); }

template<class Lhs, class Rhs> requires(Lhs::isRowVector and Rhs::isRowVector)
MatrixProduct< Lhs,MatrixTransposed<Rhs> > outer(const MatrixBase<Lhs>& lhs, const MatrixBase<Rhs>& rhs) { return MatrixProduct< Lhs,MatrixTransposed<Rhs> >(lhs, MatrixTransposed<Rhs>(rhs)); }

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_PRODUCT_HPP
