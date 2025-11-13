#ifndef FSLINALG_MATRIX_BASE_HPP
#define FSLINALG_MATRIX_BASE_HPP

#include <FSLinalg/CRTPBase.hpp>
#include <FSLinalg/Scalar.hpp>
#include <FSLinalg/misc/Logical.hpp>
#include <type_traits>

namespace FSLinalg
{

template<class Derived> struct MatrixTraits;

template<class Derived>
class MatrixBase : public CRTPBase<Derived>
{
public:
	using Base             = CRTPBase<Derived>;	
	using DerivedTraits    = MatrixTraits<Derived>;
	using Scalar           = typename DerivedTraits::Scalar;
	using RealScalar       = typename NumTraits<Scalar>::Real;
	using ReturnType       = std::conditional_t<DerivedTraits::hasWriteRandomAccess,       Scalar&, Scalar>;
	using const_ReturnType = std::conditional_t<DerivedTraits::hasWriteRandomAccess, const Scalar&, Scalar>;
	using Size             = typename DerivedTraits::Size;
	
	template<class Dst> 
	struct IsConvertibleTo : std::bool_constant<
		std::is_base_of<MatrixBase<Dst>, Dst>::value 
		and Dst::hasWriteRandomAccess
		and std::is_convertible<Scalar, typename Dst::Scalar>::value> {};
		
	template<class Src> 
	struct IsConstructibleFrom : std::bool_constant<
		std::is_base_of<MatrixBase<Src>, Src>::value 
		and DerivedTraits::hasWriteRandomAccess
		and std::is_convertible<typename Src::Scalar, Scalar>::value> {};
		
	static constexpr bool hasReadRandomAccess  = DerivedTraits::hasReadRandomAccess;
	static constexpr bool hasWriteRandomAccess = DerivedTraits::hasWriteRandomAccess;
	static constexpr bool hasFlatRandomAccess  = DerivedTraits::hasFlatRandomAccess;
	static constexpr bool causesAliasingIssues = DerivedTraits::causesAliasingIssues;
	static constexpr bool isLeaf               = DerivedTraits::isLeaf;
	static constexpr Size nRows                = DerivedTraits::nRows;
	static constexpr Size nCols                = DerivedTraits::nCols;
	static constexpr Size size                 = nRows*nCols;
	static constexpr bool isRowVector          = (nCols == 1);
	static constexpr bool isColVector          = (nRows == 1);
	
	static_assert(implies(isRowVector and hasReadRandomAccess, hasFlatRandomAccess), "Row vectors with a random access operator must have a flat random access operator");
	static_assert(implies(isColVector and hasReadRandomAccess, hasFlatRandomAccess), "Col vectors with a random access operator must have a flat random access operator");
	
	constexpr Size getRows() const { return nRows; }
	constexpr Size getCols() const { return nCols; }
	constexpr Size getSize() const { return size;  }
	
	const_ReturnType operator()(const Size& i, const Size j) const requires(hasReadRandomAccess)  { return Base::derived().getImpl(i, j); }
	      ReturnType operator()(const Size& i, const Size j)       requires(hasWriteRandomAccess) { return Base::derived().getImpl(i, j); }
	      
	const_ReturnType operator[](const Size& i) const requires(hasReadRandomAccess  and hasFlatRandomAccess) { return Base::derived().getImpl(i); }
	      ReturnType operator[](const Size& i)       requires(hasWriteRandomAccess and hasFlatRandomAccess) { return Base::derived().getImpl(i); }

	template<class Dst> bool isAliasedTo(const MatrixBase<Dst>& other) const { return Base::derived().isAliasedToImpl(other); }

	template<typename Alpha, class Dst, bool checkAliasing>
	void assignTo(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void increment(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
	
	template<typename Alpha, class Dst, bool checkAliasing>
	void decrement(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value);
};

template<typename Expr> struct IsMatrix : std::bool_constant< std::is_base_of<MatrixBase<Expr>, Expr>::value > {};

template<typename Expr> concept Matrix_concept         = IsMatrix<Expr>::value;
template<typename Expr> concept ReadableMatrix_concept = IsMatrix<Expr>::value and Expr::hasReadRandomAccess;
template<typename Expr> concept WritableMatrix_concept = IsMatrix<Expr>::value and Expr::hasWriteRandomAccess;

#define FSLINALG_DEFINE_MATRIX \
	using Base             = MatrixBase<Self>; \
	using Scalar           = typename Base::Scalar; \
	using RealScalar       = typename Base::RealScalar; \
	using ReturnType       = typename Base::ReturnType; \
	using const_ReturnType = typename Base::const_ReturnType; \
	using Size             = typename Base::Size; \
	\
	template<class Dst> using IsConvertibleTo     = typename Base::template IsConvertibleTo<Dst>; \
	template<class Src> using IsConstructibleFrom = typename Base::template IsConstructibleFrom<Src>; \
	\
	using Base::hasReadRandomAccess;\
	using Base::hasWriteRandomAccess;\
	using Base::hasFlatRandomAccess;\
	using Base::causesAliasingIssues;\
	using Base::isLeaf;\
	using Base::nRows;\
	using Base::nCols;\
    using Base::size;\
    using Base::isRowVector;\
    using Base::isColVector;\
    \

template<typename Lhs, typename Rhs> requires(Lhs::hasReadRandomAccess and Rhs::hasReadRandomAccess) bool operator==(const FSLinalg::MatrixBase<Lhs>& lhs, const FSLinalg::MatrixBase<Rhs>& rhs);
template<typename Lhs, typename Rhs> requires(Lhs::hasReadRandomAccess and Rhs::hasReadRandomAccess) bool operator!=(const FSLinalg::MatrixBase<Lhs>& lhs, const FSLinalg::MatrixBase<Rhs>& rhs) { return not(lhs == rhs); }
                
} // namespace FSLinalg

#include <FSLinalg/Matrix/Formater.hpp>

#endif // FSLINALG_MATRIX_BASE_HPP
