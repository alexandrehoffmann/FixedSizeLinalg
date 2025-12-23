#ifndef FSLINALG_MATRIX_HPP
#define FSLINALG_MATRIX_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>

#include <array>

namespace FSLinalg
{

template<typename T, unsigned int Nrows, unsigned Ncols> class Matrix;

template<typename T, unsigned int Nrows, unsigned Ncols>
struct MatrixTraits< Matrix<T, Nrows, Ncols> >
{	
	using Scalar = T;
	using Size   = unsigned int;
	
	static constexpr bool hasReadRandomAccess  = true;
	static constexpr bool hasWriteRandomAccess = true;
	static constexpr bool hasFlatRandomAccess  = true;
	static constexpr bool causesAliasingIssues = false;
	static constexpr bool isLeaf               = true;
	
	static constexpr Size nRows = Nrows;
	static constexpr Size nCols = Ncols;
};

template<typename T, unsigned int Nrows, unsigned Ncols> 
class Matrix : public MatrixBase< Matrix<T, Nrows, Ncols> >
{
public:
	using Self = Matrix<T, Nrows, Ncols>;
	FSLINALG_DEFINE_MATRIX
	
	template<class Dst>
	struct CanBeAlisaedTo : BIC::Fixed<bool,  
		    IsMatrix<Dst>::value 
		and Base::nRows == Dst::nRows
		and Base::nCols == Dst::nCols
		and std::is_same<Scalar, typename Dst::Scalar>::value > {};
	
	static constexpr bool isScalarComplex = IsComplexScalar<Scalar>::value;
	static constexpr bool isVector = isRowVector or isColVector;
	
	Matrix(const RealScalar& value = RealScalar(0))  requires(isScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] = value; } }
	Matrix(std::initializer_list< std::initializer_list<RealScalar> > values) requires(isScalarComplex);
	Matrix(std::initializer_list<RealScalar> values) requires(isScalarComplex and isVector) { std::copy(std::cbegin(values), std::cend(values), std::begin(m_data)); }
	
	Matrix(const Scalar& value = Scalar(0)) { m_data.fill(value); }
	Matrix(std::initializer_list< std::initializer_list<Scalar> > values);
	Matrix(std::initializer_list<Scalar> values) requires(isVector) { std::copy(std::cbegin(values), std::cend(values), std::begin(m_data)); }
	
	Matrix(const Matrix& other) : m_data(other.m_data) {}
	
	template<class Expr> Matrix(const MatrixBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.assignTo(BIC::fixed<bool, false>, BIC::fixed<RealScalar, RealScalar(1)>, *this); }
	
	template<class Expr> Matrix& operator= (const MatrixBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.assignTo  (BIC::fixed<bool, true>, BIC::fixed<RealScalar, RealScalar(1)>, *this); return *this; }
	template<class Expr> Matrix& operator+=(const MatrixBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.increment (BIC::fixed<bool, true>, BIC::fixed<RealScalar, RealScalar(1)>, *this); return *this; }
	template<class Expr> Matrix& operator-=(const MatrixBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.decrement (BIC::fixed<bool, true>, BIC::fixed<RealScalar, RealScalar(1)>, *this); return *this; }
	
	Matrix& operator*=(const RealScalar& alpha) requires(isScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] *= alpha; } return *this; }
	Matrix& operator/=(const RealScalar& alpha) requires(isScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] /= alpha; } return *this; }
	
	Matrix& operator*=(const Scalar& alpha) { for (Size i=0; i!=size; ++i) { m_data[i] *= alpha; } return *this; }
	Matrix& operator/=(const Scalar& alpha) { for (Size i=0; i!=size; ++i) { m_data[i] /= alpha; } return *this; }
	
	void setZero() { m_data.fill(Scalar(0)); }
	
	const_ReturnType getImpl(const Size i) const { return m_data[i]; }
	      ReturnType getImpl(const Size i)       { return m_data[i]; }
	      
	const_ReturnType getImpl(const Size i, const Size j) const { return m_data[i*nCols + j]; }
	      ReturnType getImpl(const Size i, const Size j)       { return m_data[i*nCols + j]; }
	      
	template<class Dst>           bool isAliasedToImpl(const MatrixBase<Dst>& dst) const requires(    CanBeAlisaedTo<Dst>::value) { return std::addressof(dst.derived()) == this; }
	template<class Dst> constexpr bool isAliasedToImpl(const MatrixBase<Dst>&    ) const requires(not CanBeAlisaedTo<Dst>::value) { return false; }
	
	static Matrix zero()   { return Matrix(RealScalar(0)); }
	static Matrix ones()   { return Matrix(RealScalar(1)); }
	
	static Matrix random(const RealScalar& lb = RealScalar(-1), const RealScalar& ub = RealScalar(1));
private:
	std::array<Scalar, size> m_data;
};

template<unsigned int Nrows, unsigned Ncols> using RealMatrix = Matrix<double, Nrows, Ncols>;
template<unsigned int Nrows, unsigned Ncols> using CpxMatrix  = Matrix<std::complex<double>, Nrows, Ncols>;

template<typename T, unsigned int Nrows> using RowVector     = Matrix<T, Nrows, 1>;
template<unsigned int Nrows>             using RealRowVector = Matrix<double, Nrows, 1>;
template<unsigned int Nrows>             using CpxRowVector  = Matrix<std::complex<double>, Nrows, 1>;

template<typename T, unsigned int Ncols> using ColVector     = Matrix<T, 1, Ncols>;
template<unsigned int Ncols>             using RealColVector = Matrix<double, 1, Ncols>;
template<unsigned int Ncols>             using CpxColVector  = Matrix<std::complex<double>, 1, Ncols>;

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_HPP
