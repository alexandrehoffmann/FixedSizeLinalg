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
	using Scalar     = T;
	using Size       = unsigned int;
	
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
	
	static constexpr bool IsScalarComplex = IsComplexScalar<Scalar>::value;
	
	Matrix(const RealScalar& value = RealScalar(0))  requires(IsScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] = value; } }
	Matrix(std::initializer_list< std::initializer_list<RealScalar> > values) requires(IsScalarComplex);
	
	Matrix(const Scalar& value = Scalar(0))      { m_data.fill(value); }
	Matrix(std::initializer_list< std::initializer_list<Scalar> > values);
	
	Matrix(const Matrix& other) : m_data(other.m_data) {}
	
	template<class Expr> Matrix(const MatrixBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.assignTo(1., *this, std::false_type{}); }
	
	template<class Expr> Matrix& operator= (const MatrixBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.assignTo  (1., *this,  std::true_type{}); return *this; }
	template<class Expr> Matrix& operator+=(const MatrixBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.increment (1., *this, std::true_type{}); return *this; }
	template<class Expr> Matrix& operator-=(const MatrixBase<Expr>& expr) requires(IsConstructibleFrom<Expr>::value) { expr.decrement (1., *this, std::true_type{}); return *this; }
	
	Matrix& operator*=(const RealScalar& alpha) requires(IsScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] *= alpha; } return *this; }
	Matrix& operator/=(const RealScalar& alpha) requires(IsScalarComplex) { for (Size i=0; i!=size; ++i) { m_data[i] /= alpha; } return *this; }
	
	Matrix& operator*=(const Scalar& alpha) { for (Size i=0; i!=size; ++i) { m_data[i] *= alpha; } return *this; }
	Matrix& operator/=(const Scalar& alpha) { for (Size i=0; i!=size; ++i) { m_data[i] /= alpha; } return *this; }
	
	const_ReturnType getImpl(const Size i) const { return m_data[i]; }
	      ReturnType getImpl(const Size i)       { return m_data[i]; }
	      
	const_ReturnType getImpl(const Size i, const Size j) const { return m_data[i*nCols + j]; }
	      ReturnType getImpl(const Size i, const Size j)       { return m_data[i*nCols + j]; }
	      
	template<class Dst> 
	constexpr bool isAliasedToImpl(const VectorBase<Dst>&) const { return false; }
	
	template<class Dst> 
	bool isAliasedToImpl(const MatrixBase<Dst>& dst) const 
	{ 
		if constexpr (nRows == Dst::nRows and nCols == Dst::nCols) { return std::addressof(dst.derived()) == this; }
		else                                                       { return false; }
	}
private:
	std::array<Scalar, size> m_data;
};

template<unsigned int Nrows, unsigned Ncols> using RealMatrix = Matrix<double, Nrows, Ncols>;
template<unsigned int Nrows, unsigned Ncols> using CpxMatrix  = Matrix<std::complex<double>, Nrows, Ncols>;

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_HPP
