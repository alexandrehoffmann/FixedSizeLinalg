#ifndef FSLINALG_BASIC_LINALG_GENERAL_MATRIX_MATRIX_PRODUCT_HPP
#define FSLINALG_BASIC_LINALG_GENERAL_MATRIX_MATRIX_PRODUCT_HPP

#include <FSLinalg/Scalar.hpp>
#include <FSLinalg/Matrix.hpp>

namespace FSLinalg
{
namespace BasicLinalg
{

template<bool transposeA, bool conjugateA, unsigned int nRowsA, unsigned int nColsA, bool transposeB, bool conjugateB, unsigned int nRowsB, unsigned int nColsB, bool incrDst>
struct GeneralMatrixMatrixProduct
{
	using Size = unsigned int;
	
	static constexpr Size nRowsOpA = (not transposeA) ? nRowsA : nColsA;
	static constexpr Size nColsOpA = (not transposeA) ? nColsA : nRowsA;
	static constexpr Size nRowsOpB = (not transposeB) ? nRowsB : nColsB;
	static constexpr Size nColsOpB = (not transposeB) ? nColsB : nRowsB;
	static constexpr Size nRowsY   = nRowsOpA;
	static constexpr Size nColsY   = nColsOpB;
	
	static_assert(nColsOpA == nRowsOpB, "Matrices sizes must match");
	
	template<Scalar_concept ScalarAlpha, Scalar_concept ScalarA, Scalar_concept ScalarB, Scalar_concept ScalarY>
	static void run(const ScalarAlpha& alpha, const Matrix<ScalarA,nRowsA,nColsA>& A, const Matrix<ScalarB,nRowsB,nColsB>& B, Matrix<ScalarY,nRowsY,nColsY>& Y);
	
	template<Scalar_concept ScalarAlpha, Scalar_concept ScalarA, Scalar_concept ScalarY>
	static void run(const ScalarAlpha& alpha, const Matrix<ScalarA,nRowsA,nColsA>& A, const UnitMatrix<nRowsB,nColsB>& B, Matrix<ScalarY,nRowsY,nColsY>& Y);
	
	template<Scalar_concept ScalarAlpha, Scalar_concept ScalarB, Scalar_concept ScalarY>
	static void run(const ScalarAlpha& alpha, const UnitMatrix<nRowsA,nColsA>& A, const Matrix<ScalarB,nRowsB,nColsB>& B, Matrix<ScalarY,nRowsY,nColsY>& Y);
	
	template<Scalar_concept ScalarAlpha, Scalar_concept ScalarY>
	static void run(const ScalarAlpha& alpha, const UnitMatrix<nRowsA,nColsA>& A, const UnitMatrix<nRowsB,nColsB>& B, Matrix<ScalarY,nRowsY,nColsY>& Y);
};
	
} // namespace BasicLinalg
} // namespace FSLinalg

#include <FSLinalg/BasicLinalg/GeneralMatrixMatrixProduct_impl.hpp>

#endif // FSLINALG_BASIC_LINALG_GENERAL_MATRIX_MATRIX_PRODUCT_HPP
