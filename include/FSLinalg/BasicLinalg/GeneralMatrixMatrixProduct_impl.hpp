#ifndef FSLINALG_BASIC_LINALG_GENERAL_MATRIX_MATRIX_PRODUCT_IMPL_HPP
#define FSLINALG_BASIC_LINALG_GENERAL_MATRIX_MATRIX_PRODUCT_IMPL_HPP

#include <FSLinalg/BasicLinalg/GeneralMatrixMatrixProduct.hpp>
#include <FSLinalg/BasicLinalg/TripleProduct.hpp>
#include <FSLinalg/BasicLinalg/Product.hpp>

namespace FSLinalg
{
namespace BasicLinalg
{

template<bool transposeA, bool conjugateA, unsigned int nRowsA, unsigned int nColsA, bool transposeB, bool conjugateB, unsigned int nRowsB, unsigned int nColsB, bool incrDst>
template<Scalar_concept ScalarAlpha, Scalar_concept ScalarA, Scalar_concept ScalarB, Scalar_concept ScalarY>
void GeneralMatrixMatrixProduct<transposeA,conjugateA,nRowsA,nColsA,transposeB,conjugateB,nRowsB,nColsB,incrDst>::run(
	const ScalarAlpha&                   alpha, 
	const Matrix<ScalarA,nRowsA,nColsA>& A, 
	const Matrix<ScalarB,nRowsB,nColsB>& B, 
	      Matrix<ScalarY,nRowsY,nColsY>& Y)
{
	using ScaledA = decltype(std::declval<ScalarAlpha>()*std::declval<ScalarA>());
	
	constexpr Size A_iStride = (not transposeA) ? nColsA : 1;
	constexpr Size A_kStride = (not transposeA) ?      1 : nColsA;
	constexpr Size B_kStride = (not transposeB) ? nColsB : 1;
	constexpr Size B_jStride = (not transposeB) ?      1 : nColsB;
	
	constexpr Product<false, conjugateA> prodA;
	constexpr Product<false, conjugateB> prodB;
	
	if constexpr (not incrDst) { Y.setZero(); }
	
	for (Size i=0; i!=nRowsY; ++i)
	{
		for (Size k=0; k !=nColsOpA; ++k)
		{
			ScaledA alphaA = prodA(alpha, A[i*A_iStride + k*A_kStride]);
			for (Size j=0; j!=nColsY; ++j)
			{
				Y(i,j) += prodB(alphaA, B[k*B_kStride + j*B_jStride]);
			}
		}
	}
}

template<bool transposeA, bool conjugateA, unsigned int nRowsA, unsigned int nColsA, bool transposeB, bool conjugateB, unsigned int nRowsB, unsigned int nColsB, bool incrDst>
template<Scalar_concept ScalarAlpha, Scalar_concept ScalarA, Scalar_concept ScalarY>
void GeneralMatrixMatrixProduct<transposeA,conjugateA,nRowsA,nColsA,transposeB,conjugateB,nRowsB,nColsB,incrDst>::run(
	const ScalarAlpha&                   alpha, 
	const Matrix<ScalarA,nRowsA,nColsA>& A, 
	const UnitMatrix<nRowsB,nColsB>&     B, 
	      Matrix<ScalarY,nRowsY,nColsY>& Y)
{
	constexpr Size A_iStride = (not transposeA) ? nColsA : 1;
	constexpr Size A_kStride = (not transposeA) ?      1 : nColsA;
	
	constexpr Product<false, conjugateA> prod;
	
	if constexpr (not incrDst) { Y.setZero(); }
	
	const Size k = (not transposeB) ? B.getId().i : B.getId().j;
	const Size j = (not transposeB) ? B.getId().j : B.getId().i;
	
	for (Size i=0; i!=nRowsY; ++i)
	{
		Y(i,j) += prod(alpha, A[i*A_iStride + k*A_kStride]);
	}
}

template<bool transposeA, bool conjugateA, unsigned int nRowsA, unsigned int nColsA, bool transposeB, bool conjugateB, unsigned int nRowsB, unsigned int nColsB, bool incrDst>
template<Scalar_concept ScalarAlpha, Scalar_concept ScalarB, Scalar_concept ScalarY>
void GeneralMatrixMatrixProduct<transposeA,conjugateA,nRowsA,nColsA,transposeB,conjugateB,nRowsB,nColsB,incrDst>::run(
	const ScalarAlpha&                   alpha, 
	const UnitMatrix<nRowsA,nColsA>&     A, 
	const Matrix<ScalarB,nRowsB,nColsB>& B, 
	      Matrix<ScalarY,nRowsY,nColsY>& Y)
{
	constexpr Size B_kStride = (not transposeB) ? nColsB : 1;
	constexpr Size B_jStride = (not transposeB) ?      1 : nColsB;
	
	constexpr Product<false, conjugateB> prod;
	
	if constexpr (not incrDst) { Y.setZero(); }
	
	const Size i = (not transposeA) ? A.getId().i : A.getId().j;
	const Size k = (not transposeA) ? A.getId().j : A.getId().i;
	
	
	for (Size j=0; j!=nColsY; ++j)
	{
		Y(i,j) += prod(alpha, B[k*B_kStride + j*B_jStride]);
	}
}

template<bool transposeA, bool conjugateA, unsigned int nRowsA, unsigned int nColsA, bool transposeB, bool conjugateB, unsigned int nRowsB, unsigned int nColsB, bool incrDst>
template<Scalar_concept ScalarAlpha, Scalar_concept ScalarY>
void GeneralMatrixMatrixProduct<transposeA,conjugateA,nRowsA,nColsA,transposeB,conjugateB,nRowsB,nColsB,incrDst>::run(
	const ScalarAlpha&                   alpha, 
	const UnitMatrix<nRowsA,nColsA>&     A, 
	const UnitMatrix<nRowsB,nColsB>&     B, 
	      Matrix<ScalarY,nRowsY,nColsY>& Y)
{	
	if constexpr (not incrDst) { Y.setZero(); }
	
	const Size i  = (not transposeA) ? A.getId().i : A.getId().j;
	const Size j  = (not transposeB) ? B.getId().j : B.getId().i;
	const Size k1 = (not transposeA) ? A.getId().j : A.getId().i;
	const Size k2 = (not transposeB) ? B.getId().i : B.getId().j;
	
	Y(i,j) += alpha*(k1 == k2);
}

} // namespace BasicLinalg
} // namespace FSLinalg

#endif // FSLINALG_BASIC_LINALG_GENERAL_MATRIX_MATRIX_PRODUCT_IMPL_HPP
