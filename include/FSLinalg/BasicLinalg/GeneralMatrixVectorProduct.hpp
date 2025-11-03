#ifndef FSLINALG_GENERAL_MATRIX_VECTOR_PRODUCT_HPP
#define FSLINALG_GENERAL_MATRIX_VECTOR_PRODUCT_HPP

#include <FSLinalg/Scalar.hpp>
#include <FSLinalg/Vector.hpp>
#include <FSLinalg/Matrix.hpp>

namespace FSLinalg
{
namespace BasicLinalg
{

using Size = unsigned int;

template<bool transposeA, bool conjugateA, Size nRowsA, Size nColsA, bool conjugateX, bool incrDst>
struct GeneralMatrixVectorProduct
{
	static constexpr Size nRowsX = transposeA ? nRowsA : nColsA;
	static constexpr Size nRowsY = transposeA ? nColsA : nRowsA;
	
	template<Scalar_concept ScalarAlpha, Scalar_concept ScalarA, Scalar_concept ScalarX, Scalar_concept ScalarY>
	static void run(const ScalarAlpha& alpha, const Matrix<ScalarA,nRowsA,nColsA>& A, const Vector<ScalarX,nRowsX>& x, Vector<ScalarY,nRowsY>& y)
	{			
		constexpr Size A_iStride = (not transposeA) ? nColsA : 1;
		constexpr Size A_jStride = (not transposeA) ?      1 : nColsA;
		
		if constexpr (not incrDst)
		{
			for (Size i=0; i!=nRowsY; ++i) { y[i] = 0; }
		}
		for (Size i=0; i!=nRowsY; ++i)
		{
			for (Size j=0; j!=nRowsX; ++j)
			{
				const ScalarX opXj  = (conjugateX) ? conj(x[j]) : x[j];
				const ScalarA opAij = (conjugateA) ? conj(A[i*A_iStride + j*A_jStride]) : A[i*A_iStride + j*A_jStride];
				y[i] += alpha*opAij*opXj;
			}
		}
	}
	
	template<Scalar_concept ScalarAlpha, Scalar_concept ScalarA, Scalar_concept ScalarY>
	static void run(const ScalarAlpha& alpha, const Matrix<ScalarA,nRowsA,nColsA>& A, const UnitVector<nRowsX>& x, Vector<ScalarY,nRowsY>& y)
	{			
		constexpr Size A_iStride = (not transposeA) ? nColsA : 1;
		constexpr Size A_jStride = (not transposeA) ?      1 : nColsA;
		
		if constexpr (not incrDst)
		{
			for (Size i=0; i!=nRowsY; ++i) { y[i] = 0; }
		}
		for (Size i=0; i!=nRowsY; ++i)
		{
			const Size j = x.getId();
			
			const ScalarA opAij = (conjugateA) ? conj(A[i*A_iStride + j*A_jStride]) : A[i*A_iStride + j*A_jStride];
			y[i] += alpha*opAij;
		}
	}
	
	template<Scalar_concept ScalarAlpha, Scalar_concept ScalarX, Scalar_concept ScalarY>
	static void run(const ScalarAlpha& alpha, const UnitMatrix<nRowsA,nColsA>& A, const Vector<ScalarX,nRowsX>& x, Vector<ScalarY,nRowsY>& y)
	{			
		if constexpr (not incrDst)
		{
			for (Size i=0; i!=nRowsY; ++i) { y[i] = 0; }
		}
		
		const Size i = A.getId().i;
		const Size j = A.getId().j;
		
		const ScalarX opXj = (conjugateX) ? conj(x[j]) : x[j];
		y[i] += alpha*opXj;
	}
	
	template<Scalar_concept ScalarAlpha, Scalar_concept ScalarY>
	static void run(const ScalarAlpha& alpha, const UnitMatrix<nRowsA,nColsA>& A, const UnitVector<nRowsX>& x, Vector<ScalarY,nRowsY>& y)
	{			
		if constexpr (not incrDst)
		{
			for (Size i=0; i!=nRowsY; ++i) { y[i] = 0; }
		}
		
		const Size i = A.getId().i;
		const Size j = A.getId().j;
			
		y[i] += alpha*(j == x.getId());
	}
};
	
} // namespace BasicLinalg
} // namespace FSLinalg

#endif // FSLINALG_GENERAL_MATRIX_VECTOR_PRODUCT_HPP
