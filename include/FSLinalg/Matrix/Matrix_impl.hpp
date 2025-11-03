#ifndef FSLINALG_MATRIX_IMPL_HPP
#define FSLINALG_MATRIX_IMPL_HPP

#include <FSLinalg/Matrix/Matrix.hpp>

#include <cassert>

namespace FSLinalg
{

template<typename T, unsigned int Nrows, unsigned Ncols>
Matrix<T,Nrows,Ncols>::Matrix(std::initializer_list< std::initializer_list<RealScalar> > values) requires(IsScalarComplex)
{
	using Iterator = typename std::array<Scalar, size>::iterator;
	
	assert(values.size() == Nrows);
	
	Iterator it_data = std::begin(m_data);
	for (const std::initializer_list<RealScalar>& row_values : values)
	{
		assert(row_values.size() == nCols);
		std::copy(std::cbegin(row_values), std::cend(row_values), it_data);
		it_data += nCols;
	}
}

template<typename T, unsigned int Nrows, unsigned Ncols>
Matrix<T,Nrows,Ncols>::Matrix(std::initializer_list< std::initializer_list<Scalar> > values)
{
	using Iterator = typename std::array<Scalar, size>::iterator;
	
	assert(values.size() == Nrows);
	
	Iterator it_data = std::begin(m_data);
	for (const std::initializer_list<Scalar>& row_values : values)
	{
		assert(row_values.size() == nCols);
		std::copy(std::cbegin(row_values), std::cend(row_values), it_data);
		it_data += nCols;
	}
}

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_IMPL_HPP
