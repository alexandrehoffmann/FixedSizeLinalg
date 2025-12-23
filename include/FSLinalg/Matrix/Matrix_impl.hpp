#ifndef FSLINALG_MATRIX_IMPL_HPP
#define FSLINALG_MATRIX_IMPL_HPP

#include <FSLinalg/Matrix/Matrix.hpp>

#include <cassert>
#include <random>

namespace FSLinalg
{

template<typename T, unsigned int Nrows, unsigned Ncols>
Matrix<T,Nrows,Ncols>::Matrix(std::initializer_list< std::initializer_list<RealScalar> > values) requires(isScalarComplex)
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

template<typename T, unsigned int Nrows, unsigned Ncols>
auto Matrix<T,Nrows,Ncols>::random(const RealScalar& lb, const RealScalar& ub) -> Matrix
{
	std::random_device rd; 
	std::mt19937 gen(rd()); 
	std::uniform_real_distribution<RealScalar> dist(lb, ub);
	
	Matrix ret;
	
	for (Size i=0; i!=size; ++i) { ret[i] = dist(gen); }
	
	return ret;
}

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_IMPL_HPP
