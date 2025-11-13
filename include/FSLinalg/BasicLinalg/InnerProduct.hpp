#ifndef FSLINALG_INNER_PRODUCT_HPP
#define FSLINALG_INNER_PRODUCT_HPP

#include <FSLinalg/Vector/VectorBase.hpp>
#include <FSLinalg/Matrix/MatrixBase.hpp>

namespace FSLinalg
{

template<class Lhs, class Rhs>
using InnerProductScalar = decltype(conj(std::declval<typename Lhs::Scalar>()) * std::declval<typename Rhs::Scalar>());

template<class Lhs, class Rhs> InnerProductScalar<Lhs,Rhs> inner(const MatrixBase<Lhs>& base_lhs, const MatrixBase<Rhs>& base_rhs);

	
} // namespace FSLinalg

#include <FSLinalg/BasicLinalg/InnerProduct_impl.hpp>

#endif // FSLINALG_INNER_PRODUCT_HPP
