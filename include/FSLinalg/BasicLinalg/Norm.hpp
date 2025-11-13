#ifndef FSLINALG_NORM_HPP
#define FSLINALG_NORM_HPP

#include <FSLinalg/Vector/VectorBase.hpp>
#include <FSLinalg/Matrix/MatrixBase.hpp>

namespace FSLinalg
{

template<typename Expr> typename Expr::RealScalar squaredNorm(const MatrixBase<Expr>& base_expr);
template<typename Expr> typename Expr::RealScalar norm(const MatrixBase<Expr>& base_expr) { return std::sqrt(squaredNorm(base_expr)); }

} // namespace FSLinalg

#include <FSLinalg/BasicLinalg/Norm.hpp>

#endif // FSLINALG_NORM_HPP
