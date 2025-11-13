#ifndef FSLINALG_MATRIX_BASE_IMPL_HPP
#define FSLINALG_MATRIX_BASE_IMPL_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>
#include <FSLinalg/Matrix/Matrix.hpp>

namespace FSLinalg
{

template<class Derived> template<typename Alpha, class Dst, bool checkAliasing>
void MatrixBase<Derived>::assignTo(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::nRows == nRows, "Matrix sizes must match");
	static_assert(Dst::nCols == nCols, "Matrix sizes must match");
	
	using DstScalar = typename Dst::Scalar;
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			Matrix<Scalar, nRows, nCols> tmp(*this);
			
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] = alpha*tmp[i]; }
			}
			else
			{
				for (Size i=0; i!=getRows(); ++i) { for (Size j=0; j!=getCols(); ++j) { dst(i,j) = DstScalar(alpha*tmp(i,j)); }}
			}
		}
		else
		{
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] = alpha*Base::derived().getImpl(i); }
			}
			else
			{
				for (Size i=0; i!=getRows(); ++i) { for (Size j=0; j!=getCols(); ++j) { dst(i,j) = DstScalar(alpha*Base::derived().getImpl(i,j)); }}
			}
		}
	}
	else
	{
		Base::derived().assignToImpl(alpha, dst, bc);
	}
}

template<class Derived> template<typename Alpha, class Dst, bool checkAliasing>
void MatrixBase<Derived>::increment(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::nRows == nRows, "Matrix sizes must match");
	static_assert(Dst::nCols == nCols, "Matrix sizes must match");
	
	using DstScalar = typename Dst::Scalar;
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			Matrix<Scalar, nRows, nCols> tmp(*this);
			
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] += alpha*tmp[i]; }
			}
			else
			{
				for (Size i=0; i!=getRows(); ++i) { for (Size j=0; j!=getCols(); ++j) { dst(i,j) += DstScalar(alpha*tmp(i,j)); }}
			}
		}
		else
		{
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] += alpha*Base::derived().getImpl(i); }
			}
			else
			{
				for (Size i=0; i!=getRows(); ++i) { for (Size j=0; j!=getCols(); ++j) { dst(i,j) += DstScalar(alpha*Base::derived().getImpl(i,j)); }}
			}
		}
	}
	else
	{
		Base::derived().incrementImpl(alpha, dst, bc);
	}
}

template<class Derived> template<typename Alpha, class Dst, bool checkAliasing>
void MatrixBase<Derived>::decrement(const Alpha& alpha, MatrixBase<Dst>& dst, std::bool_constant<checkAliasing> bc) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::nRows == nRows, "Matrix sizes must match");
	static_assert(Dst::nCols == nCols, "Matrix sizes must match");
	
	using DstScalar = typename Dst::Scalar;
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			Matrix<Scalar, nRows, nCols> tmp(*this);
			
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] -= alpha*tmp[i]; }
			}
			else
			{
				for (Size i=0; i!=getRows(); ++i) { for (Size j=0; j!=getCols(); ++j) { dst(i,j) -= DstScalar(alpha*tmp(i,j)); }}
			}
		}
		else
		{
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] -= alpha*Base::derived().getImpl(i); }
			}
			else
			{
				for (Size i=0; i!=getRows(); ++i) { for (Size j=0; j!=getCols(); ++j) { dst(i,j) -= DstScalar(alpha*Base::derived().getImpl(i,j)); }}
			}
		}
	}
	else
	{
		Base::derived().decrementImpl(alpha, dst, bc);
	}
}

template<typename Lhs, typename Rhs> requires(Lhs::hasReadRandomAccess and Rhs::hasReadRandomAccess) 
bool operator==(const FSLinalg::MatrixBase<Lhs>& lhs, const FSLinalg::MatrixBase<Rhs>& rhs)
{
	static_assert(Lhs::nRows == Rhs::nRows);
	static_assert(Lhs::nCols == Rhs::nCols);
	
	using Size = std::common_type_t<typename Lhs::Size, typename Rhs::Size>;
	
	bool areEquals = true;
	if constexpr (Lhs::hasFlatRandomAccess and Rhs::hasFlatRandomAccess)
	{
		for (Size i=0; i!=Lhs::size; ++i)
		{
			areEquals = areEquals and (lhs[i] == rhs[i]);
		}

	}
	else
	{
		for (Size i=0; i!=Lhs::nRows; ++i)
		{
			for (Size j=0; j!=Lhs::nCols; ++j)
			{
				areEquals = areEquals and (lhs(i,j) == rhs(i,j));
			}
		}
	}
	return areEquals;
}

} // namespace FSLinalg

#endif // FSLINALG_MATRIX_BASE_IMPL_HPP
