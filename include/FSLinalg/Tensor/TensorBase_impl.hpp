#ifndef FSLINALG_TENSOR_BASE_IMPL_HPP
#define FSLINALG_TENSOR_BASE_IMPL_HPP

#include <FSLinalg/Tensor/TensorBase.hpp>
#include <FSLinalg/Tensor/Tensor.hpp>
#include <FSLinalg/misc/NestedLoop.hpp>

namespace FSLinalg
{

template<class Derived> template<typename Bool, typename Alpha, class Dst>
void TensorBase<Derived>::assignTo(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::shape == shape, "Tensor shapes must match");
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			TensorFromShape<Scalar, shape> tmp(*this);
			
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] = alpha*tmp[i]; }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) = alpha*tmp(index);
				});
			}
		}
		else
		{
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] = alpha*CRTP::derived().getImpl(i); }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) = alpha*CRTP::derived().getImpl(index);
				});
			}
		}
	}
	else
	{
		CRTP::derived().assignToImpl(checkAliasing, alpha, dst);
	}
}

template<class Derived> template<typename Bool, typename Alpha, class Dst>
void TensorBase<Derived>::increment(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::shape == shape, "Tensor shapes must match");
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			TensorFromShape<Scalar, shape> tmp(*this);
			
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] += alpha*tmp[i]; }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) += alpha*tmp(index);
				});
			}
		}
		else
		{
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] += alpha*CRTP::derived().getImpl(i); }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) += alpha*CRTP::derived().getImpl(index);
				});
			}
		}
	}
	else
	{
		CRTP::derived().incrementImpl(checkAliasing, alpha, dst);
	}
}

template<class Derived> template<typename Bool, typename Alpha, class Dst>
void TensorBase<Derived>::decrement(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::shape == shape, "Tensor shapes must match");
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			TensorFromShape<Scalar, shape> tmp(*this);
			
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] -= alpha*tmp[i]; }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) -= alpha*tmp(index);
				});
			}
		}
		else
		{
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] -= alpha*CRTP::derived().getImpl(i); }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) -= alpha*CRTP::derived().getImpl(index);
				});
			}
		}
	}
	else
	{
		CRTP::derived().decrementImpl(checkAliasing, alpha, dst);
	}
}

template<class Derived> template<typename Bool, typename Alpha, class Dst>
void TensorBase<Derived>::multiply(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::shape == shape, "Tensor shapes must match");
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			TensorFromShape<Scalar, shape> tmp(*this);
			
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] *= alpha*tmp[i]; }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) *= alpha*tmp(index);
				});
			}
		}
		else
		{
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] *= alpha*CRTP::derived().getImpl(i); }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) *= alpha*CRTP::derived().getImpl(index);
				});
			}
		}
	}
	else
	{
		CRTP::derived().multiplyImpl(checkAliasing, alpha, dst);
	}
}

template<class Derived> template<typename Bool, typename Alpha, class Dst>
void TensorBase<Derived>::divide(const Bool checkAliasing, const Alpha& alpha, TensorBase<Dst>& dst) const requires(IsConvertibleTo<Dst>::value and IsScalar<Alpha>::value)
{
	static_assert(Dst::shape == shape, "Tensor shapes must match");
	
	if constexpr (hasReadRandomAccess)
	{
		if (checkAliasing and causesAliasingIssues and isAliasedTo(dst))
		{
			TensorFromShape<Scalar, shape> tmp(*this);
			
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] /= (alpha*tmp[i]); }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) /= (alpha*tmp(index));
				});
			}
		}
		else
		{
			if constexpr (hasFlatRandomAccess)
			{
				for (Size i=0; i!=getSize(); ++i) { dst[i] /= (alpha*CRTP::derived().getImpl(i)); }
			}
			else
			{
				misc::nestedLoop(shape, [&](const Shape& index) -> void
				{
					dst(index) /= (alpha*CRTP::derived().getImpl(index));
				});
			}
		}
	}
	else
	{
		CRTP::derived().divideImpl(checkAliasing, alpha, dst);
	}
}

template<typename Lhs, typename Rhs> requires(Lhs::hasReadRandomAccess and Rhs::hasReadRandomAccess) 
bool operator==(const FSLinalg::TensorBase<Lhs>& lhs, const FSLinalg::TensorBase<Rhs>& rhs)
{
	static_assert(Lhs::shape == Rhs::shape, "Tensor shapes must match");
	
	using Size  = std::common_type_t<typename Lhs::Size,  typename Rhs::Size>;
	using Shape = std::common_type_t<typename Lhs::Shape, typename Rhs::Shape>;
	
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
		misc::nestedLoop(Lhs::shape, [&](const Shape& index) -> void
		{
			areEquals = areEquals and (lhs(index) == rhs(index));
		});
	}
	return areEquals;
}

template<class Expr> requires(std::is_same<typename Expr::Scalar, bool>::value and Expr::hasReadRandomAccess) 
bool all(const TensorBase<Expr>& expr)
{
	using Size  = typename Expr::Size;
	using Shape = typename Expr::Shape;
	
	bool allTrue = true;
	if constexpr (Expr::hasFlatRandomAccess)
	{
		for (Size i=0; i!=Expr::size; ++i)
		{
			allTrue = allTrue and expr[i];
		}

	}
	else
	{
		misc::nestedLoop(Expr::shape, [&](const Shape& index) -> void
		{
			allTrue = allTrue and expr(index);
		});
	}
	return allTrue;
}

template<class Expr> requires(std::is_same<typename Expr::Scalar, bool>::value and Expr::hasReadRandomAccess) 
bool any(const TensorBase<Expr>& expr)
{
	using Size  = typename Expr::Size;
	using Shape = typename Expr::Shape;
	
	bool anyTrue = false;
	if constexpr (Expr::hasFlatRandomAccess)
	{
		for (Size i=0; i!=Expr::size; ++i)
		{
			anyTrue = anyTrue or expr[i];
		}

	}
	else
	{
		misc::nestedLoop(Expr::shape, [&](const Shape& index) -> void
		{
			anyTrue = anyTrue or expr(index);
		});
	}
	return anyTrue;
}

} // namespace FSLinalg

#endif // FSLINALG_TENSOR_BASE_IMPL_HPP
