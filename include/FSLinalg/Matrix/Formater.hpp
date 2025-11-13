#ifndef FSLINALG_MATRIX_FORMATER_HPP
#define FSLINALG_MATRIX_FORMATER_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>

#include <fmt/core.h>

template<FSLinalg::Matrix_concept Expr>
class fmt::formatter< Expr > : public fmt::formatter< typename FSLinalg::MatrixBase<Expr>::Scalar >
{
public:
	static_assert(Expr::hasReadRandomAccess, "Vector must have a read random access Iterator");
	
	constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }
	
	template <typename Context>
	auto format (const FSLinalg::MatrixBase<Expr>& A, Context& ctx) const 
	{
		using Size = typename FSLinalg::MatrixBase<Expr>::Size;
		
		Size rowsM1 = Size(A.getRows()-1);
		Size colsM1 = Size(A.getCols()-1);
		
		if constexpr (Expr::isRowVector or Expr::isColVector)
		{
			ctx.out() = format_to(ctx.out(), "[");
			for (Size i=0; i!=Size(A.getSize()-1); ++i)
			{
				ctx.out() = format_to(ctx.out(), "{}, ", A[i]);
			}
			if constexpr (Expr::isRowVector)
			{
				return format_to(ctx.out(), "{}]^T", A[Size(A.getSize()-1)]);
			}
			else
			{
				return format_to(ctx.out(), "{}]", A[Size(A.getSize()-1)]);
			}
		}
		else
		{
			ctx.out() = format_to(ctx.out(), "[");
			for (Size i=0; i!=rowsM1; ++i)
			{
				ctx.out() = format_to(ctx.out(), "[");
				for (Size j=0; j!=colsM1; ++j)
				{
					ctx.out() = format_to(ctx.out(), "{}, ", A(i,j));
				}
				ctx.out() = format_to(ctx.out(), "{}]\n", A(i, colsM1));
			}
			ctx.out() = format_to(ctx.out(), "[");
			for (Size j=0; j!=colsM1; ++j)
			{
				ctx.out() = format_to(ctx.out(), "{}, ", A(rowsM1,j));
			}
			return format_to(ctx.out(), "{}]]", A(rowsM1,colsM1));
		}
	}
};

#endif // FSLINALG_MATRIX_FORMATER_HPP
