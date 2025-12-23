#ifndef FSLINALG_MATRIX_FORMATER_HPP
#define FSLINALG_MATRIX_FORMATER_HPP

#include <FSLinalg/Matrix/MatrixBase.hpp>

#include <fmt/core.h>

template<FSLinalg::Matrix_concept Expr>
class fmt::formatter< Expr > : public fmt::formatter< typename FSLinalg::MatrixBase<Expr>::Scalar >
{
public:
	static_assert(Expr::hasReadRandomAccess, "Matrix must have a read random access Iterator");
	
	constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }
	
	template <typename Context>
	auto format (const FSLinalg::MatrixBase<Expr>& A, Context& ctx) const 
	{
		using Size = typename FSLinalg::MatrixBase<Expr>::Size;
		
		const Size rowsM1 = Size(A.getRows()-1);
		const Size colsM1 = Size(A.getCols()-1);
		
		auto out = ctx.out();
		
		if constexpr (Expr::isRowVector or Expr::isColVector)
		{
			out = format_to(out, "[");
			for (Size i=0; i!=Size(A.getSize()-1); ++i)
			{
				out = format_to(out, "{}, ", A[i]);
			}
			if constexpr (Expr::isRowVector)
			{
				return format_to(out, "{}]^T", A[Size(A.getSize()-1)]);
			}
			else
			{
				return format_to(out, "{}]", A[Size(A.getSize()-1)]);
			}
		}
		else
		{
			out = format_to(out, "[");
			for (Size i=0; i!=rowsM1; ++i)
			{
				out = format_to(out, "[");
				for (Size j=0; j!=colsM1; ++j)
				{
					out = format_to(out, "{}, ", A(i,j));
				}
				out = format_to(out, "{}]\n", A(i, colsM1));
			}
			out = format_to(out, "[");
			for (Size j=0; j!=colsM1; ++j)
			{
				out = format_to(out, "{}, ", A(rowsM1,j));
			}
			return format_to(out, "{}]]", A(rowsM1,colsM1));
		}
	}
};

#endif // FSLINALG_MATRIX_FORMATER_HPP
