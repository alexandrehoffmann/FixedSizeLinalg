#ifndef FSLINALG_TENSOR_FORMATER_HPP
#define FSLINALG_TENSOR_FORMATER_HPP

#include <FSLinalg/Tensor/TensorBase.hpp>
#include <FSLinalg/misc/NestedLoop.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

template<FSLinalg::Tensor_concept Expr>
class fmt::formatter< Expr > : public fmt::formatter< typename FSLinalg::TensorBase<Expr>::Scalar >
{
public:
	static_assert(Expr::hasReadRandomAccess, "Tensor must have a read random access Iterator");
	
	constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }
	
	template <typename Context>
	auto format (const FSLinalg::TensorBase<Expr>& T, Context& ctx) const 
	{
		using Size  = typename FSLinalg::TensorBase<Expr>::Size;
		using Shape = typename FSLinalg::TensorBase<Expr>::Shape;
				
		auto out = ctx.out();
		
		if constexpr (Expr::rank >= 2)
		{
			const Size rowsM1 = Size(T.shape[Expr::rank-2]-1);
			const Size colsM1 = Size(T.shape[Expr::rank-1]-1);
			
			Shape index;
			FSLinalg::misc::NestedLoop<0,Expr::rank-2>::run(Expr::shape, index, [&](const Shape& /* index */) -> void
			{
				if constexpr (Expr::rank > 2) { out = format_to(out, "[{},:,:]\n", fmt::join(std::cbegin(index), std::cend(index)-2, ",")); }
				out = format_to(out, "[");
				for (index[Expr::rank-2]=0; index[Expr::rank-2]!=rowsM1; ++index[Expr::rank-2])
				{
					out = format_to(out, "[");
					for (index[Expr::rank-1]=0; index[Expr::rank-1]!=colsM1; ++index[Expr::rank-1])
					{
						out = format_to(out, "{}, ", T(index));
					}
					index[Expr::rank-1] = colsM1;
					out = format_to(out, "{}]\n", T(index));
				}
				{
					index[Expr::rank-2] = rowsM1;
					out = format_to(out, "[");
					for (index[Expr::rank-1]=0; index[Expr::rank-1]!=colsM1; ++index[Expr::rank-1])
					{
						out = format_to(out, "{}, ", T(index));
					}
					index[Expr::rank-1] = colsM1;
					out = format_to(out, "{}]]\n", T(index));
				}
			});
			return out;
		}
		else
		{
			out = format_to(out, "[");
			for (Size i=0; i!=Size(Expr::size-1); ++i)
			{
				out = format_to(out, "{}, ", T[i]);
			}
			return format_to(out, "{}]", T[Size(Expr::size-1)]);
		}
	}
};

#endif // FSLINALG_TENSOR_FORMATER_HPP
