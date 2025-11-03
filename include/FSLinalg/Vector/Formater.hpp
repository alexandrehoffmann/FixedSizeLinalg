#ifndef FSLINALG_VECTOR_FORMATER_HPP
#define FSLINALG_VECTOR_FORMATER_HPP

#include <FSLinalg/Vector/VectorBase.hpp>

#include <fmt/core.h>

template<FSLinalg::Vector_concept Expr>
class fmt::formatter< Expr > : public fmt::formatter< typename FSLinalg::VectorBase<Expr>::Scalar >
{
public:
	static_assert(Expr::hasReadRandomAccess, "Vector must have a read random access Iterator");

	constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }

	template <typename Context>
	auto format (const FSLinalg::VectorBase<Expr>& x, Context& ctx) const 
	{
		using Size = typename FSLinalg::VectorBase<Expr>::Size;
	
		ctx.out() = format_to(ctx.out(), "[");
		for (Size i=0; i!=Size(x.getSize()-1); ++i)
		{
		ctx.out() = format_to(ctx.out(), "{}, ", x[i]);
		}
		return format_to(ctx.out(), "{}]", x[Size(x.getSize()-1)]);
	}
};

#endif // FSLINALG_VECTOR_FORMATER_HPP
