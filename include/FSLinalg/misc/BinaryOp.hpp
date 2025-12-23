#ifndef FSLINALG_BINARY_OPERATORS_HPP
#define FSLINALG_BINARY_OPERATORS_HPP

namespace FSLinalg
{
namespace BinaryOp
{
	
struct Add
{
	template<typename Lhs, typename Rhs>
	constexpr auto operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs + rhs; }
};
	
struct Sub
{
	template<typename Lhs, typename Rhs>
	constexpr auto operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs - rhs; }
};
	
struct Mul
{
	template<typename Lhs, typename Rhs>
	constexpr auto operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs * rhs; }
};
	
struct Div
{
	template<typename Lhs, typename Rhs>
	constexpr auto operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs / rhs; }
};
	
struct Mod
{
	template<typename Lhs, typename Rhs>
	constexpr auto operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs % rhs; }
};
	
struct Equal
{
	template<typename Lhs, typename Rhs>
	constexpr bool operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs == rhs; }
};
	
struct NotEqual
{
	template<typename Lhs, typename Rhs>
	constexpr bool operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs != rhs; }
};
	
struct Greater
{
	template<typename Lhs, typename Rhs>
	constexpr bool operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs > rhs; }
};
	
struct GreaterOrEqual
{
	template<typename Lhs, typename Rhs>
	constexpr bool operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs >= rhs; }
};
	
struct Lower
{
	template<typename Lhs, typename Rhs>
	constexpr bool operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs < rhs; }
};
	
struct LowerOrEqual
{
	template<typename Lhs, typename Rhs>
	constexpr bool operator() (const Lhs& lhs, const Rhs& rhs) const { return lhs <= rhs; }
};

} // namespace BinaryOp
} // namespace FSLinalg

#endif // FSLINALG_BINARY_OPERATORS_HPP
