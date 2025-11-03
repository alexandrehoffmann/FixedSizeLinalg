#include <FSLinalg/Vector.hpp>

int main()
{
	FSLinalg::UnitVector<3> e0(0);
	FSLinalg::RealVector<3> one(1);
	FSLinalg::RealVector<3> a({2, 3, 0});
	
	const auto expr = a + 3.*e0 - one/2.;
	
	fmt::print("{} + 3*{} - {}/2 = {}\n", a, e0, one, expr);
	fmt::print("expr is aliased to a = {}\n", expr.isAliasedTo(a));
	
	const auto expr2 = -2.*FSLinalg::cross(a, one);
	FSLinalg::RealVector<3> c = expr2;
	
	fmt::print("cross({}, {}) = {}\n", a, one, c);
	fmt::print("expr2 is aliased to a = {}\n", expr2.isAliasedTo(a));
	
	a = expr2;
	fmt::print("a = {}\n", a);
	
	return EXIT_SUCCESS;
}
