#include <FSLinalg/Tensor.hpp>

int main()
{	
	FSLinalg::RealTensor<3,3,3> a = FSLinalg::RealTensor<3,3,3>::random();
	FSLinalg::RealTensor<3,3,3> b = FSLinalg::RealTensor<3,3,3>::random();
	FSLinalg::RealTensor<3,3,3> c = FSLinalg::RealTensor<3,3,3>::random();
	FSLinalg::RealTensor<3,3,3> d = FSLinalg::RealTensor<3,3,3>::random();
	
	const auto expr1 = (a + b)*(c + d);
	const auto expr2 = (a - b)*(c - d);
	
	fmt::print("a =\n{}\n", a);
	fmt::print("b =\n{}\n", b);
	fmt::print("c =\n{}\n", c);
	fmt::print("d =\n{}\n", d);
	fmt::print("\n");
	fmt::print("expr1 =\n{}\n", expr1);
	fmt::print("expr2 =\n{}\n", expr2);
	fmt::print("\n");
	fmt::print("expr1 < expr2 =\n{}\n", expr1 < expr2);
	fmt::print("any(expr2) = {}\n", FSLinalg::any(expr1 < expr2));
	fmt::print("all(expr2) = {}\n", FSLinalg::all(expr1 < expr2));
	
	return EXIT_SUCCESS;
}
