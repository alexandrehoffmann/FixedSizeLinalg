#include <FSLinalg/Matrix.hpp>

int main()
{
	FSLinalg::UnitRowVector<4> e0(0);
	FSLinalg::RealRowVector<3> a({2, 3, 0});
	
	const auto expr1 = 0.5*FSLinalg::outer(e0, a)*a;
	
	fmt::print("expr = {}\n", FSLinalg::RealRowVector<4>(expr1));
	fmt::print("expr.createTemporaryMatrix = {}\n", expr1.createTemporaryLhs);
	fmt::print("expr.createTemporaryVector = {}\n", expr1.createTemporaryRhs);
	
	const auto expr2 = 0.5*FSLinalg::transpose(FSLinalg::outer(e0, a))*e0;
	
	fmt::print("expr = {}\n", FSLinalg::RealRowVector<3>(expr2));
	fmt::print("expr.createTemporaryMatrix = {}\n", expr2.createTemporaryLhs);
	fmt::print("expr.createTemporaryVector = {}\n", expr2.createTemporaryRhs);
	
	FSLinalg::RealMatrix<4, 3> A(FSLinalg::outer(e0, a));
	
	fmt::print("{}\n", A);
	
	const auto expr3 = 0.5*FSLinalg::transpose(A)*e0;
	
	fmt::print("expr = {}\n", FSLinalg::RealRowVector<3>(expr3));
	fmt::print("expr.createTemporaryMatrix = {}\n", expr3.createTemporaryLhs);
	fmt::print("expr.createTemporaryVector = {}\n", expr3.createTemporaryRhs);
	
	return EXIT_SUCCESS;
}
