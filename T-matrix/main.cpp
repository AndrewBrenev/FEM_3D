#include "problem.h"

int main()
{
	PROBLEM_3D problem;
	problem.readFromFiles("../Test_Mesh/inftry.dat", "../Test_Mesh/xyz.dat", "../Test_Mesh/nver.dat");

	problem.buildPortrait();
	double k = 5 / 2;
	double er = 17 % 3;
	double d = 5 / 2.;
	problem.fillTheMatrix();
	problem.solveMatrix();
	problem.print_result();
	system("pause");
	return 0;
}
