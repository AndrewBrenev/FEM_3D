#include "problem.h"

int main()
{
	PROBLEM_3D problem;
	problem.readFromFiles("../Test_Mesh/inftry.dat", "../Test_Mesh/xyz.dat", "../Test_Mesh/nver.dat");

	problem.buildPortrait();
	problem.fillTheMatrix();
	problem.solveMatrix();
	problem.print_result();
	system("pause");
	return 0;
}
