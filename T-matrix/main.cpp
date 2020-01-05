#include "problem.h"

int main()
{

	REGION outter, object;

	outter.X0 = 0; outter.Xn = 20; outter.Xh = 1;
	outter.Y0 = 0; outter.Yn = 20; outter.Yh = 1;
	outter.Z0 = 0; outter.Zn = 20; outter.Zh = 1;

	object.X0 = 0; object.Xn = 1; object.Xh = 1;
	object.Y0 = 0; object.Yn = 1; object.Yh = 1;
	object.Z0 = 0; object.Zn = 1; object.Zh = 1;

	PROBLEM_3D problem(outter, object);

//	problem.buildGrid();

	problem.readFromFiles("../Test_Mesh/inftry.dat", "../Test_Mesh/xyz.dat", "../Test_Mesh/nver.dat");


	problem.buildPortrait();
	problem.fillTheMatrix();
	problem.solveMatrix();
	problem.print_result();

	return 0;
}
