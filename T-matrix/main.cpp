#include "problem.h"


int main()
{
	bool t_mesh = false;

	PROBLEM_3D problem(t_mesh);
	
	problem.readGridFromFiles(mesh_inftry_path.c_str(), mesh_xyz_path.c_str(), mesh_nver_path.c_str());

	//problem.buildPortrait();
	problem.buildSetBasedPortrait();
	problem.fillTheMatrix(); 
	problem.calculateRightPart();
	problem.applyFirstEdge();
	problem.solveMatrix();
	problem.showResault();
	problem.printResult("resault.txt");

 	return 0;
}
