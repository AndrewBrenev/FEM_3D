#include "myclass.h"
#include <conio.h>

int main()
{
	GLOBAL_MATRIX my;	
	//my.kraev_2();
	//my.kraev_3();
	//my.kraev_1();
	my.MSG();
	my.print_result();
	double res;
	res = my.U(1.5, 0.5, 0.5);
	return 0;
}
