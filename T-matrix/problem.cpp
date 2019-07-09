//---------------------------------------------------------------------------
#include "problem.h"
//---------------------------------------------------------------------------
double c[EL_SIZE][EL_SIZE]={{8},{4,8},{4,2,8},{2,4,4,8},{4,2,2,1,8},
{2,4,1,2,4,8},{2,1,4,2,4,2,8},{1,2,2,4,2,4,4,8}};
double bx[EL_SIZE][EL_SIZE]={{4},{-4,4},{2,-2,4},{-2,2,-4,4},{2,-2,1,-1,4,4},
{-2,2,-1,1,-4,4},{1,-1,2,-2,2,-2,4},{-1,1,-2,2,-2,2,-4,4}};
double by[EL_SIZE][EL_SIZE]={{4},{2,4},{-4,-2,4},{-2,-4,2,4},{2,1,-2,-1,4},
{1,2,-1,-2,2,4},{-2,-1,2,1,-4,-2,4},{-1,-2,1,2,-2,-4,2,4}};
double bz[EL_SIZE][EL_SIZE]={{4},{2,4},{2,1,4},{1,2,2,4},{-4,-2,-2,-1,4},
{-2,-4,-1,-2,2,4},{-2,-1,-4,-2,2,1,4},{-1,-2,-2,-4,1,2,2,4}};
double right_vector[EL_SIZE];
double c1[EL_SIZE /2][EL_SIZE /2]={{4},{2,4},{2,1,4},{1,2,2,4}};
double rf_v[EL_SIZE][EL_SIZE]={{-4,4,-2,2,-2,2,-1,1},{-4,4,-2,2,-2,2,-1,1},
{-2,2,-4,4,-1,1,-2,2},{-2,2,-4,4,-1,1,-2,2}};

// div( lambga grad(u)) = right

double right(double *x)
{
	return (x[0] + x[1] + x[2]);

}
double u_real(double *x)
{
	return (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

}

//--------------------------------
void mul_c_vector(double *vector, double *result)
{
	for (int i = 0; i < EL_SIZE; i++) {
		result[i] = 0;

		for (int j = 0; j < EL_SIZE; j++) {
			int i1, j1;
			if (i > j) {
				i1 = i;
				j1 = j;
			}
			else {
				i1 = j;
				j1 = i;
			}
			result[i] += c[i1][j1] * vector[j];
		}
	}
}
//--------------------------------
void mul_c1_vector(double *vector, double *result)
{
	for (int i = 0; i < 4; i++) {
		result[i] = 0;

		for (int j = 0; j < 4; j++) {
			int i1, j1;
			if (i > j) {
				i1 = i;
				j1 = j;
			}
			else {
				i1 = j;
				j1 = i;
			}
			result[i] += c1[i1][j1] * vector[j];
		}
	}
}

int ELEM::count=0;
int NODE::count=0;

MATRIX::~MATRIX()
{
	di.clear();
	ig.clear();
	gg.clear();
	jg.clear();
}
 
//Считывание нумерации
void PROBLEM_3D::read_ELEM()
{
	FILE *f=fopen(FILE_ELEM_NUMBER,"r");
	fscanf(f,"%d",&ELEM::count);
	elems = new ELEM[ELEM::count];
	
	for(int i=0;i<ELEM::count;i++) {
		for(int j=0;j<EL_SIZE;j++)
			fscanf(f,"%d",&elems[i].mas[j]);
		fscanf(f,"%lf %lf",&elems[i].lambda,&elems[i].gamma);
	}

	fclose(f);
}
//----------------------------------------------------
//Считывание координат узлов
void PROBLEM_3D::read_NODE()
{
	FILE *f=fopen(FILE_NODE,"r");
	fscanf(f," %d",&NODE::count);
	n=NODE::count;	
	nodes = new NODE[n];

	for(int i=0;i<NODE::count;i++)
		for(int j=0;j<3;j++)
			fscanf(f,"%lf",&nodes[i].mas[j]);

	for(int i=0;i<NODE::count;i++)
	{
	printf("\n");
	for(int j=0;j<3;j++)
	printf("%2.1f ",nodes[i].mas[j]);
	} 
	fclose(f);
}
//------------------------------------------------------
void PROBLEM_3D::buildPortrait()
{
	int **prof;
	prof = new int*[n-1];
	for(int i=0;i<n-1;i++)		
		prof[i] = new int[i+1];

	for(int i=0;i<n-1;i++)
		for(int j=0;j<=i;j++)
			prof[i][j]=0;

	int s=0;//Кол-во ненулевых эл-тов
	for(int k=0;k<ELEM::count;k++)
		for(int i=0;i<EL_SIZE;i++)
			for(int j=i+1;j<EL_SIZE;j++) {
				int i1 = elems[k].mas[i], j1 = elems[k].mas[j], k_b;
				if(i1<j1) {
					k_b=i1;
					i1=j1;
					j1=k_b;
				}
				if(prof[i1-2][j1-1]==0) {
					s++;
					prof[i1-2][j1-1]=1;
				}
	}

	//Формирование массива ig и jg	
	jg.resize(s);
	gg.resize(s);
	ig[0]=1;
	ig[1]=1;
	for(int i=0,d=0;i<n-1;i++) {
		int k=0;
		for(int j=0;j<=i;j++)
			if(prof[i][j]==1) {
				jg[d]=j+1;
				d++;
				k++;
			}
		ig[i+2]=ig[i+1]+k;
	}
	/*
	printf("\n");
	for(int i =0; i <= n; i++)
		printf(" %d",ig[i]);
	printf("\n");
	for(int i =0; i < s; i++)
		printf(" %d",jg[i]);
		*/
	for(int i=0;i<n-1;i++)		
		delete [] prof[i];	
	delete [] prof;
}
//------------------------------------------------------
void PROBLEM_3D::add(int i,int j,double x)
{
	int k;
	for(k=ig[i]-1;k<ig[i+1]-1;k++)
		if(jg[k]==j+1)
			break;
	gg[k]+=x;
}
//------------------------------------------------------
void PROBLEM_3D::fillTheMatrix()
{
	double h[3];
	for (int i = 0; i < n; i++) {
		di[i] = 0.0;
		f[i] = 0.0;
	}
	for (int i = 0; i < ig[n] - 1; i++) {
		gg[i] = 0.0;
	}
	for (int k = 0; k < ELEM::count; k++)
	{
		//Вычисление шага
		//h----------------------------------------------------------
		for (int i = 0; i < 3; i++) { //i-по x,y,z    
			int flag = 1;
			int j;
			for (j = 1; j< EL_SIZE && flag; j++) {//1 узел фиксируем,пробегаем по остальным    
				flag = 0;
				for (int l = 0; l < 3 && !flag; l++)//проверяем,лежат ли точки на нужном ребре
					if (i != l && nodes[elems[k].mas[0] - 1].mas[l] != nodes[elems[k].mas[j] - 1].mas[l])
						flag = 1;
			}
			if (!flag)
				h[i] = fabs(nodes[elems[k].mas[0] - 1].mas[i] - nodes[elems[k].mas[j - 1] - 1].mas[i]);
		}
		//------------------------------------------------------------
		//формирование элементов матрицы
		//заполнение массива gg
		double b_k = elems[k].lambda*h[0] * h[1] * h[2] / 36.;
		double c_k = h[0] * h[1] * h[2] / 216.;
		//!с_k заменить

		//!
		//double c_k=1,b_k=0;
		double fr[EL_SIZE];

		//вектор правой части для локал. матрицы
		for (int i = 0; i< EL_SIZE; i++)
			right_vector[i] = right(nodes[elems[k].mas[i] - 1].mas);

		mul_c_vector(right_vector, fr);

		for (int i = 0; i< EL_SIZE; i++) {
			for (int j = 0; j < i; j++) {
				int i1 = elems[k].mas[i] - 1;
				int j1 = elems[k].mas[j] - 1;
				double s = elems[k].gamma*c_k*c[i][j] + b_k / (h[0] * h[0])*bx[i][j] +
					b_k / (h[1] * h[1])*by[i][j] + b_k / (h[2] * h[2])*bz[i][j];

				//добавка в gg
				if (i1 < j1)
					add(j1, i1, s);
				else
					add(i1, j1, s);
			}
			//добавка в диагональ
			di[elems[k].mas[i] - 1] += c_k * elems[k].gamma*c[i][i] + b_k / (h[0] * h[0])*bx[i][i] +
				b_k / (h[1] * h[1])*by[i][i] + b_k / (h[2] * h[2])*bz[i][i];
			//добавка в правую часть
			f[elems[k].mas[i] - 1] += c_k * fr[i];
		}
	}
}
//------------------------------------------------------
PROBLEM_3D:: PROBLEM_3D()
{	

} 

PROBLEM_3D::~PROBLEM_3D()
{
	delete [] elems;
	delete [] nodes;
	f.clear();
	x.clear();
}

//--------------------------------------
void PROBLEM_3D::print_result()
{
	FILE *f=fopen("resault.txt","w");
	for (int i = 0; i < n; i++)
	{
		double *real_f = new double[n];
		real_f[i] = right(nodes[i].mas);

		fprintf(f, "%.3lf %.3lf %.3lf %.15lf %.15lf\n", nodes[i].mas[0], nodes[i].mas[1], nodes[i].mas[2], x[i] , real_f[i]);
	}
	fclose(f);
}
//---------------------------------------
void PROBLEM_3D::kraev_1()
{
	int i;
	double q;
	long j, k;

	FILE *fl=fopen(FILE_CRAEV_1,"r");
	if (fl==NULL) return;
	int flag;
	while((flag=fscanf(fl,"%d",&i))>0) {
		q = right(nodes[i].mas);
		di[i]=1;
		f[i]=q;


		for (j = ig[i]; j < ig[i + 1]; j++)
		{
			f[jg[j-1]] -= gg[j-1] * q;
			gg[j-1] = 0.0;
		}

		for (j = i + 1; j < n; j++)
		{
			for (k = ig[j]; k < ig[j + 1]; k++)
				if (jg[k] == i)
				{
					f[j] -= gg[k] * q;
					gg[k] = 0.0;
				}
		}
	}
}
//----------------------------------------
void PROBLEM_3D::kraev_3()
{	
	int cr[4];//граничные узлы
	double teta[4];//значения потока в узлах
	double betta;
	FILE *fl=fopen(FILE_CRAEV_3,"r");
	if (fl==NULL) return;
	while(!feof(fl)) {

		for(int i = 0;i < 4;i++)
			fscanf(fl," %d",&cr[i]);

		for(int i = 0;i < 4;i++)
			fscanf(fl," %lf",&teta[i]);
		
		fscanf(fl," %lf",&betta);
					
		double hxyz = 1.0;
		double eps = 1e-12;
		double hx = 0.0, hy = 0.0;
		for(int i = 0; i < 3; i++) {
			double tmp = fabs(nodes[cr[0]-1].mas[i]-nodes[cr[1]-1].mas[i]);
			if(tmp > eps)
				hx += tmp;
			tmp = fabs(nodes[cr[0]-1].mas[i]-nodes[cr[2]-1].mas[i]);
			if(tmp > eps)
				hy += tmp;
		}
		hxyz = hx * hy;
			
		for(int i = 0; i < 4; i++) {
			di[cr[i]-1] += betta * hxyz/36.0 *c1[i][i];
			for(int j = 0; j < i; j++) {
				int ind= ig[cr[i]-1]-1;
				while(jg[ind] != cr[j])
					ind++;
				gg[ind] += betta * hxyz/36.0 *c1[i][j];
			}
		}

		double vec_tet[4];
		mul_c1_vector(teta,vec_tet);
		
		for(int i=0;i<4;i++)
			f[cr[i]-1] += hxyz/36.0 *vec_tet[i];
	}
	
}
//----------------------------------------
void PROBLEM_3D::kraev_2()
{
	int cr[4];//граничные узлы
	double teta[4];//значения потока в узлах
	FILE *fl=fopen(FILE_CRAEV_2,"r");
	if (fl==NULL) return;
	while(!feof(fl)) {

		for(int i=0;i<4;i++)
			fscanf(fl," %d",&cr[i]);

		for(int i=0;i<4;i++)
			fscanf(fl," %lf",&teta[i]);
		
		double hxyz = 1.0;
		double eps=1e-12;
		double hx = 0.0, hy = 0.0;
		for(int i = 0; i < 3; i++) {
			double tmp = fabs(nodes[cr[0]-1].mas[i]-nodes[cr[1]-1].mas[i]);
			if(tmp > eps)
				hx += tmp;
			tmp = fabs(nodes[cr[0]-1].mas[i]-nodes[cr[2]-1].mas[i]);
			if(tmp > eps)
				hy += tmp;
		}
		hxyz = hx * hy;

		double vec_tet[4];
		mul_c1_vector(teta,vec_tet);
		//=======
		/*printf("\n=========vt============\n");
		for(int i = 0; i < 4; i++)
			printf(" %lf", vec_tet[i]);
		printf("\n=====================\n");*/
		for(int i=0;i<4;i++)
			f[cr[i]-1] += hxyz/36.0 *vec_tet[i];
	}
}
//----------------------------------------
// получение значения базисной ф-ции
double PROBLEM_3D::get_psi(int num_fe, int num_basis, double x,double y, double z)
{
    int n1,n2,n3,n5;
    double hx,hy,hz;
    n1=elems[num_fe].mas[0];
    n2=elems[num_fe].mas[1];
    n3=elems[num_fe].mas[2];
    n5=elems[num_fe].mas[4];
    hx=nodes[n2-1].mas[0]-nodes[n1-1].mas[0];
    hy=nodes[n3-1].mas[1]-nodes[n1-1].mas[1];
    hz=nodes[n5-1].mas[2]-nodes[n1-1].mas[2];
    switch(num_basis) {
	case(0):
		return (nodes[n2-1].mas[0]-x)/hx*(nodes[n3-1].mas[1]-y)/hy*(nodes[n5-1].mas[2]-z)/hz;
		break;
	case(1):
		return (x-nodes[n1-1].mas[0])/hx*(nodes[n3-1].mas[1]-y)/hy*(nodes[n5-1].mas[2]-z)/hz;
		break;
	case(2):
		return (nodes[n2-1].mas[0]-x)/hx*(y-nodes[n1-1].mas[1])/hy*(nodes[n5-1].mas[2]-z)/hz;
		break;
	case(3):
		return (x-nodes[n1-1].mas[0])/hx*(y-nodes[n1-1].mas[1])/hy*(nodes[n5-1].mas[2]-z)/hz;
		break;
	case(4):
		return (nodes[n2-1].mas[0]-x)/hx*(nodes[n3-1].mas[1]-y)/hy*(z-nodes[n1-1].mas[2])/hz;
		break;
	case(5):
		return (x-nodes[n1-1].mas[0])/hx*(nodes[n3-1].mas[1]-y)/hy*(z-nodes[n1-1].mas[2])/hz;
		break;
	case(6):
		return (nodes[n2-1].mas[0]-x)/hx*(y-nodes[n1-1].mas[1])/hy*(z-nodes[n1-1].mas[2])/hz;
		break;
	case(7):
		return (x-nodes[n1-1].mas[0])/hx*(y-nodes[n1-1].mas[1])/hy*(z-nodes[n1-1].mas[2])/hz;
		break;
	}
}

//----------------------------------------
double PROBLEM_3D::U(double nx, double y, double z)
{
	bool flag = false;
	int i;
	//находим элемент
	for(i = 0; i < ELEM::count && !flag; i++) {
		int v1  = elems[i].mas[0];
		int v2  = elems[i].mas[1];
		int v3  = elems[i].mas[2];
		int v5  = elems[i].mas[4];
		if(nx >= nodes[v1-1].mas[0] && nx <= nodes[v2-1].mas[0] &&
			y >= nodes[v1-1].mas[1] && y <= nodes[v3-1].mas[1] &&
			z >= nodes[v1-1].mas[2] && z <= nodes[v5-1].mas[2])
			flag = true;
	}
	int lc = i;
	double sum = 0;
	for(i = 0; i< EL_SIZE; i++)
		sum += x[elems[lc-1].mas[i]-1] * get_psi(lc-1,i,nx,y,z);
	return sum;
}