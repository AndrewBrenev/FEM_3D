//---------------------------------------------------------------------------
#include "myclass.h"
//---------------------------------------------------------------------------
double c[M][M]={{8},{4,8},{4,2,8},{2,4,4,8},{4,2,2,1,8},
{2,4,1,2,4,8},{2,1,4,2,4,2,8},{1,2,2,4,2,4,4,8}};
double bx[M][M]={{4},{-4,4},{2,-2,4},{-2,2,-4,4},{2,-2,1,-1,4,4},
{-2,2,-1,1,-4,4},{1,-1,2,-2,2,-2,4},{-1,1,-2,2,-2,2,-4,4}};
double by[M][M]={{4},{2,4},{-4,-2,4},{-2,-4,2,4},{2,1,-2,-1,4},
{1,2,-1,-2,2,4},{-2,-1,2,1,-4,-2,4},{-1,-2,1,2,-2,-4,2,4}};
double bz[M][M]={{4},{2,4},{2,1,4},{1,2,2,4},{-4,-2,-2,-1,4},
{-2,-4,-1,-2,2,4},{-2,-1,-4,-2,2,1,4},{-1,-2,-2,-4,1,2,2,4}};
double right_vector[M];
double c1[M/2][M/2]={{4},{2,4},{2,1,4},{1,2,2,4}};
double rf_v[M][M]={{-4,4,-2,2,-2,2,-1,1},{-4,4,-2,2,-2,2,-1,1},
{-2,2,-4,4,-1,1,-2,2},{-2,2,-4,4,-1,1,-2,2}};




double right(double *x)
{
	return (x[0]* x[0] +x[1]* x[1] +x[2]*x[2]);
	/*if(x[0] >= 1.0)
		return -18.0;
	else
		return 2*(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-6.0;*/
	//return 5*cos(x[0]+x[1]);
}
//--------------------------------
void mul_c_vector(double *vector,double *result)
{
	for(int i=0;i<M;i++) {
		result[i]=0;

		for(int j=0;j<M;j++) {
			int i1,j1;
			if(i>j) {
				i1=i;
				j1=j;
			}
			else {
				i1=j;
				j1=i;
			}
			result[i]+=c[i1][j1]*vector[j];
		}
	}
}
//--------------------------------
void mul_c1_vector(double *vector,double *result)
{
	for(int i=0;i<4;i++) {
		result[i]=0;

		for(int j=0;j<4;j++) {
			int i1,j1;
			if(i>j) {
				i1=i;
				j1=j;
			}
			else {
				i1=j;
				j1=i;
			}
			result[i]+=c1[i1][j1]*vector[j];
		}
	}	
}
//----------------------------
double scalar(double *vector1,double *vector2,int n)
{
	double f=0;
	for(int i=0;i<n;i++)
		f+=vector1[i]*vector2[i];
	return f;
}
//-----------------------------
double norma(double *vector,int n)
{
	double f=scalar(vector,vector,n);
	return sqrt(f);
}
//-----------------------------
void sum(double *vector1,double *vector2,double *result,int n)
{
	for(int i=0;i<n;i++)
		result[i]=vector1[i]+vector2[i];
}
//-----------------------------
void mul(double a,double *vector,double *result,int n)
{
	for(int i=0;i<n;i++)
		result[i]=a*vector[i];
}
//-----------------------------
void mul2(double *a,double *vector,double *result,int n)
{
	for(int i=0;i<n;i++)
		result[i]=a[i]*vector[i];
}
//-----------------------------
void mov(double *vector,double *result,int n)
{
	for(int i=0;i<n;i++)
		result[i]=vector[i];
}
//-----------------------------
int local::count=0;
int cross::count=0;
//----------------------------------
void MATRIX::mul_matrix_vector(double *vector,double *result)
{
	int i;
	for(i=0;i<n;i++)
	result[i]=d[i]*vector[i];

	for(i=0;i<n;i++)
		for(int j=ig[i];j<ig[i+1];j++) {
		result[i]+=gg[j-1]*vector[jg[j-1]-1];
		result[jg[j-1]-1]+=gg[j-1]*vector[i];
		}
}
//-------------------------------------------------------
void MATRIX::solution_x_l(double *f,double *x)
{
    int i;
	for(i=0;i<n;i++)
		x[i]=f[i];
	//SY=F
	for(i=0;i<n;i++){
		for(int j=ig[i]-1;j<ig[i+1]-1;j++)
			x[i]-=gg[j]*x[jg[j]-1];
		x[i]/=d[i];
	}
}

//------------------------------
void MATRIX::solution_x_u(double*f,double *x)
{
    int i;
	for(i=0;i<n;i++)
		x[i]=f[i];

	for(i=n-1;i>=0;i--) {
		x[i]/=d[i];
		for(int j=ig[i+1]-2;j>=ig[i]-1;j--)
			x[jg[j]-1]-=gg[j]*x[i];
	}

}
//----------------------------------
void MATRIX::solution_x(double *f,double *x)
{
    int i;
	for(i=0;i<n;i++)
		x[i]=f[i];

	//SY=F 
	for(i=0;i<n;i++) {
		for(int j=ig[i]-1;j<ig[i+1]-1;j++)
			x[i]-=gg[j]*x[jg[j]-1];
		x[i]/=d[i];
	}

	//STX=Y
	for(i=n-1;i>=0;i--){
		x[i]/=d[i];
		for(int j=ig[i+1]-2;j>=ig[i]-1;j--)
			x[jg[j]-1]-=gg[j]*x[i];
	}
}


//----------------------------------
void MATRIX::SST(MATRIX *my)
{
	int i, j, k, l;

	my->n = n;
	my->d = new double[n];
	my->ig = new int[n + 1];

	for (i = 0; i <= n; i++)
		my->ig[i] = ig[i];


	my->jg = new int[ig[n] - 1];
	my->gg = new double[ig[n] - 1];
	for (i = 0; i < ig[n] - 1; i++)
		my->jg[i] = jg[i];

	for (j = 0; j < n; j++)
	{
		double sum = 0;
		for (k = ig[j] - 1; k < ig[j + 1] - 1; k++)
			sum += my->gg[k] * my->gg[k];

		my->d[j] = sqrt(fabs(d[j] - sum));

		for (i = j + 1; i < n; i++)
		{
			int number;
			int flag = 1;
			for (l = ig[i] - 1; l < ig[i + 1] - 1 && flag; l++)
				if (jg[l] == j + 1) flag = 0;

			number = l - 1;
			if (flag) continue;

			sum = 0;
			for (k = ig[i] - 1; k < ig[i + 1] - 1 && jg[k] <= j; k++)
			{
				flag = 1;
				for (l = ig[j] - 1; l < ig[i + 1] - 1 && flag&&jg[l] <= jg[k]; l++)
					if (jg[l] == jg[k]) flag = 0;

				l--;

				if (!flag)
					sum += my->gg[l] * my->gg[k];

			}

			my->gg[number] = (gg[number] - sum) / my->d[j];
		}


	}

}
MATRIX::~MATRIX()
{
	delete [] d;
	delete [] ig;
	delete [] gg;
	delete [] jg;
}
 
//Считывание нумерации
void GLOBAL_MATRIX::read_local()
{
	FILE *f=fopen(FILE_LOCAL_NUMBER,"r");
	fscanf(f,"%d",&local::count);	
	matr = new local[local::count];

	for(int i=0;i<local::count;i++) {
		for(int j=0;j<M;j++)
			fscanf(f,"%d",&matr[i].mas[j]);
		fscanf(f,"%lf %lf",&matr[i].lambda,&matr[i].gamma);
	}

	fclose(f);
}
//----------------------------------------------------
//Считывание координат узлов
void GLOBAL_MATRIX::read_cross()
{
	FILE *f=fopen(FILE_CROSS,"r");
	fscanf(f," %d",&cross::count);
	n=cross::count;	
	set = new cross[n];

	for(int i=0;i<cross::count;i++)
		for(int j=0;j<3;j++)
			fscanf(f,"%lf",&set[i].mas[j]);

	for(int i=0;i<cross::count;i++)
	{
	printf("\n");
	for(int j=0;j<3;j++)
	printf("%2.1f ",set[i].mas[j]);
	} 
	fclose(f);
}
//------------------------------------------------------
void GLOBAL_MATRIX::formier_profil()
{
	int **prof;
	prof = new int*[n-1];
	for(int i=0;i<n-1;i++)		
		prof[i] = new int[i+1];

	for(int i=0;i<n-1;i++)
		for(int j=0;j<=i;j++)
			prof[i][j]=0;

	int s=0;//Кол-во ненулевых эл-тов
	for(int k=0;k<local::count;k++)
		for(int i=0;i<M;i++)
			for(int j=i+1;j<M;j++) {
				int i1 = matr[k].mas[i], j1 = matr[k].mas[j], k_b;
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
	jg = new int[s];
	gg = new double[s];
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
void GLOBAL_MATRIX::add(int i,int j,double x)
{
	int k;
	for(k=ig[i]-1;k<ig[i+1]-1;k++)
		if(jg[k]==j+1)
			break;
	gg[k]+=x;
}
//------------------------------------------------------
void GLOBAL_MATRIX::formier_matrix()
{
	double h[3];
	for (int i = 0; i < n; i++) {
		d[i] = 0.0;
		f[i] = 0.0;
	}
	for (int i = 0; i < ig[n] - 1; i++) {
		gg[i] = 0.0;
	}
	for (int k = 0; k < local::count; k++)
	{
		//Вычисление шага
		//h----------------------------------------------------------
		for (int i = 0; i < 3; i++) { //i-по x,y,z    
			int flag = 1;
			int j;
			for (j = 1; j < M && flag; j++) {//1 узел фиксируем,пробегаем по остальным    
				flag = 0;
				for (int l = 0; l < 3 && !flag; l++)//проверяем,лежат ли точки на нужном ребре
					if (i != l && set[matr[k].mas[0] - 1].mas[l] != set[matr[k].mas[j] - 1].mas[l])
						flag = 1;
			}
			if (!flag)
				h[i] = fabs(set[matr[k].mas[0] - 1].mas[i] - set[matr[k].mas[j - 1] - 1].mas[i]);
		}
		//------------------------------------------------------------
		//формирование элементов матрицы
		//заполнение массива gg
		double b_k = matr[k].lambda*h[0] * h[1] * h[2] / 36.;
		double c_k = h[0] * h[1] * h[2] / 216.;
		//!с_k заменить

		//!
		//double c_k=1,b_k=0;
		double fr[M];

		//вектор правой части для локал. матрицы
		for (int i = 0; i < M; i++)
			right_vector[i] = right(set[matr[k].mas[i] - 1].mas);

		mul_c_vector(right_vector, fr);

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < i; j++) {
				int i1 = matr[k].mas[i] - 1;
				int j1 = matr[k].mas[j] - 1;
				double s = matr[k].gamma*c_k*c[i][j] + b_k / (h[0] * h[0])*bx[i][j] +
					b_k / (h[1] * h[1])*by[i][j] + b_k / (h[2] * h[2])*bz[i][j];

				//добавка в gg
				if (i1 < j1)
					add(j1, i1, s);
				else
					add(i1, j1, s);
			}
			//добавка в диагональ
			d[matr[k].mas[i] - 1] += c_k * matr[k].gamma*c[i][i] + b_k / (h[0] * h[0])*bx[i][i] +
				b_k / (h[1] * h[1])*by[i][i] + b_k / (h[2] * h[2])*bz[i][i];
			//добавка в правую часть
			f[matr[k].mas[i] - 1] += c_k * fr[i];
		}
	}
}
//------------------------------------------------------
GLOBAL_MATRIX::GLOBAL_MATRIX()
{
	FILE_LOCAL_NUMBER="local.txt";
	FILE_CROSS="cross.txt";
	//FILE_OUT="out.xls";
	FILE_CRAEV_1 = "boundary.txt";

	//FILE_CRAEV_1="kraev_1.txt";
	//FILE_CRAEV_2="kraev_2.txt";
	//FILE_CRAEV_3="kraev_3.txt";
	POR=1e30;
	readFromFiles("./input/inftry.dat", "./input/xyz.dat", "./input/nver.dat");
	//read_cross();
	//read_local();
	
	d = new double[n];
	f = new double[n];
	x = new double[n];
	ig = new int[n+1];

	formier_profil();
	formier_matrix();

	/*
	for(int i=0;i<ig[n]-1;i++)
		printf("%1.2f ",gg[i]);
	printf("\n=============\n");
	for(int i=0;i<n;i++)
		printf("%1.2f ",f[i]);
	printf("\n=============\n");
	for(int i=0;i<n;i++)
		printf("%1.2f ",d[i]);
	*/
} 

GLOBAL_MATRIX::~GLOBAL_MATRIX()
{
	delete [] matr;
	delete [] set;
	delete [] f;
	delete [] x;
}
//---------------------------------------------
void GLOBAL_MATRIX::MSG()
{

	double e=1e-70;
	int max=1000;
	MATRIX s;
	SST(&s);
	
	double *r=new double[n];
	double *z=new double[n];
	double *h=new double[n];
	double *h2=new double[n];
	double *h3=new double[n];

	for(int i=0;i<n;i++) 
		x[i]=0;
	mul_matrix_vector(x,h);
	mul(-1,h,h,n);
	sum(f,h,r,n);

	s.solution_x(r,z);

	double e2=norma(r,n)/norma(f,n);

	double a,b;
	int k;
	for(k=1;k<max&&e2>e;k++) {
		//a
		s.solution_x(r,h2);
		mul_matrix_vector(z,h);
		a=scalar(h2,r,n)/scalar(h,z,n);

		//r
		mov(r,h3,n);
		mul(-a,h,h,n);
		sum(r,h,r,n);

		//x
		mul(a,z,h,n);
		sum(x,h,x,n);

		//b
		b=1/scalar(h2,h3,n);
		s.solution_x(r,h2);
		b*=scalar(h2,r,n);

		//z
		mul(b,z,z,n);
		sum(z,h2,z,n);


		e2=norma(r,n)/norma(f,n);
	}
	printf("k=%d e=%10.8e\n",k-1,e2);	
	delete [] r;
	delete [] z;
	delete [] h;
	delete [] h2;
	delete [] h3;
}
//--------------------------------------
void GLOBAL_MATRIX::print_result()
{
	FILE *f=fopen("resault.txt","w");
	for (int i = 0; i < n; i++)
	{
		double *real_f = new double[n];
		real_f[i] = right(set[i].mas);


		fprintf(f, "%.3lf %.3lf %.3lf %.15lf %.15lf\n", set[i].mas[0], set[i].mas[1], set[i].mas[0], x[i] , real_f[i]);
	}
	fclose(f);
}
//---------------------------------------
void GLOBAL_MATRIX::kraev_1()
{
	int i;
	double q;
	long j, k;

	FILE *fl=fopen(FILE_CRAEV_1,"r");
	if (fl==NULL) return;
	int flag;
	while((flag=fscanf(fl,"%d",&i))>0) {
		q = right(set[i].mas);
		d[i]=1;
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
void GLOBAL_MATRIX::kraev_3()
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
			double tmp = fabs(set[cr[0]-1].mas[i]-set[cr[1]-1].mas[i]);
			if(tmp > eps)
				hx += tmp;
			tmp = fabs(set[cr[0]-1].mas[i]-set[cr[2]-1].mas[i]);
			if(tmp > eps)
				hy += tmp;
		}
		hxyz = hx * hy;
			
		for(int i = 0; i < 4; i++) {
			d[cr[i]-1] += betta * hxyz/36.0 *c1[i][i];
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
void GLOBAL_MATRIX::kraev_2()
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
			double tmp = fabs(set[cr[0]-1].mas[i]-set[cr[1]-1].mas[i]);
			if(tmp > eps)
				hx += tmp;
			tmp = fabs(set[cr[0]-1].mas[i]-set[cr[2]-1].mas[i]);
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
double GLOBAL_MATRIX::get_psi(int num_fe, int num_basis, double x,double y, double z)
{
    int n1,n2,n3,n5;
    double hx,hy,hz;
    n1=matr[num_fe].mas[0];
    n2=matr[num_fe].mas[1];
    n3=matr[num_fe].mas[2];
    n5=matr[num_fe].mas[4];
    hx=set[n2-1].mas[0]-set[n1-1].mas[0];
    hy=set[n3-1].mas[1]-set[n1-1].mas[1];
    hz=set[n5-1].mas[2]-set[n1-1].mas[2];
    switch(num_basis) {
	case(0):
		return (set[n2-1].mas[0]-x)/hx*(set[n3-1].mas[1]-y)/hy*(set[n5-1].mas[2]-z)/hz;
		break;
	case(1):
		return (x-set[n1-1].mas[0])/hx*(set[n3-1].mas[1]-y)/hy*(set[n5-1].mas[2]-z)/hz;
		break;
	case(2):
		return (set[n2-1].mas[0]-x)/hx*(y-set[n1-1].mas[1])/hy*(set[n5-1].mas[2]-z)/hz;
		break;
	case(3):
		return (x-set[n1-1].mas[0])/hx*(y-set[n1-1].mas[1])/hy*(set[n5-1].mas[2]-z)/hz;
		break;
	case(4):
		return (set[n2-1].mas[0]-x)/hx*(set[n3-1].mas[1]-y)/hy*(z-set[n1-1].mas[2])/hz;
		break;
	case(5):
		return (x-set[n1-1].mas[0])/hx*(set[n3-1].mas[1]-y)/hy*(z-set[n1-1].mas[2])/hz;
		break;
	case(6):
		return (set[n2-1].mas[0]-x)/hx*(y-set[n1-1].mas[1])/hy*(z-set[n1-1].mas[2])/hz;
		break;
	case(7):
		return (x-set[n1-1].mas[0])/hx*(y-set[n1-1].mas[1])/hy*(z-set[n1-1].mas[2])/hz;
		break;
	}
}

//----------------------------------------
double GLOBAL_MATRIX::U(double nx, double y, double z)
{
	bool flag = false;
	int i;
	//находим элемент
	for(i = 0; i < local::count && !flag; i++) {
		int v1  = matr[i].mas[0];
		int v2  = matr[i].mas[1];
		int v3  = matr[i].mas[2];
		int v5  = matr[i].mas[4];
		if(nx >= set[v1-1].mas[0] && nx <= set[v2-1].mas[0] &&
			y >= set[v1-1].mas[1] && y <= set[v3-1].mas[1] &&
			z >= set[v1-1].mas[2] && z <= set[v5-1].mas[2])
			flag = true;
	}
	int lc = i;
	double sum = 0;
	for(i = 0; i < M; i++)
		sum += x[matr[lc-1].mas[i]-1] * get_psi(lc-1,i,nx,y,z);
	return sum;
}