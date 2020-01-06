#include "Grid.h"

bool Grid::readGridFromFiles(const char* inf, const char* xyz, const char* nver)
{
	FILE* file = fopen(inf, "r");
	if (file == NULL) return false;
	int nodes_size, el_size;
	fscanf(file, " ISLAU=       0 INDKU1=       0 INDFPO=       0\nKUZLOV=%8d   KPAR=%8d    KT1=       0   KTR2=       0   KTR3=       0\nKISRS1=       0 KISRS2=       0 KISRS3=       0   KBRS=       0\n   KT7=       0   KT10=       0   KTR4=       0  KTSIM=       0\n   KT6=       0\n", &nodes_size, &el_size);
	fclose(file);

	n = nodes_size;
	
	nodes.resize(nodes_size);
	nodes.reserve(nodes_size);

	elems.reserve(el_size);
	elems.resize(el_size);

	int* tmp = new int[el_size];
	size_t t;
	FILE* fout;

	double* tmp1 = new double[3 * nodes_size];
	if (fopen_s(&fout, xyz, "rb")) return false;
	t = fread(tmp1, sizeof(double), 3 * nodes_size, fout);
	int k = feof(fout);
	fclose(fout);

	if (fopen_s(&fout, nver, "rb")) return false;
	int* tmp2 = new int[el_size * 14];
	t = fread(tmp2, sizeof(int), el_size * 14, fout);

	//ѕреобразуем считанные данные в массивы
	double a[3];
	for (int i = 0, k = 0; i < 3 * nodes_size; i += 3, k++) {
		nodes[k].value[0] = tmp1[i];
		nodes[k].value[1] = tmp1[i + 1];
		nodes[k].value[2] = tmp1[i + 2];
	}

	for (int i = 0, k = 0; k < el_size; k++, i = i + 14) {
		for (int j = 0; j < EL_SIZE; j++) {
			elems[k].value[j] = tmp2[i + j] - 1;
			elems[k].lambda = LAMBDA;
			elems[k].gamma = GAMMA;
		}
	}
	return true;
};

bool Grid::buildTestInconsistentGrid()
{
	n = 40;
	nodes.push_back(NODE(0, 3, 0, 27));
	nodes.push_back(NODE(1, 3, 0, 28));
	nodes.push_back(NODE(2, 3, 0, 29));
	nodes.push_back(NODE(0, 3, 1, 30));
	nodes.push_back(NODE(1, 3, 1, 31));
	nodes.push_back(NODE(2, 3, 1, 32));
	nodes.push_back(NODE(0, 3, 2, 33));
	nodes.push_back(NODE(1, 3, 2, 34));
	nodes.push_back(NODE(2, 3, 2, 35));
	nodes.push_back(NODE(0.5, 3, 2, 36));


	nodes.push_back(NODE(0.5, 2, 2, 37));
	nodes.push_back(NODE(0.5, 3, 1, 38));
	nodes.push_back(NODE(0.5, 2, 1, 39));
	
	ELEM frst, scnd, thrd, frth, ffth;

	frst.value[0] = 22;
	frst.value[1] = 19;
	frst.value[2] = 21;
	frst.value[3] = 18;
	frst.value[4] = 31;
	frst.value[5] = 28;
	frst.value[6] = 30;
	frst.value[7] = 27;
	frst.lambda = LAMBDA;
	frst.gamma = GAMMA;
	elems.push_back(frst);

	scnd.value[0] = 23;
	scnd.value[1] = 20;
	scnd.value[2] = 22;
	scnd.value[3] = 19;
	scnd.value[4] = 32;
	scnd.value[5] = 29;
	scnd.value[6] = 31;
	scnd.value[7] = 28;
	scnd.lambda = LAMBDA;
	scnd.gamma = GAMMA;
	elems.push_back(scnd);

	thrd.value[0] = 26;
	thrd.value[1] = 23;
	thrd.value[2] = 25;
	thrd.value[3] = 22;
	thrd.value[4] = 35;
	thrd.value[5] = 32;
	thrd.value[6] = 34;
	thrd.value[7] = 31;
	thrd.lambda = LAMBDA;
	thrd.gamma = GAMMA;
	elems.push_back(thrd);

	frth.value[0] = 25;
	frth.value[1] = 22;
	frth.value[2] = 37;
	frth.value[3] = 39;
	frth.value[4] = 34;
	frth.value[5] = 31;
	frth.value[6] = 36;
	frth.value[7] = 38;
	frth.lambda = LAMBDA;
	frth.gamma = GAMMA;
	elems.push_back(frth); 
	
	ffth.value[0] = 37;
	ffth.value[1] = 39;
	ffth.value[2] = 24;
	ffth.value[3] = 21;
	ffth.value[4] = 36;
	ffth.value[5] = 38;
	ffth.value[6] = 33;
	ffth.value[7] = 30;
	ffth.lambda = LAMBDA;
	ffth.gamma = GAMMA;
	elems.push_back(ffth);
	/*
	n = 20;
	nodes.push_back(NODE(0, 0, 0, 0));
	nodes.push_back(NODE(1, 0, 0, 1));
	nodes.push_back(NODE(3, 0, 0, 2));
	nodes.push_back(NODE(3, 3, 0, 3));
	nodes.push_back(NODE(0, 5, 0, 4));
	nodes.push_back(NODE(1, 5, 0, 5));
	nodes.push_back(NODE(3, 5, 0, 6));
	nodes.push_back(NODE(3, 0, 1, 7));
	nodes.push_back(NODE(0, 0, 4, 8));
	nodes.push_back(NODE(1, 0, 4, 9));
	nodes.push_back(NODE(3, 0, 4, 10));
	nodes.push_back(NODE(3, 3, 4, 11));
	nodes.push_back(NODE(0, 5, 4, 12));
	nodes.push_back(NODE(1, 5, 4, 13));
	nodes.push_back(NODE(3, 5, 4, 14));
	nodes.push_back(NODE(1, 3, 0, 15));
	nodes.push_back(NODE(1, 0, 1, 16));
	nodes.push_back(NODE(1, 3, 1, 17));
	nodes.push_back(NODE(3, 3, 1, 18));
	nodes.push_back(NODE(1, 3, 1, 19));


	elems.resize(4);

	elems[0].value[0] = 9;
	elems[0].value[1] = 1;
	elems[0].value[2] = 8;
	elems[0].value[3] = 0;
	elems[0].value[4] = 13;
	elems[0].value[5] = 5;
	elems[0].value[6] = 12;
	elems[0].value[7] = 4;
	elems[0].lambda = LAMBDA;
	elems[0].gamma = GAMMA;

	elems[1].value[0] = 10;
	elems[1].value[1] = 7;
	elems[1].value[2] = 9;
	elems[1].value[3] = 16;
	elems[1].value[4] = 11;
	elems[1].value[5] = 18;
	elems[1].value[6] = 19;
	elems[1].value[7] = 17;
	elems[1].lambda = LAMBDA;
	elems[1].gamma = GAMMA;

	elems[2].value[0] = 7;
	elems[2].value[1] = 2;
	elems[2].value[2] = 16;
	elems[2].value[3] = 1;
	elems[2].value[4] = 18;
	elems[2].value[5] = 3;
	elems[2].value[6] = 17;
	elems[2].value[7] = 15;
	elems[2].lambda = LAMBDA;
	elems[2].gamma = GAMMA;

	elems[3].value[0] = 11;
	elems[3].value[1] = 3;
	elems[3].value[2] = 19;
	elems[3].value[3] = 15;
	elems[3].value[4] = 14;
	elems[3].value[5] = 6;
	elems[3].value[6] = 13;
	elems[3].value[7] = 5;
	elems[3].lambda = LAMBDA;
	elems[3].gamma = GAMMA;
	*/
	return true;
};