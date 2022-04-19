#include "Polygons.h"
#include "Cell.h"
#include "Mesh.h"
#include "Functions.h"

int main()
{
	setlocale(LC_ALL, "rus");

	Mesh mesh;

	string filename = "grid.txt";

	mesh.ReadStruct(filename);

	int nCells = mesh.Get_nCells();

	Cell* cells = new Cell[nCells];

	mesh.CreateCells(cells);
	
	mesh.CreateFases();

	mesh.CellFuncs(cells);
	/*
	remove("cells.txt");
	
	for (int i = 0; i < nCells; i++) {
		cells[i].Print(i);
	}*/

	//Ввод исходных данных 

	parameters* p = new parameters[nCells];
	changes* du = new changes[nCells];

	//Initialisation

	int Nm = 1;

	Init(p, nCells, Nm);

	SetGran(mesh);

	Yw();

	double dt = 2.e-1;

	int ItMax = 50000;
	int It = 0;

	double resMin = 1.e-6;
	double res = 1.0;

	for (int i = 0; i < nCells; i++) {
		du[i].dU = new double[Nm];
	}

	Gradient* gr = new Gradient[nCells];
	for (int i = 0; i < nCells; i++) {
		gr[i].g = new Vector[Nm];
		for (int j = 0; j < Nm; j++) {
			gr[i].g[j].cx = new double[2];
		}
	}

	while ((It < ItMax) && (res > resMin)) {
		It++;
		for (int i = 0; i < nCells; i++) {
			for (int j = 0; j < Nm; j++) {
				p[i].U[j] = p[i].U1[j];
				du[i].dU[j] = 0;
			}
		}

		Gradients(cells, mesh, gr, p, Nm);


		//приращение за счет невязких потоков
		Convect(p, du, mesh, cells, It, dt);

		//приращение за счет вязкости
		Viscous(p, du, mesh, cells, dt);

		for (int i = 0; i < nCells; i++) {
			for (int j = 0; j < Nm; j++) {
				p[i].U1[j] = p[i].U[j] + du[i].dU[j];
			}
		}

		GetParams(p, nCells, Nm);

		res = 0;
		for (int i = 0; i < nCells; i++) {
			double res_ = abs(du[i].dU[0] / p[i].U1[0]);
			if (res < res_) res = res_;
		}

		int Nx = mesh.Get_Nx();
		int Ny = mesh.Get_Ny();

		Tecplot(p, cells, Nx, Ny, nCells);


	}
	cout << "It = " << It << "	res = " << res << endl;
	return 0;
}