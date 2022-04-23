#include <chrono>

#include "Polygons.h"
#include "Cell.h"
#include "Mesh.h"
#include "Functions.h"

int main()
{
	// русификация
	
	setlocale(LC_ALL, "Rus");

	Mesh mesh;

	string filename = "grid.dat";

	mesh.ReadStruct(filename);

	int nCells = mesh.Get_nCells();

	Cell* cells = new Cell[nCells];

	mesh.CreateCell(cells);

	mesh.CreateFaces();


	mesh.CellFuncs(cells);

	remove("cells.txt");


	// Ввод исходных данных
	parameters* p = new parameters[nCells];
	changes* du = new changes[nCells];

	// Initialization
	int Nm = 4;				//  NS

	Init(p, nCells, Nm);

	//exit(1);

	mesh.SetZones();

	SetGran(mesh);

	Yw(mesh, cells, nCells);

	mesh.GradCoeffs(cells);

	int ItMax = 2000;
	int It = 0;

	double resMin = 1.e-6;

	double res = 1.;

	for (int i = 0; i < nCells; i++) {
		du[i].dU = new double[Nm];
	}

	// выделение памяти под градиенты параметров
	Gradient* gr = new Gradient[nCells];
	for (int i = 0; i < nCells; i++) {
		gr[i].g = new Vector[Nm];
		for (int j = 0; j < Nm; j++)
			gr[i].g[j].cx = new double[2];
	}

	auto start = std::chrono::system_clock::now();

	while (It < ItMax && res > resMin) {

		It++;

		for (int i = 0; i < nCells; i++) {
			for (int j = 0; j < Nm; j++) {
				p[i].U[j] = p[i].U1[j];  // переприсвоение
				du[i].dU[j] = 0.;

			}
		}

		Gradients(cells, mesh, gr, p, Nm);  // NS

		double dt = 0.000001;


		// Приращение за счет невязких потоков
		ConvectNS(p, du, mesh, cells, dt, Nm);

		// Приращение за счет вязкости
		Viscous(p, du, mesh, cells, dt, gr, Nm);

		for (int i = 0; i < nCells; i++) {
			for (int j = 0; j < Nm; j++) {
				p[i].U1[j] = p[i].U[j] + du[i].dU[j];
			}
		}

		GetParams(p, nCells, Nm);

		res = 0.;
		for (int i = 0; i < nCells; i++) {
			double res_ = abs(du[i].dU[3] / p[i].U1[3]);
			if (res < res_) res = res_;
		}

		cout << "It= " << It << ", res= " << res << endl;

	}

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


	int Nx = mesh.Get_Nx();
	int Ny = mesh.Get_Ny();

	Tecplot(p, cells, Nx, Ny, nCells);
	Velocity(p, cells, Nx, Ny, nCells);
	Pressure(p, cells, Nx, Ny, nCells);
	Mach(p, cells, Nx, Ny, nCells);
	Gradplot(gr, cells, Nx, Ny, nCells);
}