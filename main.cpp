#include <chrono>

#include "Polygons.h"
#include "Cell.h"
#include "Mesh.h"
#include "Functions.h"

int main()
{
	// русификаци€
	
	setlocale(LC_ALL, "Rus");

	Mesh mesh;

	//string filename = "grid.txt";
	string filename = "grid_head_fine.dat";
	//string filename = "grid_rectangle.dat";

	mesh.ReadStruct(filename);

	int nCells = mesh.nCells;

	//cout << "nCells = " << nCells << endl;

	Cell* cells = new Cell[nCells];

	mesh.CreateCell(cells);

	mesh.CreateFaces();


	mesh.CellFuncs(cells);

	remove("cells.txt");

	//for (int i = 0; i < nCells; i++) {
	//	cells[i].Print(i);
	//}

	// ¬вод исходных данных
	parameters* p = new parameters[nCells];
	changes* du = new changes[nCells];

	// Initialization
	int Nm = 4;				//  NS

	Init(p, nCells, Nm);

	//exit(1);

	mesh.SetZones();

	//int nZones = 4;
	//Zone* z = new Zone[nZones];

	SetGran(mesh);

	Yw(mesh, cells, nCells);

	mesh.GradCoeffs(cells);

	int ItMax = 2000;	// 0000; //		50000; //0000000;	// 000;
	int It = 0;

	double resMin = 1.e-6;

	double res = 1.;

	for (int i = 0; i < nCells; i++) {
		du[i].dU = new double[Nm];
		//du[i].dU[0] = 0.;
	}

	// выделение пам€ти под градиенты параметров
	Gradient* gr = new Gradient[nCells];
	for (int i = 0; i < nCells; i++) {
		gr[i].g = new Vector[Nm];
		for (int j = 0; j < Nm; j++)
			gr[i].g[j].cx = new double[2];
	}

	auto start = std::chrono::system_clock::now();

	//ConvectNS();
	//exit(19);


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
		//if(It>1000) dt = 0.00001;

		//exit(5);

		// ѕриращение за счет нев€зких потоков
		//Convect(p, du, mesh, cells, It, dt);
		ConvectNS(p, du, mesh, cells, dt, Nm);

		// ѕриращение за счет в€зкости
		//Viscous(p, du, mesh, cells, dt);
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


	int Nx = mesh.Nx;
	int Ny = mesh.Ny;

	Tecplot(p, cells, Nx, Ny, nCells);
	Velocity(p, cells, Nx, Ny, nCells);
	Pressure(p, cells, Nx, Ny, nCells);
	Mach(p, cells, Nx, Ny, nCells);
	//#####################################################################  !!!!!!!!!!!!!!!!!!!
		// создаем поток дл€ записи
	string f = "Grad.plt";
	ofstream record(f, ios::out);
	if (record) {
		record << "VARIABLES = \"X\", \"Y\", \"hGr_x\", \"hGr_y\"" << endl;

		record << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;
		for (int i = 0; i < nCells; i++) {
			record << cells[i].Get_c().x << " " << cells[i].Get_c().y
				<< " " << gr[i].g[0].cx[0] << " " << gr[i].g[0].cx[1] << endl;
		}

	}
	record.close();
	//###################################################################### !!!!!!!!!!!!!!!!!!!!!!



}