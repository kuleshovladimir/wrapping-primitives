#pragma once
#include "Polygons.h"
#include "Mesh.h"

void Init(parameters* (&p), int nCells, int Nm);
//void Viscous(parameters* p, changes* (&du), Mesh mesh, Cell* cells, double dt);

void Viscous(parameters* p, changes* (&du), Mesh mesh, Cell* cells, double dt, Gradient* gr, int Nm);

void GetParams(parameters* (&p), int nCells, int Nm);

void Tecplot(parameters* p, Cell* cells, int Nx, int Ny, int nCells);

void Convect(parameters* p, changes* (&du), Mesh mesh, Cell* cells, int It, double dt);

void Yw(Mesh mesh, Cell* (&cells), int nCells);

double Dist(Pnt A, Pnt B, Pnt E);

void SetGran(Mesh& mesh);

void Gradients(Cell* cells, Mesh mesh, Gradient* (&gr), parameters* p, int Nm);

void Velocity(parameters* p, Cell* cells, int Nx, int Ny, int nCells);


void Matrix_Diag(double A[4][4], double L[4], double B[4][4]);
void Matrix_Matrix(double A[4][4], double B[4][4], double C[4][4]);
void Matrix_Vector(double A[4][4], double B[4], double C[4]);

void PrintMatr(double A[4][4]);

void ConvectNS();

void S_Matr(double A[4][4]);

void S_Matr(double A[4][4], double u, double v, double ro, double p, double h,
	double nx, double ny, int iMod);

void ConvectNS(parameters* p, changes* (&du), Mesh mesh, Cell* cells, double dt, int Nm);

void Pressure(parameters* p, Cell* cells, int Nx, int Ny, int nCells);
void Mach(parameters* p, Cell* cells, int Nx, int Ny, int nCells);
void Gradplot(Gradient* gr,Cell* cells, int Nx, int Ny, int nCells);
