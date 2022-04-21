#pragma once
#include "Polygons.h"
#include "Mesh.h"


struct parameters {
    // main
    double ro, p, h, H, u, v, w, T, E, e;

    // transport
    double mu, la, Pr;

    double Cp, Gm, Gam;

    // vectors
    double* U, * U1;
    double* V; // ro, u,v, h

};

struct changes
{
    double* dU;
};

struct Gradient {
    Vector* g;
};

void Init(parameters* (&p), int nCells, int Nm);

void Viscous(parameters* p, changes* du, Mesh mesh,Cell* cells, double dt);
void Viscous(parameters* p, changes* du, Mesh mesh, Cell* cells, double dt, Gradient* gr, int Nm);

void GetParams(parameters* (&p), int nCells, int Nm);

void Convect(parameters* (&p), changes* (&du), Mesh mesh, Cell* cells, int It, double dt);

void Yw(Mesh mesh, Cell* (&cells), int nCells);

void Tecplot(parameters* p, Cell* cells, int Nx, int Ny, int nCells);

double Dist(Pnt A, Pnt B, Pnt E); //расстояние от точки Е до грани АВ

void SetGran(Mesh& mesh);

void Gradients(Cell* cells, Mesh mesh, Gradient* (&gr), parameters* p, int Nm);

void Matrix_Diag(double A[4][4], double L[4], double B[4][4]);
void Matrix_Matrix(double A[4][4], double B[4][4], double C[4][4]);
void Matrix_Vector(double A[4][4], double B[4], double C[4]);

void PrintMatrix(double A[4][4]);

void ConvectNS();

void S_Matr(double A[4][4]);

void S_Matr(double A[4][4], double u, double v, double ro, double p, double h, double ns, double ny, int iMod);