#include "Polygons.h"
#include "Cell.h"
#include "Mesh.h"

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	Mesh mesh;

	string filename = "grid.txt";

	mesh.ReadStruct(filename);

	int nCells = mesh.Get_nCells();

	Cell* cells = new Cell[nCells];

	mesh.CreateCells(cells);

	
	mesh.CreateFases();

	mesh.CellFuncs(cells);

	remove("cells.txt");

	for (int i = 0; i < nCells; i++) {
		cells[i].Print(i);
	}

	return 0;
}