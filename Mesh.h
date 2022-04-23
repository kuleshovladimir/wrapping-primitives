#pragma once
#include "Cell.h"

#include "Cell.h"

struct Face {
	int nodes[2];		// индексы узлов грани	(2)
	int cr;				// индекс правой ячейки
	int cl;				// индекс левой ячейки
	bool is_boundary;	// является ли эта грань граничной
	Pnt f_centr;		// координаты центра грани
	double length;		// длина грани
	int zone;			// номер зоны, на которой устанавливаются граничные условия

	void Print(int k) {
		cout << "			face index= " << k << endl;
		cout << " nodes indexes = " << nodes[0] << ", " << nodes[1] << endl;
		cout << " cr index = " << cr << ", cl index = " << cl << endl;
		cout << " center = " << f_centr.x << ", " << f_centr.y << endl;
		cout << "  length = " << length << endl;

	}

};

struct Wall {
	int vel;	// 0 - Free slip;  1 - No slip
	int temp;	// 1 - Tw is set;  2 - qw is set
	double value;	// значение Tw или qw
};

struct Boundary {
	int* rules;
	//например, для стенки :  2 elements
		//!vel     !0 - Free slip;  1 - No slip
		//!temp    !1 - Tw is set;  2 - qw is set < ->values
	//Для сверзвукового входа : \
		//type : 1 - заданы L, U, T, p, angle < ->values(5)

	double* vals;

};

struct Zone {
	int grantype;
	//1. Wall. 2.Supersonic Inlet. 3.Symmetry
	//4. Supersonic Outlet.  5.Free boundary
	//6. Subsonic Inlet.     7. Subsonic Outlet
	Wall* wall;
	Boundary* bnd;
};

class Mesh {
private:

public:
	// переменные
	int Nx, Ny;		// число точек по x,y в структурированной сетке
	int nNodes, nFaces, nCells;
	Pnt* nodes;		// координаты узлов

	Face* faces;	// все грани сетки		***

	int nZones;
	Zone* zones;

	// функции
	Mesh();
	~Mesh();

	void ReadStruct(string filename);

	void CreateCell(Cell* (&cells));

	void CreateFaces();

	void CellFuncs(Cell* (&cells));

	void SetZones();

	void GradCoeffs(Cell* (&cells));


};

