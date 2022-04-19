#pragma once
#include "Polygons.h"

struct Vector {
	double* cx;
};


class Cell
{
	friend class Mesh;
private:
	Pnt MassC;  // координаты центра т€жести
	double S;   // площадь €чейки


	int nFaces; // число граней, окружающих €чейку
	int nNodes; // число узлов, окружающих €чейку


	int* faces; //номера граней, окружающих €чейку (nFaces)
	int* nodes; //номера узлов, окружающих €чейку (nNodes)

	bool* fType; //типы граней, окружающих €чейку (0 - внутр.грань, 1 - гранична€ грань)

	int* cells; //номера соседних €чеек (nFaces)
				//если fType = 0 - номер €чейки
				//	   fType = 1 - номер грани

	double* wk;
	Vector* ck;

public:
	Cell();
	~Cell();

	Pnt Get_MassC();
	void Set_MassC(Pnt c_);

	void Set_S(double S_);
	double Get_S();

	void Set_NFaces(int NF);
	int Get_NFaces();

	void Set_NNodes(int NN);
	int Get_NNodes();

	int Get_Node(int i);

	int Get_Face(int i);

	void Set_Face(int iFace, int fIndex);

	void Set_Nodes(int* nodes, int nNodes);

	void Set_fType(int iFace, bool ftype);

	void Set_cells(int iFace, int cell);
	int Get_Cell(int i);

	void Print(int i);

	bool Get_fType(int i);

	double Get_wk(int i);
	Vector Get_ck(int i);
};

