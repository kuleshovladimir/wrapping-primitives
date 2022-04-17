#pragma once
#include "Polygons.h"


class Cell
{
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

	void Set_Nodes(int* nodes,int nNodes);

	void Set_fType(int iFace, bool ftype);

	void Set_cells(int iFace, int cell);

	void Print(int i);
};

