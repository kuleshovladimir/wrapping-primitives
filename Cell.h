#pragma once
#include "Polygons.h"


struct Vector {
	double* cx;
};

class Cell {
private:
	Pnt c;		// ���������� �.�.
	double S;	// ������� ������

	int nFaces;		// ����� ������, ���������� ������

	int nNodes;		// ����� �����, ���������� ������



	int* nodes;    // ������ �����, ���������� ������ (nNodes)






public:

	double Yw;	// ���������� �� ������ ������ �� ������

	double* wk;  // weight coefficients [nFaces]
	Vector* ck;

	int* cells;		// ������ �������� ����� (nFaces)
					// fType = 0	-> ����� ������
					// fType = 1	-> ����� �����

	int* fType;	// ���� ������:  =0 -> ���������� �����	(nFaces)
								//   =1 -> ��������� �����
	int* faces;		// ������ ������, ���������� ������ (nFaces)


	Cell();
	~Cell();


	Pnt Get_c() { return c; };
	void Set_c(Pnt c_);

	void Set_S(double S_) { S = S_; };
	double Get_S() { return S; };

	void Set_nFaces(int nf);
	int Get_nFaces() { return nFaces; };

	void Set_nNodes(int nn);
	int Get_nNodes() { return nNodes; };

	//void Set_Nodes(int* nodes, int n);  // new

	int Get_Node(int i) { return nodes[i]; };

	int Get_Face(int i);

	void Set_Face(int iFace, int fIndex) { faces[iFace] = fIndex; }; // new

	void Set_fType(int iFace, int ftype) { fType[iFace] = ftype; };

	void Set_cells(int iFace, int cel) { cells[iFace] = cel; };

	void Print(int i);

	void Set_Nodes(int* nodes_, int nn);

	int Get_Cell(int i) { return cells[i]; };


};
