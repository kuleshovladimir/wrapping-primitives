#pragma once
#include "Cell.h"

struct  Face
{
	int nodes[2];		//������� ����� �����
	int cr;				//������ ������ ������
	int cl;				//������ ����� ������
	bool is_boundary;   //�������� �� ��� ����� ���������
	Pnt f_center;		//���������� ����� ������
	double length;		//����� �����

	void Print(int k) {
		cout << "			face index = " << k << endl;
		cout << "nodes indexes = " << nodes[0] << " " << nodes[1] << endl;
		cout << "cr index = " << cr << ", cl index = " << cl << endl;
		cout << "center = " << f_center.x << " " << f_center.y << endl;
		cout << "length = " << length << endl;
	}
};

class Mesh
{
private:
	int Nx, Ny;			//����� ����� �� x, y � ����������������� �����
	int nNodes, nCells, nFaces;
	Pnt* nodes;			//���������� �����

	Face* faces;		//��� ����� �����

public:
	Mesh();
	~Mesh();

	void Set_Nx(int n);
	int Get_Nx();

	void Set_Ny(int n);
	int Get_Ny();

	void Set_nNodes(int n);
	int Get_nNodes();

	void Set_nCells(int n);
	int Get_nCells();

	void Set_nFaces(int n);
	int Get_nFaces();

	Pnt Get_node(int i);
	void Set_node(Pnt node, int i);


	void ReadStruct(string filename);

	void CreateCells(Cell* (&cells));

	void CreateFases();

	void CellFuncs(Cell* (&cells));
};

