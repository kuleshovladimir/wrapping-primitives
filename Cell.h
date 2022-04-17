#pragma once
#include "Polygons.h"


class Cell
{
private:
	Pnt MassC;  // ���������� ������ �������
	double S;   // ������� ������
	
	
	int nFaces; // ����� ������, ���������� ������
	int nNodes; // ����� �����, ���������� ������
	

	int* faces; //������ ������, ���������� ������ (nFaces)
	int* nodes; //������ �����, ���������� ������ (nNodes)

	bool* fType; //���� ������, ���������� ������ (0 - �����.�����, 1 - ��������� �����)

	int* cells; //������ �������� ����� (nFaces)
				//���� fType = 0 - ����� ������
				//	   fType = 1 - ����� �����


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

