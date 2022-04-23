#include "Cell.h"
Cell::Cell()
{
	c.x = 0., c.y = 0.5;		// ���������� �.�.
	S = 0.;	// ������� ������

	nFaces = 0;		// ����� ������, ���������� ������

	nNodes = 0;		// ����� �����, ���������� ������

	faces = new int[nFaces];		// ������ ������, ���������� ������ (nFaces)

	nodes = new int[nNodes];    // ������ �����, ���������� ������ (nNodes)

	fType = new int[nFaces];	// ���� ������:  =0 -> ���������� �����	(nFaces)
								//   =1 -> ��������� �����

	cells = new int[nFaces];		// ������ �������� ����� (nFaces)	
}

Cell::~Cell()
{
}

void Cell::Set_c(Pnt c_)
{
	c.x = c_.x;
	c.y = c_.y;
}

void Cell::Set_nFaces(int nf)
{
	nFaces = nf;

	faces = new int[nFaces];		// ������ ������, ���������� ������ (nFaces)

	fType = new int[nFaces];	// ���� ������:  =0 -> ���������� �����	(nFaces)
								//   =1 -> ��������� �����

	cells = new int[nFaces];		// ������ �������� ����� (nFaces)

}

void Cell::Set_nNodes(int nn)
{
	nNodes = nn;

	nodes = new int[nNodes];    // ������ �����, ���������� ������ (nNodes)

}

int Cell::Get_Face(int i)
{
	if (i < nFaces)
		return faces[i];
	else
	{
		cout << "������ ����� ����� �� ����������" << endl;
		return -10;
	}
}

void Cell::Print(int m)
{
	// ������� ����� ��� ������
	string f = "cells.txt";
	ofstream record(f, ios::out | ios::app);

	if (record) {

		record << "				cell index= " << m << endl;
		record << "  center= " << c.x << ", " << c.y << endl;
		record << "  square= " << S << endl;
		record << "  number of faces= " << nFaces << endl;
		record << "  face indexes and types: " << endl;
		for (int i = 0; i < nFaces; i++) {
			record << "ind=" << faces[i] << ", type=" << fType[i] << endl;
		}

		record << "  number of nodes= " << nNodes << endl;
		record << "  node indexes: " << endl;
		for (int i = 0; i < nNodes; i++) {
			record << nodes[i];
			if (i < nNodes - 1) record << ", ";
		}
		record << endl;

		record << "  neighbours indexes: " << endl;
		for (int i = 0; i < nFaces; i++) {
			record << cells[i];
			if (i < nFaces - 1) record << ", ";
		}
		record << endl;

		record << "************************ " << endl;
	}
	else
		cout << "�� ������� ������� ���� " << f << endl;

	//record << "  square= " << S << endl;
	//record << "  square= " << S << endl;


	record.close();



}

void Cell::Set_Nodes(int* nodes_, int nn)
{
	// nn - ����� �����
	// nodes_ - ������ �������� �����, ������������� ������ ������� �������.
	nNodes = nn;

	nodes = new int[nNodes];    // ������ �����, ���������� ������ (nNodes)

	for (int i = 0; i < nn; i++) {
		nodes[i] = nodes_[i];
	}

}

