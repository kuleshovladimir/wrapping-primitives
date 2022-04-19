#include "Cell.h"

Cell::Cell()
{
	MassC.x = 0;
	MassC.y = 0;
	S = 0;   

	nFaces = 0; 
	nNodes = 0; 


	faces = new int[nFaces];
	nodes = new int[nNodes];

	fType = new bool [nFaces];

	cells = new int[nFaces];
}

Cell::~Cell()
{

}

Pnt Cell::Get_MassC()
{
	return MassC;
}

void Cell::Set_MassC(Pnt c_)
{
	MassC.x = c_.x;
	MassC.y = c_.y;
}

void Cell::Set_S(double S_)
{
	S = S_;
}

double Cell::Get_S()
{
	return S;
}

void Cell::Set_NFaces(int NF)
{
	nFaces = NF;

	faces = new int[nFaces];
	fType = new bool[nFaces];
	cells = new int[nFaces];
}

int Cell::Get_NFaces()
{
	return nFaces;
}

void Cell::Set_NNodes(int NN)
{
	nNodes = NN;

	nodes = new int[NN];
}

int Cell::Get_NNodes()
{
	return nNodes;
}

int Cell::Get_Node(int i)
{
	return nodes[i];
}

int Cell::Get_Face(int i)
{
	if (i < nFaces) {
		return faces[i];
	}
	else {
		cout << "Такого номера грани не существует" << endl;
		return (-1);
	}
}

void Cell::Set_Face(int iFace, int fIndex)
{
	faces[iFace] = fIndex;
}

void Cell::Set_Nodes(int* nodes_, int n)
{
	nNodes = n;
	nodes = new int[n]; //номера узлов, окружающих ячейку
	for (int i = 0; i < n; i++) {
		nodes[i] = nodes_[i];
	}
}

void Cell::Set_fType(int iFace, bool ftype)
{
	fType[iFace] = ftype;
}

void Cell::Set_cells(int iFace, int cell)
{
	cells[iFace] = cell;
}

int Cell::Get_Cell(int i)
{
	return cells[i];
}

void Cell::Print(int m)
{
	string f = "cells.txt";
	ofstream record(f, ios::out | ios::app);

	if (record) {

		record << "			cell index = " << m << endl;
		record << "center = " << MassC.x << ',' << MassC.y << endl;
		record << "square = " << S << endl;
		record << "number of faces = " << nFaces << endl;
		record << "face indexes and types: " << endl;
		for (int i = 0; i < nFaces; i++)
		{
			record << "ind = " << faces[i] << ", type = " << fType[i] << endl;
		}
		record << "number of nodes = " << nNodes << endl;
		record << "node indexes: " << endl;
		for (int i = 0; i < nNodes; i++) {
			record << nodes[i];
			if (i < nNodes - 1) record << ", ";
		}
		record << endl;

		record << "neighbour indexes: " << endl;
		for (int i = 0; i < nFaces; i++) {
			record << cells[i];
			if (i < nFaces - 1) record << ", ";
		}
		record << endl;
		record << "***********************" << endl;
	}
	else
		cout << "Не удалось открыть файл " << f << endl;

	record.close();
}

bool Cell::Get_fType(int i)
{
	return fType[i];
}

double Cell::Get_wk(int i)
{
	return wk[i];
}

Vector Cell::Get_ck(int i)
{
	return ck[i];
}
