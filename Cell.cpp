#include "Cell.h"
Cell::Cell()
{
	c.x = 0., c.y = 0.5;		// координаты ц.т.
	S = 0.;	// площадь €чейки

	nFaces = 0;		// число граней, окружающих €чейку

	nNodes = 0;		// число узлов, окружающих €чейку

	faces = new int[nFaces];		// номера граней, окружающих €чейку (nFaces)

	nodes = new int[nNodes];    // номера узлов, окружающих €чейку (nNodes)

	fType = new int[nFaces];	// типы граней:  =0 -> внутренн€€ грань	(nFaces)
								//   =1 -> гранична€ грань

	cells = new int[nFaces];		// номера соседних €чеек (nFaces)	
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

	faces = new int[nFaces];		// номера граней, окружающих €чейку (nFaces)

	fType = new int[nFaces];	// типы граней:  =0 -> внутренн€€ грань	(nFaces)
								//   =1 -> гранична€ грань

	cells = new int[nFaces];		// номера соседних €чеек (nFaces)

}

void Cell::Set_nNodes(int nn)
{
	nNodes = nn;

	nodes = new int[nNodes];    // номера узлов, окружающих €чейку (nNodes)

}

int Cell::Get_Face(int i)
{
	if (i < nFaces)
		return faces[i];
	else
	{
		cout << "“акого номер грани не существует" << endl;
		return -10;
	}
}

void Cell::Print(int m)
{
	// создаем поток дл€ записи
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
		cout << "Ќе удалось открыть файл " << f << endl;

	//record << "  square= " << S << endl;
	//record << "  square= " << S << endl;


	record.close();



}

void Cell::Set_Nodes(int* nodes_, int nn)
{
	// nn - число узлов
	// nodes_ - массив индексов узлов, перечисл€емых против часовой стрелки.
	nNodes = nn;

	nodes = new int[nNodes];    // номера узлов, окружающих €чейку (nNodes)

	for (int i = 0; i < nn; i++) {
		nodes[i] = nodes_[i];
	}

}

