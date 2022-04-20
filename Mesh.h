#pragma once
#include "Cell.h"


struct  Face
{
	int nodes[2];		//индексы узлов грани
	int cr;				//индекс правой ячейки
	int cl;				//индекс левой ячейки
	bool is_boundary;   //является ли эта грань граничной
	Pnt f_center;		//координаты грани центра
	double length;		//длина грани
	int zone;			//номер зоны, на которой устанавливаются граничные условия


	void Print(int k) {
		cout << "			face index = " << k << endl;
		cout << "nodes indexes = " << nodes[0] << " " << nodes[1] << endl;
		cout << "cr index = " << cr << ", cl index = " << cl << endl;
		cout << "center = " << f_center.x << " " << f_center.y << endl;
		cout << "length = " << length << endl;
	}
};

struct Wall {
	bool vel;        //0 - проскальзывание, 1 - без проскальзывания
	int temp;       //1 - задана температура, 2 - задан тепловой поток
	double value;   //значение Tw или qw
};

struct  Boudary
{
	int* rules; //определяет тип стенки и управляющие параметры
			
	// стенка : vel, temp
	// сверхзвуковой вход : type - L(размер), U(скорость), T(температура), p(давление), angle(угол потока)
	
	double* vals;

};

struct Zone
{
	int granType; //1 - стенка, 2 - сверхзвуковой вход, 3 - симметрия, 4 - сверзвуковой выход,
				  //5 - свободная границв, 6 - дозвуковой вход, 7 - дозвуковой выход
	Wall* wall;
	Boudary* bnd;
};

class Mesh
{
	friend void SetGran(Mesh& mesh);
private:
	int Nx, Ny;			//число точек по x, y в структурированной сетке
	int nNodes, nCells, nFaces;
	Pnt* nodes;			//координаты узлов

	Face* faces;		//все грани сетки

	int nZones;
	Zone* zones;

public:
	Mesh();
	~Mesh();

	void SetZones(int n);
	Zone Get_zone(int i);

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

	Face Get_face(int i);

	void ReadStruct(string filename);

	void CreateCells(Cell* (&cells));

	void CreateFases();

	void CellFuncs(Cell* (&cells));

	int Get_nZones();

	void GradCoeffs(Cell* (&cells));
};

