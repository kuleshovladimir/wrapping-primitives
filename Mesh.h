#pragma once
#include "Cell.h"

struct parameters {
	// main
	double ro, p, h, H, u, v, w, T, E, e;

	// transport
	double mu, la, Pr;

	double Cp, Gam, Gm;

	// vectors
	double* U, * U1;
	double* V;   // ro, u,v, h


};

struct changes {
	double* dU;
};

struct Gradient {
	Vector* g;
};

struct Face {
	int nodes[2];		// ������� ����� �����	(2)
	int cr;				// ������ ������ ������
	int cl;				// ������ ����� ������
	bool is_boundary;	// �������� �� ��� ����� ���������
	Pnt f_centr;		// ���������� ������ �����
	double length;		// ����� �����
	int zone;			// ����� ����, �� ������� ��������������� ��������� �������

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
	double value;	// �������� Tw ��� qw
};

struct Boundary {
	int* rules;
	//��������, ��� ������ :  2 elements
		//!vel     !0 - Free slip;  1 - No slip
		//!temp    !1 - Tw is set;  2 - qw is set < ->values
	//��� ������������� ����� : \
		//type : 1 - ������ L, U, T, p, angle < ->values(5)

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
	friend void SetGran(Mesh& mesh);
	friend void Viscous(parameters* p, changes* (&du), Mesh mesh, Cell* cells, double dt, Gradient* gr, int Nm);
private:

	int Nx, Ny;		// ����� ����� �� x,y � ����������������� �����
	int nNodes, nFaces, nCells;
	Pnt* nodes;		// ���������� �����

	Face* faces;	// ��� ����� �����		***

	int nZones;
	Zone* zones;
public:
	Mesh();
	~Mesh();

	void ReadStruct(string filename);

	void CreateCell(Cell* (&cells));

	void CreateFaces();

	void CellFuncs(Cell* (&cells));

	void SetZones();

	void GradCoeffs(Cell* (&cells));

	int Get_Nx();
	int Get_Ny();
	int Get_nFaces();
	int Get_nNodes();
	int Get_nCells();
	int Get_nZones();

	Pnt Get_node(int n);
	Face Get_face(int n);
	Zone Get_zone(int n);

};

