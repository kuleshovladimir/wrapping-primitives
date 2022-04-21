#include "Mesh.h"

Mesh::Mesh()
{
	Nx = 0;
	Ny = 0;
	nNodes = 0;
	nCells = 0;
	nFaces = 0;
	nodes = new Pnt[nNodes];
	nZones = 4;
	zones = new Zone[nZones];
}

Mesh::~Mesh()
{
}

void Mesh::SetZones(int n)
{
	nZones = n;
	zones = new Zone[nZones];
}

Zone Mesh::Get_zone(int i)
{
	return zones[i];
}

void Mesh::Set_Nx(int n)
{
	Nx = n;
}

int Mesh::Get_Nx()
{
	return Nx;
}

void Mesh::Set_Ny(int n)
{
	Ny = n;
}

int Mesh::Get_Ny()
{
	return Ny;
}

void Mesh::Set_nNodes(int n)
{
	nNodes = n;
}

int Mesh::Get_nNodes()
{
	return nNodes;
}

void Mesh::Set_nCells(int n)
{
	nCells = n;
}

int Mesh::Get_nCells()
{
	return nCells;
}

void Mesh::Set_nFaces(int n)
{
	nFaces = n;
}

int Mesh::Get_nFaces()
{
	return nFaces;
}

Pnt Mesh::Get_node(int i)
{
	return nodes[i];
}

void Mesh::Set_node(Pnt node, int i)
{
	nodes[i] = node;
}

Face Mesh::Get_face(int i)
{
	return faces[i];
}

void Mesh::ReadStruct(string filename)
{
	double tmp;

	ifstream reading(filename);
	if (reading) {
		reading >> Nx >> Ny >> tmp >> tmp;
		nNodes = Nx * Ny;
		nCells = (Nx - 1) * (Ny - 1);
		nFaces = Nx * (Ny - 1) + Ny * (Nx - 1);

		nodes = new Pnt[nNodes];

		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				reading >> tmp >> tmp >> nodes[Ny * i + j].x >> nodes[Ny * i + j].y;
			}
		}

		cout << "ƒанные считаны из файла " << filename << endl;
	}
	else {
		cout << "Ќе удалось открыть файл " << filename << endl;
	}

	

}

void Mesh::CreateCells(Cell* (&cells))
{
	for (int i = 0; i < Nx - 1; i++) {
		for (int j = 0; j < Ny - 1; j++) {
			int* nodes = new int[4];

			nodes[0] = Ny * i + j;
			nodes[1] = Ny * (i + 1) + j;
			nodes[2] = Ny * (i + 1) + (j + 1);
			nodes[3] = Ny * i + (j + 1);

			cells[(Ny - 1) * i + j].Set_Nodes(nodes, 4);
		}
	}

}

void Mesh::CreateFases()
{
	faces = new Face[nFaces];
	
	int k = 0;
	//вертикальные грани
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny - 1; j++) {
			int n1 = Ny * i + j;
			int n2 = Ny * i + j + 1;

			faces[k].nodes[0] = n1;
			faces[k].nodes[1] = n2;
			if (i == 0) {
				faces[k].is_boundary = true;
				faces[k].cr = -1;
				faces[k].cl = (Ny - 1) * i + j;
				faces[k].zone = 0;
			}
			else if (i == Nx - 1) {
				faces[k].is_boundary = true;
				faces[k].cr = (Ny - 1) * (i - 1) + j;
				faces[k].cl = -1;
				faces[k].zone = 2;
 			}
			else {
				faces[k].is_boundary = false;
				faces[k].cr = (Ny - 1) * (i - 1) + j;
				faces[k].cl = (Ny - 1) * i + j;
				faces[k].zone = -1;
			}

			faces[k].f_center.x = 0.5 * (nodes[n1].x + nodes[n2].x);
			faces[k].f_center.y = 0.5 * (nodes[n1].y + nodes[n2].y);

			double dx = nodes[n1].x - nodes[n2].x;
			double dy = nodes[n1].y - nodes[n2].y;
			faces[k].length = sqrt(dx * dx + dy * dy);

			//faces[k].Print(k);

			k++;
		}
	}
	//горизонтальные грани
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx - 1; i++) {
			int n1 = Ny * i + j;
			int n2 = Ny * (i + 1) + j;

			faces[k].nodes[0] = n1;
			faces[k].nodes[1] = n2;
			if (j == 0) {
				faces[k].is_boundary = true;
				faces[k].cl = -1;
				faces[k].cr = (Ny - 1) * i + j;
				faces[k].zone = 1;
			}
			else if (j == Ny - 1) {
				faces[k].is_boundary = true;
				faces[k].cl = (Ny - 1) * i + j - 1;
				faces[k].cr = -1;
				faces[k].zone = 3;
			}
			else {
				faces[k].is_boundary = false;
				faces[k].cr = (Ny - 1) * i + j;
				faces[k].cl = (Ny - 1) * i + j - 1;
				faces[k].zone = -1;
			}

			faces[k].f_center.x = 0.5 * (nodes[n1].x + nodes[n2].x);
			faces[k].f_center.y = 0.5 * (nodes[n1].y + nodes[n2].y);

			double dx = nodes[n1].x - nodes[n2].x;
			double dy = nodes[n1].y - nodes[n2].y;
			faces[k].length = sqrt(dx * dx + dy * dy);

			//faces[k].Print(k);

			k++;
		}
	}
}

void Mesh::CellFuncs(Cell* (&cells))
{
	for (int i = 0; i < nCells; i++) {

		//массив узлов вокруг €чейки
		int m = cells[i].Get_NNodes();

		int* nds = new int[m];
		Pnt* pnts = new Pnt[m];
		for (int j = 0; j < m; j++) {
			int n_ = cells[i].Get_Node(j);
			nds[j] = n_;

			//координаты узлов
			pnts[j].x = nodes[n_].x;
			pnts[j].y = nodes[n_].y;
 		}

		Polygons pl(pnts, m);

		Pnt center = pl.MassCenter();

		double S = pl.Square();
		
		cells[i].Set_MassC(center);
		cells[i].Set_S(S);

		cells[i].Set_NFaces(m);
	}

	for (int k = 0; k < nFaces; k++) {
		int n1 = faces[k].nodes[0];
		int n2 = faces[k].nodes[1];
		int cl = faces[k].cl;
		int cr = faces[k].cr;
		bool ftype = false;
		if (faces[k].is_boundary) ftype = true;
		
		if (cr >= 0) {
			//массив узлов €чейки c1
			int n = cells[cr].Get_NNodes(); //размер массива узлов
			int* nc1 = new int[n];			//масси индексов

			for (int i = 0; i < n; i++) {
				nc1[i] = cells[cr].Get_Node(i);
			}

			int iFace;

			for (int i = 0; i < n; i++) {
				int j = i + 1;
				if (j == n) j = 0;
				if ((nc1[i] == n1) && (nc1[j] == n2)) {
					iFace = i;
					cells[cr].Set_Face(iFace, k);
					cells[cr].Set_fType(iFace, ftype);
					cells[cr].Set_cells(iFace, cl);
				}
			}
		}

		if (cl >= 0) {
			//массив узлов €чейки c1
			int n = cells[cl].Get_NNodes(); //размер массива узлов
			int* nc1 = new int[n];			//масси индексов

			for (int i = 0; i < n; i++) {
				nc1[i] = cells[cl].Get_Node(i);
			}

			int iFace;

			for (int i = 0; i < n; i++) {
				int j = i + 1;
				if (j == n) j = 0;
				if ((nc1[i] == n2) && (nc1[j] == n1)) {
					iFace = i;
					cells[cl].Set_Face(iFace, k);
					cells[cl].Set_fType(iFace, ftype);
					cells[cl].Set_cells(iFace, cr);
				}
			}
		}
	}
}

int Mesh::Get_nZones()
{
	return nZones;
}

void Mesh::GradCoeffs(Cell* (&cells))
{
	for (int i = 0; i < nCells; i++) {
		int nFaces = cells[i].Get_NFaces();

		//массивы весовых коэффициентов wk и векторов ck
		cells[i].wk = new double[nFaces];
		cells[i].ck = new Vector[nFaces];

		int iDim = 2;

		for (int j = 0; j < nFaces; j++) {
			cells[i].ck[j].cx = new double[iDim];
		}

		//координаты центра €чейки
		Pnt xc = cells[i].Get_MassC();

		//рассто€ни€ до соседей
		double* dx = new double[nFaces];
		double* dy = new double[nFaces];

		//компоненты матрицы
		double axx = 0, axy = 0, ayy = 0;

		for (int k = 0; k < nFaces; k++) {
			int nf = cells[i].Get_Face(k);
			if (!faces[nf].is_boundary) {
				//номер соседней €чейки
				int nc = cells[i].Get_Cell(k);
				//ее центр
				Pnt xk = cells[nc].Get_MassC();

				dx[k] = xk.x - xc.x;
				dy[k] = xk.y - xc.y;

			}
			else{
				//центр грани
				Pnt xk = faces[nf].f_center;
				dx[k] = xk.x - xc.x;
				dy[k] = xk.y - xc.y;
			}

			double wk = 1 / sqrt(dx[k] * dx[k] + dy[k] * dy[k]);
			cells[i].wk[k] = wk;

			axx += wk * dx[k] * dx[k];
			axy += wk * dx[k] * dy[k];
			ayy += wk * dy[k] * dy[k];
		}

		//получение обратной матрицы
		double det = axx * ayy - axy * axy;
		double Mxx, Mxy, Myy;

		Mxx = ayy / det;
		Mxy = -axy / det;
		Myy = axx / det;

		for (int k = 0; k < nFaces; k++) {
			cells[i].ck[k].cx[0] = Mxx * dx[k] + Mxy * dy[k];
			cells[i].ck[k].cx[1] = Mxy * dx[k] + Myy * dy[k];
		}

	}
}

