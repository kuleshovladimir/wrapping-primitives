#include "Mesh.h"

Mesh::Mesh()
{
	Nx = 0;
	Ny = 0;		// число точек по x,y в структурированной сетке
	nNodes = 0;
	nCells = 0;
	nFaces = 0;
	nodes = new Pnt[nNodes];		// координаты узлов
}

Mesh::~Mesh()
{
}

void Mesh::ReadStruct(string filename)
{
	int z;

	ifstream reading(filename);

	if (reading) {

		reading >> Nx >> Ny >> z >> z;
		nNodes = Nx * Ny;
		nCells = (Nx - 1) * (Ny - 1);
		nFaces = Nx * (Ny - 1) + (Nx - 1) * Ny;

		nodes = new Pnt[nNodes];		// координаты узлов

		int k = 0;
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				reading >> z >> z >> nodes[k].x >> nodes[k].y;
				k++;
			}
		}

		k = 9;
		cout << nodes[k].x << ", " << nodes[k].y << endl;



		cout << "ƒанные считаны из файла " << filename << endl;
	}
	else
		cout << "Ќе удалось открыть файл " << filename << endl;

}

void Mesh::CreateCell(Cell* (&cells))
{
	int k = 0;		// счетчик граней
	for (int i = 0; i < Nx - 1; i++) {
		for (int j = 0; j < Ny - 1; j++) {

			int nc = (Ny - 1) * i + j;		// номер €чейки

			//cout << " nc= " << nc << endl;

			int nNodes = 4;					// число узлов вокруг €чейки

			// создаем массив индексов узлов, расположенных вокруг €чейки, проход€ узлы против часовой стрелки
			int* nodes = new int[nNodes];

			// элемент под номером 0
			int n2 = Ny * i + j;
			nodes[0] = n2;

			// элемент под номером 1
			n2 = Ny * (i + 1) + j;
			nodes[1] = n2;

			// элемент под номером 2
			n2 = Ny * (i + 1) + (j + 1);
			nodes[2] = n2;

			// элемент под номером 3
			n2 = Ny * i + (j + 1);
			nodes[3] = n2;

			cells[nc].Set_Nodes(nodes, nNodes);


		}

	}


}

void Mesh::CreateFaces()
{
	faces = new Face[nFaces];

	int k = 0;		// счетчик граней
	// вертикальные грани
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

			faces[k].f_centr.x = 0.5 * (nodes[n1].x + nodes[n2].x);
			faces[k].f_centr.y = 0.5 * (nodes[n1].y + nodes[n2].y);

			double dx = (nodes[n1].x - nodes[n2].x);
			double dy = (nodes[n1].y - nodes[n2].y);
			faces[k].length = sqrt(dx * dx + dy * dy);

			k++;
		}
	}
	//cout << "k= " << k << endl;

	// горизонтальные грани
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

			faces[k].f_centr.x = 0.5 * (nodes[n1].x + nodes[n2].x);
			faces[k].f_centr.y = 0.5 * (nodes[n1].y + nodes[n2].y);

			double dx = (nodes[n1].x - nodes[n2].x);
			double dy = (nodes[n1].y - nodes[n2].y);
			faces[k].length = sqrt(dx * dx + dy * dy);

			k++;
		}
	}


}

void Mesh::CellFuncs(Cell* (&cells))
{
	for (int i = 0; i < nCells; i++) {

		// массив точек - узлов вокруг €чейки
		int mNodes = cells[i].Get_nNodes();

		//cout << "mNodes= " << mNodes << endl;

		int* nds = new int[mNodes];
		Pnt* pnts = new Pnt[mNodes];
		for (int j = 0; j < mNodes; j++) {
			int n_ = cells[i].Get_Node(j);
			nds[j] = n_;
			//cout << "nds[j]= " << nds[j] << endl;

			//координаты узлов
			pnts[j].x = nodes[n_].x;
			pnts[j].y = nodes[n_].y;


		}

		//создаем полигон
		Polygons pl(pnts, mNodes);

		Pnt center = pl.MassCenter();

		double Sq = pl.Square();

		cells[i].Set_c(center);
		cells[i].Set_S(Sq);

		cells[i].Set_nFaces(mNodes);
	}

	//return;

	for (int k = 0; k < nFaces; k++) {

		int n1 = faces[k].nodes[0];
		int n2 = faces[k].nodes[1];
		int cl = faces[k].cl;
		int cr = faces[k].cr;
		int ftype = 0;
		if (faces[k].is_boundary) ftype = 1;

		if (cr >= 0) {
			// массив узлов €чейки cl
			int nn = cells[cr].Get_nNodes();  // размер массива узлов
			int* nncl = new int[nn];			// выделение массива под эти индексы

			for (int i = 0; i < nn; i++) {
				nncl[i] = cells[cr].Get_Node(i);			
			}

			int iFace;

			for (int i = 0; i < nn; i++) {
				int j = i + 1;
				if (j == nn) j = 0;
				if (nncl[i] == n1 && nncl[j] == n2) {
					iFace = i;
					cells[cr].Set_Face(iFace, k);
					cells[cr].Set_fType(iFace, ftype);
					cells[cr].Set_cells(iFace, cl);
				}
			}
		}

		if (cl >= 0) {
			// массив узлов €чейки cl
			int nn = cells[cl].Get_nNodes();  // размер массива узлов
			int* nncl = new int[nn];			// выделение массива под эти индексы

			//cout << " n1= " << n1 << " n2= " << n2 << endl;

			//cout << " cl: " << cl << endl;

			for (int i = 0; i < nn; i++) {
				nncl[i] = cells[cl].Get_Node(i);			//int Get_Node(int i) { return nodes[i]; };

				//cout << i << " nodes: " << nncl[i] << endl;
			}

			int iFace;
			for (int i = 0; i < nn; i++) {
				int j = i + 1;
				if (j == nn) j = 0;
				if (nncl[i] == n2 && nncl[j] == n1) {
					iFace = i;
					cells[cl].Set_Face(iFace, k);
					cells[cl].Set_fType(iFace, ftype);
					cells[cl].Set_cells(iFace, cr);
				}
			}



			//if (cl==7) exit(-6);

		}



	}


}

void Mesh::SetZones()
{
	nZones = 4;
	zones = new Zone[nZones];
}

void Mesh::GradCoeffs(Cell* (&cells))
{
	for (int i = 0; i < nCells; i++) {

		//cout << "i= " << i << endl;

		int nFaces = cells[i].Get_nFaces(); // число граней, окружающих €чейку i

		// создаем массив весовых коэффициентов wk и векторов ck
		cells[i].wk = new double[nFaces];
		cells[i].ck = new Vector[nFaces];

		int iDim = 2;

		for (int k = 0; k < nFaces; k++) {
			cells[i].ck[k].cx = new double[iDim];
		}

		// координаты центра данной €чейки
		Pnt xc = cells[i].Get_c();

		// рассто€ни€ до соседей
		double* dx = new double[nFaces];
		double* dy = new double[nFaces];

		// компоненты матрицы
		double axx = 0., axy = 0, ayy = 0.;

		//cout << "before: for (int k " << endl;

		for (int k = 0; k < nFaces; k++) {
			int nf = cells[i].Get_Face(k); // номер грани
			if (!faces[nf].is_boundary) {     // internal face

				// номер соседней €чейки
				int nc = cells[i].Get_Cell(k);
				// координаты ее центра
				Pnt xk = cells[nc].Get_c();

				dx[k] = xk.x - xc.x;
				dy[k] = xk.y - xc.y;



			}
			else {
				// координаты центра грани
				Pnt xk = faces[nf].f_centr;
				dx[k] = xk.x - xc.x;
				dy[k] = xk.y - xc.y;

			}
			//weight coefficients
			double wk = 1. / sqrt(dx[k] * dx[k] + dy[k] * dy[k]);
			cells[i].wk[k] = wk;

			axx += wk * dx[k] * dx[k];
			axy += wk * dx[k] * dy[k];
			ayy += wk * dy[k] * dy[k];

		}

		// ќбратна€ матрица M
		double det = axx * ayy - axy * axy;
		double Mxx, Mxy, Myy;
		//if (det != 0) {
		Mxx = ayy / det;
		Mxy = -axy / det;
		Myy = axx / det;

	// ¬ектор ck
		for (int k = 0; k < nFaces; k++) {

			cells[i].ck[k].cx[0] = Mxx * dx[k] + Mxy * dy[k];
			cells[i].ck[k].cx[1] = Mxy * dx[k] + Myy * dy[k];
		}


	}


}

int Mesh::Get_Nx()
{
	return Nx;
}

int Mesh::Get_Ny()
{
	return Ny;
}

int Mesh::Get_nFaces()
{
	return nFaces;
}

int Mesh::Get_nNodes()
{
	return nNodes;
}

int Mesh::Get_nCells()
{
	return nCells;
}

int Mesh::Get_nZones()
{
	return nZones;
}

Pnt Mesh::Get_node(int n)
{
	return nodes[n];
}

Face Mesh::Get_face(int n)
{
	return faces[n];
}

Zone Mesh::Get_zone(int n)
{
	return zones[n];
}
