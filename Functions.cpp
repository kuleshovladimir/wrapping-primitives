#include "Functions.h"

void Init(parameters* (&p), int nCells, int Nm)
{

	double T0, Cp, la, P0, Gm, Gam, R, ro, mu;

	double U0;

	//mesh.zones[3].bnd[0].vals[1] = 867.9; // [ m/s ] - U
	//mesh.zones[3].bnd[0].vals[2] = 75.1; // [ K ] - T
	//mesh.zones[3].bnd[0].vals[3] = 1931.; // [ Pa ] - p

	// входные данные
	T0 = 275.1;		// K
	la = 2.5658e-2;	// W/(m K)
	mu = 17.863e-6;

	U0 = 2000.;

	P0 = 1931.;	// Pa
	Gam = 1.4;
	Gm = 28.97;

	// расчет
	R = 8314.41 / Gm;  // J/(kg K)
	ro = P0 / (R * T0);
	Cp = Gam / (Gam - 1.) * R;

	//cout << "Cp= " << Cp << endl;

	//int n = nCells;

	for (int i = 0; i < nCells; i++) {   // ro, p, h, H, u, v, w, T, E;

		p[i].ro = ro;
		p[i].p = P0;
		p[i].u = U0;
		p[i].v = 0.;
		p[i].w = 0.;
		p[i].T = T0;

		p[i].Cp = Cp;
		p[i].la = la;
		p[i].mu = mu;
		p[i].Gam = Gam;
		p[i].Gm = Gm;

		p[i].Pr = p[i].mu * p[i].Cp / p[i].la;

		p[i].h = Cp * T0;

		double q2 = 0.5 * (p[i].u * p[i].u + p[i].v * p[i].v);

		p[i].H = p[i].h + q2;
		p[i].E = p[i].H - p[i].p / p[i].ro;

		p[i].e = p[i].h / p[i].Gam;

		// NS :
		// &&&&&&&&&&&&&&&&&&&

		p[i].U = new double[Nm];
		p[i].U1 = new double[Nm];

		p[i].V = new double[Nm];

		p[i].U1[0] = p[i].ro;
		p[i].U1[1] = p[i].ro * p[i].u;
		p[i].U1[2] = p[i].ro * p[i].v;
		p[i].U1[3] = p[i].ro * p[i].E;

		p[i].V[0] = p[i].ro;
		p[i].V[1] = p[i].u;
		p[i].V[2] = p[i].v;
		p[i].V[3] = p[i].h;


	}

}

//void Viscous(parameters* p, changes* (&du), Mesh mesh, Cell* cells, double dt)
//{
//	int nFaces = mesh.Get_nFaces();
//
//	for (int i = 0; i < nFaces; i++) {
//
//		Face face = mesh.faces[i];
//
//		int cr = face.cr;
//		int cl = face.cl;
//
//		if (face.is_boundary) {
//			int c = max(cr, cl);
//
//			// координаты ц.т. ячейки
//			Pnt xc = cells[c].Get_c();
//
//			// координаты центра грани
//			Pnt xf = face.f_centr;
//
//			// расстояние до стенки
//			double dx = (xc.x - xf.x);
//			double dy = (xc.y - xf.y);
//			double dl = sqrt(dx * dx + dy * dy);
//
//			double length = face.length;		// длина грани
//			double S = cells[c].Get_S();		// площадь ячейки с
//
//			int z = face.zone;
//
//			// boundary type:
//			int btype = mesh.zones[z].grantype;
//
//
//			double Tw;
//			if (face.zone == 1) Tw = 500.;
//
//			if (face.zone == 3) Tw = 200.;
//
//			if (face.zone == 1 || face.zone == 3) {
//				double hw = p[c].Cp * Tw;
//
//				// dh/dn
//				double dh_dn = (p[c].h - hw) / dl;
//
//				// mu/Pr
//				double mu_Pr = p[c].mu / p[c].Pr;
//
//				// Fv - поток через грань
//				double Fv = mu_Pr * dh_dn;
//
//				du[c].dU[0] += -Fv * length / S * dt;
//				//du[c].dU[0] = du[c].dU[0] + ( -Fv * length / S * dt );
//
//			}
//
//			if (face.zone == 0) {
//				// Fv - поток через грань
//				double Fv = 0.;
//				du[c].dU[0] += Fv * length / S * dt;
//			}
//
//			if (face.zone == 2) {
//				// Fv - поток через грань
//				double Fv = 0.;
//				du[c].dU[0] += Fv * length / S * dt;
//			}
//
//
//
//
//
//		}
//		else
//		{
//
//			// координаты правой и левой ячеек
//			Pnt xr = cells[cr].Get_c();
//			Pnt xl = cells[cl].Get_c();
//
//			// расстояние между ячейками
//			double dx = (xr.x - xl.x);
//			double dy = (xr.y - xl.y);
//			double dl = sqrt(dx * dx + dy * dy);
//
//			// dh/dn
//			double dh_dn = (p[cr].h - p[cl].h) / dl;
//
//			// среднее mu/Pr
//			double mu_Pr = 0.5 * (p[cr].mu / p[cr].Pr + p[cl].mu / p[cl].Pr);
//
//			// Fv - поток через грань
//			double Fv = mu_Pr * dh_dn;
//
//			// приращения
//			// если p[cr].h > p[cl].h => dh_dn>0 => Fv>0
//			// тепло втекает в левую ячейку => для нее в приращении Fv берется с плюсом
//			// для правой -  с минусом
//
//			double length = face.length;		// длина грани
//			double Sr = cells[cr].Get_S();		// площадь правой ячейки
//			double Sl = cells[cl].Get_S();		// площадь левой ячейки
//
//			du[cr].dU[0] += -Fv * length / Sr * dt;
//			du[cl].dU[0] += Fv * length / Sl * dt;
//
//
//		}
//
//
//	}
//
//}

void Viscous(parameters* p, changes* (&du), Mesh mesh, Cell* cells, double dt, Gradient* gr, int Nm)
{
	int nFaces = mesh.Get_nFaces();

	double* Fv = new double[Nm];

	for (int i = 0; i < nFaces; i++) {

		Face face = mesh.faces[i];

		int cr = face.cr;
		int cl = face.cl;

		double length = face.length; // длина грани

		int n1 = face.nodes[0];
		int n2 = face.nodes[1];

		Pnt x1 = mesh.nodes[n1];
		Pnt x2 = mesh.nodes[n2];

		double nx = -(x2.y - x1.y) / length;
		double ny = (x2.x - x1.x) / length;

		if (face.is_boundary) {
			int c = max(cr, cl);

			// координаты ц.т. ячейки
			Pnt xc = cells[c].Get_c();

			int z = face.zone;
			int grantype = mesh.zones[z].grantype;

			if (grantype == 1) {

				// расстояние до стенки
				double dl = cells[c].Yw;

				double du_dx, du_dy, dv_dx, dv_dy, dh_dx, dh_dy;
				double u_, v_;

				// velocity
				if (mesh.zones[z].bnd[0].rules[0] = 0) {	// Free slip
					// grad	u:
					double ux = gr[c].g[1].cx[0];
					double uy = gr[c].g[1].cx[1];
					du_dx = ny * ny * ux - nx * ny * uy;
					du_dy = -nx * ny * ux + nx * nx * uy;

					double vx = gr[c].g[2].cx[0];
					double vy = gr[c].g[2].cx[1];
					dv_dx = ny * ny * vx - nx * ny * vy;
					dv_dy = -nx * ny * vx + nx * nx * vy;

					u_ = p[c].u;
					v_ = p[c].v;

				}

				if (mesh.zones[z].bnd[0].rules[0] = 1) {	// No slip
					double du_dn = p[c].u / dl;
					double dv_dn = p[c].v / dl;
					double du_dl = 0.;
					double dv_dl = 0.;
					du_dx = du_dn * nx - du_dl * ny;
					dv_dx = dv_dn * nx - dv_dl * ny;

					du_dy = du_dn * ny + du_dl * nx;
					dv_dy = dv_dn * ny + dv_dl * nx;

					u_ = 0;
					v_ = 0;

				}

				// energy
				if (mesh.zones[z].bnd[0].rules[1] = 1) {	// Tw

					double Tw = mesh.zones[z].bnd[0].vals[0];
					double hw = p[c].Cp * Tw;
					double dh_dn = (p[c].h - hw) / dl;

					double dh_dl = 0.;

					dh_dx = dh_dn * nx - dh_dl * ny;
					dh_dy = dh_dn * ny + dh_dl * nx;


				}

				if (mesh.zones[z].bnd[0].rules[1] = 2) {	// qw = 0

					double tx = gr[c].g[3].cx[0];
					double ty = gr[c].g[3].cx[1];

					dh_dx = ny * ny * tx - nx * ny * ty;
					dh_dy = -nx * ny * tx + nx * nx * ty;


				}

				double mu = p[c].mu;
				double div = du_dx + dv_dy;
				double txx = mu * (2. * du_dx - 2. / 3. * div);
				double tyy = mu * (2. * dv_dy - 2. / 3. * div);
				double txy = mu * (du_dy + dv_dx);

				// mu/Pr
				double mu_Pr = p[c].mu / p[c].Pr;
				double qx = -mu_Pr * dh_dx;
				double qy = -mu_Pr * dh_dy;


				Fv[0] = 0.;
				Fv[1] = txx * nx + txy * ny;
				Fv[2] = txy * nx + tyy * ny;
				Fv[3] = (u_ * txx + v_ * txy - qx) * nx + (u_ * txy + v_ * tyy - qy) * ny;


			}  //if (grantype == 1) {

			if (grantype == 3) {	// symmetry

				double du_dx, du_dy, dv_dx, dv_dy, dh_dx, dh_dy;
				double u_, v_;

				// velocity

					// grad	u:
				double ux = gr[c].g[1].cx[0];
				double uy = gr[c].g[1].cx[1];
				du_dx = ny * ny * ux - nx * ny * uy;
				du_dy = -nx * ny * ux + nx * nx * uy;

				double vx = gr[c].g[2].cx[0];
				double vy = gr[c].g[2].cx[1];
				dv_dx = ny * ny * vx - nx * ny * vy;
				dv_dy = -nx * ny * vx + nx * nx * vy;

				u_ = p[c].u;
				v_ = p[c].v;

				// energy

				double tx = gr[c].g[3].cx[0];
				double ty = gr[c].g[3].cx[1];

				dh_dx = ny * ny * tx - nx * ny * ty;
				dh_dy = -nx * ny * tx + nx * nx * ty;


				double mu = p[c].mu;
				double div = du_dx + dv_dy;
				double txx = mu * (2. * du_dx - 2. / 3. * div);
				double tyy = mu * (2. * dv_dy - 2. / 3. * div);
				double txy = mu * (du_dy + dv_dx);

				// mu/Pr
				double mu_Pr = p[c].mu / p[c].Pr;
				double qx = -mu_Pr * dh_dx;
				double qy = -mu_Pr * dh_dy;


				Fv[0] = 0.;
				Fv[1] = txx * nx + txy * ny;
				Fv[2] = txy * nx + tyy * ny;
				Fv[3] = (u_ * txx + v_ * txy - qx) * nx + (u_ * txy + v_ * tyy - qy) * ny;

			} // if (grantype == 3) {	// symmetry

			if (grantype == 2 || grantype == 4) {	// Supersonic Inlet, Outlet

				for (int m = 0; m < 4; m++) {
					Fv[m] = 0.;
				}

			}

			double length = face.length;		// длина грани
			if (cr >= 0) {
				double Sr = cells[cr].Get_S();		// площадь правой ячейки
				for (int m = 0; m < Nm; m++) {
					du[cr].dU[m] += -Fv[m] * length / Sr * dt;
				}
			}

			if (cl >= 0) {
				double Sl = cells[cl].Get_S();		// площадь левой ячейки
				for (int m = 0; m < Nm; m++) {
					du[cl].dU[m] += +Fv[m] * length / Sl * dt;
				}
			}

		}
		else {		// internal face


			// координаты правой и левой ячеек
			Pnt xr = cells[cr].Get_c();
			Pnt xl = cells[cl].Get_c();

			// расстояние между ячейками
			double dx = (xr.x - xl.x);
			double dy = (xr.y - xl.y);
			double dl = sqrt(dx * dx + dy * dy);

			// градиенты в соседних ячейках
			//Vector GrR = gr[cr].g[0];
			//Vector GrL = gr[cl].g[0];
			//Vector Gr;
			//Gr.cx[0] = 0.5 * (gr[cr].g[0].cx[0] + gr[cl].g[0].cx[0]);

			double du_dx = 0.5 * (gr[cr].g[1].cx[0] + gr[cl].g[1].cx[0]);
			double du_dy = 0.5 * (gr[cr].g[1].cx[1] + gr[cl].g[1].cx[1]);

			double dv_dx = 0.5 * (gr[cr].g[2].cx[0] + gr[cl].g[2].cx[0]);
			double dv_dy = 0.5 * (gr[cr].g[2].cx[1] + gr[cl].g[2].cx[1]);

			double dh_dx = 0.5 * (gr[cr].g[3].cx[0] + gr[cl].g[3].cx[0]);
			double dh_dy = 0.5 * (gr[cr].g[3].cx[1] + gr[cl].g[3].cx[1]);

			double mu = 0.5 * (p[cr].mu + p[cl].mu);
			double div = du_dx + dv_dy;
			double txx = mu * (2. * du_dx - 2. / 3. * div);
			double tyy = mu * (2. * dv_dy - 2. / 3. * div);
			double txy = mu * (du_dy + dv_dx);

			// среднее mu/Pr
			double mu_Pr = 0.5 * (p[cr].mu / p[cr].Pr + p[cl].mu / p[cl].Pr);
			double qx = -mu_Pr * dh_dx;
			double qy = -mu_Pr * dh_dy;

			double u_ = 0.5 * (p[cr].u + p[cl].u);
			double v_ = 0.5 * (p[cr].v + p[cl].v);

			Fv[0] = 0.;
			Fv[1] = txx * nx + txy * ny;
			Fv[2] = txy * nx + tyy * ny;
			Fv[3] = (u_ * txx + v_ * txy - qx) * nx + (u_ * txy + v_ * tyy - qy) * ny;

			double length = face.length;		// длина грани
			double Sr = cells[cr].Get_S();		// площадь правой ячейки
			double Sl = cells[cl].Get_S();		// площадь левой ячейки

			for (int m = 0; m < Nm; m++) {
				du[cr].dU[m] += -Fv[m] * length / Sr * dt;
				du[cl].dU[m] += Fv[m] * length / Sl * dt;
			}


			//cout << i << ", " << cr << ", " << cl << endl;

		}

	}
}

void GetParams(parameters* (&p), int nCells, int Nm)
{

	for (int i = 0; i < nCells; i++) {

		p[i].ro = p[i].U1[0];
		p[i].u = p[i].U1[1] / p[i].ro;
		p[i].v = p[i].U1[2] / p[i].ro;
		p[i].E = p[i].U1[3] / p[i].ro;

		double q2 = 0.5 * (p[i].u * p[i].u + p[i].v * p[i].v);

		p[i].e = p[i].E - q2;

		p[i].h = p[i].e * p[i].Gam;

		p[i].T = p[i].h / p[i].Cp;

		p[i].H = p[i].h + q2;

		double R = 8314.41 / p[i].Gm;
		p[i].p = p[i].ro * R * p[i].T;

		p[i].V[0] = p[i].ro;
		p[i].V[1] = p[i].u;
		p[i].V[2] = p[i].v;
		p[i].V[3] = p[i].h;

	}

}

void Tecplot(parameters* p, Cell* cells, int Nx, int Ny, int nCells)
{
	// создаем поток для записи
	string f = "T.plt";
	ofstream record(f, ios::out);
	if (record) {
		record << "VARIABLES = \"X\", \"Y\", \"T,K\" " << endl;

		record << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;
		for (int i = 0; i < nCells; i++) {
			record << cells[i].Get_c().x << " " << cells[i].Get_c().y << " " << p[i].T << endl;
		}

	}
	record.close();
}

void Convect(parameters* p, changes* (&du), Mesh mesh, Cell* cells, int It, double dt)
{
	int nFaces = mesh.Get_nFaces();

	for (int i = 0; i < nFaces; i++) {

		Face face = mesh.Get_face(i);

		int cr = face.cr;
		int cl = face.cl;

		// nodes:
		int n1 = face.nodes[0];
		int n2 = face.nodes[1];

		double dl = face.length; // длина грани

		Pnt x1 = mesh.Get_node(n1);
		Pnt x2 = mesh.Get_node(n2);

		double nx = -(x2.y - x1.y) / dl;
		double ny = (x2.x - x1.x) / dl;

		// при таком подходе вектор n является внешним по отношению к левой ячейке и внутренним по отношению к правой

		double Fc;

		if (face.is_boundary) {
			int c = max(cr, cl);
			double S = cells[c].Get_S();		// площадь ячейки с

			if (face.zone == 1) {

				//входные параметры
				double uInlet, vInlet, T_Inlet;
				uInlet = p[c].u;
				vInlet = p[c].v;
				T_Inlet = 800.0;

				double un = uInlet * nx + vInlet * ny;

				double h = p[c].Cp * T_Inlet;
				double H = h + 0.5 * (uInlet * uInlet + vInlet * vInlet);

				double E = h / p[c].Gam + 0.5 * (uInlet * uInlet + vInlet * vInlet);

				// Fc - поток через грань
				Fc = p[c].ro * H * un;

				du[c].dU[0] += -Fc * dl / S * dt;
			}

			if (face.zone == 3) {
				double un = p[c].u * nx + p[c].v * ny;
				// Fc - поток через грань
				Fc = p[c].ro * p[c].H * un;
				du[c].dU[0] += Fc * dl / S * dt;
			}

			if (face.zone == 2 || face.zone == 4) {
				du[c].dU[0] += 0;
			}



		}
		else {		// internal face

			// средние значения
			double u_, v_, un_, H_, E_, ro_;

			u_ = 0.5 * (p[cr].u + p[cl].u);
			v_ = 0.5 * (p[cr].v + p[cl].v);
			un_ = u_ * nx + v_ * ny;

			H_ = 0.5 * (p[cr].H + p[cl].H);
			E_ = 0.5 * (p[cr].E + p[cl].E);
			ro_ = 0.5 * (p[cr].ro + p[cl].ro);

			double A = H_ / E_ * un_;
			double Apl, Amn;
			Apl = 0.5 * (A + abs(A));
			Amn = 0.5 * (A - abs(A));

			double UL, UR;
			UL = p[cl].U[0];
			UR = p[cr].U[0];

						// Fc - поток через грань
			Fc = Apl * UL + Amn * UR;



			double Sr = cells[cr].Get_S();		// площадь правой ячейки
			double Sl = cells[cl].Get_S();		// площадь левой ячейки

			du[cr].dU[0] += +Fc * dl / Sr * dt;
			du[cl].dU[0] += -Fc * dl / Sl * dt;

			//cout << i << ", " << cr << ", " << cl << endl;

		}

	}
}

void Yw(Mesh mesh, Cell* (&cells), int nCells)
{
	int nFaces = mesh.Get_nFaces();

	for (int k = 0; k < nCells; k++) {
		double z1 = 1.e10;
		Pnt E = cells[k].Get_c();

		for (int i = 0; i < nFaces; i++) {

			int nz = mesh.Get_face(i).zone;
			int grantype = mesh.Get_zone(nz).grantype;

			if (grantype == 1) {
				int n1 = mesh.Get_face(i).nodes[0];
				int n2 = mesh.Get_face(i).nodes[1];
				Pnt A = mesh.Get_node(n1);
				Pnt B = mesh.Get_node(n2);

				double z2 = Dist(A, B, E);

				if (z1 > z2) z1 = z2;

			}

		}
		cells[k].Yw = z1;

	}

	// создаем поток для записи
	int Nx = mesh.Get_Nx();
	int Ny = mesh.Get_Ny();

	string f = "Yw.plt";
	ofstream record2(f, ios::out);
	if (record2) {




		record2 << "VARIABLES = \"X\", \"Y\", \"Yw,m\"" << endl;

		record2 << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;
		for (int i = 0; i < nCells; i++) {
			record2 << cells[i].Get_c().x << " " << cells[i].Get_c().y << " " << cells[i].Yw << endl;
		}

	}
	record2.close();

}

double Dist(Pnt A, Pnt B, Pnt E)
{
	//Расстояние от точки E до грани AB

	double AB[2], BE[2], AE[2];

	double AB_BE, AB_AE, x, y;
	double x1, y1, x2, y2, mod;

	double dist;

	//vector AB
	AB[0] = B.x - A.x;
	AB[1] = B.y - A.y;

	//vector BE
	BE[0] = E.x - B.x;
	BE[1] = E.y - B.y;

	//vector AE
	AE[0] = E.x - A.x;
	AE[1] = E.y - A.y;

	//Calculating the dot product
	AB_BE = (AB[0] * BE[0] + AB[1] * BE[1]);
	AB_AE = (AB[0] * AE[0] + AB[1] * AE[1]);

	//Case 1
	if (AB_BE > 0) {
		// Finding the magnitude
		y = E.y - B.y;
		x = E.x - B.x;
		dist = sqrt(x * x + y * y);
	}
	//Case 2
	else if (AB_AE < 0) {
		y = E.y - A.y;
		x = E.x - A.x;
		dist = sqrt(x * x + y * y);
	}
	//Case 3
	else {

		//Finding the perpendicular distance
		x1 = AB[0];
		y1 = AB[1];
		x2 = AE[0];
		y2 = AE[1];
		mod = sqrt(x1 * x1 + y1 * y1);
		dist = abs(x1 * y2 - y1 * x2) / mod;

	}

	return dist;
}

void SetGran(Mesh& mesh)
{

	// NS 
	mesh.zones[0].grantype = 3;   // Left gran = Symmetry
	mesh.zones[1].grantype = 1;   // Bottom gran = Wall
	mesh.zones[2].grantype = 4;   // Right gran = Supersonic Outlet
	mesh.zones[3].grantype = 2;   // Top gran = Supersonic Inlet

	/*  1. Wall. 2.Supersonic Inlet. 3.Symmetry
		!!4. Supersonic Outlet.  5.Free boundary
		!6. Subsonic Inlet.     7. Subsonic Outlet */


	int nZones = mesh.Get_nZones();
	for (int i = 0; i < nZones; i++) {
		int tp = mesh.Get_zone(i).grantype;
		mesh.zones[i].bnd = new Boundary[1];
		if (tp == 1) {
			//mesh.zones[i].wall = new Wall[1];
			mesh.zones[i].bnd[0].rules = new int[2];
			mesh.zones[i].bnd[0].vals = new double[1];
		}
		if (tp == 2) {
			mesh.zones[i].bnd[0].rules = new int[1];
			mesh.zones[i].bnd[0].vals = new double[5];

		}

	}
	mesh.zones[1].bnd[0].rules[0] = 1;	// No slip
	mesh.zones[1].bnd[0].rules[1] = 2;	// Tw
	mesh.zones[1].bnd[0].vals[0] = 0.;

	mesh.zones[3].bnd[0].rules[0] = 1;	// Заданы: L, U, T, p, angle
	mesh.zones[3].bnd[0].vals[0] = 0.08; // [ m ] - L
	mesh.zones[3].bnd[0].vals[1] = 2000; // [ m/s ] - U
	mesh.zones[3].bnd[0].vals[2] = 275.1; // [ K ] - T
	mesh.zones[3].bnd[0].vals[3] = 1931.; // [ Pa ] - p
	mesh.zones[3].bnd[0].vals[4] = 0.; // [ m ] - angle


}

void Gradients(Cell* cells, Mesh mesh, Gradient* (&gr), parameters* p, int Nm)
{
	int nCells = mesh.Get_nCells();

	// Значение вектора Vc в центре ячейки
	double* Vc = new double[Nm];
	// Значение вектора Vk-Vc у соседа k :-)
	double* dVk = new double[Nm];

	// временный вектор градиентов для всех компонентов вектора V
	Vector* gr_ = new Vector[Nm];

	for (int m = 0; m < Nm; m++) {
		gr_[m].cx = new double[2];
	}

	for (int i = 0; i < nCells; i++) {

		int NB = cells[i].Get_nFaces();  // число окружающих граней

		for (int m = 0; m < Nm; m++) {
			gr_[m].cx[0] = 0.;
			gr_[m].cx[1] = 0.;
		}

		for (int m = 0; m < Nm; m++)
			Vc[m] = p[i].V[m];


		for (int k = 0; k < NB; k++) {   // проход по всем граням

			// тип грани
			int fType = cells[i].fType[k]; // типы граней:  =0 -> внутренняя грань	(nFaces)
											//   =1 -> граничная грань

			int nb = cells[i].cells[k];		// номер соседней ячейки   ????или грани
			int nf = cells[i].faces[k];		// номер соседней ячейки   ????или грани
				// fType = 0	-> номер ячейки
				// fType = 1	-> номер грани

			//cout << "dVk[0]= " << dVk[0] << endl;

			if (fType == 0) {  // сосед = ячейка
				for (int m = 0; m < Nm; m++) {
					dVk[m] = p[nb].V[m] - Vc[m];

				}
				//cout << "2. dVk[0]= " << dVk[0] << endl;
			}
			else {	// это граница  
				// номер зоны
				int z = mesh.Get_face(nf).zone;
				// boundary type:
				int btype = mesh.Get_zone(z).grantype;

				//cout << "3. dVk[0]= " << dVk[0] << " z= " << z << endl;

				if (btype == 1) { // WALL
					int vel = mesh.Get_zone(z).bnd[0].rules[0];		// 0 - Free slip;  1 - No slip
					int temp = mesh.Get_zone(z).bnd[0].rules[1];		// 1 - Tw is set; 2 - qw is set
					double value = mesh.Get_zone(z).bnd[0].vals[0];	// значение Tw или qw
					if (vel == 0) {
						dVk[1] = 0.;
						dVk[2] = 0.;
					}
					if (vel == 1) {
						dVk[1] = 0. - Vc[1];
						dVk[2] = 0. - Vc[2];
					}

					if (temp == 1) {
						double Tw = value;
						double hw = p[i].Cp * Tw;

						//cout << "4. dVk[0]= " << dVk[0] << endl;

						dVk[3] = hw - Vc[3];

						double ro = p[i].p * p[i].Gm / (8314.41 * Tw);	// ro
						dVk[0] = ro - Vc[0];
					}
					if (temp == 2) {
						dVk[3] = 0.;
						dVk[0] = 0.;
					}

				}

				if (btype == 2) { //Supersonic Inlet
					double u = mesh.Get_zone(z).bnd[0].vals[1];  // U
					double T = mesh.Get_zone(z).bnd[0].vals[2];  // T
					double P = mesh.Get_zone(z).bnd[0].vals[3];  // p

					double ro = P * p[i].Gm / (8314.41 * T);	// ro
					double h = p[i].Cp * T;

					dVk[0] = ro - Vc[0];
					dVk[1] = u - Vc[1];
					dVk[2] = 0. - Vc[2];
					dVk[3] = h - Vc[3];

				}

				if (btype == 3) { //Symmetry

					dVk[0] = 0;
					dVk[1] = 0;
					dVk[2] = 0. - Vc[2];
					dVk[3] = 0;

				}

				if (btype == 4) { //Supersonic Outlet

					dVk[0] = 0.;
					dVk[1] = 0.;
					dVk[2] = 0.;
					dVk[3] = 0.;

				}


			}

			for (int m = 0; m < Nm; m++) {
				gr_[m].cx[0] += cells[i].wk[k] * dVk[m] * cells[i].ck[k].cx[0];
				gr_[m].cx[1] += cells[i].wk[k] * dVk[m] * cells[i].ck[k].cx[1];
			}

		}

		for (int m = 0; m < Nm; m++) {
			gr[i].g[m].cx[0] = gr_[m].cx[0];
			gr[i].g[m].cx[1] = gr_[m].cx[1];
		}
	}

}

void Velocity(parameters* p, Cell* cells, int Nx, int Ny, int nCells)
{
	// создаем поток для записи
	string f = "V.plt";
	ofstream record(f, ios::out);
	if (record) {
		record << "VARIABLES = \"X\", \"Y\", \"u,m_s\", \"v,m_s\" " << endl;

		record << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;
		for (int i = 0; i < nCells; i++) {
			record << cells[i].Get_c().x << " " << cells[i].Get_c().y
				<< " " << p[i].u << " " << p[i].v << endl;
		}

	}
	record.close();
}

void Matrix_Diag(double A[4][4], double L[4], double B[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			B[i][j] = A[i][j] * L[j];
		}
	}
}

void Matrix_Matrix(double A[4][4], double B[4][4], double C[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			C[i][j] = 0;
			for (int k = 0; k < 4; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void Matrix_Vector(double A[4][4], double B[4], double C[4])
{
	for (int i = 0; i < 4; i++) {
		C[i] = 0;
		for (int j = 0; j < 4; j++) {
			C[i] += A[i][j] * B[j];
		}
	}

}

void PrintMatr(double A[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << setw(15) << A[i][j];
		}
		cout << endl;
	}

}

void ConvectNS()
{
	double A[4][4];

	//S_Matr(A);

	double u, v, ro, p, h, nx, ny;
	int iMod;

	u = 3000.;
	v = 1000.;

	p = 1.e5;
	double R = 8314.41 / 28.97;  // J/(kg K)
	double T = 500.;
	ro = p / (R * T);
	double Cp = 1.4 / 0.4 * R;
	h = Cp * T;
	nx = sqrt(0.5);
	ny = sqrt(0.5);

	iMod = 3;

	S_Matr(A, u, v, ro, p, h, nx, ny, iMod);

	PrintMatr(A);

}

void S_Matr(double A[4][4])
{
	A[1][1] = 10;



}

void S_Matr(double A[4][4], double u, double v, double ro, double p, double h, double nx, double ny, int iMod)
{

	double gam = h / (h - p / ro);
	double a = sqrt(gam * p / ro);

	//дополнительные коэффициенты
	double beta = gam - 1.;
	double alfa = 0.5 * (u * u + v * v);

	double Ht = h + alfa;
	double Et = Ht - p / ro;

	//касательный вектор
	double ly = nx;
	double lx = -ny;

	// скорости в системе координат n и l
	double U_ = u * nx + v * ny;
	double V_ = u * lx + v * ly;

	double S[4][4], S_[4][4];

	//Матрица S
	S[0][0] = a * a - alfa * beta;
	S[0][1] = beta * u;
	S[0][2] = beta * v;
	S[0][3] = -beta;

	S[1][0] = -V_;
	S[1][1] = lx;
	S[1][2] = ly;
	S[1][3] = 0.;

	S[2][0] = alfa * beta - U_ * a;
	S[2][1] = a * nx - beta * u;
	S[2][2] = a * ny - beta * v;
	S[2][3] = beta;

	S[3][0] = alfa * beta + U_ * a;
	S[3][1] = -a * nx - beta * u;
	S[3][2] = -a * ny - beta * v;
	S[3][3] = beta;

	//Матрица S ^ -1
	double a2 = a * a;
	S_[0][0] = 1. / a2;
	S_[0][1] = 0.;
	S_[0][2] = 1. / (2. * a2);
	S_[0][3] = 1. / (2. * a2);

	S_[1][0] = u / a2;
	S_[1][1] = lx;
	S_[1][2] = (u + a * nx) / (2. * a2);
	S_[1][3] = (u - a * nx) / (2. * a2);

	S_[2][0] = v / a2;
	S_[2][1] = ly;
	S_[2][2] = (v + a * ny) / (2. * a2);
	S_[2][3] = (v - a * ny) / (2. * a2);

	S_[3][0] = alfa / a2;
	S_[3][1] = V_;
	S_[3][2] = (Ht + a * U_) / (2. * a2);
	S_[3][3] = (Ht - a * U_) / (2. * a2);


	double L[4];

	L[0] = U_;
	L[1] = U_;
	L[2] = U_ + a;
	L[3] = U_ - a;

	double z_static = 0.5;
	double eps = z_static * (abs(U_) + a);

	if (iMod == 1) {		// сама матрица A

	}

	if (iMod == 2) {		// матрица A+
		for (int i = 0; i < 4; i++) {
			//L[i] = 0.5 * ( L[i] + sqrt( L[i] * L[i] ) );

			L[i] = 0.5 * (L[i] + sqrt(L[i] * L[i] + eps * eps));
		}
	}

	if (iMod == 3) {		// матрица A-
		for (int i = 0; i < 4; i++) {
			L[i] = 0.5 * (L[i] - sqrt(L[i] * L[i] + eps * eps));
		}
	}

	double Tmp[4][4];

	Matrix_Diag(S_, L, Tmp);

	Matrix_Matrix(Tmp, S, A);

}

void ConvectNS(parameters* p, changes* (&du), Mesh mesh, Cell* cells, double dt, int Nm)
{
	int nFaces = mesh.Get_nFaces();

	double* Fc = new double[Nm];
	double* UL = new double[Nm];
	double* UR = new double[Nm];
	double* ULR = new double[Nm];

	double A[4][4];

	double* AUl = new double[Nm];
	double* AUr = new double[Nm];

	for (int i = 0; i < nFaces; i++) {

		Face face = mesh.Get_face(i);

		int cr = face.cr;
		int cl = face.cl;

		// nodes:
		int n1 = face.nodes[0];
		int n2 = face.nodes[1];

		double dl = face.length; // длина грани

		Pnt x1 = mesh.Get_node(n1);
		Pnt x2 = mesh.Get_node(n2);

		double nx = -(x2.y - x1.y) / dl;
		double ny = (x2.x - x1.x) / dl;

		// при таком подходе вектор n является внешним по отношению к левой ячейке и внутренним по отношению к правой


		if (face.is_boundary) {
			int c = max(cr, cl);
			double S = cells[c].Get_S();		// площадь ячейки внутренней ячейки с

			int z = face.zone;

			int grantype = mesh.Get_zone(z).grantype;

			if (grantype == 1 || grantype == 3) {
				Fc[0] = 0;
				Fc[1] = p[c].p * nx;
				Fc[2] = p[c].p * ny;
				Fc[3] = 0;
			}


			if (grantype == 2) {	// Supersonic Inlet

				// скорость на входе
				double VV = mesh.Get_zone(z).bnd[0].vals[1];
				// угол потока
				double angle = mesh.Get_zone(z).bnd[0].vals[4];
				// компоненты скорости
				double u_ = VV * cos(angle);
				double v_ = VV * sin(angle);

				double T_ = mesh.Get_zone(z).bnd[0].vals[2];
				double p_ = mesh.Get_zone(z).bnd[0].vals[3];

				double Gm = 28.97;
				double Gam = 1.4;

				// расчет
				double R = 8314.41 / Gm;  // J/(kg K)
				double ro_ = p_ / (R * T_);

				//cout << "u_= " << u_ << "ro_= " << ro_ << endl;
				//exit(27);

				double h_ = p_ * Gam / ((Gam - 1.) * ro_);

				double H_ = h_ + 0.5 * (u_ * u_ + v_ * v_);

				double un = u_ * nx + v_ * ny;

				Fc[0] = ro_ * un;
				Fc[1] = ro_ * un * u_ + p_ * nx;
				Fc[2] = ro_ * un * v_ + p_ * ny;
				Fc[3] = ro_ * un * H_;

			}

			if (grantype == 4) {	// Supersonic Outlet
				double un = p[c].u * nx + p[c].v * ny;
				double ro_ = p[c].ro;

				Fc[0] = ro_ * un;
				Fc[1] = ro_ * un * p[c].u + p[c].p * nx;
				Fc[2] = ro_ * un * p[c].v + p[c].p * ny;
				Fc[3] = ro_ * un * p[c].H;

			}


			if (cr >= 0) {
				double Sr = cells[cr].Get_S();		// площадь правой ячейки
				for (int m = 0; m < Nm; m++) {
					du[cr].dU[m] += +Fc[m] * dl / Sr * dt;
				}
			}

			if (cl >= 0) {
				double Sl = cells[cl].Get_S();		// площадь левой ячейки
				for (int m = 0; m < Nm; m++) {
					du[cl].dU[m] += -Fc[m] * dl / Sl * dt;
				}
			}

		}
		else {		// internal face

			for (int m = 0; m < Nm; m++) {
				UL[m] = p[cl].U1[m];
				UR[m] = p[cr].U1[m];
				ULR[m] = 0.5 * (UL[m] + UR[m]);
			}

			// средние значения
			double u_, v_, h_, E_, ro_, p_;

			ro_ = ULR[0];
			u_ = ULR[1] / ro_;
			v_ = ULR[2] / ro_;
			E_ = ULR[3] / ro_;

			double q2 = 0.5 * (u_ * u_ + v_ * v_);

			double e_ = E_ - q2;

			double gam_ = 0.5 * (p[cr].Gam + p[cl].Gam);

			h_ = e_ * gam_;

			p_ = ro_ * (gam_ - 1) * e_;

			// Матрица A+
			int iMod = 2;

			S_Matr(A, u_, v_, ro_, p_, h_, nx, ny, iMod);

			// Умножаем A+ * UL
			Matrix_Vector(A, UL, AUl);

			// Матрица A-
			iMod = 3;

			S_Matr(A, u_, v_, ro_, p_, h_, nx, ny, iMod);

			// Умножаем A- * UR
			Matrix_Vector(A, UR, AUr);

			// Вектор невязких потоков
			//Fc[0] = 1;
			for (int m = 0; m < Nm; m++) {
				Fc[m] = AUl[m] + AUr[m];
			}

			double Sr = cells[cr].Get_S();		// площадь правой ячейки
			double Sl = cells[cl].Get_S();		// площадь левой ячейки

			// Приращения осн.вектора за счет невязких потоков
			for (int m = 0; m < Nm; m++) {
				du[cr].dU[m] += +Fc[m] * dl / Sr * dt;
				du[cl].dU[m] += -Fc[m] * dl / Sl * dt;
			}


			//cout << i << ", " << cr << ", " << cl << endl;

		}


	}

}

void Pressure(parameters* p, Cell* cells, int Nx, int Ny, int nCells)
{
	// создаем поток для записи
	string f = "p.plt";
	ofstream record(f, ios::out);
	if (record) {
		record << "VARIABLES = \"X\", \"Y\", \"p,Pa\" " << endl;

		record << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;
		for (int i = 0; i < nCells; i++) {
			record << cells[i].Get_c().x << " " << cells[i].Get_c().y
				<< " " << p[i].p << endl;
		}

	}
	record.close();
}

void Mach(parameters* p, Cell* cells, int Nx, int Ny, int nCells)
{
	// создаем поток для записи
	string f = "Mach.plt";
	ofstream record(f, ios::out);
	if (record) {
		record << "VARIABLES = \"X\", \"Y\", \"Mach\" " << endl;

		record << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;
		for (int i = 0; i < nCells; i++) {

			double vv = sqrt(p[i].u * p[i].u + p[i].v * p[i].v);
			double a = sqrt(p[i].Gam * p[i].p / p[i].ro);

			record << cells[i].Get_c().x << " " << cells[i].Get_c().y
				<< " " << vv / a << endl;
		}

	}
	record.close();
}
