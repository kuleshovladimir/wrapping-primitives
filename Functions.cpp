#include "Functions.h"

void Init(parameters* (&p), int nCells, int Nm) {

    double T0, Cp, la, P0, Gm, Gam, R, ro, mu;

    double U0;

    // входные данные 
    T0 = 293.0; //K
    la = 2.5658e-2;  // W/(m K)
    mu = 17.863e-6;

    U0 = 0.1;

    P0 = 101325.0;  //Pa
    Gam = 1.4;
    Gm = 28.97;

    //расчет 
    R = 8314.41 / Gm; //J/(kg K)
    ro = P0 / (R * T0);
    Cp = Gam / (Gam - 1.) * R;

    for (int i = 0; i < nCells; i++) {
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

        p[i].U = new double[Nm];
        p[i].U1 = new double[Nm];

        p[i].V = new double[Nm];

        p[i].U1[0] = p[i].ro * p[i].E;

        p[i].V[0] = p[i].h;
    }

}

void Viscous(parameters* p, changes* du, Mesh mesh, Cell* cells, double dt)
{
	int nFaces = mesh.Get_nFaces();

	for (int i = 0; i < nFaces; i++) {
		Face face = mesh.Get_face(i);

		int cr = face.cr;
		int cl = face.cl;

		if (face.is_boundary) {
            int c = max(cr, cl);

            Pnt xc = cells[c].Get_MassC();

            Pnt xf = face.f_center;

            double dx = xc.x - xf.x;
            double dy = xc.y - xf.y;
            double dl = sqrt(dx * dx + dy * dy);

            double length = face.length;        //длина грани
            double S = cells[c].Get_S();        //площадь ячейки с

            int z = face.zone;

            int btype = mesh.Get_zone(z).granType;

            double Tw;
            if (face.zone == 2) Tw = 500;

            if (face.zone == 4) Tw = 200;

            if (face.zone == 2 || face.zone == 4) {
                double hw = p[c].Cp * Tw;

                //dh/dn
                double dh_dn = (p[c].h - hw) / dl;

                //mu/Pr
                double mu_Pr = p[c].mu / p[c].Pr;

                //Fv - поток через грань
                double Fv = mu_Pr * dh_dn;

                du[c].dU[0] += -Fv * length / S * dt;
            }
            if (face.zone == 1) {
                double Fv = 100;
                du[c].dU[0] += Fv * length / S * dt;
            }
            
            if (face.zone == 3) {
                double Fv = -10;
                du[c].dU[0] += Fv * length / S * dt;
            }


		}
		else {
			//координаты правой и левой ячеек
			Pnt xr = cells[cr].Get_MassC();
			Pnt xl = cells[cl].Get_MassC();

			//расстояние между ячейками
			double dx = (xr.x - xl.y);
			double dy = (xr.y - xl.y);
			double dl = sqrt(dx * dx + dy * dy);

			//dh/dn
			double dh_dn = (p[cr].h - p[cl].h) / dl;

			//среднее mu/Pr
			double mu_Pr = 0.5 * (p[cr].mu / p[cr].Pr + p[cl].mu / p[cl].Pr);

			//Fv - поток через грань
			double Fv = mu_Pr * dh_dn;

			double lenght = face.length;	//длина грани
			double Sr = cells[cr].Get_S();	//площади 
			double Sl = cells[cl].Get_S();

			du[cr].dU[0] += -Fv * lenght / Sr * dt;
			du[cl].dU[0] += Fv * lenght / Sl * dt;

		}
	}
}

void GetParams(parameters* (&p), int nCells, int Nm)
{
    for (int i = 0; i < nCells; i++) {
        p[i].E = p[i].U1[0] / p[i].ro;

        double q = 0.5 * (p[i].u * p[i].u + p[i].v * p[i].v);
        double e = p[i].E - q;

        p[i].h = e * p[i].Gam;
        p[i].T = p[i].h / p[i].Cp;
        p[i].H = p[i].h + q;

        double R = 8314.41 / p[i].Gm;
        p[i].p = p[i].ro * R * p[i].T;

        p[i].V[0] = p[i].h;

    }
}

//***********
void Convect(parameters* (&p), changes* (&du), Mesh mesh, Cell* cells, int It, double dt)
{
    int nFaces = mesh.Get_nFaces();

    for (int i = 0; i < nFaces; i++) {
        Face face = mesh.Get_face(i);

        int cr = face.cr;
        int cl = face.cl;

        //nodes
        int n1 = face.nodes[0];
        int n2 = face.nodes[1];

        double dl = face.length;    //длина грани

        Pnt x1 = mesh.Get_node(n1);
        Pnt x2 = mesh.Get_node(n2);

        double nx = (x1.y - x2.y) / dl;
        double ny = (x2.x - x1.x) / dl;

        double Fc;

        if (i == -10) {
            double un = p[cl].u * nx + p[cl].v * ny;
        }

        if (face.is_boundary) {
            int c = max(cr, cl);
            double S = cells[c].Get_S();


            if (face.zone == 1) {
                double uInlet, vInlet, T_Inlet;
                uInlet = p[c].u;
                vInlet = p[c].v;
                T_Inlet = 800;

                double un = uInlet * nx + vInlet * ny;

                double h = p[c].Cp * T_Inlet;
                double H = h + 0.5 * (uInlet * uInlet + vInlet * vInlet);

                double E = h / p[c].Gam + 0.5 * (uInlet * uInlet + vInlet * vInlet);

                Fc = p[c].ro * H * un;

                du[c].dU[0] += -Fc * dl / S * dt;

            }

            if (face.zone == 3) {
                double un = p[c].u * nx + p[c].v * ny;

                Fc = p[c].ro * p[c].H * un;
                du[c].dU[0] += Fc * dl / S * dt;
            }

            if ((face.zone == 2) || (face.zone == 4)) {
                du[c].dU[0] += 0;
            }
        }
        else {
            //средние значения
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

            Fc = Apl * UL + Amn * UR;

            double Sr = cells[cr].Get_S();
            double Sl = cells[cl].Get_S();

            du[cr].dU[0] += Fc * dl / Sr * dt;
            du[cl].dU[0] -= Fc * dl / Sl * dt;


        }
    }
}

void Tecplot(parameters* p, Cell* cells, int Nx, int Ny, int nCells)
{
    // создаем поток для записи
    string f = "T.plt";
    ofstream record(f, ios::out);
    if (record) {
        record << "VARIABLES = \"X\", \"Y\", \"T,K\"" << endl;

        record << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;
        for (int i = 0; i < nCells; i++) {
            record << cells[i].Get_MassC().x << " " << cells[i].Get_MassC().y << " " << p[i].T << endl;
        }

    }
    record.close();
}



double Dist(Pnt A, Pnt B, Pnt E)
{
    double AB[2], BE[2], AE[2];

    double AB_BE, AB_AE;

    double dist;

    //векторы
    AB[0] = B.x - A.x;
    AB[1] = B.y - A.y;

    BE[0] = E.x - B.x;
    BE[1] = E.y - B.y;

    AE[0] = E.x - A.x;
    AE[1] = E.y - A.y;

    //сколярное произведение векторов
    AB_BE = (AB[0] * BE[0] + AB[1] * BE[1]);
    AB_AE = (AB[0] * AE[0] + AB[1] * AE[1]);

    if (AB_BE > 0) {
        double y = E.y - B.y;
        double x = E.x - B.x;
        dist = sqrt(x * x + y * y);
    }
    else if (AB_BE < 0) {
        double y = E.y - A.y;
        double x = E.x - A.x;
        dist = sqrt(x * x + y * y);
    }
    else {
        dist = abs(AB[0] * AE[1] - AB[1] * AE[0]) / sqrt(AB[0] * AB[0] + AB[1] * AB[1]);
    }

    return dist;
}

void SetGran(Mesh& mesh)
{
    mesh.zones[0].granType = 1;
    mesh.zones[1].granType = 1;
    mesh.zones[2].granType = 1;
    mesh.zones[3].granType = 1;

    int nZones = mesh.nZones;

    //в будущем - ввод из файла
    //тестовые значения

    for (int i = 0; i < nZones; i++) {
        int tp = mesh.zones[i].granType;
        if (tp == 1) {
            mesh.zones[i].wall = new Wall;

        }
        mesh.zones[0].wall->vel = 1;
        mesh.zones[0].wall->temp = 2;
        mesh.zones[0].wall->value = 0;

        mesh.zones[1].wall->vel = 1;
        mesh.zones[1].wall->temp = 1;
        mesh.zones[1].wall->value = 500;

        mesh.zones[2].wall->vel = 1;
        mesh.zones[2].wall->temp = 2;
        mesh.zones[2].wall->value = 0;

        mesh.zones[3].wall->vel = 1;
        mesh.zones[3].wall->temp = 1;
        mesh.zones[3].wall->value = 200;
    }
}

void Gradients(Cell* cells, Mesh mesh, Gradient* (&gr), parameters* p, int Nm)
{
    int nCells = mesh.Get_nCells();

    //значение вектора Vc в центре ячейки
    double* Vc = new double[Nm];
    //изменение вектора Vc по отношению к соседней клетке
    double* dVk = new double[Nm];

    //временный вектор градиентов для всех компнентов вектора V
    Vector* gr_ = new Vector[Nm];

    for (int i = 0; i < Nm; i++) {
        gr_[i].cx = new double[2];
    }

    for (int i = 0; i < nCells; i++) {
        int NB = cells[i].Get_NFaces();     //число окружающих граней

        for (int j = 0; j < Nm; j++) {
            gr_[j].cx[0] = 0;
            gr_[j].cx[1] = 0;
        }

        for (int j = 0; j < Nm; j++) {
            Vc[j] = p[i].V[j];
        }

        for (int j = 0; j < NB; j++) {      //проход по граням
            bool fType = cells[i].Get_fType(j);

            int nb = cells[i].Get_Cell(j);
            int nf = cells[i].Get_Face(j);

            if (!fType) {
                for (int m = 0; m < Nm; m++) {
                    dVk[m] = p[nb].V[m] - Vc[m];
                }
            }
            else {
                int z = mesh.Get_face(nf).zone;
                int btype = mesh.Get_zone(z).granType;

                if (btype == 1) {
                    int vel = mesh.Get_zone(z).wall->vel;
                    int temp = mesh.Get_zone(z).wall->temp;
                    double value = mesh.Get_zone(z).wall->value;

                    if (temp == 1) {
                        double Tw = value;
                        double hw = p[i].Cp * Tw;

                        dVk[0] = hw - Vc[0];
                    }

                    if (temp == 2) {
                        dVk[0] = 0;
                    }
                }
            }

            gr_[0].cx[0] += cells[i].Get_wk(j) * dVk[0] * cells[i].Get_ck(j).cx[0];
            gr_[0].cx[1] += cells[i].Get_wk(j) * dVk[0] * cells[i].Get_ck(j).cx[1];

        }

        for (int j = 0; j < Nm; j++) {
            gr[i].g[j].cx[0] = gr_[j].cx[0];
            gr[i].g[j].cx[1] = gr_[j].cx[1];
        }
    }
}
