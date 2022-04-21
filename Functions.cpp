#include "Functions.h"

void Init(parameters* (&p), int nCells, int Nm) {

    double T0, Cp, la, P0, Gm, Gam, R, ro, mu;

    double U0;

    // входные данные 
    T0 = 273; //температура вход€щего газа
    la = 2.5658e-2;  // коэффициент теплопроводности 
    mu = 17.863e-6;  // в€зкость

    U0 = 867.9; //скорость

    P0 = 1931;  //давление
    Gam = 1.4;  //константа адиабаты
    Gm = 28.97; //мол€рна€ масса

    //расчет 
    R = 8314.41 / Gm; 
    ro = P0 / (R * T0); //плотность
    Cp = Gam / (Gam - 1.) * R; //теплоемкость

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

        // NS: 

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
            double S = cells[c].Get_S();        //площадь €чейки с

            int z = face.zone;

            int btype = mesh.Get_zone(z).granType;

            double Tw;
            if (face.zone == 1) Tw = 500;

            if (face.zone == 3) Tw = 200;

            if (face.zone == 1 || face.zone == 3) {
                double hw = p[c].Cp * Tw;

                //dh/dn
                double dh_dn = (p[c].h - hw) / dl;

                //mu/Pr
                double mu_Pr = p[c].mu / p[c].Pr;

                //Fv - поток через грань
                double Fv = mu_Pr * dh_dn;

                du[c].dU[0] += -Fv * length / S * dt;
            }
            if (face.zone == 0) {
                double Fv = 100;
                du[c].dU[0] += Fv * length / S * dt;
            }
            
            if (face.zone == 2) {
                double Fv = -10;
                du[c].dU[0] += Fv * length / S * dt;
            }


		}
		else {
			//координаты правой и левой €чеек
			Pnt xr = cells[cr].Get_MassC();
			Pnt xl = cells[cl].Get_MassC();

			//рассто€ние между €чейками
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

void Viscous(parameters* p, changes* du, Mesh mesh, Cell* cells, double dt, Gradient* gr, int Nm)
{
    int nFaces = mesh.Get_nFaces();

    double* Fv = new double[Nm];

    for (int i = 0; i < nFaces; i++) {
        Face face = mesh.Get_face(i);

        int cr = face.cr;
        int cl = face.cl;

        double length = face.length;

        int n1 = face.nodes[0];
        int n2 = face.nodes[1];

        Pnt x1 = mesh.Get_node(n1);
        Pnt x2 = mesh.Get_node(n2);

        double nx = (x1.y - x2.y) / length;
        double ny = (x2.x - x1.x) / length;

        if (face.is_boundary) {
            int c = max(cr, cl);

            Pnt xc = cells[c].Get_MassC();

            int z = face.zone;
            int granType = mesh.Get_zone(z).granType;

            if (granType == 1) {
                double dl = cells[i].Get_Yw();

                double du_dx, du_dy, dv_dx, dv_dy, dh_dx, dh_dy;
                double  u_, v_;

                if (mesh.Get_zone(z).bnd[0].rules[0] == 0) {
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

                if (mesh.Get_zone(z).bnd[0].rules[0] == 1) {
                    double du_dn = p[c].u / dl;
                    double dv_dn = p[c].v / dl;
                    double du_dl = 0;
                    double dv_dl = 0;
                    du_dx = du_dn * nx - du_dl * ny;
                    dv_dx = dv_dn * nx - dv_dn * ny;

                    du_dy = du_dn * ny + du_dl * nx;
                    dv_dy = dv_dn * ny + dv_dl * nx;

                    u_ = 0;
                    v_ = 0;
                 }

                if (mesh.Get_zone(z).bnd[0].rules[1] == 1) {
                    double Tw = mesh.Get_zone(z).bnd[0].vals[0];
                    double hw = p[c].Cp * Tw;
                    double dh_dn = (p[c].h - hw) / dl;

                    double dh_dl = 0;

                    dh_dx = dh_dn * nx + dh_dl * ny;
                    dh_dy = dh_dn * ny + dh_dl * nx;

                }

                if (mesh.Get_zone(z).bnd[0].rules[1] == 2) {
                    double tx = gr[c].g[3].cx[0];
                    double ty = gr[c].g[3].cx[1];

                    dh_dx = ny * ny * tx - nx * nx * ty;
                    dh_dy = -nx * ny * tx + nx * nx * ty;

                }

                double mu = p[c].mu;

                double div = du_dx + dv_dy;
                double txx = mu * (2 * du_dx - 2 * div / 3);
                double tyy = mu * (2 * dv_dy - 2 * div / 3);
                double txy = mu * (du_dy + dv_dx);

                double mu_Pr = p[c].mu / p[c].Pr;

                double qx = -mu_Pr * dh_dx;
                double qy = -mu_Pr * dh_dy;

                Fv[0] = 0;
                Fv[1] = txx * nx + txy * ny;
                Fv[2] = txy * nx + tyy * ny;
                Fv[3] = (u_ * txx + v_ * txy - qx) * nx + (u_ * txy + v_ * tyy - qy) * ny;

            } 

            if (granType == 3) {
                double dl = cells[i].Get_Yw();

                double du_dx, du_dy, dv_dx, dv_dy, dh_dx, dh_dy;
                double  u_, v_;

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

                double tx = gr[c].g[3].cx[0];
                double ty = gr[c].g[3].cx[1];

                dh_dx = ny * ny * tx - nx * nx * ty;
                dh_dy = -nx * ny * tx + nx * nx * ty;

                double mu = p[c].mu;

                double div = du_dx + dv_dy;
                double txx = mu * (2 * du_dx - 2 * div / 3);
                double tyy = mu * (2 * dv_dy - 2 * div / 3);
                double txy = mu * (du_dy + dv_dx);

                double mu_Pr = p[c].mu / p[c].Pr;

                double qx = -mu_Pr * dh_dx;
                double qy = -mu_Pr * dh_dy;

                Fv[0] = 0;
                Fv[1] = txx * nx + txy * ny;
                Fv[2] = txy * nx + tyy * ny;
                Fv[3] = (u_ * txx + v_ * txy - qx) * nx + (u_ * txy + v_ * tyy - qy) * ny;

            }

            if (granType == 2 || granType == 4) {
                for (int m = 0; m < 4; m++) {
                    Fv[m] = 0;
                }
            }

            double length = face.length;
            if (cr >= 0) {
                double Sr = cells[cr].Get_S();
                for (int m = 0; m < Nm; m++) {
                    du[cr].dU[m] += -Fv[m] * length / Sr * dt;
                }
            }

            if (cl >= 0) {
                double Sl = cells[cl].Get_S();
                for (int m = 0; m < Nm; m++) {
                    du[cl].dU[m] += Fv[m] * length / Sl * dt;
                }
            }

        }
        else {


            Pnt xr = cells[cr].Get_MassC();
            Pnt xl = cells[cl].Get_MassC();

            double dx = xr.x - xl.x;
            double dy = xr.y - xl.y;
            double dl = sqrt(dx * dx + dy * dy);

            double du_dx = 0.5 * (gr[cr].g[1].cx[0] + gr[cl].g[1].cx[0]);
            double du_dy = 0.5 * (gr[cr].g[1].cx[1] + gr[cl].g[1].cx[1]);

            double dv_dx = 0.5 * (gr[cr].g[2].cx[0] + gr[cl].g[2].cx[0]);
            double dv_dy = 0.5 * (gr[cr].g[2].cx[1] + gr[cl].g[2].cx[1]);

            double dh_dx = 0.5 * (gr[cr].g[3].cx[0] + gr[cl].g[3].cx[0]);
            double dh_dy = 0.5 * (gr[cr].g[3].cx[1] + gr[cl].g[3].cx[1]);

            double mu = 0.5 * (p[cr].mu + p[cl].mu);

            double div = du_dx + dv_dy;
            double txx = mu * (2 * du_dx - 2 * div / 3);
            double tyy = mu * (2 * dv_dy - 2 * div / 3);
            double txy = mu * (du_dy + dv_dx);

            double mu_Pr = 0.5 * (p[cr].mu / p[cr].Pr + p[cl].mu / p[cl].Pr);

            double qx = -mu_Pr * dh_dx;
            double qy = -mu_Pr * dh_dy;
            
            double u_ = 0.5 * (p[cr].u + p[cl].u);
            double v_ = 0.5 * (p[cr].v + p[cl].v);

            Fv[0] = 0;
            Fv[1] = txx * nx + txy * ny;
            Fv[2] = txy * nx + tyy * ny;
            Fv[3] = (u_ * txx + v_ * txy - qx) * nx + (u_ * txy +v_ * tyy - qy) * ny;

            double lenght = face.length;
            double Sr = cells[cr].Get_S();
            double Sl = cells[cl].Get_S();

            for (int m = 0; m < Nm; m++) {
                du[cr].dU[m] += -Fv[m] * length / Sr * dt;
                du[cl].dU[m] += Fv[m] * length / Sl * dt;
            }

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

        double q = 0.5 * (p[i].u * p[i].u + p[i].v * p[i].v);

        double e = p[i].E - q;

        p[i].h = e * p[i].Gam;
        p[i].T = p[i].h / p[i].Cp;
        p[i].H = p[i].h + q;

        double R = 8314.41 / p[i].Gm;
        p[i].p = p[i].ro * R * p[i].T;

        p[i].V[0] = p[i].ro;    
        p[i].V[1] = p[i].u;
        p[i].V[2] = p[i].v;
        p[i].V[3] = p[i].h;

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
                T_Inlet = 273;

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
            //средние значени€
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

void Yw(Mesh mesh, Cell* (&cells), int nCells)
{
    int nFaces = mesh.Get_nFaces();

    for (int k = 0; k < nCells; k++) {
        double z1 = 1.e10;
        Pnt E = cells[k].Get_MassC();

        for (int i = 0; i < nFaces; i++) {
            int nz = mesh.Get_face(i).zone;
            int granType = mesh.Get_zone(nz).granType;

            if (granType == 1) {
                int n1 = mesh.Get_face(i).nodes[0];
                int n2 = mesh.Get_face(i).nodes[1];
                Pnt A = mesh.Get_node(n1);
                Pnt B = mesh.Get_node(n2);


                double z2 = Dist(A, B, E);

                if (z1 > z2) { z1 = z2; };
            }
        }
        cells[k].Set_Yw(z1);
    }

    //поток записи
    int Nx = mesh.Get_Nx();
    int Ny = mesh.Get_Ny();

    string f = "Yw.plt";
    ofstream record2(f, ios::out);
    if (record2) {
        record2 << "VARIABLES= \"X\", \"Y\", \"Yw,m\"" << endl;
        record2 << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;

        for (int i = 0; i < nCells; i++) {
            record2 << cells[i].Get_MassC().x << " " << cells[i].Get_MassC().y << " " << cells[i].Get_Yw() << endl;
        }

    }
    record2.close();
}

void Tecplot(parameters* p, Cell* cells, int Nx, int Ny, int nCells)
{
    // создаем поток дл€ записи
    string f = "T.plt";
    ofstream record(f, ios::out);
    if (record) {
        record << "VARIABLES= \"X\", \"Y\", \"T,K\"" << endl;

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

    //скол€рное произведение векторов
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
    mesh.zones[0].granType = 3;
    mesh.zones[1].granType = 1;
    mesh.zones[2].granType = 4;
    mesh.zones[3].granType = 2;

    int nZones = mesh.nZones;

    //в будущем - ввод из файла
    //тестовые значени€

    for (int i = 0; i < nZones; i++) {
        int tp = mesh.zones[i].granType;
        mesh.zones[i].bnd = new Boudary;
        if (tp == 1) {
            mesh.zones[i].bnd[0].rules = new int[2];
            mesh.zones[i].bnd[0].vals = new double;
        }
        if (tp == 2) {
            mesh.zones[i].bnd[0].rules = new int;
            mesh.zones[i].bnd[0].vals = new double[5];
        }
    }

    mesh.zones[1].bnd[0].rules[0] = 1;
    mesh.zones[1].bnd[0].rules[1] = 1;
    mesh.zones[1].bnd[0].vals[0] = 105;
    
    mesh.zones[3].bnd[0].rules[0] = 1;
    mesh.zones[3].bnd[0].vals[0] = 0.08;
    mesh.zones[3].bnd[0].vals[1] = 867.9;
    mesh.zones[3].bnd[0].vals[2] = 75.1;
    mesh.zones[3].bnd[0].vals[3] = 1931;
    mesh.zones[3].bnd[0].vals[4] = 0;
    
}

void Gradients(Cell* cells, Mesh mesh, Gradient* (&gr), parameters* p, int Nm)
{
    int nCells = mesh.Get_nCells();

    //значение вектора Vc в центре €чейки
    double* Vc = new double[Nm];
    //изменение вектора Vc по отношению к соседней клетке
    double* dVk = new double[Nm];

    //временный вектор градиентов дл€ всех компнентов вектора V
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

        for (int j = 0; j < NB; j++) {      //проход по гран€м
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
                    bool vel = mesh.Get_zone(z).bnd[0].rules[0];
                    int temp = mesh.Get_zone(z).bnd[0].rules[1];
                    double value = mesh.Get_zone(z).bnd[0].vals[0];

                    if (!vel) {
                        dVk[1] = 0;
                        dVk[2] = 0;

                    }
                    else {
                        dVk[1] = -Vc[1];
                        dVk[2] = -Vc[2];
                    }

                    if (temp == 1) {
                        double Tw = value;
                        double hw = p[i].Cp * Tw;

                        dVk[3] = hw - Vc[3];

                        double ro = p[i].p * p[i].Gm / (8314.41 * Tw);
                        dVk[0] = ro - Vc[0];
                    }

                    if (temp == 2) {
                        dVk[3] = 0;
                        dVk[0] = 0;
                    }
                }

                if (btype == 2) {
                    double u = mesh.Get_zone(z).bnd->vals[1];
                    double T = mesh.Get_zone(z).bnd->vals[2];
                    double P = mesh.Get_zone(z).bnd->vals[3];

                    double ro = P * p[i].Gm / (8314.41 * T);
                    double h = p[i].Cp * T;

                    dVk[0] = ro - Vc[0];
                    dVk[1] = u - Vc[1];
                    dVk[2] = -Vc[2];
                    dVk[3] = h - Vc[3];
                }

                if (btype == 3) {
                    dVk[0] = 0;
                    dVk[1] = 0;
                    dVk[2] = -Vc[2];
                    dVk[3] = 0;
                }

                if (btype == 4) {
                    dVk[0] = 0;
                    dVk[1] = 0;
                    dVk[2] = 0;
                    dVk[3] = 0;
                }
            }

            for (int m = 0; m < Nm; m++) {
                gr_[m].cx[0] += cells[i].Get_wk(j) * dVk[m] * cells[i].Get_ck(j).cx[0];
                gr_[m].cx[1] += cells[i].Get_wk(j) * dVk[m] * cells[i].Get_ck(j).cx[1];
            }
        }

        for (int j = 0; j < Nm; j++) {
            gr[i].g[j].cx[0] = gr_[j].cx[0];
            gr[i].g[j].cx[1] = gr_[j].cx[1];
        }
    }
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
                C[i][j] += A[i][k] + B[k][j];
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

void PrintMatrix(double A[4][4])
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << setw(15) << A[i][j];
        }
        cout << endl;
    }
}

/*void ConvectNS()
{
    double A[4][4];
    double u, double v, double ro, double p, double h, double nx, double ny;
    int iMod;
    S_Matr(A, u, v, ro, p, h, nx, ny, iMod);
}*/


void S_Matr(double A[4][4], double u, double v, double ro, double p, double h, double nx, double ny, int iMod)
{
    double gam = h / (h - p / ro);
    double a = sqrt(gam * p / ro);

    double beta = gam - 1;
    double alpha = 0.5 * (u * u + v * v);

    double Ht = h + alpha;
    double Et = Ht - p / ro;

    double ly = nx;
    double lx = -ny;

    double U_ = u * nx + v * ny;
    double V_ = u * lx + v * ly;

    double S[4][4], S_[4][4];

    S[0][0] = a * a - alpha * alpha;
    S[0][1] = beta * u;
    S[0][2] = beta * v;
    S[0][3] = -beta;

    S[1][0] = -V_;
    S[1][1] = lx;
    S[1][2] = ly;
    S[1][3] = 0;

    S[2][0] = alpha * beta - U_ * a;
    S[2][1] = a * nx - beta * u;
    S[2][2] = a * ny - beta * v;
    S[2][3] = beta;

    S[3][0] = alpha * beta + U_ * a;
    S[3][1] = -a * nx - beta * u;
    S[3][2] = -a * ny - beta * v;
    S[3][3] = beta;

    double a2 = a * a;

    S_[0][0] = 1 / a2;
    S_[0][1] = 0;
    S_[0][2] = 1 / (2 * a2);
    S_[0][3] = 1 / (2 * a2);

    S_[1][0] = u / a2;
    S_[1][1] = lx;
    S_[1][2] = (u + a * nx) / (2 * a2);
    S_[1][3] = (u - a * nx) / (2 * a2);

    S_[2][0] = v / a2;
    S_[2][1] = ly;
    S_[2][2] = (v + a * ny) / (2 * a2);
    S_[2][3] = (v - a * ny) / (2 * a2);

    S_[3][0] = alpha / a2;
    S_[3][1] = V_;
    S_[3][2] = (Ht + a * U_) / (2 * a2);
    S_[3][3] = (Ht - a * U_) / (2 * a2);

    double L[4];

    L[0] = U_;
    L[1] = U_;
    L[2] = U_ + a;
    L[3] = U_ - a;

    if (iMod == 1) { //матрица ј
        
    }

    if (iMod == 2) { //матрица B
        for (int i = 0; i < 4; i++) {
            L[0] = 0.5 * (L[i] + abs(L[i]));
        }
    }

    if (iMod == 3) { //матрица —
        for (int i = 0; i < 4; i++) {
            L[0] = 0.5 * (L[i] - abs(L[i]));
        }
    }

    double Tmp[4][4];

    Matrix_Diag(S_, L, Tmp);

    Matrix_Matrix(Tmp, S, A);
}
