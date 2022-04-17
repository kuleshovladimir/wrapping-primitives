#include "Region.h"

Region::Region() {
	gas.Lchar = 1.;
	gas.U[0] = 1.;
	gas.U[1] = 0.;
	gas.U[2] = 0.;
	gas.p = 1.e5;
	gas.T = 293.;

	iMol = 0;
	Nc = 0;
	double* Yc = new double[Nc];

	turb.turb_intens = 0.01;
	turb.turb_scale = 0.01;
	turb.MuT_Mu = 500.;

}

Region::Region(Gas gas_, int iMol_, double* Yc_, Turb turb_, int Nc_) {

	gas.Lchar = gas_.Lchar;
	gas.U[0] = gas_.U[0];
	gas.U[1] = gas_.U[1];
	gas.U[2] = gas_.U[2];
	gas.p = gas_.p;
	gas.T = gas_.T;

	iMol = iMol_;
	Nc = Nc_;
	double* Yc = new double[Nc];
	for (int i = 0; i < Nc_; i++)
		Yc[i] = Yc_[i];

	turb.turb_intens = turb_.turb_intens;
	turb.turb_scale = turb_.turb_scale;
	turb.MuT_Mu = turb_.MuT_Mu;

}

Region::~Region()
{

}

void Region::Print(int* iComps, string* comps)
{
	cout << "Хар.размер области: " << gas.Lchar << endl;
	cout << "Компоненты скорости: " << gas.U[0] << ", " << gas.U[1] << ", " << gas.U[2] << endl;
	cout << " Температура = " << gas.T << " [K]";
	cout << ", Давление = " << gas.p << " [Па]" << endl;

	if (iMol == 0)
		cout << " Массовые доли компонентов: " << endl;
	else
		cout << " Мольные доли компонентов: " << endl;

	for (int j = 0; j < Nc; j++) {
		cout << comps[iComps[j] - 1] << "=  " << Yc[j] << "; ";
	}
	cout << endl;

	cout << "Интенсивность, масштаб турбулентности, MuT/Mu " << endl;
	cout << turb.turb_intens << ", " << turb.turb_scale << ", " << turb.MuT_Mu << endl;

}

void Region::DataEntry(Gas gas_, int iMol_, double* Yc_, Turb turb_, int Nc_) {
	gas.Lchar = gas_.Lchar;
	gas.U[0] = gas_.U[0];
	gas.U[1] = gas_.U[1];
	gas.U[2] = gas_.U[2];
	gas.p = gas_.p;
	gas.T = gas_.T;

	iMol = iMol_;

	Nc = Nc_;
	Yc = new double[Nc];
	for (int i = 0; i < Nc; i++)
		Yc[i] = Yc_[i];

	turb.turb_intens = turb_.turb_intens;
	turb.turb_scale = turb_.turb_scale;
	turb.MuT_Mu = turb_.MuT_Mu;
}

Region& Region::operator=(Region d_o)
{
	this->gas.Lchar = d_o.gas.Lchar;
	this->gas.U[0] = d_o.gas.U[0];
	this->gas.U[1] = d_o.gas.U[1];
	this->gas.U[2] = d_o.gas.U[2];
	this->gas.p = d_o.gas.p;
	this->gas.T = d_o.gas.T;

	this->iMol = d_o.iMol;

	this->Nc = d_o.Nc;


	int n = d_o.Nc;

	this->Yc = new double[n];

	for (int i = 0; i < n; i++)
		this->Yc[i] = d_o.Yc[i];

	this->turb.turb_intens = d_o.turb.turb_intens;
	this->turb.turb_scale = d_o.turb.turb_scale;
	this->turb.MuT_Mu = d_o.turb.MuT_Mu;

	return *this;
}
