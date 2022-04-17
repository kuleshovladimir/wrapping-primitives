#pragma once

#include <iostream>
#include <fstream>
#include <string>


using namespace std;

struct Gas {
	double Lchar;
	double U[3];
	double T;
	double p;
};

struct Turb {
	double turb_intens, turb_scale, MuT_Mu;
};

class Region {
	//закрытая часть класса
private: 
	Gas gas;
	Turb turb;
	int iMol;
	double* Yc;
	int Nc;

	//открытая часть класса
public:
	// обязательные методы
	Region(); // конструктор по умолчанию. Создает пустой объект
	Region(Gas gas_, int iMol_, double* Yc_, Turb turb_, int Nc_); // конструктор с параметрами. Создает объект с параметрами
	~Region(); // деструктор. Удаляет объект 

	// необязательные методы (нужны для выполнения задачи)
	void Print(int* iComps, string* comps);
	void DataEntry(Gas gas_, int iMol_, double* Yc_, Turb turb_, int Nc_);

	// перегрузка оператора =
	Region& operator = (Region d_o);

	Gas GetGas() { return gas; };
	Turb GetTurb() { return turb; };
	int GetiMol() { return iMol; };
	double GetYc(int i) { return Yc[i]; };
	int GetNc() { return Nc; };


};
