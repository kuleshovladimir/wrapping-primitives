#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>      // std::setw

#include "Region.h"

using namespace std;

struct rules {
	int iDim;
	int iAxisym;
	int ItMin, ItMax;
	double ResMax;
	int iExplicit;
	int iTurb, iTrb2;
	int nInviscid;
	string cMesh;
	int iViscous;
	int i_transform;
	int iOutResults;
	int iChem;
	int iChemH2, iChemCO, iChemNO, iChemHCl;

	int Nc;		// число химических компонентов
	int* IComps;  // использумые хим.компоненты

	double CFLmax, CFLmax_Explicit;

	int iDebug;

	string comment[22];
	int nameLength[22];

};


void SplitString(string str, string& s1, string& s2, string del);
void ReadRules(rules& d, string fileName);
void ShowRules(rules r, string* comps);
void SaveRules(rules r, string fileName);
void EditRules(rules& r, int n, string* comps);

void ReadParams(Region* (&d), int& n, int Nc, string fileName);
void ShowParams(Region* z, int n, int* iComps, string* comps);
void SavingData(Region* d, int n, string fileName);
void AddData(Region* (&z), int& m, int Nc, int* iComps, string* comps);