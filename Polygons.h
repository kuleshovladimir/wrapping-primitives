#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

using namespace std;

struct Pnt {
	double x, y;
};

class Polygons{
private:
	Pnt* p;
	int n;

public:
	Polygons() ;
	~Polygons();
	Polygons(Pnt* p1, int n1) ;

	void DataEntry(string filename);
	int get_n();
	Pnt get_p(int count);
	Pnt MassCenter();
	double Square();
};

