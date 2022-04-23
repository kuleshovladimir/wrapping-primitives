#include "Polygons.h"

Polygons::Polygons()
{
	n = 0;
	p = new Pnt[0];
}

Polygons::Polygons(Pnt* p_in, int n_in)
{
	n = n_in;
	p = new Pnt[n];

	for (int i = 0; i < n; i++) {
		p[i] = p_in[i];
	}
}

void Polygons::DataEntry(string filename)
{
	ifstream reading(filename);

	if (reading) {
		cout << "Файл points открыт\n";

		reading >> n;
		p = new Pnt[n];
		for (int i = 0; i < n; i++) {
			reading >> p[i].x >> p[i].y;
		}
	}
	else {
		cout << "Не удалось открыть файл points\n";
	}
}

int Polygons::get_n()
{
	return n;
}

Pnt Polygons::get_p(int count)
{
	return p[count];
}

Pnt Polygons::MassCenter()
{
	Pnt* xn = new Pnt[n + 1];

	for (int i = 0; i < n; i++) {
		xn[i].x = p[i].x;
		xn[i].y = p[i].y;
	}
	xn[n].x = xn[0].x;
	xn[n].y = xn[0].y;

	double A = 0;
	double Cx = 0, Cy = 0;

	for (int i = 0; i < n; i++) {
		A = A + xn[i].x * xn[i + 1].y - xn[i + 1].x * xn[i].y;
		Cx = Cx + (xn[i].x + xn[i + 1].x) * (xn[i].x * xn[i + 1].y - xn[i + 1].x * xn[i].y);
		Cy = Cy + (xn[i].y + xn[i + 1].y) * (xn[i].x * xn[i + 1].y - xn[i + 1].x * xn[i].y);
	}
	A = A / 2;

	Pnt ret;
	ret.x = Cx / (6 * A);
	ret.y = Cy / (6 * A);

	
	return ret;
}

double Polygons::Square()
{
	Pnt* xn = new Pnt[n + 1];

	for (int i = 0; i < n; i++) {
		xn[i].x = p[i].x;
		xn[i].y = p[i].y;
	}
	xn[n].x = xn[0].x;
	xn[n].y = xn[0].y;

	double A = 0;
	double Cx = 0, Cy = 0;

	for (int i = 0; i < n; i++) {
		A = A + xn[i].x * xn[i + 1].y - xn[i + 1].x * xn[i].y;
	}
	A = fabs(A / 2);

	return (A);
}

Polygons::~Polygons()
{

}