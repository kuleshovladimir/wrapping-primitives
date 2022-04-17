#include "Functions.h"

void SplitString(string str, string& s1, string& s2, string del) {
	auto pos = str.find(del);
	if (pos != string::npos)
	{
		s1 = str.substr(0, pos);
		s2 = str.substr(pos + 1);
	}
}

void GetName(string str, string& s1, int& n, string& s2)
{
	string del;
	del = "! ";
	
	auto pos = str.find(del);
	if (pos != string::npos)
	{
		s1 = str.substr(0, pos);
		s2 = str.substr(pos + 1);
	}

	n = s1.size();

}

void Split(string str, string* (&parts), int n) {
	const int N = 256;      //Максимальная длина строки
	char word[N] = {};          //Буфер для считывания строки
	stringstream x; //Создание потоковой переменной
	x << str;                //Перенос строки в поток

	int i = 0;
	while (x >> word) {
		if (i < n) parts[i] = word;
		i++;
	}
}


void ReadRules(rules& d, string fileName) {

	string f = "Input/" + fileName;
	ifstream reading(f);

	if (reading) {
		string line;
		string name[23];

		for (int i = 0; i < 22; i++) {
			getline(reading, line);
			GetName(line, name[i], d.nameLength[i], d.comment[i]);

			//cout << "name[i]  " << name[i] << endl;
			//cout << "d.comment[i]  " << d.comment[i] << endl;
			//cout << "d.nameLength[i]  " << d.nameLength[i] << endl;

		}
		d.iDim = atoi(name[0].c_str());
		d.iAxisym = atoi(name[1].c_str());
		d.ItMin = atoi(name[2].c_str());
		d.ItMax = atoi(name[3].c_str());
		d.ResMax = atof(name[4].c_str());
		d.iExplicit = atoi(name[5].c_str());
		d.iTurb = atoi(name[6].c_str());
		d.iTrb2 = atoi(name[7].c_str());
		d.nInviscid = atoi(name[8].c_str());
		d.cMesh = name[9];
		d.iViscous = atoi(name[10].c_str());
		d.i_transform = atoi(name[11].c_str());
		d.iOutResults = atoi(name[12].c_str());
		d.iChem = atoi(name[13].c_str());
		d.iChemH2 = atoi(name[14].c_str());
		d.iChemCO = atoi(name[15].c_str());
		d.iChemNO = atoi(name[16].c_str());
		d.iChemHCl = atoi(name[17].c_str());

		d.CFLmax = atof(name[18].c_str());
		d.CFLmax_Explicit = atof(name[19].c_str());
		d.iDebug = atoi(name[20].c_str());
		d.Nc = atoi(name[21].c_str());

		d.IComps = new int[d.Nc];

		/*
				string* Comps = new string[d.Nc];
				Split(name[22], Comps, d.Nc);
				for (int i = 0; i < d.Nc; i++) {
					d.IComps[i] = atoi(Comps[i].c_str());
					cout << "d.IComps[i]= " << d.IComps[i] << endl;
				}
		*/
		for (int i = 0; i < d.Nc; i++) {
			reading >> d.IComps[i];

			//cout << "d.IComps[i]= " << d.IComps[i] << endl;
		}


		
	}
	else
		cout << "Ошибка открытия файла в ReadRules!" << endl;


	reading.close();
}


void ShowRules(rules r, string* comps) {
	string a, c, b;

	system("clear"); //очистка консоли

	cout << "ПРАВИЛА: " << endl;

	cout << "1. Размерность задачи: " << r.iDim << endl;

	c = "Нет";
	if (r.iAxisym == 1) c = "Да";
	cout << "2. Осесимметричность: " << c << endl;

	cout << "3. Мин. и макс.число итераций:: " << r.ItMin <<", " << r.ItMax << endl;

	cout << "4. Макс. невязка: " << r.ResMax << endl;

	c = "Неявная";
	if (r.iExplicit == 1) c = "Явная";
	cout << "5. Тип численной схемы: " << c << endl;

	string s[5] = { "Ламинарная","k-epsilon","k-omega","k-g","SST" };
	int i = r.iTurb;
	c = "без учета сжимаемости";
	if (r.iTrb2 == 2) c = "с учетом сжимаемости";
	cout << "6. Модель турбулентности: " << s[i] <<  ", " << c << endl;

	cout << "7. Порядок представления невязких потоков: " << r.nInviscid << endl;

	cout << "8. Сеточный файл: " << r.cMesh << endl;

	c = "Нет";
	if (r.iViscous == 1) c = "Да";
	cout << "9. Учитываются ли вязкие потоки: " << c << endl;

	a = "Задаются в программе (файл \"InitAll.txt\")";
	if (r.i_transform == 2) a = "Берутся из предыдущего расчета из файла \"Results.txt\"";
	cout << "10. Способ задания начальных условий: " << a << endl;

	c = "Нет";
	if (r.iOutResults == 1) c = "Да";
	cout << "11. Создавать файл \"Results.txt\"?: " << c << endl;

	c = "Нет";
	if (r.iChem == 1) c = "Да";
	cout << "12. Учитываются химические реакции: " << c << endl;

	string d;
	a = ""; b = ""; c = ""; d = "";
	
	if (r.iChemH2 == 1) a = "H2+O2 упрощенная, ";
	if (r.iChemH2 == 2) a = "H2+O2 полная, ";
	if (r.iChemCO == 1) b = "CO+O2 упрощенная, ";
	if (r.iChemCO == 2) b = "CO+O2 полная, ";
	if (r.iChemNO == 1) c = "N2+O2 упрощенная, ";
	if (r.iChemNO == 2) c = "N2+O2 + e, ";
	if (r.iChemNO == 3) c = "N2+O2 полная, ";
	if (r.iChemHCl == 1) c = " реакции хлора ";
	cout << "13. Хим.кинетика: " << a << b << c << d <<  endl;

	int q = r.IComps[0] - 1;
	c = comps[q];

	for (int i = 1; i < r.Nc; i++) {
		q = r.IComps[i] - 1;
		c += ", " + comps[q];
	}

	cout << "14. Химические компоненты: " << c << endl;

	cout << "+ В систему при расчете будут добавлены компоненты, " ;
	cout << "входящие в указанные выше реакции" << endl;


	cout << "15. Максимальное значение числа CFL: " << r.CFLmax << endl;
	cout << "16. Максимальное значение числа CFL в явной схеме: " << r.CFLmax_Explicit << endl;

	c = "Нет";
	if (r.iDebug == 1) c = "Да";
	cout << "17. Режим отладки: " << c << endl;
	
	cout << "------------------------------------------- "  << endl << endl;



}

void SaveRules(rules r, string fileName)
{
	// создаем поток для записи
	string f = "Input/" + fileName;
	ofstream record(f, ios::out);

	if (record) {
		record << setw(r.nameLength[0]) << left << r.iDim << "!" << r.comment[0]  << endl;
			record << setw(r.nameLength[1]) << left << r.iAxisym << "!" << r.comment[1] << endl;
			record << setw(r.nameLength[2]) << left << r.ItMin << "!" << r.comment[2] << endl;
			record << setw(r.nameLength[3]) << left << r.ItMax << "!" << r.comment[3] << endl;

			record.setf(ios::scientific);
			record.precision(2); // две цифры после запятой
			record << setw(r.nameLength[4]) << left << r.ResMax << "!" << r.comment[4] << endl;
			record.unsetf(ios::scientific);

			record << setw(r.nameLength[5]) << left << r.iExplicit << "!" << r.comment[5] << endl;
			record << setw(r.nameLength[6]) << left << r.iTurb << "!" << r.comment[6] << endl;
			record << setw(r.nameLength[7]) << left << r.iTrb2 << "!" << r.comment[7] << endl;
			record << setw(r.nameLength[8]) << left << r.nInviscid << "!" << r.comment[8] << endl;
			record << setw(r.nameLength[9]) << left << r.cMesh << "!" << r.comment[9] << endl;
			record << setw(r.nameLength[10]) << left << r.iViscous << "!" << r.comment[10] << endl;
			record << setw(r.nameLength[11]) << left << r.i_transform << "!" << r.comment[11] << endl;
			record << setw(r.nameLength[12]) << left << r.iOutResults << "!" << r.comment[12] << endl;
			record << setw(r.nameLength[13]) << left << r.iChem << "!" << r.comment[13] << endl;
			record << setw(r.nameLength[14]) << left << r.iChemH2 << "!" << r.comment[14] << endl;
			record << setw(r.nameLength[15]) << left << r.iChemCO << "!" << r.comment[15] << endl;
			record << setw(r.nameLength[16]) << left << r.iChemNO << "!" << r.comment[16] << endl;
			record << setw(r.nameLength[17]) << left << r.iChemHCl << "!" << r.comment[17] << endl;

			record.setf(ios::scientific);
			record.precision(2); // две цифры после запятой
			record << setw(r.nameLength[18]) << left << r.CFLmax << "!" << r.comment[18] << endl;
			record << setw(r.nameLength[19]) << left << r.CFLmax_Explicit << "!" << r.comment[19] << endl;
			record.unsetf(ios::scientific);

			record << setw(r.nameLength[20]) << left << r.iDebug << "!" << r.comment[20] << endl;
			record << setw(r.nameLength[21]) << left << r.Nc << "!" << r.comment[21] << endl;

			for (int i = 0; i < r.Nc; i++) {
				record << r.IComps[i] << " ";
			}
	}
	else
		cout << "Ошибка открытия файла!" << endl;

	record.close();
}

void EditRules(rules& r, int n, string* comps)
{
	int i;
	string s;

	switch (n)
	{
	case 1:
		cout << "Введите размерность задачи (2 или 3): " << endl;
		cin >> i;
		if (i == 2 || i == 3)
			r.iDim = i;
		else
			cout << "Ошибка ввода!" << endl << endl;
		break;

	case 2:
		cout << "Задача осесимметричная? (1 - Да; любое другое число - Нет): " << endl;
		cin >> i;
		if (i == 1)
			r.iAxisym = 1;
		else
			r.iAxisym = 0;
		break;

	case 3:
		cout << "Мин. и максимальное число итераций (любые числа > 1): " << endl;
		cin >> r.ItMin >> r.ItMax;
		break;

	case 4:
		cout << " Максимальная невязка: " << endl;
		cin >> r.ResMax;
		break;

	case 5:
		cout << "Какая схема используется? (1 - Явная; любое другое число - Неявная): " << endl;
		cin >> i;
		if (i == 1)
			r.iExplicit = 1;
		else
			r.iExplicit = 0;
		break;

	case 6:
		//string s[5] = { "Ламинарная","","k-omega","k-g","SST" };

		cout << "Какая Модель турбулентности используется? " << endl;
		cout << "(0 - лам.течение; 1 - k-epsilon; 2 - k-omega; 3 - k-g; любое другое число - SST): " << endl;
		cin >> i;
		if (i < 4)
			r.iTurb = i;
		else
			r.iTurb = 4;

		cout << "Учитывается сжимаемость? (1 - Да; любое другое число - Нет): " << endl;
		cin >> i;
		if (i == 1)
			r.iTrb2 = 2;
		else
			r.iTrb2 = 1;

		break;

	case 7:
		cout << "Порядок представления невязких потоков? (1 - первый; любое другое число - 2-ой): " << endl;
		cin >> i;
		if (i == 1)
			r.nInviscid = 1;
		else
			r.nInviscid = 2;
		break;

	case 8:
		cout << "Имя сеточного файла " << endl;
		cin >> s;
		if (s != "")  r.cMesh = s;

		break;

	case 9:
		cout << "Учитываются вязкие потоки? (1 - Да; любое другое число - Нет): " << endl;
		cin >> i;
		if (i == 1)
			r.iViscous = 1;
		else
			r.iViscous = 0;

		break;

	case 10:
		cout << "Какой способ задания начальных условий?" << endl;
		cout << "(0 - Задаются в программе; любое другое число - Берутся из предыдущего расчета из файла \"Results.txt\": " << endl;
		cin >> i;
		if (i == 0)
			r.i_transform = i;
		else
			r.i_transform = 2;

		break;

	case 11:
		cout << "Создавать файл \"Results.txt\"? (1 - Да; любое другое число - Нет) : " << endl;
		cin >> i;
		if (i == 1)
			r.iOutResults = 1;
		else
			r.iOutResults = 0;

		break;

	case 12:
		cout << "Учитывать химические реакции? (1 - Да; любое другое число - Нет) : " << endl;
		cin >> i;
		if (i == 1)
			r.iChem = 1;
		else
			r.iChem = 0;

		break;

	case 13:

		cout << "Какая схема для H2+O2?  : " << endl;
		cin >> r.iChemH2;
		cout << "Какая схема для CO+O2?  : " << endl;
		cin >> r.iChemCO;
		cout << "Какая схема для N2+O2?  : " << endl;
		cin >> r.iChemNO;
		cout << "Какая схема для HCl?  : " << endl;
		cin >> r.iChemHCl;

		break;

	case 14:
		cout << "Число хим. компонентов  : " << endl;
		cin >> r.Nc;
		delete(r.IComps);
		r.IComps = new int[r.Nc];

		cout << "Выберите номера хим. компонентов из списка  : " << endl;
		for (int i = 0; i < 30; i++)
		{
			if (comps[i] != "ZZ")
				cout << i + 1 << " -> " << comps[i] << endl;
		}
		int j;
		for (int i = 0; i < r.Nc; i++)
		{
			cin >> j;
			r.IComps[i] = j;
		}

		break;

	case 15:
		cout << "Максимальное значение числа CFL: " << endl;
		cin >> r.CFLmax;
		break;

	case 16:
		cout << "Максимальное значение числа CFL в явной схеме: " << endl;
		cin >> r.CFLmax_Explicit;
		break;

	case 17:
		cout << "Режим отладки? (1 - Да; любое другое число - Нет) : " << endl;
		cin >> i;
		if (i == 1)
			r.iDebug = 1;
		else
			r.iDebug = 0;
		break;

	default:
		break;
	}

}

void ReadParams(Region* (&d), int& n, int Nc, string fileName)
{
	string reg;
	string f = "Input/" + fileName;
	ifstream reading(f);


	if (reading) {

		reading >> n;

		// временные переменные
		Gas gas;
		Turb turb;
		int iMol;
		double* Yc = new double[Nc];
		//cout << "Nc = " << Nc << endl;
		//  выделяем память
		d = new Region[n];

		for (int i = 0; i < n; i++) {
			reading >> reg;
			reading >> gas.Lchar >> gas.U[0] >> gas.U[1] >> gas.U[2] >> gas.T >> gas.p;
			reading >> iMol;

			for (int j = 0; j < Nc; j++)
				reading >> Yc[j];

			//cout << "Yc[0] = " << Yc[0] << ", Yc[1] = " << Yc[1] << endl;

			reading >> turb.turb_intens >> turb.turb_scale >> turb.MuT_Mu;

			d[i].DataEntry(gas, iMol, Yc, turb, Nc);
		}

		// cout << "Данные считаны in ReadParams!" << endl;
		//cout << endl;
	}
	else {
		cout << "Ошибка открытия файла в ReadParams" << endl;
	}

	reading.close();
}

void ShowParams(Region* z, int n, int* iComps, string* comps)
{
	cout << "*****************************" << endl;
	
	for (int i = 0; i < n; i++) {
		cout << "       Участок № " << i + 1 << ":" << endl;

		z[i].Print(iComps, comps);
		cout << "___________________________" << endl;

	}
	cout << " " << endl;

}

void SavingData(Region* d, int n, string fileName)
{
	// создаем поток для записи
	string f = "Input/" + fileName;
	ofstream record(f, ios::out);

	if (record) {
		record << n << endl;

		for (int i = 0; i < n; i++) {
			record << "region" << i + 1 << endl;
			record << d[i].GetGas().Lchar << " " << d[i].GetGas().U[0] << " " << d[i].GetGas().U[1]
				<< " " << d[i].GetGas().U[2] << " " << d[i].GetGas().T << " " << d[i].GetGas().p << endl;
			record << d[i].GetiMol() << endl;

			//Yc
			for (int j = 0; j < d[i].GetNc(); j++) {
				record << d[i].GetYc(j) << " ";
			}
			record << endl;

			record << d[i].GetTurb().turb_intens << " " << d[i].GetTurb().turb_scale
				<< " " << d[i].GetTurb().MuT_Mu << endl;
		}

		cout << "Данные записаны!" << endl;
	}
	else
		cout << "Ошибка открытия файла!" << endl;

	record.close();
}

void Copy(Region* (&d_n), Region* (&d_o), int n) {
	for (int i = 0; i < n; i++) {
		d_n[i] = d_o[i];
	}
}

void AddData(Region* (&z), int& m, int Nc, int* iComps, string* comps)
{

	// временные переменные
	Gas gas;
	Turb turb;
	int iMol;

	double* Yc = new double[Nc];

	// временный массив
	Region* buf = new Region[m];


	Copy(buf, z, m);

	//// Выделяем новую память
	z = new Region[m + 1];

	//z[m].GetYc() = new double[Nc];

	Copy(z, buf, m);
	cout << "Введите характерный размер [м]:  ";
	cin >> gas.Lchar;

	cout << "Введите компоненты скорости u, v, w [м/с]:  ";
	cin >> gas.U[0];
	cin >> gas.U[1];
	cin >> gas.U[2];
	cout << "Введите температуру [K]:  ";
	cin >> gas.T;
	cout << "Введите давление [Па]:  ";
	cin >> gas.p;

	int f;
	cout << "Используются мольные (1) или массовые доли (0) компонентов?  ";
	cin >> f;
	if (f == 0)
		iMol = 0;
	else
		iMol = 1;

	for (int j = 0; j < Nc; j++) {
		if (iMol == 0)
			cout << "Введите массовую долю " << j + 1 << "-го компонента:   ";
		else
			cout << "Введите мольную долю " << j + 1 << "-го компонента:   ";

		cin >> Yc[j];
	}

	cout << "Введите интенсивность турбулентности (в долях от скорости):  ";
	cin >> turb.turb_intens;

	cout << "Введите масштаб турбулентности (в долях от характ.размера):  ";
	cin >> turb.turb_scale;

	cout << "Введите MuT/Mu:  ";
	cin >> turb.MuT_Mu;

	z[m].DataEntry(gas, iMol, Yc, turb, Nc);


	m++;

	system("clear");
	delete[] buf;

	cout << "Данные добавлены!" << endl;
	system("pause");
}
