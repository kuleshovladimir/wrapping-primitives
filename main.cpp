#include "Functions.h"

int _stateMenu, _stateMenu1, _stateMenu2;

void Menu() {
    cout << "Выберите действие: " << endl
        << "(0) Выход из программы" << endl
        << "(1) Работа с файлом \"rules.txt\"" << endl
        << "(2) Работа с файлом \"InitAll.txt\"" << endl
        << "(3) Работа с файлом \"zones.txt\"" << endl
        << "(4) Работа с файлом \"InitCond.txt\" (нач.условия)" << endl
        << "Ваш выбор: ";
    cin >> _stateMenu;
}

void Menu1() {
    cout << "Выберите действие: " << endl
        << "  Продолжение (0)" << endl
        << "  Редактирование этого файла (укажите № пункта)" << endl
        << "  Повторный просмотр этого файла (20)" << endl
        << "  Запись этого файла (30)" << endl
        << "Ваш выбор: ";
    cin >> _stateMenu1;
}

void Menu2() {
    cout << "Выберите действие: " << endl
        << "  Продолжение (0)" << endl
        << "  Редактирование этого файла (укажите № участка)" << endl
        << "  Добавить участок (11)" << endl
        << "  Удалить участок (22)" << endl
        << "  Повторный просмотр этого файла (20)" << endl
        << "  Запись этого файла (30)" << endl
        << "Ваш выбор: ";
    cin >> _stateMenu2;
}

int main()
{
    setlocale(LC_ALL, "Russian");

    string* comps = new string[30]{
        "H",     "O",  "OH",  "H2",  "O2",   "H2O", "ZZ", "N2", "CO", "CO2",
        "ZZ",    "Ar", "ZZ",  "HO2", "H2O2", "HCl", "ZZ", "ZZ", "NO", "N",
        "C7H16", "C",  "N2+", "O2+", "NO+",  "N+",  "O+", "e",  "Fe", "FeO"
    };

    rules r;
    string fileName;

    ReadRules(r, "rules.txt");

    int  amountOfRegions = 0;
    Region* reg = new Region[amountOfRegions];

    ReadParams(reg, amountOfRegions, r.Nc, "InitAll.txt");


    //cout << r.CFLmax << endl;

    //exit(9);

    //cout << "размерность = " << r.iDim << endl;
    //cout << "iChemH2 = " << r.iChemH2 << endl;

    Menu();

    //cout << "_stateMenu = " << _stateMenu << endl;
    while (_stateMenu != 0) {



        switch (_stateMenu)
        {
        case 1:   //  Работа с файлом "rules.txt" 
                //  ????
            ShowRules(r, comps);

            Menu1();
            while (_stateMenu1 != 0) {
                if (_stateMenu1 == 0) {

                }
                else
                if (_stateMenu1 == 20) {
                    system("cls");
                    ShowRules(r, comps);
                    Menu1();
                }
                else
                if (_stateMenu1 == 30) {
                    int i;
                    system("cls");
                    cout << "Записать в старый файл (1) или создать новый (любое число)? ";
                    cin >> i;

                    if (i == 1) {
                        fileName = "rules.txt";
                    }
                    else {
                        cout << "Введите имя файла: ";
                        cin >> fileName;
                    }
                    SaveRules(r, fileName);

                    ShowRules(r, comps);
                    Menu1();
                }
                else
                if (_stateMenu1 >= 1 && _stateMenu1 <= 17) {
                        EditRules(r, _stateMenu1, comps);
                        system("cls");
                        ShowRules(r, comps);
                        Menu1();
                }

            }

            system("pause"); // задержка консоли

            system("cls"); //очистка консоли
            Menu();
            break;
        case 2:   //  Работа с файлом "InitAll.txt" 
                //  ????
            system("cls"); //очистка консоли

            ShowParams(reg, amountOfRegions, r.IComps, comps);

            Menu2();
            while (_stateMenu2 != 0)
            {
                if (_stateMenu2 == 0) {

                }
                else
                    if (_stateMenu2 == 20) {
                        system("cls"); //очистка консоли
                        ShowParams(reg, amountOfRegions, r.IComps, comps);
                        Menu2();
                    }
                    else

                    if (_stateMenu2 == 30) {
                            system("cls"); //очистка консоли
                            int i;
                            cout << "Записать в старый файл (1) или создать новый (любое число)? ";
                            cin >> i;

                            if (i == 1) {
                                fileName = "InitAll.txt";
                            }
                            else {
                                cout << "Введите имя файла: ";
                                cin >> fileName;
                            }

                            SavingData(reg, amountOfRegions, fileName);
                            system("pause"); // задержка консоли
                            system("cls"); //очистка консоли
                            ShowParams(reg, amountOfRegions, r.IComps, comps);
                            Menu2();

                     }
                     else
                     if (_stateMenu2 == 11) {
                                AddData(reg, amountOfRegions, r.Nc, r.IComps, comps);
                                system("cls"); //очистка консоли
                                ShowParams(reg, amountOfRegions, r.IComps, comps);
                                Menu2();
                      }
                     else

                     if (_stateMenu2 == 22) {
                                    //DeleteData(reg, amountOfRegions);
                                    system("pause"); // задержка консоли
                                    system("cls"); //очистка консоли

                                    ShowParams(reg, amountOfRegions, r.IComps, comps);
                                    Menu2();
                      }
                     else
                     {
                                    if (_stateMenu2 > amountOfRegions) {

                                        system("cls"); //очистка консоли
                                        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                                        cout << " Неверно введен номер действия!!" << endl;
                                        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                                        ShowParams(reg, amountOfRegions, r.IComps, comps);
                                        Menu2();
                                    }
                                    else {
                                        //EditRegions(reg, _stateMenu2);
                                        system("cls"); //очистка консоли
                                        ShowParams(reg, amountOfRegions, r.IComps, comps);
                                        Menu2();
                                    }

                     }

            }

            Menu();
            break;
        default:
            system("cls"); //очистка консоли
            Menu();
            break;
        }

    }


}


