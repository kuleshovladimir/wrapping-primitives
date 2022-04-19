void Tecplot(parameters* p, Cell* cells, int Nx, int Ny, int nCells) 
{
    // создаем поток для записи
    string f = "Tecplot/T.plt";
    ofstream record(f, ios::out);
    if (record) {
        record << "VARIABLES = \"X\", \"Y\", \"T,K\"" << endl;

        record << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;
        for (int i = 0; i < nCells; i++) {
            record << cells[i].Get_c().x << " " << cells[i].Get_c().y << " " << p[i].T << endl;
        }

    }
    record.close();
} 