// вынесли из-под скобок
if (face.zone == 1) {
    //Fv - поток через грань 
    double Fv = 100.;
    du[c].du[0] += Fv * length / S *dt;
}
// вынесли из-под скобок

if (fase.zone == 3) {
    //Fv - поток через грань 
    double Fv = -100.;
    du[c].du[0] += Fv * length / S *dt;
}

