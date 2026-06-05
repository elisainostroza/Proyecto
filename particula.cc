#include "particula.h"

// carga elemental
const double e_carga = 1.602176634e-19; 

Particula crear_proton() { //creamos un protón con su masa y carga
    Particula p;
    p.nombre = "Proton";
    p.masa = 1.6726219e-27; 
    p.carga = e_carga;
    return p;
}

Particula crear_pion_mas() { //creamos un pion positivo con su masa y carga
    Particula p;
    p.nombre = "Pion+";
    p.masa = 2.48832e-28;   
    p.carga = e_carga;
    return p;
}

Particula crear_kaon_mas() { //creamos un kaon positivo con su masa y carga
    Particula p;
    p.nombre = "Kaon+";
    p.masa = 8.8033e-28;   
    p.carga = e_carga;
    return p;
}