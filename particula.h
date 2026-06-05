#ifndef PARTICULA_H
#define PARTICULA_H
#include <string>

using namespace std;

// estructura para representar partículas comunes con su masa y carga
struct Particula {
    string nombre;
    double masa;   //kg
    double carga;  //c
};

// funciones para crear partículas comunes
Particula crear_proton();
Particula crear_pion_mas(); 
Particula crear_kaon_mas();

#endif