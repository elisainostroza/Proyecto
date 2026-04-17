#ifndef DETECTOR_H
#define DETECTOR_H

#include <vector>
#include "estructuras.h"

using namespace std;

Derivada calcular_derivada(const Estado& estado_actual, double m, double q, const Vec& E, const Vec& B);

Estado rk4(const Estado& estado_actual, double dt, double m, double q, const Vec& E, const Vec& B);

vector<Vec> detectar_intersecciones(Estado estado, double dt, int N, double m, double q, const Vec& E, const Vec& B, const vector<double>& capas, vector<Estado>& trayectoria_completa);

#endif