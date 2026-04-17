#ifndef DETECTOR_H
#define DETECTOR_H

#include <vector>
#include "estructuras.h"

using namespace std;

//función para calcular las derivadas de posición y momento de la partícula según las ecuaciones de movimiento relativistas
Derivada calcular_derivada(const Estado& estado_actual, double m, double q, const Vec& E, const Vec& B); 

//función para avanzar el estado de la partícula usando el método de Runge-Kutta de cuarto orden
Estado rk4(const Estado& estado_actual, double dt, double m, double q, const Vec& E, const Vec& B); 

//función para detectar las intersecciones de la partícula con las capas del detector, además guarda la trayectoria completa
vector<Vec> detectar_intersecciones(Estado estado, double dt, int N, double m, double q, const Vec& E, const Vec& B, const vector<double>& capas, vector<Estado>& trayectoria_completa); 

#endif