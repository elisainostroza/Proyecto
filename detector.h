#ifndef DETECTOR_H
#define DETECTOR_H

#include <vector>
#include "estructuras.h"
#include <TRandom3.h>

using namespace std;

//función para calcular las derivadas de posición y momento de la partícula según las ecuaciones de movimiento relativistas
Derivada calcular_derivada(const Estado& estado_actual, double m, double q, const Vec& E, const Vec& B); 

//función para avanzar el estado de la partícula usando el método de Runge-Kutta de cuarto orden
Estado rk4(const Estado& estado_actual, double dt, double m, double q, const Vec& E, const Vec& B); 

//función para detectar las intersecciones de la partícula con las capas del detector, además guarda la trayectoria completa
vector<Vec> detectar_intersecciones(Estado estado, double dt, int N, double m, double q, const Vec& E, const Vec& B, const vector<double>& capas, vector<Estado>& trayectoria_completa); 

//función para aplicar ruido gaussiano a los hits. Recibe el generador por referencia
//para que su estado avance entre llamadas y no repita siempre la misma secuencia aleatoria
vector<Vec> ruido(const vector<Vec>& hits_exactos, double sigma_xy, double sigma_z, TRandom3& generador);

//función que ajusta la hélice y calcula el momento de la particula usando mínimos cuadrados
TrackAjustado ajustar_helice(const vector<Vec>& hits_ruidosos, double B, double q);

#endif