#include "detector.h"
#include <cmath>

const double c = 299792458.0; //velocidad de la luz en el vacío

Derivada calcular_derivada(const Estado& estado_actual, double m, double q, const Vec& E, const Vec& B){
    Derivada d;
    double p2 = norma2(estado_actual.p);
    double gamma = sqrt(1.0 + p2 / (m*m*c*c));
    Vec v = estado_actual.p * (1.0/(gamma*m));

    d.drdt = v; 
    d.dpdt = (E + cruz(v, B)) * q;
    return d;
}

Estado rk4(const Estado& estado_actual, double dt, double m, double q, const Vec& E, const Vec& B) { //implementamos el método de Runge-Kutta de cuarto orden para avanzar el estado de la partícula
    Derivada k1 = calcular_derivada(estado_actual, m, q, E, B);

    Estado s2;
    s2.r = estado_actual.r + k1.drdt * (dt/2.0);
    s2.p = estado_actual.p + k1.dpdt * (dt/2.0);
    Derivada k2 = calcular_derivada(s2, m, q, E, B);

    Estado s3;
    s3.r = estado_actual.r + k2.drdt * (dt/2.0);
    s3.p = estado_actual.p + k2.dpdt * (dt/2.0);
    Derivada k3 = calcular_derivada(s3, m, q, E, B);

    Estado s4;
    s4.r = estado_actual.r + k3.drdt * dt;
    s4.p = estado_actual.p + k3.dpdt * dt;
    Derivada k4 = calcular_derivada(s4, m, q, E, B);

    Estado siguiente;
    siguiente.r = estado_actual.r + (dt/6.0) * (k1.drdt + 2.0*k2.drdt + 2.0*k3.drdt + k4.drdt);
    siguiente.p = estado_actual.p + (dt/6.0) * (k1.dpdt + 2.0*k2.dpdt + 2.0*k3.dpdt + k4.dpdt);

    return siguiente;
}

vector<Vec> detectar_intersecciones(Estado estado, double dt, int N, double m, double q, const Vec& E, const Vec& B, const vector<double>& capas, vector<Estado>& trayectoria_completa) {
    vector<Vec> intersecciones;
    vector<bool> cruzada(capas.size(), false); //lista para no guardar múltiples intersecciones con la misma capa

    double rho_antes = sqrt(estado.r.x*estado.r.x + estado.r.y*estado.r.y); //calcula donde está la partícula antes de moverse
    
    trayectoria_completa.push_back(estado);

    for(int i = 0; i < N; i++){
        estado = rk4(estado, dt, m, q, E, B);
        trayectoria_completa.push_back(estado); //guardamos la trayectoria de la partícula
        
        double rho = sqrt(estado.r.x*estado.r.x + estado.r.y*estado.r.y);
        for(int j = 0; j < capas.size(); j++){

            if(!cruzada[j] && rho_antes < capas[j] && rho >= capas[j]){ //verificamos si en este paso se cruzó la capa j y si no se había cruzado antes
                intersecciones.push_back(estado.r); //guardamos la posición de la intersección
                cruzada[j] = true; //guardamos que esta capa ya fue cruzada para no contarla de nuevo
            }
        }
        rho_antes = rho;
    }
    return intersecciones;
}