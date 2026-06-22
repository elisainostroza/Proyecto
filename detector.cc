#include "detector.h"
#include <cmath>
#include <TRandom3.h>

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

        if (std::abs(estado.r.z) > 1.5) break; //si la partícula sale del detector, terminamos la simulación

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


vector<Vec> ruido(const vector<Vec>& hits_exactos, double sigma_xy, double sigma_z, TRandom3& generador) { //aplicamos ruido gaussiano a los hits exactos para simular la resolución del detector
    vector<Vec> hits_gauss;

    for (const auto& hit : hits_exactos) {
        Vec hit_ruido;
        //agregamos el error aleatorio en cada coordenada
        hit_ruido.x = hit.x + generador.Gaus(0, sigma_xy);
        hit_ruido.y = hit.y + generador.Gaus(0, sigma_xy);
        hit_ruido.z = hit.z + generador.Gaus(0, sigma_z);
        hits_gauss.push_back(hit_ruido);
    }
    return hits_gauss;
}

TrackAjustado ajustar_helice(const vector<Vec>& hits, double B, double q) { //ajustamos la trayectoria a una hélice usando mínimos cuadrados para reconstruir el momento de la partícula
    TrackAjustado track;
    int n = hits.size();
    
    if (n < 3) { // si hay menos de 3 hits, no se puede ajustar un círculo
        track.pT = 0; track.pz = 0; track.p_total = 0; 
        return track;
    }

    double Sxx=0, Sxy=0, Sx=0, Syy=0, Sy=0, Sxw=0, Syw=0, Sw=0; 
    
    for(const auto& h : hits) { //ajuste de círculo en el plano XY para obtener el radio y centro de la hélice
        double x = h.x;
        double y = h.y;
        double w = x*x + y*y;
        
        Sxx += x*x; Sxy += x*y; Sx += x;
        Syy += y*y; Sy += y;
        Sxw += x*w; Syw += y*w; Sw += w;
    }

    //resolvemos el sistema usando Cramer para obtener el círculo ajustado
    double detM = Sxx*(Syy*n - Sy*Sy) - Sxy*(Sxy*n - Sy*Sx) + Sx*(Sxy*Sy - Syy*Sx);
    double detA = Sxw*(Syy*n - Sy*Sy) - Sxy*(Syw*n - Sw*Sy) + Sx*(Syw*Sy - Syy*Sw);
    double detB = Sxx*(Syw*n - Sw*Sy) - Sxw*(Sxy*n - Sy*Sx) + Sx*(Sxy*Sw - Syw*Sx);
    double detC = Sxx*(Syy*Sw - Sy*Syw) - Sxy*(Sxy*Sw - Syw*Sx) + Sxw*(Sxy*Sy - Syy*Sx);

    double coefA = detA / detM;
    double coefB = detB / detM;
    double coefC = detC / detM;

    track.centro_x = coefA / 2.0;
    track.centro_y = coefB / 2.0;
    track.radio = sqrt(coefC + track.centro_x*track.centro_x + track.centro_y*track.centro_y);

    track.pT = std::abs(q) * std::abs(B) * track.radio; // calculamos el momento transversal usando pT = |q||B|r

    // ajuste Lineal en el plano Z vs S
    double sum_s = 0, sum_z = 0, sum_sz = 0, sum_s2 = 0;
    double s_acumulado = 0.0;
    
    for(int i = 0; i < n; i++) { 
        if(i > 0) {
            double dx = hits[i].x - hits[i-1].x; //calculamos la distancia entre hits consecutivos para obtener s
            double dy = hits[i].y - hits[i-1].y; 
            s_acumulado += sqrt(dx*dx + dy*dy);
        }
        
        double z = hits[i].z; //coordenada z del hit
        sum_s += s_acumulado;
        sum_z += z;
        sum_sz += s_acumulado * z;
        sum_s2 += s_acumulado * s_acumulado;
    }

    // pendiente dz/ds
    double m_pendiente = (n * sum_sz - sum_s * sum_z) / (n * sum_s2 - sum_s * sum_s);
    
    // calculamos pz y p_total
    track.pz = track.pT * m_pendiente;
    track.p_total = sqrt(track.pT*track.pT + track.pz*track.pz);

    return track;
}