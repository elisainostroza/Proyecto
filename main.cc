#include <iostream>
#include <vector>

#include "estructuras.h"
#include "detector.h"
#include "particula.h"

#include <TFile.h>
#include <TTree.h>

using namespace std;

int main() {

    Particula mi_particula = crear_proton(); //probamos con un proton, pero se pueden probar con otras partículas cambiando esta línea por crear_pion_mas() o crear_kaon_mas()

    Vec E = {0.0, 0.0, 0.0};
    Vec B = {0.0, 0.0, 3.8}; //definimos un campo magnético uniforme en la dirección z
    
    Estado estado0;
    estado0.r = {0.0, 0.0, 0.0}; //generamos la partícula en el origen
    estado0.p = {1e-19, 0.0, 2e-20}; 

    vector<double> capas = {0.05, 0.10, 0.20, 0.40, 0.60}; //radios del detector en metros
    double dt = 1e-11;
    int N = 20000;

    vector<Estado> trayectoria; //vector para guardar la trayectoria completa de la partícula
    
    vector<Vec> intersecciones = detectar_intersecciones(estado0, dt, N, mi_particula.masa, mi_particula.carga, E, B, capas, trayectoria); // ejecutamos la simuación

    vector<Vec> hits_reales = ruido(intersecciones, 0.0001, 0.0001); //agregamos ruido gaussiano a los hits exactos para simular la resolución del detector, con una sigma de 100 micrómetros que es la resolución típica de  los sensores de silicio de ATLAS

    TrackAjustado track = ajustar_helice(hits_reales, B.z, mi_particula.carga);//reconstruimos el momento

./

    // Exportamos los datos a root
    TFile *archivo_salida = new TFile("datos.root", "RECREATE");

    // trayectoria teórica completa
    TTree *arbol_trayectoria = new TTree("Trayectoria", "Puntos continuos del RK4");
    double x_out, y_out, z_out;
    
    arbol_trayectoria->Branch("x", &x_out, "x/D");
    arbol_trayectoria->Branch("y", &y_out, "y/D");
    arbol_trayectoria->Branch("z", &z_out, "z/D");
    
    for(const auto& est : trayectoria){
        x_out = est.r.x;
        y_out = est.r.y;
        z_out = est.r.z;
        arbol_trayectoria->Fill(); 
    }

    // HITS
    TTree *arbol_hits = new TTree("Hits", "Intersecciones con las capas del detector");
    double x_ex, y_ex, z_ex;       // exactos (sin ruido)
    double x_ruido, y_ruido, z_ruido; // con resolucion gaussiana
    
    arbol_hits->Branch("x_exacto", &x_ex, "x_exacto/D");
    arbol_hits->Branch("y_exacto", &y_ex, "y_exacto/D");
    arbol_hits->Branch("z_exacto", &z_ex, "z_exacto/D");
    
    arbol_hits->Branch("x_ruido", &x_ruido, "x_ruido/D");
    arbol_hits->Branch("y_ruido", &y_ruido, "y_ruido/D");
    arbol_hits->Branch("z_ruido", &z_ruido, "z_ruido/D");
    
    for(size_t i = 0; i < intersecciones.size(); i++){ 
        x_ex = intersecciones[i].x; y_ex = intersecciones[i].y; z_ex = intersecciones[i].z;
        x_ruido = hits_reales[i].x; y_ruido = hits_reales[i].y; z_ruido = hits_reales[i].z;
        arbol_hits->Fill();
    }

    // reconstrucción de la hélice
    TTree *arbol_reco = new TTree("Reconstruccion", "Parametros ajustados de la helice");
    double pT_real, pz_real;
    double pT_ajustado, pz_ajustado, radio_ajustado;
    
    arbol_reco->Branch("pT_real", &pT_real, "pT_real/D");
    arbol_reco->Branch("pz_real", &pz_real, "pz_real/D");
    arbol_reco->Branch("pT_ajustado", &pT_ajustado, "pT_ajustado/D");
    arbol_reco->Branch("pz_ajustado", &pz_ajustado, "pz_ajustado/D");
    arbol_reco->Branch("radio_ajustado", &radio_ajustado, "radio_ajustado/D");
    
    // valores reales de la partícula
    pT_real = sqrt(estado0.p.x*estado0.p.x + estado0.p.y*estado0.p.y);
    pz_real = estado0.p.z;
    
    // valores ajustados de la reconstrucción
    pT_ajustado = track.pT;
    pz_ajustado = track.pz;
    radio_ajustado = track.radio;
    
    arbol_reco->Fill();

    archivo_salida->Write();
    archivo_salida->Close();
    
    cout << "Simulación completada" << endl;

    return 0;
}