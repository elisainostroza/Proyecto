#include <iostream>
#include <vector>

#include "estructuras.h"
#include "detector.h"

#include <TFile.h>
#include <TTree.h>

using namespace std;

int main() {
    // Definimos un protón para verificar que todo funcione correctamente
    double masa_proton = 1.6726219e-27;     
    double carga_proton = 1.602176634e-19;   

    Vec E = {0.0, 0.0, 0.0};
    Vec B = {0.0, 0.0, 3.8}; //Definimos un campo magnético uniforme en la dirección z
    
    Estado estado0;
    estado0.r = {0.0, 0.0, 0.0}; //generamos la partícula en el origen
    estado0.p = {1e-19, 0.0, 2e-20}; 

    vector<double> capas = {0.05, 0.10, 0.20, 0.40, 0.60}; // radios de las capas del detector en metros
    double dt = 1e-11;
    int N = 20000;

    vector<Estado> trayectoria; //vector para guardar la trayectoria completa de la partícula
    
    // Ejecutamos la simulación
    vector<Vec> intersecciones = detectar_intersecciones(estado0, dt, N, masa_proton, carga_proton, E, B, capas, trayectoria);

    // Exportamos a archivos ROOT
    TFile *archivo_salida = new TFile("datos.root", "RECREATE");
    TTree *arbol = new TTree("TrackTree", "Trayectoria y cruces cilindricos");
    
    double x_out, y_out, z_out;
    bool es_hit;
    
    arbol->Branch("x", &x_out, "x/D");
    arbol->Branch("y", &y_out, "y/D");
    arbol->Branch("z", &z_out, "z/D");
    arbol->Branch("es_hit", &es_hit, "es_hit/O");
    
    for(const auto& est : trayectoria){
        x_out = est.r.x;
        y_out = est.r.y;
        z_out = est.r.z;
        es_hit = false; 
        
        for(const auto& inter : intersecciones){
            if(est.r.x == inter.x && est.r.y == inter.y && est.r.z == inter.z){
                es_hit = true;
                break;
            }
        }
        arbol->Fill(); 
    }
    
    archivo_salida->Write();
    archivo_salida->Close();

    return 0;
}