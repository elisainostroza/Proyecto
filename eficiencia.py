import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# abrimos el archivo root y accedemos al árbol "Eficiencia" que contiene los datos necesarios para el mapa de eficiencia
file = uproot.open("datos.root")
tree = file["Eficiencia"]

# cargamos los datos
df = tree.arrays(library="pd")

#definimos los bins, que deben ser los rangos definidos en main.cc para pT y eta
pt_bins = np.linspace(1.0e-19, 5.5e-18, 40)
eta_bins = np.linspace(-3.0, 3.0, 40)

total, xedges, yedges = np.histogram2d(df['pT'], df['eta'], bins=[pt_bins, eta_bins])

encontrado_hist, _, _ = np.histogram2d(df[df['encontrado']]['pT'], df[df['encontrado']]['eta'], bins=[pt_bins, eta_bins]) #criterio geométrico
reco_hist, _, _ = np.histogram2d(df[df['reconstruido']]['pT'], df[df['reconstruido']]['eta'], bins=[pt_bins, eta_bins]) #criterio de pT ajustado


eff_encontrado = np.divide(encontrado_hist, total, out=np.zeros_like(encontrado_hist), where=total > 0) #de las partículas lanzadas, cuántas cruzan al menos 3 capas (criterio geométrico/aceptancia)
eff_reconstruido = np.divide(reco_hist, total, out=np.zeros_like(reco_hist), where=total > 0) #de las partículas lanzadas, cuántas terminan con un ajuste de pT correcto (criterio de calidad)

eff_ajuste_condicional = np.divide(reco_hist, encontrado_hist, out=np.zeros_like(reco_hist), where=encontrado_hist > 0) #de las que si cruzaron al menos 3 capas, cuántas terminan con un ajuste de pT correcto 


fig, axes = plt.subplots(1, 3, figsize=(20, 6.5), sharey=True)

mapas = [
    (eff_encontrado, "Eficiencia geométrica\n"),
    (eff_ajuste_condicional, "Eficiencia del ajuste\n"),
    (eff_reconstruido, "Eficiencia total de reconstrucción\n"),
]

for ax, (mapa, titulo) in zip(axes, mapas):
    im = ax.imshow(mapa.T, origin='lower',
                    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                    aspect='auto', cmap='plasma', vmin=0, vmax=1)
    ax.set_title(titulo, fontsize=11)
    ax.set_xlabel('Momento Transversal ($p_T$) [kg·m/s]')

axes[0].set_ylabel(r'Pseudorapidez ($\eta$)')
fig.colorbar(im, ax=axes, label='Eficiencia', fraction=0.025, pad=0.02)

fig.suptitle(r'Mapas de Eficiencia: $p_T$ vs $\eta$ (Protón)', fontsize=14)
plt.show()
