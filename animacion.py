import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

archivo_root = "datos.root"

capas = [0.05, 0.10, 0.20, 0.40, 0.60] #radios de las capas

salto = 20        
intervalo = 10
fps = 30

archivo = uproot.open(archivo_root) #cargamos datos de la trayectoria y los hits

tray = archivo["Trayectoria"]
x = tray["x"].array(library="np")[::salto]
y = tray["y"].array(library="np")[::salto]
z = tray["z"].array(library="np")[::salto]


hits = archivo["Hits"] #guardamos los hits con ruido para graficarlos en la animación

hx_ruido = hits["x_ruido"].array(library="np")
hy_ruido = hits["y_ruido"].array(library="np")
hz_ruido = hits["z_ruido"].array(library="np")

n_frames = len(x)


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

margen_xy = 1.15 * max(np.max(np.abs(x)), np.max(np.abs(y)), max(capas)) #margen para que la animacion se vea bien sin importar la escala de la trayectoria o las capas del detector
z_min = min(0.0, np.min(z))
z_max = 1.15 * np.max(z) if np.max(z) > 0 else 1.0

ax.set_xlim(-margen_xy, margen_xy)
ax.set_ylim(-margen_xy, margen_xy)
ax.set_zlim(z_min, z_max)
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')
ax.set_title('Trayectoria Helicoidal de Protón', fontsize=14)

# añadimos las capas del detector y los hits con ruido para dar contexto visual 
theta = np.linspace(0, 2 * np.pi, 60)
z_cil = np.linspace(z_min, z_max, 2)
theta_grid, z_grid = np.meshgrid(theta, z_cil)
for r in capas:
    x_cil = r * np.cos(theta_grid)
    y_cil = r * np.sin(theta_grid)
    ax.plot_surface(x_cil, y_cil, z_grid, alpha=0.06, color='gray', linewidth=0)

ax.scatter(hx_ruido, hy_ruido, hz_ruido, color='blue', s=30, marker='o', label='Hits con ruido', depthshade=False)


linea, = ax.plot([], [], [], lw=2, color='magenta', label='Trayectoria RK4') #la trayectoria reconstruida con RK4 se actualiza en cada frame
punto, = ax.plot([], [], [], 'o', color='green', markersize=7) #el punto que marca la posición de la partícula actual también se actualiza en cada frame

ax.legend(loc='upper right', fontsize=8)


def inicializar(): #inicializamos la línea y el punto vacíos para el primer frame
    linea.set_data([], [])
    linea.set_3d_properties([])
    punto.set_data([], [])
    punto.set_3d_properties([])
    return linea, punto


def actualizar(num): #actualizamos la línea y el punto para el frame actual, mostrando la trayectoria hasta ese punto y la posición actual de la partícula
    linea.set_data(x[:num], y[:num])
    linea.set_3d_properties(z[:num])
    punto.set_data([x[num]], [y[num]])
    punto.set_3d_properties([z[num]])
    return linea, punto


anim = animation.FuncAnimation(fig, actualizar, frames=n_frames, init_func=inicializar, interval=intervalo, blit=True, repeat=True)

plt.show()