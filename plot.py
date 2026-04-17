import uproot
import numpy as np
import matplotlib.pyplot as plt

# Abrir archivo ROOT
file = uproot.open("avance_puntos1y2.root")
tree = file["TrackTree"]

# Leer ramas
x = tree["x"].array(library="np")
y = tree["y"].array(library="np")
z = tree["z"].array(library="np")
es_hit = tree["es_hit"].array(library="np")

# ----------------- 1) XY -----------------
plt.figure()
plt.plot(x, y, label="Trayectoria", color="pink")
plt.scatter(x[es_hit], y[es_hit], s=60, label="Hits", color="magenta", zorder=2)
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.grid(True)
plt.legend()
plt.show()
