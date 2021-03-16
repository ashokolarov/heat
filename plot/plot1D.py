import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

with open('data/time.txt') as tfile:
    l = tfile.readlines()[0].split()
    t = [float(x) for x in l]

with open('data/coordinates.txt') as xfile:
    l = xfile.readlines()[0].split()
    x = [float(x) for x in l]

temp = []

with open('data/temp.txt') as sfile:
    for l in sfile.readlines():
        te = [float(x) for x in l.split()]
        temp.append(te)


x, t = np.meshgrid(x,t)
temp = np.array(temp)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


# Plot the surface.
surf = ax.plot_surface(x, t, temp, cmap=cm.jet,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

ax.set_xlabel('x-axis', fontsize=10)
ax.set_ylabel('Time [s]', fontsize=10)
ax.set_zlabel('Temperature [K]', fontsize=10)

ax.set_title(r"1D Laplace equation, y(0)=10sin(x), $\alpha$=0.1")
plt.show()
