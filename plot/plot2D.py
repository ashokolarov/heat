import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import matplotlib

plt.rcParams['animation.html'] = 'html5'

with open('data/time.txt') as tfile:
    l = tfile.readlines()[0].split()
    t = [float(x) for x in l]

with open('data/coordinates.txt') as xfile:
    l = xfile.readlines()[0].split()
    x = [float(x) for x in l]

x = np.unique(x)
x,y = np.meshgrid(x,x)

temp = []

with open('data/temp.txt') as sfile:
    for l in sfile.readlines():
        te = [float(x) for x in l.split()]
        temp.append(te)

temp = np.array(temp)
size = x.shape[0]


fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

t=np.array(t)
t = t[t<0.6]
N = len(t)
dt = t[1] - t[0]
fps = 1/dt

wframe = None
Z = temp[0].reshape(size,size)
Zmax = np.max(temp)

ax.set_xlabel('x-axis', fontsize=10)
ax.set_ylabel('y-axis', fontsize=10)
ax.set_zlabel('Temperature [K]', fontsize=10)

def update(idx):
    global wframe
    ax.set_title(fr"Laplace 2D, $\alpha$=0.1, t={t[idx]}[s]")
    if wframe:
        ax.collections.remove(wframe)
    ax.set_zlim(0, Zmax)

    # Plot the new wireframe and pause briefly before continuing.
    wframe = ax.plot_surface(x, y, temp[idx].reshape(size,size), cmap="jet",
                              linewidth=0, antialiased=False, vmin=0, vmax=Zmax)
    
    cb = fig.colorbar(wframe, shrink=0.5, aspect=5)
    plt.pause(.001)
    cb.remove()
    
ani = animation.FuncAnimation(fig, update, N, interval=dt)

fn = 'Laplace 2d'
ani.save(fn+'.gif',writer='imagemagick',fps=10)


