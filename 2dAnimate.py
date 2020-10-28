import pandas
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import cv2
import os
import os.path

name = 'dens'



N = pandas.read_csv("exports/"+name+"_0.csv", header=None, skiprows=None, nrows=1).to_numpy()
N = N[0]
Nx = N[0]
Ny = N[1]
Nz = N[2]

dim = (int)(Nx>1) + (int)(Ny>1) + (int)(Nz>1)
print("Dimensions:",dim)

if dim != 2:
    sys.exit("There is not just one dimension")

d = pandas.read_csv("exports/"+name+"_0.csv", header=None, skiprows=1, nrows=1).to_numpy()
d = d[0]
dx = d[0]
dy = d[1]
dz = d[2]

print(f"Resolution: {Nx} x {Ny} x {Nz}")
print(f"Step sizes: {dx} & {dy} & {dz}")

p0 = pandas.read_csv("exports/"+name+"_0.csv", header=None, skiprows=2, nrows=1).to_numpy()
p0 = p0[0]
x0 = p0[0]
y0 = p0[1]
z0 = p0[2]

xn = x0+(Nx-1)*dx
yn = y0+(Ny-1)*dy
zn = z0+(Nz-1)*dz

print(f"Domain range:[{x0}; {xn}] x [{y0}; {yn}] x [{z0}; {zn}]")

files = []

print('\nMaking all the plots - this may take a while')

#fig, ax = plt.subplots(figsize=(5, 5))

def plot(i, var):
    path = "exports/"+var+"_"+str(i)+".csv"
    if not os.path.isfile(path):
        return False
    data = pandas.read_csv(path, header=None, skiprows=3, nrows=1).to_numpy()
    data = data[0]
    print("Rendering data at ",path)

    print("Minimum:",np.min(data))
    print("Maximum:",np.max(data))
    
    x1 = 0
    x2 = 0
    if Nx > 1 and Ny > 1:
        x1 = np.linspace(x0, xn, Nx)
        x2 = np.linspace(y0, yn, Ny)
        data = np.reshape(data, (-1,Nx))
    if Nx > 1 and Nz > 1:
        x1 = np.linspace(x0, xn, Nx)
        x2 = np.linspace(z0, zn, Nz)
        data = np.reshape(data, (-1,Nx))
    if Ny > 1 and Nz > 1:
        x1 = np.linspace(y0, yn, Ny)
        x2 = np.linspace(z0, zn, Nz)
        data = np.reshape(data, (-1,Ny))

    fig, ax = plt.subplots()
    #cont = ax.contourf(x1, x2, data, 300)
    cont = ax.pcolormesh(x1, x2, data)

    if Nx > 1 and Ny > 1:
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    if Nx > 1 and Nz > 1:
        ax.set_xlabel('x')
        ax.set_ylabel('z')
    if Ny > 1 and Nz > 1:
        ax.set_xlabel('y')
        ax.set_ylabel('z')
    
    ax.set_title(path)
    ax.set_aspect(1)

    cbar = fig.colorbar(cont)
##    plt.contourf(x1, x2, data, 300)
##    fig.colorbar()
##    plt.Axes.set_aspect(aspect=2)

    fname = 'exports/_tmp%03d.png' % i
    files.append(fname)
    plt.savefig(fname, dpi=300)

    return True


for i in range(50000):
    if not plot(i, name):
        break

    

print('Making movie 2dmovie.avi - this may take a while')

img_array = []
path = 'exports/*.png'
for filename in glob.glob(path):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)

out = cv2.VideoWriter(name+'.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)

for i in range(len(img_array)):
    out.write(img_array[i])
out.release()


# cleanup
for fname in files:
    os.remove(fname)




