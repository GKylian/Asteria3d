import pandas
import numpy as np
import matplotlib.pyplot as plt
import sys



name = 'dens'
nid = 21



N = pandas.read_csv("../exports/"+name+"_0.csv", header=None, skiprows=None, nrows=1).to_numpy()
N = N[0]
Nx = N[0]
Ny = N[1]
Nz = N[2]

dim = (int)(Nx>1) + (int)(Ny>1) + (int)(Nz>1)
print("Dimensions:",dim)

if dim > 1:
    sys.exit("There is more than one dimension")

d = pandas.read_csv("../exports/"+name+"_0.csv", header=None, skiprows=1, nrows=1).to_numpy()
d = d[0]
dx = d[0]
dy = d[1]
dz = d[2]

print(f"Resolution: {Nx} x {Ny} x {Nz}")
print(f"Step sizes: {dx} & {dy} & {dz}")

p0 = pandas.read_csv("../exports/"+name+"_0.csv", header=None, skiprows=2, nrows=1).to_numpy()
p0 = p0[0]
x0 = p0[0]
y0 = p0[1]
z0 = p0[2]

xn = x0+(Nx-1)*dx
yn = y0+(Ny-1)*dy
zn = z0+(Nz-1)*dz

print(f"Domain range:[{x0}; {xn}] x [{y0}; {yn}] x [{z0}; {zn}]")


def plot(i, var):
    path = "../exports/"+var+"_"+str(i)+".csv"
    data = pandas.read_csv(path, header=None, skiprows=3, nrows=1).to_numpy()
    data = data[0]
    print("Rendering data at ",path)

    print("Minimum:",np.min(data))
    print("Maximum:",np.max(data))
    
    x = 0
    if Nx > 1:
        x = np.linspace(x0, xn, Nx)
    if Ny > 1:
        x = np.linspace(y0, yn, Ny)
    if Nz > 1:
        x = np.linspace(z0, zn, Nz)

    #plt.xlim(0.1,0.3)
    #plt.ylim(0.25,0.45)
    plt.plot(x, data)
    plt.show()
    #plt.savefig(var+".png", dpi=800)

plot(nid, name)
