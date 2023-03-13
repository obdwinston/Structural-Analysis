import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mesh = 'wrench.su2'  # mesh file

E = 200e9  # Young's modulus [Pa]
nu = .25  # Poisson's ratio
t = .01  # thickness [m]

Np = [895, 903, 1089, 1327, 1268, 1569, 345, 428, 1870, 1869, 863, 1865, 1824, 1819, 1818, 1789, 1695, 1672, 1521, 1662, 369, 1553, 1326, 1087, 1306, 1578]  # point force node [Np1, Np2, ...]
Fp = []  # point force [[Fpx, Fpy], ...]

p = 10e6  # surface force
L = .06  # surface length
Fs = p*L/(len(Np) - 1)*t/2  # point force equivalent
for i in range(len(Np)):
    if i == 0 or i == len(Np) - 1:
        Fp.append([0., -Fs])
    else:
        Fp.append([0., -2*Fs])

print(len(Fp))
print(Fp)

# program start

file = open(mesh, 'r')
lines = file.readlines()

cells = []
nodes = []

n_cells = int(lines[1].split()[1])
n_original_nodes = int(lines[2 + n_cells].split()[1])

original_nodes = []  # coordinates
original_nodes_list = []  # indices
for i in range(n_original_nodes):
    nx = float(lines[3 + n_cells + i].split()[0])
    ny = float(lines[3 + n_cells + i].split()[1])
    original_nodes.append([nx, ny])

for i in range(n_cells):
    n1 = int(lines[2 + i].split()[1])
    n2 = int(lines[2 + i].split()[2])
    n3 = int(lines[2 + i].split()[3])
    
    if n1 not in original_nodes_list:
        original_nodes_list.append(n1)
        nodes.append([original_nodes[n1][0], original_nodes[n1][1]])
    if n2 not in original_nodes_list:
        original_nodes_list.append(n2)
        nodes.append([original_nodes[n2][0], original_nodes[n2][1]])
    if n3 not in original_nodes_list:
        original_nodes_list.append(n3)
        nodes.append([original_nodes[n3][0], original_nodes[n3][1]])
n_nodes = len(nodes)

for i in range(n_cells):
    nodes_list = []

    x1 = original_nodes[int(lines[2 + i].split()[1])][0]
    y1 = original_nodes[int(lines[2 + i].split()[1])][1]
    for j in range(n_nodes):
        if x1 == nodes[j][0] and y1 == nodes[j][1]:
            nodes_list.append(j)

    x2 = original_nodes[int(lines[2 + i].split()[2])][0]
    y2 = original_nodes[int(lines[2 + i].split()[2])][1]
    for j in range(n_nodes):
        if x2 == nodes[j][0] and y2 == nodes[j][1]:
            nodes_list.append(j)

    x3 = original_nodes[int(lines[2 + i].split()[3])][0]
    y3 = original_nodes[int(lines[2 + i].split()[3])][1]
    for j in range(n_nodes):
        if x3 == nodes[j][0] and y3 == nodes[j][1]:
            nodes_list.append(j)

    cells.append(nodes_list)

original_fixed_nodes = []
for line in lines[6 + n_cells + n_original_nodes:]:
    n1 = int(line.split()[1])
    n2 = int(line.split()[2])
    if n1 not in original_fixed_nodes:
        original_fixed_nodes.append(n1)
    if n2 not in original_fixed_nodes:
        original_fixed_nodes.append(n2)
n_original_fixed_nodes = len(original_fixed_nodes)

fixed_nodes = []
for i in range(n_original_fixed_nodes):
    x1 = original_nodes[original_fixed_nodes[i]][0]
    y1 = original_nodes[original_fixed_nodes[i]][1]
    for j in range(n_nodes):
        if x1 == nodes[j][0] and y1 == nodes[j][1]:
            fixed_nodes.append(j)
n_fixed_nodes = len(fixed_nodes)

free_nodes = list(range(n_nodes))
for node in fixed_nodes:
    free_nodes.remove(node)
n_free_nodes = len(free_nodes)

A = []
for i in range(n_cells):
    x1 = nodes[cells[i][0]][0]
    y1 = nodes[cells[i][0]][1]
    x2 = nodes[cells[i][1]][0]
    y2 = nodes[cells[i][1]][1]
    x3 = nodes[cells[i][2]][0]
    y3 = nodes[cells[i][2]][1]

    A.append(.5*abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)))

# check mesh
plt.figure(figsize=(12, 6))

for i in range(n_nodes):
    px = nodes[i][0]
    py = nodes[i][1]
    plt.text(px, py, i, c='r', horizontalalignment='center', verticalalignment='center')

for i in range(n_cells):
        x1 = nodes[cells[i][0]][0]
        y1 = nodes[cells[i][0]][1]
        x2 = nodes[cells[i][1]][0]
        y2 = nodes[cells[i][1]][1]
        x3 = nodes[cells[i][2]][0]
        y3 = nodes[cells[i][2]][1]
        plt.plot([x1, x2, x3, x1], [y1, y2, y3, y1], c='k', linewidth=0.5)

        # ccx = (x1 + x2 + x3)/3
        # ccy = (y1 + y2 + y3)/3
        # plt.text(ccx, ccy, i, c='b', horizontalalignment='center', verticalalignment='center')

plt.axis('equal')
plt.grid('on')
plt.show()

# global stiffness matrix
B = np.zeros((3, 6))
C = E/(1 - nu*nu)*np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu)/2]])
K = np.zeros((2*n_nodes, 2*n_nodes))
for i in range(n_cells):
    x1 = nodes[cells[i][0]][0]
    y1 = nodes[cells[i][0]][1]
    x2 = nodes[cells[i][1]][0]
    y2 = nodes[cells[i][1]][1]
    x3 = nodes[cells[i][2]][0]
    y3 = nodes[cells[i][2]][1]

    b1 = 1/(2*A[i])*(y2 - y3)
    b2 = 1/(2*A[i])*(y3 - y1)
    b3 = 1/(2*A[i])*(y1 - y2)
    c1 = 1/(2*A[i])*(x3 - x2)
    c2 = 1/(2*A[i])*(x1 - x3)
    c3 = 1/(2*A[i])*(x2 - x1)

    B[0, 0] = b1
    B[1, 1] = c1
    B[2, 0] = c1
    B[2, 1] = b1

    B[0, 2] = b2
    B[1, 3] = c2
    B[2, 2] = c2
    B[2, 3] = b2

    B[0, 4] = b3
    B[1, 5] = c3
    B[2, 4] = c3
    B[2, 5] = b3

    k = t*A[i]*np.matmul(np.transpose(B), np.matmul(C, B))

    n1 = cells[i][0]
    n2 = cells[i][1]
    n3 = cells[i][2]

    K[2*n1:2*n1 + 2, 2*n1:2*n1 + 2] += k[0:2, 0:2]
    K[2*n2:2*n2 + 2, 2*n1:2*n1 + 2] += k[2:4, 0:2]
    K[2*n3:2*n3 + 2, 2*n1:2*n1 + 2] += k[4:6, 0:2]

    K[2*n1:2*n1 + 2, 2*n2:2*n2 + 2] += k[0:2, 2:4]
    K[2*n2:2*n2 + 2, 2*n2:2*n2 + 2] += k[2:4, 2:4]
    K[2*n3:2*n3 + 2, 2*n2:2*n2 + 2] += k[4:6, 2:4]

    K[2*n1:2*n1 + 2, 2*n3:2*n3 + 2] += k[0:2, 4:6]
    K[2*n2:2*n2 + 2, 2*n3:2*n3 + 2] += k[2:4, 4:6]
    K[2*n3:2*n3 + 2, 2*n3:2*n3 + 2] += k[4:6, 4:6]

# global force vector
F = np.zeros(2*n_nodes)
for i in range(len(Np)):
    F[2*Np[i]] = Fp[i][0]
    F[2*Np[i] + 1] = Fp[i][1]

# solve displacements
delete = []
for node in fixed_nodes:
    delete.append(2*node)
    delete.append(2*node + 1)
F1 = np.delete(F, delete, 0)
K1 = np.delete(K, delete, 0)
K2 = np.delete(K1, delete, 1)
D1 = np.linalg.solve(K2, F1)

D = np.zeros(2*n_nodes)
for i in range(n_free_nodes):
    ni = free_nodes[i]
    D[2*ni] = D1[2*i]
    D[2*ni + 1] = D1[2*i + 1]

u = np.zeros(n_nodes)
v = np.zeros(n_nodes)
for i in range(n_nodes):
    u[i] = D[2*i]
    v[i] = D[2*i + 1]

# solve stresses
Sa = np.zeros(n_cells)  # average stress
Sv = np.zeros(n_cells)  # von Mises stress
S1 = np.zeros(n_cells)  # maximum principal stress
S2 = np.zeros(n_cells)  # minimum principal stress
for i in range(n_cells):
    x1 = nodes[cells[i][0]][0]
    y1 = nodes[cells[i][0]][1]
    x2 = nodes[cells[i][1]][0]
    y2 = nodes[cells[i][1]][1]
    x3 = nodes[cells[i][2]][0]
    y3 = nodes[cells[i][2]][1]

    b1 = 1/(2*A[i])*(y2 - y3)
    b2 = 1/(2*A[i])*(y3 - y1)
    b3 = 1/(2*A[i])*(y1 - y2)
    c1 = 1/(2*A[i])*(x3 - x2)
    c2 = 1/(2*A[i])*(x1 - x3)
    c3 = 1/(2*A[i])*(x2 - x1)

    B[0, 0] = b1
    B[1, 1] = c1
    B[2, 0] = c1
    B[2, 1] = b1

    B[0, 2] = b2
    B[1, 3] = c2
    B[2, 2] = c2
    B[2, 3] = b2

    B[0, 4] = b3
    B[1, 5] = c3
    B[2, 4] = c3
    B[2, 5] = b3

    n1 = cells[i][0]
    n2 = cells[i][1]
    n3 = cells[i][2]

    d1 = D[2*n1]
    d2 = D[2*n1 + 1]
    d3 = D[2*n2]
    d4 = D[2*n2 + 1]
    d5 = D[2*n3]
    d6 = D[2*n3 + 1]

    d = np.array([d1, d2, d3, d4, d5, d6])
    s = np.matmul(C, np.matmul(B, d))
    sx = s[0]
    sy = s[1]
    txy = s[2]

    Sa[i] = .5*(sx + sy)
    Sv[i] = m.sqrt(sx**2 + sy**2 - sx*sy + 3*txy**2)
    S1[i] = .5*(sx + sy) + m.sqrt((.5*(sx - sy))**2 + txy**2)
    S2[i] = .5*(sx + sy) - m.sqrt((.5*(sx - sy))**2 + txy**2)

# show results
print('Maximum x-deflection = %e m' % max(u))
print('Minimum x-deflection = %e m' % min(u))
print('Maximum y-deflection = %e m' % max(v))
print('Minimum y-deflection = %e m' % min(v))

print('Maximum principal stress = %e Pa' % max(S1))
print('Minimum principal stress = %e Pa' % min(S2))
print('Maximum von Mises stress = %e Pa' % max(Sv))
print('Minimum von Mises stress = %e Pa' % min(Sv))

fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 6), height_ratios=[20, 1])

S = Sv  # change stress
cmap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=min(S), vmax=max(S))
ticks = np.linspace(min(S), max(S), 10)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax2, orientation='horizontal', label='Stress [Pa]', ticks=ticks)

for i in range(n_cells):
    x1 = nodes[cells[i][0]][0]
    y1 = nodes[cells[i][0]][1]
    x2 = nodes[cells[i][1]][0]
    y2 = nodes[cells[i][1]][1]
    x3 = nodes[cells[i][2]][0]
    y3 = nodes[cells[i][2]][1]

    n1 = cells[i][0]
    n2 = cells[i][1]
    n3 = cells[i][2]

    x1n = x1 + D[2*n1]
    y1n = y1 + D[2*n1 + 1]
    x2n = x2 + D[2*n2]
    y2n = y2 + D[2*n2 + 1]
    x3n = x3 + D[2*n3]
    y3n = y3 + D[2*n3 + 1]

    color = cmap(norm(S[i]))
    ax1.plot([x1, x2, x3, x1], [y1, y2, y3, y1], c='k', linewidth=0.5, zorder=1)
    ax1.fill([x1n, x2n, x3n, x1n], [y1n, y2n, y3n, y1n], facecolor=color, edgecolor=None, zorder=2)

ax1.axis('equal')
ax1.grid('on')
plt.show()
