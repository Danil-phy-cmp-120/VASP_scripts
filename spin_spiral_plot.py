# -*- coding: utf-8 -*-
import numpy as np
import pyvista as pv
import ase.io.vasp
import ase.spacegroup
from ase.visualize.plot import plot_atoms



def get_magmom():
    f = open('OUTCAR', "r")
    outcar = f.readlines()
    f.close()

    for i in np.arange(len(outcar)-1000, len(outcar)):
        inp = outcar[i].split()
        if len(inp) > 1 and inp[0] == 'magnetization' and inp[1] == '(x)':
            n = i
   
    magmom_x = []; magmom_y = []; magmom_z = []
    for line in outcar[n+3:n+len(cell_0.positions)+4]:
        inp = line.split()
        if len(inp) >1:
            magmom_x += [float(inp[4])]   
    for line in outcar[n+len(cell_0.positions)+13:n+2*len(cell_0.positions)+13]:
        inp = line.split()
        if len(inp) >1:
            magmom_y += [float(inp[4])]  
    for line in outcar[n+2*len(cell_0.positions)+22:n+3*len(cell_0.positions)+22]:
        inp = line.split()
        if len(inp) >1:
            magmom_z += [float(inp[4])]

    return np.column_stack((magmom_x, magmom_y, magmom_z))
   
   
cell_0 = ase.io.vasp.read_vasp("POSCAR")

magmom = get_magmom()


q = 0.25
T = 1
a = q*2*np.pi


cell = cell_0*(1,1,T)

rotation_matrix = np.array([[np.cos(a), -np.sin(a), 0],
                           [np.sin(a), np.cos(a), 0],
                           [0, 0, 1]])

print(a*180/np.pi, np.dot(np.array([4,0,0]), rotation_matrix))


n = 0
for i in range(T):
    for i in range(len(cell_0)): 
        if cell[i+n].symbol == 'Mn':
            cell[i+n].magmom = magmom[i,:]
        if cell[i+n].symbol == 'Au':
            cell[i+n].magmom = np.array([0,0,0])
    n += len(cell_0)
    #magmom = np.dot(magmom, rotation_matrix)
    

colors = {'Mn':'#f24c00', 'Au':'#e7e7e7'}
opacities = {'Mn':0.6, 'Au':0.1}
p = pv.Plotter()

### axes origin abc ###
ax = np.dot(np.array([1, 0, 0]), cell.cell)
ay = np.dot(np.array([0, 1, 0]), cell.cell)
az = np.dot(np.array([0, 0, 1]), cell.cell)
p.add_arrows(np.array([-1, -1, -2]), ax/np.linalg.norm(ax), mag=4, color='r')
p.add_arrows(np.array([-1, -1, -2]), ay/np.linalg.norm(ay), mag=4, color='g')
p.add_arrows(np.array([-1, -1, -2]), az/np.linalg.norm(az), mag=4, color='b')

points = np.array([[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0], [0, 0, 0], [0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0], [0, 1, 0]])
points = np.dot(points, cell.cell)  
p.add_lines(points, color='grey', width=3)
for i in range(1, T):
    points = np.array([[0, 0, i*1/T], [1, 0, i*1/T], [1, 1, i*1/T], [0, 1, i*1/T], [0, 0, i*1/T]])
    points = np.dot(points, cell.cell)  
    p.add_lines(points, color='grey', width=2)


for c in cell:   
    for i in np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0], [1,1,0]]):  
        new_position = c.position + np.dot(i, cell.cell)
        if new_position[0] <= np.dot(np.array([1,0,0]), cell.cell)[0] and new_position[1] <= np.dot(np.array([0,1,0]), cell.cell)[1] and new_position[2] <= np.dot(np.array([0,0,1]), cell.cell)[2]:
            mesh = pv.Sphere(radius=0.7, center=(new_position))
            p.add_mesh(mesh, smooth_shading=True, opacity=opacities[c.symbol], color = colors[c.symbol])
            p.add_arrows(c.position, c.magmom, mag=(0.15*(c.magmom[0]**2+c.magmom[1]**2+c.magmom[2]**2)**0.5), color = 'green', scalar_bar_args={'title': r'', 'color': 'k', 'label_font_size': 26, 'title_font_size': 26})

   
mesh = pv.Sphere(radius=0.7, center=(np.dot(np.array([0, 0, 1]), cell.cell)))
p.add_mesh(mesh, smooth_shading=True, opacity=opacities[cell[0].symbol], color = colors[cell[0].symbol])
mesh = pv.Sphere(radius=0.7, center=(np.dot(np.array([1, 0, 1]), cell.cell)))
p.add_mesh(mesh, smooth_shading=True, opacity=opacities[cell[0].symbol], color = colors[cell[0].symbol])
mesh = pv.Sphere(radius=0.7, center=(np.dot(np.array([0, 1, 1]), cell.cell)))
p.add_mesh(mesh, smooth_shading=True, opacity=opacities[cell[0].symbol], color = colors[cell[0].symbol])
mesh = pv.Sphere(radius=0.7, center=(np.dot(np.array([1, 1, 1]), cell.cell)))
p.add_mesh(mesh, smooth_shading=True, opacity=opacities[cell[0].symbol], color = colors[cell[0].symbol])
p.add_arrows(np.dot(np.array([0, 0, 1]), cell.cell), cell[0].magmom, mag=(0.15*(cell[0].magmom[0]**2+cell[0].magmom[1]**2+cell[0].magmom[2]**2)**0.5), color = 'green', scalar_bar_args={'title': r'', 'color': 'k', 'label_font_size': 26, 'title_font_size': 26})

p.set_position([-14, 8, 14], reset=True)
p.camera.zoom(1.8)

p.set_background('w')
p.show(screenshot = True, full_screen=True)
p.screenshot(filename='swd.png', transparent_background=True)
