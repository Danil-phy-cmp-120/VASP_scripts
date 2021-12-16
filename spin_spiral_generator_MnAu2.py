# -*- coding: utf-8 -*-
import numpy as np
import pyvista as pv
import ase.io.vasp

theta = 90
q_z = 0.1

cell_0 = ase.io.vasp.read_vasp("POSCAR")
T = int(1/q_z)
cell = cell_0 * (1,1,T)

q = (2*np.pi/cell_0.cell[2,2])*np.array([0, 0, q_z])


for c in cell: 
    if c.symbol == 'Mn':
    
        mx = np.cos(np.dot(q, c.position)) * np.sin((theta*np.pi)/180)
        my = np.sin(np.dot(q, c.position)) * np.sin((theta*np.pi)/180)
        mz = np.cos((theta*np.pi)/180)
        print(mx, my, mz) 
        c.magmom = np.array([mx, my, mz])
    if c.symbol == 'Au':
        c.magmom = np.array([0,0,0])
    

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
    for i in np.array([[0,0,0]]): #, [0,0,1], [0,1,0], [1,0,0], [1,1,0]]):  
        new_position = c.position + np.dot(i, cell.cell)
        if new_position[0] <= np.dot(np.array([1,0,0]), cell.cell)[0] and new_position[1] <= np.dot(np.array([0,1,0]), cell.cell)[1] and new_position[2] <= np.dot(np.array([0,0,1]), cell.cell)[2]:
            mesh = pv.Sphere(radius=0.7, center=(new_position))
            p.add_mesh(mesh, smooth_shading=True, opacity=opacities[c.symbol], color = colors[c.symbol])
            p.add_arrows(c.position, c.magmom, mag=(2*(c.magmom[0]**2+c.magmom[1]**2+c.magmom[2]**2)**0.5), scalar_bar_args={'title': r'', 'color': 'k', 'label_font_size': 26, 'title_font_size': 26})

p.set_position([-14, 7, 7], reset=True)
p.camera.zoom(1.8)

p.set_background('w')
p.show(screenshot = True, full_screen=True)
p.screenshot(filename='swd.png', transparent_background=True)
