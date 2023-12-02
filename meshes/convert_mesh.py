import sys
import numpy as np
import meshio

filename = sys.argv[1]
mesh = meshio.read(filename)
mesh.write('mesh.mesh')

print(mesh)

with open('mesh.mesh', 'r') as file :
    data = file.read()
    
    print("\nSPHERE: 0 - OUTLET, 1 - WALL, 2 - SYMMETRYY, 3 - SYMMETRYZ\n")
    
    print('SPHERE')
    print('-5 ', data.count(' -5\n'))
    print('-6 ', data.count(' -6\n'), '-> 0, OUTLET')
    print('-7 ', data.count(' -7\n'), '-> 1, WALL')
    print('-8 ', data.count(' -8\n'), '-> 2, SYMMETRYY')
    print('-9 ', data.count(' -9\n'), '-> 3, SYMMETRYZ')
    print('-10', data.count(' -10\n'))

    data = data.replace(' -6\n', ' 0\n')
    data = data.replace(' -7\n', ' 1\n')
    data = data.replace(' -8\n', ' 2\n')
    data = data.replace(' -9\n', ' 3\n')

with open('mesh.mesh', 'w') as file:
    file.write(data)
