import numpy as np
import matplotlib.pyplot as plt


height = 200
width = 500

cylinder_diameter = 50
cylinder_radius = cylinder_diameter / 2
cylinder_center = [(height / 2), 100]


def in_circle(x, y): 
    """ Determine if parameters x and y are in the cylinder's 2d circle profile. """
    dist = np.sqrt((x - cylinder_center[0])**2 + (y - cylinder_center[1])**2)
    if (dist <= cylinder_radius):
        return True
    return False

solid_body2 = []
rows = []
cols = []
solid_body3 = np.zeros((width, height))

def solid_body_setup(width, height):
    for i in range(width):
        for j in range(height): 
            if in_circle(i, j):
                # solid_body3[i, j] = 1
                rows.append(i)
                cols.append(j)
            # dist = np.sqrt((i - cylinder_center[0])**2 + (j - cylinder_center[1])**2)
            # if (dist <= cylinder_radius):
            #     solid_body2.append((i, j))

solid_body_setup(width, height)

solid_body3[rows, cols] = 1

figNum = 1
fig = plt.figure(figNum)
plt.axes().set_aspect('equal')
data_graphable = np.flipud(np.rot90(solid_body3))

plt.pcolor(data_graphable)

plt.show()
