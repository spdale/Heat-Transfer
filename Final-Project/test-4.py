import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

import cv2
import os

import time
start_time = time.time()


fileName = "Final-Project"

height = 200
width = 500

num_time_steps = 10

cylinder_diameter = 50
cylinder_radius = cylinder_diameter / 2
cylinder_center = [(height / 2), 100]

error_limit = 0.01                  # 1% maximum change for convergence

U_inf = 2                           # m/s uniform inflow
F = 1.9                             # over-relaxation factor
free_lid = U_inf * (height / 2)     # free-lid streamfunction constant

Re_D = 200                          # Given Reynolds number

T_surface = 400                     # K
T_boundary = 300                    # K
T_init = min(T_surface, T_boundary) # Bulk fluid initial temp

# Constants picked for air around room temp
rho = 3000              # kg/m^3
c = 840                 # J/(kg*C)
#h = 28                  # W/(m^2*C)  Convective Heat Transfer Coefficient
k = 5.2                 # W/(m*C)    Thermal Conductivity
alpha = k / (rho * c)   # m^2/s      Thermal Diffusivity
nu = 1.48 * 10**(-5)    # m^s/s      kinematic viscosity

h_1 = 9 * nu / U_inf #10 * nu / U_inf
h_2 = 9 * alpha / U_inf #10 * alpha / U_inf
h = min(h_1, h_2)       # grid spacing

dt = (h / U_inf) / 2


omega = np.zeros((width, height))           # vorticity
psi = np.zeros((width, height))             # streamfunction
temp = np.zeros((width, height)) + T_init   # temperature

solid_body = np.zeros((width, height), dtype=bool)

def circle_setup(width, height):
    for i in range(width):
        for j in range(height): 
            dist = np.sqrt((i - cylinder_center[0])**2 + (j - cylinder_center[1])**2)
            if (dist <= cylinder_radius):
                solid_body[i][j] = True

circle_setup(width, height)


# solid_body2 = []

# def solid_body_setup(width, height):
#     for i in range(width):
#         for j in range(height): 
#             dist = np.sqrt((i - cylinder_center[0])**2 + (j - cylinder_center[1])**2)
#             if (dist <= cylinder_radius):
#                 solid_body2.append([i, j])

# solid_body_setup(width, height)





def in_circle(x, y): 
    """ Determine if parameters x and y are in the cylinder's 2d circle profile. """
    dist = np.sqrt((x - cylinder_center[0])**2 + (y - cylinder_center[1])**2)
    if (dist <= cylinder_radius):
        return True
    return False




solid_rows = []
solid_cols = []

def solid_body_setup(width, height):
    for i in range(width):
        for j in range(height): 
            dist = np.sqrt((i - cylinder_center[0])**2 + (j - cylinder_center[1])**2)
            if (dist <= cylinder_radius):
                solid_rows.append(i)
                solid_cols.append(j)

solid_body_setup(width, height)


# print(solid_body2)

solid_body3 = np.zeros((width, height))
solid_body3[solid_rows][solid_cols] = 1


figNum = 1
fig = plt.figure(figNum)
plt.axes().set_aspect('equal')
data_graphable = np.flipud(np.rot90(solid_body3))

plt.pcolor(data_graphable)

plt.show()


def is_boundary_condition(x, y):
    """ @return: True if (x,y) is a position determined by boundary conditions. """
    if (y == 0):
        return True
    if (y == (height - 1)):
        return True
    if (x == 0):
        return True
    if (x == (width - 1)):
        return True
    if in_circle(x, y):
        return True
    return False

def solid_boundary(x, y, checked = False):
    """ 
    Determine if (x,y) is a point along the wall of a solid body.

    @param checked: False if checking whether the point is a solid body.
                    True if determining what the adjacent free (outside of solid body) point is. 
    
    @return: if checked = False: True if a wall point (adjacent to fluid) of body.
             if checked = True: returns (x, y) values of the adjacent point in the fluid (i.e. outside of the solid body).
    """
    if not checked: 
        if in_circle(x, y) and (not in_circle(x - 1, y) or not in_circle(x + 1, y) or not in_circle(x, y - 1) or not in_circle(x, y + 1)): 
            return True
    if checked: 
        if in_circle(x, y):
            if not in_circle(x - 1, y):
                return (x - 1, y)
            if not in_circle(x + 1, y):
                return (x + 1, y)
            if not in_circle(x, y - 1):
                return (x, y - 1)
            if not in_circle(x, y + 1):
                return (x, y + 1)
    return False

def is_outflow_boundary(x, y):
    """ Determine if (x,y) is at outflow boundary. """
    if (x == (width - 1)) and (y != 0) and (y != (height - 1)):
        return True
    return False

def gauss_seidel_iteration(data, error_limit, omega = "", initial = False):
    """ 
    Perform Gauss-Seidel Iteration 

    @param data: 2D array (width, height) of values to be relaxed by Gauss-Seidel Iteration
    @param omega (necesary not necessary for initial): omega array used to calculate new psi.
    @param initial: Determines which Poisson/Laplacian equation will be used.

    @return: Data array post-relaxation iteration (same dimensions as @param data). 
    """
    error_flag = True
    while error_flag:
        data_old = data.copy()

        data[1:-1, 1:-1] = data[1:-1, 1:-1] + (1 / 4) * (data[0:-2, 1:-1] + data[2:, 1:-1] + data[1:-1,0:-2] + data[1:-1, 2:] - 4 * data[1:-1, 1:-1])
        
        if not initial:
            data +=  h * h * omega  # Multiply by F

        data[0, :] = data[3,:]
        data[width - 1, :] = data[width - 2, :]

        
        data[solid_body2] = 0


        temp = np.absolute(data - data_old)
        error_array = np.divide(temp, data_old, out = np.zeros_like(data), where = ((temp != 0) & (data_old != 0)))
        error_array[error_array == np.inf] = 0
        error_term = np.amax(error_array)

        if (error_term <= error_limit):
            error_flag = False





        # large_error_term_found = False

        # for i in range(width):
        #     for j in range(height): 
        #         datum_old = data[i, j]

        #         if not is_boundary_condition(i, j):
        #             datum_old = data[i, j]
        #             if initial:
        #                 data[i, j] = data[i, j] + (F / 4) * (data[i + 1, j] + data[i - 1, j] + data[i, j + 1] + data[i, j - 1] - 4 * data[i, j])
        #             else:
        #                 data[i, j] = data[i, j] + (F / 4) * (data[i + 1, j] + data[i - 1, j] + data[i, j + 1] + data[i, j - 1] + 4 * (h ** 2) * omega[i, j] - 4 * psi[i, j]) 
                    
        #             if not large_error_term_found:
        #                 error_term = abs((data[i, j] - datum_old) / datum_old)
        #                 if (error_term <= error_limit):
        #                     error_flag = False
        #                 else:
        #                     error_flag = True
        #                     large_error_term_found = True
                
        #         if (i == 1):
        #             data[0, j] = data[i + 2, j]
        #         elif (i == (width - 2)):
        #             data[i + 1, j] = data[i, j]
    return data



###############################################################
#  Initial Conditions
###############################################################
for i in range(width):
    for j in range(height):
        if in_circle(i, j):
            psi[i, j] = 0
            temp[i, j] = T_surface
            
        # establish initial uniform gradient in psi
        if not is_boundary_condition(i, j):
            psi[i, j] = U_inf * j - free_lid
        
        # # Boundary Conditions for psi
        # if (i == 1):
        #         psi[0, j] = psi[i + 2, j]
        # elif (i == (width - 2)):
        #     psi[i + 1, j] = psi[i, j]

    psi[i, cylinder_center[1]] = 0 
    psi[i, 0] = -free_lid
    psi[i, (height - 1)] = free_lid

    omega[i, 0] = T_boundary
    omega[i, (height - 1)] = T_boundary

psi = gauss_seidel_iteration(psi, error_limit, initial = True)





print("--- %s seconds ---" % (time.time() - start_time))






figNum = 1
fig = plt.figure(figNum)
plt.axes().set_aspect('equal')
data_graphable = np.flipud(np.rot90(psi))


num_streamlines = 31
max_streamline = np.max(data_graphable)
min_streamline = np.min(data_graphable)
contours_before = np.linspace(min_streamline, max_streamline, num=(num_streamlines + 3))
contours = contours_before[(contours_before != 0) & (contours_before != min_streamline) & (contours_before != max_streamline)]

plt.contour(data_graphable, levels = contours, colors = 'black', linestyles = 'solid')


plt.xlim(0, width)
plt.ylim(0, height)
plt.xticks(np.arange(0, width + 1, 50))
plt.yticks(np.arange(0, height + 1, 20))
plt.tick_params(top=True, right=True)

plt.style.use('grayscale')
heatmap = plt.pcolor(data_graphable)
plt.clim(np.amin(data_graphable), np.amax(data_graphable))

plt.savefig(fileName + "/images/" + fileName + "-Figure.png")
plt.show()
