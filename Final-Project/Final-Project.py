import numpy as np
import matplotlib.pyplot as plt

fileName = "Final-Project"

# Grid points:
height = 200
width = 500

num_time_steps = 10

cylinder_diameter = 50
cylinder_radius = cylinder_diameter / 2
cylinder_center = [(height / 2), 100]

error_limit = 0.01      # 1% maximum change for convergence

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
h = 28                  # W/(m^2*C)  Convective Heat Transfer Coefficient
k = 5.2                 # W/(m*C)    Thermal Conductivity
alpha = k / (rho * c)   # m^2/s      Thermal Diffusivity
nu = 1.48 * 10**(-5)    # m^s/s      kinematic viscosity

h_1 = 9 * nu / U_inf #10 * nu / U_inf
h_2 = 9 * alpha / U_inf #10 * alpha / U_inf
h = min(h_1, h_2)       # grid spacing

dt = (h / U_inf) / 2


# Create array and initialize to T-initial
omega = np.zeros((width, height))           # vorticity
psi = np.zeros((width, height))             # streamfunction
temp = np.zeros((width, height)) + T_init   # temperature


def in_circle(x, y): 
    """ Determine if parameters x and y are in the cylinder's 2d circle profile. """
    dist = np.sqrt((x - cylinder_center[0])**2 + (y - cylinder_center[1])**2)
    if (dist <= cylinder_radius):
        return True
    return False

def is_boundary_condition(x, y):
    """ @return: True if (x,y) is a position determined by boundary conditions. """
    if (y == 0):
        return True
    if (y == (height - 1)):
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

def set_initial_conditions():
    for i in range(width):
        for j in range(height):
            if in_circle(i, j):
                psi[i, j] = 0
                
            # establish initial uniform gradient in psi
            if not is_boundary_condition(i, j):
                psi[i, j] = U_inf * j - free_lid

            if solid_boundary(i, j):
                temp[i, j] = T_surface

        psi[i, cylinder_center[1]] = 0 
        psi[i, 0] = -free_lid
        psi[i, (height - 1)] = free_lid

        omega[i, 0] = T_boundary
        omega[i, (height - 1)] = T_boundary
set_initial_conditions()

def gauss_seidel_iteration(data, datatype, omega = ""):
    """ 
    Perform Gauss-Seidel Iteration 

    @param data: 2D array (width, height) of values to be relaxed by Gauss-Seidel Iteration
    @param datatype: either "psi" or "omega". Determines which Poisson/Laplacian equation will be used.
    @param omega (necesary only for datatype "psi"): omega array used to calculate new psi.

    @return: Data array post-relaxation iteration (same dimensions as @param data). 
    """
    error_flag = True
    while error_flag:
        large_error_term_found = False

        # Gauss-Seidel Iteration
        for i in range(width):
            for j in range(height): 
                datum_old = data[i, j]

                if not is_boundary_condition(i, j):
                    datum_old = data[i, j]
                    if datatype == "psi":
                        data[i, j] = data[i, j] + (F / 4) * (data[i + 1, j] + data[i - 1, j] + data[i, j + 1] + data[i, j - 1] - 4 * data[i, j])
                    elif datatype == "omega":
                        data[i, j] = data[i, j] + (F / 4) * (data[i + 1, j] + data[i - 1, j] + data[i, j + 1] + data[i, j - 1] + 4 * (h ** 2) * omega - 4 * psi[i, j]) 

                    if not large_error_term_found:
                        error_term = abs((data[i, j] - datum_old) / datum_old)
                        if (error_term <= error_limit):
                            error_flag = False
                        else:
                            error_flag = True
                            large_error_term_found = True
    return data

psi = gauss_seidel_iteration(psi, "psi")



###############################################################
#  Initial Conditions established, now "turn on" vorticity
###############################################################
u = np.zeros((width, height))
v = np.zeros((width, height))

omega_history = [omega.copy()]
psi_history = [psi.copy()]
temp_history = [temp.copy()]


for n in range(1, num_time_steps):
    for i in range(width):
        for j in range(height):
            if not is_boundary_condition(i, j):
                # Bulk Fluid
                u[i, j] = (psi[i, j + 1] - psi[i, j - 1]) / (2 * h)
                v[i, j] = (psi[i - 1, j] - psi[i + 1, j]) / (2 * h)

                delta_u_omega = 0
                if (u[i, j] < 0):
                    delta_u_omega = u[i + 1, j] * omega[i + 1, j] - u[i, j] * omega[i, j]
                elif (u[i, j] > 0):
                    delta_u_omega = u[i, j] * omega[i, j] - u[i - 1, j]
                    
                delta_v_omega = 0
                if (v[i, j] < 0):
                    delta_v_omega = v[i, j + 1] * omega[i, j + 1] - v[i, j] * omega[i, j]
                elif (v[i, j] > 0):
                    delta_v_omega = v[i, j] * omega[i, j] - v[i, j - 1] * omega[i, j - 1]

                vorticity_laplacian = (omega[i + 1, j] + omega[i - 1, j] + omega[i, j + 1] + omega[i, j - 1] - 4 * omega[i, j]) / (h ** 2)

                omega[i, j] = omega[i, j] + dt * (- delta_u_omega / h - delta_v_omega / h + v * vorticity_laplacian)



                u_delta_T = 0
                if (u[i, j] < 0):
                    u_delta_T = u[i, j] * (temp[i + 1, j] - temp[i, j])
                elif (u[i, j] < 0): 
                    u_delta_T = u[i, j] * (temp[i, j] - temp[i - 1, j])

                v_delta_T = 0
                if (v[i, j] < 0):
                    v_delta_T = v[i, j] * (temp[i, j + 1] - temp[i, j])
                elif (v[i, j] > 0):
                    v_delta_T = v[i, j] * (temp[i, j] - temp[i, j - 1])
                
                temp_laplacian = (temp[i + 1, j] + temp[i - 1, j] + temp[i, j + 1] + temp[i, j - 1] - 4 * temp[i, j]) / (h ** 2)

                temp[i, j] = temp[i, j] + dt * (- u_delta_T / h - v_delta_T / h + alpha * temp_laplacian)

            if solid_boundary(i, j, checked = False):
                (adj_ext_pt_x, adj_ext_pt_y) = solid_boundary(i, j, checked = True)
                omega[i, j] = -2 * (psi[adj_ext_pt_x, adj_ext_pt_y] - psi[i, j]) / (h ** 2)

            if is_outflow_boundary(i, j):
                psi[i, j] = 2 * psi[i - 1, j] - psi[i - 2, j]
                omega[i, j] = omega[i - 1, j]
            

    omega_history.append(omega.copy())
    psi_history.append(psi.copy())
    temp_history.append(temp.copy())



###############################################################
#  Graphs and Plots
###############################################################

figNum = 0

def print_data_in_console():
    """ Print the data in the console (readable format) """
    print(np.rot90(psi))
# print_data_in_console()


def plot_streamlines():
    data_graphable = np.flipud(np.rot90(psi))

    figNum += 1
    fig = plt.figure(figNum)
    plt.axes().set_aspect('equal')


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

    plt.savefig(fileName + "/images/" + fileName + "-Figure-" + str(figNum) + ".png")
    plt.show()
# plot_streamlines()


def plot_temperatures():
    data_graphable = np.flipud(np.rot90(omega))

    figNum += 1
    plt.figure(figNum)
    plt.axes().set_aspect('equal')
    plt.style.use('classic')

    heatmap = plt.pcolor(data_graphable)

    plt.axis("off")

    plt.xlim(0, width)
    plt.ylim(0, height)

    cbar = plt.colorbar(heatmap)
    cbar.set_label("Temperature (\N{DEGREE SIGN}C)")
    plt.clim(np.amin(data_graphable), np.amax(data_graphable))

    plt.savefig(fileName + "/images/" + fileName + "-Figure-" + str(figNum) + ".png")
    plt.show()
plot_temperatures()