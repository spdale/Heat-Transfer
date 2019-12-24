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

num_time_steps = 2

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

def test_setup():
    solid_body3 = np.zeros((width, height))
    solid_body3[solid_rows, solid_cols] = 1


    figNum = 1
    fig = plt.figure(figNum)
    plt.axes().set_aspect('equal')
    data_graphable = np.flipud(np.rot90(solid_body3))

    plt.pcolor(data_graphable)

    plt.show()
# test_setup()


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

        # data[i, j] = data[i, j] + (F / 4) * (data[i + 1, j] + data[i - 1, j] + data[i, j + 1] + data[i, j - 1] - 4 * data[i, j])

        data[1:-1, 1:-1] = data[1:-1, 1:-1] + (1 / 4) * (data[0:-2, 1:-1] + data[2:, 1:-1] + data[1:-1,0:-2] + data[1:-1, 2:] - 4 * data[1:-1, 1:-1])
        if not initial:
            data += h * h * omega   # Multiply by F

        data[0, :] = data[3,:]
        data[width - 1, :] = data[width - 2, :]

        
        data[solid_rows, solid_cols] = 0


        temp = np.absolute(data - data_old)
        error_array = np.divide(temp, data_old, out = np.zeros_like(data), where = ((temp != 0) & (data_old != 0)))
        error_array[error_array == np.inf] = 0
        error_term = np.amax(error_array)

        if (error_term <= error_limit):
            error_flag = False

    return data



###############################################################
#  Initial Conditions
###############################################################
psi[solid_rows, solid_cols] = 0
temp[solid_rows, solid_cols] = T_surface

psi[:, cylinder_center[1]] = 0
psi[:, 0] = -free_lid
psi[:, (height - 1)] = free_lid

omega[:, 0] = T_boundary
omega[:, (height - 1)] = T_boundary

for i in range(width):
    for j in range(height):
        if not is_boundary_condition(i, j):
            psi[i, j] = U_inf * j - free_lid

psi = gauss_seidel_iteration(psi, error_limit, initial = True)


print("--- %s seconds ---" % (time.time() - start_time))



###############################################################
#  Initial Conditions established, now "turn on" vorticity
###############################################################
u = np.zeros((width, height))
v = np.zeros((width, height))

omega_history = [omega.copy()]
psi_history = [psi.copy()]
temp_history = [temp.copy()]




# data[i, j] = data[i, j] + (F / 4) * (data[i + 1, j] + data[i - 1, j] + data[i, j + 1] + data[i, j - 1] - 4 * data[i, j])
# data[1:-1, 1:-1] = data[1:-1, 1:-1] + (1 / 4) * (data[0:-2, 1:-1] + data[2:, 1:-1] + data[1:-1, 0:-2] + data[1:-1, 2:] - 4 * data[1:-1, 1:-1])
# data[1:-1, 1:-1]  current
# data[0:-2, 1:-1]  left
# data[2:,   1:-1]  right
# data[1:-1, 0:-2]  down
# data[1:-1, 2:  ]    up

delta_u_omega = np.zeros((width, height))
delta_v_omega = np.zeros((width, height))

for n in range(1, num_time_steps):
    u[1:-1, 1:-1] = (psi[1:-1, 2:  ] - psi[1:-1,0:-2]) / (2 * h)
    v[1:-1, 1:-1] = (psi[0:-2, 1:-1] - psi[2:,  1:-1]) / (2 * h)

    delta_u_omega = 0
    delta_v_omega = 0

    delta_u_omega[1:-1, 1:-1] = u[2:,   1:-1] * omega[2:,   1:-1] - u[1:-1, 1:-1] * omega[1:-1, 1:-1]
    delta_u_omega[1:-1, 1:-1] = u[1:-1, 1:-1] * omega[1:-1, 1:-1] - u[] * omega[]

    for i in range(width):
        for j in range(height):
            if not is_boundary_condition(i, j):
                # Bulk Fluid
                delta_u_omega = 0
                if (u[i, j] < 0):
                    delta_u_omega = u[i + 1, j] * omega[i + 1, j] - u[i, j] * omega[i, j]
                elif (u[i, j] > 0):
                    delta_u_omega = u[i, j] * omega[i, j] - u[i - 1, j] * omega[i - 1, j]
                    
                delta_v_omega = 0
                if (v[i, j] < 0):
                    delta_v_omega = v[i, j + 1] * omega[i, j + 1] - v[i, j] * omega[i, j]
                elif (v[i, j] > 0):
                    delta_v_omega = v[i, j] * omega[i, j] - v[i, j - 1] * omega[i, j - 1]

                vorticity_laplacian = (omega[i + 1, j] + omega[i - 1, j] + omega[i, j + 1] + omega[i, j - 1] - 4 * omega[i, j]) / (h ** 2)

                omega[i, j] = omega[i, j] + dt * (-delta_u_omega / h - delta_v_omega / h + nu * vorticity_laplacian)


                psi = gauss_seidel_iteration(psi, error_limit, omega)


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

# for n in range(1, num_time_steps):
#     for i in range(width):
#         for j in range(height):
#             if not is_boundary_condition(i, j):
#                 # Bulk Fluid
#                 u[i, j] = (psi[i, j + 1] - psi[i, j - 1]) / (2 * h)
#                 v[i, j] = (psi[i - 1, j] - psi[i + 1, j]) / (2 * h)

#                 delta_u_omega = 0
#                 if (u[i, j] < 0):
#                     delta_u_omega = u[i + 1, j] * omega[i + 1, j] - u[i, j] * omega[i, j]
#                 elif (u[i, j] > 0):
#                     delta_u_omega = u[i, j] * omega[i, j] - u[i - 1, j] * omega[i - 1, j]
                    
#                 delta_v_omega = 0
#                 if (v[i, j] < 0):
#                     delta_v_omega = v[i, j + 1] * omega[i, j + 1] - v[i, j] * omega[i, j]
#                 elif (v[i, j] > 0):
#                     delta_v_omega = v[i, j] * omega[i, j] - v[i, j - 1] * omega[i, j - 1]

#                 vorticity_laplacian = (omega[i + 1, j] + omega[i - 1, j] + omega[i, j + 1] + omega[i, j - 1] - 4 * omega[i, j]) / (h ** 2)

#                 omega[i, j] = omega[i, j] + dt * (-delta_u_omega / h - delta_v_omega / h + nu * vorticity_laplacian)


#                 psi = gauss_seidel_iteration(psi, error_limit, omega)


#                 u_delta_T = 0
#                 if (u[i, j] < 0):
#                     u_delta_T = u[i, j] * (temp[i + 1, j] - temp[i, j])
#                 elif (u[i, j] < 0): 
#                     u_delta_T = u[i, j] * (temp[i, j] - temp[i - 1, j])

#                 v_delta_T = 0
#                 if (v[i, j] < 0):
#                     v_delta_T = v[i, j] * (temp[i, j + 1] - temp[i, j])
#                 elif (v[i, j] > 0):
#                     v_delta_T = v[i, j] * (temp[i, j] - temp[i, j - 1])
                
#                 temp_laplacian = (temp[i + 1, j] + temp[i - 1, j] + temp[i, j + 1] + temp[i, j - 1] - 4 * temp[i, j]) / (h ** 2)

#                 temp[i, j] = temp[i, j] + dt * (- u_delta_T / h - v_delta_T / h + alpha * temp_laplacian)

#             if solid_boundary(i, j, checked = False):
#                 (adj_ext_pt_x, adj_ext_pt_y) = solid_boundary(i, j, checked = True)
#                 omega[i, j] = -2 * (psi[adj_ext_pt_x, adj_ext_pt_y] - psi[i, j]) / (h ** 2)

#             if is_outflow_boundary(i, j):
#                 psi[i, j] = 2 * psi[i - 1, j] - psi[i - 2, j]
#                 omega[i, j] = omega[i - 1, j]
            

#     omega_history.append(omega.copy())
#     psi_history.append(psi.copy())
#     temp_history.append(temp.copy())


def print_data_in_console():
    """ Print the data in the console (readable format) """
    print(np.rot90(psi))
# print_data_in_console()


###############################################################
#  Graphs and Plots
###############################################################
for plot_index in range(num_time_steps):
    figNum = plot_index
    fig = plt.figure(figNum, figsize=(10, 10.5))



    ###############################################################
    #  Streamfunction Plot
    ###############################################################
    sub1 = plt.subplot(3, 1, 1, aspect = 'equal')
    data_graphable = np.flipud(np.rot90(psi_history[plot_index]))

    plt.title("Streamfunction")

    num_streamlines = 31
    max_streamline = np.max(data_graphable)
    min_streamline = np.min(data_graphable)
    contours_before = np.linspace(min_streamline, max_streamline, num=(num_streamlines + 3))
    contours = contours_before[(contours_before != 0) & (contours_before != min_streamline) & (contours_before != max_streamline)]

    plt.contour(data_graphable, levels = contours, colors = 'black', linestyles = 'solid')


    plt.style.use('grayscale')
    plt.xticks(np.arange(0, width + 1, 50))
    plt.yticks(np.arange(0, height + 1, 20))
    plt.tick_params(top=True, right=True)

    plt.pcolor(data_graphable)

    # Color bar: 
    divider1 = make_axes_locatable(sub1)
    cax1 = divider1.append_axes('right', size = '3%', pad = 0.3)
    im = sub1.imshow(data_graphable, origin = 'lower', aspect = 'equal', interpolation = 'none')
    fig.colorbar(im, cax = cax1, orientation = 'vertical')




    ###############################################################
    #  Vorticity Plot
    ###############################################################
    sub2 = plt.subplot(3, 1, 2, aspect = 'equal')
    data_graphable = np.flipud(np.rot90(omega_history[plot_index]))

    plt.title("Vorticity")
    plt.axis("off")

    plt.pcolor(data_graphable)

    # Color bar: 
    divider2 = make_axes_locatable(sub2)
    cax2 = divider2.append_axes('right', size = '3%', pad = 0.3)
    im = sub2.imshow(data_graphable, origin = 'lower', aspect = 'equal', interpolation = 'none')
    fig.colorbar(im, cax = cax2, orientation = 'vertical')




    ###############################################################
    #  Temperature Plot
    ###############################################################
    sub3 = fig.add_subplot(3, 1, 3, aspect = 'equal')
    data_graphable = np.flipud(np.rot90(temp_history[plot_index]))
    plt.title("Temperature")

    plt.style.use('classic')
    plt.axis("off")

    # Color bar: 
    divider3 = make_axes_locatable(sub3)
    cax3 = divider3.append_axes('right', size = '3%', pad = 0.3)
    im = sub3.imshow(data_graphable, origin = 'lower', aspect = 'equal', interpolation = 'none')
    cbar = fig.colorbar(im, cax=cax3, orientation = 'vertical')
    cbar.set_label("Temperature (\N{DEGREE SIGN}C)")




    plt.savefig(fileName + "/images/" + fileName + "-Figure-" + str(figNum) + ".png")
    plt.show()




###############################################################
#  Generate Video
###############################################################
def generate_video():
    image_folder = "C:\\Users\\samda\\Documents\\GitHub\\Heat-Transfer\\Final-Project\\images"
    output_folder = "C:\\Users\\samda\\Documents\\GitHub\\Heat-Transfer\\Final-Project\\"
    video_name = output_folder + "final-project.mp4"
    fps = 5

    images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))

    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))

    cv2.destroyAllWindows()
    video.release()
# generate_video()