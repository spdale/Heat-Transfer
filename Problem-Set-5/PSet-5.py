import numpy as np
import matplotlib.pyplot as plt

fileName = "Problem-Set-5"

# Grid points:
height = 30
width = 30

# Constants
rho = 3000  # kg/m^3
c = 840     # J/(kg*C)
h = 28      # W/(m^2*C)  Convective Heat Transfer Coefficient
k = 5.2     # W/(m*C)    Thermal Conductivity

# Make set temperatures on fixed positions
T_alpha   = 0    # (bottom boundary temperature)
T_bravo   = 40   # (left boundary temperature)
T_charlie = 100  # (top boundary temperature)
T_delta   = 100  # (right boundary temperature)

T_initial = 15
T_right = 38


# Initialize matrix of zeros for that size
# Note: index 0,0 is bottom left
#default_temp = (max(T_alpha, T_bravo, T_charlie, T_delta) + min(T_alpha, T_bravo, T_charlie, T_delta)) / 2
#data = np.zeros((height, width)) + T_initial

data = np.zeros((width, height)) + T_initial

for j in range(height):
    data[(width - 1), j] = T_right


# data = [[0, 1, 2, 3, 4],
#         [5, 6, 7, 8, 9]]

# # for j in range(height):
# #     for i in range(width):
# #         print(data[i,j], end = ' ')
# #     print("")

# # print(np.flip(data))
# print(data)
# print(np.rot90(data))

# data = np.rot90(data)

# Set boundary conditions
# for i in range(width):
#     data[i, 0] = T_alpha
#     data[i, (height - 1)] = T_charlie
# for j in range(1, (height - 1)):
#     data[0, j] = T_bravo
#     data[(width - 1), j] = T_delta


# error_flag = True
# error_limit = 1e-4
# while error_flag:
#     large_error_term_found = False

#     # Gauss-Seidel Iteration
#     for n in range(1, (height - 1)):
#         for m in range(1, (width - 1)): 
#             data_old = data[n, m]
#             data[n, m] = 0.25 * (data[(n + 1), m] + data[(n - 1), m] + data[n, (m + 1)] + data[n, (m - 1)])

#             if not large_error_term_found:
#                 error_term = abs(data[n, m] - data_old) / data_old
#                 if (error_term <= error_limit):
#                     error_flag = False
#                 else:
#                     error_flag = True
#                     large_error_term_found = True

print(np.rot90(data))

data_printable = np.flipud(np.rot90(data))

figNum = 1
plt.figure(figNum)
plt.axes().set_aspect('equal')
plt.style.use('classic')
heatmap = plt.pcolor(data_printable)

plt.text(0.5, -0.02, "T = " + str(T_alpha) + "\N{DEGREE SIGN}C",
         horizontalalignment='center',
         verticalalignment='top',
         rotation=0,
         clip_on=False,
         transform=plt.gca().transAxes)
plt.text(0, 0.5, "T = " + str(T_bravo) + "\N{DEGREE SIGN}C",
         horizontalalignment='right',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=plt.gca().transAxes)
plt.text(0.5, 1, "T = " + str(T_charlie) + "\N{DEGREE SIGN}C",
         horizontalalignment='center',
         verticalalignment='bottom',
         rotation=0,
         clip_on=False,
         transform=plt.gca().transAxes)
plt.text(1, 0.5, "T = " + str(T_delta) + "\N{DEGREE SIGN}C",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=270,
         clip_on=False,
         transform=plt.gca().transAxes)

plt.axis("off")

plt.xlim(0, width)
plt.ylim(0, height)

cbar = plt.colorbar(heatmap)
cbar.set_label("Temperature (\N{DEGREE SIGN}C)")
plt.clim(0, 100)

plt.savefig(fileName + "/images/" + fileName + "-" + str(figNum) + ".png")
plt.show()
