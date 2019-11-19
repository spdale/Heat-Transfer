import numpy as np
import matplotlib.pyplot as plt

fileName = "Problem-Set-5"

num_time_steps = 100
delta = 1

# Grid points:
height = 30
width = 30

# Constants
rho = 3000  # kg/m^3
c = 840     # J/(kg*C)
h = 28      # W/(m^2*C)  Convective Heat Transfer Coefficient
k = 5.2     # W/(m*C)    Thermal Conductivity
alpha = k / (rho * c)

#dt = k * (delta ^ 2) / (2 * h * delta + 4 * k)  # Characteristic time
dt = (delta ^ 2) / (4 * alpha)
Fo = alpha * dt / (delta ^ 2)   # Fourier Number
Bi = h * delta / k              # Biot Number

T_initial = 10
T_right = 38
T_inf = 0

# Create array and initialize to T-initial
data = np.zeros((width, height)) + T_initial

# Set the right boundary to T_right
for j in range(height):
    data[(width - 1), j] = T_right


for t in range(num_time_steps):
    data_old = data.copy()

    # Internal Nodes
    for m in range(1, width - 1):
        for n in range(1, height - 1):
            data[m, n] = (data_old[m + 1, n] + data_old[m - 1, n] + data_old[m, n + 1] + data_old[m, n - 1]) / 4
    
    # Convective Boundary Nodes (Left)
    for n in range(1, height - 1):
        data[m, n] = Fo * (2 * Bi * (T_inf - data_old[m, n]) + 2 * data_old[m + 1, n] + data_old[m, n + 1] + data_old[m, n - 1] - 4 * data_old[m, n]) + data_old[m, n]

    # Insulated Boundary Noes (Top)
    for m in range(1, width - 1):
        data[m, 0] = Fo * (2 * data_old[m, n - 1] + data_old[m - 1, n] + data_old[m + 1, n]) + (1 - 4 * Fo) * data_old[m, n]

    



#print(np.rot90(data))

data_printable = np.rot90(data) #np.flipud(np.rot90(data))

figNum = 1
plt.figure(figNum)
plt.axes().set_aspect('equal')
plt.style.use('classic')
heatmap = plt.pcolor(data_printable)

plt.text(0.5, -0.02, "T = " + str(T_initial) + "\N{DEGREE SIGN}C",
         horizontalalignment='center',
         verticalalignment='top',
         rotation=0,
         clip_on=False,
         transform=plt.gca().transAxes)
plt.text(0, 0.5, "Convective Boundary",
         horizontalalignment='right',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=plt.gca().transAxes)
plt.text(0.5, 1, "Insulated Surface",
         horizontalalignment='center',
         verticalalignment='bottom',
         rotation=0,
         clip_on=False,
         transform=plt.gca().transAxes)
plt.text(1, 0.5, "T = " + str(T_right) + "\N{DEGREE SIGN}C",
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
plt.clim(0, 50)

plt.savefig(fileName + "/images/" + fileName + "-" + str(figNum) + ".png")
plt.show()
