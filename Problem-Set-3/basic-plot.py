import numpy as np
import matplotlib.pyplot as plt 

plt.axes().set_aspect('equal')
plt.style.use('classic')

# Grid squares values:
height = 21
width = 31

# Initialize matrix of zeros for that size
data = np.zeros((height, width)) + 75

# Make set temperatures on fixed positions
T_alpha   = 0    # (bottom boundary temperature)
T_bravo   = 40   # (left boundary temperature)
T_charlie = 100  # (top boundary temperature)
T_delta   = 100  # (right boundary temperature)

# Note: index 0,0 is bottom left

for i in range(width):
    data[0, i] = T_alpha
    data[(height - 1), i] = T_charlie
for j in range(1, (height - 1)):
    data[j, 0] = T_bravo
    data[j, (width - 1)] = T_delta

#print(data)

heatmap = plt.pcolor(data)

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

plt.show()
