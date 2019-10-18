import numpy as np
import matplotlib.pyplot as plt 

plt.style.use('classic')

# Grid squares values:
height = 21
width = 31

# Initialize matrix of zeros for that size
data = np.zeros((height, width)) + 75
#data = np.random.randint(100, size=(width, height))

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

fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
#ax = plt.axes()
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.yaxis.set_major_locator(plt.NullLocator())
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.set_xlabel("T = " + str(T_alpha) + "\N{DEGREE SIGN}C")
ax.set_ylabel("T = " + str(T_bravo) + "\N{DEGREE SIGN}C")

#ax2 = ax.twinx()
ax2.yaxis.set_label_position('right')
ax2.yaxis.set_ticks_position('none')
ax2.yaxis.set_major_locator(plt.NullLocator())
ax2.set_ylabel("T = " + str(T_delta) + "\N{DEGREE SIGN}C")

# ax3 = ax.twiny()
# ax3.xaxis.set_label_position('top')
# ax3.xaxis.set_ticks_position('none')
# ax3.xaxis.set_major_formatter(plt.NullFormatter())
# ax3.set_xlabel("T = " + str(T_charlie) + "\N{DEGREE SIGN}C")

# ax2.set_ylim(ax.get_ylim())
# #ax3.set_xlim(ax.get_xlim())
# ax.set_adjustable('box')
# ax2.set_adjustable('box')

plt.xlim(0, width)
plt.ylim(0, height)

cbar = plt.colorbar(heatmap)
cbar.set_label("Temperature (\N{DEGREE SIGN}C)")
plt.clim(0, 100)

plt.show()
