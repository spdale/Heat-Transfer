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

plt.text(0.5, -0.05, "T = " + str(T_alpha) + "\N{DEGREE SIGN}C",
         horizontalalignment='center',
         verticalalignment='center',
         rotation=0,
         clip_on=False,
         transform=plt.gca().transAxes)
plt.text(-0.05, 0.5, "T = " + str(T_bravo) + "\N{DEGREE SIGN}C",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=plt.gca().transAxes)
plt.text(0.5, 1.02, "T = " + str(T_charlie) + "\N{DEGREE SIGN}C",
         horizontalalignment='center',
         verticalalignment='center',
         rotation=0,
         clip_on=False,
         transform=plt.gca().transAxes)
plt.text(1, 0.5, "T = " + str(T_delta) + "\N{DEGREE SIGN}C",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=270,
         clip_on=False,
         transform=plt.gca().transAxes)

# plt.text(1.02, 0.5, r"\bf{level set} $\phi$", {'color': 'C2', 'fontsize': 20},
#          horizontalalignment='left',
#          verticalalignment='center',
#          rotation=90,
#          clip_on=False,
#          transform=plt.gca().transAxes)

plt.axis("off")
# ax = plt.axes()
# ax.xaxis.set_ticks_position('none')
# ax.yaxis.set_ticks_position('none')
# ax.yaxis.set_major_locator(plt.NullLocator())
# ax.xaxis.set_major_formatter(plt.NullFormatter())
# ax.set_xlabel("T = " + str(T_alpha) + "\N{DEGREE SIGN}C")
# ax.set_ylabel("T = " + str(T_bravo) + "\N{DEGREE SIGN}C")

# ax2 = ax.twinx()
# ax2.yaxis.set_label_position('right')
# ax2.yaxis.set_ticks_position('none')
# ax2.yaxis.set_major_locator(plt.NullLocator())
# ax2.set_ylabel("T = " + str(T_delta) + "\N{DEGREE SIGN}C")

# ax3 = ax.twiny()
# ax3.xaxis.set_label_position('top')
# ax3.xaxis.set_ticks_position('none')
# ax3.xaxis.set_major_formatter(plt.NullFormatter())
# ax3.set_xlabel("T = " + str(T_charlie) + "\N{DEGREE SIGN}C")

plt.xlim(0, width)
plt.ylim(0, height)

cbar = plt.colorbar(heatmap)
cbar.set_label("Temperature (\N{DEGREE SIGN}C)")
plt.clim(0, 100)

plt.show()
