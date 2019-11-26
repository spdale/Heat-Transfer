import numpy as np
import matplotlib.pyplot as plt

fileName = "Problem-Set-6"

# Grid points:
height = 200
width = 500

cylinder_diameter = 50
cylinder_radius = cylinder_diameter / 2
cylinder_center = [(height / 2), 100]

error = 0.01    # 1% maximum change for convergence

# Create array and initialize to T-initial
omega = np.zeros((width, height))
psi = np.zeros((width, height))


# Set the right boundary to T_right
for i in range(width):
    for j in range(height):
        dist = np.sqrt((i - cylinder_center[0])**2 + (j - cylinder_center[1])**2)
        if (dist <= cylinder_radius):
            omega[i, j] = 100



# Print the data in the console (readable format)
#print(np.rot90(data))


figNum = 1
plt.figure(figNum)
plt.axes().set_aspect('equal')
plt.style.use('classic')
data_graphable = np.flipud(np.rot90(omega))
heatmap = plt.pcolor(data_graphable)

plt.axis("off")

plt.xlim(0, width)
plt.ylim(0, height)

plt.clim(np.amin(omega), np.amax(omega))

plt.savefig(fileName + "/images/" + fileName + "-Figure.png")
plt.show()
