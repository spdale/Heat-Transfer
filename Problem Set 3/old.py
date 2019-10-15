#%%
import numpy as np
import matplotlib.pyplot as plt 

plt.style.use('classic')

harvest = np.array([[0.8, 2.4, 2.5, 3.9, 0.0, 4.0, 0.0],
                    [2.4, 0.0, 4.0, 1.0, 2.7, 0.0, 0.0],
                    [1.1, 2.4, 0.8, 4.3, 1.9, 4.4, 0.0],
                    [0.6, 0.0, 0.3, 0.0, 3.1, 0.0, 0.0],
                    [0.7, 1.7, 0.6, 2.6, 2.2, 6.2, 0.0],
                    [1.3, 1.2, 0.0, 0.0, 0.0, 75, 5.1],
                    [0.1, 2.0, 0.0, 1.4, 0.0, 1.9, 6.3]])

height = 21
width = 31

delta = dx = dy = 0.1

x = np.linspace(0, delta, width)
x = np.zeros(width)
y = np.zeros(height)

fig, ax = plt.subplots()


plt.imshow(harvest)
ax.set_xticks(np.arange(7))
ax.set_yticks(np.arange(7))
#plt.colorbar()
plt.colorbar().set_label("Temperature (\N{DEGREE SIGN}C)")
# cbar = plt.colorbar()
# cbar.ax.get_yaxis().labelpad = 15
# cbar.ax.set_ylabel('# of contacts', rotation=270)
#plt.clim(0, 100)

plt.show()

# fig = plt.figure()
# im = ax.imshow(u.copy(), cmap=plt.get_cmap('hot'), vmin=Tcool,vmax=Thot)
# cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
# cbar_ax.set_xlabel('$T$ / K', labelpad=20)
# fig.colorbar(im, cax=cbar_ax)
# plt.show()

#%%
