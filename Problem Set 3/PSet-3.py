import numpy as np
import matplotlib.pyplot as plt 

plt.style.use('classic')

height = 21
width = 31

delta = dx = dy = 0.1

x = np.linspace(0, delta, width)
x = np.zeros()

I = np.sin(x) * np.cos(x[:, np.newaxis])

plt.imshow(I)
#plt.colorbar()
plt.colorbar().set_label("Temperature (\N{DEGREE SIGN}C)")
# cbar = plt.colorbar()
# cbar.ax.get_yaxis().labelpad = 15
# cbar.ax.set_ylabel('# of contacts', rotation=270)
plt.clim(0, 100)

plt.show()

# fig = plt.figure()
# im = ax.imshow(u.copy(), cmap=plt.get_cmap('hot'), vmin=Tcool,vmax=Thot)
# cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
# cbar_ax.set_xlabel('$T$ / K', labelpad=20)
# fig.colorbar(im, cax=cbar_ax)
# plt.show()