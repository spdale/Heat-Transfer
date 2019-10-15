import matplotlib.pyplot as plt
import numpy as np

plt.style.use('classic')

data = np.random.rand(5, 4)
heatmap = plt.pcolor(data)

# for y in range(data.shape[0]):
#     for x in range(data.shape[1]):
#         plt.text(x + 0.5, y + 0.5, '%.4f' % data[y, x],
#                  horizontalalignment='center',
#                  verticalalignment='center',
#                  )

plt.colorbar(heatmap)

plt.show()