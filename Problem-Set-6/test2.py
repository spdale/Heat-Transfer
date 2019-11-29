import numpy as np
import matplotlib.pyplot as plt

data_graphable = np.linspace(-20, 20)

num_streamlines = 15
max_streamline = np.max(data_graphable)
min_streamline = np.min(data_graphable)
contours_before = np.linspace(min_streamline, max_streamline, num=(num_streamlines + 3))
contours = contours_before[(contours_before != 0) & (contours_before != min_streamline) & (contours_before != max_streamline)]

print(contours)