import numpy as np

a = np.zeros(3) + 1
b = np.arange(1, 4)

error_term = max(np.absolute(a - b) / b)
print(a)
print(b)
print(error_term)