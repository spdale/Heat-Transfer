import numpy as np

a = np.array([1, 1, 0, 1, 1])

b = np.zeros(5) + 2

print(a)
print(b)

b[1:-1, 1:-1] = 1

print(b)