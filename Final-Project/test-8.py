import numpy as np

a = np.zeros((5, 5))

count = 1
for i in range(5):
    for j in range(5):
        a[i][j] = count
        count += 1
print(a)
print()

indices = np.nonzero((a % 3 == 0) & (a % 5 != 0))

# print(indices)
# print()

# print(a[indices])
# print()


# indices[:][1][:] += 1
# print(indices)
# print()
indices_right = indices[:]
indices_right[:][1][:] += 1



print(a[indices_right])


# print("\nAlso Desired:")
# print("[ 4. 7. 13. 19. 21. 25.]")
# print("Using a[indices] of some form")