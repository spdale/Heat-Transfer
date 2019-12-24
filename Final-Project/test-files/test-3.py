import numpy as np
import matplotlib.pyplot as plt

t = np.arange(0.0, 20.0, 1)
s = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

# fig = plt.figure()

fig, (sub1, sub2) = plt.subplots(2, 1)
# sub1 = fig.add_subplot(211)
sub1.plot(t, s)
sub1.set_aspect('equal')
# sub1.aspect('equal')

# sub2 = fig.add_subplot(212)
sub2.plot(t, s)


# plt.axes().set_aspect('equal')
plt.show()