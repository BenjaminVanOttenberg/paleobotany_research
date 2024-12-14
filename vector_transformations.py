from matplotlib import pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection = "3d")

v1 = [1,2,3]
ax.plot3d(0, 0, 0, v1[0], v1[1], v1[2], color="r")
ax.set_xlim([-3, 3])
ax.set_ylim([-3 ,3])
ax.set_zlim([-3, 3])

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.title("3D Vector Plot")
plt.show()


