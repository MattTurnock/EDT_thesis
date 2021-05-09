# import numpy as np
# import math
# from matplotlib import patches as patches_module
# # from matplotlib import pyplot as plt
# # from matplotlib import animation
# #
# # # Create figure
# # fig = plt.figure()
# # ax = fig.gca()
# #
# # # Axes labels and title are established
# # ax = fig.gca()
# # ax.set_xlabel('x')
# # ax.set_ylabel('y')
# #
# # ax.set_ylim(-2,2)
# # ax.set_xlim(-2,2)
# # plt.gca().set_aspect('equal', adjustable='box')
# #
# x = np.linspace(-1,1,20)
# y  = np.linspace(-1,1,20)
# dx = np.zeros(len(x))
# dy = np.zeros(len(y))
#
# for i in range(len(x)):
#     dx[i] = math.sin(x[i])
#     dy[i] = math.cos(y[i])
# arrowPatch = patches_module.Arrow(x[0], y[0], dx[0], dy[0] )
# #
# #
# # def init():
# #     ax.add_patch(patch)
# #     return patch,
# #
# # def animate(t):
# #     global patch
# #
# #     t %= 20 # get only 0-19 to loop animation and get color t/20 as 0.0-1.0
# #
# #     ax.patches.remove(patch)
# #
# #     patch = patches.Arrow(x[t], y[t], dx[t], dy[t])
# #
# #     patch.update({'facecolor': (t/20,t/20,t/20,1.0)})
# #
# #     ax.add_patch(patch)
# #
# #     return patch,
# #
# # anim = animation.FuncAnimation(fig, animate,
# #                                init_func=init,
# #                                interval=20,
# #                                blit=False)
# #
# # plt.show()
#
#
# import matplotlib
# matplotlib.use('Qt5Agg') #use Qt5 as backend, comment this line for default backend
#
# from matplotlib import pyplot as plt
# from matplotlib import animation
#
# fig = plt.figure()
#
# ax = plt.axes(xlim=(0, 2), ylim=(0, 100))
#
# N = 4
# lines = [plt.plot([], [])[0] for _ in range(N)] #lines to animate
#
# rectangles = plt.bar([0.5,1,1.5],[50,40,90],width=0.1) #rectangles to animate
#
# patches = lines + list(rectangles) + [arrowPatch] #things to animate
# # patches.append(arrowPatch)
#
# def init():
#     #init lines
#     for line in lines:
#         line.set_data([], [])
#
#     #init rectangles
#     for rectangle in rectangles:
#         rectangle.set_height(0)
#
#     patches = lines + list(rectangles) + [arrowPatch]
#     # print(patches)
#
#     return patches #return everything that must be updated
#
# def animate(i):
#     #animate lines
#     for j,line in enumerate(lines):
#         line.set_data([0, 2], [10 * j,i])
#
#     #animate rectangles
#     for j,rectangle in enumerate(rectangles):
#         rectangle.set_height(i/(j+1))
#
#
#     # arrowPatch = plt.Arrow(x[i], y[i], dx[i], dy[i] )
#     # ax.add_patch(arrowPatch)
#
#
#     # global arrowPatch
#     #
#     # i %= 20 # get only 0-19 to loop animation and get color t/20 as 0.0-1.0
#     #
#     # ax.patches.remove(arrowPatch)
#     #
#     # arrowPatch = patches_module.Arrow(x[i], y[i], dx[i], dy[i])
#     #
#     patches = lines + list(rectangles) + [arrowPatch]
#     #
#     # arrowPatch.update({'facecolor': (i/20,i/20,i/20,1.0)})
#     #
#     # ax.add_patch(arrowPatch)
#
#     return patches #return everything that must be updated
#
# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                frames=20, interval=20, blit=True)
#
# plt.show()












# constantX = [1]
# xCoords = np.array([constantX,]*Nframes)




import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.cm import get_cmap
from mpl_toolkits.mplot3d import Axes3D

# Use matplotlib ggplot stylesheet if available
try:
    plt.style.use('ggplot')
except OSError:
    pass

# Set which type of animation will be plotted. One of:
# quiver, 3d_contour, polar, scatter, fill
animation_type = 'quiver'

Nframes = 120



# Create an invisible axis
fig, ax = plt.subplots(figsize=(4, 3))
xlim = 4
ax.set(xlim=(-xlim, xlim), ylim=(-3, 3),
       xticklabels=[], yticklabels=[])

# x = np.linspace(0, 2*np.pi, 1, endpoint=False)
t = np.linspace(0, 2*np.pi, Nframes, endpoint=False)
# Arrays Qx and Qy are the positions of the base
# points of the arrow with size (Nframes, Narrows) where
# Narrows is the number of arrows
xCoords = np.linspace(-1,1,Nframes)
yCoords = np.linspace(-1,1,Nframes)
Qx = xCoords[:, None]
Qy = yCoords[:, None]

# Arrow vectors are the same shape as Qx and Qy
UCoords = np.linspace(-1,1,Nframes)
VCoords = np.repeat(np.linspace(1,-1,12), 10)
U = UCoords[:, None]
V = VCoords[:, None]

# For frame 1, plot the 0th set of arrows
s = np.s_[0, :]
qax = ax.quiver(Qx[s], Qy[s], U[s], V[s],
                facecolor="red", scale=10)

def animate(i):
    # Update to frame i
    s = np.s_[i, :]
    # Change direction of arrows
    qax.set_UVC(U[s], V[s])
    # Change base position of arrows
    qax.set_offsets(np.c_[Qx[s].flatten(), Qy[s].flatten()])



# ----------------------------------------------------------------------------
# Save the animation
anim = FuncAnimation(
    fig, animate, interval=50, frames=Nframes, repeat=False)
fig.show()
# anim.save(animation_type + '.gif', writer='imagemagick')
plt.show()