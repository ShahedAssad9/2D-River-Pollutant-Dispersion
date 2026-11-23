import numpy as np
import matplotlib.pyplot as plt

# Defining our problem
height = 100    # m
length = 200  # m
time = 900     # 15 minutes in seconds
dx = 2 # m
dy = 2 # m
dt = 0.1 # seconds
x_cordinate = 20 # x injection coordinate (m)
y_cordinate = 50 # y injection coordinate (m)
u = 1 # longitudinal velocity m/s
v = 0 # lateral velocity m/s
E = 5 # dispersion coefficient m^2/s

# Defining the equation coefficients (eps1, eps2, eps3, eps4)
eps1 = dt * ((1 / dt) - (u / dx) - 2 * E * ((1 / dx**2) + (1 / dy**2)))
eps2 = dt * ((u / dx) + (E / dx**2))
eps3 = (dt * E) / dx**2
eps4 = (dt * E) / dy**2

# Initialization
nodes_x = int(length / dx)
nodes_y = int(height / dy)
t_nodes = int(time / dt)

River_grid = np.zeros((nodes_x, nodes_y))  # river concentration is initially zero at all points

# Boundary Conditions
x_inj_idx = int(x_cordinate / dx)
y_inj_idx = int(y_cordinate / dy)
River_grid[x_inj_idx, y_inj_idx] = 100 # the concentration at injection point is 100



# Visualizing with a plot

x_length = np.linspace(0, length, nodes_x)  # 100 points along length (0 to 200)
y_width = np.linspace(0, height, nodes_y)   # 50 points along width (0 to 100)

fig, axis = plt.subplots(figsize=(10, 5))
pcm = axis.pcolormesh(x_length, y_width, River_grid.T, cmap=plt.cm.jet, vmin=0, vmax=100, shading='auto')
axis.set_aspect('equal')
axis.set_xlabel("River length [m]")
axis.set_ylabel("River width [m]")
plt.colorbar(pcm, ax=axis, label="Concentration [%]")

# Simulating
counter = 0

while counter < time:

    w = River_grid.copy()

    for i in range(1, nodes_x - 1):
        for j in range(1, nodes_y - 1):

            # Explicit  equation 
            River_grid[i, j] = (
                eps1 * w[i, j]
                + eps2 * w[i - 1, j]
                + eps3 * w[i + 1, j]
                + eps4 * (w[i, j - 1] + w[i, j + 1])
            )
    # Boundary Conditions:
    # 1. Neumann zero-gradient boundary conditions (dC/dx = 0 and dC/dy = 0) on all 4 outer walls
    River_grid[0, :] = River_grid[1, :]      # Upstream boundary:   C(0, j, t) = C(1, j, t)
    River_grid[-1, :] = River_grid[-2, :]    # Downstream boundary: C(i, j, t) = C(i-1, j, t)
    River_grid[:, 0] = River_grid[:, 1]      # Bottom bank:         C(i, -2, t) = C(i, -1, t)
    River_grid[:, -1] = River_grid[:, -2]    # Top bank:            C(i, 2, t) = C(i, 1, t)

    # 2. Dirichlet boundary condition at the continuous injection point (C = 100%)
    River_grid[x_inj_idx, y_inj_idx] = 100

    counter += dt

    print("t: {:.3f} [s], Average concentration: {:.2f} %".format(counter, np.average(River_grid)))

    # Updating the plot
    pcm.set_array(River_grid.T)
    axis.set_title("Distribution at t: {:.3f} [s].".format(counter))
    plt.pause(0.01)

# After the while loop finishes (at t = 900 s), plot Fig. 15.6 contour lines:
fig_contour, ax_c = plt.subplots(figsize=(10, 5))

levels = [35, 40, 45, 50, 55, 60, 70, 80, 90]
cs = ax_c.contour(x_length, y_width, River_grid.T, levels=levels, colors='black', linewidths=1.0)
ax_c.clabel(cs, inline=True, fontsize=10, fmt='%d')

# Mark the injection point
ax_c.plot(x_cordinate * dx, y_cordinate * dy, 'ko', markersize=5)
ax_c.text(x_cordinate * dx, y_cordinate * dy - 8, 'injection\npoint', ha='center')

ax_c.set_aspect('equal')
ax_c.set_xlabel("x (m)")
ax_c.set_ylabel("y (m)")
ax_c.set_title(f"Isoconcentration lines at t = {counter:.1f} s ({counter/3600:.2f} h)")
plt.show()


plt.show()