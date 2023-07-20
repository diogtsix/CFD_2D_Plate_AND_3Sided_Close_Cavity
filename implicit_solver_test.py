import numpy as np
from preprocessor import Preprocessor
from solver_implicit import Solver_implicit
import matplotlib.pyplot as plt

# Create a preprocessor instance and initialize the velocity field
preprocessor = Preprocessor(dx=0.1, dy=0.1, plate_length=0.5, grid_x_size=1, grid_y_size=1)
preprocessor.initialize_velocity_field()

# Create a solver instance and solve the system
solver = Solver_implicit(preprocessor.grid_nodes_x, preprocessor.grid_nodes_y,
                         preprocessor.grid_u_velocity, preprocessor.grid_v_velocity,
                         preprocessor.dx, preprocessor.dy)
solver.u_i_plus_1_system_eq()

# Visualize the velocity fields
fig, axs = plt.subplots(2, 1, figsize=(6, 8))
sc1 = axs[0].scatter(preprocessor.grid_nodes_x.flatten(), preprocessor.grid_nodes_y.flatten(),
                     c=solver.grid_u_velocity.flatten(), cmap='coolwarm')
axs[0].set_title('Grid U-Velocity')
axs[0].set_xlabel('X')
axs[0].set_ylabel('Y')
fig.colorbar(sc1, ax=axs[0])

sc2 = axs[1].scatter(preprocessor.grid_nodes_x.flatten(), preprocessor.grid_nodes_y.flatten(),
                     c=solver.grid_v_velocity.flatten(), cmap='coolwarm')
axs[1].set_title('Grid V-Velocity')
axs[1].set_xlabel('X')
axs[1].set_ylabel('Y')
fig.colorbar(sc2, ax=axs[1])

plt.tight_layout()
plt.show()