import numpy as np
from solver_explicit import Solver_explicit
from preprocessor import Preprocessor
import matplotlib.pyplot as plt

PreProcess = Preprocessor()
PreProcess.create_grid()
PreProcess.initialize_velocity_field()

Result = Solver_explicit(PreProcess.grid_nodes_x,
                         PreProcess.grid_nodes_y,
                         PreProcess.grid_u_velocity,
                         PreProcess.grid_v_velocity,
                         PreProcess.dx, 
                         PreProcess.dy)


def visualize_matrix(matrix):
    shape = np.shape(matrix)
    rows, cols = shape[0], shape[1]

    # Flatten the matrix to get a 1D array of values
    values = matrix.flatten()

    # Generate coordinates for each point in the scatter plot
    x_coords, y_coords = np.meshgrid(range(cols), range(rows))

    # Plot the scatter plot
    plt.scatter(x_coords, y_coords, c=values, cmap='jet')
    plt.colorbar()

    # Set the axis limits and labels
    plt.xlim(0, cols-1)
    plt.ylim(0, rows-1)
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')

    # Show the plot
    plt.show()

Result.solve()

plt.quiver(PreProcess.grid_nodes_x, PreProcess.grid_nodes_y,
           Result.grid_u_velocity, Result.grid_v_velocity)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Velocity Field')
plt.show()

visualize_matrix(Result.grid_u_velocity)