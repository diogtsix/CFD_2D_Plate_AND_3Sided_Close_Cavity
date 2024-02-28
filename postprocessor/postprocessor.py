import numpy as np
import matplotlib.pyplot as plt

class PostProcessor:
    def __init__(self, grid_nodes_x, grid_nodes_y, grid_u_velocity, 
                 solver_type):
        self.grid_nodes_x = grid_nodes_x
        self.grid_nodes_y = grid_nodes_y
        self.grid_u_velocity = grid_u_velocity
        self.solver_type = solver_type

    def plot_colored_velocity_field(self):
        plt.figure(figsize=(8, 6))
        
        mask = (self.grid_nodes_x >= 0.5) & (self.grid_nodes_x <= 10)
        
        # Apply mask to u_velocity, setting values outside the range to NaN so they won't be colored
        filtered_u_velocity = np.where(mask, self.grid_u_velocity, np.nan)
    
        plt.pcolormesh(self.grid_nodes_x, self.grid_nodes_y, filtered_u_velocity, cmap='viridis')

        plt.colorbar(label='u-velocity')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Colored Velocity Field (u-velocity)')
        plt.xlim(0, 10)  # Set the x-axis range to 0-11 meters
        plt.grid()
        plt.show()

    def plot_u_velocity_profiles(self, x_values):
        plt.figure(figsize=(8, 6))
        for x in x_values:
            idx = np.abs(self.grid_nodes_x[0, :] - x).argmin()
            
            if self.solver_type == "explicit":
                plt.plot(self.grid_u_velocity[:, idx], self.grid_nodes_y[:, idx], label=f'x = {x}')
            elif self.solver_type == "implicit":
                plt.plot( self.grid_nodes_y[:, idx],self.grid_u_velocity[:, idx], label=f'x = {x}')

        plt.xlabel('u-velocity')
        plt.ylabel('y')
        plt.title('u-velocity Profiles for Different x')
        plt.legend()
        plt.grid()
        plt.show()