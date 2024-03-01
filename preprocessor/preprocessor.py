import numpy as np

class Preprocessor:
    
    
    def __init__(self, dx=0.01, dy=0.01, plate_length=10, grid_x_size=12, grid_y_size=1, free_flow_velocity = 1):
        if grid_x_size < plate_length:
            raise ValueError("Grid x length must be greater than or equal to the plate's length.")

        self.dx = dx
        self.dy = dy
        self.plate_length = plate_length
        self.grid_x_size = grid_x_size
        self.grid_y_size = grid_y_size
        self.grid_nodes_x = None
        self.grid_nodes_y = None
        self.grid_u_velocity = None
        self.grid_v_velocity = None
        self.free_flow_velocity = free_flow_velocity

    def create_grid(self):
        x_nodes = int(self.grid_x_size / self.dx) + 1
        y_nodes = int(self.grid_y_size / self.dy) + 1

        x_coords = np.arange(0, self.grid_x_size + self.dx, self.dx)
        y_coords = np.arange(0, self.grid_y_size + self.dy, self.dy)

        self.grid_nodes_x, self.grid_nodes_y = np.meshgrid(x_coords, y_coords)

        return self.grid_nodes_x, self.grid_nodes_y
    
        
    def initialize_velocity_field(self):
        self.create_grid()   

    # Initialize velocity fields with zeros
        self.grid_u_velocity = np.ones_like(self.grid_nodes_x)
        self.grid_v_velocity = np.zeros_like(self.grid_nodes_y) * self.free_flow_velocity
        

    # Set velocity to 0 for points on the top boundary where grid_nodes_x is between 0 and 10
        for i in range(self.grid_u_velocity.shape[1]):
            mask = (self.grid_nodes_x[0, i] >= 1 - self.dx) & (self.grid_nodes_x[0, i] <= 10 +1)
            self.grid_u_velocity[0, i, mask] = 0

    # Set velocity to free_flow_velocity for the last row (bottom boundary)
        self.grid_u_velocity[-1, :] = self.free_flow_velocity
        
      
