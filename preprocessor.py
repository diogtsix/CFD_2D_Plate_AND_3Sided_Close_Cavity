import numpy as np

class Preprocessor:
    
    
    def __init__(self, dx=0.001, dy=0.001, plate_length=10, grid_x_size=12, grid_y_size=1):
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
        self.grid_u_velocity = np.zeros_like(self.grid_nodes_x)
        self.grid_v_velocity = np.zeros_like(self.grid_nodes_y)
        
        # Set boundary conditions
        self.grid_u_velocity[:, 0:1] = 1.0  # First 2 columns, all rows
        self.grid_u_velocity[0, :] = 0.0  # First row, all columns
        