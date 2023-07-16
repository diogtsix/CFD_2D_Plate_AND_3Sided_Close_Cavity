class Preprocessor:
    
    def __init__(self,dx,dy,plate_length = 10,grid_x_size = 12,grid_y_size = 10):
        
        
        if grid_x_size < plate_length:
            raise ValueError("Grid x length must be greater than or equal to the plate's length.")
        
        self.dx = dx
        self.dy = dy
        self.__plate_length = plate_length
        self.grid_x_size = grid_x_size
        self.grid_y_size = grid_y_size


