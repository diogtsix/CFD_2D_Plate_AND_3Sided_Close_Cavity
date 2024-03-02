import numpy as np

class cavity_preprocessor:
    
    
    def __init__(self, wall_length, dx, dy, dt, relax_factor, wall_vel, 
                 reynolds, error):
        
        
        self.wall_length = wall_length
        self.dx = dx
        self.dy = dy
        self.dt = dt
        self.relax_factor = relax_factor
        self.wall_vel = wall_vel
        self.reynolds = reynolds
        self.error = error
        
        self.nodx = None 
        self.nody = None 
        self.xCoords = None
        self.yCoords = None
        self.u_vel = None
        self.v_vel = None
        self.psi1 = None
        self.psi2 = None
        self.psi3 = None
        self.zhta1 = None
        self.zhta2 = None
        self.zhta3 = None
        
        
    def create_grid(self):
        
        # Initialize matrices and constants
        self.nodx = round(self.wall_length / self.dx) + 1
        self.nody = round(self.wall_length / self.dy) + 1
        
        self.xCoords = np.linspace(0, self.wall_length  ,self.nodx) # coordinates as a vector because of symmetry
        self.yCoords = np.linspace(0, self.wall_length  ,self.nody) # coordinates as a vector because of symmetry
        
        self.u_vel = np.zeros((self.nody, self.nodx))
        self.v_vel = np.zeros((self.nody, self.nodx))
        self.psi1 = np.zeros((self.nody, self.nodx))
        self.psi2 = np.zeros((self.nody, self.nodx))
        self.psi3 = np.zeros((self.nody, self.nodx))
        self.zhta1 = np.zeros((self.nody, self.nodx))
        self.zhta2 = np.zeros((self.nody, self.nodx))
        self.zhta3 = np.zeros((self.nody, self.nodx))
    
    
    def initialize_velocity_field(self):
        
        self.create_grid()   
        # Set top wall velocity 
        self.u_vel[-1,:] = self.wall_vel
        